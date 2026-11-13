"""Indel lengths: the counting (phylobot.indels) and the two pages that show it."""

import json

import pytest
from django.core.cache import cache
from django.template.loader import render_to_string
from django.test import Client

from phylobot import gaps, indel_views, indels, jobs
from phylobot.models import Job
from phylobot.modes import local_user
from tests.conftest import FakeStorage
from tests.test_views import viewer_job


@pytest.fixture(autouse=True)
def empty_cache() -> None:
    cache.clear()


def presence(text: str) -> list[bool]:
    return [ch == "x" for ch in text]


def lengths(parent: str, child: str) -> list[tuple[int, bool]]:
    found = indels.branch_events(presence(parent), presence(child))
    return [(e.length, e.terminal) for e in found]


def test_a_run_on_a_branch_is_one_event() -> None:
    assert lengths("xxxxxxxx", "xx---xxx") == [(3, False)]
    assert lengths("xx---xxx", "xxxxxxxx") == [(3, False)]  # gained, same length
    assert lengths("xxxxxxxx", "xxxxxxxx") == []


def test_a_change_of_direction_is_two_events() -> None:
    assert lengths("xxx--xxxx", "xxxxx--xx") == [(2, False), (2, False)]
    assert lengths("xxxx--xxx", "xx--xxxxx") == [(2, False), (2, False)]


def test_columns_absent_at_both_ends_do_not_split_or_end_a_run() -> None:
    assert lengths("xx-xx-xxx", "x----xxxx")[0] == (3, False)  # bridges one absent column
    assert lengths("--xxxxx--", "--xx-xx--") == [(1, False)]  # not an end of the alignment


def test_runs_reaching_either_end_are_terminal() -> None:
    assert lengths("xxxxxxxx", "---xxxxx") == [(3, True)]
    assert lengths("xxxxxxxx", "xxxxx---") == [(3, True)]
    assert lengths("--xxxxxx", "----xxxx") == [(2, True)]  # the end of what the branch holds
    assert lengths("xxxxxxxx", "-xx--xx-") == [(1, True), (2, False), (1, True)]


def test_on_the_tree_a_shared_gap_is_one_event() -> None:
    rows = {"a": "MK--KL", "b": "MK--KL", "c": "MKRRKL", "d": "MKRRKL"}
    found = indels.on_tree("((a:1,b:1)N2:1,(c:1,d:1)N3:1)N1;", rows)
    assert found.counts == {2: 1} and found.end_gaps == 0
    assert found.on_tree and found.rate == pytest.approx(1 / found.exposure)


def test_end_gaps_are_counted_apart() -> None:
    rows = {"a": "---KLM", "b": "MKRKLM", "c": "MKRKLM", "d": "MKRKLM"}
    found = indels.on_tree("((a:1,b:1)N2:1,(c:1,d:1)N3:1)N1;", rows)
    assert found.counts == {} and found.end_gaps == 1


def test_an_unlabelled_root_still_counts_its_branches() -> None:
    rows = {"a": "MK--KL", "b": "MK--KL", "c": "MKRRKL", "d": "MKRRKL"}
    assert indels.on_tree("((a:1,b:1):1,(c:1,d:1):1);", rows).counts == {2: 1}


def test_a_tip_missing_from_the_alignment_is_a_key_error() -> None:
    with pytest.raises(KeyError):
        indels.on_tree("(a:1,b:1,z:1)N1;", {"a": "MK", "b": "M-"})


def test_without_a_tree_each_distinct_interior_gap_counts_once() -> None:
    rows = {"a": "MK--KL", "b": "MK--KL", "c": "M---KL", "d": "--RRKL", "e": "MKRR--"}
    found = indels.in_alignment(rows)
    assert found.counts == {2: 1, 3: 1} and found.end_gaps == 2 and not found.on_tree


def test_summary_figures() -> None:
    found = indels.Lengths({1: 3, 2: 1, 10: 1}, 0)
    assert (found.events, found.mean, found.median, found.longest) == (5, 3.0, 1, 10)
    assert found.binned(2) == [3, 1, 1]
    assert indels.Lengths({}, 0).mean == 0 and indels.Lengths({}, 0).median == 0


def test_the_geometric_sums_to_one_and_has_its_mean() -> None:
    bins = indels.geometric(3.0, 40)
    assert sum(bins) == pytest.approx(1)
    assert sum(k * p for k, p in enumerate(bins[:-1], 1)) == pytest.approx(3.0, abs=0.01)
    assert indels.geometric(1.0, 3) == [1.0, 0.0, 0.0, 0.0]


def test_fitch_presence_reports_labelled_nodes_root_first() -> None:
    rows = {"a": "MK-", "b": "M-K", "c": "MKK"}
    found = gaps.fitch_presence("((a,b)Node2,c)Node1;", rows)
    assert list(found) == ["Node1", "Node2"]
    assert found == {"Node1": [True] * 3, "Node2": [True] * 3}


def measured(aligner: str, counts: dict[int, int]) -> indel_views.Measured:
    return indel_views.Measured(aligner, aligner, indels.Lengths(counts, 0, 1.0), "x/y/z")


def test_every_thumbnail_has_the_same_axes() -> None:
    """Whatever the lengths in each alignment, and on any job."""
    short = measured("short", {1: 9, 2: 1})
    long = measured("long", {3: 1, 19: 2, 300: 7})
    empty = indel_views.Measured("none", "none", indels.Lengths({}, 0), "")
    row = indel_views.thumbnails([short, long, empty])
    other_job = indel_views.thumbnails([measured("other", {7: 1})])
    ticks = {tuple(c.x_ticks) for c in [*row, *other_job]}
    assert len(ticks) == 1
    assert [text for _, text in next(iter(ticks))] == ["1", "5", "10", "15", ">20"]
    assert len({tuple(c.y_ticks) for c in row}) == 1  # a row shares its y scale too


def test_lengths_past_the_axis_share_the_last_bar() -> None:
    chart = indel_views.thumbnails([measured("a", {3: 1, 21: 1, 300: 2})])[0]
    assert [b.title for b in chart.bars] == ["Length 3: 1 (25%)", "Longer than 20: 3 (75%)"]
    assert [b.longer for b in chart.bars] == [False, True]


def test_bars_stay_inside_the_plot() -> None:
    for chart in indel_views.thumbnails([measured("a", {1: 1}), measured("b", {2: 3, 30: 1})]):
        for bar in chart.bars:
            assert chart.top <= bar.y and bar.y + bar.height <= chart.bottom + 1e-6
            assert chart.left <= bar.x and bar.x + bar.width <= chart.right


def test_the_large_chart_draws_the_geometric() -> None:
    chart = indel_views.large(measured("a", {1: 5, 2: 3, 4: 1}))
    assert chart.expected and len(chart.expected.split()) == indel_views.CHART_LAST + 1
    assert not indel_views.large(measured("b", {1: 5})).expected  # mean 1: nothing to compare
    svg = render_to_string("phylobot/_indel_histogram.html", {"c": chart, "small": False})
    assert 'stroke-dasharray="5 4"' in svg and "Indel length (residues)" in svg


def test_tsv_keeps_every_length_past_the_axis() -> None:
    text = indel_views.tsv(measured("a", {2: 1, 300: 1}))
    assert "300\t1" in text and "on tree x/y/z" in text


def test_the_best_tree_is_ranked_by_bic_and_rooted_when_it_can_be() -> None:
    combos = [
        {"key": "m/a/LG", "aligner": "m", "tree": "a", "model": "LG", "rank_bic": 2},
        {"key": "m/b/WAG", "aligner": "m", "tree": "b", "model": "WAG", "rank_bic": 1},
        {"key": "o/a/LG", "aligner": "o", "tree": "a", "model": "LG", "rank_bic": 1},
    ]
    keys = {"m/a/LG/tree.nwk", "m/b/WAG/tree.nwk", "m/b/WAG/tree.rooted.nwk", "o/a/LG/tree.nwk"}
    summary = {"combinations": combos}
    assert indel_views.best_tree(summary, keys, "m") == ("m/b/WAG", "m/b/WAG/tree.rooted.nwk")
    assert indel_views.best_tree(summary, keys, "o") == ("o/a/LG", "o/a/LG/tree.nwk")
    assert indel_views.best_tree(summary, keys, "none") == ("", "")


def test_a_long_end_gap_stays_out_of_the_mean() -> None:
    """Why SR238's mean is 2.0, not 8.9."""
    rows = {"a": "-----MKLA", "b": "MKRSTMKLA", "c": "MKR--MKLA", "d": "MKRSTMKLA"}
    found = indels.on_tree("((a:1,b:1)N2:1,(c:1,d:1)N3:1)N1;", rows)
    assert found.mean == 2 and found.end_gaps == 1


@pytest.mark.django_db
class TestPages:
    def test_job_page_shows_a_thumbnail_per_alignment(
        self, client: Client, fake_storage: FakeStorage
    ) -> None:
        job = viewer_job(fake_storage)
        page = client.get(f"/jobs/{job.id}/").content.decode()
        assert page.count('class="plot indel-thumb"') == 2
        assert page.index("?aligner=mafft#indels") < page.index("?aligner=muscle#indels")
        assert "Indel lengths" in page
        assert f"/jobs/{job.id}/alignment/?aligner=mafft#indels" in page.replace("&amp;", "&")

    def test_viewer_shows_the_chart_and_its_counts(
        self, client: Client, fake_storage: FakeStorage
    ) -> None:
        job = viewer_job(fake_storage)
        page = client.get(f"/jobs/{job.id}/alignment/").content.decode()
        assert 'id="indels"' in page and "end gap" in page
        reply = client.get(f"/jobs/{job.id}/alignment/?export=indels")
        assert reply["Content-Disposition"].endswith('.indel-lengths.tsv"')
        assert b"length\tevents" in reply.content

    def test_the_chart_swaps_with_the_alignment(
        self, client: Client, fake_storage: FakeStorage
    ) -> None:
        job = viewer_job(fake_storage)
        page = client.get(
            f"/jobs/{job.id}/alignment/?aligner=muscle", headers={"HX-Request": "true"}
        ).content.decode()
        assert 'id="indels"' in page and "<html" not in page

    def test_an_alignments_only_job_counts_gaps_and_says_so(
        self, client: Client, fake_storage: FakeStorage
    ) -> None:
        job = Job.objects.create(owner=local_user(), params={"stop_after": "alignments"})
        assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
        prefix = f"jobs/{job.id}/results/"
        entry = {"key": "mafft", "label": "MAFFT", "path": "mafft/alignment.fasta"}
        summary = {"alignments": [entry], "combinations": []}
        fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
        fake_storage.objects[prefix + "mafft/alignment.fasta"] = b">a\nMK--KL\n>b\nMKRRKL\n"
        page = client.get(f"/jobs/{job.id}/").content.decode()
        assert page.count('class="plot indel-thumb"') == 1 and "no tree" in page
        assert "each distinct gap counts once" not in page  # the caption's "no tree" says it

    def test_an_unreadable_tree_never_breaks_the_page(
        self, client: Client, fake_storage: FakeStorage
    ) -> None:
        job = viewer_job(fake_storage)
        prefix = f"jobs/{job.id}/results/"
        fake_storage.objects[prefix + "mafft/iqtree/modelfinder/tree.nwk"] = b"not a tree"
        assert client.get(f"/jobs/{job.id}/").status_code == 200
