"""The same ancestor across combinations (phylobot.analogs, phylobot.analog_views): the
matching, the tree viewer's popovers, the ancestor viewer's comparison and "Add all"."""

import json

import pytest
from django.test import Client

from phylobot import analogs, jobs
from phylobot.models import Job, Library
from phylobot.modes import local_user
from phylobot.pipeline import ancestor_map
from tests.conftest import FakeStorage

LG, WAG, MUSCLE = "mafft/iqtree/LG", "mafft/iqtree/WAG", "muscle/iqtree/LG"
# Outgroup e. MAFFT's two trees agree, named differently; MUSCLE's resolves the ingroup
# as ((a,b),c),d: the same ingroup clade (Node9) by another history, and no (c,d) clade.
TREES = {
    LG: ("(e:0.1,(a:0.1,b:0.1)Node2:0.1,(c:0.1,d:0.1)Node3:0.1)Node1;",
         "(e:0.05,((a:0.1,b:0.1)Node2:0.1,(c:0.1,d:0.1)Node3:0.1)Node1:0.05);"),
    WAG: ("(e:0.1,(c:0.1,d:0.1)Node7:0.1,(a:0.1,b:0.1)Node8:0.1)Node6;",
          "(e:0.05,((c:0.1,d:0.1)Node7:0.1,(a:0.1,b:0.1)Node8:0.1)Node6:0.05);"),
    MUSCLE: ("(e:0.1,((a:0.1,b:0.1)Node5:0.1,c:0.1)Node4:0.1,d:0.1)Node9;",
             "(e:0.05,(((a:0.1,b:0.1)Node5:0.1,c:0.1)Node4:0.1,d:0.1)Node9:0.05);"),
}  # fmt: skip
ALIGNMENTS = {
    "mafft": {"a": "MKL", "b": "MKL", "c": "MRL", "d": "MRL", "e": "MKV"},
    "muscle": {"a": "MK-L", "b": "MK-L", "c": "MR-L", "d": "MR-L", "e": "MKV-"},
}
ANCESTORS = {
    LG: {"Node1": "MKL", "Node2": "MKL", "Node3": "MRL"},
    WAG: {"Node6": "MKL", "Node7": "MRL", "Node8": "MKL"},
    MUSCLE: {"Node9": "MR-L", "Node4": "MK-L", "Node5": "MK-L"},
}


def places() -> dict[str, dict[str, analogs.Place]]:
    return {key: analogs.places(rooted) for key, (_, rooted) in TREES.items()}


def test_matched_by_clade_flagging_a_different_history() -> None:
    found = analogs.index(places(), LG, by_clade=True)
    assert found["Node1"] == [
        analogs.Analog(WAG, "Node6", alike=True),
        analogs.Analog(MUSCLE, "Node9", alike=False),  # same clade, other branching inside
    ]
    assert found["Node2"] == [
        analogs.Analog(WAG, "Node8", True),
        analogs.Analog(MUSCLE, "Node5", True),
    ]
    assert found["Node3"] == [analogs.Analog(WAG, "Node7", True)]
    assert analogs.missing(places(), LG, found["Node3"]) == [MUSCLE]
    assert analogs.find(places(), LG, "NodeX", by_clade=True) == []


def test_unrooted_trees_match_by_partition_wherever_their_files_root_them() -> None:
    unrooted = {key: analogs.places(newick) for key, (newick, _) in TREES.items()}
    # The same tree rooted elsewhere: its nodes still match, by the subtrees meeting there.
    unrooted["other"] = analogs.places(
        "((a:0.1,b:0.1)NodeB:0.1,c:0.1,(d:0.1,e:0.1)NodeA:0.1)NodeC;"
    )
    found = analogs.index(unrooted, LG, by_clade=False)
    assert analogs.Analog("other", "NodeB", True) in found["Node2"]
    # Node1 (a,b | c,d | e) has no counterpart in MUSCLE's tree: it needs the same
    # branching on every side, which is the price of having no root.
    assert [a.combination for a in found["Node1"]] == [WAG]


def map_rows(rooted: bool) -> list[analogs.MapRow]:
    """ancestor_map.tsv for TREES, as ANCESTOR_MAP writes and the viewers read it."""
    lines = [
        f"{key}\t{'rooted' if rooted else 'unrooted'}\t{pair[1] if rooted else pair[0]}"
        for key, pair in TREES.items()
    ]
    return analogs.read_map(ancestor_map.build("\n".join(lines)))


@pytest.mark.parametrize("by_clade", [True, False])
def test_the_map_matches_exactly_as_the_trees_do(by_clade: bool) -> None:
    """The IDs stand in for the leaf sets: every node's analogs, alike flags included, are
    those of matching on the parsed trees."""
    trees = {
        key: analogs.places(rooted if by_clade else newick)
        for key, (newick, rooted) in TREES.items()
    }
    mapped = analogs.placed(map_rows(by_clade), by_clade)
    for combination in TREES:
        assert analogs.index(mapped, combination, by_clade) == analogs.index(
            trees, combination, by_clade
        )


def test_map_ids_name_the_sequences_not_their_order() -> None:
    """The same tree written another way round (input reordered, "Duplicate and run
    again") gives the same IDs; another clade gives another."""
    one = analogs.map_rows("x", "(e,((a,b)N2,(c,d)N3)N1);", rooted=True)
    two = analogs.map_rows("y", "(((d,c)M3,(b,a)M2)M1,e);", rooted=True)
    assert {r.place for r in one} == {r.place for r in two}
    assert analogs.clade_id(["a", "b"]) == analogs.clade_id(["b", "a"]) != analogs.clade_id(["a"])
    assert {r.place.leaves for r in one} == {2, 4}
    unrooted = analogs.map_rows("x", "(e,(a,b)N2,(c,d)N3)N1;", rooted=False)
    assert all(r.place.below == "" for r in unrooted)  # no clades without a root


def test_the_map_reads_back_what_it_wrote() -> None:
    rows = map_rows(rooted=True)
    text = analogs.write_map(rows)
    assert text.splitlines()[0] == "combination\tnode\tclade\tparts\tleaves"
    assert sorted(analogs.read_map(text), key=str) == sorted(rows, key=str)
    with pytest.raises(ValueError):
        analogs.read_map("node\tclade\n")


def test_identity_counts_against_the_longer_ancestor() -> None:
    assert analogs.identity("MKL", "MKL", 3, 3).percent == 100.0
    assert analogs.identity("MK-L", "MR-L", 3, 3).percent == pytest.approx(66.7)
    # A residue the frame could not draw still lengthens the ancestor it belongs to.
    assert analogs.identity("MKL", "MKL", 3, 4).fraction == 0.75
    assert analogs.identity("", "", 0, 0).fraction == 0.0


def finished_job(fake_storage: FakeStorage) -> Job:
    job = Job.objects.create(
        owner=local_user(), name="analogs", params={"outgroup": ["e"], "gap_methods": []}
    )
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combos = []
    for key, (newick, rooted) in TREES.items():
        aligner, tree, model = key.split("/")
        combos.append(
            {"key": key, "aligner": aligner, "tree": tree, "model": model,
             "aligner_label": aligner.upper(), "model_label": model, "rank_bic": 1,
             "label": f"{aligner.upper()} · IQ-TREE · {model}"}
        )  # fmt: skip
        rows = ANCESTORS[key]
        posteriors = "node\tsite\tK\tL\tM\tR\n" + "".join(
            f"{node}\t{site}\t" + "\t".join("1.0" if s == c else "0.0" for s in "KLMR") + "\n"
            for node, seq in rows.items()
            for site, c in enumerate(seq, start=1)
        )
        files = {
            "tree.nwk": newick,
            "tree.rooted.nwk": rooted,
            "alignment.fasta": "".join(f">{n}\n{s}\n" for n, s in ALIGNMENTS[aligner].items()),
            "ancestors.fasta": "".join(f">{n}\n{s}\n" for n, s in rows.items()),
            "ancestors_posteriors.tsv": posteriors,
        }
        for name, text in files.items():
            fake_storage.objects[f"{prefix}{key}/{name}"] = text.encode()
    summary = {"primary": LG, "combinations": combos}
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    fake_storage.map_ancestors(job)
    return job


@pytest.mark.django_db
def test_tree_viewer_popovers_list_the_same_ancestor_elsewhere(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    page = client.get(f"/jobs/{job.id}/tree/", {"combination": LG}).content.decode()
    assert 'class="tree-node" data-node="Node1"' in page and "phylobot/tree.js" in page
    popover = page.split('class="node-popover" data-node="Node1"', 1)[1].split("</div>", 1)[0]
    assert "The same ancestor in 2 of 2 other combinations:" in popover
    assert "<strong>Node6</strong>" in popover and "<strong>Node9</strong>" in popover
    assert popover.count('class="pop-flag"') == 1  # Node9: same clade, another history
    assert "#analogs" not in popover  # no "Compare them" link (owner, 2026-10-30)
    lacking = page.split('class="node-popover" data-node="Node3"', 1)[1].split("</div>", 1)[0]
    assert "1 of 2" in lacking and "Not in MUSCLE · LG: its tree lacks this clade." in lacking


@pytest.mark.django_db
def test_ancestor_viewer_compares_and_measures_identity(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    page = client.get(
        f"/jobs/{job.id}/ancestors/", {"combination": LG, "node": "Node1"}
    ).content.decode()
    section = page.split('id="analogs"', 1)[1].split("</section>", 1)[0]
    rows = section.split("<tbody>", 1)[1].split("</tbody>", 1)[0].split("<tr")[1:]
    assert len(rows) == 3 and "<strong>Node1</strong>" in rows[0]
    assert "Node6" in rows[1] and "100.0%" in rows[1]
    # MUSCLE's Node9 differs at one residue of three, placed through the extant residues.
    assert "Node9" in rows[2] and "66.7%" in rows[2] and 'class="pop-flag"' in rows[2]
    assert ">Add all 3 to a new library</button>" in section
    assert 'form="library-add-analogs"' in section and 'id="library-add-analogs"' in page
    data = json.loads(
        section.split('id="aln-data-c0" type="application/json">', 1)[1].split("<", 1)[0]
    )
    # One tree source, so rows name the model alone; ◂ marks the ancestor shown.
    assert data["names"][-3:] == ["MAFFT Node1 · LG ◂", "MAFFT Node6 · WAG", "MUSCLE Node9 · LG"]
    # Another layout comes from the menu; anything else falls back to this alignment.
    each = client.get(
        f"/jobs/{job.id}/ancestors/", {"combination": LG, "node": "Node1", "compare_on": "own"}
    ).content.decode()
    assert 'id="aln-data-c1"' in each
    junk = client.get(
        f"/jobs/{job.id}/ancestors/", {"combination": LG, "node": "Node1", "compare_on": "x"}
    ).content.decode()
    assert 'id="aln-data-c1"' not in junk and "On MAFFT columns" in junk


@pytest.mark.django_db
def test_switching_combination_keeps_the_ancestor(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    page = client.get(
        f"/jobs/{job.id}/ancestors/",
        {"aligner": "muscle", "combination": MUSCLE, "node": "Node1", "shown": LG},
    ).content.decode()
    assert '<option value="Node9" selected>' in page


@pytest.mark.django_db
def test_add_all_puts_every_analog_in_the_library(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    url = f"/jobs/{job.id}/library/add-analogs/"
    data = {"combination": LG, "gaps": "original", "node": "Node1"}
    response = client.post(url, data | {"next": f"/jobs/{job.id}/ancestors/"})
    assert response.status_code == 302
    targets = Library.objects.get().targets
    assert [(t["combination"], t["node"]) for t in targets] == [
        (LG, "Node1"), (WAG, "Node6"), (MUSCLE, "Node9"),
    ]  # fmt: skip
    assert all(t["version"] == "original" for t in targets)
    client.post(url, data)  # again: nothing new
    assert len(Library.objects.get().targets) == 3
    # Only the job's own values; the analogs are found on the server, never posted.
    assert client.post(url, data | {"node": "Node99"}).status_code == 404
    assert client.post(url, data | {"combination": "x/y/z"}).status_code == 404
    assert client.post(url, data | {"gaps": "fitch"}).status_code == 404
