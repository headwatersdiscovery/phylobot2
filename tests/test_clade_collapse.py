"""Folding clades into triangles in the tree viewer (phylobot.clade_collapse): folded once,
the triangle, its label and its style show on every combination that has the clade, and the
fold persists. User text and colours are cleaned and never reach a pipeline file."""

from typing import Any

import pytest
from django.contrib.auth.models import User
from django.test import Client

from phylobot import clade_collapse
from phylobot.models import CladeCollapse, Job
from tests.conftest import FakeStorage
from tests.test_extend import split
from tests.test_job_import import download, edit_manifest, refused, upload
from tests.test_views import viewer_job

pytestmark = pytest.mark.django_db

WAG, MUSCLE = "mafft/iqtree/WAG", "muscle/iqtree/modelfinder"
# viewer_job's trees: MAFFT's ((a,b)Node2,c)Node1, MUSCLE's the same numbered the other way
# round, so MAFFT's Node2 clade is MUSCLE's Node1.


def collapse(
    client: Client, job: Job, node: str, combination: str = WAG,
    action: str = "collapse", **more: str,
) -> Any:  # fmt: skip
    data = {"combination": combination, "aligner": combination.split("/")[0], "gaps": "original",
            "node": node, "action": action, **more}  # fmt: skip
    return client.post(f"/jobs/{job.id}/tree/collapse/", data)


def tree(client: Client, job: Job, combination: str = MUSCLE, **extra: str) -> str:
    query = {"combination": combination, "aligner": combination.split("/")[0], **extra}
    return client.get(f"/jobs/{job.id}/tree/", query).content.decode()


# ── the value cleaners ──────────────────────────────────────────────────────────────────


def test_labels_are_kept_printable_and_bounded() -> None:
    assert clade_collapse.clean_label("  Fungal IME2s ") == "Fungal IME2s"  # spaces allowed
    assert clade_collapse.clean_label("a\x00b<c>") == "abc"  # control + angle brackets out
    assert clade_collapse.clean_label("x" * 80) == "x" * clade_collapse.LABEL_MAX


def test_colours_are_validated() -> None:
    fill = clade_collapse.DEFAULT_FILL
    assert clade_collapse.clean_colour("#AbC123", "#000000") == "#abc123"
    assert clade_collapse.clean_colour("red", fill) == fill  # not a hex: the fallback
    assert clade_collapse.clean_colour("none", "#000000", allow_none=True) == "none"
    assert clade_collapse.clean_colour("none", "#000000") == "#000000"  # not allowed here


def test_border_width_snaps_to_an_offered_value() -> None:
    assert clade_collapse.clean_width("2.4") == 2.5
    assert clade_collapse.clean_width("nonsense") == clade_collapse.DEFAULT_BORDER_WIDTH


# ── the web flow ────────────────────────────────────────────────────────────────────────


def test_folded_once_shown_on_every_combination(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    response = collapse(client, job, "Node2", label="AB clade", fill="#dfe7d6")
    assert response.status_code == 302
    row = CladeCollapse.objects.get()
    assert row.job_id == job.id and row.label == "AB clade" and row.fill == "#dfe7d6"
    # MUSCLE numbers the same clade Node1: its triangle is drawn there too.
    page = tree(client, job)
    assert 'class="tree-clade" data-node="Node1"' in page
    assert "<polygon" in page and 'fill="#dfe7d6"' in page and ">AB clade</tspan>" in page
    assert 'clade-popover" data-node="Node1"' in page  # the expand/restyle popover
    # The member tip "a" is hidden, so its leaf label is gone from the drawing.
    assert ">a</text>" not in page


def test_a_named_ancestor_shows_on_its_triangle(client: Client, fake_storage: FakeStorage) -> None:
    """A folded clade's ancestor is not drawn as a node, so its user name would vanish: the
    triangle carries it at an apex dot (where the ancestor sits)."""
    job = viewer_job(fake_storage)
    client.post(
        f"/jobs/{job.id}/ancestors/name/",
        {
            "combination": WAG,
            "aligner": "mafft",
            "gaps": "original",
            "node": "Node2",
            "ancestor_name": "AncX",
        },  # fmt: skip
    )
    collapse(client, job, "Node2")
    page = tree(client, job)  # MUSCLE numbers the same clade Node1
    group = page.split('class="tree-clade" data-node="Node1"', 1)[1].split("</g>", 1)[0]
    assert '<tspan font-weight="bold">AncX</tspan>' in group  # the name, at the apex
    assert 'class="dot"' in group  # the apex dot marking the ancestor
    # Hiding node labels hides the name but keeps the dot.
    hidden = tree(client, job, node_labels="none")
    hidden_group = hidden.split('class="tree-clade" data-node="Node1"', 1)[1].split("</g>", 1)[0]
    assert '<tspan font-weight="bold">AncX</tspan>' not in hidden_group  # label hidden
    assert 'class="dot"' in hidden_group  # dot kept; the title tooltip still names it
    assert "<title>AncX —" in hidden_group


def test_restyle_then_expand(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    collapse(client, job, "Node2")
    collapse(client, job, "Node2", action="update", label="renamed", border_width="4.0")
    row = CladeCollapse.objects.get()
    assert row.label == "renamed" and row.border_width == 4.0
    collapse(client, job, "Node1", combination=MUSCLE, action="expand")  # expand elsewhere
    assert not CladeCollapse.objects.exists()
    assert 'class="tree-clade"' not in tree(client, job)


def test_the_node_popover_offers_collapse(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    page = client.get(
        f"/jobs/{job.id}/tree/", {"combination": WAG, "aligner": "mafft"}
    ).content.decode()
    popover = page.split('class="node-popover" data-node="Node2"', 1)[1].split("</div>", 1)[0]
    assert "data-pop-collapse" in popover and "Collapse clade" in popover
    assert 'id="tree-clade-collapse"' in page  # the POST form it goes through


def test_only_the_owner_folds_real_nodes(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    assert collapse(client, job, "Node99").status_code == 404
    Job.objects.filter(id=job.id).update(owner=User.objects.create(username="other"))
    assert collapse(client, job, "Node2").status_code == 404
    assert client.get(f"/jobs/{job.id}/tree/collapse/").status_code in (404, 405)
    assert not CladeCollapse.objects.exists()


def test_an_extension_shares_its_parents_folds(client: Client, fake_storage: FakeStorage) -> None:
    parent, child = split(fake_storage)
    collapse(client, child, "Node2")  # folded from the extension's tree viewer
    assert CladeCollapse.objects.get().job_id == parent.id  # held by the parent
    assert 'class="tree-clade"' in tree(client, child, MUSCLE)


@pytest.mark.parametrize(
    "where",
    ["https://evil.example/", "/jobs/", "/jobs/{other}/tree/", "/jobs/{job}/ancestors/?x=1"],
)  # fmt: skip
def test_it_returns_only_to_this_jobs_tree(
    client: Client, fake_storage: FakeStorage, where: str
) -> None:
    job = viewer_job(fake_storage)
    asked = where.format(job=job.id, other="0" * 8 + "-0000-0000-0000-" + "0" * 12)
    response = collapse(client, job, "Node2", return_to=asked)
    assert response["Location"] == f"/jobs/{job.id}/tree/"


def test_folds_travel_in_the_zip(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    collapse(client, job, "Node2", label="AB clade", fill="#dfe7d6", border_width="2.5")
    data = download(client, job)
    Job.objects.filter(id=job.id).update(name="original")  # no name clash on import
    assert upload(client, data).status_code == 302
    imported = Job.objects.get(imported_at__isnull=False)
    row = imported.clade_collapses.get()
    assert row.label == "AB clade" and row.fill == "#dfe7d6" and row.border_width == 2.5


def test_bad_folds_in_a_zip_are_refused(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    bad = [{"clade": "zz", "parts": "0" * 16, "label": "x"}]  # clade not a well-formed ID
    data = edit_manifest(download(client, job), lambda m: m["job"].update(clade_collapses=bad))
    Job.objects.filter(id=job.id).update(name="original")
    assert "collapsed clades" in refused(client, data)
    assert not CladeCollapse.objects.exists()
