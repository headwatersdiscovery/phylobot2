import json
import re
from typing import Any

import pytest
from django.core.cache import cache
from django.test import Client

from phylobot import jobs
from phylobot.models import Job
from phylobot.modes import local_user
from tests.conftest import FakeStorage
from tests.test_views import finished, other_users_job, viewer_job

pytestmark = pytest.mark.django_db


@pytest.fixture(autouse=True)
def empty_cache() -> None:
    cache.clear()  # composition tracks are cached per job


def data(page: str) -> dict[str, Any]:
    match = re.search(r'<script id="aln-data" type="application/json">(.*?)</script>', page, re.S)
    assert match, "no viewer data on the page"
    return dict(json.loads(match.group(1)))


def test_defaults_and_links(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    page = client.get(f"/jobs/{job.id}/alignment/").content.decode()
    d = data(page)
    assert d["names"] == ["a", "b", "c<i>"] and d["rows"] == ["MK-", "M-K", "MKK"]
    assert not any(d["ancestor"]) and d["scheme"] == "clustalx" and not d["codon"]
    assert [t["key"] for t in d["tracks"]] == ["jsd"]  # owner decision, 2026-11-11
    assert "c&lt;i&gt;" not in page or "c<i>" not in page.split('id="aln-data"')[0]
    # Tree and model choose the ancestors and the tree order: hidden until either is used
    assert 'name="combination" value="mafft/iqtree/modelfinder"' in page  # kept, hidden
    assert "Show ancestors" in page and '<select name="combination">' not in page
    shown = client.get(f"/jobs/{job.id}/alignment/?ancestors=on").content.decode()
    assert '<option value="mafft/iqtree/modelfinder" selected>ModelFinder → LG+G4' in shown
    assert "3 sequences × 3 columns" in page
    assert d["ancestorViewer"].startswith(f"/jobs/{job.id}/ancestors/?aligner=mafft")
    job_page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Alignment Viewer →" in job_page
    assert f"/jobs/{job.id}/alignment/?aligner=muscle" in job_page
    tree = client.get(f"/jobs/{job.id}/tree/").content.decode()
    assert "alignment in tree order" in tree
    ancestors = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
    assert "show in the alignment" in ancestors


def test_tree_order_with_ancestors(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    url = f"/jobs/{job.id}/alignment/?order=tree&ancestors=on"
    page = client.get(url).content.decode()
    d = data(page)
    assert d["names"] == ["Node1", "c<i>", "Node2", "a", "b"]
    assert d["ancestor"] == [True, False, True, False, False]
    assert d["rows"][0] == "MK-"  # the job's gap method (Fitch) by default
    assert '<option value="fitch" selected>' in page
    d = data(client.get(url + "&gaps=original&aligner=muscle").content.decode())
    assert d["rows"][0] == "MRK" and d["names"][0] == "Node2"  # muscle's root is Node2
    # alignment order: ancestors below the extant sequences, root first
    d = data(client.get(f"/jobs/{job.id}/alignment/?ancestors=on").content.decode())
    assert d["names"] == ["a", "b", "c<i>", "Node1", "Node2"]


def test_tracks_are_edited_in_the_url(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    base = f"/jobs/{job.id}/alignment/"
    response = client.get(base + "?add_track=changes")
    assert response.status_code == 302
    assert "track=jsd&track=changes" in response["Location"]
    page = client.get(response["Location"]).content.decode()
    changes = data(page)["tracks"][1]
    assert changes["key"] == "changes"
    # column 2: Node1 K -> Node2 R -> a K; Node2 -> b is a gap (not a substitution)
    assert changes["series"][0]["values"] == [0, 2, 0]
    location = client.get(base + "?tracks_set=1&track=jsd&track=entropy&track_up=1")["Location"]
    assert "track=entropy&track=jsd" in location
    location = client.get(base + "?tracks_set=1&track=jsd&remove_track=0")["Location"]
    assert "track=" not in location and "tracks_set=1" in location
    assert data(client.get(location).content.decode())["tracks"] == []
    # unknown or unavailable tracks are dropped
    d = data(client.get(base + "?tracks_set=1&track=evil&track=rate&track=helix").content.decode())
    assert [t["key"] for t in d["tracks"]] == ["helix"]
    # HTMX gets the fragment, and the address bar the canonical view
    response = client.get(base + "?add_track=helix", headers={"HX-Request": "true"})
    assert response.status_code == 200 and "track=helix" in response["HX-Push-Url"]
    assert "<h1>" not in response.content.decode()


def test_fallbacks_for_unknown_values(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    page = client.get(
        f"/jobs/{job.id}/alignment/?aligner=evil&scheme=evil&zoom=evil&reference=evil"
        "&combination=evil&order=evil&gaps=evil&states=codons"
    ).content.decode()
    d = data(page)
    assert d["scheme"] == "clustalx" and d["reference"] == -1 and not d["codon"]
    assert d["cell"] == 11.0  # the default zoom
    assert d["names"] == ["a", "b", "c<i>"]
    d = data(
        client.get(f"/jobs/{job.id}/alignment/?scheme=differences&reference=b").content.decode()
    )
    assert d["reference"] == 1


def test_zoom(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)

    def cell(zoom: str) -> object:
        page = client.get(f"/jobs/{job.id}/alignment/", {"zoom": zoom}).content.decode()
        return data(page)["cell"]

    assert cell("fit") == "fit" and cell("3.5") == 3.5
    assert cell("999") == 20.0 and cell("0") == 0.05  # the slider's ends
    assert cell("small") == 6.0  # links saved from the old menu
    for junk in ("nan", "evil", ""):
        assert cell(junk) == 11.0
    page = client.get(f"/jobs/{job.id}/alignment/", {"zoom": "3.5"}).content.decode()
    assert 'name="zoom" value="3.5"' in page and "zoom-slider" in page


def test_exports(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    base = f"/jobs/{job.id}/alignment/?order=tree&ancestors=on&gaps=original"
    fasta = client.get(base + "&export=fasta")
    assert fasta["Content-Disposition"] == 'attachment; filename="alignment.mafft.fasta"'
    assert fasta.content.decode().startswith(">Node1\nMKK\n>c<i>\nMKK\n>Node2\n")
    clustal = client.get(base + "&export=clustal").content.decode()
    assert clustal.startswith("CLUSTAL W") and "Node2      MRK" in clustal
    tracks = client.get(base + "&export=tracks").content.decode().splitlines()
    header = tracks[0].split("\t")
    assert header[0] == "column" and "Conservation (Jensen-Shannon)" in header
    assert "Substitutions on the tree" in header and "Charge" in header  # every track
    assert len(tracks) == 4


def test_private_and_unfinished(client: Client, fake_storage: FakeStorage) -> None:
    assert client.get(f"/jobs/{other_users_job().id}/alignment/").status_code == 404
    job = Job.objects.create(owner=local_user())
    assert jobs.mark(job.id, Job.Status.RUNNING)
    page = client.get(f"/jobs/{job.id}/alignment/").content.decode()
    assert "the job has not finished successfully" in page and "aln-data" not in page


def test_codon_alignment(client: Client, fake_storage: FakeStorage) -> None:
    job = finished({"seqtype": "codon", "genetic_code": 1})
    prefix = f"jobs/{job.id}/results/"
    mafft = {"key": "mafft", "label": "MAFFT", "path": "mafft/alignment.fasta"}
    summary = {"stop_after": "alignments", "combinations": [], "alignments": [mafft]}
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    fake_storage.objects[prefix + "mafft/alignment.fasta"] = b">a\nMK\n>b\nM-\n"
    fake_storage.objects[prefix + "mafft/alignment.cds.fasta"] = b">a\nATGAAA\n>b\nATG---\n"
    page = client.get(f"/jobs/{job.id}/alignment/").content.decode()
    assert not data(page)["codon"] and "Synonymous differences" not in page
    page = client.get(f"/jobs/{job.id}/alignment/?states=codons&scheme=synonymous").content.decode()
    d = data(page)
    assert d["codon"] and d["rows"] == ["ATGAAA", "ATG---"] and d["amino"] == ["MK", "M-"]
    assert d["scheme"] == "synonymous"
    assert "2 sequences × 2 columns (codons)" in page
    job_page = client.get(f"/jobs/{job.id}/").content.decode()
    assert f"/jobs/{job.id}/alignment/?aligner=mafft&amp;states=codons" in job_page


def test_branch_site_job_shows_its_codon_alignment(
    client: Client, fake_storage: FakeStorage
) -> None:
    parent = finished({"seqtype": "codon"})
    job = finished({"branch_site": True, "genetic_code": 1, "branch_site_of": {
        "job": str(parent.id), "aligner": "mafft", "node": "Node3", "scope": "clade"}})  # fmt: skip
    prefix = f"jobs/{job.id}/"
    fake_storage.objects[prefix + "results/summary.json"] = json.dumps(
        {"branch_site": True, "stop_after": "branch_site", "alignments": []}
    ).encode()
    fake_storage.objects[prefix + "inputs/alignment1.fasta"] = b">a\nATGAAA\n>b\nATGAAG\n"
    fake_storage.objects[prefix + "results/branch_site/sites.tsv"] = (
        b"site\tresidue\tp0\tp1\tp2a\tp2b\tp_positive\n1\tM\t1\t0\t0\t0\t0.0\n2\tK\t0\t0\t.5\t.4\t0.9\n"
    )
    page = client.get(f"/jobs/{job.id}/alignment/?tracks_set=1&track=branch_site").content.decode()
    d = data(page)
    assert d["codon"] and d["amino"] == ["MK", "MK"]
    assert d["tracks"][0]["series"] == [{"label": "Node3 (clade)", "values": [0.0, 0.9]}]
    assert d["ancestorViewer"] == "" and "click a column" not in page
    # the parent job's viewer layers its tests (here: one, for this aligner)
    parent_prefix = f"jobs/{parent.id}/results/"
    fake_storage.objects[parent_prefix + "summary.json"] = json.dumps({"alignments": [
        {"key": "mafft", "label": "MAFFT", "path": "mafft/alignment.fasta"}]}).encode()  # fmt: skip
    fake_storage.objects[parent_prefix + "mafft/alignment.fasta"] = b">a\nMK\n>b\nMK\n"
    fake_storage.objects[parent_prefix + "mafft/alignment.cds.fasta"] = b">a\nATGAAA\n>b\nATGAAG\n"
    d = data(client.get(f"/jobs/{parent.id}/alignment/?add_track=branch_site",
                        follow=True).content.decode())  # fmt: skip
    assert d["tracks"][-1]["series"][0]["values"] == [0.0, 0.9]


def test_site_rate_and_likelihood_tracks(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    folder = f"jobs/{job.id}/results/mafft/iqtree/modelfinder/"
    base = f"/jobs/{job.id}/alignment/?tracks_set=1&track=rate&track=site_lh&track=gap_calls"
    assert data(client.get(base).content.decode())["tracks"][0]["key"] == "gap_calls"
    fake_storage.objects[folder + "site_rates.tsv"] = (
        b"site\trate\tcategory\tcategory_rate\n1\t0.2\t1\t0.1\n2\t1.5\t3\t1.4\n3\t0.9\t2\t0.8\n"
    )
    fake_storage.objects[folder + "site_lh.tsv"] = b"site\tlog_likelihood\n1\t-1.5\n2\t-7\n3\t-3\n"
    tracks = data(client.get(base).content.decode())["tracks"]
    assert [t["key"] for t in tracks] == ["rate", "site_lh", "gap_calls"]
    assert tracks[0]["series"][0]["values"] == [0.2, 1.5, 0.9]
    assert tracks[1]["series"][0]["values"] == [-1.5, -7, -3]
    # one series per gap method, layered in one plot
    assert [s["label"] for s in tracks[2]["series"]] == ["Fitch", "ML"]
    assert tracks[2]["series"][0]["values"] == [0, 0.5, 0.5]  # Fitch: Node1 MK-, Node2 M-K


def test_viewer_height(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    page = client.get(f"/jobs/{job.id}/alignment/").content.decode()
    assert data(page)["height"] == "fit"  # the window, scrollbar on screen
    assert '<option value="fit" selected>Fit the window' in page
    assert (
        data(client.get(f"/jobs/{job.id}/alignment/?height=20").content.decode())["height"] == "20"
    )
    assert (
        data(client.get(f"/jobs/{job.id}/alignment/?height=9999").content.decode())["height"]
        == "fit"
    )


def test_the_tracks_chosen_are_remembered_per_alignment(
    client: Client, fake_storage: FakeStorage
) -> None:
    """Adding, removing or reordering a track is remembered for that alignment: the viewer
    opens it with those tracks next time, and "none" is remembered as well as any."""
    job = viewer_job(fake_storage)
    base = f"/jobs/{job.id}/alignment/"
    shown = lambda url: [t["key"] for t in data(client.get(url).content.decode())["tracks"]]  # noqa: E731
    assert shown(base) == ["jsd"]
    client.get(base + "?add_track=entropy")
    assert shown(base) == ["jsd", "entropy"]  # a plain visit: what was last chosen
    client.get(base + "?tracks_set=1&track=jsd&track=entropy&track_up=1")
    assert shown(base) == ["entropy", "jsd"]
    client.get(base + "?tracks_set=1&track=entropy&track=jsd&remove_track=0")
    client.get(base + "?tracks_set=1&track=jsd&remove_track=0")
    assert shown(base) == []  # every track removed: still none on the next visit
    # A URL that names its tracks shows those, and looking is not choosing.
    assert shown(base + "?tracks_set=1&track=jsd") == ["jsd"]
    assert shown(base) == []
