"""My alignments and the alignment editor (phylobot.alignment_editor_views), and their
use on the New Job page."""

import json
import re
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.core.cache import cache
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import Client

from phylobot import alignments, jobs
from phylobot.models import Alignment, Job
from phylobot.modes import local_user
from tests.conftest import FakeStorage
from tests.test_views import ALIGNED_FOUR, finished, post

pytestmark = pytest.mark.django_db

UPLOAD = ">seq2 second\nMKA-L\n>seq10\nMRA-L\n>Alpha\nM-G-V\n"


@pytest.fixture(autouse=True)
def empty_cache() -> None:
    cache.clear()


def data(page: str) -> dict[str, Any]:
    match = re.search(r'<script id="aln-data" type="application/json">(.*?)</script>', page, re.S)
    assert match, "no viewer data on the page"
    return dict(json.loads(match.group(1)))


def uploaded(client: Client, text: str = UPLOAD, **fields: Any) -> Alignment:
    response = client.post(
        "/alignments/upload/", {"file": SimpleUploadedFile("my aln.fasta", text.encode()), **fields}
    )
    assert response.status_code == 302, response.content.decode()[-3000:]
    record = Alignment.objects.get(id=response["Location"].split("/")[2])
    return record


def job_with_alignments(fake_storage: FakeStorage, name: str = "Ime2 run") -> Job:
    job = finished({"stop_after": "alignments", "aligners": ["mafft"]})
    Job.objects.filter(id=job.id).update(name=name)  # not save(): `job` still says queued
    job.refresh_from_db()
    prefix = f"jobs/{job.id}/results/"
    summary = {"stop_after": "alignments", "combinations": [], "alignments": [
        {"key": "mafft", "label": "MAFFT", "path": "mafft/alignment.fasta"},
        {"key": "provided1", "label": "Curated", "path": "provided1/alignment.fasta"},
    ]}  # fmt: skip
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    for key in ("mafft", "provided1"):
        fake_storage.objects[f"{prefix}{key}/alignment.fasta"] = ALIGNED_FOUR.encode()
    return job


def test_upload_opens_the_editor(client: Client, fake_storage: FakeStorage) -> None:
    record = uploaded(client)
    assert (
        record.name == "my aln" and record.origin == "upload" and record.filename == "my aln.fasta"
    )
    assert (record.sequences, record.columns, record.seqtype) == (3, 5, "protein")
    assert fake_storage.objects[record.storage_key].decode() == UPLOAD
    page = client.get(f"/alignments/{record.id}/").content.decode()
    d = data(page)
    assert d["names"] == ["seq2", "seq10", "Alpha"] and d["rows"][2] == "M-G-V"
    assert d["editable"] and d["sort"] == "" and not d["codon"]
    assert "Remove 1 all-gap column" in page
    listing = client.get("/alignments/").content.decode()
    assert "my aln" in listing and f"/alignments/{record.id}/download/" in listing


def test_bad_uploads_are_refused(client: Client, fake_storage: FakeStorage) -> None:
    response = client.post("/alignments/upload/", {"text": ">a\nMKL\n>b\nMK\n"})
    assert response.status_code == 400
    assert "different lengths" in response.content.decode()
    response = client.post("/alignments/upload/", {"text": ""})
    assert "Choose a file or paste an alignment" in response.content.decode()
    assert not Alignment.objects.exists() and not fake_storage.objects


def test_a_name_used_twice_is_resolved_on_upload(client: Client, fake_storage: FakeStorage) -> None:
    # Two different rows named b, and a repeat of a (the same sequence, used once).
    text = ">a\nMK-LV\n>b\nMKALV\n>c\nMKALI\n>d\nMR-LV\n>b\nMKSLV\n>a\nMK-LV\n"

    def check(**fields: Any) -> str:
        upload = SimpleUploadedFile("dup.fasta", text.encode())
        response = client.post("/alignments/upload/check/", {"file": upload, **fields})
        return re.sub(r"\s+", " ", response.content.decode())

    page = check()
    assert "1 sequence name is used for different rows" in page and "Not saved" not in page
    assert "<code>b.2</code>" in page and "Cancel this upload" in page
    assert "1 row repeats a sequence already in the alignment" in page
    assert "✓ 5 protein sequences × 5 columns: ready to save" in check(duplicates="rename")
    response = client.post(
        "/alignments/upload/", {"file": SimpleUploadedFile("dup.fasta", text.encode())}
    )
    assert response.status_code == 400 and not Alignment.objects.exists()
    renamed = uploaded(client, text, duplicates="rename")

    def names(record: Alignment) -> list[str]:
        stored = fake_storage.objects[record.storage_key].decode()
        return [r.name for r in alignments.read_rows(stored)]

    assert names(renamed) == ["a", "b", "c", "d", "b.2"]
    first = uploaded(client, text, duplicates="first")
    assert names(first) == ["a", "b", "c", "d"]
    nothing = client.post("/alignments/upload/check/", {"text": ""})
    assert nothing.status_code == 200 and not nothing.content.decode().strip()


PHYLIP_AS_FASTA = ">seq2\nMKA-LV\n>seq10\nMRA-LV\n>Alpha\nM-G-VV\n"


@pytest.mark.parametrize(
    "text",
    [
        "3 6\nseq2 MKA-LV\nseq10 MRA-LV\nAlpha M-G-VV\n",  # relaxed, one line each
        " 3  6 \n\nseq2\tMKA -LV\nseq10  MRA-LV\nAlpha M-G-VV",  # tabs, spaces, blanks
        "3 6\nseq2 MKA\n-LV\nseq10 MR\nA-LV\nAlpha M-G-VV\n",  # sequential, wrapped
        "3 6\nseq2 MKA\nseq10 MRA\nAlpha M-G\n\n-LV\n-LV\n-VV\n",  # interleaved
        "3 6\nseq2      MKA-LV\nseq10     MRA-LV\nAlpha     M-G-VV\n",  # strict, padded
        "3 6\nseq2______MKA-LV\nseq10_____MRA-LV\nAlpha_____M-G-VV\n".replace("_", " "),
    ],
)
def test_phylip_layouts(text: str) -> None:
    assert alignments.phylip_to_fasta(text) == PHYLIP_AS_FASTA


def test_phylip_strict_names_without_a_space() -> None:
    text = "2 4\nHomo sapieMKLV\nMus muscluMK-V\n"
    assert alignments.phylip_to_fasta(text) == ">Homo_sapie\nMKLV\n>Mus_musclu\nMK-V\n"


def test_phylip_refusals() -> None:
    assert alignments.phylip_to_fasta(">a\nMK\n>b\nMK\n") is None
    assert alignments.phylip_to_fasta("") is None
    with pytest.raises(alignments.AlignmentError, match="do not fit"):
        alignments.phylip_to_fasta("2 5\na MKLV\nb MKLV\n")
    with pytest.raises(alignments.AlignmentError, match="do not fit"):
        alignments.phylip_to_fasta("3 4\na MKLV\nb MKLV\n")
    with pytest.raises(alignments.AlignmentError, match="header says 0"):
        alignments.phylip_to_fasta("0 4\n")


def test_phylip_paste_and_upload(client: Client, fake_storage: FakeStorage) -> None:
    text = "3 6\nseq2 MKA\nseq10 MRA\nAlpha M-G\n\n-LV\n-LV\n-VV\n"
    check = client.post("/alignments/upload/check/", {"text": text}).content.decode()
    assert "✓ 3 protein sequences × 6 columns: ready to save" in check
    response = client.post("/alignments/upload/", {"text": text})
    assert response.status_code == 302
    record = Alignment.objects.get(id=response["Location"].split("/")[2])
    assert record.name == "Pasted alignment" and (record.sequences, record.columns) == (3, 6)
    assert fake_storage.objects[record.storage_key].decode() == PHYLIP_AS_FASTA
    record = uploaded(client, text)
    assert fake_storage.objects[record.storage_key].decode() == PHYLIP_AS_FASTA
    # "?" (missing data) is stored as a gap, in PHYLIP as in FASTA.
    record = uploaded(client, "3 6\nseq2 MKA?LV\nseq10 MRA-LV\nAlpha M?G.VV\n")
    assert fake_storage.objects[record.storage_key].decode() == PHYLIP_AS_FASTA
    bad = client.post("/alignments/upload/", {"text": "3 7\nseq2 MKA-LV\n"})
    assert bad.status_code == 400 and "PHYLIP header says 3 sequences" in bad.content.decode()


def test_codon_upload(client: Client, fake_storage: FakeStorage) -> None:
    record = uploaded(client, ">a\nATG---AAA\n>b\nATGCCCAAA\n", genetic_code="1")
    assert record.seqtype == "codon" and record.columns == 3
    d = data(client.get(f"/alignments/{record.id}/").content.decode())
    assert d["codon"] and d["amino"] == ["M-K", "MPK"]
    response = client.get(f"/alignments/{record.id}/download/?states=amino")
    assert response.content.decode() == ">a\nM-K\n>b\nMPK\n"
    # Cutting codon column 2 cuts three characters.
    response = client.get(f"/alignments/{record.id}/?cut_columns=2")
    d = data(client.get(response["Location"]).content.decode())
    assert d["rows"] == ["ATGAAA", "ATGAAA"]


def test_edits_go_through_the_url(client: Client, fake_storage: FakeStorage) -> None:
    record = uploaded(client)
    base = f"/alignments/{record.id}/"
    response = client.get(base, {"sort_toggle": "1"})
    assert response.status_code == 302 and "sort=name" in response["Location"]
    response = client.get(response["Location"] + "&cut_columns=1-2")  # shown: Alpha, seq2, seq10
    location = response["Location"]
    assert "cut=1-2" in location
    response = client.get(location + "&drop_rows=1", HTTP_HX_REQUEST="true")
    assert response.status_code == 200
    pushed = response["HX-Push-Url"]
    assert "drop=3" in pushed and "cut=1-2" in pushed and "sort=name" in pushed
    d = data(response.content.decode())
    assert d["names"] == ["seq2", "seq10"] and d["rows"] == ["A-L", "A-L"]
    assert d["origin"] == [3, 4, 5]
    page = response.content.decode()
    assert "Remove 1 all-gap column" in page and "Start over" in page
    assert "Back to file order" in page and ">Sort by name" not in page  # the header sorts
    assert page.index('class="aln-save"') < page.index('id="alignment-viewer"')  # Save on top
    exported = client.get(pushed + "&export=fasta")
    assert exported.content.decode() == ">seq2 second\nA-L\n>seq10\nA-L\n"
    assert 'filename="my_aln.edited.fasta"' in exported["Content-Disposition"]
    # The stored alignment never changes.
    assert fake_storage.objects[record.storage_key].decode() == UPLOAD


def test_save_makes_a_new_alignment(client: Client, fake_storage: FakeStorage) -> None:
    record = uploaded(client)
    posted = {"cut": "4", "drop": "3", "sort": "name", "name": "  Trimmed   v1 "}
    response = client.post(f"/alignments/{record.id}/save/", posted)
    new = Alignment.objects.exclude(id=record.id).get()
    assert response["Location"] == f"/alignments/{new.id}/"
    assert new.name == "Trimmed v1" and new.origin == "edit" and new.derived_from == record
    assert new.edits == {"parent_sites": 5, "kept_sites": "1-3,5", "removed_sequences": ["Alpha"],
                         "sort": "name", "parent_version": 1}  # fmt: skip
    assert fake_storage.objects[new.storage_key].decode() == ">seq2 second\nMKAL\n>seq10\nMRAL\n"
    assert (new.sequences, new.columns) == (2, 4)
    page = client.get(f"/alignments/{new.id}/").content.decode()
    assert "kept columns 1-3,5 of 5" in page and f"/alignments/{record.id}/" in page


def test_rename_and_delete(client: Client, fake_storage: FakeStorage) -> None:
    record = uploaded(client)
    client.post(f"/alignments/{record.id}/rename/", {"name": "Kinase domain"})
    record.refresh_from_db()
    assert record.name == "Kinase domain"
    client.post(f"/alignments/{record.id}/rename/", {"name": "   "})  # blank: unchanged
    record.refresh_from_db()
    assert record.name == "Kinase domain"
    # The editor's title pencil posts with the page (and its edits) to come back to.
    page = client.get(f"/alignments/{record.id}/", {"cut": "2"}).content.decode()
    assert "data-rename" in page and 'aria-label="Alignment name"' in page
    back = f"/alignments/{record.id}/?cut=2"
    response = client.post(f"/alignments/{record.id}/rename/", {"name": "Kinase", "next": back})
    assert response["Location"] == back
    record.refresh_from_db()
    assert record.name == "Kinase"
    response = client.post(f"/alignments/{record.id}/delete/")
    assert response.status_code == 302
    record.refresh_from_db()
    assert record.deleted_at is not None and not fake_storage.objects  # the record stays
    assert client.get(f"/alignments/{record.id}/").status_code == 404


def test_others_alignments_are_not_found(client: Client, fake_storage: FakeStorage) -> None:
    other = User.objects.create(username="other")
    record = alignments.store(other, "theirs", alignments.read_rows(UPLOAD), "protein", 1,
                              origin=Alignment.Origin.UPLOAD)  # fmt: skip
    for url in ("", "download/"):
        assert client.get(f"/alignments/{record.id}/{url}").status_code == 404
    for url in ("rename/", "delete/", "save/", "duplicate/"):
        assert client.post(f"/alignments/{record.id}/{url}", {"name": "x"}).status_code == 404
    assert "theirs" not in client.get("/alignments/").content.decode()


def test_jobs_register_their_alignments(client: Client, fake_storage: FakeStorage) -> None:
    job = job_with_alignments(fake_storage)
    page = client.get("/alignments/").content.decode()
    assert not Alignment.objects.filter(job=job).exists()  # the page itself only links the section
    assert 'hx-get="/alignments/jobs/"' in page
    page = client.get("/alignments/jobs/").content.decode()  # older jobs: registered on first visit
    registered = list(Alignment.objects.filter(job=job).order_by("job_key"))
    assert [(a.job_key, a.name, a.origin, a.path) for a in registered] == [
        ("mafft", "MAFFT", "job", "mafft/alignment.fasta"),
        ("provided1", "Curated", "job", "provided1/alignment.fasta"),
    ]
    assert alignments.register_job(job) == registered  # idempotent
    assert "Ime2 run" in page and "Edit a copy" in page
    mafft = registered[0]
    assert mafft.title == "Ime2 run · MAFFT"
    assert f"/alignments/{mafft.id}/rename/" not in page  # read-only
    assert client.post(f"/alignments/{mafft.id}/rename/", {"name": "x"}).status_code == 404
    assert client.post(f"/alignments/{mafft.id}/delete/").status_code == 404
    editor = client.get(f"/alignments/{mafft.id}/").content.decode()
    assert "Made by the job" in editor and data(editor)["names"][0] == "seq0"
    # Saving an edit of a job's alignment gives the user their own.
    client.post(f"/alignments/{mafft.id}/save/", {"cut": "10", "name": ""})
    copy = Alignment.objects.get(origin="edit")
    assert copy.name == "Ime2 run · MAFFT (edited)" and copy.job is None and copy.columns == 10
    # The job's viewer links to the editor.
    viewer = client.get(f"/jobs/{job.id}/alignment/?scheme=chemistry&zoom=5").content.decode()
    head = viewer.split('id="aln-viewer-actions"')[1].split("</div>")[0]
    assert f'href="/alignments/{mafft.id}/?scheme=chemistry&amp;zoom=5.0' in head
    assert "Edit a copy" in head and "track=jsd" in head
    # It follows the alignment chosen (an out-of-band swap with the viewer).
    swap = client.get(f"/jobs/{job.id}/alignment/?aligner=provided1", HTTP_HX_REQUEST="true")
    curated = Alignment.objects.get(job=job, job_key="provided1")
    oob = swap.content.decode().split('id="aln-viewer-actions" hx-swap-oob="true"')[1]
    assert f"/alignments/{curated.id}/" in oob.split("</div>")[0]
    # Deleting the job takes its alignments out of the lists; the copy stays.
    jobs.delete([job])
    assert client.get(f"/alignments/{mafft.id}/").status_code == 404
    page = client.get("/alignments/").content.decode()
    assert "Ime2 run · MAFFT (edited)" in page
    assert "Edit a copy" not in client.get("/alignments/jobs/").content.decode()


def test_a_branch_site_test_registers_nothing(fake_storage: FakeStorage) -> None:
    job = finished({"branch_site": True})
    assert alignments.register_job(job) == []
    running = Job.objects.create(owner=local_user(), params={})
    assert alignments.register_job(running) == []


def test_new_job_from_my_alignments(client: Client, fake_storage: FakeStorage) -> None:
    job = job_with_alignments(fake_storage)
    record = uploaded(client, ALIGNED_FOUR.replace("seq", "s"), name="Mine")
    # "Use in a new job": that alignment is chosen in step 3 (and shown as loaded in step 2
    # by the page's script), and its names fill the outgroup.
    page = client.get(f"/jobs/new/?alignment={record.id}").content.decode()
    assert f'value="{record.id}" checked' in page and 'id="loaded-alignments"' in page
    assert 'value="s0"' in page  # the outgroup checklist
    other = alignments.store(User.objects.create(username="other"), "theirs",
                             alignments.read_rows(ALIGNED_FOUR), "protein", 1,
                             origin=Alignment.Origin.UPLOAD)  # fmt: skip
    page = client.get(f"/jobs/new/?alignment={other.id}&alignment=nonsense").content.decode()
    assert "theirs" not in page and 'name="saved_alignments" value' not in page
    client.get("/alignments/pick/")  # opening the picker registers older jobs' alignments
    mafft = Alignment.objects.get(job=job, job_key="mafft")
    assert mafft.last_used_at is None
    response = post(client, aligners=[], run_trees="", run_ancestors="",
                    saved_alignments=[str(mafft.id)])  # fmt: skip
    assert response.status_code == 302, response.content.decode()[-2000:]
    new = Job.objects.exclude(id=job.id).get()
    assert new.params["alignment_labels"] == ["Ime2 run · MAFFT"]
    assert fake_storage.objects[f"jobs/{new.id}/inputs/alignment1.fasta"] == ALIGNED_FOUR.encode()
    assert f"jobs/{new.id}/inputs/sequences.fasta" not in fake_storage.objects
    mafft.refresh_from_db()
    assert mafft.last_used_at is not None  # first in the picker now


def test_the_picker_lists_recently_used_first(client: Client, fake_storage: FakeStorage) -> None:
    job = job_with_alignments(fake_storage)
    kinase = uploaded(client, name="Kinase domain")
    lysin = uploaded(client, name="Lysin trimmed")

    def rows(query: str = "") -> list[str]:
        page = client.get("/alignments/pick/" + query).content.decode()
        return re.findall(r'data-title="([^"]*)"', page)

    # Never used: newest first (the job finished before these were uploaded).
    assert rows() == ["Lysin trimmed", "Kinase domain", "Ime2 run · MAFFT", "Ime2 run · Curated"]
    client.get(f"/alignments/{kinase.id}/")  # opening one uses it
    assert rows()[0] == "Kinase domain"
    client.get(f"/alignments/{kinase.id}/?sort_toggle=1", HTTP_HX_REQUEST="true")
    client.get(f"/alignments/{lysin.id}/download/")
    assert rows()[:2] == ["Lysin trimmed", "Kinase domain"]
    assert rows("?q=ime2+cur") == ["Ime2 run · Curated"]  # every word, name or job
    assert "No alignment matches" in client.get("/alignments/pick/?q=zzz").content.decode()
    assert job.alignments.count() == 2


def test_the_picker_pages_as_it_scrolls(
    client: Client, fake_storage: FakeStorage, monkeypatch: pytest.MonkeyPatch
) -> None:
    from phylobot import alignment_editor_views

    monkeypatch.setattr(alignment_editor_views, "PICKER_PAGE", 2)
    for n in range(5):
        uploaded(client, name=f"aln {n}")
    page = client.get("/alignments/pick/").content.decode()
    assert re.findall(r'data-title="([^"]*)"', page) == ["aln 4", "aln 3"]
    assert 'hx-trigger="revealed"' in page and "page=1" in page
    last = client.get("/alignments/pick/?page=2").content.decode()
    assert re.findall(r'data-title="([^"]*)"', last) == ["aln 0"] and "revealed" not in last


def test_new_job_labels_and_limits(client: Client, fake_storage: FakeStorage) -> None:
    job = job_with_alignments(fake_storage)
    client.get("/alignments/jobs/")  # registers the job's alignments
    mafft = Alignment.objects.get(job=job, job_key="mafft")
    response = post(
        client, aligners=[], run_trees="", run_ancestors="", saved_alignments=[str(mafft.id)],
        alignment_files=[SimpleUploadedFile("Ime2 run · MAFFT.fasta", ALIGNED_FOUR.encode())],
    )  # fmt: skip
    assert response.status_code == 302, response.content.decode()[-2000:]
    new = Job.objects.exclude(id=job.id).get()
    assert new.params["alignment_labels"] == ["Ime2 run · MAFFT", "Ime2 run · MAFFT (2)"]
    # Another user's alignment is not a choice.
    other = alignments.store(User.objects.create(username="other"), "theirs",
                             alignments.read_rows(ALIGNED_FOUR), "protein", 1,
                             origin=Alignment.Origin.UPLOAD)  # fmt: skip
    response = post(client, aligners=[], run_trees="", run_ancestors="",
                    saved_alignments=[str(other.id)])  # fmt: skip
    assert response.status_code == 200 and "is not one of the available choices" in (
        response.content.decode()
    )
    assert Job.objects.count() == 2


def test_save_overwrites_an_uploaded_alignment(client: Client, fake_storage: FakeStorage) -> None:
    record = uploaded(client)
    page = client.get(f"/alignments/{record.id}/").content.decode()
    assert "Save (overwrite)" in page and "disabled" in page.split("Save (overwrite)")[0][-200:]
    copy_response = client.post(f"/alignments/{record.id}/save/", {"cut": "1", "name": "Child"})
    child = Alignment.objects.get(name="Child")
    assert copy_response.status_code == 302 and child.edits["parent_version"] == 1
    # Overwrite: cut column 4 (all gaps) and remove Alpha, then sort.
    posted = {"mode": "overwrite", "version": "1", "cut": "4", "drop": "3", "sort": "name"}
    response = client.post(f"/alignments/{record.id}/save/", posted)
    assert response["Location"] == f"/alignments/{record.id}/"
    record.refresh_from_db()
    assert fake_storage.objects[record.storage_key].decode() == ">seq2 second\nMKAL\n>seq10\nMRAL\n"
    assert (record.version, record.sequences, record.columns) == (2, 2, 4)
    assert record.edits == {"parent_sites": 5, "kept_sites": "1-3,5",
                            "removed_sequences": ["Alpha"], "sort": "name"}  # fmt: skip
    assert Alignment.objects.count() == 2  # no new alignment
    page = client.get(f"/alignments/{record.id}/").content.decode()
    assert (
        "overwritten 1 time since" in " ".join(page.split())
        and "columns 1-3,5 of the uploaded file's 5" in page
    )
    assert data(page)["rows"] == ["MKAL", "MRAL"] and "|2|" in data(page)["stateKey"]
    # The copy made before now names an earlier version of it.
    assert "an earlier version of" in client.get(f"/alignments/{child.id}/").content.decode()
    # A page opened before the overwrite can't overwrite (its positions name other rows).
    stale = client.post(f"/alignments/{record.id}/save/", posted)
    record.refresh_from_db()
    assert stale.status_code == 302 and record.version == 2
    # Nothing to overwrite without edits.
    client.post(f"/alignments/{record.id}/save/", {"mode": "overwrite", "version": "2"})
    record.refresh_from_db()
    assert record.version == 2


def test_a_jobs_alignment_is_never_overwritten(client: Client, fake_storage: FakeStorage) -> None:
    job = job_with_alignments(fake_storage)
    client.get("/alignments/jobs/")
    mafft = Alignment.objects.get(job=job, job_key="mafft")
    page = client.get(f"/alignments/{mafft.id}/?cut=1").content.decode()
    assert "Save (overwrite)" not in page and "Save a copy" in page
    response = client.post(f"/alignments/{mafft.id}/save/",
                           {"mode": "overwrite", "version": "1", "cut": "1"})  # fmt: skip
    assert response.status_code == 404
    assert fake_storage.objects[mafft.storage_key] == ALIGNED_FOUR.encode()


def test_duplicate(client: Client, fake_storage: FakeStorage) -> None:
    record = uploaded(client, name="Kinase")
    response = client.post(f"/alignments/{record.id}/duplicate/")
    assert response["Location"] == "/alignments/"
    copy = Alignment.objects.exclude(id=record.id).get()
    assert copy.name == "Kinase (copy)" and copy.derived_from == record and copy.edits["duplicate"]
    assert fake_storage.objects[copy.storage_key] == fake_storage.objects[record.storage_key]
    assert (copy.sequences, copy.columns) == (record.sequences, record.columns)
    assert "Copied from" in client.get("/alignments/").content.decode()
    assert "copied from" in client.get(f"/alignments/{copy.id}/").content.decode()
    long = uploaded(client, name="x" * alignments.MAX_NAME)
    client.post(f"/alignments/{long.id}/duplicate/")
    longer = Alignment.objects.get(derived_from=long)
    assert longer.name.endswith("x (copy)") and len(longer.name) == alignments.MAX_NAME
    job = job_with_alignments(fake_storage)
    alignments.ensure_registered(job.owner)
    theirs = Alignment.objects.filter(job=job).first()
    assert theirs is not None
    assert client.post(f"/alignments/{theirs.id}/duplicate/").status_code == 404
