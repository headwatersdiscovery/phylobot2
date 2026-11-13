"""Job import (phylobot.job_import) from the ZIPs phylobot.bundle writes."""

import hashlib
import io
import json
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.http import QueryDict
from django.test import Client, override_settings

from phylobot import bundle, job_import, jobs, manifest, usage
from phylobot.models import Job, Library, LibraryExport, SiteSettings
from phylobot.modes import local_user
from phylobot.settings import middleware_for
from tests.conftest import FakeStorage

pytestmark = pytest.mark.django_db
FINISHED = datetime(2026, 10, 1, 12, 0, tzinfo=UTC)


def make_job(fake_storage: FakeStorage, name: str = "Ime2 fungi", **fields: Any) -> Job:
    fields = {
        "status": Job.Status.SUCCEEDED, "finished_at": FINISHED, "started_at": FINISHED,
        "params": {"aligners": ["mafft"], "seed": 7}, "cpu_seconds": 321.0,
        "billable_cpu_seconds": 321.0, "work_units": 4.5, "progress": {"done": 3, "total": 3},
        **fields,
    }  # fmt: skip
    job = Job.objects.create(owner=local_user(), name=name, **fields)
    prefix = job.storage_prefix
    fake_storage.objects[prefix + "inputs/sequences.fasta"] = b">a\nMK\n>b\nMR\n"
    fake_storage.objects[prefix + "results/mafft/iqtree/LG/tree.nwk"] = b"(a,b);\n"
    fake_storage.objects[prefix + "results/summary.json"] = b'{"combinations": []}'
    return job


def download(client: Client, job: Job) -> bytes:
    response = client.get(f"/jobs/{job.id}/download/")
    return b"".join(response.streaming_content)  # type: ignore[attr-defined]


def upload(client: Client, data: bytes) -> Any:
    return client.post("/jobs/import/", {"archive": SimpleUploadedFile("job.zip", data)})


def rezip(data: bytes, change: Any) -> bytes:
    """The ZIP with its members passed through `change(name, bytes) -> {name: bytes}`."""
    out = io.BytesIO()
    with zipfile.ZipFile(io.BytesIO(data)) as src, zipfile.ZipFile(out, "w") as dst:
        for name in src.namelist():
            for new_name, content in change(name, src.read(name)).items():
                dst.writestr(new_name, content)
    return out.getvalue()


def edit_manifest(data: bytes, edit: Any) -> bytes:
    def change(name: str, content: bytes) -> dict[str, bytes]:
        if name.endswith(bundle.MANIFEST):
            manifest = json.loads(content)
            edit(manifest)
            content = json.dumps(manifest).encode()
        return {name: content}

    return rezip(data, change)


def test_download_has_the_manifest(client: Client, fake_storage: FakeStorage) -> None:
    job = make_job(fake_storage)
    job.add_tags(["ime2", "draft"])
    archive = zipfile.ZipFile(io.BytesIO(download(client, job)))
    folder = bundle.folder_name(job)
    manifest = json.loads(archive.read(f"{folder}/{bundle.MANIFEST}"))
    assert manifest["format"] == "phylobot-job" and manifest["format_version"] == 2
    assert manifest["job"]["name"] == "Ime2 fungi" and manifest["job"]["tags"] == ["draft", "ime2"]
    assert manifest["job"]["params"] == {"aligners": ["mafft"], "seed": 7}
    tree = b"(a,b);\n"
    assert manifest["files"]["results/mafft/iqtree/LG/tree.nwk"] == {
        "sha256": hashlib.sha256(tree).hexdigest(), "size": len(tree),
    }  # fmt: skip
    assert set(manifest["files"]) == {
        "inputs/sequences.fasta", "results/mafft/iqtree/LG/tree.nwk", "results/summary.json",
    }  # fmt: skip


def test_round_trip(client: Client, fake_storage: FakeStorage) -> None:
    job = make_job(fake_storage)
    job.add_tags(["ime2"])
    record = Library.objects.create(
        owner=local_user(), job=job, name="Library 1", code="IME", spec={"seed": 3},
        notes="Order these first.",
    )  # fmt: skip
    LibraryExport.objects.create(
        library=record, number=1, spec={"seed": 3}, spec_version=1, member_count=2,
        phylobot_version="0.1",
    )  # fmt: skip
    fake_storage.objects[f"libraries/{record.id}/exports/1/library.fasta"] = b">m1\nMK\n"
    data = download(client, job)
    jobs.delete([job])  # the original is gone: no name clash, and its files are purged
    assert not fake_storage.list_keys(job.storage_prefix)
    assert not manifest.paths(job)  # and so are their manifest rows

    response = upload(client, data)
    imported = Job.objects.get(imported_at__isnull=False)
    assert response.status_code == 302 and response["Location"] == f"/jobs/{imported.id}/"
    assert imported.id != job.id and imported.name == "Ime2 fungi"
    assert imported.status == Job.Status.SUCCEEDED and imported.finished_at == FINISHED
    assert imported.params == job.params and imported.work_units == 4.5
    # The CPU it used is provenance only, never ours.
    assert imported.cpu_seconds is None and imported.billable_cpu_seconds is None
    assert imported.imported_from["cpu_seconds"] == 321.0
    assert imported.imported_from["id"] == str(job.id)
    assert list(imported.tags.values_list("name", flat=True)) == ["ime2"]
    prefix = imported.storage_prefix
    assert fake_storage.objects[prefix + "results/mafft/iqtree/LG/tree.nwk"] == b"(a,b);\n"
    # Its results are indexed from the ZIP (the files were stored before the Job row).
    stored = {k[len(prefix + "results/") :] for k in fake_storage.list_keys(prefix + "results/")}
    assert manifest.paths(imported) == stored and imported.results_indexed_at is not None
    assert imported.storage_bytes == sum(
        len(fake_storage.objects[k]) for k in fake_storage.list_keys(prefix)
    )
    (copy,) = imported.libraries.all()
    assert (copy.name, copy.code, copy.spec) == ("Library 1", "IME", {"seed": 3})
    assert copy.notes == "Order these first."
    (receipt,) = copy.exports.all()
    assert receipt.member_count == 2
    assert fake_storage.objects[f"libraries/{copy.id}/exports/1/library.fasta"] == b">m1\nMK\n"
    assert not jobs.resumable(imported)

    page = client.get(f"/jobs/{imported.id}/").content.decode()
    assert "imported-banner" in page and "not counted here" in page
    listing = client.get("/jobs/").content.decode()
    assert 'class="badge imported"' in listing and 'id="import-job-open"' in listing

    # Downloaded again, it keeps where it first came from.
    again = json.loads(
        zipfile.ZipFile(io.BytesIO(download(client, imported))).read(
            f"{bundle.folder_name(imported)}/{bundle.MANIFEST}"
        )
    )
    assert again["job"]["imported_from"]["id"] == str(job.id)


def test_a_name_clash_asks_for_a_name(client: Client, fake_storage: FakeStorage) -> None:
    job = make_job(fake_storage)
    data = download(client, job)
    response = upload(client, data)
    assert response.status_code == 302 and "/jobs/import/" in response["Location"]
    assert Job.objects.count() == 1  # nothing imported yet; the upload waits
    page = client.get(response["Location"]).content.decode()
    assert 'value="Ime2 fungi (copy)"' in page

    taken = client.post(response["Location"], {"name": "ime2 FUNGI"}).content.decode()
    assert "You already have a job named" in taken and Job.objects.count() == 1
    done = client.post(response["Location"], {"name": "Ime2 fungi (copy)"})
    imported = Job.objects.get(imported_at__isnull=False)
    assert done["Location"] == f"/jobs/{imported.id}/" and imported.name == "Ime2 fungi (copy)"
    assert not fake_storage.list_keys(job_import.STAGED)
    # The next clash suggests the next free name.
    second = upload(client, data)
    assert 'value="Ime2 fungi (copy 2)"' in client.get(second["Location"]).content.decode()
    cancel = client.post(second["Location"] + "cancel/")
    assert cancel["Location"] == "/jobs/" and not fake_storage.list_keys(job_import.STAGED)
    assert client.get(second["Location"])["Location"] == "/jobs/"  # no longer waiting


def test_another_users_upload_is_not_theirs(client: Client, fake_storage: FakeStorage) -> None:
    bob = User.objects.create_user("bob")
    fake_storage.objects[job_import.staged_key(bob, "a" * 32)] = b"PK"
    assert client.get(f"/jobs/import/{'a' * 32}/")["Location"] == "/jobs/"
    assert client.get("/jobs/import/..%2f..%2fx/").status_code == 404


def refused(client: Client, data: bytes) -> str:
    response = upload(client, data)
    assert response.status_code == 302 and response["Location"] == "/jobs/"
    assert not Job.objects.filter(imported_at__isnull=False).exists()
    page = client.get("/jobs/").content.decode()
    assert "Could not import that ZIP." in page
    return page


@pytest.mark.parametrize(
    ("change", "says"),
    [
        # Paths that would escape the job's folder, or that a job cannot have.
        (lambda n, c: {n.replace("inputs/sequences.fasta", "inputs/../../x"): c}, "unacceptable"),
        (lambda n, c: {n.replace("inputs/sequences.fasta", "inputs/a\\b"): c}, "unacceptable"),
        (lambda n, c: {n.replace("inputs/sequences.fasta", "work/x.sh"): c}, "cannot have"),
        # A file the manifest does not list, or one changed after export.
        (lambda n, c: {n: c, n.rsplit("/", 1)[0] + "/extra.txt": b"x"}
         if n.endswith("summary.json") else {n: c}, "does not match"),
        (lambda n, c: {n: b"(a,c);\n" if n.endswith("tree.nwk") else c}, "changed or damaged"),
        (lambda n, c: {} if n.endswith(bundle.MANIFEST) else {n: c}, "no PhyloBot record"),
    ],
)  # fmt: skip
def test_hostile_zips_are_refused(
    client: Client, fake_storage: FakeStorage, change: Any, says: str
) -> None:
    data = rezip(download(client, make_job(fake_storage)), change)
    jobs.delete(list(Job.objects.all()))
    before = set(fake_storage.objects)
    assert says in refused(client, data)
    assert set(fake_storage.objects) <= before  # nothing left behind


@pytest.mark.parametrize(
    ("edit", "says"),
    [
        (lambda m: m["job"].update(status="running"), "had not finished"),
        (lambda m: m.update(format_version=3), "newer PhyloBot"),
        # Format 1 is the layout before the large-job rework: nothing reads it any more.
        (lambda m: m.update(format_version=1), "older PhyloBot"),
        (lambda m: m.update(format="other"), "not a PhyloBot job"),
        (lambda m: m["job"].update(params=["x"]), "unreadable (settings)"),
        (lambda m: m["job"].update(finished_at="yesterday"), "unreadable (finished at)"),
        (lambda m: m["job"].update(cpu_seconds=-1), "unreadable (cpu seconds)"),
        (lambda m: m.update(libraries="x"), "libraries is unreadable"),
    ],
)
def test_bad_records_are_refused(
    client: Client, fake_storage: FakeStorage, edit: Any, says: str
) -> None:
    data = edit_manifest(download(client, make_job(fake_storage)), edit)
    jobs.delete(list(Job.objects.all()))
    assert says in refused(client, data)


def test_not_a_zip(client: Client, fake_storage: FakeStorage) -> None:
    assert "not a ZIP" in refused(client, b"hello")


def test_bad_tags_are_dropped(client: Client, fake_storage: FakeStorage) -> None:
    data = edit_manifest(
        download(client, make_job(fake_storage)),
        lambda m: m["job"].update(tags=["ok", "Not OK", 3, "<b>"]),
    )
    jobs.delete(list(Job.objects.all()))
    upload(client, data)
    imported = Job.objects.get(imported_at__isnull=False)
    assert list(imported.tags.values_list("name", flat=True)) == ["ok"]


def test_links_to_other_jobs(client: Client, fake_storage: FakeStorage) -> None:
    """A replicate names its simulation: imported first, its link waits for it."""
    simulation = make_job(fake_storage, name="Sim", params={"simulate": True})
    replicate = make_job(fake_storage, name="Rep")
    replicate.params = {"simulated": True, "simulated_from": {"job": str(simulation.id),
                        "seed": 2, "name": "Sim"}}  # fmt: skip
    replicate.save()
    sim_zip, rep_zip = download(client, simulation), download(client, replicate)
    jobs.delete([simulation, replicate])

    upload(client, rep_zip)
    rep = Job.objects.get(imported_at__isnull=False, name="Rep")
    assert "job" not in rep.params["simulated_from"]
    assert rep.params["simulated_from"][job_import.UNLINKED] == str(simulation.id)
    assert client.get(f"/jobs/{rep.id}/").status_code == 200  # names it, unlinked
    assert client.get(f"/jobs/new/?from={rep.id}").status_code == 200

    upload(client, sim_zip)
    sim = Job.objects.get(imported_at__isnull=False, name="Sim")
    rep.refresh_from_db()
    assert rep.params["simulated_from"] == {"job": str(sim.id), "seed": 2, "name": "Sim"}


def test_usage_counts_only_their_storage() -> None:
    now = datetime(2026, 10, 16, 12, tzinfo=UTC)
    owner = local_user()
    imported = Job(
        owner=owner, created_at=datetime(2026, 10, 10, tzinfo=UTC), finished_at=FINISHED,
        storage_bytes=700, imported_at=datetime(2026, 10, 14, tzinfo=UTC),
        imported_from={"cpu_seconds": 999},
    )  # fmt: skip
    window = usage.Filters.parse(QueryDict("range=7d"), now.date(), now.date())
    result = usage.compute([imported], window, now)
    tiles = {t.label: t for t in result.tiles}
    assert tiles["Jobs submitted"].value == "0" and tiles["Storage now"].value == "700 bytes"
    storage = {p.label: p.value for p in result.storage.points}
    assert storage["Tue 13 Oct 2026"] == "0 bytes" and storage["Wed 14 Oct 2026"] == "700 bytes"


@override_settings(PHYLOBOT_MODE="cloud", MIDDLEWARE=middleware_for("cloud"))
def test_cloud_mode_has_no_import(client: Client, fake_storage: FakeStorage) -> None:
    user = User.objects.create_user("ann", password="pw")
    client.force_login(user)
    assert 'id="import-job-open"' not in client.get("/jobs/").content.decode()
    assert upload(client, b"PK").status_code == 404
    assert client.get(f"/jobs/import/{'a' * 32}/").status_code == 404


def test_read_without_the_web(tmp_path: Path, fake_storage: FakeStorage) -> None:
    """`read` checks a file on disk; staging keeps one upload per user."""
    path = tmp_path / "x.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("a/phylobot-job.json", "[]")
    with pytest.raises(job_import.Refused, match="not a PhyloBot job"):
        job_import.read(path)
    owner = local_user()
    first = job_import.stage(owner, path)
    second = job_import.stage(owner, path)
    assert fake_storage.list_keys(job_import.STAGED) == [job_import.staged_key(owner, second)]
    assert not job_import.fetch(owner, first, tmp_path / "y.zip")


@pytest.mark.parametrize(
    ("field", "value", "says"),
    [
        ("import_max_megabytes", 1, "larger than 1 MB"),
        ("import_max_unpacked_megabytes", 1, "add up to more than 1 MB"),
        ("import_max_files", 10, "more than 10 files"),
        ("import_max_record_megabytes", 1, "record is larger than 1 MB"),
        ("import_max_libraries", 0, "more than 0 libraries"),
        ("import_max_exports", 0, "more than 0 exports"),
    ],
)
def test_limits_come_from_the_settings_page(
    client: Client, fake_storage: FakeStorage, field: str, value: int, says: str
) -> None:
    job = make_job(fake_storage)
    record = Library.objects.create(owner=local_user(), job=job, name="L", code="L", spec={})
    LibraryExport.objects.create(
        library=record, number=1, spec={}, spec_version=1, member_count=0,
        phylobot_version="0.1",
    )  # fmt: skip
    big = b"x" * 1_100_000  # past 1 MB uncompressed, and as a stored member past 1 MB
    data = download(client, job)
    if field in {"import_max_megabytes", "import_max_unpacked_megabytes"}:
        name = f"{bundle.folder_name(job)}/results/big.txt"
        data = edit_manifest(
            rezip(data, lambda n, c: {n: c, name: big} if n.endswith("summary.json") else {n: c}),
            lambda m: m["files"].update({"results/big.txt": {
                "sha256": hashlib.sha256(big).hexdigest(), "size": len(big)}}),
        )  # fmt: skip
    if field == "import_max_record_megabytes":
        data = edit_manifest(data, lambda m: m["job"]["params"].update(padding="x" * 1_100_000))
    if field == "import_max_files":
        stub = {"sha256": hashlib.sha256(b"").hexdigest(), "size": 0}
        names = [f"{bundle.folder_name(job)}/results/f{i}.txt" for i in range(10)]
        data = edit_manifest(
            rezip(data, lambda n, c: {n: c} | (dict.fromkeys(names, b"") if n.endswith(".json")
                                               and "/results/" in n else {})),
            lambda m: m["files"].update({f"results/f{i}.txt": stub for i in range(10)}),
        )  # fmt: skip
    jobs.delete([job])
    SiteSettings.objects.create(pk=1, **{field: value})
    assert says in refused(client, data)
    SiteSettings.objects.all().delete()  # the defaults take it
    assert upload(client, data)["Location"].startswith("/jobs/")
    assert Job.objects.filter(imported_at__isnull=False).exists()
    assert "Up to 2 GB" in client.get("/jobs/").content.decode()
