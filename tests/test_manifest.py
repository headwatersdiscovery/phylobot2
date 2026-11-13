"""Each job's result files in Postgres (phylobot.manifest): written where results are
stored, read by every page instead of a storage listing."""

import json
from pathlib import Path
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.test import Client
from django.utils import timezone

from phylobot import filelist, jobs, manifest, storage
from phylobot.models import Job, ResultFile
from phylobot.modes import local_user
from tests.conftest import FakeStorage

pytestmark = pytest.mark.django_db


def finished_job(**fields: Any) -> Job:
    job = Job.objects.create(owner=local_user(), **fields)
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    job.refresh_from_db()
    return job


def test_directory_sizes_name_files_as_uploaded(tmp_path: Path) -> None:
    (tmp_path / "mafft" / "iqtree").mkdir(parents=True)
    (tmp_path / "summary.json").write_text("{}")
    (tmp_path / "mafft" / "iqtree" / "tree.nwk").write_text("(a,b);\n")
    assert manifest.directory_sizes(tmp_path) == {"mafft/iqtree/tree.nwk": 7, "summary.json": 2}


def test_an_upload_adds_and_resizes_never_removes() -> None:
    """A resumed attempt uploads over the first one's files and never deletes any, so
    neither does its manifest."""
    job = finished_job()
    manifest.indexed(job.id, {"a.tsv": 1, "logs/pipeline.log": 5})
    manifest.unindex(job.id)
    job.refresh_from_db()
    assert job.results_indexed_at is None
    manifest.indexed(job.id, {"logs/pipeline.log": 9, "b.tsv": 2})
    job.refresh_from_db()
    assert job.results_indexed_at is not None
    assert manifest.sizes(job) == {"a.tsv": 1, "logs/pipeline.log": 9, "b.tsv": 2}
    assert manifest.paths(job, "logs/") == {"logs/pipeline.log"}
    assert manifest.has(job, "b.tsv") and not manifest.has(job, "c.tsv")


def test_maintenance_indexes_older_jobs_from_storage(fake_storage: FakeStorage) -> None:
    """Jobs from before the manifest, or whose upload broke off, are listed once."""
    fake_storage.objects.indexing = False  # stored with no manifest
    old = finished_job()
    running = Job.objects.create(owner=local_user())
    assert jobs.mark(running.id, Job.Status.RUNNING)
    deleted = finished_job(deleted_at=timezone.now())
    for job in (old, running, deleted):
        fake_storage.objects[f"{job.storage_prefix}results/mafft/alignment.fasta"] = b">a\nMK\n"
    fake_storage.objects[f"{old.storage_prefix}inputs/sequences.fasta"] = b">a\nMK\n"

    assert manifest.index_unindexed() == 1
    assert manifest.sizes(old) == {"mafft/alignment.fasta": 6}  # results/ only
    assert not manifest.paths(running) and not manifest.paths(deleted)
    assert manifest.index_unindexed() == 0  # once


def test_purging_forgets_the_files(fake_storage: FakeStorage) -> None:
    job = finished_job()
    fake_storage.objects[f"{job.storage_prefix}results/summary.json"] = b"{}"
    assert manifest.paths(job) == {"summary.json"}
    jobs.delete([job])
    assert not ResultFile.objects.filter(job=job).exists()


def no_listing(monkeypatch: pytest.MonkeyPatch, fake_storage: FakeStorage) -> None:
    """Listing a job's results/ fails: pages must read the manifest."""
    listed = fake_storage.list_keys

    def list_keys(prefix: str) -> list[str]:
        assert "/results" not in prefix, f"listed {prefix}"
        return listed(prefix)

    monkeypatch.setattr(storage, "list_keys", list_keys)
    monkeypatch.setattr(storage, "list_sizes", lambda prefix: pytest.fail(f"listed {prefix}"))


def large_job(fake_storage: FakeStorage, combinations: int = 3, nodes: int = 60) -> Job:
    """A job with more files than the page lists: properties per node, as on GFP."""
    job = finished_job()
    prefix = f"{job.storage_prefix}results/"
    combos = [
        {"key": f"mafft/iqtree/M{n}", "aligner": "mafft", "tree": "iqtree", "model": f"M{n}",
         "label": f"MAFFT · IQ-TREE · M{n}"}
        for n in range(combinations)
    ]  # fmt: skip
    summary = {"alignments": [], "combinations": combos}
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    fake_storage.objects[prefix + "mafft/alignment.fasta"] = b">a\nMK\n"
    for combo in combos:
        fake_storage.objects[f"{prefix}{combo['key']}/tree.nwk"] = b"(a,b);\n"
        for node in range(nodes):
            for name in ("draws.tsv", "summary.json"):
                key = f"{prefix}{combo['key']}/properties/Node{node}/parsimony/{name}"
                fake_storage.objects[key] = b"x" * 10
    return job


def test_a_large_job_page_shows_folders_without_listing_storage(
    client: Client, fake_storage: FakeStorage, monkeypatch: pytest.MonkeyPatch
) -> None:
    job = large_job(fake_storage, combinations=4)  # 4 x (60 x 2 + 1) + 2 = 486 files
    assert len(manifest.paths(job)) > filelist.INLINE_FILES
    no_listing(monkeypatch, fake_storage)
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert f'href="/jobs/{job.id}/folders/mafft/"' in page
    assert "4 combinations, 485 files" in page
    assert "Node0" not in page  # no file of the 482 under mafft/ is listed

    # Down a level at a time, until the files are few enough to list.
    page = client.get(f"/jobs/{job.id}/folders/mafft/").content.decode()
    assert f'href="/jobs/{job.id}/folders/mafft/iqtree/"' in page
    assert f'href="/jobs/{job.id}/files/mafft/alignment.fasta"' in page
    page = client.get(f"/jobs/{job.id}/folders/mafft/iqtree/M1/").content.decode()
    assert "MAFFT · IQ-TREE · M1" in page and "properties/Node59/parsimony/draws.tsv" in page
    assert client.get(f"/jobs/{job.id}/folders/mafft/nothing/").status_code == 404
    assert client.get(f"/jobs/{job.id}/folders/../").status_code == 404


def test_a_small_job_lists_every_file(client: Client, fake_storage: FakeStorage) -> None:
    job = large_job(fake_storage, combinations=1, nodes=2)
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "folders/" not in page
    assert f'href="/jobs/{job.id}/files/mafft/iqtree/M0/tree.nwk"' in page


def test_others_folders_are_not_found(client: Client, fake_storage: FakeStorage) -> None:
    job = large_job(fake_storage, combinations=1, nodes=2)
    other = User.objects.create(username="other")
    Job.objects.filter(id=job.id).update(owner=other)
    assert client.get(f"/jobs/{job.id}/folders/mafft/").status_code == 404


def test_a_finished_job_not_indexed_yet_says_so(client: Client, fake_storage: FakeStorage) -> None:
    fake_storage.objects.indexing = False
    job = finished_job()
    fake_storage.objects[f"{job.storage_prefix}results/summary.json"] = b"{}"
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Listing this job's files" in page


def test_a_running_job_shows_its_live_log(
    client: Client, fake_storage: FakeStorage, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Uploaded during the run, so not in the manifest: found by its known path."""
    job = Job.objects.create(owner=local_user())
    assert jobs.mark(job.id, Job.Status.RUNNING)
    fake_storage.objects.indexing = False
    fake_storage.objects[f"{job.storage_prefix}results/logs/pipeline.log"] = b"step one\n"
    no_listing(monkeypatch, fake_storage)
    page = client.get(f"/jobs/{job.id}/status/").content.decode()
    assert "step one" in page
