"""Failed jobs' work folders archived to object storage (phylobot.work_archive): packed
with their file times, restored at the same path, links kept inside the job's folder,
and removed once nothing needs them."""

import os
import shutil
import tarfile
from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest
from django.contrib.auth.models import User

from phylobot import executors, jobs, work_archive
from phylobot.models import Job, SiteSettings
from tests.conftest import FakeStorage

pytestmark = pytest.mark.django_db
STEP_TIME = 1_600_000_000  # a modification time Nextflow's cache keys would include


def failed_job(root: Path, attempt: int = 1) -> tuple[Job, Path]:
    job = Job.objects.create(
        owner=User.objects.get_or_create(username="alice")[0], status=Job.Status.FAILED,
        attempt=attempt, finished_at=datetime.now(UTC), resume_cache=True,
    )  # fmt: skip
    folder = root / "jobs" / str(job.id)
    step = folder / "work" / "ab" / "123"
    step.mkdir(parents=True)
    (step / "out.fasta").write_text(">a\nMK\n")
    os.utime(step / "out.fasta", (STEP_TIME, STEP_TIME))
    later = folder / "work" / "cd" / "456"
    later.mkdir(parents=True)
    (later / "in.fasta").symlink_to(step / "out.fasta")  # Nextflow stages absolute links
    (folder / "params.json").write_text("{}")
    (folder / ".nextflow").mkdir()
    return job, folder


def switch_on() -> None:
    SiteSettings.objects.update_or_create(pk=1, defaults={"archive_failed_work": True})


def test_an_archive_restores_only_its_own_job(tmp_path: Path) -> None:
    job, folder = failed_job(tmp_path)
    archive = tmp_path / "a.tar.gz"
    work_archive.pack(folder, archive)
    with pytest.raises(tarfile.FilterError):
        work_archive.unpack(archive, tmp_path / "elsewhere", "other-name")


def test_restored_at_the_same_path(tmp_path: Path) -> None:
    job, folder = failed_job(tmp_path)
    archive = tmp_path / "a.tar.gz"
    work_archive.pack(folder, archive)
    shutil.rmtree(folder)
    work_archive.unpack(archive, tmp_path / "jobs", str(job.id))
    out = folder / "work" / "ab" / "123" / "out.fasta"
    link = folder / "work" / "cd" / "456" / "in.fasta"
    assert out.stat().st_mtime == STEP_TIME
    assert link.is_symlink() and not os.path.isabs(os.readlink(link))
    assert link.resolve() == out.resolve() and (folder / ".nextflow").is_dir()


def test_links_outside_the_folder_are_dropped_and_escapes_refused(tmp_path: Path) -> None:
    job, folder = failed_job(tmp_path)
    (folder / "outside").symlink_to("/etc/passwd")
    archive = tmp_path / "a.tar.gz"
    work_archive.pack(folder, archive)
    shutil.rmtree(folder)
    work_archive.unpack(archive, tmp_path / "jobs", str(job.id))
    assert not (folder / "outside").exists() and not (folder / "outside").is_symlink()
    evil = tmp_path / "evil.tar.gz"
    with tarfile.open(evil, "w:gz") as tar:
        info = tarfile.TarInfo(f"{job.id}/../../escaped")
        tar.addfile(info)
    with pytest.raises(tarfile.FilterError):
        work_archive.unpack(evil, tmp_path / "jobs", str(job.id))
    assert not (tmp_path / "escaped").exists()


def test_archived_only_when_switched_on_then_restored(
    tmp_path: Path, fake_storage: FakeStorage
) -> None:
    job, folder = failed_job(tmp_path)
    assert work_archive.archive_pending(tmp_path, set()) == 0  # off by default
    switch_on()
    assert work_archive.archive_pending(tmp_path, {str(job.id)}) == 0  # running: not yet
    assert work_archive.archive_pending(tmp_path, set()) == 1
    assert work_archive.archive_pending(tmp_path, set()) == 0  # this attempt is done
    job.refresh_from_db()
    assert job.work_archive_attempt == 1 and job.work_archive_bytes
    assert list(fake_storage.objects) == [f"work-archives/{job.id}/1.tar.gz"]
    shutil.rmtree(folder)  # the worker's disk is lost
    assert jobs.resume(job)
    job.refresh_from_db()
    assert job.work_archive_attempt == 1  # kept through the resume, to restore from
    assert work_archive.restore(tmp_path, job)
    assert (folder / "work" / "ab" / "123" / "out.fasta").stat().st_mtime == STEP_TIME


def test_a_later_attempt_replaces_the_archive(tmp_path: Path, fake_storage: FakeStorage) -> None:
    switch_on()
    job, _ = failed_job(tmp_path)
    assert work_archive.archive(tmp_path, job.id, 1)
    Job.objects.filter(id=job.id).update(attempt=2)  # resumed, failed again
    assert work_archive.archive(tmp_path, job.id, 2)
    assert list(fake_storage.objects) == [f"work-archives/{job.id}/2.tar.gz"]
    # Resumed while archiving: the archive is not recorded, and is dropped.
    Job.objects.filter(id=job.id).update(attempt=3, status=Job.Status.RUNNING)
    assert not work_archive.archive(tmp_path, job.id, 2)


def test_restore_failure_starts_from_the_beginning(
    tmp_path: Path, fake_storage: FakeStorage
) -> None:
    job, folder = failed_job(tmp_path)
    shutil.rmtree(folder)
    job.work_archive_attempt = 1  # recorded, but the object is missing
    assert not work_archive.restore(tmp_path, job)
    assert not folder.exists()


def test_archives_go_when_nothing_needs_them(tmp_path: Path, fake_storage: FakeStorage) -> None:
    switch_on()
    done, _ = failed_job(tmp_path)
    expired, _ = failed_job(tmp_path)
    kept, _ = failed_job(tmp_path)
    for job in (done, expired, kept):
        assert work_archive.archive(tmp_path, job.id, 1)
    Job.objects.filter(id=done.id).update(status=Job.Status.SUCCEEDED)  # resumed and done
    assert work_archive.tidy() == 1
    SiteSettings.objects.filter(pk=1).update(keep_failed_work_days=7)
    Job.objects.filter(id=expired.id).update(finished_at=datetime.now(UTC) - timedelta(days=8))
    shutil.rmtree(tmp_path / "jobs" / str(expired.id))  # only the archive is left
    assert executors.expire_work_folders(tmp_path, set()) == 1
    assert sorted(fake_storage.objects) == [f"work-archives/{kept.id}/1.tar.gz"]
    expired.refresh_from_db()
    assert expired.work_expired_at and expired.work_archive_attempt is None
    Job.objects.filter(id=kept.id).update(deleted_at=datetime.now(UTC))
    kept.refresh_from_db()
    assert jobs.purge(kept) and fake_storage.objects == {}
