"""Failed jobs' work folders, archived to object storage so a resume survives losing the
worker's disk (Garage locally, S3 in the cloud; site_settings.archive_failed_work).

A failed job's local folder (jobs/<id>/: Nextflow's resume cache and every step's files)
is packed into one tar per attempt, `work-archives/<id>/<attempt>.tar.gz`, by the
worker's archiver thread. A resume on a machine without the folder unpacks it first.

Nextflow reuses a step only if its inputs keep their path, size and modification time:
a tar keeps the times (copying files one by one to S3 would not), and the folder comes
back at the same path, so the worker must keep the checkout at a fixed place. Restored
elsewhere, the job still runs, but every step again.

The archive is ours, but treat it as untrusted when unpacking: members must stay inside
the job's own folder, and Nextflow's absolute links (to other steps' files) are turned
into relative ones inside it, or dropped.
"""

import logging
import os
import shutil
import tarfile
import uuid
from pathlib import Path

from django.db.models import F

from phylobot import jobs, site_settings, storage
from phylobot.models import Job

logger = logging.getLogger("phylobot.worker")
PREFIX = "work-archives/"


def key(job_id: object, attempt: int) -> str:
    return f"{PREFIX}{job_id}/{attempt}.tar.gz"


def pack(folder: Path, dest: Path) -> int:
    """Pack `folder` into `dest` (members named <folder name>/...); its size in bytes.
    Links are stored as links, never followed."""
    with tarfile.open(dest, "w:gz", compresslevel=1) as tar:
        tar.add(folder, arcname=folder.name)
    return dest.stat().st_size


def unpack(archive: Path, jobs_dir: Path, name: str) -> None:
    """Unpack an archive made by `pack` into jobs_dir/<name>, keeping file times."""
    folder = (jobs_dir / name).resolve()

    def keep(member: tarfile.TarInfo, path: str) -> tarfile.TarInfo | None:
        if member.name != name and not member.name.startswith(f"{name}/"):
            raise tarfile.FilterError(f"{member.name} is outside the job's folder")
        if member.issym() and os.path.isabs(member.linkname):
            target = Path(member.linkname)
            if not target.is_relative_to(folder):
                return None  # restored at another path: the step runs again
            start = (Path(path) / member.name).parent
            member = member.replace(linkname=os.path.relpath(target, start), deep=False)
        return tarfile.data_filter(member, path)

    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(jobs_dir, filter=keep)


def archive(root: Path, job_id: uuid.UUID | str, attempt: int) -> bool:
    """Archive a failed attempt's folder, if it is still that attempt's; older archives
    of the job are removed once the new one is recorded."""
    folder = root / "jobs" / str(job_id)
    temp = root / "jobs" / f".archive-{job_id}.tar.gz"
    try:
        size = pack(folder, temp)
        storage.put_file(key(job_id, attempt), temp)
    except Exception:
        logger.warning("job %s: could not archive its work folder", job_id, exc_info=True)
        return False
    finally:
        temp.unlink(missing_ok=True)
    recorded = Job.objects.filter(
        id=job_id, attempt=attempt, status=jobs.Status.FAILED, deleted_at__isnull=True,
        work_expired_at__isnull=True,
    ).update(work_archive_attempt=attempt, work_archive_bytes=size)  # fmt: skip
    for stored in storage.list_keys(f"{PREFIX}{job_id}/"):
        # Resumed meanwhile, the new one may hold a half-changed folder: drop it.
        if stored != key(job_id, attempt) or not recorded:
            storage.delete_key(stored)
    if recorded:
        logger.info("job %s: archived attempt %d's work folder (%d bytes)", job_id, attempt, size)
    return bool(recorded)


def archive_pending(root: Path, running: set[str]) -> int:
    """Archive failed jobs' folders whose latest attempt has no archive yet, if the
    Settings page says so (worker maintenance, in its own thread: big folders take
    minutes)."""
    if not site_settings.archive_failed_work():
        return 0
    folders = {p.name for p in (root / "jobs").glob("*") if p.is_dir()}
    pending = (
        Job.objects.filter(
            status=jobs.Status.FAILED,
            deleted_at__isnull=True,
            work_expired_at__isnull=True,
            id__in=[n for n in folders - running if jobs.is_job_id(n)],
        )  # fmt: skip
        .exclude(work_archive_attempt=F("attempt"))
        .values_list("id", "attempt")
    )
    return sum(archive(root, job_id, attempt) for job_id, attempt in pending)


def restore(root: Path, job: Job) -> bool:
    """Bring back a resumed job's folder from its archive, before its next attempt."""
    if job.work_archive_attempt is None:
        return False
    jobs_dir = root / "jobs"
    folder = jobs_dir / str(job.id)
    temp = jobs_dir / f".restore-{job.id}.tar.gz"
    try:
        storage.get_file(key(job.id, job.work_archive_attempt), temp)
        shutil.rmtree(folder, ignore_errors=True)
        unpack(temp, jobs_dir, str(job.id))
    except Exception:
        logger.warning("job %s: could not restore its work folder", job.id, exc_info=True)
        shutil.rmtree(folder, ignore_errors=True)  # start from the beginning
        return False
    finally:
        temp.unlink(missing_ok=True)
    logger.info("job %s: restored its folder from attempt %d", job.id, job.work_archive_attempt)
    return True


def remove(job_id: object) -> None:
    storage.delete_prefix(f"{PREFIX}{job_id}/")


def tidy() -> int:
    """Remove the archives of jobs that no longer need them: succeeded or cancelled
    after a resume, or deleted (worker maintenance)."""
    done = Job.objects.filter(work_archive_attempt__isnull=False).exclude(
        status__in=[jobs.Status.FAILED, jobs.Status.QUEUED, jobs.Status.RUNNING],
        deleted_at__isnull=True,
    )
    count = 0
    for job_id in done.values_list("id", flat=True):
        remove(job_id)
        Job.objects.filter(id=job_id).update(work_archive_attempt=None, work_archive_bytes=None)
        count += 1
    return count
