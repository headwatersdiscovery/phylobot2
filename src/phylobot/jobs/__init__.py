"""The job queue. All job status changes go through this module.

Backend: the Job table itself is the queue (Postgres, SELECT ... FOR UPDATE SKIP
LOCKED). A later SQS or AWS Batch backend replaces the functions here, not callers.
Rules that hold for any backend:
- a queue message is only a job ID; data always comes from the Job row;
- status only moves forward (see TRANSITIONS), so repeated or out-of-order updates
  are harmless: `mark` returns False and changes nothing;
- whoever runs a job refreshes `heartbeat_at`; `sweep_stale` fails jobs whose runner died;
- deleting a job cancels it if unfinished and removes its files at once, but keeps the
  record (`deleted_at`, `purged_at`, runtimes, `cpu_seconds`) for accounting. A running
  job's worker stops at its next heartbeat (refused), uploads nothing more, and purges
  again in case a live-log upload was in flight;
- stopping a job (`stop`, from the job page or My jobs) cancels it and keeps its files;
  its runner stops the containers at the next heartbeat. Cancelled is terminal, so a
  stopped job cannot be resumed -- only duplicated and run again;
- "duplicate and run again" never rewinds a job's status: it starts a new job from the
  same inputs and settings and leaves the old one as it is;
- resuming a failed job (`resume`) is the one way back to queued: it starts the job's
  next attempt (Job.attempt, a JobAttempt row each), which reuses the steps that
  finished. Within an attempt status still only moves forward, and every write from a
  runner (`mark`, `heartbeat`, `record_usage`) names its attempt, so a superseded runner
  (one that slept through being declared stopped, say) changes nothing.
"""

import logging
import uuid

from collections.abc import Collection
from datetime import timedelta
from typing import Any

from django.db import transaction
from django.db.models import Count, F, Q, Sum
from django.utils import timezone

from phylobot import manifest, site_settings, storage
from phylobot.models import Job, JobAttempt, limits_for

logger = logging.getLogger(__name__)

Status = Job.Status
TERMINAL = {Status.SUCCEEDED, Status.FAILED, Status.CANCELLED}
STOPPED = "The job stopped unexpectedly. You can resume it."  # its worker went silent
INPUT_ERROR = "Input error:"  # VALIDATE's messages (executors.user_error): resuming repeats them
TRANSITIONS: dict[str, set[str]] = {
    Status.QUEUED: {Status.RUNNING, Status.FAILED, Status.CANCELLED},
    Status.RUNNING: {Status.SUCCEEDED, Status.FAILED, Status.CANCELLED},
    Status.SUCCEEDED: set(),
    Status.FAILED: set(),
    Status.CANCELLED: set(),
}


def enqueue(job: Job) -> None:
    """Make a saved, queued job visible to workers. With the table as the queue, the
    committed row is the message; other backends send the job ID here (after commit)."""
    assert job.status == Status.QUEUED


def claim_next(worker: str) -> Job | None:
    """Atomically take the oldest queued job whose owner is under their concurrency
    limit, and mark it running. Safe with any number of concurrent workers."""
    with transaction.atomic():
        candidates = (
            Job.objects.select_for_update(skip_locked=True, of=("self",))
            .filter(status=Status.QUEUED)
            .select_related("owner")
            .order_by("created_at")[:50]
        )
        for job in candidates:
            running = Job.objects.filter(owner_id=job.owner_id, status=Status.RUNNING).count()
            if running >= limits_for(job.owner).max_concurrent_jobs:
                continue
            now = timezone.now()
            job.status = Status.RUNNING
            job.worker = worker
            job.started_at = job.heartbeat_at = now
            job.save(update_fields=["status", "worker", "started_at", "heartbeat_at"])
            JobAttempt.objects.update_or_create(
                job=job, number=job.attempt, defaults={"worker": worker, "started_at": now}
            )
            return job
    return None


def mark(job_id: Any, status: str, error: str = "", attempt: int | None = None) -> bool:
    """Move a job to `status` if that is a legal forward step; returns whether it moved.
    A runner passes its `attempt`, and nothing moves once the job has gone on to another;
    deleting and sweeping act on whichever attempt is current."""
    allowed_from = [source for source, targets in TRANSITIONS.items() if status in targets]
    now = timezone.now()
    fields: dict[str, Any] = {"status": status}
    if status in TERMINAL:
        fields["finished_at"] = now
    if error:
        fields["error"] = error
    with transaction.atomic():
        rows = Job.objects.select_for_update().filter(id=job_id, status__in=allowed_from)
        if attempt is not None:
            rows = rows.filter(attempt=attempt)
        current = rows.values_list("attempt", flat=True).first()
        if current is None:
            return False
        Job.objects.filter(id=job_id).update(**fields)
        if status in TERMINAL:
            JobAttempt.objects.filter(job_id=job_id, number=current).update(
                status=status, finished_at=now, error=error
            )
    return True


def heartbeat(
    job_id: Any, progress: dict[str, Any] | None = None, attempt: int | None = None
) -> bool:
    """Record that a running job's runner is alive (and optionally its progress). False
    once the job is no longer running, or no longer on the runner's `attempt`."""
    fields: dict[str, Any] = {"heartbeat_at": timezone.now()}
    if progress is not None:
        fields["progress"] = progress
    rows = Job.objects.filter(id=job_id, status=Status.RUNNING)
    if attempt is not None:
        rows = rows.filter(attempt=attempt)
    return rows.update(**fields) == 1


def sweep_stale(timeout: timedelta = timedelta(minutes=10), alive: Collection[str] = ()) -> int:
    """Fail running jobs whose runner has not sent a heartbeat within `timeout`. `alive`:
    the sweeping worker's own running jobs, never swept (after the machine sleeps, the
    sweep can come round before their next heartbeat)."""
    cutoff = timezone.now() - timeout
    stale = Job.objects.filter(status=Status.RUNNING).filter(
        Q(heartbeat_at__lt=cutoff) | Q(heartbeat_at__isnull=True)
    )
    count = 0
    for job_id in stale.values_list("id", flat=True):
        if str(job_id) not in alive:
            count += mark(job_id, Status.FAILED, STOPPED)
    return count


def record_usage(
    job_id: Any,
    attempt: int,
    completed: float,
    failed: float,
    run: int | None = None,
    reused: int | None = None,
) -> None:
    """Record an attempt's compute (CPU seconds of steps completed and of steps failed or
    stopped), billing failed steps at the admin's current rate, and total the job's."""
    rate = site_settings.failed_task_rate_percent()
    JobAttempt.objects.update_or_create(
        job_id=job_id,
        number=attempt,
        defaults={
            "completed_cpu_seconds": completed,
            "failed_cpu_seconds": failed,
            "tasks_run": run,
            "tasks_reused": reused,
            "failed_rate_percent": rate,
            "billable_cpu_seconds": completed + failed * rate / 100,
        },
    )
    totals = JobAttempt.objects.filter(job_id=job_id, completed_cpu_seconds__isnull=False)
    figures = totals.aggregate(
        used=Sum(F("completed_cpu_seconds") + F("failed_cpu_seconds")),
        billable=Sum("billable_cpu_seconds"),
    )
    Job.objects.filter(id=job_id).update(
        cpu_seconds=figures["used"], billable_cpu_seconds=figures["billable"]
    )


def is_job_id(name: str) -> bool:
    """Whether a folder or key name is a job id (a UUID)."""
    try:
        uuid.UUID(name)
    except ValueError:
        return False
    return True


def resumable(job: Job) -> bool:
    """Whether the job can be resumed: failed, not deleted, not over its input (which
    would only fail again), and run here (an imported job has no work folder here, and
    its settings were never checked by this server)."""
    return (
        job.status == Status.FAILED
        and job.deleted_at is None
        and job.imported_at is None
        and not job.error.startswith(INPUT_ERROR)
    )


STOPPED_BY_OWNER = "Stopped by its owner."


def stop(job: Job) -> bool:
    """Stop a queued or running job at its owner's request. Returns whether it moved;
    refuses quietly if it had already finished (or someone got there first).

    The job's own runner notices at its next heartbeat and stops its containers, so the
    machine is freed within `settings.PHYLOBOT_HEARTBEAT_SECONDS`, not at once. A stopped
    job is `cancelled`, and `resumable` accepts only `failed`, so **this cannot be
    continued afterwards** -- say so wherever it is offered. Its files are kept (nothing
    is deleted here); to carry on from its inputs, duplicate it and run again.
    """
    return mark(job.id, Status.CANCELLED, STOPPED_BY_OWNER)


def resume(job: Job) -> bool:
    """Queue a failed job's next attempt. The worker runs it in the job's work directory
    with Nextflow's -resume, so steps that finished are reused (and billed only once).
    Returns whether it was queued; refuses quietly if someone got there first."""
    with transaction.atomic():
        locked = Job.objects.select_for_update().filter(id=job.id).first()
        if locked is None or not resumable(locked):
            return False
        locked.status = Status.QUEUED
        locked.attempt += 1
        locked.error = ""
        locked.progress = {}
        locked.worker = ""
        locked.finished_at = locked.heartbeat_at = None
        locked.work_bytes = locked.work_expired_at = None
        locked.save(
            update_fields=[
                "status", "attempt", "error", "progress", "worker", "finished_at",
                "heartbeat_at", "work_bytes", "work_expired_at",
            ]
        )  # fmt: skip
    job.refresh_from_db()
    enqueue(job)
    return True


def measure_storage(job: Job) -> bool:
    """Record how much object storage a finished job's files take (`storage_bytes`)."""
    try:
        size = storage.prefix_size(job.storage_prefix)
    except Exception:
        logger.warning("job %s: could not measure its storage", job.id, exc_info=True)
        return False
    # Not once purged: the files are gone, and the size measured before is the record.
    return bool(Job.objects.filter(id=job.id, purged_at__isnull=True).update(storage_bytes=size))


def measure_unmeasured(limit: int = 50) -> int:
    """Measure finished jobs not measured yet (older jobs; worker maintenance)."""
    pending = Job.objects.filter(
        storage_bytes__isnull=True, deleted_at__isnull=True, status__in=TERMINAL
    ).order_by("-created_at")[:limit]
    return sum(measure_storage(job) for job in pending)


def delete(jobs: list[Job]) -> int:
    """Delete jobs for their owner: cancel unfinished ones, mark all deleted, and remove
    their files."""
    now = timezone.now()
    for job in jobs:
        if job.status not in TERMINAL:
            mark(job.id, Status.CANCELLED, "Deleted by its owner.")
        Job.objects.filter(id=job.id, deleted_at__isnull=True).update(deleted_at=now)
    for job in Job.objects.filter(id__in=[j.id for j in jobs], status__in=TERMINAL):
        purge(job)
    return len(jobs)


MAX_NAME = 100  # Job.name


def rename(job: Job, name: str) -> None:
    """Rename a job (its owner's label; blank is allowed and shows as "Job")."""
    job.name = " ".join(name.split())[:MAX_NAME]
    Job.objects.filter(id=job.id).update(name=job.name)


def purge(job: Job) -> bool:
    """Remove a deleted, finished job's files from storage; record when. Idempotent.

    A library's exports live under their own prefix, not the job's, so they have to be
    named here: the ancestors they were drawn from are going, and an export that outlived
    them would promise a provenance nothing could check. The `Library` and `LibraryExport`
    records stay, as deleted jobs' own records do.
    """
    if job.deleted_at is None or job.status not in TERMINAL:
        return False
    try:
        storage.delete_prefix(job.storage_prefix)
        storage.delete_prefix(f"work-archives/{job.id}/")  # phylobot.work_archive
        for prefix in job.libraries.values_list("id", flat=True):
            storage.delete_prefix(f"libraries/{prefix}/")
    except Exception:
        logger.exception("job %s: could not delete files; will retry", job.id)
        return False
    manifest.forget(job.id)
    Job.objects.filter(id=job.id).update(purged_at=timezone.now())
    return True


def purge_deleted() -> int:
    """Retry file removal for deleted jobs whose purge failed (worker maintenance)."""
    pending = Job.objects.filter(
        deleted_at__isnull=False, purged_at__isnull=True, status__in=TERMINAL
    )
    return sum(purge(job) for job in pending)


def running_counts() -> dict[int, int]:
    """Running jobs per owner id (for admin and diagnostics)."""
    rows = Job.objects.filter(status=Status.RUNNING).values("owner_id").annotate(n=Count("id"))
    return {row["owner_id"]: row["n"] for row in rows}
