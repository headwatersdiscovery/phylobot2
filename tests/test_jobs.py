from datetime import timedelta

import pytest
from django.contrib.auth import get_user_model
from django.contrib.auth.models import AnonymousUser, User
from django.test import override_settings
from django.utils import timezone

from phylobot import jobs
from phylobot.models import Job, SiteSettings, UserLimits, limits_for

Status = Job.Status
pytestmark = pytest.mark.django_db


def user(name: str) -> User:
    created: User = get_user_model().objects.create_user(name)
    return created


def test_limits_default_then_override() -> None:
    alice = user("alice")
    with override_settings(
        PHYLOBOT_DEFAULT_LIMITS={
            "max_sequences": 10,
            "max_length": 20,
            "max_upload_bytes": 30,
            "max_concurrent_jobs": 1,
        }
    ):
        assert limits_for(alice).max_sequences == 10
        UserLimits.objects.create(user=alice, max_sequences=99)  # other fields blank
        limits = limits_for(alice)
        assert (limits.max_sequences, limits.max_length) == (99, 20)


def test_jobs_are_private() -> None:
    alice, bob = user("alice"), user("bob")
    job = Job.objects.create(owner=alice)
    assert list(Job.objects.for_user(alice)) == [job]
    assert not Job.objects.for_user(bob).exists()
    assert not Job.objects.for_user(AnonymousUser()).exists()


def test_claim_takes_oldest_and_never_twice() -> None:
    alice, bob = user("alice"), user("bob")
    first = Job.objects.create(owner=alice)
    second = Job.objects.create(owner=bob)
    claimed = jobs.claim_next("w1")
    assert claimed is not None and claimed.id == first.id
    assert claimed.status == Status.RUNNING and claimed.worker == "w1"
    again = jobs.claim_next("w2")
    assert again is not None and again.id == second.id
    assert jobs.claim_next("w3") is None


def test_claim_respects_concurrency_limit() -> None:
    alice, bob = user("alice"), user("bob")
    UserLimits.objects.create(user=alice, max_concurrent_jobs=1)
    a1 = Job.objects.create(owner=alice)
    a2 = Job.objects.create(owner=alice)
    b1 = Job.objects.create(owner=bob)
    assert jobs.claim_next("w").id == a1.id  # type: ignore[union-attr]
    assert jobs.claim_next("w").id == b1.id  # type: ignore[union-attr]  # a2 waits
    assert jobs.claim_next("w") is None
    assert jobs.mark(a1.id, Status.SUCCEEDED)
    assert jobs.claim_next("w").id == a2.id  # type: ignore[union-attr]


def test_mark_only_moves_forward_and_is_idempotent() -> None:
    job = Job.objects.create(owner=user("alice"))
    assert not jobs.mark(job.id, Status.SUCCEEDED)  # queued -> succeeded skips running
    assert jobs.mark(job.id, Status.RUNNING)
    assert jobs.mark(job.id, Status.SUCCEEDED)
    assert not jobs.mark(job.id, Status.SUCCEEDED)  # repeated: no change
    assert not jobs.mark(job.id, Status.RUNNING)  # late, out-of-order event: ignored
    job.refresh_from_db()
    assert job.status == Status.SUCCEEDED and job.finished_at is not None


def test_failure_message_recorded() -> None:
    job = Job.objects.create(owner=user("alice"))
    assert jobs.mark(job.id, Status.FAILED, "Input error: bad FASTA")
    job.refresh_from_db()
    assert job.error == "Input error: bad FASTA"


def test_heartbeat_and_stale_sweep() -> None:
    alice = user("alice")
    live, dead = Job.objects.create(owner=alice), Job.objects.create(owner=alice)
    for job in (live, dead):
        assert jobs.mark(job.id, Status.RUNNING)
    assert jobs.heartbeat(live.id, {"done": 3, "total": 12})
    Job.objects.filter(id=dead.id).update(heartbeat_at=timezone.now() - timedelta(hours=1))
    assert jobs.sweep_stale(timedelta(minutes=10)) == 1
    live.refresh_from_db()
    dead.refresh_from_db()
    assert live.status == Status.RUNNING and live.progress == {"done": 3, "total": 12}
    assert dead.status == Status.FAILED and "stopped unexpectedly" in dead.error
    assert not jobs.heartbeat(dead.id)  # finished jobs take no heartbeats


def test_purge_retries_failed_deletes(monkeypatch: pytest.MonkeyPatch) -> None:
    from phylobot import storage

    job = Job.objects.create(owner=user("alice"))
    assert jobs.mark(job.id, Status.FAILED)

    def broken(prefix: str) -> int:
        raise OSError("storage down")

    monkeypatch.setattr(storage, "delete_prefix", broken)
    jobs.delete([job])
    job.refresh_from_db()
    assert job.deleted_at is not None and job.purged_at is None  # recorded; files remain
    monkeypatch.setattr(storage, "delete_prefix", lambda prefix: 0)
    assert jobs.purge_deleted() == 1  # the worker's sweep finishes the job
    job.refresh_from_db()
    assert job.purged_at is not None


def test_storage_refuses_broad_deletes() -> None:
    from phylobot import storage

    for prefix in ("", "jobs/", "jobs", "/"):
        with pytest.raises(ValueError):
            storage.delete_prefix(prefix)


def test_worker_removes_local_folders_of_deleted_jobs(tmp_path: object) -> None:
    from pathlib import Path

    from phylobot.management.commands.run_worker import Command

    root = Path(str(tmp_path))
    kept = Job.objects.create(owner=user("alice"))
    gone = Job.objects.create(owner=kept.owner)
    Job.objects.filter(id=gone.id).update(deleted_at=timezone.now())
    for job in (kept, gone):
        (root / "jobs" / str(job.id)).mkdir(parents=True)
    (root / "jobs" / "not-a-job").mkdir()
    Command.remove_local_folders(root, running=set())
    remaining = sorted(p.name for p in (root / "jobs").iterdir())
    assert remaining == sorted([str(kept.id), "not-a-job"])


def claimed(owner: User) -> Job:
    job = Job.objects.create(owner=owner)
    found = jobs.claim_next("w1")
    assert found is not None and found.id == job.id
    return found


def test_claim_opens_an_attempt() -> None:
    job = claimed(user("alice"))
    attempt = job.attempts.get()
    assert (attempt.number, attempt.worker, attempt.status) == (1, "w1", "")
    assert attempt.started_at == job.started_at
    assert jobs.mark(job.id, Status.FAILED, "boom", attempt=1)
    attempt.refresh_from_db()
    assert attempt.status == Status.FAILED and attempt.error == "boom" and attempt.finished_at


def test_resume_starts_the_next_attempt() -> None:
    alice = user("alice")
    job = claimed(alice)
    assert not jobs.resume(job)  # running: nothing to resume
    assert jobs.mark(job.id, Status.FAILED, jobs.STOPPED, attempt=1)
    job.refresh_from_db()
    assert jobs.resumable(job) and jobs.resume(job)
    assert (job.status, job.attempt, job.error, job.finished_at) == (Status.QUEUED, 2, "", None)
    assert not jobs.resume(job)  # already queued: a second click does nothing
    again = jobs.claim_next("w2")
    assert again is not None and again.attempt == 2
    assert [(a.number, a.worker) for a in job.attempts.all()] == [(1, "w1"), (2, "w2")]
    assert job.attempts.get(number=1).status == Status.FAILED  # the first stays as it ended


def test_superseded_runner_changes_nothing() -> None:
    job = claimed(user("alice"))
    Job.objects.filter(id=job.id).update(heartbeat_at=timezone.now() - timedelta(hours=1))
    assert jobs.sweep_stale(timedelta(minutes=10)) == 1  # its machine slept, say
    job.refresh_from_db()
    assert jobs.resume(job) and jobs.claim_next("w2") is not None
    # The first runner wakes up: its heartbeat and its outcome are refused.
    assert not jobs.heartbeat(job.id, {"done": 1}, attempt=1)
    assert not jobs.mark(job.id, Status.SUCCEEDED, attempt=1)
    assert jobs.heartbeat(job.id, {"done": 2}, attempt=2)
    job.refresh_from_db()
    assert job.status == Status.RUNNING and job.progress == {"done": 2}


def test_sweep_spares_the_sweepers_own_jobs() -> None:
    job = claimed(user("alice"))
    Job.objects.filter(id=job.id).update(heartbeat_at=timezone.now() - timedelta(hours=1))
    assert jobs.sweep_stale(timedelta(minutes=10), alive={str(job.id)}) == 0
    assert jobs.sweep_stale(timedelta(minutes=10)) == 1


def test_not_resumable() -> None:
    alice = user("alice")
    bad_input = claimed(alice)
    assert jobs.mark(bad_input.id, Status.FAILED, "Input error: bad FASTA")
    deleted = claimed(alice)
    assert jobs.mark(deleted.id, Status.FAILED, jobs.STOPPED)
    Job.objects.filter(id=deleted.id).update(deleted_at=timezone.now())
    done = claimed(alice)
    assert jobs.mark(done.id, Status.SUCCEEDED)
    for job in (bad_input, deleted, done):
        job.refresh_from_db()
        assert not jobs.resumable(job) and not jobs.resume(job)


def test_usage_bills_failed_steps_at_the_rate_of_the_day() -> None:
    job = claimed(user("alice"))
    SiteSettings.objects.create(pk=1, failed_task_rate_percent=50)
    jobs.record_usage(job.id, 1, completed=100.0, failed=40.0, run=5, reused=0)
    assert jobs.mark(job.id, Status.FAILED, jobs.STOPPED, attempt=1)
    job.refresh_from_db()
    assert job.cpu_seconds == 140 and job.billable_cpu_seconds == 120
    SiteSettings.objects.filter(pk=1).update(failed_task_rate_percent=0)
    assert jobs.resume(job) and jobs.claim_next("w2") is not None
    jobs.record_usage(job.id, 2, completed=30.0, failed=10.0, run=2, reused=5)
    job.refresh_from_db()
    first, second = job.attempts.all()
    assert (first.failed_rate_percent, first.billable_cpu_seconds) == (50, 120)  # kept
    assert (second.failed_rate_percent, second.billable_cpu_seconds) == (0, 30)
    assert job.cpu_seconds == 180 and job.billable_cpu_seconds == 150


def test_stop_cancels_a_live_job_and_cannot_be_resumed() -> None:
    alice = user("alice")
    queued = Job.objects.create(owner=alice)
    assert jobs.stop(queued)  # before a worker takes it
    queued.refresh_from_db()
    assert (queued.status, queued.error) == (Status.CANCELLED, jobs.STOPPED_BY_OWNER)
    assert queued.finished_at is not None
    # Cancelled is terminal: only a failed job resumes, so a stopped one never continues.
    assert not jobs.resumable(queued) and not jobs.resume(queued)

    running = claimed(alice)
    assert jobs.stop(running)
    running.refresh_from_db()
    assert running.status == Status.CANCELLED
    assert not jobs.stop(running)  # a second click changes nothing
    assert running.deleted_at is None  # stopping keeps the files; only delete removes them


def test_stop_does_nothing_to_a_finished_job() -> None:
    alice = user("alice")
    done = claimed(alice)
    assert jobs.mark(done.id, Status.SUCCEEDED)
    assert not jobs.stop(done)
    done.refresh_from_db()
    assert done.status == Status.SUCCEEDED and done.error == ""
