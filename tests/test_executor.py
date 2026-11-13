import re
import threading
from collections.abc import Callable
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

import pytest
from django.conf import settings
from django.contrib.auth.models import User
from django.test import override_settings

from phylobot import executors, jobs, manifest, storage
from phylobot.executors import (
    GENERIC_FAILURE,
    LocalDockerExecutor,
    build_params,
    expected_tasks,
    progress_from_trace,
    read_trace,
    run_job,
    staged_inputs,
    submitted_tasks,
    user_error,
)
from phylobot.models import Job, JobAttempt, UserLimits

ROOT = Path(__file__).resolve().parent.parent
TRACE_HEADER = "task_id\thash\tname\tstatus\texit\trealtime\tworkdir\tcpus\n"
pytestmark = pytest.mark.django_db


def make_job(params: dict[str, Any] | None = None) -> Job:
    owner = User.objects.create_user("alice")
    job = Job.objects.create(owner=owner, params=params or {})
    assert jobs.mark(job.id, Job.Status.RUNNING)
    return job


def test_runner_image_matches_compose() -> None:
    compose = (ROOT / "docker-compose.yml").read_text()
    runner = re.search(r"runner:\n(?:.*\n)*?\s+image: (\S+)", compose)
    assert runner and runner.group(1) == settings.PHYLOBOT_NEXTFLOW_IMAGE


def test_only_allowlisted_params_reach_nextflow(tmp_path: Path) -> None:
    job = make_job(
        {"aligners": "mafft", "models": ["LG+G4"], "outdir": "/etc", "input": "/etc/passwd",
         "max_sequences": 10**9, "-resume": True}
    )  # fmt: skip
    UserLimits.objects.create(user=job.owner, max_sequences=42)
    (tmp_path / "sequences.fasta").write_text(">a\nMK\n")
    params = build_params(job, tmp_path, tmp_path / "out")
    assert params["input"] == str(tmp_path / "sequences.fasta")
    assert params["outdir"] == str(tmp_path / "out")
    assert params["max_sequences"] == 42
    assert "-resume" not in params and params["models"] == ["LG+G4"]


def test_kept_branch_lengths_reach_nextflow(tmp_path: Path) -> None:
    """The oracle with true branch lengths is a job param: dropped by the
    allowlist, the job would silently refit them."""
    job = make_job(
        {"tree_labels": ["True tree", "True tree"], "tree_branch_lengths": ["fit", "keep"]}
    )
    params = build_params(job, tmp_path, tmp_path / "out")
    assert params["tree_branch_lengths"] == ["fit", "keep"]


def test_staged_inputs_in_numeric_order(tmp_path: Path) -> None:
    for name in ("alignment10.fasta", "alignment2.fasta", "tree1.nwk", "notes.txt"):
        (tmp_path / name).write_text("x")
    params = staged_inputs(tmp_path)
    assert [Path(p).name for p in params["alignment"]] == ["alignment2.fasta", "alignment10.fasta"]
    assert [Path(p).name for p in params["tree"]] == ["tree1.nwk"]
    assert "input" not in params


@pytest.mark.parametrize(
    ("params", "tasks"),
    [
        # VALIDATE, MAFFT, EXPORT_ALIGNMENT, MODELFINDER, TREE, ANCESTRAL, EXPORT, SUMMARY,
        # ANCESTOR_MAP (with combinations: + 1 in each case below that has trees)
        ({}, 9),
        # 3 alignments (each exported, each through ModelFinder) x 2 trees = 6 combinations
        ({"alignment": ["a1", "a2"], "tree": ["t1"]}, 1 + 1 + 3 + 3 + 6 * 3 + 1 + 1),
        ({"models": ["LG+G4", "WAG"], "fixed_models": True}, 1 + 1 + 1 + 1 + 3 * 3 + 1 + 1),
        # RAxML-NG: per combination a search, then IQ-TREE's fit (and the usual steps)
        ({"tree_builders": ["iqtree", "raxml"]}, 1 + 1 + 1 + 1 + 1 + 2 * 3 + 1 + 1),
        ({"tree_builders": ["raxml"], "stop_after": "trees"}, 1 + 1 + 1 + 1 + 1 + 2 + 1 + 1),
        # fixed models only: no MODELFINDER
        ({"models": ["LG+G4"], "fixed_models": True, "modelfinder": False}, 8),
        ({"aligners": "none", "alignment": ["a1"], "tree_builders": "none", "tree": ["t"]}, 8),
        ({"gap_methods": ["parsimony"]}, 9),  # parsimony runs inside EXPORT
        ({"gap_methods": ["parsimony", "ml"]}, 10),  # ml adds GAPS_ML per combination
        ({"gap_methods": ["parsimony", "ml", "prank", "grasp"]}, 12),  # one step each
        ({"outgroup": ["s1"]}, 10),  # ROOT_TREE per combination
        # stopping early: alignments only (VALIDATE, 2 aligners, 3 exports, SUMMARY) ...
        (
            {
                "stop_after": "alignments",
                "aligners": ["mafft", "muscle"],
                "alignment_labels": ["x"],
            },
            1 + 2 + 3 + 1,
        ),
        # ... or trees without ancestors: a tree and EXPORT_TREE per combination
        (
            {"stop_after": "trees", "gap_methods": ["ml"], "outgroup": ["s1"]},
            1 + 1 + 1 + 1 + 2 + 1 + 1,
        ),
    ],
)
def test_expected_tasks(params: dict[str, Any], tasks: int) -> None:
    assert expected_tasks(params) == tasks


def trace_rows(*rows: tuple[str, str, str]) -> list[dict[str, str]]:
    return [{"name": n, "status": s, "workdir": w, "realtime": "1s"} for n, s, w in rows]


def test_progress_from_trace() -> None:
    rows = trace_rows(("VALIDATE (1)", "COMPLETED", ""), ("ALIGN_MAFFT (mafft)", "COMPLETED", ""))
    rows[0] |= {"hash": "ab/123456", "start": "2026-09-29 14:32:05.123", "complete": "-"}
    running = {"ab/123456": ("VALIDATE (1)", "t0"), "cd/654321": ("TREE_IQTREE (mafft)", "t1")}
    progress = progress_from_trace(rows, 6, running)
    assert (progress["done"], progress["total"], progress["last"]) == (2, 6, "ALIGN_MAFFT (mafft)")
    assert progress["tasks"][0]["start"] == "2026-09-29T14:32:05.123000+00:00"
    assert progress["tasks"][0]["complete"] == ""
    # started tasks not yet in the trace are running
    assert progress["tasks"][2:] == [
        {"name": "TREE_IQTREE (mafft)", "status": "RUNNING", "start": "t1"}
    ]


def test_submitted_tasks(tmp_path: Path) -> None:
    log = tmp_path / "pipeline.log"
    assert submitted_tasks(log) == {}
    log.write_text("N E X T F L O W\n[e0/fe27ec] Submitted process > PHYLOBOT:VALIDATE (1)\n"
                   "[88/91f962] Submitted process > ALIGN_MAFFT (mafft)\nother\n")  # fmt: skip
    assert submitted_tasks(log) == {"e0/fe27ec": "VALIDATE (1)", "88/91f962": "ALIGN_MAFFT (mafft)"}


def test_trace_names_lose_the_workflow(tmp_path: Path) -> None:
    # main.nf's entry workflow runs PHYLOBOT or BRANCH_SITE: task names carry it.
    trace = tmp_path / "trace.tsv"
    trace.write_text(TRACE_HEADER + "1\tab/cdef12\tBRANCH_SITE:VALIDATE (1)\tFAILED\t1\t1s\t"
                     f"{tmp_path}\t1\n")  # fmt: skip
    rows = read_trace(trace)
    assert rows[0]["name"] == "VALIDATE (1)"
    (tmp_path / ".command.err").write_text("Input error: tree 1: no foreground branch\n")
    assert user_error(rows) == "Input error: tree 1: no foreground branch"


def test_user_error_only_from_validate(tmp_path: Path) -> None:
    (tmp_path / ".command.err").write_text("noise\nInput error: need at least 4 sequences\n")
    assert user_error(trace_rows(("VALIDATE (1)", "FAILED", str(tmp_path)))) == (
        "Input error: need at least 4 sequences"
    )
    assert user_error(trace_rows(("TREE_IQTREE (mafft)", "FAILED", str(tmp_path)))) == (
        GENERIC_FAILURE
    )


FakeRun = Callable[[list[str], Path, Callable[[], bool]], int]


def fake_pipeline(exit_code: int, trace: str, validate_err: str = "") -> FakeRun:
    """Stands in for `docker run ... nextflow`: writes a trace and results, ticks once."""

    def run(command: list[str], log: Path, tick: Callable[[], bool]) -> int:
        workdir = log.parent
        task = workdir / "work" / "ab" / "cd"
        task.mkdir(parents=True)
        (task / ".command.err").write_text(validate_err)
        info = workdir / "results" / "pipeline_info"
        info.mkdir(parents=True)
        (info / "trace.tsv").write_text(TRACE_HEADER + trace.format(task=task))
        if exit_code == 0:
            (workdir / "results" / "summary.json").write_text("{}")
        log.write_text("N E X T F L O W\n")
        tick()
        return exit_code

    return run


@pytest.fixture
def fake_storage(monkeypatch: pytest.MonkeyPatch) -> dict[str, list[str]]:
    uploaded: dict[str, list[str]] = {}

    def download_prefix(prefix: str, dest: Path) -> list[Path]:
        (dest / "sequences.fasta").write_text(">a\nMK\n")
        return [dest / "sequences.fasta"]

    def upload_tree(directory: Path, prefix: str, alive: Any = None, every: float = 30.0) -> int:
        files = [p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file()]
        uploaded[prefix] = sorted(files)
        if alive is not None:
            uploaded["alive"] = [str(alive())]  # the heartbeat the worker passes
        return len(files)

    monkeypatch.setattr(storage, "download_prefix", download_prefix)
    monkeypatch.setattr(storage, "upload_tree", upload_tree)
    return uploaded


def test_successful_run(tmp_path: Path, fake_storage: dict[str, list[str]]) -> None:
    job = make_job({"aligners": "mafft"})
    trace = (
        "1\tab\tVALIDATE (1)\tCOMPLETED\t0\t1s\t{task}\t1\n"
        "2\tcd\tTREE_IQTREE (mafft/modelfinder)\tCOMPLETED\t0\t1m 2s\t{task}\t4\n"
    )
    JobAttempt.objects.create(job=job, number=1)  # as jobs.claim_next makes it
    with override_settings(PHYLOBOT_PROJECT_ROOT=str(tmp_path)):
        executor = LocalDockerExecutor(run=fake_pipeline(0, trace))
        command = executor.nextflow_command(job, tmp_path / "jobs" / str(job.id))
        run_job(job, executor, threading.Event())
    job.refresh_from_db()
    assert job.status == Job.Status.SUCCEEDED and job.error == ""
    assert job.progress["done"] == 2 and job.progress["total"] == 9
    assert job.cpu_seconds == 1 + 62 * 4
    files = fake_storage[f"jobs/{job.id}/results/"]
    assert {"summary.json", "logs/pipeline.log", "logs/params.json"} <= set(files)
    assert fake_storage["alive"] == ["True"]  # it heartbeats while uploading
    # every uploaded file in the manifest, which pages read instead of listing storage
    assert manifest.paths(job) == set(files) and job.results_indexed_at is not None
    assert manifest.sizes(job)["summary.json"] == 2
    assert not (tmp_path / "jobs" / str(job.id)).exists()  # cleaned on success
    assert f"{tmp_path}:{tmp_path}" in command and command[-6:-4] == ["-params-file", command[-5]]
    # the attempt's share of the machine, written before launch (phylobot.machine)
    workdir = tmp_path / "jobs" / str(job.id)
    assert command[-2:] == ["-c", str(workdir / "resources.config")]
    attempt = job.attempts.get()
    assert (attempt.cpus, attempt.memory_bytes) == (8 // 1, 16 * 1024**3 - 2 * 1024**3 - 1024**3)


def test_nextflow_gets_the_cpu_setting(tmp_path: Path) -> None:
    """nextflow.config turns PHYLOBOT_DOCKER_CPUSET into each task's --cpuset-cpus."""
    job = make_job()
    for value in ("0-6", ""):
        with override_settings(PHYLOBOT_PROJECT_ROOT=str(tmp_path), PHYLOBOT_DOCKER_CPUSET=value):
            command = LocalDockerExecutor().nextflow_command(job, tmp_path)
        assert f"PHYLOBOT_DOCKER_CPUSET={value}" in command


def test_input_error_reaches_user(tmp_path: Path, fake_storage: dict[str, list[str]]) -> None:
    job = make_job()
    trace = "1\tab\tVALIDATE (1)\tFAILED\t1\t1s\t{task}\t1\n"
    err = "Input error: need at least 4 sequences, found 2\n"
    with override_settings(PHYLOBOT_PROJECT_ROOT=str(tmp_path)):
        run_job(job, LocalDockerExecutor(run=fake_pipeline(1, trace, err)), threading.Event())
    job.refresh_from_db()
    assert job.status == Job.Status.FAILED
    assert job.error == "Input error: need at least 4 sequences, found 2"
    assert "logs/tasks/VALIDATE_1.err" in fake_storage[f"jobs/{job.id}/results/"]
    assert (tmp_path / "jobs" / str(job.id) / "work").exists()  # kept for debugging


class RecordingExecutor(LocalDockerExecutor):
    """Records container stops instead of calling docker."""

    stopped: list[str]

    def stop(self, job: Job) -> None:
        self.stopped.append(str(job.id))

    def stop_earlier(self, job: Job) -> None:
        self.stopped.append(f"before {job.id}-{job.attempt}")


def test_timeout_and_crash(tmp_path: Path, fake_storage: dict[str, list[str]]) -> None:
    def too_slow(command: list[str], log: Path, tick: Callable[[], bool]) -> int:
        raise TimeoutError

    def crash(command: list[str], log: Path, tick: Callable[[], bool]) -> int:
        raise RuntimeError("docker not found")

    with override_settings(PHYLOBOT_PROJECT_ROOT=str(tmp_path)):
        slow = make_job()
        executor = RecordingExecutor(run=too_slow)
        executor.stopped = []
        run_job(slow, executor, threading.Event())
        slow.refresh_from_db()
        assert slow.status == Job.Status.FAILED and "time limit" in slow.error
        assert executor.stopped == [str(slow.id)]  # the container is stopped, not orphaned

        broken = Job.objects.create(owner=slow.owner)
        assert jobs.mark(broken.id, Job.Status.RUNNING)
        run_job(broken, LocalDockerExecutor(run=crash), threading.Event())
        broken.refresh_from_db()
        assert broken.status == Job.Status.FAILED and broken.error == GENERIC_FAILURE


def record(ticks: list[int]) -> bool:
    ticks.append(1)
    return True


def test_run_command_heartbeats_and_times_out(tmp_path: Path) -> None:
    ticks: list[int] = []
    with override_settings(PHYLOBOT_HEARTBEAT_SECONDS=0.2, PHYLOBOT_JOB_TIME_LIMIT=timedelta(0)):
        with pytest.raises(TimeoutError):
            executors.run_command(["sleep", "5"], tmp_path / "log", lambda: record(ticks))
    assert ticks
    with override_settings(PHYLOBOT_HEARTBEAT_SECONDS=0.2):
        with pytest.raises(executors.Cancelled):
            executors.run_command(["sleep", "5"], tmp_path / "log", lambda: False)
    with override_settings(PHYLOBOT_HEARTBEAT_SECONDS=5):
        command = ["sh", "-c", "echo hi; exit 3"]
        assert executors.run_command(command, tmp_path / "log", lambda: True) == 3
    assert (tmp_path / "log").read_text() == "hi\n"


def test_deleted_while_running_is_stopped_and_purged(
    tmp_path: Path, fake_storage: dict[str, list[str]], monkeypatch: pytest.MonkeyPatch
) -> None:
    job = make_job()
    purged: list[str] = []

    def delete_prefix(prefix: str) -> int:
        purged.append(prefix)
        return 0

    monkeypatch.setattr(storage, "delete_prefix", delete_prefix)

    def deleted_midway(command: list[str], log: Path, tick: Callable[[], bool]) -> int:
        assert tick()  # running: heartbeat accepted
        jobs.delete([Job.objects.get(id=job.id)])  # the owner deletes it now
        if not tick():
            raise executors.Cancelled
        return 0

    with override_settings(PHYLOBOT_PROJECT_ROOT=str(tmp_path)):
        executor = RecordingExecutor(run=deleted_midway)
        executor.stopped = []
        run_job(job, executor, threading.Event())
    job.refresh_from_db()
    assert job.status == Job.Status.CANCELLED and job.deleted_at is not None
    assert executor.stopped == [str(job.id)]
    assert f"jobs/{job.id}/results/" not in fake_storage  # no results uploaded
    # purged when deleted, and again by the worker after stopping (catches in-flight uploads)
    assert purged == [f"jobs/{job.id}/", f"work-archives/{job.id}/"] * 2
    assert job.purged_at is not None
    assert not (tmp_path / "jobs" / str(job.id)).exists()


def test_species_files_become_reconciliation_params(tmp_path: Path) -> None:
    for name in ("sequences.fasta", "species.tsv", "species_tree.nwk"):
        (tmp_path / name).write_text("x")
    params = staged_inputs(tmp_path)
    assert params["species_map"] == str(tmp_path / "species.tsv")
    assert params["species_tree"] == str(tmp_path / "species_tree.nwk")
    # The worker counts tasks from these params, which lack the web's `reconcile` flag.
    assert expected_tasks(params) == expected_tasks({"reconcile": True})


def test_usage_splits_completed_failed_and_reused() -> None:
    rows = [
        {"status": "COMPLETED", "realtime": "10s", "cpus": "2"},
        {"status": "CACHED", "realtime": "1h", "cpus": "1"},
        {"status": "FAILED", "realtime": "5s", "cpus": "1"},
        {"status": "ABORTED", "realtime": "", "cpus": "1"},
    ]
    assert executors.usage(rows) == executors.Usage(completed=20, failed=5, run=3, reused=1)
    assert executors.cpu_seconds(rows) == 25  # reused steps were counted where they ran
    # A step stopped part-way has no runtime: estimated from its start to the run's end.
    rows[3] |= {"start": "2026-10-19 01:31:37.000", "cpus": "2"}
    ended = datetime(2026, 10, 19, 1, 32, 7, tzinfo=UTC)
    assert executors.usage(rows, ended).failed == 5 + 30 * 2


def resumable_pipeline(runs: list[list[str]], traces: list[str]) -> FakeRun:
    """Like fake_pipeline, but runs again in the same folder: run n writes traces[n]."""

    def run(command: list[str], log: Path, tick: Callable[[], bool]) -> int:
        workdir = log.parent
        runs.append([*command, (workdir / "params.json").read_text()])
        task = workdir / "work" / "ab" / "cd"
        task.mkdir(parents=True, exist_ok=True)
        info = workdir / "results" / "pipeline_info"
        info.mkdir(parents=True, exist_ok=True)
        trace = traces[len(runs) - 1]
        (info / "trace.tsv").write_text(TRACE_HEADER + trace.format(task=task))
        log.write_text(f"N E X T F L O W run {len(runs)}\n")
        assert tick()
        return 0 if "FAILED" not in trace else 1

    return run


def test_resume_reuses_the_folder_and_bills_only_new_steps(
    tmp_path: Path, fake_storage: dict[str, list[str]], monkeypatch: pytest.MonkeyPatch
) -> None:
    downloads: list[str] = []
    real_download = storage.download_prefix

    def counting(prefix: str, dest: Path) -> list[Path]:
        downloads.append(prefix)
        return real_download(prefix, dest)

    monkeypatch.setattr(storage, "download_prefix", counting)
    owner = User.objects.create_user("alice")
    Job.objects.create(owner=owner, params={"aligners": "mafft"})
    job = jobs.claim_next("w1")
    assert job is not None
    runs: list[list[str]] = []
    traces = [
        "1\tab\tVALIDATE (1)\tCOMPLETED\t0\t1s\t{task}\t1\n"
        "2\tcd\tTREE_IQTREE (mafft/modelfinder)\tFAILED\t137\t1m\t{task}\t4\n",
        "1\tab\tVALIDATE (1)\tCACHED\t0\t1s\t{task}\t1\n"
        "2\tef\tTREE_IQTREE (mafft/modelfinder)\tCOMPLETED\t0\t2m\t{task}\t4\n",
    ]
    folder = tmp_path / "jobs" / str(job.id)
    with override_settings(PHYLOBOT_PROJECT_ROOT=str(tmp_path)):
        executor = RecordingExecutor(run=resumable_pipeline(runs, traces))
        executor.stopped = []
        run_job(job, executor, threading.Event())
        job.refresh_from_db()
        assert job.status == Job.Status.FAILED and folder.is_dir() and job.resume_cache
        (folder / ".nextflow").mkdir()  # what Nextflow leaves in the folder it ran from
        assert jobs.resume(job)
        again = jobs.claim_next("w1")
        assert again is not None
        run_job(again, executor, threading.Event())
    job.refresh_from_db()
    assert job.status == Job.Status.SUCCEEDED and job.attempt == 2
    assert executor.stopped == [f"before {job.id}-2"]  # no earlier container left running
    assert downloads == [f"jobs/{job.id}/inputs/"]  # inputs fetched once, kept for the resume
    (*first, params1), (*second, params2) = runs
    assert params1 == params2  # not rebuilt: a changed limit would invalidate every step
    assert "-resume" not in first and second[-1] == "-resume"
    assert second[second.index("-w") + 1] == str(folder)  # its own resume cache
    assert "--name" in second and second[second.index("--name") + 1].endswith(f"{job.id}-2")
    one, two = job.attempts.all()
    assert (one.completed_cpu_seconds, one.failed_cpu_seconds, one.status) == (1, 240, "failed")
    assert (two.completed_cpu_seconds, two.tasks_run, two.tasks_reused) == (480, 1, 1)
    assert job.cpu_seconds == 1 + 240 + 480 and job.billable_cpu_seconds == job.cpu_seconds
    files = fake_storage[f"jobs/{job.id}/results/"]
    assert {"logs/attempts/1/pipeline.log", "logs/attempts/1/trace.tsv"} <= set(files)
    assert "logs/pipeline.log" in files and not folder.exists()  # cleaned on success


def test_superseded_attempt_leaves_the_folder(
    tmp_path: Path, fake_storage: dict[str, list[str]]
) -> None:
    job = make_job()
    folder = tmp_path / "jobs" / str(job.id)

    def overtaken(command: list[str], log: Path, tick: Callable[[], bool]) -> int:
        assert tick()
        Job.objects.filter(id=job.id).update(attempt=2)  # declared stopped and resumed
        if not tick():
            raise executors.Cancelled
        return 0

    with override_settings(PHYLOBOT_PROJECT_ROOT=str(tmp_path)):
        executor = RecordingExecutor(run=overtaken)
        executor.stopped = []
        run_job(job, executor, threading.Event())
    job.refresh_from_db()
    assert job.status == Job.Status.RUNNING  # the outcome of attempt 1 is refused
    assert folder.is_dir()  # attempt 2 is using it


def failed_with_folder(root: Path, ended: datetime, size: int = 0) -> tuple[Job, Path]:
    job = Job.objects.create(
        owner=User.objects.get_or_create(username="alice")[0], status=Job.Status.FAILED,
        finished_at=ended, resume_cache=True,
    )  # fmt: skip
    folder = root / "jobs" / str(job.id)
    (folder / "work" / "ab").mkdir(parents=True)
    (folder / "work" / "ab" / "out.txt").write_bytes(b"x" * size)
    return job, folder


def test_work_folders_kept_until_the_settings_period_ends(tmp_path: Path) -> None:
    from phylobot.models import SiteSettings

    now = datetime.now(UTC)
    old, old_folder = failed_with_folder(tmp_path, now - timedelta(days=10))
    new, new_folder = failed_with_folder(tmp_path, now - timedelta(days=2))
    busy, busy_folder = failed_with_folder(tmp_path, now - timedelta(days=30))
    assert executors.expire_work_folders(tmp_path, set(), now) == 0  # unset: kept forever
    SiteSettings.objects.create(pk=1, keep_failed_work_days=7)
    assert executors.expire_work_folders(tmp_path, {str(busy.id)}, now) == 1
    assert not old_folder.exists() and new_folder.is_dir() and busy_folder.is_dir()
    old.refresh_from_db()
    assert old.work_expired_at and not old.resume_cache
    assert executors.expire_work_folders(tmp_path, {str(busy.id)}, now) == 0
    # A resume starts afresh, with a cache of its own again.
    assert jobs.resume(old)
    old.refresh_from_db()
    assert old.work_expired_at is None and old.attempt == 2


def test_work_folders_measured_once_without_following_links(tmp_path: Path) -> None:
    job, folder = failed_with_folder(tmp_path, datetime.now(UTC), size=1000)
    (folder / "work" / "cd").mkdir()
    (folder / "work" / "cd" / "in.txt").symlink_to(folder / "work" / "ab" / "out.txt")
    assert executors.measure_work_folders(tmp_path, set()) == 1
    job.refresh_from_db()
    assert job.work_bytes is not None and 1000 <= job.work_bytes < 1200  # once, plus the link
    assert executors.measure_work_folders(tmp_path, set()) == 0
    assert jobs.resume(job)  # the folder is the next attempt's now: measured again later
    job.refresh_from_db()
    assert job.work_bytes is None


def test_record_unrecorded_reads_a_stopped_attempts_trace(tmp_path: Path) -> None:
    owner = User.objects.create_user("alice")
    Job.objects.create(owner=owner)
    job = jobs.claim_next("w1")
    assert job is not None
    info = tmp_path / "jobs" / str(job.id) / "results" / "pipeline_info"
    info.mkdir(parents=True)
    (info / "trace.tsv").write_text(
        TRACE_HEADER + "1\tab\tVALIDATE (1)\tCOMPLETED\t0\t30s\t/x\t2\n"
    )
    assert jobs.mark(job.id, Job.Status.FAILED, jobs.STOPPED)  # swept: nothing recorded
    assert executors.record_unrecorded(tmp_path) == 1
    assert executors.record_unrecorded(tmp_path) == 0
    job.refresh_from_db()
    assert job.cpu_seconds == 60 and job.attempts.get().completed_cpu_seconds == 60


def test_orphans_are_only_the_jobs_own_tasks() -> None:
    folder = Path("/p/jobs/abc")
    inspected = (
        "aaa /p/jobs/abc/work/12/3456\n"
        "bbb /p/jobs/abcdef/work/12/3456\n"  # another job whose id starts the same
        "ccc /p/results/golden/ime2/work/ab/cd\n"
        "ddd \n"
    )
    assert executors.orphans(inspected, folder) == ["aaa"]


def test_upload_heartbeats_and_stops_when_the_job_is_gone(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A 524k-file upload took 1.5 h without a heartbeat. It now calls `alive` as it goes
    and stops when the job is no longer this attempt's."""
    for n in range(5):
        (tmp_path / f"{n}.tsv").write_text("x")
    put: list[str] = []
    monkeypatch.setattr(storage, "put_file", lambda key, path: put.append(key))
    beats: list[int] = []

    def alive() -> bool:
        beats.append(len(put))
        return len(put) < 3

    with pytest.raises(storage.UploadStopped):
        storage.upload_tree(tmp_path, "jobs/x/results/", alive=alive, every=0)
    assert beats == [1, 2, 3] and len(put) == 3


def test_a_job_deleted_during_its_upload_is_cancelled(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    job = make_job({"aligners": "mafft"})
    JobAttempt.objects.create(job=job, number=1)

    def upload_tree(directory: Path, prefix: str, alive: Any = None, every: float = 30.0) -> int:
        Job.objects.filter(id=job.id).update(attempt=2)  # resumed elsewhere meanwhile
        assert alive is not None and not alive()
        raise storage.UploadStopped

    monkeypatch.setattr(storage, "download_prefix", lambda prefix, dest: [])
    monkeypatch.setattr(storage, "upload_tree", upload_tree)
    trace = "1\tab\tVALIDATE (1)\tCOMPLETED\t0\t1s\t{task}\t1\n"
    with override_settings(PHYLOBOT_PROJECT_ROOT=str(tmp_path)):
        outcome = LocalDockerExecutor(run=fake_pipeline(0, trace)).run(job)
    assert outcome.status == Job.Status.CANCELLED
    job.refresh_from_db()
    assert job.results_indexed_at is None  # left for maintenance, never half-indexed
