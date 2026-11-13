"""Work records: each finished job's estimate beside the CPU time it used, stage by stage,
for the admin (Settings -> Work records). What `workload.calibration` learns from, and
what a refit of the step costs would start from.

Read from Job rows only (the trace's rows are kept in `progress`), never storage. A
stage's CPU time counts each step once, as the estimate does: steps that ran to the end
(COMPLETED) or were reused from an earlier attempt (CACHED). Failed steps and retries
are the job's total's, not any stage's.
"""

import csv
import io
from collections import defaultdict
from dataclasses import dataclass, field
from statistics import median

from phylobot import workload
from phylobot.durations import seconds
from phylobot.models import Job, JobAttempt
from phylobot.usage import stage_of

DONE = {"COMPLETED", "CACHED"}
DEFAULT_LIMIT = 100
MAX_LIMIT = 2000


@dataclass
class Stage:
    key: str
    label: str
    tasks: int = 0
    estimated_units: float | None = None  # None: the job's input size was not recorded
    cpu_seconds: float = 0.0
    peak_rss_bytes: int | None = None  # largest of its steps, where the trace has it

    def ratio(self, seconds_per_unit: float) -> float | None:
        if not self.estimated_units or not self.cpu_seconds:
            return None
        return self.cpu_seconds / (self.estimated_units * seconds_per_unit)


@dataclass
class Record:
    job: Job
    size: workload.Size | None
    cpus: int | None  # the share its last attempt was given (None before it was recorded)
    memory_bytes: int | None
    stages: list[Stage] = field(default_factory=list)
    counted: bool = False  # among the jobs this server's calibration learns from

    def ratio(self, seconds_per_unit: float) -> float | None:
        if not self.job.work_units or not self.job.cpu_seconds:
            return None
        return self.job.cpu_seconds / (self.job.work_units * seconds_per_unit)


def memory(text: str) -> int | None:
    """A trace's "1.2 GB" / "830.5 MB" / "12 KB" as bytes."""
    units = {"B": 1, "KB": 1024, "MB": 1024**2, "GB": 1024**3, "TB": 1024**4}
    number, _, unit = text.strip().partition(" ")
    try:
        return int(float(number) * units[unit])
    except KeyError, ValueError:
        return None


def stages(job: Job) -> list[Stage]:
    """The job's stages in pipeline order: estimated (at its recorded size) and used."""
    by_key = {k: Stage(k, label) for k, label in workload.GROUP_LABELS.items()}
    if (size := workload.recorded_size(job.input_size)) is not None:
        try:
            planned = workload.tasks(job.params or {}, size)
        except Exception:  # an older job's params that today's estimate no longer reads
            planned = []
        for kind, cost in planned:
            stage = by_key[workload.group(kind)]
            stage.estimated_units = (stage.estimated_units or 0) + cost / workload.SECONDS_PER_UNIT
    for task in (job.progress or {}).get("tasks", []):
        if task.get("status") not in DONE:
            continue
        stage = by_key[stage_of(task.get("name", ""))]
        stage.tasks += 1
        if (runtime := seconds(task.get("realtime", ""))) is not None:
            stage.cpu_seconds += runtime * int(task.get("cpus") or 1)
        if (peak := memory(task.get("peak_rss") or "")) is not None:
            stage.peak_rss_bytes = max(stage.peak_rss_bytes or 0, peak)
    return [s for s in by_key.values() if s.tasks or s.estimated_units]


def records(limit: int = DEFAULT_LIMIT) -> list[Record]:
    """The most recent finished jobs with recorded CPU time, everyone's, newest first."""
    finished = list(
        Job.objects.filter(cpu_seconds__gt=0, imported_at__isnull=True)
        .exclude(status__in=[Job.Status.QUEUED, Job.Status.RUNNING])
        .select_related("owner")
        .order_by("-finished_at")[:limit]
    )
    counted = set(
        Job.objects.filter(status=Job.Status.SUCCEEDED, cpu_seconds__gt=0)
        .order_by("-finished_at")
        .values_list("id", flat=True)[: workload.HISTORY]
    )
    shares = {
        a.job_id: (a.cpus, a.memory_bytes)
        for a in JobAttempt.objects.filter(job__in=finished).order_by("number")
    }  # the last attempt's wins
    return [
        Record(
            job=job,
            size=workload.recorded_size(job.input_size),
            cpus=shares.get(job.id, (None, None))[0],
            memory_bytes=shares.get(job.id, (None, None))[1],
            stages=stages(job),
            counted=job.id in counted,
        )
        for job in finished
    ]


@dataclass(frozen=True)
class StageSummary:
    label: str
    jobs: int  # with both an estimate and CPU time for this stage
    median_ratio: float | None
    cpu_seconds: float


def summary(rows: list[Record], seconds_per_unit: float) -> list[StageSummary]:
    """Per stage, over successful jobs: how far its estimates are from what it used."""
    ratios: dict[str, list[float]] = defaultdict(list)
    used: dict[str, float] = defaultdict(float)
    labels: dict[str, str] = {}
    for row in rows:
        if row.job.status != Job.Status.SUCCEEDED:
            continue
        for stage in row.stages:
            labels[stage.key] = stage.label
            used[stage.key] += stage.cpu_seconds
            if (ratio := stage.ratio(seconds_per_unit)) is not None:
                ratios[stage.key].append(ratio)
    return [
        StageSummary(
            labels[k],
            len(ratios[k]),
            median(ratios[k]) if ratios[k] else None,
            used[k],
        )
        for k in workload.GROUP_LABELS
        if k in labels
    ]


CSV_COLUMNS = [
    "job", "name", "owner", "status", "finished", "sequences", "length", "cpus_given",
    "memory_given_gb", "job_work_units", "job_cpu_seconds", "counted_in_calibration",
    "stage", "stage_tasks", "stage_estimated_units", "stage_cpu_seconds",
    "stage_peak_rss_mb",
]  # fmt: skip


def to_csv(rows: list[Record]) -> str:
    """One line per job and stage: the table, for refitting elsewhere."""
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(CSV_COLUMNS)
    for row in rows:
        job = row.job
        common = [
            job.id, job.name, job.owner.get_username(), job.status,
            job.finished_at.isoformat() if job.finished_at else "",
            row.size.sequences if row.size else "",
            f"{row.size.length:.1f}" if row.size else "",
            row.cpus or "",
            f"{row.memory_bytes / 1024**3:.2f}" if row.memory_bytes else "",
            f"{job.work_units:.3f}" if job.work_units else "",
            f"{job.cpu_seconds:.1f}" if job.cpu_seconds else "",
            int(row.counted),
        ]  # fmt: skip
        for stage in row.stages or [Stage("", "")]:
            writer.writerow(common + [
                stage.key, stage.tasks,
                f"{stage.estimated_units:.3f}" if stage.estimated_units is not None else "",
                f"{stage.cpu_seconds:.1f}",
                f"{stage.peak_rss_bytes / 1024**2:.1f}" if stage.peak_rss_bytes else "",
            ])  # fmt: skip
    return out.getvalue()
