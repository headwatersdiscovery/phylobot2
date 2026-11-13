"""Rough time-remaining estimates from how long each kind of step took in past jobs.

Kinds are process names (TREE_IQTREE, ANCESTRAL_IQTREE, ...). The estimate sums the
median past duration of each remaining step, i.e. as if steps ran one at a time, so
it errs long when combinations run in parallel. It does not yet scale with sequence
count or length.
"""

from statistics import median

from phylobot.durations import seconds
from phylobot.models import Job
from phylobot.stages import expected_kinds, kind

HISTORY = 50  # most recent successful jobs to learn from


def typical_durations(history: int = HISTORY) -> dict[str, float]:
    samples: dict[str, list[float]] = {}
    recent = Job.objects.filter(status=Job.Status.SUCCEEDED).order_by("-finished_at")[:history]
    for progress in recent.values_list("progress", flat=True):
        for task in (progress or {}).get("tasks", []):
            if (duration := seconds(task.get("realtime", ""))) is not None:
                samples.setdefault(kind(task["name"]), []).append(duration)
    return {name: median(values) for name, values in samples.items()}


def remaining_seconds(job: Job, typical: dict[str, float] | None = None) -> float | None:
    """Estimated seconds left for a queued or running job; None if there is no history
    for some remaining kind of step."""
    typical = typical_durations() if typical is None else typical
    remaining = expected_kinds(job.params)
    for task in job.progress.get("tasks", []):
        if task.get("status") in {"COMPLETED", "CACHED"}:
            remaining[kind(task["name"])] -= 1
    total = 0.0
    for name, count in (+remaining).items():
        if name not in typical:
            return None
        total += count * typical[name]
    return total


def describe(total_seconds: float | None, *, short: bool = False) -> str:
    """Human-readable remaining time; `short` for table cells (no long explanation)."""
    if total_seconds is None:
        return "estimating…" if short else "Estimate available after a few jobs have completed"
    if total_seconds < 60:
        return "less than a minute"
    minutes = round(total_seconds / 60)
    if minutes < 90:
        return f"about {minutes} min"
    return f"about {total_seconds / 3600:.1f} h"
