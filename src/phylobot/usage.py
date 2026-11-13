"""My usage (`phylobot.usage_views`): what a user's jobs used, CPU time and storage, over
a chosen window, as figures, charts and tables.

Everything comes from fields the Job records already hold (`cpu_seconds`,
`storage_bytes`, the steps' run times in `progress`, `deleted_at`), never from listing
storage, so the page costs one query. Deleted jobs count by default: their compute was
used, and their storage was held until they were deleted. A job belongs to the window by
its submission date. Charts are geometry for server-drawn SVG (as `treeview`); their
colours are the validated categorical slots `--viz-1`.. in site.css, one per job kind.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta
from typing import Any

from django.http import QueryDict
from django.utils import timezone

from phylobot import job_types
from phylobot.durations import seconds
from phylobot.models import Job

# Each job counted once, by the furthest thing it did; slot = its colour (--viz-<slot>).
# In pipeline order, which is the order the chart and the table read in. Not every
# `job_types` key belongs here: "reconciled" is an attribute of an ancestors job rather
# than a further stage, so `kind_of` passes over it (see there).
KINDS = {
    "alignments": ("Alignments only", 1),
    "trees": ("Trees", 2),
    "ancestors": ("Ancestors", 3),
    "properties": ("Ancestral properties", 6),
    "uncertainty": ("Ancestral uncertainty", 7),
    # Never beside Simulations: slot 8 is too close to slot 5 in the dark skin (site.css).
    "esr": ("ESR", 8),
    "dnds": ("dN/dS", 4),
    "simulation": ("Simulations", 5),
}
KIND_HELP = (
    "Each job counts once, by the furthest it went: alignments only, up to trees, up to "
    "ancestors, measuring properties or uncertainty on a finished job's ancestors, an "
    "extant-sequence reconstruction (ESR), codeml dN/dS (branch-site tests too), or a "
    "simulation."
)
PRESETS = {
    "7d": ("Last 7 days", 7),
    "30d": ("Last 30 days", 30),
    "90d": ("Last 90 days", 90),
    "1y": ("Last 12 months", 365),
    "all": ("All time", None),
}
DEFAULT_PRESET = "90d"
FAILED = (Job.Status.FAILED, Job.Status.CANCELLED)
# Pipeline stages for "where the CPU time went", in pipeline order (workload's groups).
STAGES = {
    "setup": "Checks and exports",
    "align": "Alignment",
    "modelfinder": "ModelFinder",
    "trees": "Tree searches and fits",
    "ancestors": "Ancestral reconstruction",
    "properties": "Ancestral properties",
    "codeml": "codeml (dN/dS)",
    "simulation": "Simulation (AliSim)",
}


def kind_of(params: dict[str, Any]) -> str:
    """Which KIND counts this job: the furthest of its types that is one.

    `job_types.types_of` also returns types that are not stages a job can stop at --
    "reconciled" is a property of an ancestors job, not something past ancestors -- and
    taking its last entry blindly returned a key this page has no slot for, so the job
    was dropped from the breakdown without trace (every reconciled and every
    properties-only job: 14 of 148 locally when this was found).
    """
    keys = [t.key for t in job_types.types_of(params) if t.key in KINDS]
    return keys[-1] if keys else "ancestors"


def stage_of(task: str) -> str:
    """A trace task name ("PHYLOBOT:ALIGN_MAFFT (mafft)") -> its workload group."""
    from phylobot.workload import group  # workload imports executors; keep this module light

    name = task.rsplit(":", 1)[-1].split(" (", 1)[0].strip()
    return group(name)


# --- The filters -------------------------------------------------------------------


@dataclass(frozen=True)
class Window:
    start: date
    end: date  # inclusive
    preset: str  # a PRESETS key, or "custom"

    @property
    def days(self) -> int:
        return (self.end - self.start).days + 1

    @property
    def label(self) -> str:
        if self.preset in PRESETS:
            return PRESETS[self.preset][0].lower()
        return f"{self.start:%Y-%m-%d} to {self.end:%Y-%m-%d}"

    def previous(self) -> Window | None:
        """The window of the same length just before this one (for the change figure)."""
        if self.preset == "all":
            return None
        return Window(self.start - timedelta(days=self.days), self.start - timedelta(days=1), "")


def _date(value: str | None) -> date | None:
    try:
        return date.fromisoformat(value) if value else None
    except ValueError:
        return None


@dataclass
class Filters:
    window: Window
    kinds: list[str] = field(default_factory=list)  # empty = every kind
    include_deleted: bool = True

    @classmethod
    def parse(cls, query: QueryDict, today: date, first_day: date) -> Filters:
        preset = query.get("range", DEFAULT_PRESET)
        start, end = _date(query.get("from")), _date(query.get("to"))
        if preset == "custom" or (preset not in PRESETS and (start or end)):
            end = min(end or today, today)
            start = min(start or first_day, end)
            window = Window(start, end, "custom")
        else:
            preset = preset if preset in PRESETS else DEFAULT_PRESET
            days = PRESETS[preset][1]
            begin = min(first_day, today) if days is None else today - timedelta(days=days - 1)
            window = Window(begin, today, preset)
        return cls(
            window=window,
            kinds=[k for k in query.getlist("kind") if k in KINDS],
            include_deleted=query.get("deleted", "1") != "0",
        )

    def wants(self, job: Job) -> bool:
        return (not self.kinds or kind_of(job.params) in self.kinds) and (
            self.include_deleted or job.deleted_at is None
        )


# --- Formatting --------------------------------------------------------------------


def duration(total: float) -> str:
    """CPU time as "45 s", "3 min", "2 h 10 min" or "31 h" (no seconds past an hour)."""
    total = int(round(total))
    if total < 60:
        return f"{total} s"
    minutes = total // 60
    if minutes < 60:
        return f"{minutes} min"
    hours, minutes = divmod(minutes, 60)
    return f"{hours} h {minutes} min" if minutes and hours < 10 else f"{hours:,} h"


def size(total: float) -> str:
    """Bytes in the units My Jobs uses (Django's filesizeformat: 1024-based)."""
    for unit in ("bytes", "KB", "MB", "GB", "TB"):
        if total < 1024 or unit == "TB":
            return f"{total:.0f} {unit}" if unit == "bytes" else f"{total:.1f} {unit}"
        total /= 1024
    raise AssertionError


def percent(part: float, whole: float) -> str:
    if not whole:
        return "0%"
    share = 100 * part / whole
    return "<1%" if 0 < share < 1 else f"{share:.0f}%"


# --- Time buckets ------------------------------------------------------------------


@dataclass(frozen=True)
class Bucket:
    start: date
    end: date  # inclusive

    def label(self, unit: str) -> str:
        if unit == "day":
            return f"{self.start:%a %d %b %Y}"
        if unit == "week":
            return f"Week of {self.start:%d %b %Y}"
        return f"{self.start:%B %Y}"


def unit_for(window: Window) -> str:
    """Days up to two months, weeks up to two years, then months (about 60 points at most)."""
    return "day" if window.days <= 62 else "week" if window.days <= 730 else "month"


def buckets(window: Window, unit: str) -> list[Bucket]:
    out = []
    day = window.start
    if unit == "week":
        day -= timedelta(days=day.weekday())  # Mondays
    elif unit == "month":
        day = day.replace(day=1)
    while day <= window.end:
        if unit == "day":
            following = day + timedelta(days=1)
        elif unit == "week":
            following = day + timedelta(days=7)
        else:
            following = (day.replace(day=28) + timedelta(days=4)).replace(day=1)
        out.append(Bucket(max(day, window.start), min(following - timedelta(days=1), window.end)))
        day = following
    return out


def local_day(moment: datetime) -> date:
    return timezone.localtime(moment).date()


# --- Chart geometry ----------------------------------------------------------------

WIDTH = 560  # about a card's width at two cards a row, so text stays near its set size


def nice_ceiling(value: float) -> tuple[float, float]:
    """(axis top, tick step) for values up to `value`: 4 ticks at 1, 2, 2.5 or 5 x 10^n."""
    if value <= 0:
        return 1.0, 0.25
    raw = value / 4
    power = 10 ** math.floor(math.log10(raw))
    step = next(m * power for m in (1, 2, 2.5, 5, 10) if m * power >= raw)
    return step * math.ceil(value / step - 1e-9), step


@dataclass
class Point:
    x: float
    y: float
    label: str  # the bucket
    value: str  # formatted for the tooltip and table


@dataclass
class LineChart:
    points: list[Point]
    line: str  # SVG path
    area: str
    yticks: list[tuple[float, str]]
    xticks: list[tuple[float, str]]
    left: float
    right: float
    top: float
    bottom: float
    end_label: str  # the last value, labelled at the line's end
    width: int = WIDTH
    height: int = 230
    hit_width: float = 0.0  # each point's hover column

    @property
    def empty(self) -> bool:
        return not self.points

    @property
    def plot_height(self) -> float:
        return self.bottom - self.top

    @property
    def tick_x(self) -> float:
        return self.left - 8

    @property
    def xlabel_y(self) -> float:
        return self.bottom + 20

    @property
    def last(self) -> Point:
        return self.points[-1]

    @property
    def end_label_y(self) -> float:
        """Above every point under the end label (right-aligned, ~8 px per character), so
        the label never sits on the line."""
        reach = self.last.x - 9 * len(self.end_label) - 14
        return min(p.y for p in self.points if p.x >= reach) - 12

    @property
    def hits(self) -> list[tuple[Point, float]]:
        """Each point with the left edge of its hover column (columns tile the plot)."""
        half = self.hit_width / 2
        return [(p, max(self.left - half, p.x - half)) for p in self.points]


def axis_scale(top: float, kind: str) -> tuple[float, str]:
    """(divisor, unit) for an axis whose largest value is `top` (seconds or bytes)."""
    if kind == "cpu":
        if top < 120:
            return 1, "s"
        if top < 2 * 3600:
            return 60, "min"
        return 3600, "h"
    scale = 1.0
    for unit in ("bytes", "KB", "MB", "GB", "TB"):
        if top < 1024 * scale or unit == "TB":
            return scale, unit
        scale *= 1024
    raise AssertionError


def tick_text(value: float, unit: str) -> str:
    text = f"{value:,.0f}" if value == int(value) else f"{value:,.1f}"
    return f"{text} {unit}" if value else "0"


def line_chart(
    values: list[float], spans: list[Bucket], unit: str, kind: str, fmt: Any
) -> LineChart:
    left, right, top, bottom = 64.0, WIDTH - 20.0, 16.0, 230 - 30.0
    if not values:
        return LineChart([], "", "", [], [], left, right, top, bottom, "")
    divisor, axis_unit = axis_scale(max(values), kind)
    ceiling, step = nice_ceiling(max(values) / divisor)
    span = right - left
    n = len(values)
    xs = [left + span / 2] if n == 1 else [left + i * span / (n - 1) for i in range(n)]

    def y_of(v: float) -> float:
        return bottom - (v / divisor) / ceiling * (bottom - top)

    points = [
        Point(round(x, 1), round(y_of(v), 1), b.label(unit), fmt(v))
        for x, v, b in zip(xs, values, spans, strict=True)
    ]
    line = "M" + " L".join(f"{p.x},{p.y}" for p in points)
    area = f"{line} L{points[-1].x},{bottom} L{points[0].x},{bottom} Z"
    yticks = [
        (round(y_of(i * step * divisor), 1), tick_text(i * step, axis_unit))
        for i in range(int(round(ceiling / step)) + 1)
    ]
    every = max(1, math.ceil(n / 6))
    xticks = [(xs[i], _short(spans[i], unit)) for i in range(0, n, every)]
    return LineChart(
        points, line, area, yticks, xticks, left, right, top, bottom, fmt(values[-1]),
        hit_width=span / max(n - 1, 1),
    )  # fmt: skip


def _short(bucket: Bucket, unit: str) -> str:
    if unit == "month":
        return f"{bucket.start:%b %Y}" if bucket.start.month == 1 else f"{bucket.start:%b}"
    return f"{bucket.start:%d %b}".lstrip("0")


@dataclass
class Segment:
    key: str
    label: str
    slot: int
    x: float
    width: float
    value: str
    share: str
    path: str
    text: str  # shown inside when it fits, else "" (the legend, tooltip and table carry it)

    @property
    def text_x(self) -> float:
        return self.x + 7


@dataclass
class StackRow:
    label: str
    y: float
    total: str
    segments: list[Segment]

    @property
    def mid_y(self) -> float:
        return self.y + BAR / 2


BAR = 22.0  # bar thickness (<= 24)
BAR_LEFT = 196.0  # bar charts: labels to the left of this


def rounded_end(x: float, y: float, w: float, h: float, r: float = 4.0) -> str:
    """A bar from a square baseline (left) to a rounded data end (right)."""
    r = max(0.0, min(r, w / 2, h / 2))
    return (
        f"M{x:.1f},{y:.1f} H{x + w - r:.1f} Q{x + w:.1f},{y:.1f} {x + w:.1f},{y + r:.1f} "
        f"V{y + h - r:.1f} Q{x + w:.1f},{y + h:.1f} {x + w - r:.1f},{y + h:.1f} H{x:.1f} Z"
    )


def stack_rows(rows: list[tuple[str, dict[str, float], Any]]) -> list[StackRow]:
    """100% stacked bars, one per (label, {kind: value}, formatter), 2px gaps between."""
    left, right = 84.0, WIDTH - 4.0
    out = []
    for i, (label, values, fmt) in enumerate(rows):
        y = 8 + i * (BAR + 30)
        total = sum(values.values())
        present = [(k, v) for k, v in values.items() if v > 0]
        gaps = 2.0 * max(len(present) - 1, 0)
        x = left
        segments = []
        for j, (key, value) in enumerate(present):
            width = (right - left - gaps) * value / total
            text = f"{KINDS[key][0]} {percent(value, total)}"
            last = j == len(present) - 1
            path = (
                rounded_end(x, y, width, BAR)
                if last
                else f"M{x:.1f},{y:.1f} h{width:.1f} v{BAR} h{-width:.1f} Z"
            )
            fits = len(text) * 7.4 + 14 <= width
            segments.append(
                Segment(key, KINDS[key][0], KINDS[key][1], round(x, 1), round(width, 1),
                        fmt(value), percent(value, total), path, text if fits else "")
            )  # fmt: skip
            x += width + 2
        out.append(StackRow(label, y, fmt(total), segments))
    return out


@dataclass
class Bar:
    label: str
    y: float
    width: float
    path: str
    value: str
    share: str

    @property
    def tip_x(self) -> float:
        return BAR_LEFT + self.width + 6

    @property
    def mid_y(self) -> float:
        return self.y + BAR / 2


def bar_rows(values: list[tuple[str, float]], fmt: Any) -> tuple[list[Bar], float]:
    """Horizontal bars from one baseline, values at the tips; (bars, chart height)."""
    left, right = BAR_LEFT, WIDTH - 90.0
    top = max((v for _, v in values), default=0) or 1
    total = sum(v for _, v in values)
    bars = []
    for i, (label, value) in enumerate(values):
        y = 6 + i * (BAR + 12)
        width = max((right - left) * value / top, 2.0)
        bars.append(Bar(label, y, width, rounded_end(left, y, width, BAR), fmt(value),
                        percent(value, total)))  # fmt: skip
    return bars, 6 + len(values) * (BAR + 12)


# --- The page's numbers ------------------------------------------------------------


@dataclass
class Tile:
    label: str
    value: str
    detail: str = ""
    change: str = ""  # vs the previous window, e.g. "+12% vs the 30 days before"


@dataclass
class Usage:
    filters: Filters
    unit: str
    tiles: list[Tile]
    cpu: LineChart
    storage: LineChart
    kinds: list[StackRow]
    kind_table: list[tuple[str, int, str]]  # (kind, jobs, CPU time), for the table view
    stages: list[Bar]
    stages_height: float
    largest_cpu: list[Job]
    largest_storage: list[Job]
    jobs_in_window: int
    unmeasured: int  # deleted jobs whose storage was never measured (deleted before 2026-10-16)
    present_kinds: list[str] = field(default_factory=list)
    # For the downloads (export_table): each bucket with its CPU seconds and the bytes held
    # at its end; the jobs behind "Most CPU time" and "Most storage now", every one of them.
    series: list[tuple[Bucket, float, float]] = field(default_factory=list)
    export_jobs: list[Job] = field(default_factory=list)
    window_ids: set[Any] = field(default_factory=set)


def compute(jobs: list[Job], filters: Filters, now: datetime) -> Usage:
    window = filters.window
    today = local_day(now)
    chosen = [job for job in jobs if filters.wants(job)]
    # Imported jobs ran elsewhere: they count only as storage, from when they arrived.
    ran = [j for j in chosen if j.imported_at is None]
    inside = [j for j in ran if window.start <= local_day(j.created_at) <= window.end]
    unit = unit_for(window)
    spans = buckets(window, unit)

    cpu_by_bucket = [0.0] * len(spans)
    index = {}
    for i, bucket in enumerate(spans):
        day = bucket.start
        while day <= bucket.end:
            index[day] = i
            day += timedelta(days=1)
    for job in inside:
        cpu_by_bucket[index[local_day(job.created_at)]] += job.cpu_seconds or 0

    # Storage held at the end of each bucket (now, for the current one): measured jobs
    # from when they finished (imported ones, from their import) until they were deleted.
    held = [
        (j.storage_bytes, since, j.deleted_at)
        for j in chosen
        if j.storage_bytes and (since := j.imported_at or j.finished_at)
    ]
    storage_by_bucket = []
    for bucket in spans:
        moment = now if bucket.end >= today else _end_of(bucket.end)
        storage_by_bucket.append(
            float(
                sum(
                    b for b, start, end in held if start <= moment and (end is None or end > moment)
                )
            )
        )

    total_cpu = sum(j.cpu_seconds or 0 for j in inside)
    failed_cpu = sum(j.cpu_seconds or 0 for j in inside if j.status in FAILED)
    live = [j for j in chosen if j.deleted_at is None]
    stored_now = sum(j.storage_bytes or 0 for j in live)
    freed = sum(j.storage_bytes or 0 for j in chosen if j.deleted_at is not None)
    succeeded = sum(j.status == Job.Status.SUCCEEDED for j in inside)
    failed = sum(j.status in FAILED for j in inside)

    change = ""
    if (before := window.previous()) is not None:
        earlier = sum(
            j.cpu_seconds or 0 for j in ran
            if before.start <= local_day(j.created_at) <= before.end
        )  # fmt: skip
        if earlier:
            delta = 100 * (total_cpu - earlier) / earlier
            change = f"{delta:+.0f}% vs the {window.days} days before"
        elif total_cpu:
            change = f"none in the {window.days} days before"

    tiles = [
        Tile("CPU time", duration(total_cpu), window.label.capitalize(), change),
        Tile("Storage now", size(stored_now),
             f"{len(live)} job{'s' if len(live) != 1 else ''}"
             + (f"; {size(freed)} freed by deleting jobs" if freed else "")),
        Tile("Jobs submitted", f"{len(inside):,}",
             f"{succeeded} succeeded, {failed} failed or cancelled"),
        Tile("CPU on failed jobs", duration(failed_cpu),
             f"{percent(failed_cpu, total_cpu)} of the CPU time" if total_cpu else window.label),
    ]  # fmt: skip

    counts: dict[str, tuple[int, float]] = {}
    for key in KINDS:
        members = [j for j in inside if kind_of(j.params) == key]
        if members:
            counts[key] = (len(members), sum(j.cpu_seconds or 0 for j in members))
    kinds = (
        stack_rows([
            ("Jobs", {k: float(n) for k, (n, _) in counts.items()}, lambda v: f"{v:,.0f}"),
            ("CPU time", {k: c for k, (_, c) in counts.items()}, duration),
        ])
        if inside else []
    )  # fmt: skip
    if kinds and not any(c for _, c in counts.values()):
        kinds = kinds[:1]  # no CPU recorded yet: the jobs row only

    by_stage: dict[str, float] = defaultdict(float)
    for job in inside:
        for task in (job.progress or {}).get("tasks", []):
            runtime = seconds(task.get("realtime", ""))
            if runtime is not None:
                by_stage[stage_of(task.get("name", ""))] += runtime * int(task.get("cpus") or 1)
    stage_values = [(STAGES[k], by_stage[k]) for k in STAGES if by_stage.get(k)]
    stages, stages_height = bar_rows(stage_values, duration)

    return Usage(
        filters=filters,
        unit=unit,
        tiles=tiles,
        cpu=line_chart(cpu_by_bucket, spans, unit, "cpu", duration),
        storage=line_chart(storage_by_bucket, spans, unit, "bytes", size),
        kinds=kinds,
        kind_table=[(KINDS[k][0], n, duration(c)) for k, (n, c) in counts.items()],
        stages=stages,
        stages_height=stages_height,
        largest_cpu=sorted(
            (j for j in inside if j.cpu_seconds), key=lambda j: -(j.cpu_seconds or 0)
        )[:8],
        largest_storage=sorted(
            (j for j in live if j.storage_bytes), key=lambda j: -(j.storage_bytes or 0)
        )[:8],
        jobs_in_window=len(inside),
        unmeasured=sum(1 for j in chosen if j.deleted_at is not None and j.storage_bytes is None),
        present_kinds=list(counts),
        series=list(zip(spans, cpu_by_bucket, storage_by_bucket, strict=True)),
        export_jobs=sorted(
            {j.id: j for j in [*inside, *live]}.values(),
            key=lambda j: (-(j.cpu_seconds or 0), -(j.storage_bytes or 0)),
        ),
        window_ids={j.id for j in inside},
    )  # fmt: skip


# The page's downloads (?export=<key>&format=csv|xlsx): "jobs" serves both job lists.
EXPORTS = {"cpu": "CPU time", "storage": "Storage", "jobs": "Jobs"}


def export_table(data: Usage, what: str) -> tuple[list[str], list[list[str]]]:
    """One panel's figures as (columns, rows), numbers unrounded (seconds, bytes)."""
    if what == "cpu":
        columns = [f"{data.unit.capitalize()} from", "To", "CPU seconds", "CPU time"]
        rows = [[f"{b.start:%Y-%m-%d}", f"{b.end:%Y-%m-%d}", f"{cpu:g}", duration(cpu)]
                for b, cpu, _ in data.series]  # fmt: skip
    elif what == "storage":
        columns = [f"End of {data.unit}", "Storage bytes", "Storage"]
        rows = [[f"{b.end:%Y-%m-%d}", f"{held:.0f}", size(held)] for b, _, held in data.series]
    elif what == "jobs":
        columns = ["Job ID", "Name", "Kind", "Status", "Submitted", "Finished", "Deleted",
                   "CPU seconds", "CPU time", "Storage bytes now", "Storage now",
                   "Submitted in this window"]  # fmt: skip
        rows = []
        for job in data.export_jobs:
            stored = (job.storage_bytes or 0) if job.deleted_at is None else 0
            rows.append([
                str(job.id), job.name or "", KINDS.get(kind_of(job.params), ("",))[0],
                str(job.status), _stamp(job.created_at), _stamp(job.finished_at),
                _stamp(job.deleted_at),
                f"{job.cpu_seconds:g}" if job.cpu_seconds is not None else "",
                duration(job.cpu_seconds or 0), str(stored), size(stored),
                "yes" if job.id in data.window_ids else "no",
            ])  # fmt: skip
    else:
        raise ValueError(what)
    return columns, rows


def _stamp(moment: datetime | None) -> str:
    return f"{timezone.localtime(moment):%Y-%m-%d %H:%M}" if moment else ""


def _end_of(day: date) -> datetime:
    start = datetime.combine(day + timedelta(days=1), datetime.min.time())
    return timezone.make_aware(start) - timedelta(microseconds=1)
