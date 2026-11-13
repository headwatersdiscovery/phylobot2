"""Histograms beside the homolog filters, so the user can see where a threshold falls.

One histogram per percentage filter, of the hits that reach that filter (the filters
run in order): query coverage for the hits within the E-value cut-off, then identity
for those that also pass coverage. Each 5% bin is split into the hits the current
settings keep (blue) and remove (grey); the kept range is shaded and each threshold
marked. Server-drawn SVG like the ancestor viewer's plots (no chart library); styles
are SVG attributes so a downloaded image matches the page.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass

from phylobot.ancestors import Frame, nice_step
from phylobot.homologs.filters import Settings
from phylobot.homologs.hits import Hit

SIZE = (440, 236)
MARGIN = (52, 14, 30, 44)  # left, right, top (threshold labels), bottom
BIN = 5  # % per bar
KEPT, REMOVED = "#2457a6", "#b3bac5"
GAP = 2  # px of background between a bar's two segments, and between bars


@dataclass(frozen=True)
class Segment:
    y: float
    height: float


@dataclass(frozen=True)
class Bar:
    x: float
    width: float
    kept: Segment | None
    removed: Segment | None
    label: str  # the bin, e.g. "55-60%"
    kept_count: int
    removed_count: int


@dataclass(frozen=True)
class Marker:
    x: float
    label: str
    anchor: str  # text-anchor for the label, so labels at the edges stay inside


@dataclass(frozen=True)
class Histogram:
    key: str
    title: str
    x_label: str
    frame: Frame
    bars: list[Bar]
    x_ticks: list[tuple[float, str]]
    y_ticks: list[tuple[float, str]]
    band: tuple[float, float]  # the kept range: x and width
    markers: list[Marker]
    total: int  # hits reaching this filter
    kept: int


def histogram(
    key: str,
    title: str,
    x_label: str,
    hits: list[Hit],
    value: Callable[[Hit], float],
    low: float,
    high: float,
) -> Histogram:
    """Bin `value` of each hit; hits with low <= value <= high are kept."""
    frame = Frame.of(SIZE, MARGIN)
    bins = 100 // BIN
    kept = [0] * bins
    removed = [0] * bins
    for hit in hits:
        v = value(hit)
        slot = min(max(int(v // BIN), 0), bins - 1)
        (kept if low <= v <= high else removed)[slot] += 1

    step = nice_step(max([k + r for k, r in zip(kept, removed, strict=True)] + [1]), 4)
    top = step * math.ceil(max([k + r for k, r in zip(kept, removed, strict=True)] + [1]) / step)

    def y(count: float) -> float:
        return round(frame.bottom - count / top * frame.plot_height, 1)

    def x(percent: float) -> float:
        return round(frame.left + percent / 100 * frame.plot_width, 1)

    slot_width = frame.plot_width / bins
    bars = []
    for i in range(bins):
        k, r = kept[i], removed[i]
        kept_segment = Segment(y(k), round(frame.bottom - y(k), 1)) if k else None
        removed_segment = None
        if r:
            base = y(k) - (GAP if k else 0)
            height = frame.bottom - y(r)
            removed_segment = Segment(round(base - height, 1), round(height, 1))
        bars.append(
            Bar(
                x=round(frame.left + i * slot_width + GAP / 2, 1),
                width=round(slot_width - GAP, 1),
                kept=kept_segment,
                removed=removed_segment,
                label=f"{i * BIN}-{(i + 1) * BIN}%",
                kept_count=k,
                removed_count=r,
            )
        )
    markers = []
    for threshold in (low, high):
        if 0 < threshold < 100:
            anchor = "start" if threshold < 15 else "end" if threshold > 85 else "middle"
            markers.append(Marker(x(threshold), f"{threshold:g}%", anchor))
    return Histogram(
        key=key,
        title=title,
        x_label=x_label,
        frame=frame,
        bars=bars,
        x_ticks=[(x(p), f"{p}") for p in range(0, 101, 20)],
        y_ticks=[(y(c), f"{c:,}") for c in range(0, top + 1, step)],
        band=(x(low), round(x(high) - x(low), 1)),
        markers=markers,
        total=len(hits),
        kept=sum(kept),
    )


def filter_histograms(hits: list[Hit], settings: Settings) -> list[Histogram]:
    """The coverage and identity histograms for these hits under these settings."""
    s = settings
    reach_coverage = [h for h in hits if h.evalue <= s.max_evalue]
    reach_identity = [h for h in reach_coverage if h.coverage >= s.min_coverage]
    return [
        histogram(
            "coverage",
            "Query coverage",
            "% of the query aligned",
            reach_coverage,
            lambda h: h.coverage,
            s.min_coverage,
            100,
        ),
        histogram(
            "identity",
            "Identity to the query",
            "% identical in the alignment",
            reach_identity,
            lambda h: h.identity,
            s.min_identity,
            s.max_identity,
        ),
    ]
