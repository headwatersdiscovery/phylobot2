"""The ancestor viewer's property panel.

PROPERTIES writes, per combination, node and gap version, every property's draws plus a
summary (`phylobot.pipeline.properties`). This module turns those files into what the
panel draws: per property a **histogram of the draws**, with the posterior mean, f(ML)
and f(AltAll) marked as vertical lines, as server-side SVG geometry in the `treeview` /
`ancestors.node_plots` idiom.

Pure but for the storage reads in `read` and `read_histograms`: everything below takes
parsed data.

The draws are read as well as the summaries — a histogram is the point of the panel, and
a summary cannot be rebinned. One node's draws are a few ranged reads of the
combination's draws.parquet (`phylobot.property_store`). `read_histograms` tolerates
their absence: a node without draws still gets its numbers and markers, without bars.

f(AltAll) is computed per request rather than read, because AltAll depends on the
viewer's threshold slider, which the pipeline never saw.
"""

from __future__ import annotations

import functools
import io
from dataclasses import dataclass
from typing import Any

from phylobot import property_store, storage
from phylobot.property_compare import LANE_HEIGHT, lanes

#: The plot's drawing area, in viewBox units: the SVG scales to its column, so these set
#: the aspect ratio and the label spacing, not a pixel size.
PLOT_WIDTH = 460
PLOT_HEIGHT = 150
MARGIN_LEFT = 8.0
MARGIN_RIGHT = 8.0
MARGIN_TOP = 10.0
#: Room under the axis for the end labels.
MARGIN_BOTTOM = 22.0
#: Bars wider than this get a gap between them; narrower ones are drawn touching.
BAR_GAP = 1.0
#: How many bins a continuous property gets. A count property is binned per integer
#: instead (`bins`), since half-integer bars would misrepresent a discrete quantity.
BINS = 24
#: At most this many distinct values and a property counts as discrete.
DISCRETE_VALUES = 15
#: How much of the span to leave either side, so a marker at an end is not clipped.
PAD = 0.04


SUMMARIES_CACHED = 64  # summaries.json files kept per process (finished jobs' never change)


def has(keys: set[str], combination: str) -> bool:
    """Whether a combination's properties were measured here (`keys`: results files)."""
    return property_store.summaries_path(combination) in keys


@functools.lru_cache(maxsize=SUMMARIES_CACHED)
def _summaries(key: str) -> str:
    """One combination's summaries.json: read once per process, for every node of it."""
    return storage.read_bytes(key).decode()


def read(prefix: str, combination: str, node: str, version: str) -> dict[str, Any] | None:
    """The stored summary, or None when this job did not measure properties here (the
    box was unticked, or not this node). Unreadable or malformed JSON is also None: the
    panel is an extra, and must never break the viewer it sits in."""
    try:
        text = _summaries(prefix + property_store.summaries_path(combination))
    except Exception:
        return None
    return property_store.summary(text, node, version)


def read_histograms(
    prefix: str, combination: str, node: str, version: str, keys: list[str]
) -> dict[str, list[float]]:
    """Each property's draws, by ranged reads of the combination's draws.parquet. A
    property without draws is left out, and its row is then drawn without bars."""
    key = prefix + property_store.draws_path(combination)
    try:
        ranged = property_store.RangedFile(
            storage.size(key), lambda start, n: storage.read_range(key, start, n)
        )
        return property_store.read_draws(io.BufferedReader(ranged, READ_SIZE), node, version, keys)
    except Exception:
        return {}


READ_SIZE = 256 * 1024  # bytes per ranged read of draws.parquet


def number(value: float) -> str:
    """A property value for display. Properties range from counts (3) through charge
    (-3.24) to the hydrophobic moment (0.0356), so significant figures, not decimals —
    but integers stay integers rather than becoming 3.00."""
    if value == int(value) and abs(value) < 1e6:
        return str(int(value))
    text = f"{value:.3g}"
    return "0" if text in {"-0", "-0.0"} else text


@dataclass(frozen=True)
class Bar:
    """One histogram bar, in viewBox coordinates."""

    x: float
    y: float
    width: float
    height: float
    count: int
    share: float
    label: str  # the bin, for the tooltip


@dataclass(frozen=True)
class Marker:
    """A vertical line over the histogram: the posterior mean, f(ML) or f(AltAll)."""

    kind: str  # "mean" | "ml" | "altall"
    label: str
    x: float
    value: str
    inside: bool  # within the drawn axis, so the template can say when it is not
    #: The label drawn above the line, so a reader need not look back to a legend
    #: (owner, 2026-10-21). Stacked into lanes exactly as the comparison page's mean
    #: labels are, because these three cluster for the same reason.
    label_x: float = 0.0
    label_y: float = 0.0


@dataclass(frozen=True)
class Plot:
    """One property's histogram and its markers."""

    bars: list[Bar]
    markers: list[Marker]
    low_label: str
    high_label: str
    #: True when the property took few enough distinct values to draw one bar each,
    #: rather than binning a continuum.
    discrete: bool
    #: Room the marker labels need above the plot, in SVG units. The labels are drawn
    #: at negative y and the viewBox is extended upwards to show them, so every bar and
    #: the baseline keep the coordinates they always had.
    headroom: float = 0.0

    @property
    def svg_height(self) -> float:
        return PLOT_HEIGHT + self.headroom

    @property
    def view_box(self) -> str:
        return f"0 {-self.headroom:g} {PLOT_WIDTH} {self.svg_height:g}"


def bins(values: list[float]) -> tuple[float, float, int]:
    """(low, high, count) for the histogram's bins.

    A property with few distinct values — `hydropathy_regions` is a count, and on one
    ancestor takes perhaps six values — gets one bin per integer, offset by a half so
    each bar sits over its value. Anything else gets `BINS` equal bins over its range.
    """
    low, high = min(values), max(values)
    if len(set(values)) <= DISCRETE_VALUES and all(v == int(v) for v in values):
        return low - 0.5, high + 0.5, int(high - low) + 1
    if high <= low:  # every draw identical, and not an integer
        return low - 0.5, high + 0.5, 1
    return low, high, BINS


def plot(values: list[float], marks: list[tuple[str, str, float]]) -> Plot:
    """The histogram of `values`, with `marks` ((kind, label, value), ...) over it.

    The axis spans the draws *and* every mark, because f(ML) or f(AltAll) can fall
    outside the drawn range; a marker that silently vanished would be the worst outcome,
    since a reader would take its absence for agreement.
    """
    low, high, count = bins(values)
    discrete = count != BINS
    span_low = min([low, *(v for _, _, v in marks)])
    span_high = max([high, *(v for _, _, v in marks)])
    pad = (span_high - span_low) * PAD or 0.5
    axis_low, axis_high = span_low - pad, span_high + pad
    axis_span = axis_high - axis_low
    plot_width = PLOT_WIDTH - MARGIN_LEFT - MARGIN_RIGHT
    plot_height = PLOT_HEIGHT - MARGIN_TOP - MARGIN_BOTTOM

    def x(value: float) -> float:
        return MARGIN_LEFT + (value - axis_low) / axis_span * plot_width

    counts = [0] * count
    width = (high - low) / count
    for value in values:
        index = int((value - low) / width) if width else 0
        counts[min(max(index, 0), count - 1)] += 1
    tallest = max(counts) or 1
    bars = []
    for i, n in enumerate(counts):
        left, right = x(low + i * width), x(low + (i + 1) * width)
        height = n / tallest * plot_height
        gap = BAR_GAP if right - left > 3 else 0.0
        bars.append(Bar(
            x=round(left, 1), y=round(MARGIN_TOP + plot_height - height, 1),
            width=round(max(right - left - gap, 0.6), 1), height=round(height, 1),
            count=n, share=n / len(values),
            label=f"{number(low + i * width)} to {number(low + (i + 1) * width)}",
        ))  # fmt: skip

    # Lay the labels out before the markers: the number of lanes they need is how far
    # the topmost sits above the plot, and `headroom` is what the template leaves for it.
    placed = lanes(
        [(i, x(value), label) for i, (_, label, value) in enumerate(marks)],
        left=MARGIN_LEFT,
        right=MARGIN_RIGHT,
        width=PLOT_WIDTH,
    )
    used = max((lane for *_, lane in placed), default=-1) + 1
    by_index = {i: (label_x, lane) for i, _, _, label_x, lane in placed}
    markers = []
    for i, (kind, label, value) in enumerate(marks):
        label_x, lane = by_index[i]
        markers.append(Marker(
            kind, label, round(x(value), 1), number(value),
            axis_low <= value <= axis_high,
            label_x=label_x,
            # Lane 0 sits just above the plot; higher lanes stack upwards from it.
            # Negative: above the plot, lane 0 nearest it.
            label_y=round(-(lane * LANE_HEIGHT) - 5, 1),
        ))  # fmt: skip
    return Plot(bars, markers, number(low), number(high), discrete, headroom=used * LANE_HEIGHT)


@dataclass(frozen=True)
class Row:
    """One property in the panel."""

    key: str
    label: str
    units: str
    mean: str
    sd: str
    low: str
    high: str
    ml: str
    altall: str  # "" when this reconstruction offers no AltAll
    #: f(ML) - mean, in standard deviations, formatted with its sign. Substantial in
    #: magnitude but random in sign, so it is shown as a distance and never called a bias.
    shift: str
    shift_notable: bool
    #: How the number is computed and what it is not, with citations: an `about.Method`,
    #: or None for a property with no entry there.
    method: Any
    plot: Plot | None
    #: How many draws this property was measured on. Per property, not per panel: a
    #: later run may have used another number (`merge`).
    draws: int = 0
    #: Where this number came from, when it was not the job's own run: "measured later,
    #: on <date>". Empty for the job's own measurement, which is
    #: the ordinary case and needs no note.
    source: str = ""


#: Where the offset stops being a rounding detail. It sits below the typical value on
#: purpose: it is a "look at this" threshold, not a significance test.
NOTABLE_SHIFT = 0.25


def merge(summaries: list[dict[str, Any] | None]) -> dict[str, Any] | None:
    """Several stored summaries of the same ancestor, read in order, into one.

    A property measured after the job finished lives under its own run's prefix, so the
    panel is assembled from more than one file. Later
    summaries win per property: asking for a property again is how a reader replaces a
    measurement, and the newest answer is the one the registry's current code produced.

    The surrounding fields (draws, sampling, level) come from the last summary that
    carries each property, so the caption describes the run the numbers came from.
    """
    kept: dict[str, Any] = {}
    merged: dict[str, Any] = {}
    for summary in summaries:
        if not summary or not isinstance(summary.get("properties"), dict):
            continue
        if not merged:
            merged = {k: v for k, v in summary.items() if k != "properties"}
        for key, entry in summary["properties"].items():
            # Each property keeps the draw count of the run that measured it: runs may
            # have been asked for different numbers, and the caption must not claim one
            # of them for all three.
            kept[key] = {"draws": summary.get("draws"), **entry}
    if not kept:
        return None
    return merged | {"properties": kept}


def rows(
    summary: dict[str, Any],
    histograms: dict[str, list[float]] | None = None,
    altall: dict[str, float] | None = None,
    methods: dict[str, Any] | None = None,
    sources: dict[str, str] | None = None,
) -> list[Row]:
    """The panel's rows, in the registry's order as PROPERTIES wrote them.

    A property whose stored entry is missing a number is skipped rather than drawn with
    zeros: the summaries are written by our own pipeline, so a gap means a version
    mismatch, and a silently wrong figure is worse than an absent one.
    """
    histograms, altall, methods = histograms or {}, altall or {}, methods or {}
    sources = sources or {}
    # A summary read straight from one file carries its draw count once, for every
    # property; `merge` stamps each entry with its own run's count instead.
    everyones = summary.get("draws")
    out = []
    for key, entry in summary.get("properties", {}).items():
        try:
            low, high = float(entry["low"]), float(entry["high"])
            mean, sd, ml = float(entry["mean"]), float(entry["sd"]), float(entry["ml"])
            standardised = float(entry["standardised_shift"])
        except KeyError, TypeError, ValueError:
            continue
        marks = [("mean", "posterior mean", mean), ("ml", "f(ML)", ml)]
        if key in altall:
            marks.append(("altall", "f(AltAll)", altall[key]))
        values = histograms.get(key)
        out.append(
            Row(
                key=key,
                label=str(entry.get("label") or key),
                units=str(entry.get("units") or ""),
                mean=number(mean),
                sd=number(sd),
                low=number(low),
                high=number(high),
                ml=number(ml),
                altall=number(altall[key]) if key in altall else "",
                shift=f"{standardised:+.2f}",
                shift_notable=abs(standardised) >= NOTABLE_SHIFT,
                method=methods.get(key),
                draws=int(entry.get("draws") or everyones or 0),
                plot=plot(values, marks) if values else None,
                source=sources.get(key, ""),
            )  # fmt: skip
        )
    return out


@dataclass(frozen=True)
class Panel:
    """Everything the property panel needs, or None from `of` when there is nothing."""

    rows: list[Row]
    #: The draw count every property shares, or 0 when they differ — then the caption
    #: says so and each row carries its own (`Row.draws`).
    draws: int
    level: str  # the credible level as a percentage, e.g. "95"
    #: The sampling caveat, as the pipeline stored it. Shown verbatim, because a figure
    #: made from these numbers has to carry it.
    sampling: str
    weighted: bool
    #: The AltAll threshold f(AltAll) was computed at, for the caption.
    altall_threshold: str = ""
    width: int = PLOT_WIDTH
    height: int = PLOT_HEIGHT

    # Geometry the template would otherwise have to compute; Django templates cannot do
    # arithmetic, and `add` on floats is a trap.
    @property
    def left(self) -> float:
        return MARGIN_LEFT

    @property
    def right(self) -> float:
        return self.width - MARGIN_RIGHT

    @property
    def top(self) -> float:
        return MARGIN_TOP

    @property
    def baseline(self) -> float:
        return self.height - MARGIN_BOTTOM

    @property
    def label_y(self) -> float:
        return self.baseline + 14


def of(
    summary: dict[str, Any] | None,
    histograms: dict[str, list[float]] | None = None,
    altall: dict[str, float] | None = None,
    methods: dict[str, Any] | None = None,
    altall_threshold: str = "",
    sources: dict[str, str] | None = None,
) -> Panel | None:
    """The panel for one stored summary, or None when there is nothing to draw."""
    if not summary:
        return None
    drawn = rows(summary, histograms, altall, methods, sources)
    if not drawn:
        return None
    first: dict[str, Any] = next(iter(summary["properties"].values()), {})
    try:
        level = float(first.get("level") or 0.95)
    except TypeError, ValueError:
        level = 0.95
    counts = {row.draws for row in drawn}
    draws = counts.pop() if len(counts) == 1 else 0
    return Panel(
        rows=drawn,
        draws=draws,
        level=number(level * 100),
        sampling=str(summary.get("sampling") or ""),
        weighted=bool(summary.get("weighted")),
        altall_threshold=altall_threshold if any(r.altall for r in drawn) else "",
    )
