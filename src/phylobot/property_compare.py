"""One ancestor's properties across every combination that recovered it.

The same biological node is reconstructed once per alignment x tree x model combination.
This module takes each one's draws for a property and lays them out for comparison: a
**shared x-axis**, an overlay of every distribution, a small multiple per combination,
and a table of the numbers.

Pure: it takes already-read values and returns SVG geometry, in the `treeview` /
`property_views` idiom. Matching which nodes *are* the same ancestor is `phylobot.analogs`
and happens before this.

**The axis never moves when series are toggled.** Both the bins and the y-scale are
computed over every candidate, shown or not, so hiding a combination cannot rescale the
plot underneath the ones left. Otherwise two readings of the same page would not be
comparable — which is the one thing this page exists to make possible.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

#: The plot's drawing area, in viewBox units.
WIDTH = 760
HEIGHT = 240
LEFT = 42.0
RIGHT = 12.0
TOP = 10.0
BOTTOM = 30.0
#: Bins for a continuous property; a discrete one gets a bin per integer (`bins`).
BINS = 28
DISCRETE_VALUES = 20
#: How many series colours there are (`.viz-line.s1` ... `.s5` in site.css, a palette
#: already validated for colour blindness). Colours repeat beyond this.
COLOURS = 5


#: The histogram palettes offered on the comparison page. Each is five hues that stay
#: distinguishable when laid over each other with alpha; the colours themselves live in
#: site.css under `.compare[data-palette="<key>"]`, so a skin can still adjust them.
PALETTES: dict[str, tuple[str, str]] = {
    "bright": (
        "Bright",
        "PhyloBot's default: five saturated hues, the widest separation of the three.",
    ),
    "okabe": (
        "Colour-blind safe",
        "Okabe & Ito's qualitative palette, designed to stay distinct under the common "
        "forms of colour blindness. The safest choice for a figure others will read.",
    ),
    "ink": (
        "Print / greyscale",
        "One hue at five lightnesses: the series stay separable when a figure is "
        "printed in black and white or photocopied, where hue is lost.",
    ),
}
DEFAULT_PALETTE = "bright"


def palette_of(asked: object) -> str:
    """The chosen palette, or the default for anything unknown (raw request text)."""
    return asked if isinstance(asked, str) and asked in PALETTES else DEFAULT_PALETTE


# How many combinations per alignment the page compares, best by BIC first (owner,
# 2026-10-30): 285 overlaid distributions are unreadable and slow to read and draw, so
# the default is a few; more than that warns before it loads (compare.js).
TOP_CHOICES = (1, 2, 3, 5, 10, 20)
DEFAULT_TOP = 3
ALL = "all"


def top_of(asked: object) -> int | None:
    """Combinations per alignment to compare, None for all (raw request text)."""
    if asked == ALL:
        return None
    try:
        n = int(str(asked))
    except ValueError:
        return DEFAULT_TOP
    return n if n in TOP_CHOICES else DEFAULT_TOP


def within_top(combo: dict[str, Any], top: int | None) -> bool:
    """Whether a combination is among its alignment's `top` by BIC (`rank_bic` ranks
    within each alignment and data type; a combination without a fit is last)."""
    if top is None:
        return True
    rank = combo.get("rank_bic")
    return isinstance(rank, int) and rank <= top


def number(value: float) -> str:
    """A property value for display, as `property_views.number`."""
    if value == int(value) and abs(value) < 1e6:
        return str(int(value))
    text = f"{value:.3g}"
    return "0" if text in {"-0", "-0.0"} else text


def bins(values: list[float]) -> tuple[float, float, int]:
    """(low, high, count) over the pooled values of every candidate."""
    low, high = min(values), max(values)
    if len(set(values)) <= DISCRETE_VALUES and all(v == int(v) for v in values):
        return low - 0.5, high + 0.5, int(high - low) + 1
    if high <= low:
        return low - 0.5, high + 0.5, 1
    return low, high, BINS


def counts(values: list[float], low: float, high: float, count: int) -> list[float]:
    """`values` as proportions per bin — proportions, not counts, so series with
    different numbers of draws can be laid over each other."""
    width = (high - low) / count
    out = [0.0] * count
    for value in values:
        index = int((value - low) / width) if width else 0
        out[min(max(index, 0), count - 1)] += 1
    return [c / len(values) for c in out] if values else out


@dataclass(frozen=True)
class Axis:
    """The shared x-axis: the same for the overlay and every small multiple."""

    low: float
    high: float
    count: int
    height: float
    ticks: list[tuple[float, str]]  # (x, label)
    #: Space reserved above the plot for the mean labels, in viewBox units.
    headroom: float = 0.0

    @property
    def plot_width(self) -> float:
        return WIDTH - LEFT - RIGHT

    @property
    def left(self) -> float:
        return LEFT

    @property
    def right(self) -> float:
        return WIDTH - RIGHT

    @property
    def top(self) -> float:
        return TOP + self.headroom

    @property
    def baseline(self) -> float:
        return self.height - BOTTOM

    @property
    def width(self) -> float:
        return WIDTH

    def x(self, value: float) -> float:
        span = self.high - self.low
        return LEFT + ((value - self.low) / span if span else 0.5) * self.plot_width


def axis_of(
    low: float, high: float, count: int, height: float = HEIGHT, headroom: float = 0.0
) -> Axis:
    """Five labelled ticks across the axis, in the property's own units."""
    span = high - low
    ticks = []
    for i in range(5):
        value = low + span * i / 4
        x = LEFT + (i / 4) * (WIDTH - LEFT - RIGHT)
        ticks.append((round(x, 1), number(value)))
    return Axis(low, high, count, height + headroom, ticks, headroom)


#: Label geometry for the mean lines.
LANE_HEIGHT = 16.0
LABEL_FONT = 10.5
#: Rough width of a character at LABEL_FONT, for collision tests. An estimate is enough:
#: the lanes only have to stop labels from overlapping, and erring wide is harmless.
CHAR_WIDTH = 5.6
LABEL_GAP = 4.0


@dataclass(frozen=True)
class MeanMark:
    """One series' posterior mean, as a line over the histogram with a label above it."""

    colour: int
    x: float
    text: str
    label_x: float  # clamped so the text stays inside the plot
    lane: int  # 0 is nearest the plot; higher lanes sit further above
    y: float  # the label's baseline


def lanes(
    marks: list[tuple[int, float, str]],
    left: float = LEFT,
    right: float = RIGHT,
    width: float = WIDTH,
) -> list[tuple[int, float, str, float, int]]:
    """Place each label in the lowest lane where it does not touch one already there.

    Means cluster — that is the usual finding, not an edge case — so labels at nearly the
    same x are the normal case, and nudging them sideways would detach a label from the
    line it names. Stacking instead keeps every label directly over its own line.

    The bounds are parameters because the ancestor viewer's own plots label their markers
    the same way (`property_views.plot`) on a narrower axis.
    """
    placed: list[list[tuple[float, float]]] = []  # per lane: (left, right) taken
    out = []
    for colour, x, text in sorted(marks, key=lambda m: m[1]):
        half = len(text) * CHAR_WIDTH / 2
        label_x = min(max(x, left + half), width - right - half)
        starts, ends = label_x - half - LABEL_GAP, label_x + half + LABEL_GAP
        lane = 0
        while True:
            if lane == len(placed):
                placed.append([])
            if all(ends <= a or starts >= b for a, b in placed[lane]):
                placed[lane].append((starts, ends))
                break
            lane += 1
        out.append((colour, x, text, round(label_x, 1), lane))
    return out


@dataclass(frozen=True)
class Series:
    """One combination's distribution of one property."""

    key: str  # the combination
    node: str
    label: str
    colour: int  # 1..COLOURS, the .sN class
    shown: bool
    outline: str  # SVG polyline points: the histogram as a step
    area: str  # SVG polygon points: the same, closed to the baseline, drawn with alpha
    mean_x: float
    ml_x: float
    #: Numbers for the table, already formatted.
    mean: str
    sd: str
    low: str
    high: str
    ml: str
    shift: str
    draws: int


@dataclass(frozen=True)
class Panel:
    """One property, every candidate."""

    key: str
    label: str
    units: str
    axis: Axis
    series: list[Series]
    level: str
    #: The posterior means as vertical lines, labelled above the plot. Empty when the
    #: reader has switched them off.
    means: list[MeanMark]

    @property
    def drawn(self) -> list[Series]:
        """Only the series the reader has left switched on."""
        return [s for s in self.series if s.shown]

    @property
    def spread(self) -> str:
        """The range of the shown series' means — the headline number: how far method
        choice moves this property, in the property's own units."""
        means = [float(s.mean) for s in self.drawn] if self.drawn else []
        if len(means) < 2:
            return ""
        return f"{number(min(means))} to {number(max(means))}"


@dataclass(frozen=True)
class Candidate:
    """One combination that recovered this ancestor, with its draws and summary."""

    key: str
    node: str
    label: str
    aligner: str
    shown: bool
    values: dict[str, list[float]]  # property key -> draws
    summary: dict[str, Any]


def step(axis: Axis, shares: list[float], top: float) -> tuple[str, str]:
    """A histogram as a step outline and the same closed to the baseline.

    Steps rather than bars: a dozen bar series over one axis is unreadable, while
    outlines stay legible and the filled areas show where they agree.
    """
    span = (axis.high - axis.low) / axis.count
    plot_height = axis.baseline - axis.top
    points = []
    for i, share in enumerate(shares):
        y = axis.baseline - (share / top if top else 0) * plot_height
        points.append((axis.x(axis.low + i * span), round(y, 1)))
        points.append((axis.x(axis.low + (i + 1) * span), round(y, 1)))
    outline = " ".join(f"{round(x, 1)},{y}" for x, y in points)
    closed = f"{round(points[0][0], 1)},{axis.baseline} {outline} "
    closed += f"{round(points[-1][0], 1)},{axis.baseline}"
    return outline, closed


def panel(
    key: str,
    label: str,
    units: str,
    candidates: list[Candidate],
    level: str = "95",
    show_means: bool = True,
) -> Panel | None:
    """One property laid out over every candidate that measured it."""
    having = [c for c in candidates if c.values.get(key)]
    if not having:
        return None
    pooled = [v for c in having for v in c.values[key]]
    low, high, count = bins(pooled)
    shares = {c.key: counts(c.values[key], low, high, count) for c in having}
    # The tallest bar of ANY candidate, so switching one off does not stretch the rest.
    top = max((max(s) for s in shares.values() if s), default=1.0) or 1.0

    # The mean labels are laid out first: how many lanes they need sets how much room the
    # plot has to leave above itself, and the lane count is not known until they are
    # placed. x depends only on the value and the axis's ends, not on its height, so
    # this is not circular.
    plain = axis_of(low, high, count)
    wanted = []
    if show_means:
        for index, c in enumerate(having):
            entry = c.summary.get("properties", {}).get(key, {})
            try:
                value = float(entry["mean"])
            except KeyError, TypeError, ValueError:
                continue
            if c.shown:
                # The combination's name, not the number: the question these lines
                # answer is which combination sits where relative to the others, and
                # every value is in the table below.
                wanted.append((index % COLOURS + 1, plain.x(value), c.label))
    laid = lanes(wanted)
    used = max((lane for *_, lane in laid), default=-1) + 1
    headroom = used * LANE_HEIGHT

    axis = axis_of(low, high, count, HEIGHT, headroom)
    means = [
        MeanMark(
            colour=colour,
            x=round(x, 1),
            text=text,
            label_x=label_x,
            lane=lane,
            # Lane 0 sits just above the plot; later lanes stack upwards.
            y=round(TOP + headroom - lane * LANE_HEIGHT - 3, 1),
        )  # fmt: skip
        for colour, x, text, label_x, lane in laid
    ]
    series = []
    for index, c in enumerate(having):
        entry = c.summary.get("properties", {}).get(key, {})

        def figure(name: str, entry: dict[str, Any] = entry) -> str:
            try:
                return number(float(entry[name]))
            except KeyError, TypeError, ValueError:
                return ""

        outline, area = step(axis, shares[c.key], top)
        try:
            standardised = float(entry["standardised_shift"])
            shift = f"{standardised:+.2f}"
        except KeyError, TypeError, ValueError:
            shift = ""
        mean = figure("mean")
        ml = figure("ml")
        series.append(
            Series(
                key=c.key,
                node=c.node,
                label=c.label,
                colour=index % COLOURS + 1,
                shown=c.shown,
                outline=outline,
                area=area,
                mean_x=round(axis.x(float(mean)), 1) if mean else axis.left,
                ml_x=round(axis.x(float(ml)), 1) if ml else axis.left,
                mean=mean,
                sd=figure("sd"),
                low=figure("low"),
                high=figure("high"),
                ml=ml,
                shift=shift,
                draws=len(c.values[key]),
            )  # fmt: skip
        )
    return Panel(key, label, units, axis, series, level, means)
