"""A mutational path drawn as a straight line.

On the tree a path turns corners through whatever lies between its two ends. Here it is
pulled straight, like a line on a transit map: a station for each ancestor (and the
sequence it may end at), the stretches between them as long as their branches, or all
the same length when the real lengths crowd the short ones together. Each stretch is a
band coloured by how widely its changes are shared with the job's other combinations
(`phylobot.path_compare`), so the line shows at a glance where the path is robust.

Pure geometry over numbers, as `phylobot.treeview` is; the template draws the SVG.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from phylobot.path_compare import LEVELS

WIDTH = 1000
PAD = 28  # either side of the line
MIN_STRETCH = 14  # px: a branch of no length still separates its two stations
BAND = 10  # the line's thickness where a branch has changes
THIN = 3  # and where it has none
LANE = 15  # px per row of station names
MAX_LANES = 8
CHAR = 6.3  # px per character of a station name at 11.5px
TOP = 10

# One colour per level (path_compare.LEVELS), light enough and far enough apart to tell
# without relying on hue alone: a blue ramp from every combination down to a few, a warm
# colour for a change no other combination makes, grey where nothing could be compared.
LEVEL_COLOURS = {
    "all": "#1f66c1",
    "most": "#6aa6e8",
    "few": "#bcd7f3",
    "none": "#eb6834",
    "na": "#b3bac5",
}


@dataclass(frozen=True)
class Band:
    x: float
    width: float
    level: str
    count: int

    @property
    def colour(self) -> str:
        return LEVEL_COLOURS[self.level]


@dataclass(frozen=True)
class Stretch:
    """One branch of the path."""

    index: int
    x1: float
    x2: float
    length: float  # the tree's, substitutions per site
    total: int  # its changes
    bands: list[Band]
    title: str

    @property
    def middle(self) -> float:
        return (self.x1 + self.x2) / 2

    @property
    def width(self) -> float:
        return round(self.x2 - self.x1, 2)

    @property
    def roomy(self) -> bool:
        """Wide enough to print the number of changes under it."""
        return self.x2 - self.x1 >= 8 + 7 * len(str(self.total))


@dataclass(frozen=True)
class Station:
    name: str
    label: str
    x: float
    extant: bool
    url: str
    lane: int = 0
    label_x: float = 0.0
    label_y: float = 0.0
    anchor: str = "middle"

    @property
    def stem_top(self) -> float:
        return self.label_y + 4


@dataclass
class Line:
    width: int
    height: float
    y: float  # the line's
    stations: list[Station]
    stretches: list[Stretch]
    equal: bool
    scale_width: float = 0.0  # px of the scale bar (0: none, the lengths are not real)
    scale_label: str = ""
    stretched: int = 0  # branches drawn longer than their length, to keep stations apart
    levels: list[tuple[str, str]] = field(default_factory=list)  # (level, colour) present

    @property
    def band_top(self) -> float:
        return self.y - BAND / 2

    @property
    def thin_top(self) -> float:
        return self.y - THIN / 2

    @property
    def stem_y(self) -> float:
        return self.y - 9

    @property
    def count_y(self) -> float:
        return self.y + 22

    @property
    def scale_y(self) -> float:
        return self.y + 42

    @property
    def scale_end(self) -> float:
        return PAD + self.scale_width

    @property
    def scale_text_x(self) -> float:
        return self.scale_end + 8

    @property
    def scale_text_y(self) -> float:
        return self.scale_y + 4


def _positions(lengths: list[float], equal: bool) -> tuple[list[float], float, int]:
    """Each station's x, the px per unit length (0 when equal), and how many branches
    were drawn at the minimum instead of their own length."""
    room = WIDTH - 2 * PAD
    n = len(lengths)
    if not n:
        return [WIDTH / 2], 0.0, 0
    floor = min(MIN_STRETCH, room / n)
    if equal or sum(lengths) <= 0:
        return [PAD + room * i / n for i in range(n + 1)], 0.0, 0

    def drawn(scale: float) -> list[float]:
        return [max(length * scale, floor) for length in lengths]

    lo, hi = 0.0, room / max(min(v for v in lengths if v > 0), 1e-12)
    for _ in range(60):  # the scale at which the stretches fill the room
        mid = (lo + hi) / 2
        if sum(drawn(mid)) > room:
            hi = mid
        else:
            lo = mid
    widths = drawn(lo)
    xs = [float(PAD)]
    for width in widths:
        xs.append(xs[-1] + width)
    return xs, lo, sum(1 for length in lengths if length * lo < floor - 1e-9)


def _scale_bar(scale: float) -> tuple[float, str]:
    """A round length whose bar is 50-130 px wide."""
    if scale <= 0:
        return 0.0, ""
    target = 90 / scale
    power = 10 ** math.floor(math.log10(target))
    value = min((1, 2, 5, 10), key=lambda m: abs(m * power - target)) * power
    return value * scale, f"{value:g} substitutions per site"


def _lanes(labels: list[tuple[float, str]]) -> list[tuple[int, float, str]]:
    """A row for each station's name so that none overlaps its neighbours: the lowest
    free row, names kept within the drawing. Returns (lane, x, text-anchor) per name."""
    ends: list[float] = []  # the right edge of the last name in each lane
    placed = []
    for x, text in labels:
        width = len(text) * CHAR
        left = min(max(x - width / 2, 2.0), WIDTH - 2.0 - width)
        lane = next((i for i, end in enumerate(ends) if end + 8 <= left), len(ends))
        lane = min(lane, MAX_LANES - 1)
        if lane == len(ends):
            ends.append(left + width)
        else:
            ends[lane] = max(ends[lane], left + width)
        placed.append((lane, left + width / 2, "middle"))
    return placed


def line(
    names: list[str],
    labels: list[str],
    extant: set[str],
    urls: list[str],
    lengths: list[float],
    levels: list[dict[str, int]],
    titles: list[str],
    equal: bool = False,
) -> Line:
    """The drawing. `names` is the path, oldest first; `labels` what each station is
    called; `urls` where each leads ("" for none); `lengths` the tree's branch lengths
    and `levels` each branch's changes counted by level (`Comparison.levels`), both one
    per branch; `titles` each branch's tooltip."""
    xs, scale, stretched = _positions(lengths, equal)
    placed = _lanes(list(zip(xs, labels, strict=True)))
    lanes = max((lane for lane, _, _ in placed), default=0) + 1
    y = TOP + lanes * LANE + 14
    stations = [
        Station(
            name, label, round(x, 1), name in extant, url, lane,
            round(label_x, 1), round(y - 14 - lane * LANE, 1), anchor,
        )
        for name, label, x, url, (lane, label_x, anchor) in zip(
            names, labels, xs, urls, placed, strict=True
        )
    ]  # fmt: skip
    stretches = []
    present: dict[str, None] = {}
    for i, (counts, title) in enumerate(zip(levels, titles, strict=True)):
        x1, x2 = xs[i], xs[i + 1]
        total = sum(counts.values())
        bands = []
        at = x1
        for level in LEVELS:
            count = counts.get(level, 0)
            if not count:
                continue
            present[level] = None
            width = (x2 - x1) * count / total
            bands.append(Band(round(at, 2), round(width, 2), level, count))
            at += width
        stretches.append(
            Stretch(
                i, round(x1, 1), round(x2, 1), lengths[i] if lengths else 0.0, total, bands, title
            )
        )
    bar, text = (0.0, "") if equal else _scale_bar(scale)
    return Line(
        WIDTH, y + (58 if bar else 34), y, stations, stretches, equal or scale <= 0,
        round(bar, 1), text, stretched,
        [(level, LEVEL_COLOURS[level]) for level in LEVELS if level in present],
    )  # fmt: skip
