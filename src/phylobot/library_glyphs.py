"""Glyphs for the library page (phylobot.library_views): server-drawn, no chart library.

Three views of the same thing — where a library's ancestors are uncertain, since that is
what the library has to cover:

- `overview`: one thin strip per target, laid out by alignment column and coloured by
  the ML state's posterior probability, so uncertainty that piles up along a lineage
  shows as a column of warm cells.
- `sparkline`: per target, the ML state's probability along the ancestor, with a tick at
  each ambiguous site.
- `site_map`: the ancestral sequence itself, as tiles coloured like the ancestor viewer's,
  with the alternative states stacked beneath every ambiguous site.
- `length_chart`: the extant sequences' lengths as a histogram, with every ancestor's
  length marked on it, so one far longer than any real protein stands out (and
  `ExtantRange.check` flags such lengths in the tables).

Pure Python over the same distributions as the viewers (`ancestors.distributions`).
Colours come from `ancestors.PALETTES`; the SVG's own greys are ones the Dark skin remaps.
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass, field
from typing import Any

from phylobot import alignment_map, ancestors, library

GAP = "-"
LINE = 60  # residues per site-map line, the usual width of a sequence listing
NO_POSTERIORS = "#b3bac5"  # PRANK and GRASP reconstructions: residues without probabilities
# Alternatives are kept down to the lowest threshold a user may pick, so the page can show
# any threshold without asking the server again.
FLOOR = round(library.MIN_THRESHOLD * 100)  # as a percentage


def percent(pp: float) -> int:
    """A probability as a whole percentage, rounded down: the page's colour class, and what
    its script compares with the threshold. Thresholds are on a 0.01 grid, so
    `percent(pp) >= 100 * t` exactly when `pp >= t`, and the server's counts and the
    slider's agree. (The epsilon: 0.29 * 100 is 28.999... in floating point.)"""
    return max(0, min(int(pp * 100 + 1e-9), 100))


def threshold_percent(threshold: float) -> int:
    return round(threshold * 100)


@dataclass(frozen=True)
class Alternative:
    state: str
    pp: float

    @property
    def pct(self) -> int:
        return percent(self.pp)


@dataclass(frozen=True)
class Tile:
    residue: int  # 1-based position in the ungapped ancestor
    column: int  # 1-based alignment column (the number in member names, K300R)
    state: str
    pp: float | None  # None for a reconstruction without posteriors
    alternatives: list[Alternative] = field(default_factory=list)  # down to FLOOR %
    shown: list[Alternative] = field(default_factory=list)  # those at the page's threshold

    @property
    def pct(self) -> int:
        return percent(self.pp) if self.pp is not None else 0

    @property
    def top(self) -> int:
        """The strongest alternative, as a percentage (0 without one): the site is
        uncertain at every threshold up to this, which is what the slider compares."""
        return max((a.pct for a in self.alternatives), default=0)

    def uncertain(self, threshold: float) -> bool:
        return self.pp is not None and self.top >= max(threshold_percent(threshold), FLOOR)

    @property
    def code(self) -> str:
        """Every alternative, compactly, for the page's script: "R35 Q25"."""
        return " ".join(f"{a.state}{a.pct:02d}" for a in self.alternatives)

    @property
    def title(self) -> str:
        """The hover read-out's fixed part (the script appends the alternatives)."""
        where = f"residue {self.residue} · column {self.column}"
        if self.pp is None:
            return f"{where} · {self.state} (no posterior probabilities)"
        return f"{where} · {self.state} {self.pp:.2f}"


@dataclass(frozen=True)
class SiteMap:
    lines: list[list[Tile]]
    length: int
    threshold: float
    mean_pp: float | None
    sequence: str  # ungapped, for copying

    @property
    def tiles(self) -> list[Tile]:
        return [t for line in self.lines for t in line]

    @property
    def ambiguous(self) -> int:
        """Sites with an alternative at or above the threshold."""
        return sum(1 for t in self.tiles if t.uncertain(self.threshold))

    def fasta(self, name: str) -> str:
        wrapped = "\n".join(self.sequence[i : i + LINE] for i in range(0, len(self.sequence), LINE))
        return f">{name}\n{wrapped}\n"


def colour_classes(palette: str) -> list[tuple[int, str]]:
    """The page's tile colours, one class per whole percentage (`q0` ... `q100`): a tile
    names its class instead of carrying a style, which keeps a page of several thousand
    residues small."""
    return [(i, ancestors.colour(i / 100, palette)) for i in range(101)]


def ml_pp(state: str, dist: list[tuple[str, float]]) -> float:
    return next((pp for s, pp in dist if s == state), 0.0)


def site_map(
    template: str,
    dists: list[list[tuple[str, float]]] | None,
    threshold: float,
    per_line: int = LINE,
) -> SiteMap:
    """The ancestor (an aligned row; gap columns are skipped) as tiles. Every state other
    than the reconstructed one whose probability reaches `threshold` is shown beneath it
    — not only the second best, since a third plausible state is exactly what ML and
    AltAll alone cannot cover.

    Alternatives are kept down to FLOOR (%), the lowest threshold a user may pick, so the page
    can move its threshold without asking the server again: `threshold` only decides
    which are shown, and counted, first."""
    if dists is not None and len(dists) != len(template):
        raise ValueError(f"columns differ: ancestor {len(template)}, posteriors {len(dists)}")
    least = max(threshold_percent(threshold), FLOOR)
    tiles: list[Tile] = []
    for column, state in enumerate(template, start=1):
        if state == GAP:
            continue
        residue = len(tiles) + 1
        if dists is None:
            tiles.append(Tile(residue, column, state, None))
            continue
        dist = dists[column - 1]
        alternatives = [Alternative(s, p) for s, p in dist if s != state and percent(p) >= FLOOR]
        shown = [a for a in alternatives if a.pct >= least]
        tiles.append(Tile(residue, column, state, ml_pp(state, dist), alternatives, shown))
    known = [t.pp for t in tiles if t.pp is not None]
    return SiteMap(
        lines=[tiles[i : i + per_line] for i in range(0, len(tiles), per_line)],
        length=len(tiles),
        threshold=threshold,
        mean_pp=sum(known) / len(known) if known else None,
        sequence="".join(t.state for t in tiles),
    )


@dataclass(frozen=True)
class Sparkline:
    width: int
    height: int
    points: str  # an SVG polyline
    ticks: list[tuple[float, int]]  # (x, strongest alternative %) per site that has one
    baseline: float  # y of pp 0

    @property
    def empty(self) -> bool:
        return not self.points

    @property
    def code(self) -> str:
        """Every tick for the page's script ("x:percentage ..."), which redraws the ticks
        and recounts the uncertain sites as the threshold moves."""
        return " ".join(f"{x:g}:{top}" for x, top in self.ticks)

    def tick_path(self, threshold_pct: int) -> str:
        """The ticks at a threshold as one SVG path (a line per tick is ten times larger)."""
        rise = self.baseline - 6
        return "".join(
            f"M{x:g} {self.baseline:g}V{rise:g}" for x, top in self.ticks if top >= threshold_pct
        )

    def count(self, threshold_pct: int) -> int:
        return sum(1 for _, top in self.ticks if top >= threshold_pct)


def sparkline(tiles: list[Tile], width: int = 220, height: int = 30) -> Sparkline:
    """The ML state's probability along the ancestor (y up = certain), as a sliding-window
    mean, so it shows *regions* of uncertainty — single sites are the ticks beneath. The
    window is about a fortieth of the ancestor (at least 5 residues), then the curve is
    sampled at most once per pixel."""
    pad = 2.0
    known = [t.pp for t in tiles if t.pp is not None]
    if not known:
        return Sparkline(width, height, "", [], height - pad)
    n = len(known)
    half = max(5, n // 40) // 2
    running = [0.0]
    for pp in known:
        running.append(running[-1] + pp)
    smooth = [
        (running[min(i + half + 1, n)] - running[max(i - half, 0)])
        / (min(i + half + 1, n) - max(i - half, 0))
        for i in range(n)
    ]
    samples = min(n, width)
    picked = [
        smooth[round(j * (n - 1) / (samples - 1))] if samples > 1 else smooth[0]
        for j in range(samples)
    ]

    def x(i: float, count: int) -> float:
        return pad + (width - 2 * pad) * (i / (count - 1) if count > 1 else 0.5)

    def y(pp: float) -> float:
        return pad + (height - 2 * pad) * (1 - pp)

    points = " ".join(f"{x(j, samples):.1f},{y(p):.1f}" for j, p in enumerate(picked))
    measured = [t for t in tiles if t.pp is not None]
    ticks = [(round(x(i, n), 1), t.top) for i, t in enumerate(measured) if t.alternatives]
    return Sparkline(width, height, points, ticks, y(0.0))


# --- Length against the extant sequences ------------------------------------------------
# An ancestor carrying every lineage's insertions (Ime2's unmasked original: 1052 residues
# against 369-800) is the failure the gap methods exist to prevent, so the tables flag any
# length outside the extant range, and `length_chart` draws them all on one histogram.


@dataclass(frozen=True)
class LengthWarning:
    beyond: int  # 1 longer than every extant sequence, -1 shorter
    text: str


@dataclass(frozen=True)
class ExtantRange:
    shortest: int
    longest: int
    count: int

    def check(self, length: int) -> LengthWarning | None:
        """A warning when `length` is outside the extant range, else None."""
        extant = f"The extant sequences are {self.shortest}–{self.longest} residues ({self.count})."
        if length > self.longest:
            return LengthWarning(1, (
                f"{length - self.longest} residues longer than the longest extant sequence. "
                f"{extant} An ancestor this long usually carries insertions from several "
                "lineages at once: use a gap version (Parsimony, ML) rather than the original "
                "reconstruction."
            ))  # fmt: skip
        if length < self.shortest:
            return LengthWarning(-1, (
                f"{self.shortest - length} residues shorter than the shortest extant sequence. "
                f"{extant} Gap placement may have removed columns this ancestor had: compare "
                "it with the other gap versions."
            ))  # fmt: skip
        return None


def extant_range(extant: list[int]) -> ExtantRange | None:
    return ExtantRange(min(extant), max(extant), len(extant)) if extant else None


# The same comparison as one figure: a single histogram of the extant lengths with every
# ancestor of the library drawn on it as a line, labelled above the plot (hover for its
# full name and length), and the extant lengths' summary in the corner. Its axis stretches
# to take every ancestor, so nothing is pinned; a far outlier squeezes the histogram,
# which is itself the message.
LC_W = 640
LC_LEFT, LC_RIGHT, LC_BOTTOM = 44, 12, 40
LC_HEADER = 26  # the title and the summary line
LC_PLOT_H = 146
LABEL_ROW = 13  # px per row of line labels, stacked upwards from the plot
LABEL_CHAR = 6.2  # px per character at the labels' 11 px: enough to keep them apart
LABEL_ROWS = 6  # past these, labels may overlap rather than the band growing further


@dataclass(frozen=True)
class ChartBar:
    x: float
    y: float
    width: float
    height: float
    title: str


@dataclass(frozen=True)
class ChartMark:
    """Every ancestor at one x (to the pixel): one line, one label, one tooltip."""

    x: float
    beyond: bool  # outside the extant range
    names: list[str]
    title: str
    label: str  # the short names, written above the line
    label_x: float  # where the label starts (it runs rightwards from its line)
    label_y: float
    top: float  # where the line starts: just under its label


@dataclass(frozen=True)
class LengthChart:
    width: int
    height: int
    left: int
    right: float
    top: float
    bottom: float
    bars: list[ChartBar]
    marks: list[ChartMark]
    x_ticks: list[tuple[float, int]]
    y_ticks: list[tuple[float, int]]
    stats: str  # "Extant (n = 7): mean 598 · median 601 · sd 140"
    count: int

    @property
    def empty(self) -> bool:
        return self.count == 0

    @property
    def center(self) -> float:
        return (self.left + self.right) / 2

    @property
    def middle(self) -> float:
        return (self.top + self.bottom) / 2


def length_stats(extant: list[int]) -> str:
    """The extant lengths' mean, median and (sample) standard deviation, as one line."""
    n = len(extant)
    sd = statistics.stdev(extant) if n > 1 else 0.0
    return (
        f"Extant (n = {n}): mean {statistics.fmean(extant):.0f} · "
        f"median {statistics.median(extant):g} · sd {sd:.0f}"
    )


def label_clashes(
    span: tuple[float, float, float],
    row: int,
    placed: list[tuple[float, float, float]],
    rows: list[int],
) -> int:
    """How many placed labels a label (x, start, end) would clash with in `row`."""
    x, start, end = span
    n = 0
    for (ox, ostart, oend), orow in zip(placed, rows, strict=False):
        if orow == row and start < oend + 4 and ostart < end + 4:
            n += 1  # labels overlap
        elif orow > row and start - 2 <= ox <= end + 2:
            n += 1  # this label sits on a line rising past it
        elif row > orow and ostart - 2 <= x <= oend + 2:
            n += 1  # this line rises through that label
    return n


def label_rows(spans: list[tuple[float, float, float]]) -> list[int]:
    """A row for each label (x, start, end), row 0 nearest the plot. Each line rises to
    its own label, so two labels clash when they overlap in one row, or when a label
    covers the line of another that rises past it. Labels start at their line and run
    right, so placing them right to left, a line never rises through a label already
    placed; only labels side by side can clash, and they stack. When no row is free,
    the least crowded is taken."""
    order = sorted(range(len(spans)), key=lambda i: -spans[i][0])
    placed: list[tuple[float, float, float]] = []
    placed_rows: list[int] = []
    rows = [0] * len(spans)
    for i in order:
        counts = [label_clashes(spans[i], row, placed, placed_rows) for row in range(LABEL_ROWS)]
        rows[i] = counts.index(min(counts))
        placed.append(spans[i])
        placed_rows.append(rows[i])
    return rows


def length_chart(extant: list[int], marks_at: list[tuple[str, str, int]]) -> LengthChart:
    """`marks_at` are the ancestors' (short name, full name, length): each becomes a line
    over the histogram, red when it is shorter or longer than every extant sequence,
    labelled above the plot with its short name and hovered for the full one."""
    right = LC_W - LC_RIGHT
    if not extant:
        top = LC_HEADER + LABEL_ROW
        bottom = top + LC_PLOT_H
        return LengthChart(LC_W, round(bottom + LC_BOTTOM), LC_LEFT, right, top, bottom,
                           [], [], [], [], "", 0)  # fmt: skip
    shortest, longest = min(extant), max(extant)
    everything = extant + [n for _, _, n in marks_at]
    low_n, high_n = min(everything), max(everything)
    pad = max((high_n - low_n) * 0.04, 5)
    low, high = low_n - pad, high_n + pad

    def x_at(length: float) -> float:
        return round(LC_LEFT + (right - LC_LEFT) * (length - low) / (high - low), 1)

    # The labels first: how many rows they take decides where the plot starts.
    grouped: dict[float, list[tuple[str, str, int]]] = {}
    for short, label, n in marks_at:
        grouped.setdefault(x_at(n), []).append((short, label, n))
    placed = []
    for x, group in sorted(grouped.items()):
        text = ", ".join(dict.fromkeys(short for short, _, _ in group))
        width = len(text) * LABEL_CHAR
        start = min(x - 2, LC_W - 2 - width)  # from the line rightwards, kept on the figure
        placed.append((x, group, text, start, start + width))
    rows = label_rows([(x, start, end) for x, _, _, start, end in placed])
    used = max(rows, default=-1) + 1
    top = LC_HEADER + max(used, 1) * LABEL_ROW + 2
    bottom = top + LC_PLOT_H

    marks = []
    for (x, group, text, start, _), row in zip(placed, rows, strict=True):
        lines = []
        for _, label, n in group:
            if n > longest:
                lines.append(f"{label}: {n} residues, {n - longest} longer than the longest extant")
            elif n < shortest:
                lines.append(
                    f"{label}: {n} residues, {shortest - n} shorter than the shortest extant"
                )
            else:
                lines.append(f"{label}: {n} residues")
        beyond = any(not shortest <= n <= longest for _, _, n in group)
        baseline = top - 4 - row * LABEL_ROW
        marks.append(ChartMark(
            x, beyond, [label for _, label, _ in group], "\n".join(lines),
            text, round(start, 1), baseline, baseline + 3,
        ))  # fmt: skip

    # Bins of a round width over the extant range: about 2 sqrt(n) of them, 10 to 30.
    wanted = min(30, max(10, round(2 * len(extant) ** 0.5)))
    size = round_width(max(longest - shortest, 1) / wanted)
    first = shortest // size * size
    counts: dict[int, int] = {}
    for n in extant:
        start_bin = first + (n - first) // size * size
        counts[start_bin] = counts.get(start_bin, 0) + 1
    tallest = max(counts.values())
    y_step = ancestors.nice_step(tallest, most=4)
    y_max = -(-tallest // y_step) * y_step

    def y_at(count: float) -> float:
        return round(bottom - (bottom - top) * count / y_max, 1)

    bars = [
        ChartBar(
            x_at(b), y_at(n), max(x_at(b + size) - x_at(b) - 1, 1), round(bottom - y_at(n), 1),
            f"{b}–{b + size - 1} residues: {n} extant sequence{'s' if n != 1 else ''}",
        )
        for b, n in sorted(counts.items())
    ]  # fmt: skip
    x_step = ancestors.nice_step(round(high - low), most=7)
    x_ticks = [
        (x_at(t), t)
        for t in range(-(-round(low) // x_step) * x_step, round(high) + 1, x_step)
        if low <= t <= high
    ]
    y_ticks = [(y_at(c), c) for c in range(0, y_max + 1, y_step)]
    return LengthChart(LC_W, round(bottom + LC_BOTTOM), LC_LEFT, right, top, bottom, bars,
                       marks, x_ticks, y_ticks, length_stats(extant), len(extant))  # fmt: skip


def round_width(width: float) -> int:
    """A round bin width (1, 2, 5, 10, 20, 50, ...) at least `width` residues."""
    step = 1
    while True:
        for factor in (1, 2, 5):
            if factor * step >= width:
                return factor * step
        step *= 10


def extant_lengths(alignments: list[dict[str, str]]) -> list[int]:
    """Each extant sequence's length, gaps removed, once per name: every alignment of a
    job holds the same proteins, so a library mixing aligners still has one histogram."""
    seen: dict[str, int] = {}
    for rows in alignments:
        for name, row in rows.items():
            seen.setdefault(name, len(row.replace(GAP, "")))
    return list(seen.values())


# --- Layouts: which frame each target is drawn on ---------------------------------------


@dataclass(frozen=True)
class Entry:
    """A target with what laying it out needs: its alignment's extant rows (same columns
    as `template`) and which alignment that is."""

    label: str
    anchor: str
    title: str
    template: str
    dists: list[list[tuple[str, float]]] | None
    aligner: str
    aligner_label: str
    alignment: dict[str, str]
    # A library member drawn on its ancestor's columns (L5c): the columns where it differs
    # from that ancestor's ML sequence, outlined instead of the uncertain sites; `cut` when
    # the budget left it out, which dims the row.
    member: bool = False
    marked: frozenset[int] = frozenset()  # 0-based columns of `template`
    cut: bool = False

    @property
    def fits(self) -> bool:
        """Its alignment really has the template's columns (always, in practice)."""
        return bool(self.alignment) and all(
            len(row) == len(self.template) for row in self.alignment.values()
        )


# Colour schemes the library page offers. Posterior is its own (the ancestors' certainty
# is the thing this page is about); the rest colour a residue on its own, so they can be
# drawn both in the viewer frames and in the site maps beneath each ancestor. Clustal X,
# percentage identity and differences are left out on purpose: they judge a residue by its
# column, which a site map has no notion of, and a scheme that coloured only half the page
# would be worse than not offering it.
SCHEMES = {
    "posterior": "Posterior probability",
    "chemistry": "Chemistry",
    "hydropathy": "Hydropathy (Kyte-Doolittle)",
    "structure": "Secondary-structure propensity (Chou-Fasman)",
    "none": "None",
}
DEFAULT_SCHEME = "posterior"
BY_POSTERIOR = "posterior"  # the one scheme the palette applies to


def parse_scheme(text: object) -> str:
    key = str(text)
    return key if key in SCHEMES else DEFAULT_SCHEME


OWN = "own"  # each alignment on its own columns
ON_SEQUENCE = "sequence"  # no alignment: one extant sequence's residues


@dataclass(frozen=True)
class Layout:
    """What the library's viewer is laid out on: an alignment (OWN, an aligner, or
    ON_SEQUENCE) and the extant sequences — shown beside the ancestors on an alignment,
    or the one sequence that is the frame on ON_SEQUENCE."""

    alignment: str
    extant: tuple[str, ...] = ()


def aligners(entries: list[Entry]) -> dict[str, str]:
    """{aligner: label}, in the order the targets first use them."""
    return {e.aligner: e.aligner_label for e in entries}


def alignment_options(entries: list[Entry]) -> list[tuple[str, str]]:
    used = aligners(entries)
    options = [(OWN, "Each alignment separately")] if len(used) > 1 else []
    options += [(key, f"{label}") for key, label in used.items()]
    return options + [(ON_SEQUENCE, "Single sequence")]


def extant_options(entries: list[Entry]) -> list[str]:
    return next((list(e.alignment) for e in entries if e.fits), [])


def resolve_layout(entries: list[Entry], alignment: str, extant: list[str]) -> Layout:
    """The layout asked for, with anything unknown replaced by the defaults: each
    alignment separately when there are several (never one frame by default), and on a
    sequence, exactly one sequence (the first, unless one is chosen)."""
    keys = {key for key, _ in alignment_options(entries)}
    used = list(aligners(entries))
    if alignment not in keys:
        alignment = OWN if len(used) > 1 else used[0] if used else OWN
    names = extant_options(entries)
    chosen = tuple(dict.fromkeys(n for n in extant if n in names))
    if alignment == ON_SEQUENCE:
        chosen = chosen[:1] or tuple(names[:1])
    return Layout(alignment, chosen)


@dataclass(frozen=True)
class Row:
    """One row of a frame, one entry per frame unit (column or residue)."""

    name: str
    ancestor: bool
    states: str  # "-" where the row has nothing in that unit
    pp: list[int]  # whole percentage of the shown state, -1 for none
    top: list[int]  # the strongest alternative's percentage, 0 for none
    alts: list[str]  # every alternative down to FLOOR, "R35 Q25"
    faded: list[int] = field(default_factory=list)  # units whose column only partly maps
    anchor: str = ""  # the table row an ancestor (or a member) opens
    member: bool = False
    marked: list[int] = field(default_factory=list)  # units outlined: a member's changes
    dim: bool = False  # a member the budget left out


@dataclass(frozen=True)
class Frame:
    """One alignment-viewer instance on the library page."""

    title: str
    note: str
    unit: str  # "column" or "residue"
    units: int
    rows: list[Row]
    tracks: list[dict[str, Any]] = field(default_factory=list)
    index: int = 0  # its place on the page
    group: str = ""  # "" for the ancestors' viewers, "m" for the members'

    @property
    def json_id(self) -> str:
        """The id of its data's script element (alignment.js finds "aln-data..." ids)."""
        return f"aln-data-{self.group}{self.index}"

    @property
    def members(self) -> int:
        return sum(1 for r in self.rows if r.member)

    def data(
        self,
        palette: str,
        threshold_pct: int,
        cell: float | str = "fit",
        scheme: str = DEFAULT_SCHEME,
    ) -> dict[str, Any]:
        """The alignment viewer's data (alignment.js), with the posterior scheme's extras."""
        return {
            "names": [r.name for r in self.rows],
            "ancestor": [r.ancestor for r in self.rows],
            "anchors": [r.anchor for r in self.rows],
            "rows": [r.states for r in self.rows],
            "amino": [],
            "codon": False,
            "extant": sum(1 for r in self.rows if not r.ancestor),
            "reference": -1,
            "scheme": scheme,
            "cell": cell,
            "row": 16,
            "height": len(self.rows),
            "tracks": self.tracks,
            "ancestorViewer": None,
            "rowEvents": True,
            "pp": [r.pp for r in self.rows],
            "top": [r.top for r in self.rows],
            "alts": [r.alts for r in self.rows],
            "faded": [r.faded for r in self.rows],
            "member": [r.member for r in self.rows],
            "marked": [r.marked for r in self.rows],
            "dim": [r.dim for r in self.rows],
            "palette": palette_stops(palette),
            "threshold": threshold_pct,
        }


def palette_stops(palette: str) -> list[list[Any]]:
    """A palette as the page's script interpolates it: [[position, [r, g, b]], ...]."""
    stops = ancestors.PALETTES.get(palette, ancestors.PALETTES[ancestors.DEFAULT_PALETTE])[1]
    return [[position, list(rgb)] for position, rgb in stops]


Cellinfo = tuple[str, int, int, str]  # state, pp %, top %, alternatives code


def column_info(template: str, dists: list[list[tuple[str, float]]] | None) -> list[Cellinfo]:
    """Per alignment column of an ancestor: its state, the state's probability, its
    strongest alternative and all alternatives down to FLOOR (percentages)."""
    out: list[Cellinfo] = []
    for column, state in enumerate(template):
        if state == GAP or dists is None:
            out.append((state, -1, 0, ""))
            continue
        alternatives = [
            Alternative(s, p) for s, p in dists[column] if s != state and percent(p) >= FLOOR
        ]
        top = max((a.pct for a in alternatives), default=0)
        code = " ".join(f"{a.state}{a.pct:02d}" for a in alternatives)
        out.append((state, percent(ml_pp(state, dists[column])), top, code))
    return out


def _row(
    e: Entry, name: str, cells: list[Cellinfo | None], faded: list[int], marked: list[int]
) -> Row:
    """An ancestor's or a member's row. A member's outline is its changes from the ML
    sequence, so its uncertain sites (`top`) are not outlined as well."""
    filled = [c if c is not None else (GAP, -1, 0, "") for c in cells]
    top = [0] * len(filled) if e.member else [c[2] for c in filled]
    return Row(
        name, True, "".join(c[0] for c in filled), [c[1] for c in filled],
        top, [c[3] for c in filled], faded, e.anchor, e.member, marked, e.cut,
    )  # fmt: skip


def _extant(name: str, states: str) -> Row:
    n = len(states)
    return Row(name, False, states, [-1] * n, [0] * n, [""] * n)


def _named(e: Entry, entries: list[Entry]) -> str:
    """A row's name, naming its alignment when the frame holds several combinations (the
    same node — and its members — can come from each, so the aligner tells them apart)."""
    several = len({x.aligner for x in entries}) > 1
    named = several and e.aligner_label
    return f"{e.aligner_label} {e.label}" if named else e.label


def agreement_track(values: list[float], unit: str) -> dict[str, Any]:
    return {
        "key": "agreement",
        "label": "Aligners agree",
        "description": f"Per {unit}: 1 where every alignment builds the same column.",
        "series": [{"label": "agreement", "values": [round(v, 3) for v in values]}],
        "low": 0,
        "high": 1,
        "text": "",
    }


def disagreement_track(values: list[float], unit: str) -> dict[str, Any]:
    """Columns where an ancestor's combinations reconstruct different ML residues
    (phylobot.library_disagreement): the sites the Disputed column counts, on the sequence."""
    return {
        "key": "disagreement",
        "label": "Combinations disagree",
        "description": f"Per {unit}: marked where your combinations give this ancestor a "
        "different ML residue (the Disputed column, on the sequence).",
        "series": [{"label": "disagreement", "values": [round(v, 3) for v in values]}],
        "low": 0,
        "high": 1,
        "text": "",
    }


def frames(entries: list[Entry], layout: Layout) -> list[Frame]:
    """The library's viewers for a layout. Targets that cannot be placed on the chosen
    frame — another sequence set, or a sequence their alignment lacks — keep their own
    columns in frames of their own: never squeezed onto a frame they do not share."""
    if layout.alignment == ON_SEQUENCE and layout.extant:
        return _on_sequence(entries, layout.extant[0])
    if layout.alignment not in (OWN, ON_SEQUENCE):
        reference = next((e for e in entries if e.aligner == layout.alignment and e.fits), None)
        if reference is not None:
            return _on_alignment(entries, reference, layout.extant)
    return _own(entries, layout.extant)


def _own(entries: list[Entry], extant: tuple[str, ...], note: str = "") -> list[Frame]:
    out = []
    for key, label in aligners(entries).items():
        mine = [e for e in entries if e.aligner == key]
        width = len(mine[0].template)
        rows = [_extant(n, mine[0].alignment[n]) for n in extant if n in mine[0].alignment]
        rows += [
            _row(e, e.label, list(column_info(e.template, e.dists)), [], sorted(e.marked))
            for e in mine
        ]
        out.append(Frame(f"{label} alignment · {width} columns", note, "column", width, rows))
    return out


def _on_alignment(entries: list[Entry], reference: Entry, extant: tuple[str, ...]) -> list[Frame]:
    frame = reference.alignment
    width = len(reference.template)
    placed = [e for e in entries if e.fits and alignment_map.same_sequences(e.alignment, frame)]
    rows = [_extant(n, frame[n]) for n in extant if n in frame]
    maps: dict[str, alignment_map.ColumnMap] = {}
    hidden = 0
    for e in placed:
        info = column_info(e.template, e.dists)
        if e.aligner == reference.aligner:
            rows.append(_row(e, _named(e, placed), list(info), [], sorted(e.marked)))
            continue
        if e.aligner not in maps:
            maps[e.aligner] = alignment_map.map_columns(e.alignment, frame)
        mapped = maps[e.aligner]
        # Each frame column shows the ancestor's residue from the column mapped there that
        # agrees best (the first such, on a tie); residues sharing a frame column with a
        # better one, or in a column holding no extant residue, are not shown.
        best: dict[int, int] = {}
        for column, state in enumerate(e.template):
            if state == GAP:
                continue
            if column not in mapped.target:  # no extant residue to place it by
                hidden += 1
                continue
            m = mapped.target[column]
            if m not in best:
                best[m] = column
                continue
            hidden += 1  # two of its residues meet in one frame column: one is not shown
            if mapped.agreement[column] > mapped.agreement[best[m]]:
                best[m] = column
        cells: list[Cellinfo | None] = [info[best[m]] if m in best else None for m in range(width)]
        faded = [m for m, column in sorted(best.items()) if not mapped.exact(column)]
        marked = [m for m, column in sorted(best.items()) if column in e.marked]
        rows.append(_row(e, _named(e, placed), cells, faded, marked))
    tracks: list[dict[str, Any]] = []
    note = ""
    if maps:
        agreement = [
            min(m.reference_agreement.get(column, 0.0) for m in maps.values())
            for column in range(width)
        ]
        tracks.append(agreement_track(agreement, "column"))
        others = ", ".join(dict.fromkeys(e.aligner_label for e in placed if e.aligner in maps))
        note = (
            f"{others} ancestors are placed by their extant residues, faded where the "
            "alignments disagree about which residues belong together."
        )
        if hidden:
            note += f" {hidden} of their residues share a column with another and are not shown."
    title = f"On {reference.aligner_label} columns · {width}"
    out = [Frame(title, note, "column", width, rows, tracks)]
    rest = [e for e in entries if e not in placed]
    return out + _own(rest, extant, "Not the same sequences: shown on their own columns.")


def _on_sequence(entries: list[Entry], name: str) -> list[Frame]:
    placed = [e for e in entries if e.fits and name in e.alignment]
    if not placed:
        return _own(entries, ())
    length = len(placed[0].alignment[name].replace(GAP, ""))
    rows = [_extant(name, placed[0].alignment[name].replace(GAP, ""))]
    inserted = [0] * length
    alignments: dict[str, dict[str, str]] = {}
    for e in placed:
        alignments.setdefault(e.aligner, e.alignment)
        frame = alignment_map.sequence_frame(e.alignment[name])
        info = column_info(e.template, e.dists)
        cells: list[Cellinfo | None] = [None] * length
        extra = [0] * (length + 1)
        marked = []
        for column, residue in enumerate(frame.residue):
            if residue is not None:
                cells[residue - 1] = info[column]
                if column in e.marked:
                    marked.append(residue - 1)
            elif e.template[column] != GAP:
                extra[frame.slot[column]] += 1
        for k in range(length):  # residues inserted after residue k + 1 (and before the first)
            inserted[k] = max(inserted[k], extra[k + 1] + (extra[0] if k == 0 else 0))
        rows.append(_row(e, _named(e, placed), cells, [], marked))
    tracks: list[dict[str, Any]] = [
        {
            "key": "inserted",
            "label": "Inserted residues",
            "description": (
                f"Per residue: the most residues any ancestor has after it where {name} has none."
            ),
            "series": [{"label": "inserted", "values": inserted}],
            "low": 0,
            "high": None,
            "text": "",
        }
    ]
    if len(alignments) > 1:
        tracks.insert(
            0,
            agreement_track(
                alignment_map.residue_agreement(list(alignments.values()), name), "residue"
            ),
        )
    note = (
        f"Each ancestor on {name}'s residues, through its own alignment. Residues an ancestor "
        f"has where {name} has none are not drawn; the track counts them."
    )
    out = [Frame(f"On {name} · {length} residues", note, "residue", length, rows, tracks)]
    rest = [e for e in entries if e not in placed]
    return out + _own(
        rest, (), f"Their alignments do not contain {name}: shown on their own columns."
    )


MIN_CELL, MAX_CELL = 0.05, 20.0  # the zoom slider's ends, in pixels per column


def parse_zoom(text: object) -> float | str:
    """The viewers' zoom from request text: "fit" (the whole alignment in view, the
    default) or a cell width in pixels, clamped to the slider's range."""
    if str(text) == "fit":
        return "fit"
    try:
        return round(min(max(float(str(text)), MIN_CELL), MAX_CELL), 2)
    except ValueError:
        return "fit"


# The Design section's coverage chart: wells across, the share of plausible states those
# wells carry up the side, one line per strategy. Pure geometry, as the tree viewer's is:
# the template draws it, so the page needs no chart library.
CHART_W, CHART_H = 440, 150
CHART_PAD = 34  # room for the axis labels on the left and under


@dataclass(frozen=True)
class Line:
    """One strategy's curve: its points, and what it reaches at the chosen budget."""

    key: str
    label: str
    points: str
    at_budget: float  # the share covered at the budget, 0-1
    covered: int

    @property
    def percent(self) -> float:
        """The share as a percentage, to one decimal: at a useful budget the strategies
        often differ by a fraction of a point, which a whole number would hide."""
        return round(self.at_budget * 100, 1)


@dataclass(frozen=True)
class CoverageChart:
    """Coverage against budget, for every strategy at once."""

    width: int
    height: int
    pad: int
    lines: list[Line]
    total: int
    budget: int
    max_wells: int
    budget_x: float
    ticks: list[tuple[float, int]]  # x, wells
    grid: list[tuple[float, str]]  # y, label ("50%")
    middle_x: float = 0  # the plot's centre, where the axis is named

    @property
    def empty(self) -> bool:
        return self.total == 0 or self.max_wells < 1

    @property
    def tied(self) -> bool:
        """Whether every strategy reaches the same states at this budget — common once
        the budget is generous, and worth saying rather than leaving three lines on top
        of one another to be puzzled over."""
        return len({line.covered for line in self.lines}) == 1 and len(self.lines) > 1


def coverage_chart(
    curves: dict[str, tuple[str, list[int]]], total: int, budget: int
) -> CoverageChart:
    """`curves` maps a strategy key to its label and its coverage curve (index = wells).

    Every curve is drawn over the same axes, so the three strategies can be read against
    each other, and the chosen budget is a line down the chart: where it crosses each
    curve is what that choice would buy.
    """
    max_wells = max((len(c) - 1 for _, c in curves.values()), default=0)
    inner_w = CHART_W - CHART_PAD - 6
    inner_h = CHART_H - CHART_PAD - 6
    if total <= 0 or max_wells < 1:
        return CoverageChart(CHART_W, CHART_H, CHART_PAD, [], total, budget, max_wells, 0, [], [])

    def x_at(wells: float) -> float:
        return CHART_PAD + inner_w * (wells / max_wells)

    def y_at(share: float) -> float:
        return 6 + inner_h * (1 - min(max(share, 0.0), 1.0))

    lines = []
    for key, (label, curve) in curves.items():
        points = " ".join(f"{x_at(i):.1f},{y_at(v / total):.1f}" for i, v in enumerate(curve))
        at = curve[min(budget, len(curve) - 1)]
        lines.append(Line(key, label, points, at / total, at))
    step = max(1, round(max_wells / 4))
    ticks = [(x_at(w), w) for w in range(0, max_wells + 1, step)]
    grid = [(y_at(share), f"{round(share * 100)}%") for share in (0.0, 0.5, 1.0)]
    return CoverageChart(
        CHART_W, CHART_H, CHART_PAD, lines, total, budget, max_wells,
        x_at(min(budget, max_wells)), ticks, grid, x_at(max_wells / 2),
    )  # fmt: skip
