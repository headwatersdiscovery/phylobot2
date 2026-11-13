"""The ancestor viewer: one ancestral node's sites and posterior distributions.

As in the 2016 PhyloBot (Hanson-Smith & Johnson, Fig. 5C): one row per alignment
column, giving the column, the residue in each extant sequence the user
picks, and the ancestral states as tiles, most probable first, coloured by posterior
probability (cool = well supported, warm = uncertain).

The residues and probabilities are IQ-TREE's (ancestors_posteriors.tsv). A gap-masked
version (parsimony, fitch, ml) shows a gap tile where it masked the column and IQ-TREE's
distribution elsewhere. PRANK and GRASP reconstruct ancestors themselves and report no
posterior probabilities, so their versions show their own residue, uncoloured.
"""

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from phylobot.gaps import MASKING
from phylobot.trees import load_tree

GAP = "-"
MIN_SHOWN = 0.01  # alternative states below this probability get no tile (the TSV has all)
ORIGINAL = "original"

# Gap versions: key -> (menu label, file). "original" is IQ-TREE's unmasked reconstruction.
VERSIONS = {
    ORIGINAL: ("None (IQ-TREE's reconstruction, a residue at every column)", "ancestors.fasta"),
    "parsimony": ("Parsimony (shared presence)", "ancestors.gaps-parsimony.fasta"),
    "fitch": ("Parsimony (Fitch)", "ancestors.gaps-fitch.fasta"),
    "ml": ("Maximum likelihood", "ancestors.gaps-ml.fasta"),
    "prank": ("PRANK's reconstruction (no posteriors)", "ancestors.gaps-prank.fasta"),
    "grasp": ("GRASP's reconstruction (no posteriors)", "ancestors.gaps-grasp.fasta"),
}
# The same versions named for a table cell, where the menu's explanatory label would wrap
# a row to three lines. The full label stays as the cell's tooltip.
SHORT_VERSIONS = {
    ORIGINAL: "None",
    "parsimony": "Parsimony (shared)",
    "fitch": "Parsimony (Fitch)",
    "ml": "Maximum likelihood",
    "prank": "PRANK",
    "grasp": "GRASP",
}
# PAML codeml's reconstructions (coding DNA): "codeml-<model>", files in codeml/<model>/,
# with codon posteriors summed into amino acids like IQ-TREE's codon reconstructions.
CODEML_PREFIX = "codeml-"


def codeml_model(version: str) -> str | None:
    return version[len(CODEML_PREFIX) :] if version.startswith(CODEML_PREFIX) else None


def version_label(version: str) -> str:
    if model := codeml_model(version):
        return f"PAML codeml {model} (its own reconstruction)"
    return VERSIONS[version][0]


def short_version_label(version: str) -> str:
    """The version named in a few words, for a table cell."""
    if model := codeml_model(version):
        return f"codeml {model}"
    return SHORT_VERSIONS[version]


def version_folder(version: str) -> str:
    """Where a version's ancestor files are, inside its combination's folder."""
    model = codeml_model(version)
    return f"codeml/{model}/" if model else ""


CODON_POSTERIORS = "ancestors_codon_posteriors.tsv"  # codon reconstructions only


def coding_file(version: str) -> str:
    """The version's coding ancestors (codon reconstructions), inside its combination."""
    if model := codeml_model(version):
        return f"codeml/{model}/ancestors.cds.fasta"
    if version == ORIGINAL:
        return "ancestors.cds.fasta"
    return f"ancestors.gaps-{version}.cds.fasta"


def codon_distributions(tsv: str, node: str) -> list[list[tuple[str, float]]]:
    """Per alignment column, `node`'s codons and probabilities (those kept: 0.01 or
    more), most probable first, from ancestors_codon_posteriors.tsv."""
    _header, *rows = posterior_lines(tsv, node)
    by_site: dict[int, list[tuple[str, float]]] = {}
    for row in rows:
        _node, site, codon, p = row.split("\t")
        by_site.setdefault(int(site), []).append((codon, float(p)))
    return [sorted(by_site[s], key=lambda cp: (-cp[1], cp[0])) for s in sorted(by_site)]


def has_posteriors(version: str) -> bool:
    """IQ-TREE's reconstruction and its masked copies, and codeml's, have posteriors;
    PRANK's and GRASP's don't."""
    return version in (ORIGINAL, *MASKING) or codeml_model(version) is not None


# Heatmap palettes: posterior probability -> colour, by linear interpolation between
# stops. Each runs from warm/uncertain (pp 0) to cool/certain (pp 1) (greys: light to
# dark); stops are packed
# near the top, where most sites and all the decisions are.
Stop = tuple[float, tuple[int, int, int]]
PALETTES: dict[str, tuple[str, list[Stop]]] = {
    "coolwarm": ("Red → blue", [
        (0.0, (215, 48, 39)), (0.4, (252, 141, 89)), (0.6, (254, 224, 144)),
        (0.8, (145, 191, 219)), (1.0, (69, 117, 180)),
    ]),
    "spectral": ("Spectral", [
        (0.0, (213, 62, 79)), (0.4, (252, 141, 89)), (0.6, (254, 224, 139)),
        (0.8, (153, 213, 148)), (1.0, (50, 136, 189)),
    ]),
    "viridis": ("Viridis (colour-blind safe)", [
        (0.0, (253, 231, 37)), (0.4, (94, 201, 98)), (0.6, (33, 145, 140)),
        (0.8, (59, 82, 139)), (1.0, (68, 1, 84)),
    ]),
    "greys": ("Greys (dark = certain)", [(0.0, (245, 245, 245)), (1.0, (70, 70, 70))]),
}  # fmt: skip
DEFAULT_PALETTE = "coolwarm"


def parse_fasta(text: str) -> dict[str, str]:
    """{name: sequence} from FASTA text; names end at the first whitespace."""
    records: dict[str, str] = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0] if line[1:].strip() else ""
            records[name] = ""
        elif name is not None:
            records[name] += line.strip()
    return records


def gap_probabilities(tsv: str, node: str) -> list[float]:
    """`node`'s P(gap) per column from gaps-ml.tsv (node, site, p_gap)."""
    rows = [line.split("\t") for line in posterior_lines(tsv, node)[1:]]
    return [float(p) for _, _, p in sorted(rows, key=lambda r: int(r[1]))]


def posterior_lines(tsv: str, node: str) -> list[str]:
    """The header and `node`'s rows of ancestors_posteriors.tsv (or ancestors_pp.tsv)."""
    lines = tsv.splitlines()
    return lines[:1] + [line for line in lines[1:] if line.split("\t", 1)[0] == node]


def distributions(tsv: str, node: str) -> list[list[tuple[str, float]]]:
    """Per alignment column, `node`'s states and probabilities, most probable first."""
    header, *rows = posterior_lines(tsv, node)
    states = header.split("\t")[2:]
    by_site: dict[int, list[tuple[str, float]]] = {}
    for row in rows:
        fields = row.split("\t")
        pairs = [(s, float(p)) for s, p in zip(states, fields[2:], strict=True)]
        by_site[int(fields[1])] = sorted(pairs, key=lambda sp: (-sp[1], sp[0]))
    return [by_site[site] for site in sorted(by_site)]


def rgb(pp: float, palette: str) -> tuple[int, int, int]:
    """The palette's colour for a posterior probability."""
    stops = PALETTES.get(palette, PALETTES[DEFAULT_PALETTE])[1]
    pp = min(max(pp, 0.0), 1.0)
    i = next(i for i in range(1, len(stops)) if pp <= stops[i][0])
    (low, a), (high, b) = stops[i - 1], stops[i]
    t = (pp - low) / (high - low)
    r, g, b_ = (round(x + (y - x) * t) for x, y in zip(a, b, strict=True))
    return r, g, b_


def text_colour(background: tuple[int, int, int]) -> str:
    """Dark or light text, whichever reads better on `background`."""
    luminance = (0.299 * background[0] + 0.587 * background[1] + 0.114 * background[2]) / 255
    return "#111111" if luminance > 0.55 else "#ffffff"


def colour(pp: float, palette: str) -> str:
    """Inline style for a tile: background from the palette, text dark or light to suit."""
    rgb_ = rgb(pp, palette)
    return f"background-color: rgb{rgb_}; color: {text_colour(rgb_)}"


@dataclass(frozen=True)
class Tile:
    state: str
    pp: float | None  # None: no posterior (a gap, or a PRANK/GRASP residue)
    style: str = ""
    kind: str = "state"  # "state", "gap", "external"
    note: str = ""  # codon tiles: the amino acid the codon encodes


@dataclass(frozen=True)
class SiteRow:
    column: int  # 1-based alignment column
    extants: list[tuple[str, int | None]]  # per chosen extant sequence: its residue
    # there ("–" at a gap) and that residue's 1-based number (None at a gap)
    tiles: list[Tile]
    temperature: float | None = None  # calibrated views: the column's T


def residue_numbers(aligned: str) -> list[int | None]:
    """Per alignment column, the 1-based residue number in this row (None at a gap)."""
    numbers: list[int | None] = []
    count = 0
    for char in aligned:
        if char == GAP:
            numbers.append(None)
        else:
            count += 1
            numbers.append(count)
    return numbers


def site_rows(
    version: str,
    ancestor: str,
    dists: list[list[tuple[str, float]]],
    extants: list[str],
    palette: str = DEFAULT_PALETTE,
    p_gap: list[float] | None = None,
    notes: dict[str, str] | None = None,
) -> list[SiteRow]:
    """The viewer's rows. `ancestor` is the node's row in the chosen version's FASTA;
    `dists` its IQ-TREE distributions (per column); `extants` the chosen sequences' rows.

    `p_gap` (ml gap placement: the node's P(gap) per column, from gaps-ml.tsv) adds the
    gap as a scored tile: first where the column is masked, otherwise after the residues
    when P(gap) >= MIN_SHOWN. Residue probabilities are IQ-TREE's, given a residue."""
    lengths = {len(ancestor), len(dists), *(len(row) for row in extants)}
    if len(lengths) != 1:
        raise ValueError(
            f"columns differ: ancestor {len(ancestor)}, posteriors {len(dists)}, "
            f"extant {[len(row) for row in extants]}"
        )
    if p_gap is not None and len(p_gap) != len(dists):
        raise ValueError(f"columns differ: gap probabilities {len(p_gap)}, posteriors {len(dists)}")
    numbers = [residue_numbers(row) for row in extants]
    rows = []
    for column, (char, dist) in enumerate(zip(ancestor, dists, strict=True)):
        gap = p_gap[column] if p_gap is not None else None
        if char == GAP:
            tiles = [Tile("–", gap, colour(gap, palette) if gap is not None else "", "gap")]
        elif not has_posteriors(version):
            tiles = [Tile(char.upper(), None, kind="external")]
        else:
            shown = [(s, p) for i, (s, p) in enumerate(dist) if i == 0 or p >= MIN_SHOWN]
            notes = notes or {}
            kind = "codon" if notes else "state"
            tiles = [
                Tile(s, p, colour(p, palette), kind, notes.get(s, "")) if notes
                else Tile(s, p, colour(p, palette), kind)
                for s, p in shown
            ]  # fmt: skip
            if gap is not None and gap >= MIN_SHOWN:
                tiles.append(Tile("–", gap, colour(gap, palette), "gap"))
        shown_extants = [
            (row[column].upper() if row[column] != GAP else "–", n[column])
            for row, n in zip(extants, numbers, strict=True)
        ]
        rows.append(SiteRow(column + 1, shown_extants, tiles))
    return rows


# Plots of a node's support (the 2016 paper's Fig. 5B, plus a histogram). SVG user units;
# the SVGs scale uniformly with the page, so text stays in proportion.
TRACE_SIZE, HIST_SIZE = (1000, 260), (420, 260)
# Margins (left, right, top, bottom) hold the tick labels and axis titles; the trace's
# are larger because its text is set larger (it is drawn at a smaller scale).
HIST_MARGIN, TRACE_MARGIN = (48, 12, 10, 42), (66, 14, 12, 58)
Y_TICKS = (0.0, 0.25, 0.5, 0.75, 1.0)
TRACE_RANKS = ("ML residue", "2nd best", "3rd best")

# Line styles for the trace, for readers who need more contrast, thicker lines, or no
# reliance on colour: key -> (menu label, [(colour, width, dash pattern)] per rank).
LineStyle = tuple[str, float, str]
LINE_STYLES: dict[str, tuple[str, list[LineStyle]]] = {
    "standard": ("Blue, orange, grey", [
        ("#1f4e8c", 1.8, ""), ("#e08214", 1.1, ""), ("#8c8c8c", 0.9, ""),
    ]),
    "okabe-ito": ("Colour-blind safe (Okabe–Ito)", [
        ("#0072b2", 2.0, ""), ("#e69f00", 1.4, ""), ("#009e73", 1.2, ""),
    ]),
    "dashes": ("Black, dashed, dotted", [
        ("#111111", 1.8, ""), ("#444444", 1.4, "6 4"), ("#777777", 1.4, "1.5 3"),
    ]),
    "bold": ("Thick, high contrast", [
        ("#000000", 3.2, ""), ("#d55e00", 2.4, "8 4"), ("#0072b2", 2.2, "2 3"),
    ]),
    "greys": ("Greys", [("#111111", 2.0, ""), ("#6e6e6e", 1.4, ""), ("#b0b0b0", 1.2, "")]),
}  # fmt: skip
DEFAULT_LINES = "standard"


@dataclass(frozen=True)
class Trace:
    rank: str
    points: str  # SVG polyline points
    colour: str
    width: float
    dash: str  # SVG stroke-dasharray ("" = solid)


@dataclass(frozen=True)
class Frame:
    """A plot's drawing area inside its viewBox."""

    width: int
    height: int
    left: float
    top: float
    plot_width: float
    plot_height: float

    @classmethod
    def of(cls, size: tuple[int, int], margin: tuple[int, int, int, int]) -> Frame:
        (w, h), (left, right, top, bottom) = size, margin
        return cls(w, h, left, top, w - left - right, h - top - bottom)

    @property
    def right(self) -> float:
        return self.left + self.plot_width

    @property
    def bottom(self) -> float:
        return self.top + self.plot_height

    @property
    def middle(self) -> float:
        return self.top + self.plot_height / 2

    @property
    def center(self) -> float:
        return self.left + self.plot_width / 2

    def y(self, p: float) -> float:
        return self.top + (1 - p) * self.plot_height

    @property
    def y_ticks(self) -> list[tuple[float, str]]:
        return [(round(self.y(t), 1), f"{t:g}") for t in Y_TICKS]


@dataclass(frozen=True)
class Bar:
    x: float
    y: float
    width: float
    height: float
    proportion: float
    label: str  # the bin, e.g. "0.9-1.0"
    fill: str  # the palette's colour at the bin's middle


@dataclass(frozen=True)
class NodePlots:
    sites: int
    mean: float  # of the ML residue's posterior probability, over the ancestor's sites
    sd: float
    histogram: list[Bar]
    hist_frame: Frame
    hist_ticks: list[tuple[float, str]]  # x positions of bin edges 0, 0.2, ... 1
    traces: list[Trace]  # ML residue first
    trace_frame: Frame
    trace_ticks: list[tuple[float, int]]  # (x, site number)
    tooltip: list[list[object]]  # per site: number, column, then state/pp pairs (top 3)


def nice_step(n: int, most: int = 10) -> int:
    """A 1/2/5 x 10^k step giving at most `most` ticks over 1..n."""
    magnitude = 1
    while True:
        for factor in (1, 2, 5):
            if n / (factor * magnitude) <= most:
                return factor * magnitude
        magnitude *= 10


def node_plots(
    version: str,
    ancestor: str,
    dists: list[list[tuple[str, float]]],
    palette: str,
    lines: str = DEFAULT_LINES,
    unit: str = "residue",  # "codon" when the viewer shows codons
) -> NodePlots | None:
    """Support plots over the ancestor's sites (its non-gap columns in this version).
    None for versions without posteriors (PRANK, GRASP) or an ancestor with no sites."""
    if not has_posteriors(version):
        return None
    columns = [c for c, char in enumerate(ancestor) if char != GAP]
    if not columns:
        return None
    top = [(dists[c] + [("", 0.0)] * 3)[:3] for c in columns]
    ml = [site[0][1] for site in top]
    n = len(ml)
    mean = sum(ml) / n
    sd = (sum((p - mean) ** 2 for p in ml) / n) ** 0.5

    hist = Frame.of(HIST_SIZE, HIST_MARGIN)
    counts = [0] * 10
    for p in ml:
        counts[min(int(p * 10), 9)] += 1
    slot = hist.plot_width / 10
    bars = []
    for i, count in enumerate(counts):
        share = count / n
        height = share * hist.plot_height
        bars.append(Bar(
            x=round(hist.left + i * slot + 1, 1), y=round(hist.bottom - height, 1),
            width=round(slot - 2, 1), height=round(height, 1), proportion=share,
            label=f"{i / 10:.1f}–{(i + 1) / 10:.1f}", fill=f"rgb{rgb((i + 0.5) / 10, palette)}",
        ))  # fmt: skip
    hist_ticks = [(round(hist.left + k * 2 * slot, 1), f"{k * 0.2:.1f}") for k in range(6)]

    trace = Frame.of(TRACE_SIZE, TRACE_MARGIN)

    def x(i: int) -> float:
        return trace.left + (i / (n - 1) if n > 1 else 0.5) * trace.plot_width

    styles = LINE_STYLES.get(lines, LINE_STYLES[DEFAULT_LINES])[1]
    traces = [
        Trace(
            rank,
            " ".join(f"{x(i):.1f},{trace.y(site[r][1]):.1f}" for i, site in enumerate(top)),
            *styles[r],
        )
        for r, rank in enumerate((f"ML {unit}", *TRACE_RANKS[1:]))
    ]
    step = nice_step(n)
    numbers = [1, *range(step, n + 1, step)]
    trace_ticks = [(round(x(k - 1), 1), k) for k in numbers]
    tooltip: list[list[object]] = [
        [i + 1, c + 1, *(v for state, pp in site for v in (state, round(pp, 3)))]
        for i, (c, site) in enumerate(zip(columns, top, strict=True))
    ]
    return NodePlots(n, mean, sd, bars, hist, hist_ticks, traces, trace, trace_ticks, tooltip)


@dataclass(frozen=True)
class AltAllPlot:
    """How many sites AltAll would switch, by the probability of the alternative.

    Deliberately a histogram of the **next-most-probable state's** probability, not of
    the ML residue's: the AltAll threshold is a cut on the alternative, so drawing it on
    the ML histogram beside it would put a line at a value that histogram does not plot.
    Sites to the right of the threshold are the ones AltAll switches, which is the
    question the slider asks.
    """

    bars: list[Bar]
    frame: Frame
    ticks: list[tuple[float, str]]  # x positions of bin edges 0, 0.2, ... 1
    sites: int  # sites with any alternative at all


def altall_plot(pps: Sequence[float], palette: str) -> AltAllPlot | None:
    """The alternatives' probabilities as ten bins over 0-1, on the same frame as the
    ML histogram so the two read as a pair. None when nothing has an alternative."""
    if not pps:
        return None
    frame = Frame.of(HIST_SIZE, HIST_MARGIN)
    counts = [0] * 10
    for pp in pps:
        counts[min(int(pp * 10), 9)] += 1
    slot = frame.plot_width / 10
    total = len(pps)
    bars = [
        Bar(
            x=round(frame.left + i * slot + 1, 1),
            y=round(frame.bottom - count / total * frame.plot_height, 1),
            width=round(slot - 2, 1),
            height=round(count / total * frame.plot_height, 1),
            proportion=count / total,
            label=f"{i / 10:.1f}–{(i + 1) / 10:.1f}",
            # The palette means "warm = uncertain". Here a *low* alternative probability
            # is a *certain* site, the opposite of the ML histogram beside it, so the
            # scale is reversed to keep the colours meaning the same thing on both.
            fill=f"rgb{rgb(1 - (i + 0.5) / 10, palette)}",
        )
        for i, count in enumerate(counts)
    ]
    ticks = [(round(frame.left + k * 2 * slot, 1), f"{k * 0.2:.1f}") for k in range(6)]
    return AltAllPlot(bars, frame, ticks, total)


@dataclass(frozen=True)
class NodeChoice:
    name: str
    label: str
    clade: frozenset[str] = frozenset()  # the extant sequences below it


def node_choices(
    newick: str,
    names: Iterable[str],
    support: dict[str, str],
    called: dict[str, Any] | None = None,
) -> list[NodeChoice]:
    """Ancestral nodes in `names`, root first (preorder), each described by the extant
    sequences below it, and by its user-given name where `called` has one
    (`ancestor_names.for_selection`). Pass the outgroup-rooted tree when there is one:
    descendants are only meaningful on a rooted tree."""
    wanted = set(names)
    called = called or {}
    tree = load_tree(newick)
    choices = []
    for node in tree.preorder_internal_node_iter():
        if node.label not in wanted:
            continue
        leaves = [leaf.taxon.label for leaf in node.leaf_iter()]
        shown = f"{called[node.label].name} ({node.label})" if node.label in called else node.label
        label = f"{shown} · {len(leaves)} sequences"
        if support.get(node.label, "") not in ("", "-"):
            label += f" · support {support[node.label]}"
        choices.append(NodeChoice(node.label, label, frozenset(leaves)))
    # Nodes the tree does not label (should not happen) still get listed.
    listed = {c.name for c in choices}
    choices += [NodeChoice(name, name) for name in names if name not in listed]
    return choices


# EXPORT writes a byte-offset index beside each posteriors table, because the table is
# ~170 bytes a row and the viewers want one node (at the pilot limits hundreds of
# megabytes per combination, read on every request). A node's rows are then one ranged
# read of the table itself, with no copy of it (a per-node copy of every table was a
# quarter of a 285-combination job's bytes and 31,320 of its objects).
INDEX_FIELDS = ("node", "start", "length")


def index_name(name: str) -> str:
    """The index of a posteriors table, beside it: ancestors_posteriors.index.tsv."""
    return name.removesuffix(".tsv") + ".index.tsv"


@dataclass(frozen=True)
class ByteIndex:
    header: int  # bytes of the table's header line, which every node's rows need
    nodes: dict[str, list[tuple[int, int]]]  # node -> (start, length) of its runs of rows

    def text(self, node: str, read_range: Callable[[int, int], bytes]) -> str | None:
        """`node`'s rows with the header, as a table of its own (None if it has none)."""
        runs = self.nodes.get(node)
        if not runs:
            return None
        parts = [read_range(0, self.header)] + [read_range(start, n) for start, n in runs]
        return b"".join(parts).decode()


def index_table(path: Path) -> ByteIndex:
    """The byte offsets of each node's rows in a posteriors table (node in the first
    column). Streams; a node whose rows are not contiguous gets several runs."""
    nodes: dict[str, list[tuple[int, int]]] = {}
    with open(path, "rb") as handle:
        header = len(handle.readline())
        offset, current = header, None
        for line in handle:
            node = line.split(b"\t", 1)[0].decode()
            runs = nodes.setdefault(node, [])
            if node == current:
                start, n = runs[-1]
                runs[-1] = (start, n + len(line))
            else:
                runs.append((offset, len(line)))
                current = node
            offset += len(line)
    return ByteIndex(header, nodes)


def write_index(index: ByteIndex) -> str:
    rows = ["\t".join(INDEX_FIELDS), f"\t0\t{index.header}"]  # the header: no node
    rows += [f"{node}\t{start}\t{n}" for node, runs in index.nodes.items() for start, n in runs]
    return "\n".join(rows) + "\n"


def read_index(text: str) -> ByteIndex:
    lines = [line.split("\t") for line in text.splitlines() if line]
    if not lines or tuple(lines[0]) != INDEX_FIELDS or len(lines) < 2 or lines[1][0]:
        raise ValueError("not a posteriors index")
    nodes: dict[str, list[tuple[int, int]]] = {}
    for node, start, n in lines[2:]:
        nodes.setdefault(node, []).append((int(start), int(n)))
    return ByteIndex(int(lines[1][2]), nodes)


# An uncertainty run's own outputs (sites/, drivers/, pairs/) are a file per node with an
# index of node -> file name; file names are numbered, so no node name reaches them.
SPLIT_INDEX = "index.tsv"


def split_index(text: str) -> dict[str, str]:
    """{node: file name} from a per-node folder's index.tsv."""
    rows = (line.split("\t") for line in text.splitlines()[1:] if line)
    return {node: file for node, file in rows}


def file_range(path: Path) -> Callable[[int, int], bytes]:
    """read_range for a local table."""

    def read_range(start: int, n: int) -> bytes:
        with open(path, "rb") as handle:
            handle.seek(start)
            return handle.read(n)

    return read_range
