"""Tree viewer layout: a rectangular tree growing left to right, extant sequences on the
right, drawn by the template as SVG. Internal nodes become links to the ancestor viewer.

Three ways to place nodes horizontally (BRANCH_MODES):
- "adjusted" (default): branch lengths, with each branch clamped to
  [SHORT_FLOOR, LONG_CAP] x the tree's median branch length, so very short branches stay
  visible and a very long one does not squeeze the rest of the tree. Clamped branches are
  flagged (long ones drawn with a break mark); tooltips give every true length.
- "cladogram": no branch lengths; extant sequences aligned at the right.
- "true": branch lengths as inferred, with a scale bar.

Vertical placement is the same in every mode: extant sequences evenly spaced in
ladderized order (smaller clades first), each internal node midway between its first
and last child.

Branch display (Display): colour and thickness each by nothing, a branch support measure
(the support of the split the branch defines, i.e. of the node below it) or codeml's
free-ratio dN/dS; support labels at nodes (any measures, as numbers or stars); dN/dS
labels on branches. Every branch's tooltip lists what is known about it. Legends are
drawn inside the SVG, so downloaded figures keep them.
"""

import math
import statistics
from dataclasses import dataclass, field, replace
from typing import Any

from phylobot.ancestors import DEFAULT_PALETTE, PALETTES, rgb
from phylobot.paths import BranchRate
from phylobot.reconcile import Overlap
from phylobot.trees import Clade, canonical_split, load_tree

# Branch support measures (supports.tsv columns, then "agree") and their display names.
# "agree" is not a statistic of this tree: the share of the job's other alignment x model
# combinations (same tree source) whose trees have the same branch (views.agreement),
# labelled "k/n".
AGREE = "agree"
SUPPORT_MEASURES = (("sh_alrt", "SH-aLRT"), ("ufboot", "UFBoot"), (AGREE, "Agreement"))
# Stars for support labels: one per threshold reached (%). UFBoot: 95 is the usual
# "well supported" cut-off; SH-aLRT: 80; agreement: half, three quarters, 90%, all.
STARS = {"ufboot": (80, 90, 95, 99), "sh_alrt": (70, 80, 90, 99), AGREE: (50, 75, 90, 100)}
DNDS = "dnds"
NEUTRAL = "#1d2330"  # branches drawn without a value, and with no styling
MISSING = "#b3b8c2"  # styled by a value this branch lacks


@dataclass(frozen=True)
class Display:
    """How branches are drawn: colour and width each "" (plain), a support measure key
    or "dnds"; the support measures labelled at nodes; stars instead of numbers; dN/dS
    numbers on branches."""

    colour: str = ""
    width: str = ""
    labels: tuple[str, ...] = ("sh_alrt", "ufboot")
    stars: bool = False
    dnds_labels: bool = False
    palette: str = DEFAULT_PALETTE  # ancestors.PALETTES, for branch colours


def stars(measure: str, value: float) -> str:
    return "*" * sum(value >= t for t in STARS[measure]) or "·"


def _hex(colour: tuple[int, int, int]) -> str:
    return "#" + "".join(f"{c:02x}" for c in colour)


def support_colour(value: float, palette: str = DEFAULT_PALETTE) -> str:
    """Support (%) on the ancestor viewer's scale for posterior probabilities: warm for
    weak, cool for strong (support / 100)."""
    return _hex(rgb(value / 100, palette if palette in PALETTES else DEFAULT_PALETTE))


def dnds_position(omega: float) -> float:
    """dN/dS on a log scale, 0 at 0.1 or less, 0.5 at 1, 1 at 10 or more."""
    return min(max((math.log10(max(omega, 1e-6)) + 1) / 2, 0.0), 1.0)


def dnds_colour(omega: float, palette: str = DEFAULT_PALETTE) -> str:
    """dN/dS on the palette reversed: the cool end at 0.1 or less (purifying selection),
    the warm end at 10 or more, and 1 at the palette's 0.6, the pale middle of the
    red-blue and spectral palettes (log scale either side)."""
    palette = palette if palette in PALETTES else DEFAULT_PALETTE
    position = dnds_position(omega)
    t = 1 - 0.8 * position if position <= 0.5 else 0.6 - 1.2 * (position - 0.5)
    return _hex(rgb(t, palette))


def support_width(value: float) -> float:
    return round(1.0 + 5.0 * min(max((value - 50) / 50, 0.0), 1.0), 2)


def dnds_width(omega: float) -> float:
    return round(1.0 + 6.0 * dnds_position(omega), 2)


BRANCH_MODES = {
    "adjusted": "Branch lengths, adjusted (short ones visible, long ones cut)",
    "cladogram": "Cladogram (no branch lengths)",
    "true": "True branch lengths",
}
DEFAULT_MODE = "adjusted"
SHORT_FLOOR, LONG_CAP = 0.25, 3.0  # x the median branch length, in "adjusted" mode
TIP_SPACING = (0.3, 2.5)  # the tip-spacing slider's bounds, x the default row
# More drawn nodes than this and NodeN beside every one collides with its neighbours, so
# the viewer then labels only the ancestors the user named, unless asked for all.
CROWDED_NODES = 60

WIDTH = 1000  # SVG user units; the SVG scales to the page width
ROW = 22  # vertical space per extant sequence
MARGIN = {"left": 16, "right": 12, "top": 16, "bottom": 44}
CHAR_WIDTH = 7.4  # approximate width of a 13-unit label character
NODE_CHAR_WIDTH = 6.4  # of an 11-unit node label character
LABEL_GAP = 8
LEGEND_ROW = 26
LEGEND_PAD = 8  # above the first legend row
# A folded clade (collapse_clade) is drawn as a short triangle whatever its size, to save
# vertical space: its band is COLLAPSE_BASE rows tall, growing slowly with the clade (log),
# never past COLLAPSE_MAX_ROWS. Kept under two rows for the smallest clade, so folding any
# clade of two or more tips is always shorter than the tips it replaces.
COLLAPSE_BASE, COLLAPSE_PER, COLLAPSE_MAX_ROWS = 1.2, 0.5, 4.0


@dataclass(frozen=True)
class Branch:
    x1: float  # parent's x
    x2: float  # child's x
    y: float
    length: float | None  # as inferred (None: no length in the file)
    clamped: str = ""  # "short" or "long" in adjusted mode
    on_path: bool = False  # part of a traced mutational path
    stroke: str = NEUTRAL
    stroke_width: float = 1.6
    dashed: bool = False  # styled by dN/dS where it is undefined
    label: str = ""  # dN/dS text drawn on the branch
    tooltip: str = ""  # everything known about the branch, one item per line
    foreground: bool = False  # chosen as foreground for the branch-site test
    child: str = ""  # the node or sequence the branch leads to (path picking, tree.js)
    parent: str = ""  # the node it comes from ("" for an unlabelled root)
    parent_y: float | None = None  # the parent's y: a traced path draws the elbow too

    @property
    def middle(self) -> float:
        """Where a cut branch gets its break mark."""
        return round((self.x1 + self.x2) / 2, 1)

    @property
    def plain(self) -> bool:
        """Drawn without colour or thickness from a value."""
        return self.stroke == NEUTRAL and self.stroke_width == 1.6 and not self.dashed

    @property
    def drawn_stroke(self) -> str:
        """A plain branch on a traced path is orange; a styled one keeps its colour and
        gets an orange halo instead."""
        return "#e08214" if self.on_path and self.plain else self.stroke

    @property
    def drawn_width(self) -> float:
        return 4.0 if self.on_path and self.plain else self.stroke_width

    @property
    def halo_width(self) -> float:
        return round(self.stroke_width + 6, 2)


@dataclass(frozen=True)
class Connector:
    x: float
    y1: float
    y2: float


@dataclass(frozen=True)
class Event:
    """A reconciliation event at a node (phylobot.reconcile): "duplication" or "transfer"
    (speciations are not marked), the species-tree node where it happened, and for a
    duplication how many species its copies share (`reconcile.Overlap`): a weak one is
    drawn dashed."""

    kind: str
    species: str
    overlap: Overlap | None = None

    @property
    def weak(self) -> bool:
        return bool(self.overlap is not None and self.overlap.weak)


@dataclass(frozen=True)
class DrawnNode:
    name: str  # internal node label ("" for an unlabelled root)
    x: float
    y: float
    leaves: int
    support: str = ""
    link: bool = False  # has an ancestor to view
    href: str = ""  # its ancestor viewer URL (set by the view)
    called: str = ""  # the user's name for this ancestor (phylobot.ancestor_names)
    called_alike: bool = True  # False: same clade as where it was named, other branching
    event: Event | None = None  # a duplication or transfer, from the reconciliation


@dataclass(frozen=True)
class Leaf:
    name: str
    x: float
    y: float
    outgroup: bool = False
    note: str = ""  # drawn after the name, e.g. its species (reconciliation)


@dataclass(frozen=True)
class CollapsedClade:
    """A clade the user folded into a triangle (phylobot.clade_collapse): apex at the
    ancestral node, base at the deepest tip below it; its member tips are hidden, so the
    band is short. Clickable (data-node = name) to expand or restyle it."""

    name: str  # the folded ancestral node's label
    apex_x: float
    base_x: float
    y: float  # the band's vertical centre
    height: float
    leaves: int  # sequences folded away
    label: str = ""  # the user's text, drawn beside the triangle
    fill: str = "#cfd8e6"
    border_colour: str = "#2457a6"
    border_width: float = 1.6
    alike: bool = True  # same subtrees meet here as where it was folded (else approximate)
    called: str = ""  # the user's name for the clade's ancestor, drawn at the apex dot
    called_alike: bool = True  # that name's clade matches but the branching differs

    @property
    def top(self) -> float:
        return round(self.y - self.height / 2, 1)

    @property
    def bottom(self) -> float:
        return round(self.y + self.height / 2, 1)

    @property
    def points(self) -> str:
        """The polygon: apex on the left, vertical base on the right."""
        return f"{self.apex_x},{self.y} {self.base_x},{self.top} {self.base_x},{self.bottom}"


@dataclass
class TreeDrawing:
    width: int
    height: int
    mode: str
    branches: list[Branch] = field(default_factory=list)
    connectors: list[Connector] = field(default_factory=list)
    nodes: list[DrawnNode] = field(default_factory=list)
    leaves: list[Leaf] = field(default_factory=list)
    collapsed: list[CollapsedClade] = field(default_factory=list)
    scale_bar: tuple[float, float, float, str] | None = None  # x1, x2, y, label
    note: str = ""  # how lengths were adjusted, if they were
    legends: list[Legend] = field(default_factory=list)
    root_y: float = 0.0
    root_x: float = 0.0
    stem_from: float | None = None  # x where a stem above an added root starts
    # Legends sit above the tree: `top` is their height, and the tree (with its scale bar)
    # is drawn shifted down by it; legend y values are from the top of the figure.
    top: int = 0
    events_y: float | None = None  # the reconciliation legend row, if there is one
    doubtful_y: float | None = None  # its second row, for doubtful duplications

    @property
    def weak_events(self) -> bool:
        """Any duplication drawn as doubtful (its legend row is shown only then)."""
        return any(n.event is not None and n.event.weak for n in self.nodes)


@dataclass(frozen=True)
class LegendItem:
    label: str
    stroke: str
    width: float
    dashed: bool = False
    x: float = 0.0  # set by layout


@dataclass(frozen=True)
class Legend:
    title: str
    items: list[LegendItem]
    y: float = 0.0  # set by layout


def legend(
    kind: str, styled_by: str, names: dict[str, str], palette: str = DEFAULT_PALETTE
) -> Legend:
    """The legend for branch colour or width ("colour"/"width") styled by a support
    measure or dN/dS."""
    colour = kind == "colour"
    name = "colour" if colour else "thickness"
    if styled_by == DNDS:
        ratios = (0.1, 0.3, 1.0, 3.0, 10.0)
        items = [
            LegendItem(f"{v:g}", dnds_colour(v, palette) if colour else NEUTRAL,
                       2.6 if colour else dnds_width(v))
            for v in ratios
        ]  # fmt: skip
        items.append(LegendItem("undefined", MISSING, 1.6, dashed=True))
        title = f"Branch {name}: dN/dS (codeml free ratios)"
    else:
        values = (
            (20.0, 50.0, 70.0, 80.0, 90.0, 95.0, 100.0)
            if colour
            else (50.0, 70.0, 80.0, 90.0, 95.0, 100.0)
        )
        items = [
            LegendItem(f"{v:g}{'' if v > 50 or colour else ' or less'}",
                       support_colour(v, palette) if colour else NEUTRAL,
                       2.6 if colour else support_width(v))
            for v in values
        ]  # fmt: skip
        items.append(LegendItem("none", MISSING, 1.6))
        title = f"Branch {name}: {names.get(styled_by, styled_by)} support (%)"
        if styled_by == AGREE:
            title = f"Branch {name}: share of the other combinations with this branch (%)"
    return Legend(title, items)


def nice_length(target: float) -> float:
    """A 1/2/5 x 10^k length no larger than `target` (for the scale bar)."""
    magnitude = 10.0 ** int(f"{target:e}".split("e")[1])
    for factor in (5, 2, 1):
        if factor * magnitude <= target:
            return factor * magnitude
    return magnitude / 2


def layout(
    newick: str,
    mode: str = DEFAULT_MODE,
    ancestors: set[str] | None = None,
    support: dict[str, str] | None = None,
    outgroup: set[str] | None = None,
    path: list[str] | None = None,
    supports: dict[str, dict[str, float]] | None = None,
    rates: dict[frozenset[str], BranchRate] | None = None,
    display: Display | None = None,
    foreground: set[str] | None = None,
    events: dict[str, Event] | None = None,
    losses: dict[Clade, list[str]] | None = None,
    agreement: dict[str, tuple[int, int]] | None = None,
    tip_notes: dict[str, str] | None = None,
    collapsed: dict[str, Any] | None = None,
) -> TreeDrawing:
    """Coordinates for drawing `newick`. `ancestors` are the node names with a
    reconstruction (they become links); `support` maps node names to support text
    (overrides `supports` for labels when given); `supports` maps node names to
    {measure: %} (supports.tsv); `rates` codeml's free-ratio estimates by branch ends
    (paths.branch_rates); `display` how to style branches and labels; `path` (names,
    oldest first) marks the branches of a traced path; `foreground` names the nodes
    whose branches are the branch-site test's foreground (phylobot.branch_site_views);
    `events` the reconciliation's duplications and transfers by node name, and `losses`
    the species that lost a copy on each branch, by its unrooted split
    (`trees.canonical_split`), so they land on the right branch however it is drawn.
    `collapsed` maps node names the user folded to a style object (label, fill,
    border_colour, border_width): each becomes a triangle and its tips are hidden."""
    ancestors = ancestors or set()
    events = events or {}
    losses = losses or {}
    outgroup = outgroup or set()
    supports = supports or {}
    rates = rates or {}
    display = display or Display()
    names = dict(SUPPORT_MEASURES)
    measures = [m for m, _ in SUPPORT_MEASURES if any(m in v for v in supports.values())]

    agreement = agreement or {}

    def label_text(node_name: str) -> str:
        values = supports.get(node_name, {})
        shown = [m for m in measures if m in display.labels and m != AGREE]
        parts = []
        if shown and any(m in values for m in shown):
            if display.stars:
                parts.append("/".join(stars(m, values[m]) if m in values else "-" for m in shown))
            else:
                parts.append("/".join(f"{values[m]:g}" if m in values else "-" for m in shown))
        if AGREE in display.labels and node_name in agreement:
            k, n = agreement[node_name]
            parts.append(stars(AGREE, values[AGREE]) if display.stars else f"{k}/{n}")
        return " · ".join(parts)

    if support is None:
        support = {n: text for n in supports if (text := label_text(n))}
    if mode not in BRANCH_MODES:
        mode = DEFAULT_MODE
    tree = load_tree(newick)
    tree.ladderize(ascending=True)
    root = tree.seed_node
    below = [n for n in tree.preorder_node_iter() if n is not root]
    lengths = [n.edge.length for n in below if n.edge.length is not None]
    positive = [x for x in lengths if x > 0]
    median = statistics.median(positive) if positive else 1.0
    floor, cap = SHORT_FLOOR * median, LONG_CAP * median

    # Horizontal position: distance from the root under the chosen mode.
    drawn_length: dict[int, float] = {}
    clamped: dict[int, str] = {}
    for node in below:
        length = node.edge.length
        if mode == "true":
            drawn_length[id(node)] = length or 0.0
        elif mode == "adjusted":
            value = length or 0.0
            if value < floor:
                drawn_length[id(node)], clamped[id(node)] = floor, "short"
            elif value > cap:
                drawn_length[id(node)], clamped[id(node)] = cap, "long"
            else:
                drawn_length[id(node)] = value
    depth: dict[int, float] = {}
    if mode == "cladogram":
        height: dict[int, int] = {}
        for node in tree.postorder_node_iter():
            kids = node.child_nodes()
            height[id(node)] = 0 if not kids else 1 + max(height[id(k)] for k in kids)
        top = height[id(root)]
        for node in tree.preorder_node_iter():
            depth[id(node)] = float(top - height[id(node)])
    else:
        for node in tree.preorder_node_iter():
            parent = node.parent_node
            depth[id(node)] = 0.0 if parent is None else depth[id(parent)] + drawn_length[id(node)]

    leaves = list(tree.leaf_node_iter())
    everyone = frozenset(leaf.taxon.label for leaf in leaves)
    tip_notes = tip_notes or {}

    # Clades the user folded into triangles (phylobot.clade_collapse): map the requested
    # node labels to tree nodes (internal, labelled, never the root), then find for every
    # node the outermost fold it sits inside. A node strictly inside a fold is hidden; the
    # fold itself is drawn as a triangle only where no other fold encloses it.
    collapsed = collapsed or {}
    by_label = {n.label: n for n in tree.preorder_internal_node_iter() if n.label and n is not root}
    folded = {by_label[name]: style for name, style in collapsed.items() if name in by_label}
    folded_set = set(folded)
    top_fold: dict[int, Any] = {}
    for node in tree.preorder_node_iter():
        parent = node.parent_node
        inherited = top_fold.get(id(parent)) if parent is not None else None
        top_fold[id(node)] = inherited or (node if node in folded_set else None)

    def hidden(node: Any) -> bool:
        fold = top_fold.get(id(node))
        return fold is not None and fold is not node

    drawn_folds = [n for n in folded_set if not hidden(n)]  # the outermost folds (drawn)
    visible_leaves = [leaf for leaf in leaves if top_fold.get(id(leaf)) is None]

    longest = max(
        [len(leaf.taxon.label) + (len(tip_notes[leaf.taxon.label]) + 1
                                  if leaf.taxon.label in tip_notes else 0)
         for leaf in visible_leaves]
        + [len(folded[n].label) for n in drawn_folds if folded[n].label],
        default=4,
    )  # fmt: skip
    label_space = min(360.0, longest * CHAR_WIDTH + LABEL_GAP + 4)
    # A labelled root sits at the left edge with its label to its left: make room.
    root_text = f"{root.label} {support.get(root.label, '')}".strip() if root.label else ""
    left = MARGIN["left"] + (len(root_text) * NODE_CHAR_WIDTH + 6 if root_text else 0)
    plot_width = WIDTH - left - MARGIN["right"] - label_space
    deepest = max(depth.values()) or 1.0
    stem = 0.0 if root.label else 12.0  # a short stem above an added (unlabelled) root

    def x(node: object) -> float:
        return left + stem + depth[id(node)] / deepest * (plot_width - stem)

    def clade_height(count: int) -> float:
        return ROW * min(COLLAPSE_MAX_ROWS, COLLAPSE_BASE + COLLAPSE_PER * math.log2(max(count, 2)))

    # Vertical order of the drawn terminals: each visible leaf (one ROW), and each folded
    # clade once in place of its hidden leaves (a short band), so the figure stands shorter
    # than the sequences it represents.
    terminals: list[Any] = []
    seen: set[int] = set()
    for leaf in leaves:
        fold = top_fold.get(id(leaf))
        if fold is None:
            terminals.append(leaf)
        elif id(fold) not in seen:
            seen.add(id(fold))
            terminals.append(fold)

    y: dict[int, float] = {}
    cursor = float(MARGIN["top"])
    for terminal in terminals:
        span = ROW if terminal.is_leaf() else clade_height(len(terminal.leaf_nodes()))
        y[id(terminal)] = cursor + span / 2
        cursor += span
    for node in tree.postorder_internal_node_iter():
        if node in folded_set or hidden(node):
            continue  # a folded node already has its band y; a hidden node is not drawn
        kids = node.child_nodes()
        y[id(node)] = (y[id(kids[0])] + y[id(kids[-1])]) / 2

    drawing = TreeDrawing(
        width=WIDTH,
        height=int(cursor + MARGIN["bottom"]),
        mode=mode,
        root_x=x(root),
        root_y=y[id(root)],
        stem_from=left if stem else None,
    )
    on_path = set((path or [])[1:])  # branches leading to the path's later members
    fg = foreground or set()

    def colour_of(by: str, values: dict[str, float], rate: BranchRate | None) -> str:
        if by == DNDS:
            if rate is not None and rate.defined:
                return dnds_colour(rate.omega, display.palette)
            return MISSING
        return support_colour(values[by], display.palette) if by in values else MISSING

    def width_of(by: str, values: dict[str, float], rate: BranchRate | None) -> float:
        if by == DNDS:
            return dnds_width(rate.omega) if rate is not None and rate.defined else 1.6
        return support_width(values[by]) if by in values else 1.6

    for node in below:
        if hidden(node):
            continue  # a branch inside a folded clade is not drawn
        length = node.edge.length
        name = node.taxon.label if node.is_leaf() else (node.label or "")
        parent = node.parent_node.label or ""
        values = supports.get(name, {}) if not node.is_leaf() else {}
        rate = rates.get(frozenset((parent, name))) if parent else None
        stroke = colour_of(display.colour, values, rate) if display.colour else NEUTRAL
        width = width_of(display.width, values, rate) if display.width else 1.6
        # dN/dS styling where it is undefined: dashed
        undefined = rate is None or not rate.defined
        dashed = DNDS in (display.colour, display.width) and undefined
        tip = [f"{parent or 'root'} → {name}" if name else "branch"]
        tip.append(f"length {length:.4g}" if length is not None else "length not given")
        tip += [f"{names[m]} {values[m]:g}%" for m in measures if m in values]
        if rate is not None:
            omega = f"{rate.omega:.3g}" if rate.defined else "undefined (< 1 synonymous change)"
            tip.append(f"dN/dS {omega}; dN {rate.dn:.4g}, dS {rate.ds:.4g}; expected changes "
                       f"≈ {rate.n_dn:.1f} non-synonymous, {rate.s_ds:.1f} synonymous")  # fmt: skip
        split = canonical_split(frozenset(n.taxon.label for n in node.leaf_iter()), everyone)
        if lost := losses.get(split):
            copies = "copy" if len(lost) == 1 else "copies"
            tip.append(f"{len(lost)} {copies} lost here (in {'; '.join(lost)})")
        label = ""
        if display.dnds_labels and rate is not None:
            label = f"{rate.omega:.2f}" if rate.defined else "–"
        drawing.branches.append(Branch(
            round(x(node.parent_node), 1), round(x(node), 1), round(y[id(node)], 1),
            length, clamped.get(id(node), ""), name in on_path, stroke, width, dashed, label,
            "\n".join(tip + (["foreground (branch-site test)"] if name in fg else [])),
            name in fg, child=name, parent=parent, parent_y=round(y[id(node.parent_node)], 1),
        ))  # fmt: skip
    for node in tree.preorder_internal_node_iter():
        if node in folded_set or hidden(node):
            continue  # a folded node is drawn as a triangle below, not as a dot
        kids = node.child_nodes()
        drawing.connectors.append(
            Connector(round(x(node), 1), round(y[id(kids[0])], 1), round(y[id(kids[-1])], 1))
        )
        name = node.label or ""
        drawing.nodes.append(DrawnNode(
            name, round(x(node), 1), round(y[id(node)], 1), len(node.leaf_nodes()),
            support.get(name, "") if support.get(name, "") not in ("", "-") else "",
            name in ancestors, event=events.get(name),
        ))  # fmt: skip
    for leaf in visible_leaves:
        name = leaf.taxon.label
        position = round(x(leaf), 1), round(y[id(leaf)], 1)
        drawing.leaves.append(
            Leaf(name, *position, outgroup=name in outgroup, note=tip_notes.get(name, ""))
        )
    for node in drawn_folds:
        style = folded[node]
        apex = x(node)
        base = max((x(leaf) for leaf in node.leaf_iter()), default=apex)
        count = len(node.leaf_nodes())
        drawing.collapsed.append(CollapsedClade(
            name=node.label or "",
            apex_x=round(apex, 1),
            base_x=round(max(base, apex + 6), 1),
            y=round(y[id(node)], 1),
            height=round(clade_height(count), 1),
            leaves=count,
            label=style.label,
            fill=style.fill,
            border_colour=style.border_colour,
            border_width=style.border_width,
            alike=getattr(style, "alike", True),
        ))  # fmt: skip

    # Legends above the tree, each one row: title, then samples.
    for kind, by in (("colour", display.colour), ("width", display.width)):
        if by:
            drawing.legends.append(legend(kind, by, names, display.palette))
    legend_top = drawing.height  # the tree's own height: its scale bar sits below it
    placed = []
    for i, entry in enumerate(drawing.legends):
        row_y = LEGEND_PAD + LEGEND_ROW * i + LEGEND_ROW / 2
        x_at = MARGIN["left"] + len(entry.title) * NODE_CHAR_WIDTH + 14
        items = []
        for item in entry.items:
            items.append(
                LegendItem(item.label, item.stroke, item.width, item.dashed, round(x_at, 1))
            )
            x_at += 30 + len(item.label) * NODE_CHAR_WIDTH + 14
        placed.append(Legend(entry.title, items, round(row_y, 1)))
    drawing.legends = placed
    rows = len(placed)
    if any(n.event for n in drawing.nodes) or losses:
        drawing.events_y = round(LEGEND_PAD + LEGEND_ROW * rows + LEGEND_ROW / 2, 1)
        rows += 1
        if drawing.weak_events:
            drawing.doubtful_y = round(LEGEND_PAD + LEGEND_ROW * rows + LEGEND_ROW / 2, 1)
            rows += 1
    top = LEGEND_PAD + LEGEND_ROW * rows if rows else 0
    drawing.top = top
    drawing.height += top

    bar_y = legend_top - MARGIN["bottom"] / 2
    if mode == "true" and deepest > 0:
        length = nice_length(deepest / 5)
        units = length / deepest * (plot_width - stem)
        drawing.scale_bar = (
            left, round(left + units, 1), bar_y,
            f"{length:g} substitutions per site",
        )  # fmt: skip
    elif mode == "adjusted":
        short = sum(1 for c in clamped.values() if c == "short")
        long = sum(1 for c in clamped.values() if c == "long")
        parts = []
        if short:
            parts.append(f"{short} short branch{'es' if short != 1 else ''} lengthened")
        if long:
            parts.append(f"{long} long branch{'es' if long != 1 else ''} cut (marked //)")
        drawing.note = (
            f"Branch lengths adjusted to {SHORT_FLOOR:g}–{LONG_CAP:g} × the median "
            f"({median:.3g}): "
            + ("; ".join(parts) if parts else "none needed")
            + ". Hover a branch for its true length."
        )
    return drawing


def slider(value: Any, bounds: tuple[float, float], default: float = 1.0) -> float:
    """A slider's value from the query, clamped to its bounds; the default if unreadable."""
    try:
        number = float(value)
    except TypeError, ValueError:
        return default
    if not math.isfinite(number):
        return default
    return min(max(number, bounds[0]), bounds[1])


def stretch(drawing: TreeDrawing, factor: float, top: float | None = None,
            bottom: float | None = None) -> TreeDrawing:  # fmt: skip
    """The drawing with its rows `factor` times as far apart (the tip-spacing slider): y
    measured from `top` (the top margin); `bottom` the margin kept under the last row."""
    if factor == 1.0:
        return drawing
    top = MARGIN["top"] if top is None else top
    bottom = MARGIN["bottom"] if bottom is None else bottom

    def y(value: float) -> float:
        return round(top + (value - top) * factor, 2)

    drawing.branches = [
        replace(b, y=y(b.y), parent_y=None if b.parent_y is None else y(b.parent_y))
        for b in drawing.branches
    ]
    drawing.connectors = [replace(c, y1=y(c.y1), y2=y(c.y2)) for c in drawing.connectors]
    drawing.nodes = [replace(n, y=y(n.y)) for n in drawing.nodes]
    drawing.leaves = [replace(lf, y=y(lf.y)) for lf in drawing.leaves]
    drawing.collapsed = [
        replace(c, y=y(c.y), height=round(c.height * factor, 2)) for c in drawing.collapsed
    ]
    drawing.root_y = y(drawing.root_y)
    lowest = max([lf.y for lf in drawing.leaves] + [c.bottom for c in drawing.collapsed],
                 default=top)  # fmt: skip
    if drawing.scale_bar:
        x1, x2, _, label = drawing.scale_bar
        drawing.scale_bar = (x1, x2, lowest + 24, label)
    drawing.height = int(drawing.top + lowest + bottom)
    return drawing
