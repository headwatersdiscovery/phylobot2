"""The Uncertainty page's overview figures (uncertainty_views lays out the rest): the
explainer at the top, a map of the job's reconstructions, and a pie per ancestor on the
reference tree. Pure: they read a run's summary and per-ancestor parts, never storage.

What "two reconstructions chosen at random" means, which the explainer draws: the job
reconstructed each ancestor once per cell (alignment x tree x model x gap version). Pick
two cells at random by their weights (possibly the same cell twice), sample a sequence
from each cell's posterior, and count the positions where the two differ. The expected
count splits into the posterior's part (what two samples of one cell already disagree
on) and each choice's part (Shapley values).

The map: `summary["pairs"][mode][i][j]` is that expected count for cells i and j, averaged
over the ancestors. With D_ii a cell's own disagreement, D_ij - (D_ii + D_jj) / 2 is the
squared distance between the two cells' posterior means (sum over columns of
||p_i - p_j||^2), so classical MDS places the cells exactly as far as the 2-D projection
allows. A halo of radius sqrt(D_ii / 2), the root-mean-square distance of one sample from
its cell's mean, shows the posterior's own spread on the same scale.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

# The --viz slots (uncertainty_views.FILLS), repeated to keep this module free of Django.
FILLS = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a", 4: "#eda100", 5: "#e87ba4", 6: "#7d3c98"}
INK, MUTED, GRID = "#1d2330", "#5d6675", "#d8dde5"


def fmt(value: float) -> str:
    return f"{value:.1f}" if abs(value) < 100 else f"{value:.0f}"


# --- The explainer -------------------------------------------------------------------


@dataclass(frozen=True)
class Explainer:
    cells: int
    design: str  # "5 alignments × 2 trees × 2 models × 2 gap versions"
    grid: list[tuple[float, float, bool]]  # x, y, picked: the cells as small squares
    total: str
    bar: list[tuple[float, float, str, str]]  # x, width, fill, label (the split, 0..W)
    posterior_end: float  # where the posterior's part ends on the bar
    choices_label: str  # "alignment 26.1 · model 0.4 · …", the larger choices first


GRID_COLUMNS = 10
BAR_WIDTH = 250.0

PLURALS = {"Alignment": "alignments", "Tree": "trees", "Model": "models",
           "Gap placement": "gap versions"}  # fmt: skip


def explainer(
    summary: dict[str, Any], parts: dict[str, dict[str, float]], components: tuple[Any, ...]
) -> Explainer:
    """The three-step picture: the job's cells, two picked and sampled, the split bar."""
    cells = len(summary.get("cells", []))
    factors = [f for f in summary.get("factors", []) if f.get("levels")]
    sizes = [
        f"{len(f['levels'])} {PLURALS.get(f['label'], f['label'].lower() + 's')}"
        if len(f["levels"]) != 1
        else f"1 {PLURALS.get(f['label'], f['label'].lower() + 's').rstrip('s')}"
        for f in factors
    ]
    shown = min(cells, 60)
    picked = {0, min(shown - 1, GRID_COLUMNS + 3)} if shown > 1 else {0}
    grid = [
        (float((i % GRID_COLUMNS) * 13), float((i // GRID_COLUMNS) * 13), i in picked)
        for i in range(shown)
    ]
    whole = parts["all"]
    total = whole["total"] or 1.0
    bar, x = [], 0.0
    for key, label, colour in components:
        value = whole.get(key, 0.0)
        if value <= 0:
            continue
        width = BAR_WIDTH * value / total
        bar.append((x, width, FILLS[colour], f"{label}: {fmt(value)} positions"))
        x += width
    posterior_end = BAR_WIDTH * whole.get("within", 0.0) / total
    choices = sorted(
        ((label, whole.get(key, 0.0)) for key, label, _ in components if key != "within"),
        key=lambda kv: -kv[1],
    )
    choices_label = " · ".join(f"{label.lower()} {fmt(v)}" for label, v in choices if v >= 0.05)
    return Explainer(
        cells=cells,
        design=" × ".join(sizes),
        grid=grid,
        total=fmt(whole["total"]),
        bar=bar,
        posterior_end=posterior_end,
        choices_label=choices_label or "none",
    )


# --- The map of reconstructions --------------------------------------------------------


@dataclass(frozen=True)
class MapPoint:
    x: float
    y: float
    fill: str  # its level's colour (of the factor the dots are coloured by)
    level: str  # that level's label
    tip: str


@dataclass(frozen=True)
class MapLabel:
    x: float
    y: float
    text: str


@dataclass(frozen=True)
class MapChoice:
    """One factor's checkboxes: (value "factor:level", label, ticked)."""

    label: str
    levels: list[tuple[str, str, bool]]


@dataclass(frozen=True)
class ReconstructionMap:
    points: list[MapPoint]
    halos: list[tuple[float, float, float]]  # x, y, r: one per spatial cluster
    labels: list[MapLabel]  # a level's name where its dots sit together
    legend: list[tuple[str, str, int]]  # fill, level, dots
    by: str  # the factor the dots are coloured by ("alignment")
    colour_by: str  # its key
    colour_choices: list[tuple[str, str]]  # varied factors: key, label
    choices: list[MapChoice]  # every varied factor's levels, to show or hide
    hidden: int  # reconstructions left out
    size: int
    shown: int  # percent of the spread the two axes show
    largest: str  # the farthest pair shown, in words
    noise: str  # the mean posterior spread, in positions
    note: str = ""  # why nothing is drawn (too few left)


MAP_SIZE = 420
MAP_PAD = 46
MERGE = 34.0  # px: labels closer than this join ("MAFFT, MUSCLE")
CLUSTER = 24.0  # px: dots closer than this share one halo
COMPACT = 40.0  # px: a level is labelled on the map only if its dots sit within this
# Level colours: the --viz slots (remapped in dark.css), then ink.
LEVEL_FILLS = [FILLS[2], FILLS[1], FILLS[3], FILLS[4], FILLS[5], FILLS[6], INK]


def map_shown(params: Any) -> set[str] | None:
    """The "factor:level" values ticked in the map's form, or None (all) before it is used.
    An unticked box sends nothing, so the form carries map_set=1 beside the boxes."""
    if not params.get("map_set"):
        return None
    return set(params.getlist("show"))


def reconstruction_map(
    summary: dict[str, Any],
    mode: str,
    parts: dict[str, dict[str, float]],
    components: tuple[Any, ...],
    colour_by: str = "",
    shown: set[str] | None = None,
) -> ReconstructionMap | None:
    """The cells shown as dots, placed by classical MDS on the distances between their
    means; laid out afresh on whichever cells are shown."""
    cells = summary.get("cells", [])
    matrix = (summary.get("pairs") or {}).get(mode)
    n = len(cells)
    if n < 3 or not matrix or len(matrix) != n:
        return None
    factors = summary.get("factors", [])
    varied = [f for f in factors if f.get("varied") and len(f.get("levels", [])) > 1]
    if not varied:
        return None
    # Colour by the asked factor, else the varied choice that adds the most.
    whole = parts["all"]
    keys = [f["key"] for f in varied]
    if colour_by not in keys:
        colour_by = max(keys, key=lambda k: whole.get(k, 0.0))
    position = {f["key"]: i for i, f in enumerate(factors)}
    by = next(f for f in varied if f["key"] == colour_by)
    names = {lv["key"]: lv["label"] for lv in by["levels"]}
    fills = {lv["key"]: LEVEL_FILLS[i % len(LEVEL_FILLS)] for i, lv in enumerate(by["levels"])}

    def level_of(cell: dict[str, Any], key: str) -> str:
        levels, i = cell.get("levels", []), position.get(key, -1)
        return levels[i] if 0 <= i < len(levels) else ""

    every = {f"{f['key']}:{lv['key']}" for f in varied for lv in f["levels"]}
    ticked = every if shown is None else shown & every
    keep = [i for i, c in enumerate(cells)
            if all(f"{f['key']}:{level_of(c, f['key'])}" in ticked for f in varied)]  # fmt: skip
    choices = [
        MapChoice(f["label"], [(f"{f['key']}:{lv['key']}", lv["label"],
                                f"{f['key']}:{lv['key']}" in ticked) for lv in f["levels"]])
        for f in varied
    ]  # fmt: skip
    base = dict(by=by["label"].lower(), colour_by=colour_by,
                colour_choices=[(f["key"], f["label"]) for f in varied], choices=choices,
                hidden=n - len(keep), size=MAP_SIZE)  # fmt: skip
    if len(keep) < 3:
        return ReconstructionMap(points=[], halos=[], labels=[], legend=[], shown=0, largest="",
                                 noise="", note="Fewer than three reconstructions are ticked: "
                                 "tick more to draw the map.", **base)  # fmt: skip
    d = np.asarray(matrix, dtype=float)[np.ix_(keep, keep)]
    m = len(keep)
    own = np.diag(d).copy()
    squared = np.clip(d - (own[:, None] + own[None, :]) / 2, 0.0, None)
    np.fill_diagonal(squared, 0.0)
    centre = np.eye(m) - 1.0 / m
    values, vectors = np.linalg.eigh(-0.5 * centre @ squared @ centre)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    positive = values.clip(min=0)
    if positive.sum() <= 1e-9:
        same = "The reconstructions ticked are all the same: there is nothing to place."
        return ReconstructionMap(points=[], halos=[], labels=[], legend=[], shown=0, largest="",
                                 noise=fmt(float(own.mean())), note=same, **base)  # fmt: skip
    coords = vectors[:, :2] * np.sqrt(positive[:2])
    explained = round(100 * positive[:2].sum() / positive.sum())
    spread = np.sqrt(own / 2) if mode == "posterior" else np.zeros(m)
    # One scale for both axes, fitting the dots and their posterior spread.
    span = float(max(
        (coords[:, 0] + spread).max() - (coords[:, 0] - spread).min(),
        (coords[:, 1] + spread).max() - (coords[:, 1] - spread).min(),
        1e-9,
    ))  # fmt: skip
    k = (MAP_SIZE - 2 * MAP_PAD) / span
    mid = (coords.max(axis=0) + coords.min(axis=0)) / 2
    xy = (coords - mid) * k + MAP_SIZE / 2
    points = []
    for row, i in enumerate(keep):
        cell = cells[i]
        level = level_of(cell, colour_by)
        version = cell.get("version", "")
        label = cell.get("label", cell.get("key", ""))
        tip = (f"{label}{' · ' + version if version else ''}: "
               f"own posterior {fmt(own[row])} positions")  # fmt: skip
        points.append(MapPoint(float(xy[row, 0]), float(xy[row, 1]), fills.get(level, INK),
                               names.get(level, level), tip))  # fmt: skip
    halos = [
        (x, y, r * k)
        for x, y, r in clusters_of([(p.x, p.y) for p in points], spread.tolist(), CLUSTER)
        if mode == "posterior"
    ]
    points = fanned(points)
    counts: dict[str, int] = {}
    for p in points:
        counts[p.level] = counts.get(p.level, 0) + 1
    legend = [(fills[lv["key"]], lv["label"], counts[lv["label"]])
              for lv in by["levels"] if lv["label"] in counts]  # fmt: skip
    far = np.unravel_index(int(np.argmax(squared)), squared.shape)
    i, j = int(far[0]), int(far[1])
    largest = (
        f"{cells[keep[i]].get('label', '')} and {cells[keep[j]].get('label', '')}: "
        f"{fmt(float(d[i, j]))} positions apart"
        if squared.max() > 0
        else ""
    )
    return ReconstructionMap(
        points=points,
        halos=halos,
        labels=map_labels(points, [h for h in halos]),
        legend=legend,
        shown=int(explained),
        largest=largest,
        noise=fmt(float(own.mean())),
        **base,
    )


FAN = 4.0  # px: dots closer than this are fanned out on a small ring, so none hides
FAN_RING = 6.0


def fanned(points: list[MapPoint]) -> list[MapPoint]:
    """Dots that would sit on top of each other moved apart on a ring around their centre
    (a few pixels: they barely differ, which the hint says), so every colour shows."""
    out = list(points)
    for x, y, _ in clusters_of([(p.x, p.y) for p in points], [0.0] * len(points), FAN):
        members = [i for i, p in enumerate(points) if math.dist((p.x, p.y), (x, y)) < FAN]
        if len(members) < 2:
            continue
        ring = FAN_RING * max(1.0, len(members) / 6)
        for n, i in enumerate(members):
            angle = 2 * math.pi * n / len(members) - math.pi / 2
            p = points[i]
            out[i] = MapPoint(x + ring * math.cos(angle), y + ring * math.sin(angle),
                              p.fill, p.level, p.tip)  # fmt: skip
    return out


def clusters_of(
    xy: list[tuple[float, float]], radius: list[float], gap: float
) -> list[tuple[float, float, float]]:
    """Single-linkage clusters of points closer than `gap`: centroid and mean radius."""
    parent = list(range(len(xy)))

    def root(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a in range(len(xy)):
        for b in range(a + 1, len(xy)):
            if math.dist(xy[a], xy[b]) < gap:
                parent[root(a)] = root(b)
    groups: dict[int, list[int]] = {}
    for a in range(len(xy)):
        groups.setdefault(root(a), []).append(a)
    return [
        (sum(xy[a][0] for a in g) / len(g), sum(xy[a][1] for a in g) / len(g),
         sum(radius[a] for a in g) / len(g))
        for g in groups.values()
    ]  # fmt: skip


def map_labels(points: list[MapPoint], halos: list[tuple[float, float, float]]) -> list[MapLabel]:
    """A level's name above its dots where they sit together (within COMPACT); names whose
    places fall closer than MERGE join into one."""
    placed: list[tuple[float, float, list[str], float]] = []  # x, y, names, top
    for level in dict.fromkeys(p.level for p in points):
        members = [p for p in points if p.level == level]
        x = sum(p.x for p in members) / len(members)
        y = sum(p.y for p in members) / len(members)
        if max(math.dist((p.x, p.y), (x, y)) for p in members) > COMPACT:
            continue  # spread over the map: the legend names it
        ring = max((r for hx, hy, r in halos if math.dist((hx, hy), (x, y)) < CLUSTER), default=0.0)
        top = min(p.y for p in members) - max(ring - (y - min(p.y for p in members)), 0) - 12
        for n, (px, py, names, ptop) in enumerate(placed):
            if math.dist((px, py), (x, y)) < MERGE:
                placed[n] = (px, py, [*names, level], min(ptop, top))
                break
        else:
            placed.append((x, y, [level], top))
    return [MapLabel(x, max(top, 14.0), ", ".join(names)) for x, _, names, top in placed]


# --- Pies on the tree ------------------------------------------------------------------


@dataclass(frozen=True)
class Pie:
    node: str
    x: float
    y: float
    r: float
    slices: list[tuple[str, str, str]]  # path d, fill (colour or hatch url), tip
    href: str
    tip: str


PIE_MAX = 10.0  # px: the radius of the most uncertain ancestor's pie


def pies(
    places: dict[str, tuple[float, float]],
    node_parts: dict[str, dict[str, dict[str, float]]],
    hrefs: dict[str, str],
    components: tuple[Any, ...],
    prefix: str,
    scale: float = 1.0,
) -> list[Pie]:
    """A pie at each ancestor drawn on the tree: its area in proportion to the ancestor's
    total (its bar's length), its slices the bar's segments in the bar's order, solid for
    which residue and hatched for residue or gap."""
    top = max((node_parts[n]["all"]["total"] for n in places if n in node_parts), default=0.0)
    out = []
    for node, (x, y) in places.items():
        if node not in node_parts or top <= 0:
            continue
        p = node_parts[node]
        total = p["all"]["total"]
        if total <= 0:
            continue
        r = max(PIE_MAX * scale * (total / top) ** 0.5, 2.0)
        pieces = []
        for key, label, colour in components:
            for part, hatched, what in (("residue", False, "which residue"),
                                        ("presence", True, "residue or gap")):  # fmt: skip
                v = p.get(part, {}).get(key, 0.0)
                if v > 0:
                    fill = f"url(#{prefix}-hatch-{colour})" if hatched else FILLS[colour]
                    pieces.append((v, fill, f"{node} · {label}, {what}: {fmt(v)} positions"))
        whole = sum(v for v, _, _ in pieces) or 1.0
        slices, start = [], 0.0
        for v, fill, tip in pieces:
            end = start + v / whole
            slices.append((arc(x, y, r, start, end), fill, tip))
            start = end
        out.append(Pie(node, x, y, r, slices, hrefs.get(node, ""),
                       f"{node}: {fmt(total)} positions"))  # fmt: skip
    return out


def arc(x: float, y: float, r: float, start: float, end: float) -> str:
    """The SVG path of a pie slice from `start` to `end` (fractions of a turn, clockwise
    from 12 o'clock); a whole turn is a circle."""
    if end - start >= 0.9999:
        return (f"M {x - r:.2f} {y:.2f} a {r:.2f} {r:.2f} 0 1 0 {2 * r:.2f} 0 "
                f"a {r:.2f} {r:.2f} 0 1 0 {-2 * r:.2f} 0 Z")  # fmt: skip

    def point(t: float) -> tuple[float, float]:
        angle = 2 * math.pi * t - math.pi / 2
        return x + r * math.cos(angle), y + r * math.sin(angle)

    (x0, y0), (x1, y1) = point(start), point(end)
    large = 1 if end - start > 0.5 else 0
    return (f"M {x:.2f} {y:.2f} L {x0:.2f} {y0:.2f} "
            f"A {r:.2f} {r:.2f} 0 {large} 1 {x1:.2f} {y1:.2f} Z")  # fmt: skip


# --- Spacing the tree ------------------------------------------------------------------

PIE_SIZES = (0.5, 3.0)  # the pie-size slider's range (x PIE_MAX)
TREE_SPACING = (0.3, 2.5)  # the tip-spacing slider's range (treeview.TIP_SPACING)


def slider(value: Any, bounds: tuple[float, float], default: float = 1.0) -> float:
    """A slider's value from the query (treeview.slider)."""
    from phylobot import treeview  # dendropy: only when a tree is drawn

    return treeview.slider(value, bounds, default)


def stretch(drawing: Any, factor: float, top: float, bottom: float) -> Any:
    """The tree drawing with its rows `factor` times as far apart (treeview.stretch)."""
    from phylobot import treeview

    return treeview.stretch(drawing, factor, top, bottom)


# --- What it means ----------------------------------------------------------------------


def percent(part: float, whole: float) -> str:
    return f"{round(100 * part / whole)}%" if whole > 0 else "0%"


def takeaways(
    summary: dict[str, Any],
    parts: dict[str, dict[str, float]],
    totals: dict[str, float],
    components: tuple[Any, ...],
    mode: str,
) -> list[str]:
    """The page's "so what", in sentences: what the number is, over what, how much is the
    posterior's and how much the choices', which choice matters most, and the spread over
    the ancestors (`totals`: each ancestor's total)."""
    whole = parts["all"]
    total, within = whole["total"], whole.get("within", 0.0)
    choices = total - within
    nodes, columns = summary.get("nodes", len(totals)), summary.get("columns", 0)
    out = [
        f"This job has {nodes} ancestors (every internal node of the reference tree, each "
        f"matched to the same clade in the other combinations). Take any two of an "
        f"ancestor's reconstructions and a sequence from each: on average they differ at "
        f"{fmt(total)} of its {columns} alignment columns. That is not one sampled pair but "
        f"the exact average over every ancestor, every pair of reconstructions and every "
        f"sequence their posteriors allow."
    ]
    if total <= 0:
        return out
    if mode == "posterior":
        out.append(
            f"{fmt(within)} of those positions ({percent(within, total)}) are the posterior's "
            f"own uncertainty: they would remain with every choice fixed, and only more or "
            f"better data, or a better model of evolution, narrows them."
        )
    ranked = sorted(
        ((label, whole.get(key, 0.0)) for key, label, _ in components
         if key not in ("within", "other")),
        key=lambda kv: -kv[1],
    )  # fmt: skip
    if ranked and ranked[0][1] > 0:
        label, value = ranked[0]
        rest = [f"{lb.lower()} {fmt(v)}" for lb, v in ranked[1:] if v >= 0.05]
        out.append(
            f"{fmt(choices)} ({percent(choices, total)}) come from your choices, most of them "
            f"from the {label.lower()} ({fmt(value)} positions"
            + (f"; then {', '.join(rest)}" if rest else "")
            + f"). That is how far your ancestors would move under another {label.lower()}: "
            f"before ordering one, compare its versions under each {label.lower()} (the "
            f"library's Disputed column does this site by site)."
        )
        key = next(k for k, lb, _ in components if lb == label)
        gap = parts.get("presence", {}).get(key, 0.0)
        if value > 0 and gap / value > 0.5:
            out.append(
                f"Most of the {label.lower()}'s part ({fmt(gap)} of {fmt(value)}) is about "
                f"whether a residue is there at all, not which residue: where the ancestor "
                f"has insertions and deletions, rather than what it is made of."
            )
    if totals:
        ordered = sorted(totals.values())
        most = max(totals, key=lambda n: totals[n])
        typical = ordered[len(ordered) // 2]
        high = ordered[min(len(ordered) - 1, int(0.9 * len(ordered)))]
        out.append(
            f"Ancestors vary: half differ at {fmt(typical)} positions or fewer, nine in ten "
            f"at {fmt(high)} or fewer, and the most uncertain, {most}, at {fmt(totals[most])}."
        )
    return out


@dataclass(frozen=True)
class Tile:
    """One figure of the page's opening summary: a number, what it is, and the explanation
    its "?" opens. `colour` is a component's --viz slot (None: no swatch)."""

    key: str
    label: str
    value: str
    unit: str
    help: str
    colour: int | None = None
    node: str = ""  # the ancestor the tile names, for its link


@dataclass(frozen=True)
class Glance:
    """The page's answer at a glance: the headline number and the few figures behind it,
    the same numbers `takeaways` puts in sentences."""

    total: str
    columns: int
    nodes: int
    share: str  # the total as a share of the alignment's columns
    tiles: list[Tile]


def glance(
    summary: dict[str, Any],
    parts: dict[str, dict[str, float]],
    totals: dict[str, float],
    components: tuple[Any, ...],
    mode: str,
) -> Glance:
    whole = parts["all"]
    total, within = whole["total"], whole.get("within", 0.0)
    choices = total - within
    nodes, columns = summary.get("nodes", len(totals)), summary.get("columns", 0)
    tiles: list[Tile] = []
    if total > 0:
        colours = {key: colour for key, _, colour in components}
        if mode == "posterior":
            tiles.append(
                Tile(
                    "posterior", "The posterior's own", fmt(within),
                    f"positions · {percent(within, total)}",
                    "Uncertainty each reconstruction already admits in its posterior "
                    "probabilities. It would remain with every choice fixed: only more or "
                    "better data, or a better model of evolution, narrows it.",
                    colours.get("within"),
                )
            )  # fmt: skip
        ranked = sorted(
            ((key, label, whole.get(key, 0.0)) for key, label, _ in components
             if key not in ("within", "other")),
            key=lambda item: -item[2],
        )  # fmt: skip
        tiles.append(
            Tile(
                "choices", "Your choices", fmt(choices),
                f"positions · {percent(choices, total)}",
                "How far the ancestors move between this job's aligners, trees, models and "
                "gap placements: "
                + ", ".join(f"{label.lower()} {fmt(value)}" for _, label, value in ranked)
                + ". These describe this job's choices, not the methods in general.",
            )
        )  # fmt: skip
        if ranked and ranked[0][2] > 0:
            key, label, value = ranked[0]
            tiles.append(
                Tile(
                    "largest", f"Mostly the {label.lower()}", fmt(value),
                    f"positions · {percent(value, choices)} of the choices",
                    f"The choice that moves the ancestors most. Before ordering an ancestor, "
                    f"compare its versions under each {label.lower()}: the library's "
                    f"Disputed column does this site by site.",
                    colours.get(key),
                )
            )  # fmt: skip
        presence = parts.get("presence", {}).get("total", 0.0)
        tiles.append(
            Tile(
                "presence", "Residue or gap", fmt(presence),
                f"positions · {percent(presence, total)}",
                "Disagreement over whether the ancestor has a residue at a position at all "
                "(its insertions and deletions), rather than which residue it is. Hatched "
                "in every bar; the solid part is which residue.",
            )
        )  # fmt: skip
    if totals:
        ordered = sorted(totals.values())
        most = max(totals, key=lambda n: totals[n])
        typical = ordered[len(ordered) // 2]
        high = ordered[min(len(ordered) - 1, int(0.9 * len(ordered)))]
        tiles.append(
            Tile(
                "most", f"Most uncertain: {most}", fmt(totals[most]), "positions",
                f"Ancestors vary: half differ at {fmt(typical)} positions or fewer, nine in "
                f"ten at {fmt(high)} or fewer.",
                node=most,
            )
        )  # fmt: skip
    return Glance(fmt(total), columns, nodes, percent(total, columns), tiles)


# --- Uncertainty with depth -------------------------------------------------------------


def nearest_tip(newick: str) -> dict[str, float]:
    """{internal node label: path length to its nearest descendant tip}."""
    from phylobot.trees import load_tree  # dendropy: only when the scatter is drawn

    tree = load_tree(newick)
    out: dict[str, float] = {}
    below: dict[int, float] = {}
    for node in tree.postorder_node_iter():
        children = node.child_nodes()
        if not children:
            below[id(node)] = 0.0
            continue
        below[id(node)] = min((c.edge_length or 0.0) + below[id(c)] for c in children)
        if node.label:
            out[node.label] = below[id(node)]
    return out


@dataclass(frozen=True)
class Scatter:
    dots: list[tuple[float, float, str, str]]  # x, y, tip, href
    trend: str  # path d of the binned medians
    x_ticks: list[tuple[float, str]]
    y_ticks: list[tuple[float, str]]
    x_label: str
    y_label: str
    fill: str
    width: int
    height: int
    left: int
    top: int
    plot_width: int
    plot_height: int
    note: str


SCATTER_W, SCATTER_H = 640, 300
SCATTER_LEFT, SCATTER_TOP, SCATTER_BOTTOM = 58, 12, 46
X_AXES = {"tip": "Distance to the nearest descendant tip (substitutions per site)",
          "below": "Sequences below (log scale)"}  # fmt: skip
Y_SCALES = {"share": "Share of the ancestor's total", "positions": "Positions"}


def ticks(low: float, high: float, count: int = 5) -> list[float]:
    """Round tick values covering low..high."""
    span = high - low or 1.0
    raw = span / count
    step = 10 ** math.floor(math.log10(raw))
    step = next(m * step for m in (1, 2, 5, 10) if m * step >= raw)
    first = math.floor(low / step) * step
    return [
        round(first + i * step, 10) for i in range(count + 3) if first + i * step <= high + 1e-9
    ]


def depth_scatter(
    rows: list[tuple[str, float, int, dict[str, dict[str, float]]]],
    component: str,
    x_axis: str,
    y_scale: str,
    components: tuple[Any, ...],
    hrefs: dict[str, str],
) -> Scatter | None:
    """One dot per ancestor: x its depth (nearest tip, or sequences below), y one part of
    its uncertainty (a share of its total, or positions); a line through the medians of
    eight bins of x shows the trend. `rows`: (node, nearest tip, sequences below, parts)."""
    colours = {key: (label, colour) for key, label, colour in components}
    if component == "total":
        name, fill, y_scale = "Total", INK, "positions"
    elif component in colours:
        name, fill = colours[component][0], FILLS[colours[component][1]]
    else:
        return None
    points = []
    for node, tip, below, p in rows:
        whole = p["all"]
        value = whole["total"] if component == "total" else whole.get(component, 0.0)
        if y_scale == "share":
            if whole["total"] <= 0:
                continue
            value /= whole["total"]
        x = tip if x_axis == "tip" else math.log10(max(below, 1))
        points.append((node, x, value))
    if len(points) < 3:
        return None
    xs = [x for _, x, _ in points]
    x_low, x_high = min(xs), max(xs) or 1.0
    y_high = 1.0 if y_scale == "share" else max(v for _, _, v in points) or 1.0
    x_tick_values = ticks(0 if x_axis == "tip" else x_low, x_high)
    x_low = x_tick_values[0]
    x_high = max(x_high, x_tick_values[-1])
    y_tick_values = ticks(0, y_high)
    y_high = max(y_high, y_tick_values[-1])
    pw, ph = SCATTER_W - SCATTER_LEFT - 16, SCATTER_H - SCATTER_TOP - SCATTER_BOTTOM

    def px(x: float) -> float:
        return round(SCATTER_LEFT + pw * (x - x_low) / ((x_high - x_low) or 1.0), 2)

    def py(y: float) -> float:
        return round(SCATTER_TOP + ph * (1 - y / (y_high or 1.0)), 2)

    unit = "" if y_scale == "share" else " positions"
    dots = [
        (px(x), py(v), f"{node}: {fmt(100 * v) + '%' if y_scale == 'share' else fmt(v)}"
         f"{unit} {name.lower()}", hrefs.get(node, ""))
        for node, x, v in points
    ]  # fmt: skip
    # The trend: medians of eight bins with equal numbers of ancestors.
    ordered = sorted(points, key=lambda t: t[1])
    size = max(len(ordered) // 8, 1)
    trend = []
    for start in range(0, len(ordered), size):
        chunk = ordered[start : start + size]
        if len(chunk) < max(size // 2, 1):
            break
        xm = sorted(x for _, x, _ in chunk)[len(chunk) // 2]
        ym = sorted(v for _, _, v in chunk)[len(chunk) // 2]
        trend.append(f"{px(xm)},{py(ym)}")
    x_label = X_AXES[x_axis]
    x_ticks = [(px(t), fmt_tick(t) if x_axis == "tip" else f"{10**t:.0f}") for t in x_tick_values]
    y_ticks = [(py(t), f"{round(100 * t)}%" if y_scale == "share" else fmt_tick(t))
               for t in y_tick_values]  # fmt: skip
    y_label = (f"{name}: share of the total" if y_scale == "share"
               else f"{name}: positions")  # fmt: skip
    return Scatter(
        dots=dots,
        trend="M " + " L ".join(trend) if len(trend) > 1 else "",
        x_ticks=x_ticks,
        y_ticks=y_ticks,
        x_label=x_label,
        y_label=y_label,
        fill=fill,
        width=SCATTER_W,
        height=SCATTER_H,
        left=SCATTER_LEFT,
        top=SCATTER_TOP,
        plot_width=pw,
        plot_height=ph,
        note=trend_note(points, name, y_scale),
    )


def fmt_tick(value: float) -> str:
    return f"{value:g}"


def trend_note(points: list[tuple[str, float, float]], name: str, y_scale: str) -> str:
    """The trend in words: the median in the shallowest and the deepest third."""
    ordered = sorted(points, key=lambda t: t[1])
    third = max(len(ordered) // 3, 1)

    def median(chunk: list[tuple[str, float, float]]) -> float:
        values = sorted(v for _, _, v in chunk)
        return values[len(values) // 2]

    near, deep = median(ordered[:third]), median(ordered[-third:])
    show = (lambda v: f"{round(100 * v)}%") if y_scale == "share" else fmt
    what = "of the total" if y_scale == "share" else "positions"
    return (
        f"{name}: a median {show(near)} {what} in the third of ancestors nearest the tips, "
        f"{show(deep)} in the deepest third."
    )
