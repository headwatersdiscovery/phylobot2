"""The Uncertainty page (`/jobs/<id>/uncertainty/`).

Reads an uncertainty run's files (`uncertainty_run`), never computes them: the run is a
follow-up job, and its files never change, so each one is parsed once per process. The
figures are server SVG, laid out here, drawn by the template: a headline bar, a bar per
ancestor, and for one ancestor the per-column track, the sites that vary most with what
each choice says there, and the differences between every pair of cells.

Components keep one colour each wherever they appear (the --viz slots, validated for
colour blindness): posterior 1, alignment 2, tree 3, model 4, gap placement 5, other 6.
The presence part of each is hatched, the residue part solid (§1.3), so the split never
rests on colour alone. Colours are SVG attributes (the light slots), not classes, so a
downloaded figure looks as shown (`plots.js`); dark.css remaps each one.
"""

from __future__ import annotations

import functools
import json
import logging
import re
from dataclasses import dataclass, field
from typing import Any
from urllib.parse import urlencode

from django.http import Http404, HttpRequest, HttpResponse, QueryDict
from django.shortcuts import get_object_or_404, render
from django.urls import reverse

from phylobot import ancestors, joined, storage, treeview
from phylobot import uncertainty as u
from phylobot import uncertainty_overview as overview
from phylobot import uncertainty_run as run_mod
from phylobot import uncertainty_submit as submit_mod
from phylobot import uncertainty_tables as tables
from phylobot.gaps import MASKING
from phylobot.models import Job

logger = logging.getLogger(__name__)

COMPONENTS = (
    ("within", "Posterior", 1),
    ("alignment", "Alignment", 2),
    ("tree", "Tree", 3),
    ("model", "Model", 4),
    ("gaps", "Gap placement", 5),
    ("other", "Other", 6),
)
LABELS = {key: label for key, label, _ in COMPONENTS}
# The --viz slots of site.css (light); dark.css remaps these exact values in SVG.
FILLS = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a", 4: "#eda100", 5: "#e87ba4", 6: "#7d3c98"}
INK, MUTED, GRID, WEAK = "#1d2330", "#5d6675", "#d8dde5", "#b3bac5"


def fill(colour: int, hatched: bool, prefix: str) -> str:
    return f"url(#{prefix}-hatch-{colour})" if hatched else FILLS[colour]


SORTS = ("tree", "total")
SCALES = ("residues", "share")
TOP_SITES = 20
FILES_CACHED = 256


@functools.lru_cache(maxsize=FILES_CACHED)
def _text(key: str) -> str:
    """A run's file, read once (a finished job's files never change)."""
    return storage.read_bytes(key).decode()


def run_file(run: Job, relative: str) -> str | None:
    try:
        return _text(f"{run.storage_prefix}results/{run_mod.OUT}/{relative}")
    except Exception:
        return None


def rows_of(text: str | None) -> list[dict[str, str]]:
    if not text:
        return []
    header, *lines = text.splitlines()
    names = header.split("\t")
    return [dict(zip(names, line.split("\t"), strict=True)) for line in lines if line]


def fmt(value: float) -> str:
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


# --- Stacked bars ------------------------------------------------------------------


@dataclass(frozen=True)
class Segment:
    x: float
    width: float
    fill: str  # a slot's colour, or its hatch pattern (the presence part)
    title: str


@dataclass
class Bar:
    label: str
    total: float
    segments: list[Segment]
    href: str = ""
    present: str = ""  # "k/n" when the node is missing from some combinations
    note: str = ""
    y: float = 0.0
    end: float = 0.0  # where the bar stops, for its total
    values: dict[str, float] = field(default_factory=dict)  # component -> residues

    @property
    def text(self) -> str:
        return fmt(self.total)


GAP_PX = 2.0  # the surface gap between neighbouring fills


def segments(
    values: dict[str, dict[str, float]], scale: float, left: float, share: bool, prefix: str
) -> list[Segment]:
    """One bar's segments: per component, its residue part then its presence part.
    `values` is {part: {component: residues}} for "presence" and "residue"."""
    total = sum(values["presence"][k] + values["residue"][k] for k, _, _ in COMPONENTS)
    unit = (1 / total if total > 0 else 0) if share else 1.0
    out = []
    x = left
    for key, label, colour in COMPONENTS:
        for part, hatched in (("residue", False), ("presence", True)):
            amount = values[part][key] * unit
            width = amount * scale
            if width < 0.5:
                continue
            what = "residue identity" if part == "residue" else "residue or gap"
            value = f"{100 * amount:.1f}%" if share else f"{fmt(amount)} positions"
            drawn = max(width - GAP_PX, 0.5)
            title = f"{label}, {what}: {value}"
            out.append(Segment(round(x, 2), round(drawn, 2), fill(colour, hatched, prefix), title))
            x += width
    return out


@dataclass
class BarChart:
    bars: list[Bar]
    width: int
    height: int
    left: int
    plot: int
    ticks: list[tuple[float, str]]
    row: int
    patterns: str  # a prefix for the hatch pattern ids, unique per chart
    label: str = ""

    @property
    def baseline(self) -> int:
        return self.height - 22


def tick(fraction: float, top: float, share: bool) -> str:
    if share:
        return f"{100 * fraction:.0f}%"
    value = top * fraction
    return f"{value:g}" if top >= 10 else f"{value:.2g}"


def nice_max(value: float) -> float:
    if value <= 0:
        return 1.0
    step = 10 ** max(len(str(int(value))) - 1, 0)
    for m in (1, 2, 2.5, 3, 4, 5, 6, 8, 10):
        if value <= m * step:
            return float(m * step)
    return float(10 * step)


def bar_chart(
    rows: list[tuple[str, dict[str, dict[str, float]], str, str]],
    share: bool,
    patterns: str,
    left: int = 150,
    plot: int = 620,
    row: int = 26,
) -> BarChart:
    """rows: (label, {part: {component: value}}, href, present)."""
    totals = [
        sum(v["presence"][k] + v["residue"][k] for k, _, _ in COMPONENTS) for _, v, _, _ in rows
    ]
    top = 1.0 if share else nice_max(max(totals, default=0))
    scale = plot / top
    bars = []
    for i, ((label, values, href, present), total) in enumerate(zip(rows, totals, strict=True)):
        drawn = segments(values, scale, left, share, patterns)
        end = drawn[-1].x + drawn[-1].width if drawn else left
        bar = Bar(label, total, drawn, href, present, y=i * row + 4, end=round(end, 2))
        bar.values = {k: values["presence"][k] + values["residue"][k] for k, _, _ in COMPONENTS}
        bars.append(bar)
    ticks = [(left + plot * f, tick(f, top, share)) for f in (0, 0.25, 0.5, 0.75, 1)]
    height = len(rows) * row + 30
    return BarChart(bars, left + plot + 70, height, left, plot, ticks, row, patterns)


# --- The run's numbers -------------------------------------------------------------


def parts_of(rows: list[dict[str, str]]) -> dict[str, dict[str, float]]:
    """{part: {component: value}} from the rows of one node and mode (nodes.tsv), with
    each part's total."""
    keys = ["total", *(k for k, _, _ in COMPONENTS)]
    return {r["part"]: {k: float(r[k]) for k in keys} for r in rows}


def job_parts(summary: dict[str, Any], mode: str) -> dict[str, dict[str, float]]:
    keys = ["total", *(k for k, _, _ in COMPONENTS)]
    job = summary.get("job", {}).get(mode, {})
    return {part: {k: float(job.get(part, {}).get(k, 0.0)) for k in keys} for part in run_mod.PARTS}


def sentence(parts: dict[str, dict[str, float]], mode: str, factors: list[dict[str, Any]]) -> str:
    """The headline in words (§7): descriptive of this job, never a claim about methods."""
    total = parts["all"]["total"]
    if total <= 0:
        return "Every reconstruction of these ancestors is the same at every position."
    varied = {f["key"] for f in factors if f.get("varied")}
    pieces = []
    for key, label, _ in COMPONENTS:
        value = parts["all"][key]
        if key == "other" and value < 0.005:
            continue
        if key not in ("within", "other") and key not in varied:
            continue
        name = "the posterior itself" if key == "within" else f"the {label.lower()}"
        pieces.append(f"{fmt(value)} from {name}")
    what = (
        "the most probable sequences of two reconstructions picked at random"
        if mode == "map"
        else "a sample from each of two reconstructions picked at random"
    )
    return (
        f"On average over the ancestors, {what} differ at "
        f"{fmt(total)} positions: " + ", ".join(pieces) + "."
    )


@dataclass(frozen=True)
class Track:
    """The per-column stacked track of one node: one path per component and part."""

    width: int
    height: int
    paths: list[tuple[str, str, str]]  # fill, d, title
    agreement: str  # path of the low-agreement marks under the track
    columns: int
    ticks: list[tuple[float, str]]


def site_parts(row: dict[str, str]) -> dict[str, dict[str, float]]:
    """{part: {value: number}} from one sites/<n>.tsv row, "all" the sum of the two."""
    keys = ["total", *(k for k, _, _ in COMPONENTS)]
    parts = {part: {k: float(row[f"{part}_{k}"]) for k in keys} for part in ("presence", "residue")}
    parts["all"] = {k: parts["presence"][k] + parts["residue"][k] for k in keys}
    return parts


TRACK_WIDTH = 880  # drawing units: the figure scales to the card, so its text must not


def column_step(columns: int, most: int = 10) -> int:
    """A 1-2-5 tick step giving at most `most` labels along `columns`."""
    step = 1
    while columns / step > most:
        for m in (2, 5, 10):
            if columns / (step * m) <= most:
                return step * m
        step *= 10
    return step


def track(rows: list[dict[str, str]], mode: str) -> Track:
    """Laid out at TRACK_WIDTH whatever the frame's length (a column is TRACK_WIDTH / n
    units), since `.plot > svg` scales the figure to the card: drawn one unit per column,
    a short frame was blown up several times over, tick labels with it."""
    by_column: dict[int, dict[str, dict[str, float]]] = {}
    agreement: dict[int, float] = {}
    for r in rows:
        if r["mode"] != mode:
            continue
        column = int(r["column"])
        by_column[column] = site_parts(r)
        agreement[column] = float(r["agreement"])
    columns = max(by_column, default=0)
    unit = TRACK_WIDTH / max(columns, 1)
    height = 120
    base = 100.0
    paths = []
    stacks = dict.fromkeys(by_column, 0.0)
    for key, label, colour in COMPONENTS:
        for part, hatched in (("residue", False), ("presence", True)):
            d = []
            for column, parts in sorted(by_column.items()):
                value = parts.get(part, {}).get(key, 0.0)
                if value * base < 0.05:
                    continue
                y = base - (stacks[column] + value) * base
                x = (column - 1) * unit
                d.append(f"M{x:.2f} {y:.2f}h{unit:.3f}v{value * base:.2f}h-{unit:.3f}z")
                stacks[column] += value
            if d:
                part_label = "presence" if hatched else "residue"
                paths.append((fill(colour, hatched, "tr"), "".join(d), f"{label} ({part_label})"))
    weak = "".join(
        f"M{(c - 1) * unit:.2f} {base + 6:.1f}h{unit:.3f}v{(1 - a) * 12:.2f}h-{unit:.3f}z"
        for c, a in sorted(agreement.items())
        if a < 1
    )
    step = column_step(columns)
    ticks = [((c - 0.5) * unit, str(c)) for c in range(step, columns + 1, step)]
    return Track(TRACK_WIDTH, height, paths, weak, columns, ticks)


@dataclass
class SiteRow:
    column: int
    total: float
    within: float
    between: float
    shares: list[tuple[str, int, float]]  # label, colour, value
    says: list[tuple[str, list[str]]]  # factor label -> "level: state (p)"
    agreement: float
    viewer: str = ""


def top_sites(
    site_rows: list[dict[str, str]],
    drivers: list[dict[str, str]],
    mode: str,
    level_labels: dict[tuple[str, str], str],
    varied: list[str],
) -> list[SiteRow]:
    # (column, agreement, the whole site's values)
    alls = [
        (int(r["column"]), float(r["agreement"]), site_parts(r)["all"])
        for r in site_rows
        if r["mode"] == mode
    ]
    ranked = sorted(alls, key=lambda a: -(a[2]["total"] - a[2]["within"]))[:TOP_SITES]
    said: dict[int, dict[str, list[str]]] = {}
    for d in drivers:
        said.setdefault(int(d["column"]), {}).setdefault(d["factor"], []).append(
            f"{level_labels.get((d['factor'], d['level']), d['level'])}: "
            f"{d['state']} ({float(d['p']):.2f})"
        )
    out = []
    for column, agreement, values in ranked:
        between = values["total"] - values["within"]
        if between < 0.005:
            break
        shares = [
            (label, colour, values[key])
            for key, label, colour in COMPONENTS
            if key != "within" and values[key] >= 0.005
        ]
        here = said.get(column, {})
        says = [(LABELS[f], here[f]) for f in varied if here.get(f)]
        out.append(
            SiteRow(column, values["total"], values["within"], between, shares, says, agreement)
        )
    return out


# The heatmap's sequential palettes: near (0) to far (the node's largest pair), one hue
# or a perceptual ramp, light to dark, never a rainbow. Interpolated colours are not
# remapped by dark.css, so each ramp is chosen to read on both backgrounds.
Stop = tuple[float, tuple[int, int, int]]
HEAT_PALETTES: dict[str, tuple[str, list[Stop]]] = {
    "blues": ("Blues", [(0.0, (239, 243, 255)), (0.5, (107, 174, 214)), (1.0, (8, 69, 148))]),
    "viridis": ("Viridis (colour-blind safe)", [
        (0.0, (253, 231, 37)), (0.25, (94, 201, 98)), (0.5, (33, 145, 140)),
        (0.75, (59, 82, 139)), (1.0, (68, 1, 84)),
    ]),
    "cividis": ("Cividis (colour-blind safe)", [
        (0.0, (255, 234, 70)), (0.5, (124, 123, 120)), (1.0, (0, 32, 77)),
    ]),
    "oranges": ("Oranges", [(0.0, (254, 237, 222)), (0.5, (253, 141, 60)), (1.0, (166, 54, 3))]),
    "greys": ("Greys", [(0.0, (240, 240, 240)), (0.5, (150, 150, 150)), (1.0, (37, 37, 37))]),
}  # fmt: skip
DEFAULT_HEAT = "blues"
LEGEND_STOPS = 11  # gradient stops drawn in the legend


def heat_colour(fraction: float, palette: str) -> str:
    """The palette's colour at a fraction of the largest distance, as #rrggbb."""
    stops = HEAT_PALETTES.get(palette, HEAT_PALETTES[DEFAULT_HEAT])[1]
    fraction = min(max(fraction, 0.0), 1.0)
    i = next(i for i in range(1, len(stops)) if fraction <= stops[i][0])
    (low, a), (high, b) = stops[i - 1], stops[i]
    t = (fraction - low) / (high - low)
    return "#" + "".join(f"{round(x + (y - x) * t):02x}" for x, y in zip(a, b, strict=True))


@dataclass(frozen=True)
class Heatmap:
    labels: list[tuple[float, str]]  # y, text
    cells: list[tuple[float, float, str, str]]  # x, y, fill, title
    width: int
    height: int
    left: int
    step: int
    legend_y: float  # the colour key, under the grid
    legend_width: float
    gradient: list[tuple[str, str]]  # offset (%), colour
    ticks: list[tuple[float, str]]  # x, label


def heatmap(
    pairs: list[dict[str, str]],
    mode: str,
    cells: list[dict[str, Any]],
    factors: list[dict[str, Any]],
    palette: str = DEFAULT_HEAT,
) -> Heatmap | None:
    """The pairs of one node (pairs/<n>.tsv: cells by their position in summary.json),
    coloured from 0 to the node's largest distance, with that scale drawn under it."""
    keys = [c["key"] for c in cells]
    if len(keys) < 2:
        return None
    values: dict[tuple[int, int], float] = {}
    for p in pairs:
        a, b = int(p["a"]), int(p["b"])
        if p["mode"] == mode and 0 <= a < len(keys) and 0 <= b < len(keys):
            values[(a, b)] = float(p["differences"])
    if not values:
        return None
    top = max(values.values()) or 1.0
    step = 30 if len(keys) <= 16 else 22 if len(keys) <= 32 else 12
    labels = [cell_label(c, factors) for c in cells]
    # 11px text; at least room for the colour key's label, which sits in the same margin
    left = max(12 + round(6.4 * max(len(label) for label in labels)), 190)
    drawn = [
        (float(left + j * step), float(8 + i * step), heat_colour(v / top, palette),
         f"{labels[i]} vs {labels[j]}: {fmt(v)} positions")
        for (i, j), v in sorted(values.items())
    ]  # fmt: skip
    placed = [(8 + i * step + step / 2 + 4, text) for i, text in enumerate(labels)]
    size = len(keys) * step
    legend_y = 8 + size + 18
    legend_width = float(max(size, 220))
    gradient = [
        (f"{100 * k / (LEGEND_STOPS - 1):.0f}", heat_colour(k / (LEGEND_STOPS - 1), palette))
        for k in range(LEGEND_STOPS)
    ]
    ticks = [(left + legend_width * f, fmt(top * f) if f else "0") for f in (0, 0.5, 1)]
    width = int(left + max(size, legend_width) + 30)
    return Heatmap(
        placed,
        drawn,
        width,
        int(legend_y + 52),
        left,
        step,
        legend_y,
        legend_width,
        gradient,
        ticks,
    )


def cell_label(cell: dict[str, Any], factors: list[dict[str, Any]]) -> str:
    """A cell named by the choices the job varied (the rest are the same in every cell)."""
    names = []
    for i, factor in enumerate(factors):
        if not factor.get("varied"):
            continue
        labels = {lv["key"]: lv["label"] for lv in factor.get("levels", [])}
        level = cell["levels"][i]
        names.append(labels.get(level, level))
    return " · ".join(names) or str(cell.get("label", cell["combination"]))


# --- The view ----------------------------------------------------------------------


def heat_of(params: QueryDict) -> str:
    asked = params.get("heat", DEFAULT_HEAT)
    return asked if asked in HEAT_PALETTES else DEFAULT_HEAT


def downloads(job: Job, run: Job, mode: str, table: str, node: str = "") -> list[tuple[str, str]]:
    """(label, URL) for a table's CSV and Excel files, for the run and mode shown (the
    run by its id, so a run finishing meanwhile does not change what is downloaded)."""
    base = reverse("uncertainty_download", args=[job.id])
    out = []
    for label, kind in (("CSV", "csv"), ("Excel", "xlsx")):
        fields = {"run": str(run.id), "mode": mode, "table": table, "format": kind}
        if node:
            fields["node"] = node
        out.append((label, f"{base}?{urlencode(fields)}"))
    return out


def chosen_run(runs: list[Job], asked: str) -> Job | None:
    by_id = {str(r.id): r for r in runs}
    if asked in by_id:
        return by_id[asked]
    return next((r for r in runs if r.status == Job.Status.SUCCEEDED), runs[0] if runs else None)


def query(base: QueryDict, **changes: str) -> str:
    q = base.copy()
    for key, value in changes.items():
        if value:
            q[key] = value
        else:
            q.pop(key, None)
    return "?" + urlencode(sorted((k, v) for k, v in q.items()))


def form_values(run: Job | None, offered: submit_mod.Choices) -> dict[str, Any]:
    """The recompute form, filled from the run shown (or the defaults)."""
    settings = (
        run_mod.Settings.parse((run.params.get("uncertainty_of") or {}).get("settings") or {})
        if run
        else run_mod.Settings()
    )
    versions = (
        list(settings.versions)
        or [v for v in offered.versions if v in MASKING]
        or [ancestors.ORIGINAL]
    )
    factors = []
    for factor in u.FACTORS:
        levels = offered.levels.get(factor, [])
        factors.append({
            "key": factor,
            "label": LABELS[factor],
            "levels": [
                {"name": f"weight-{factor}-{i}", "label": lv.label,
                 "value": f"{settings.multipliers.get(factor, {}).get(lv.key, 1.0):g}"}
                for i, lv in enumerate(levels)
            ],
        })  # fmt: skip
    return {
        "weighting": settings.weighting,
        "fills": FILLS,
        "reference": settings.reference or offered.reference,
        "versions": versions,
        "factors": factors,
        "version_choices": [
            (v, submit_mod.version_label(v), v in run_mod.OPT_IN) for v in offered.versions
        ],
    }


def page(request: HttpRequest, job_id: str) -> HttpResponse:
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if job.status != Job.Status.SUCCEEDED or job.params.get("uncertainty_only"):
        raise Http404
    files = submit_mod.published(job)
    offered = submit_mod.choices(files)
    runs = submit_mod.runs_of(job)
    run = chosen_run(runs, request.GET.get("run", ""))
    mode = request.GET.get("mode", "posterior")
    mode = mode if mode in run_mod.MODES else "posterior"
    sort = request.GET.get("sort", "tree")
    sort = sort if sort in SORTS else "tree"
    scale = request.GET.get("scale", "residues")
    scale = scale if scale in SCALES else "residues"
    share = scale == "share"
    context: dict[str, Any] = {
        "job": job, "runs": runs, "run": run, "mode": mode, "sort": sort, "scale": scale,
        "modes": [(m, run_mod.MODE_LABELS[m]) for m in run_mod.MODES],
        "offered": offered, "form": form_values(run, offered), "components": COMPONENTS,
        "can_run": bool(offered.combinations),
        "fills": sorted(FILLS.items()),
        "legend": [(label, FILLS[colour]) for _, label, colour in COMPONENTS],
    }  # fmt: skip
    summary_text = (
        run_file(run, "summary.json") if run and run.status == Job.Status.SUCCEEDED else None
    )
    if run is None or summary_text is None:
        return render(request, "phylobot/uncertainty.html", context)
    summary = json.loads(summary_text)
    factors = summary.get("factors", [])
    parts = job_parts(summary, mode)
    context |= {
        "summary": summary,
        "factors": factors,
        "varied": [f for f in factors if f.get("varied")],
        "not_varied": [f for f in factors if not f.get("varied")],
        "sentence": sentence(parts, mode, factors),
        "headline": bar_chart([("Mean ancestor", parts, "", "")], share, "hl", row=34),
        "headline_rows": [
            (label, parts[part]["total"], [parts[part][k] for k, _, _ in COMPONENTS])
            for part, label in (
                ("all", "All"),
                ("presence", "Residue or gap"),
                ("residue", "Which residue"),
            )
        ],
        "point_estimates": summary.get("point_estimates", []),
        "explainer": overview.explainer(summary, parts, COMPONENTS),
        "map_reset": query(request.GET, show="", map_set="", colour_by=""),
        "recon_map": overview.reconstruction_map(
            summary,
            mode,
            parts,
            COMPONENTS,
            request.GET.get("colour_by", ""),
            overview.map_shown(request.GET),
        ),
        "weighting": (run.params.get("uncertainty_of") or {}).get("settings", {}),
    }
    by_node: dict[str, list[dict[str, str]]] = {}
    present: dict[str, tuple[str, str, str]] = {}
    for r in rows_of(run_file(run, "nodes.tsv")):
        if r["mode"] == mode:
            by_node.setdefault(r["node"], []).append(r)
            present[r["node"]] = (r["present"], r["of"], r["leaves"])
    order = list(by_node)
    node_parts = {n: parts_of(rows) for n, rows in by_node.items()}
    if sort == "total":
        order.sort(key=lambda n: -node_parts[n]["all"]["total"])
    node = request.GET.get("node", "")
    rows = []
    for n in order:
        k, of, _ = present[n]
        rows.append(
            (n, node_parts[n], query(request.GET, node=n) + "#node", f"{k}/{of}" if k != of else "")
        )
    chart = bar_chart(rows, share, "nd")
    for bar in chart.bars:
        bar.note = f"{present[bar.label][2]} sequences below"
    context["nodes"] = chart
    tree = pie_tree(job, summary, request.GET, node_parts)
    newick = tree.pop("newick", "")
    context |= tree
    context |= depth_view(request.GET, newick, order, node_parts, present)
    context["takeaways"] = overview.takeaways(
        summary, parts, {n: node_parts[n]["all"]["total"] for n in order}, COMPONENTS, mode
    )
    glance = overview.glance(
        summary, parts, {n: node_parts[n]["all"]["total"] for n in order}, COMPONENTS, mode
    )
    context["glance"] = glance
    # Each tile with its swatch colour and, where it names an ancestor, that ancestor's link.
    context["glance_tiles"] = [
        (
            tile,
            FILLS.get(tile.colour or 0, ""),
            query(request.GET, node=tile.node) + "#node" if tile.node else "",
        )
        for tile in glance.tiles
    ]
    context["ancestor_downloads"] = downloads(job, run, mode, "ancestors")
    context["node_table"] = [
        {"node": n, "present": present[n][0], "of": present[n][1], "leaves": present[n][2],
         "values": [(key, node_parts[n]["all"][key]) for key, _, _ in COMPONENTS],
         "total": node_parts[n]["all"]["total"], "presence": node_parts[n]["presence"]["total"],
         "residue": node_parts[n]["residue"]["total"], "href": query(request.GET, node=n) + "#node"}
        for n in order
    ]  # fmt: skip
    if node in by_node:
        context |= node_detail(run, summary, node, mode, request.GET, job)
    return render(request, "phylobot/uncertainty.html", context)


def pie_tree(
    job: Job, summary: dict[str, Any], params: QueryDict, node_parts: dict[str, Any]
) -> dict[str, Any]:
    """The reference combination's tree (as the tree viewer draws it, in the asked branch
    lengths) with a pie at each ancestor (uncertainty_overview.pies); absent without it."""
    lengths = params.get("lengths", treeview.DEFAULT_MODE)
    lengths = lengths if lengths in treeview.BRANCH_MODES else treeview.DEFAULT_MODE
    reference = summary.get("reference", "")
    files = joined.cached(job)
    name = next((f"{reference}/{n}" for n in ("tree.rooted.nwk", "tree.nwk")
                 if f"{reference}/{n}" in files.keys), "")  # fmt: skip
    if not reference or not name:
        return {}
    spacing = overview.slider(params.get("spacing"), overview.TREE_SPACING)
    pie_size = overview.slider(params.get("pie_size"), overview.PIE_SIZES)
    try:
        newick = files.read(name)
        drawing = treeview.layout(
            newick, lengths, ancestors=set(node_parts),
            outgroup=set(job.params.get("outgroup") or []),
        )  # fmt: skip
    except Exception:
        logger.exception("job %s: could not draw %s", job.id, name)
        return {}
    drawing = overview.stretch(drawing, spacing, treeview.MARGIN["top"], treeview.MARGIN["bottom"])
    places = {n.name: (n.x, n.y) for n in drawing.nodes if n.name}
    hrefs = {n: query(params, node=n) + "#node" for n in node_parts}
    return {
        "newick": newick,
        "pie_tree": drawing,
        "pies": overview.pies(places, node_parts, hrefs, COMPONENTS, "pt", pie_size),
        "tree_lengths": lengths,
        "tree_modes": list(treeview.BRANCH_MODES.items()),
        "spacing": f"{spacing:g}",
        "pie_size": f"{pie_size:g}",
        "spacing_bounds": overview.TREE_SPACING,
        "pie_bounds": overview.PIE_SIZES,
    }


def depth_view(
    params: QueryDict,
    newick: str,
    order: list[str],
    node_parts: dict[str, Any],
    present: dict[str, tuple[str, str, str]],
) -> dict[str, Any]:
    """The Uncertainty with Depth scatter: each ancestor's depth against one part."""
    if not newick:
        return {}
    try:
        tips = overview.nearest_tip(newick)
    except Exception:
        logger.exception("could not measure depths")
        return {}
    component = params.get("depth_y", "within")
    keys = [k for k, _, _ in COMPONENTS] + ["total"]
    component = component if component in keys else "within"
    x_axis = params.get("depth_x", "tip")
    x_axis = x_axis if x_axis in overview.X_AXES else "tip"
    y_scale = params.get("depth_scale", "share")
    y_scale = y_scale if y_scale in overview.Y_SCALES else "share"
    rows = [
        (n, tips[n], int(present[n][2] or 0), node_parts[n])
        for n in order if n in tips
    ]  # fmt: skip
    hrefs = {n: query(params, node=n) + "#node" for n in order}
    scatter = overview.depth_scatter(rows, component, x_axis, y_scale, COMPONENTS, hrefs)
    if scatter is None:
        return {}
    return {
        "scatter": scatter,
        "depth_y": component,
        "depth_x": x_axis,
        "depth_scale": y_scale,
        "depth_components": [(k, label) for k, label, _ in COMPONENTS if k != "other"]
        + [("total", "Total")],
        "depth_axes": list(overview.X_AXES.items()),
        "depth_scales": list(overview.Y_SCALES.items()),
    }


def node_detail(
    run: Job, summary: dict[str, Any], node: str, mode: str, params: QueryDict, job: Job
) -> dict[str, Any]:
    index = ancestors.split_index(run_file(run, f"sites/{ancestors.SPLIT_INDEX}") or "")
    name = index.get(node)
    if name is None:
        return {}
    site_rows = rows_of(run_file(run, f"sites/{name}"))
    drivers = rows_of(run_file(run, f"drivers/{name}"))
    pairs = rows_of(run_file(run, f"pairs/{name}"))
    level_labels = {
        (f["key"], lv["key"]): lv["label"]
        for f in summary.get("factors", [])
        for lv in f.get("levels", [])
    }
    varied = [f["key"] for f in summary.get("factors", []) if f.get("varied")]
    sites = top_sites(site_rows, drivers, mode, level_labels, varied)
    viewer = reverse("ancestor_viewer", args=[job.id])
    reference = summary.get("reference", "")
    for s in sites:
        s.viewer = f"{viewer}?{urlencode({'combination': reference, 'node': node})}#site-{s.column}"
    return {
        "node": node,
        "track": track(site_rows, mode),
        "sites": sites,
        "heatmap": heatmap(
            pairs, mode, summary.get("cells", []), summary.get("factors", []), heat_of(params)
        ),
        "heat": heat_of(params),
        "heat_palettes": [(k, label) for k, (label, _) in HEAT_PALETTES.items()],
        "downloads": downloads(job, run, mode, "sites", node),
        "viewer_url": f"{viewer}?{urlencode({'combination': reference, 'node': node})}",
        "close_url": query(params, node=""),
    }


XLSX = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"


def download(request: HttpRequest, job_id: str) -> HttpResponse:
    """A run's table as CSV or Excel: ?table=sites&node=Node3 (every column of the
    reference frame, not only the ones the page lists) or ?table=ancestors; ?run= and
    ?mode= as on the page."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if job.status != Job.Status.SUCCEEDED or job.params.get("uncertainty_only"):
        raise Http404
    run = chosen_run(submit_mod.runs_of(job), request.GET.get("run", ""))
    summary_text = (
        run_file(run, "summary.json") if run and run.status == Job.Status.SUCCEEDED else None
    )
    if run is None or summary_text is None:
        raise Http404
    summary = json.loads(summary_text)
    mode = request.GET.get("mode", "posterior")
    mode = mode if mode in run_mod.MODES else "posterior"
    kind = request.GET.get("format", "csv")
    what = request.GET.get("table", "")
    if kind not in ("csv", "xlsx") or what not in ("sites", "ancestors"):
        raise Http404
    if what == "sites":
        node = request.GET.get("node", "")
        name = ancestors.split_index(run_file(run, f"sites/{ancestors.SPLIT_INDEX}") or "").get(
            node
        )
        if name is None:
            raise Http404
        table = tables.site_table(
            rows_of(run_file(run, f"sites/{name}")),
            rows_of(run_file(run, f"drivers/{name}")),
            mode,
            summary.get("factors", []),
        )
        stem = f"{node}.sites"
    else:
        table = tables.ancestor_table(rows_of(run_file(run, "nodes.tsv")), mode)
        stem = "ancestors"
    stem = re.sub(r"[^A-Za-z0-9._+-]", "_", f"uncertainty.{mode}.{stem}")
    if kind == "csv":
        response = HttpResponse(tables.as_csv(table), content_type="text/csv; charset=utf-8")
    else:
        settings = (run.params.get("uncertainty_of") or {}).get("settings") or {}
        notes = tables.about(job.name or str(job.id), summary, mode, settings)
        if what == "sites":
            notes.insert(1, ("Ancestor", request.GET.get("node", "")))
        response = HttpResponse(tables.as_xlsx(table, notes), content_type=XLSX)
    response["Content-Disposition"] = f'attachment; filename="{stem}.{kind}"'
    return response
