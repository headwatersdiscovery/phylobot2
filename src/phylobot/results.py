"""Read a finished job's results for display: one results table per alignment.

Each alignment gets a table with one row per combination (tree source x model): its
final-tree fit statistics, BIC rank, and links to its tree and ancestors. ModelFinder's
candidate ranking is shown beneath, as a secondary table.

Colours are relative within each table: blue = better, red = worse, pale = middling.
Information criteria (and log-likelihood, as 2*dlogL) are coloured by their difference
from the table's best row, capped at DELTA_CAP: differences of 10 or more mean
essentially no support, and a cap keeps the contrast where decisions are made instead
of spreading it over the thousands of units a long ModelFinder table can span.
Weights are coloured by their share of the table's largest weight, so a solid runner-up
is not painted as bad just because support is split. Columns with no better/worse
direction (parameter counts) are not coloured.
"""

from collections.abc import Collection
from dataclasses import dataclass, field
from typing import Any

from phylobot import catalog

DELTA_CAP = 10.0
MODELFINDER_SHOWN = 15  # rows shown before "show all"
ICON_LINKS = frozenset({"view", "download"})  # links a Cell may show as icons (_icon.html)
# Rows a results table shows before it scrolls inside its own box (its header stays put):
# a 285-combination job has ~60 models per alignment.
SCROLL_ROWS = 20

# Column -> how "better" is judged: "max"/"min" (difference-from-best, capped),
# "weight" (share of the table's largest weight), or None (not coloured).
DIRECTIONS: dict[str, str | None] = {
    "log_likelihood": "max",
    "aic": "min",
    "aicc": "min",
    "bic": "min",
    "delta_bic": "min",
    "w_aic": "weight",
    "w_aicc": "weight",
    "w_bic": "weight",
    "free_parameters": None,
}

WORSE, MIDDLE, BETTER = (214, 96, 77), (247, 247, 247), (67, 147, 195)
TINT = 0.5  # blend with white so text stays readable


@dataclass
class Cell:
    text: str
    style: str = ""
    href: str = ""
    label: bool = False  # a text column (left-aligned), not a number
    links: list[tuple[str, str]] = field(default_factory=list)  # (text, href), shown after text
    badge: str = ""  # e.g. a rooting warning
    title: str = ""  # a tooltip for the cell
    icons: bool = False  # "view" and "download" links shown as icons, their words tooltips

    @property
    def link_items(self) -> list[tuple[str, str, bool, str]]:
        """(text, href, as an icon, separator before it) for the template: " · " between
        words, a space beside an icon."""
        out: list[tuple[str, str, bool, str]] = []
        for i, (text, href) in enumerate(self.links):
            icon = self.icons and text in ICON_LINKS
            before = out[-1][2] if out else False
            sep = "" if i == 0 else (" " if icon or before else " · ")
            out.append((text, href, icon, sep))
        return out


@dataclass
class Table:
    title: str
    subtitle: str
    columns: list[str]
    rows: list[list[Cell]]
    hidden_rows: list[list[Cell]] = field(default_factory=list)
    sortable: frozenset[str] = frozenset()  # columns the reader can sort by (sort.js)
    # The table as a CSV or Excel download gives it (`results_export`): values as
    # model_comparison.tsv has them, unrounded, and the combination's key; no links.
    export_columns: list[str] = field(default_factory=list)
    export_rows: list[list[str]] = field(default_factory=list)

    @property
    def total(self) -> int:
        return len(self.rows) + len(self.hidden_rows)

    @property
    def scrolls(self) -> bool:
        return len(self.rows) > SCROLL_ROWS


def read_tsv(text: str) -> list[dict[str, str]]:
    header, *lines = [line for line in text.splitlines() if line.strip()]
    keys = header.split("\t")
    return [dict(zip(keys, line.split("\t"), strict=True)) for line in lines]


def tint(goodness: float) -> str:
    """Background for goodness in [0, 1] (1 = best): red -> pale -> blue, lightened."""
    g = min(max(goodness, 0.0), 1.0)
    low, high, t = (WORSE, MIDDLE, g / 0.5) if g < 0.5 else (MIDDLE, BETTER, (g - 0.5) / 0.5)
    rgb = [
        round((a + (b - a) * t) * (1 - TINT) + 255 * TINT) for a, b in zip(low, high, strict=True)
    ]
    return f"background-color: rgb({rgb[0]}, {rgb[1]}, {rgb[2]})"


def goodness(column: str, value: float, values: list[float]) -> float | None:
    direction = DIRECTIONS.get(column)
    if direction is None:
        return None
    if direction == "weight":
        return value / max(values) if max(values) > 0 else 0.0
    if direction == "max":
        delta = 2 * (max(values) - value)  # 2*dlogL, on the information-criterion scale
    else:
        delta = value - min(values)
    return 1 - min(delta, DELTA_CAP) / DELTA_CAP


def number(column: str, text: str) -> str:
    value = float(text)
    if column.startswith("w_"):
        return f"{value:.3g}"
    if column == "free_parameters":
        return str(int(value))
    return f"{value:.2f}" if column == "delta_bic" else f"{value:.1f}"


def heat_cells(rows: list[dict[str, str]], columns: list[str]) -> list[dict[str, Cell]]:
    values = {c: [float(r[c]) for r in rows] for c in columns}
    out = []
    for row in rows:
        cells = {}
        for column in columns:
            g = goodness(column, float(row[column]), values[column])
            style = tint(g) if g is not None else ""
            cells[column] = Cell(number(column, row[column]), style)
        out.append(cells)
    return out


MODELFINDER_COLUMNS = ["log_likelihood", "aic", "w_aic", "aicc", "w_aicc", "bic", "w_bic"]
COMBINATION_COLUMNS = [
    "log_likelihood", "free_parameters", "aic", "aicc", "bic", "delta_bic", "w_bic",
]  # fmt: skip
HEADINGS = {
    "log_likelihood": "log L", "free_parameters": "params", "aic": "AIC", "w_aic": "AIC weight",
    "aicc": "AICc", "w_aicc": "AICc weight", "bic": "BIC", "w_bic": "BIC weight",
    "delta_bic": "ΔBIC",
}  # fmt: skip


def modelfinder_table(text: str, trees: list[str] | None = None) -> Table:
    """ModelFinder's ranking; `trees` names the tree sources it belongs to, when an
    alignment's sources have different rankings (jobs from before MODELFINDER ran once per
    alignment, when provided trees ran their own)."""
    rows = read_tsv(text)
    heated = heat_cells(rows, MODELFINDER_COLUMNS)
    body = []
    for row, cells in zip(rows, heated, strict=True):
        marker = " ✓" if row.get("bic_in_95") == "yes" else ""
        body.append(
            [Cell(row["model"] + marker, label=True), *(cells[c] for c in MODELFINDER_COLUMNS)]
        )
    subtitle = (
        "Candidates ranked by BIC, scored on ModelFinder's initial tree. "
        "✓ = in the 95% BIC confidence set."
    )
    title = "ModelFinder candidates"
    if trees:
        title += f" · tree{'s' if len(trees) > 1 else ''}: {', '.join(trees)}"
    return Table(
        title=title,
        subtitle=subtitle,
        columns=["Model", *(HEADINGS[c] for c in MODELFINDER_COLUMNS)],
        rows=body[:MODELFINDER_SHOWN],
        hidden_rows=body[MODELFINDER_SHOWN:],
    )


@dataclass(frozen=True)
class AlignmentStats:
    sequences: int
    columns: int
    gap_percent: float  # gap characters as a share of all cells (sequences x columns)


def alignment_stats(fasta: str) -> AlignmentStats:
    rows: list[str] = []
    for line in fasta.splitlines():
        if line.startswith(">"):
            rows.append("")
        elif rows:
            rows[-1] += line.strip()
    columns = max((len(row) for row in rows), default=0)
    cells = len(rows) * columns
    gaps = sum(row.count("-") for row in rows)
    return AlignmentStats(len(rows), columns, 100 * gaps / cells if cells else 0.0)


@dataclass
class AlignmentResult:
    """Everything shown for one alignment: its combinations, and ModelFinder's candidates."""

    aligner: str
    label: str
    alignment_href: str
    stats: AlignmentStats | None
    table: Table
    modelfinder: list[Table]
    codon_href: str | None = None  # the codon alignment, for coding-DNA jobs
    data: str = "amino"  # this table's models: "amino" or "codon" (ranked apart)
    first: bool = True  # the alignment's first table (the overview lists each once)
    best: str = ""  # the top-ranked combination's key
    comparison: str = ""  # codon tables: their best tree against the amino-acid one
    codeml: list[CodemlResult] = field(default_factory=list)  # codon tables

    @property
    def anchor(self) -> str:
        return f"alignment-{self.aligner}" + ("-codon" if self.data == "codon" else "")


def ancestor_files(key: str, available: Collection[str]) -> list[str]:
    """A combination's ancestor files, as its row's ZIP holds them: the reconstructions
    (original, gap versions, coding DNA), posteriors, gap probabilities and property
    measurements. EXPORT's byte index is for the viewers, not for readers."""
    folder = f"{key}/"
    out = []
    for path in available:
        if not path.startswith(folder):
            continue
        name = path[len(folder) :]
        if name.startswith("properties/") or name == "gaps-ml.tsv":
            out.append(path)
        elif "/" not in name and name.startswith("ancestors") and not name.endswith(".index.tsv"):
            out.append(path)
    return sorted(out)


def ancestor_contents(key: str, available: Collection[str]) -> list[str]:
    """What a row's ZIP holds, in words, for the download's tooltip."""
    out = ["original"]
    out += [name for method, name in GAP_VERSIONS
            if f"{key}/ancestors.gaps-{method}.fasta" in available]  # fmt: skip
    if f"{key}/ancestors.cds.fasta" in available:
        out.append("coding DNA")
    out.append("site probabilities")
    if f"{key}/gaps-ml.tsv" in available:
        out.append("gap probabilities")
    # Property posteriors: written only when the job asked for any.
    if f"{key}/properties/index.json" in available:
        out.append("properties")
    return out


def data_type(combo: dict[str, Any]) -> str:
    """'codon' for a codon model's combination (summary.json's `data`; older jobs had
    amino-acid models only)."""
    if combo.get("data"):
        return str(combo["data"])
    model = combo["model"]
    codon = model == catalog.CODON_MODELFINDER or catalog.is_codon_model(model)
    return "codon" if codon else "amino"


# Gap-masked ancestor versions, linked when present (ancestors.gaps-<method>.fasta).
GAP_VERSIONS = [
    ("parsimony", "shared presence"),
    ("fitch", "Fitch"),
    ("ml", "ML"),
    ("prank", "PRANK"),
    ("grasp", "GRASP"),
]

RESULT_COLUMNS = [
    "log_likelihood", "free_parameters", "aic", "aicc", "bic", "delta_bic", "w_bic",
]  # fmt: skip
# The alignment tables sort by these in the browser; the server's order is by rank.
SORTABLE = frozenset(["Rank", "Model", *(HEADINGS[c] for c in RESULT_COLUMNS)])


def combination_order(primary: str) -> Any:
    """Sort key: the primary's alignment first, then by aligner; ModelFinder first."""

    def order(combo: dict[str, Any]) -> tuple[bool, str, bool, str, str]:
        return (
            combo["aligner"] != primary.split("/")[0],
            combo["aligner"],
            not catalog.is_modelfinder(combo["model"]),
            combo["tree"],
            combo["model"],
        )

    return order


def fitted_model(combo: dict[str, Any], stats: dict[str, dict[str, str]]) -> str:
    row = stats.get(combo["key"])
    return str(row["model"] if row else combo.get("fit", {}).get("model", ""))


def merge_modelfinder_duplicates(
    group: list[dict[str, Any]], stats: dict[str, dict[str, str]]
) -> tuple[list[dict[str, Any]], dict[str, list[str]]]:
    """Drop fixed-model rows that repeat ModelFinder's pick on the same tree source.

    They are the same fit computed twice, so one row shows it: ModelFinder's. Returns
    the remaining combinations and {ModelFinder key: keys of the rows merged into it}.
    """
    merged: dict[str, list[str]] = {}
    dropped: set[str] = set()
    for finder in group:
        picked = fitted_model(finder, stats)
        if not catalog.is_modelfinder(finder["model"]) or not picked:
            continue
        for combo in group:
            if (
                not catalog.is_modelfinder(combo["model"])
                and combo["tree"] == finder["tree"]
                and fitted_model(combo, stats) == picked
            ):
                merged.setdefault(finder["key"], []).append(combo["key"])
                dropped.add(combo["key"])
    return [c for c in group if c["key"] not in dropped], merged


def alignment_results(
    summary: dict[str, Any],
    comparison: str | None,
    modelfinder: dict[str, str],
    href: Any,
    alignments: dict[str, str] | None = None,
    available: set[str] | None = None,
    rooting: dict[str, dict[str, Any]] | None = None,
    viewer: Any = None,
    tree_viewer: Any = None,
    reconciliation: dict[str, dict[str, Any]] | None = None,
    topology: dict[str, dict[str, str]] | None = None,
    download: Any = None,
) -> list[AlignmentResult]:
    """One AlignmentResult per alignment, in display order.

    `comparison` is model_comparison.tsv (None for jobs from before it existed: rows
    then show links without statistics); `modelfinder` maps combination keys to their
    modelfinder.tsv; `href(relative_path)` gives a result file's URL; `alignments` maps
    aligner keys to their alignment FASTA text, for the heading's statistics;
    `available` (result file paths) decides which gap-masked ancestor versions to link;
    `rooting` maps combination keys to their rooting.json (outgroup-rooted jobs);
    `viewer(key)` and `tree_viewer(key)`, if given, link a combination's ancestors and
    tree in their viewers; `reconciliation` maps combination keys to their
    reconciliation.json (jobs reconciled with a species tree); `topology` maps them to
    their row of the AU test (`topology_rows`; jobs with GeneRax's reconciled tree).
    """
    rooting = rooting or {}
    reconciliation = reconciliation or {}
    topology = topology or {}
    alignments = alignments or {}
    available = available or set()
    with_ancestors = summary.get("stop_after", "ancestors") == "ancestors"
    primary = summary.get("primary", "")
    combos = sorted(summary.get("combinations", []), key=combination_order(primary))
    stats = {
        f"{r['aligner']}/{r['tree']}/{r['model_key']}": r
        for r in (read_tsv(comparison) if comparison else [])
    }
    several_trees = len({c["tree"] for c in combos}) > 1

    def codon_href(aligner: str) -> str | None:
        path = f"{aligner}/alignment.cds.fasta"
        return href(path) if path in available else None

    with_codons = any(data_type(c) == "codon" for c in combos)
    tables = [
        (aligner, data)
        for aligner in dict.fromkeys(c["aligner"] for c in combos)
        for data in ("amino", "codon")
        if any(c["aligner"] == aligner and data_type(c) == data for c in combos)
    ]
    out: list[AlignmentResult] = []
    for aligner, data in tables:
        group = [c for c in combos if c["aligner"] == aligner and data_type(c) == data]
        # Rank order when statistics exist; otherwise display order.
        if all(c["key"] in stats for c in group):
            group.sort(
                key=lambda c: (
                    float(stats[c["key"]]["bic"]),
                    not catalog.is_modelfinder(c["model"]),
                )
            )
        group, merged = merge_modelfinder_duplicates(group, stats)
        rows = [stats[c["key"]] for c in group if c["key"] in stats]
        # A single row has nothing to be compared against: no tints.
        tints: dict[str, dict[str, Cell]] = {}
        if len(rows) > 1:
            for stat, heated in zip(rows, heat_cells(rows, RESULT_COLUMNS), strict=True):
                tints[f"{stat['aligner']}/{stat['tree']}/{stat['model_key']}"] = heated
        reconciled = any(c["key"] in reconciliation for c in group)
        tested = any(c["key"] in topology for c in group)
        body = []
        exported = []
        for combo in group:
            key = combo["key"]
            row = stats.get(key)
            fitted = row["model"] if row else combo.get("fit", {}).get("model", "")
            model = combo.get("model_label", combo["model"])
            if catalog.is_modelfinder(combo["model"]) and fitted:
                model = f"ModelFinder\u00a0→ {fitted}"  # wraps after the arrow only
                if key in merged:
                    model += " (also ticked)"
            cells = [Cell(row["rank_bic"] if row else "–")]
            if several_trees:
                cells.append(Cell(combo.get("tree_label", combo["tree"]), label=True))
            cells.append(Cell(model, label=True))
            for column in RESULT_COLUMNS:
                if row is None:
                    cells.append(Cell("–"))
                elif key in tints:
                    cells.append(tints[key][column])
                else:
                    cells.append(Cell(number(column, row[column])))
            if reconciled:
                events = reconciliation.get(key)
                if events is None:
                    cells.append(Cell("–"))
                else:
                    text = f"{events.get('duplications', 0)} / {events.get('losses', 0)}"
                    if weak := events.get("weak"):
                        text += f" · {weak} doubtful"
                    cells.append(Cell(text, href=href(f"{key}/reconciliation.tsv")))
            if tested:
                cells.append(au_cell(topology.get(key), href))
            exported.append(export_row(combo, row, fitted, key in merged, several_trees,
                                       reconciliation.get(key, {}) if reconciled else None,
                                       topology.get(key, {}) if tested else None))  # fmt: skip
            tree_links = [("view", tree_viewer(key))] if tree_viewer else []
            tree_links.append(("tree.nwk", href(f"{key}/tree.nwk")))
            if f"{key}/tree.rooted.nwk" in available:
                tree_links.append(("rooted", href(f"{key}/tree.rooted.nwk")))
            report = rooting.get(key)
            warning = ""
            if report and not report.get("monophyletic", True):
                inside = ", ".join(report.get("intruders", []))
                warning = f"⚠ outgroup is not one group here (also inside: {inside})"
            cells.append(Cell("", label=True, links=tree_links, badge=warning, icons=True))
            if not with_ancestors:
                body.append(cells)
                continue
            # View, and every ancestor file as one ZIP (`download(key)`), as two icons: a
            # link per file made the column the table's widest. The cell's tooltip says
            # what the ZIP holds.
            links = [("view", viewer(key))] if viewer else []
            links += [("download", download(key))] if download else []
            holds = "ZIP of " + ", ".join(ancestor_contents(key, available))
            cells.append(Cell("", label=True, links=links, title=holds, icons=True))
            body.append(cells)
        columns = ["Rank", *(["Tree source"] if several_trees else []), "Model"]
        columns += [HEADINGS[c] for c in RESULT_COLUMNS]
        columns += ["Duplications / losses"] if reconciled else []
        columns += ["AU\u00a0p"] if tested else []
        columns += ["Tree"]
        columns += ["Ancestors"] if with_ancestors else []
        first = group[0]
        label = first.get("aligner_label", aligner)
        kind = {"amino": " · amino-acid models", "codon": " · codon models"}[data]
        title = f"Alignment: {label}" + (kind if with_codons else "")
        subtitle = (
            "One row per combination. Scores are "
            "comparable only within this alignment"
            + (" and data type (codon models fit the codon alignment, amino-acid "
               "models its translation)." if with_codons else ".")
        )  # fmt: skip
        if any(c["tree"] == "generax" for c in group):
            subtitle += (
                " GeneRax's reconciled tree trades fit to the sequences for fit to the "
                "species tree, so it ranks low by BIC; the AU test says whether the "
                "sequences alone reject a tree (p < 0.05)."
            )
        if any((reconciliation.get(c["key"]) or {}).get("weak") for c in group):
            subtitle += (
                " A doubtful duplication's two copies share under 10% of their species; a real "
                "duplication leaves both copies in many of the same species, so this usually "
                "means the tree grouped the sequences wrongly."
            )
        # One ModelFinder run per alignment, shared by every tree source: one table.
        rankings: dict[str, list[str]] = {}
        for c in group:
            if c["key"] in modelfinder:
                trees = rankings.setdefault(modelfinder[c["key"]], [])
                trees.append(c.get("tree_label", c["tree"]))
        finder = [
            modelfinder_table(text, trees if len(rankings) > 1 else None)
            for text, trees in rankings.items()
        ]
        out.append(AlignmentResult(
            aligner=aligner,
            label=label,
            alignment_href=href(f"{first['key']}/alignment.fasta"),
            codon_href=codon_href(aligner),
            data=data,
            first=all(r.aligner != aligner for r in out),
            best=first["key"],
            stats=alignment_stats(alignments[aligner]) if aligner in alignments else None,
            table=Table(title=title, subtitle=subtitle, columns=columns, rows=body,
                        sortable=SORTABLE,
                        export_columns=export_columns(several_trees, reconciled, tested),
                        export_rows=exported),
            modelfinder=finder,
        ))  # fmt: skip
    # Alignments without combinations (a job that stopped after alignments).
    shown = {result.aligner for result in out}
    for entry in summary.get("alignments", []):
        if entry["key"] in shown:
            continue
        out.append(AlignmentResult(
            aligner=entry["key"],
            label=entry.get("label", entry["key"]),
            alignment_href=href(entry["path"]),
            codon_href=codon_href(entry["key"]),
            stats=alignment_stats(alignments[entry["key"]]) if entry["key"] in alignments else None,
            table=Table(title=f"Alignment: {entry.get('label', entry['key'])}", subtitle="",
                        columns=[], rows=[]),
            modelfinder=[],
        ))  # fmt: skip
    return out


def topology_rows(aligner: str, model: str, tsv: str) -> dict[str, dict[str, str]]:
    """An AU test's rows (<aligner>/topology_tests/<model>.tsv) by combination key; the
    row keeps the file's path for its link."""
    path = f"{aligner}/topology_tests/{model}.tsv"
    return {
        f"{aligner}/{row['tree']}/{model}": row | {"path": path}
        for row in read_tsv(tsv)
        if row.get("tree")
    }


EXPORT_STATISTICS = [
    "log_likelihood", "free_parameters", "aic", "w_aic", "aicc", "w_aicc", "bic", "delta_bic",
    "w_bic",
]  # fmt: skip
EXPORT_HEADINGS = {**HEADINGS, "log_likelihood": "log L"}


def export_columns(several_trees: bool, reconciled: bool, tested: bool) -> list[str]:
    columns = ["Rank", "Combination", *(["Tree source"] if several_trees else [])]
    columns += ["Model", "Fitted model", *(EXPORT_HEADINGS[c] for c in EXPORT_STATISTICS)]
    columns += ["Duplications", "Losses", "Doubtful duplications"] if reconciled else []
    columns += ["AU p"] if tested else []
    return columns


def export_row(
    combo: dict[str, Any],
    row: dict[str, str] | None,
    fitted: str,
    merged: bool,
    several_trees: bool,
    events: dict[str, Any] | None,
    au: dict[str, str] | None,
) -> list[str]:
    """One combination as the downloads give it (export_columns): unrounded. `events`
    and `au` are None when the table has no such columns, {} when this row lacks them."""
    model = combo.get("model_label", combo["model"])
    if catalog.is_modelfinder(combo["model"]):
        model = "ModelFinder" + (" (also ticked)" if merged else "")
    out = [row["rank_bic"] if row else "", combo["key"]]
    out += [combo.get("tree_label", combo["tree"])] if several_trees else []
    out += [model, fitted]
    out += [row.get(c, "") if row else "" for c in EXPORT_STATISTICS]
    if events is not None:
        out += [str(events.get(k, "")) for k in ("duplications", "losses", "weak")]
    if au is not None:
        out.append(au.get("p_au", ""))
    return out


def au_cell(row: dict[str, str] | None, href: Any) -> Cell:
    """The AU p of one combination's tree: "0.33", or "0.012 · rejected" below 0.05."""
    try:
        p = float(row["p_au"]) if row else None
    except KeyError, ValueError:
        p = None
    if p is None or row is None:
        return Cell("–")
    text = f"{p:.2g}" + (" · rejected" if p < 0.05 else "")
    return Cell(text, href=href(row["path"]))


@dataclass
class CodemlResult:
    """PAML codeml on one codon combination: a row per site model, the tests for
    positive selection, and the sites BEB places under it."""

    combination: str  # key
    label: str  # tree source x model, as in the codon table
    table: Table
    tests: list[str]  # one sentence each
    sites: list[str]  # one sentence per model with BEB sites


def codeml_result(
    combo: dict[str, Any],
    codeml_tsv: str,
    selection_tsv: str,
    beb_tsv: str,
    href: Any,
    viewer: Any = None,
    available: set[str] | None = None,
) -> CodemlResult:
    """The codeml section for one combination, from its codeml/ files."""
    available = available or set()
    key = combo["key"]
    rows = []
    for row in read_tsv(codeml_tsv):
        model = row["model"]
        lnl, k = float(row["log_likelihood"]), int(row["free_parameters"])
        folder = f"{key}/codeml/{model}"
        links = [("view", viewer(key, model))] if viewer else []
        links += [("protein", href(f"{folder}/ancestors.fasta"))]
        if f"{folder}/ancestors.cds.fasta" in available or not available:
            links.append(("coding DNA", href(f"{folder}/ancestors.cds.fasta")))
        links.append(("site probabilities", href(f"{folder}/ancestors_pp.tsv")))
        if f"{folder}/branches.tsv" in available:  # free ratios: dN/dS per branch
            links.append(("per branch", href(f"{folder}/branches.tsv")))
        rows.append([
            Cell(model, label=True), Cell(row["estimates"], label=True),
            Cell(f"{lnl:.1f}"), Cell(str(k)), Cell(f"{2 * k - 2 * lnl:.1f}"),
            Cell("", label=True, links=links),
        ])  # fmt: skip
    tests = []
    for row in read_tsv(selection_tsv):
        p = float(row["p_value"])
        verdict = "positive selection supported" if p < 0.05 else "not significant"
        tests.append(
            f"{row['test']}: 2ΔlnL = {float(row['statistic']):.1f}, p = {p:.2g} "
            f"(χ², 2 df): {verdict} at 5%."
        )
    sites: dict[str, list[str]] = {}
    for row in read_tsv(beb_tsv):
        if float(row["p_positive"]) > 0.95:
            sites.setdefault(row["model"], []).append(f"{row['site']}{row['residue']}")
    label = " · ".join(
        part for part in (combo.get("tree_label", combo["tree"]), combo.get("model_label"))
        if part
    )  # fmt: skip
    return CodemlResult(
        combination=key,
        label=label,
        table=Table(
            title=f"PAML codeml on {label}",
            subtitle=(
                "Site models (and the free-ratio branch model, if chosen) on this "
                "combination's tree; log-likelihoods and AIC compare the models with each "
                "other (same alignment and tree)."
            ),
            columns=["Model", "dN/dS estimates", "log L", "params", "AIC", "Ancestors"],
            rows=rows,
        ),
        tests=tests,
        sites=[
            f"{model}: {len(found)} site{'s' * (len(found) > 1)} with P(dN/dS > 1) > 0.95 "
            f"(Bayes empirical Bayes): {', '.join(found)} (numbered by alignment column; "
            "residues of the first sequence)."
            for model, found in sites.items()
        ],
    )
