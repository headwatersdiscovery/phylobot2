"""The pipeline's stages as a small graph, for the job page (and the ETA).

A stage is one Nextflow process (ALIGN_MAFFT, TREE_IQTREE, ...); it runs once per
alignment, combination, etc. `expected_kinds` says which stages a job runs and how many
tasks each; `graph` lays them out left to right with their status and times, drawn by
`_job_stages.html` as SVG.

DEPENDS lists the data each process takes from others, read from pipeline/main.nf; keep
the two in step. The drawing shows only the stages a job runs, and only the edges not
implied by a longer path (a transitive reduction), so the picture stays sparse.
"""

from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any

from django.utils import timezone

from phylobot import workload
from phylobot.durations import seconds
from phylobot.executors import (
    alignment_count,
    as_list,
    branch_site,
    branch_site_runs,
    codeml_combinations,
    codeml_runs,
    coding,
    codon_amino_trees,
    codon_model_count,
    codon_modelfinder,
    codon_search,
    combination_count,
    esr,
    esr_tips,
    model_count,
    modelfinder_runs,
    planned_combinations,
    properties_only,
    property_targets,
    reconciles,
    simulation,
    simulation_seeds,
    task_name,
    uncertainty_only,
)

ALIGNERS = {
    "MAFFT": "mafft",
    "MUSCLE": "muscle",
    "Clustal Ω": "clustalo",
    "PRANK": "prank",
    "PRANK codon": "prank_codon",
}

LABELS = {
    "VALIDATE": "Check input",
    **{f"ALIGN_{key.upper()}": f"Align: {name}" for name, key in ALIGNERS.items()},
    "CODON_ALIGNMENT": "Codon alignment",
    "EXPORT_ALIGNMENT": "Save alignment",
    "MODELFINDER": "ModelFinder",
    "MODELFINDER_CODON": "ModelFinder (codon)",
    "TREE_IQTREE": "IQ-TREE search",
    "TREE_RAXML": "RAxML-NG search",
    "FIT_IQTREE": "IQ-TREE fit",  # provided and RAxML-NG topologies
    "FIT_CODON_AA": "Codon fit on aa tree",  # codon models on amino-acid trees
    "TREE_GENERAX": "GeneRax search",
    "FIT_GENERAX": "IQ-TREE fit (GeneRax)",
    "TOPOLOGY_TEST": "Topology test (AU)",
    "ANCESTRAL_IQTREE": "Ancestors",
    "ROOT_TREE": "Root tree",
    "GAPS_ML": "Gaps: ML",
    "PROPERTIES": "Properties",
    "UNCERTAINTY": "Uncertainty",
    "GAPS_PRANK": "Gaps: PRANK",
    "GAPS_GRASP": "Gaps: GRASP",
    "EXPORT": "Export results",
    "CODEML": "codeml",
    "CODEML_EXPORT": "codeml results",
    "BRANCH_SITE_EXPORT": "Branch-site test",
    "EXPORT_TREE": "Export trees",
    "SPECIES_TREE": "Species tree",
    "RECONCILE": "Reconcile (GeneRax)",
    "RECONCILE_EXPORT": "Reconciliation",
    "SUMMARY": "Summary",
    "ANCESTOR_MAP": "Ancestor map",
    "SIM_PREPARE": "Prepare",
    "SIM_ALISIM": "AliSim",
    "SIM_COLLECT": "Collect truth",
    "ESR_FIT": "Fit the model",
    "ESR_TIP": "Hold out a sequence",
    "ESR_COLLECT": "Score",
}
ALIGN = [f"ALIGN_{key.upper()}" for key in ALIGNERS.values()]
TREES = ["TREE_IQTREE", "TREE_RAXML", "FIT_IQTREE", "FIT_CODON_AA", "FIT_GENERAX"]
# The trees that guide SPECIES_TREE: every source but GeneRax, which needs its output.
SEARCHED = TREES[:-1]
# Where alignments come from: the aligners, and for coding DNA CODON_ALIGNMENT (both
# forms of each aligner's alignment; the protein form feeds every later step).
ALIGNED = [*ALIGN, "CODON_ALIGNMENT"]

# stage -> the stages whose outputs it reads (pipeline/main.nf)
DEPENDS: dict[str, list[str]] = {
    "VALIDATE": [],
    **{stage: ["VALIDATE"] for stage in ALIGN},
    "CODON_ALIGNMENT": ["VALIDATE", *ALIGN],
    "EXPORT_ALIGNMENT": ["VALIDATE", *ALIGNED],
    "MODELFINDER": ["VALIDATE", *ALIGNED],
    "MODELFINDER_CODON": ["VALIDATE", *ALIGNED],
    "TREE_IQTREE": ["VALIDATE", *ALIGNED, "MODELFINDER", "MODELFINDER_CODON"],
    "TREE_RAXML": ["VALIDATE", *ALIGNED, "MODELFINDER"],
    "FIT_IQTREE": ["VALIDATE", *ALIGNED, "MODELFINDER", "MODELFINDER_CODON", "TREE_RAXML"],
    "FIT_CODON_AA": [*ALIGNED, "MODELFINDER_CODON", "TREE_IQTREE", "FIT_IQTREE"],
    "SPECIES_TREE": ["VALIDATE", *SEARCHED],
    "TREE_GENERAX": [*ALIGNED, "TREE_IQTREE", "FIT_IQTREE", "SPECIES_TREE"],
    "FIT_GENERAX": [*ALIGNED, "MODELFINDER", "TREE_GENERAX"],
    "TOPOLOGY_TEST": [*ALIGNED, *TREES],
    "ANCESTRAL_IQTREE": [*ALIGNED, *TREES],
    "ROOT_TREE": ["VALIDATE", "ANCESTRAL_IQTREE"],
    "GAPS_ML": ["ANCESTRAL_IQTREE"],
    "GAPS_PRANK": ["ROOT_TREE"],
    "GAPS_GRASP": ["ROOT_TREE"],
    "EXPORT": [
        "VALIDATE",
        "ANCESTRAL_IQTREE",
        "MODELFINDER",
        "MODELFINDER_CODON",
        *TREES,
        "ROOT_TREE",
        "GAPS_ML",
        "GAPS_PRANK",
        "GAPS_GRASP",
    ],  # fmt: skip
    "PROPERTIES": ["EXPORT"],
    "UNCERTAINTY": ["EXPORT"],  # a follow-up: a finished job's published ancestors
    "CODEML": ["VALIDATE", *ALIGNED, "ANCESTRAL_IQTREE", "EXPORT"],
    "CODEML_EXPORT": ["CODEML"],
    "BRANCH_SITE_EXPORT": ["VALIDATE", "CODEML"],
    "EXPORT_TREE": ["VALIDATE", "MODELFINDER", "MODELFINDER_CODON", *TREES],
    "RECONCILE": [*ALIGNED, *TREES, "SPECIES_TREE"],
    "RECONCILE_EXPORT": ["VALIDATE", "ANCESTRAL_IQTREE", "SPECIES_TREE", "RECONCILE"],
    "SUMMARY": ["VALIDATE", *ALIGNED, "EXPORT", "EXPORT_TREE"],
    "ANCESTOR_MAP": ["EXPORT", "EXPORT_TREE"],
    "SIM_PREPARE": [],
    "SIM_ALISIM": ["SIM_PREPARE"],
    "SIM_COLLECT": ["SIM_ALISIM"],
    "ESR_FIT": [],
    "ESR_TIP": ["ESR_FIT"],
    "ESR_COLLECT": ["ESR_FIT", "ESR_TIP"],
}
ORDER = list(DEPENDS)  # a topological order; breaks ties within a column
POSITION = {name: i for i, name in enumerate(ORDER)}

# Drawing, in SVG user units (the SVG scales down to the page width).
NODE_W, NODE_H = 118, 50
COL_GAP, ROW_GAP = 22, 10
PAD = 2


def kind(name: str) -> str:
    """ALIGN_MAFFT from "ALIGN_MAFFT (mafft)" (or "PHYLOBOT:ALIGN_MAFFT (mafft)", as
    recorded by workers before they dropped the workflow prefix)."""
    return task_name(name).split(" ")[0]


def expected_kinds(params: dict[str, Any]) -> Counter[str]:
    """How many tasks of each kind a job with these params runs (their sum is
    executors.expected_tasks)."""
    if branch_site(params):
        return Counter({"VALIDATE": 1, "CODEML": len(branch_site_runs()), "BRANCH_SITE_EXPORT": 1})
    if properties_only(params):
        return Counter({"PROPERTIES": property_targets(params)})
    if uncertainty_only(params):
        return Counter({"UNCERTAINTY": 1})
    if simulation(params):
        n = len(simulation_seeds(params))
        return +Counter({"SIM_PREPARE": n, "SIM_ALISIM": n, "SIM_COLLECT": n})
    if esr(params):
        return Counter({"ESR_FIT": 1, "ESR_TIP": esr_tips(params), "ESR_COLLECT": 1})
    aligners = [a for a in as_list(params.get("aligners", "mafft")) if a != "none"]
    counts: Counter[str] = Counter({"VALIDATE": 1, "SUMMARY": 1})
    for aligner in aligners:
        counts[f"ALIGN_{aligner.upper()}"] += 1
    if coding(params):
        counts["CODON_ALIGNMENT"] = len(aligners)
    counts["EXPORT_ALIGNMENT"] = alignment_count(params)
    combos = combination_count(params)
    planned = planned_combinations(params)
    if combos and planned is not None:  # an extension: only what it runs (main.nf skips)
        counts["MODELFINDER"] = modelfinder_runs(params, planned)
        counts["TREE_IQTREE"] = sum(t == "iqtree" for _a, t, _m in planned)
        counts["TREE_RAXML"] = sum(t == "raxml" for _a, t, _m in planned)
        counts["FIT_IQTREE"] = sum(t != "iqtree" for _a, t, _m in planned)
    elif combos:
        builders = [t for t in as_list(params.get("tree_builders", "iqtree")) if t != "none"]
        alignments, models = alignment_count(params), model_count(params)
        if params.get("modelfinder", True):
            counts["MODELFINDER"] = alignments
        if codon_modelfinder(params):
            counts["MODELFINDER_CODON"] = alignments
        codon = codon_model_count(params)  # no RAxML-NG for codon models
        searches = models * ("iqtree" in builders) + codon * codon_search(params)
        counts["TREE_IQTREE"] = alignments * searches
        counts["TREE_RAXML"] = alignments * ("raxml" in builders) * models
        # IQ-TREE fits each provided tree and each RAxML-NG tree
        trees = params.get("tree") or params.get("tree_labels") or []  # provided trees
        raxml = "raxml" in builders
        counts["FIT_IQTREE"] = alignments * (len(trees) * (models + codon) + raxml * models)
        counts["FIT_CODON_AA"] = alignments * codon * codon_amino_trees(params)
        if "generax" in builders:  # amino-acid models only
            for step in ("TREE_GENERAX", "FIT_GENERAX", "TOPOLOGY_TEST"):
                counts[step] = alignments * models
    if combos:
        counts["ANCESTOR_MAP"] = 1  # job-level, from every exported tree
    if combos and params.get("stop_after", "ancestors") == "ancestors":
        counts["ANCESTRAL_IQTREE"] = combos
        counts["EXPORT"] = combos
        if as_list(params.get("outgroup")):
            counts["ROOT_TREE"] = combos
        for method in as_list(params.get("gap_methods")):
            if method in ("ml", "prank", "grasp"):
                counts[f"GAPS_{method.upper()}"] = combos
        if as_list(params.get("properties")):
            counts["PROPERTIES"] = combos
    elif combos:
        counts["EXPORT_TREE"] = combos
    if chosen := codeml_combinations(params):
        counts["CODEML"] = chosen * len(codeml_runs(params))
        counts["CODEML_EXPORT"] = chosen
    if reconciles(params):
        counts["SPECIES_TREE"] = 1
        counts["RECONCILE"] = counts["RECONCILE_EXPORT"] = combos
    return +counts


def edges(stages: list[str]) -> list[tuple[str, str]]:
    """Dependencies among `stages`, without those implied by a longer path."""
    present = set(stages)
    parents = {s: [p for p in DEPENDS.get(s, []) if p in present] for s in stages}

    def ancestors(stage: str) -> set[str]:
        found: set[str] = set()
        todo = list(parents[stage])
        while todo:
            if (p := todo.pop()) not in found:
                found.add(p)
                todo.extend(parents[p])
        return found

    reduced = []
    for stage in stages:
        for parent in parents[stage]:
            others = set(parents[stage]) - {parent}
            if not any(parent in ancestors(other) for other in others):
                reduced.append((parent, stage))
    return reduced


@dataclass
class Stage:
    name: str  # the Nextflow process
    label: str
    expected: int
    tasks: list[dict[str, Any]] = field(default_factory=list)
    status: str = "waiting"  # waiting, running, succeeded, failed, stopped, not run
    done: int = 0
    elapsed: float | None = None  # seconds
    started: datetime | None = None
    completed: datetime | None = None
    x: float = 0
    y: float = 0

    @property
    def group(self) -> str:
        """The step's pipeline stage (align, trees, ...), for skins that colour by it: the
        work estimate's group, except that rooting is tree work (the estimate counts it
        with the ancestors that follow it)."""
        return "trees" if self.name == "ROOT_TREE" else workload.group(self.name)

    @property
    def count(self) -> str:
        """'2/3' for a stage of several tasks, while not all have finished."""
        total = max(self.expected, len(self.tasks))
        return f"{self.done}/{total}" if total > 1 and self.done < total else ""

    @property
    def tooltip(self) -> str:
        lines = [f"{self.label} ({self.name})"]
        for task in self.tasks:
            lines.append(f"{task['name']}: {task['status'].lower()}"
                         + (f", {task['realtime']}" if task.get("realtime") else ""))  # fmt: skip
        return "\n".join(lines)


def when(text: str | None) -> datetime | None:
    try:
        return datetime.fromisoformat(text) if text else None
    except ValueError:
        return None


def summarise(stage: Stage, finished: bool, now: datetime) -> None:
    """Status and times of a stage from its tasks."""
    statuses = [task["status"] for task in stage.tasks]
    stage.done = sum(s in ("COMPLETED", "CACHED") for s in statuses)
    running = "RUNNING" in statuses
    starts = [t for task in stage.tasks if (t := when(task.get("start")))]
    completes = [t for task in stage.tasks if (t := when(task.get("complete")))]
    stage.started = min(starts) if starts else None
    if any(s in ("FAILED", "ABORTED") for s in statuses):
        stage.status = "failed"
    elif running and not finished:
        stage.status = "running"
    elif stage.tasks and stage.done >= stage.expected:
        stage.status = "succeeded"
    elif stage.tasks:
        stage.status = "stopped" if finished else "running"
    else:
        stage.status = "not run" if finished else "waiting"
    if stage.status == "running" and stage.started:
        stage.elapsed = (now - stage.started).total_seconds()
    elif stage.status != "running":
        if stage.started and completes and not running:
            stage.completed = max(completes)
            stage.elapsed = (stage.completed - stage.started).total_seconds()
        elif len(stage.tasks) == 1:  # jobs from before the trace had timestamps
            stage.elapsed = seconds(stage.tasks[0].get("realtime", ""))


def layout(stages: list[Stage], links: list[tuple[str, str]]) -> tuple[float, float]:
    """Place stages in columns, the main line along the top row. A stage goes in the
    column just before its first consumer (so a stage feeding a later one sits next to
    it rather than beside its own inputs), a stage nothing reads (e.g. Save alignment)
    just after its inputs. Within a column, stages whose inputs are nearest come first,
    so edges that skip columns run below the top row. Returns width and height."""
    parents = {s.name: [p for p, c in links if c == s.name] for s in stages}
    children = {s.name: [c for p, c in links if p == s.name] for s in stages}
    earliest: dict[str, int] = {}
    for stage in stages:  # stages come in ORDER, a topological order
        earliest[stage.name] = 1 + max((earliest[p] for p in parents[stage.name]), default=-1)
    rank: dict[str, int] = {}
    for stage in reversed(stages):
        kids = children[stage.name]
        rank[stage.name] = min(rank[c] for c in kids) - 1 if kids else earliest[stage.name]
    columns: dict[int, list[Stage]] = {}
    for stage in stages:
        columns.setdefault(rank[stage.name], []).append(stage)
    row: dict[str, int] = {}

    def key(stage: Stage) -> tuple[float, int, bool, int]:
        ups = parents[stage.name]
        centre = sum(row[p] for p in ups) / len(ups) if ups else 0
        reach = max((rank[stage.name] - rank[p] for p in ups), default=1)
        return centre, reach, not children[stage.name], POSITION.get(stage.name, len(ORDER))

    for col in sorted(columns):
        columns[col].sort(key=key)
        for i, stage in enumerate(columns[col]):
            row[stage.name] = i
            stage.x = PAD + col * (NODE_W + COL_GAP)
            stage.y = PAD + i * (NODE_H + ROW_GAP)
    rows = max(len(column) for column in columns.values())
    cols = max(columns) + 1
    width = cols * NODE_W + (cols - 1) * COL_GAP
    return width + 2 * PAD, rows * NODE_H + (rows - 1) * ROW_GAP + 2 * PAD


@dataclass(frozen=True)
class Link:
    path: str  # SVG path data
    state: str  # "done" (upstream finished), "flowing" (into a running stage) or ""


def graph(params: dict[str, Any], progress: dict[str, Any], finished: bool) -> dict[str, Any]:
    """The job's stages, laid out, for `_job_stages.html`."""
    expected = expected_kinds(params)
    tasks = (progress or {}).get("tasks", [])
    for task in tasks:  # a stage we did not predict still gets drawn
        expected.setdefault(kind(task["name"]), 1)
    names = [s for s in ORDER if s in expected] + sorted(set(expected) - set(ORDER))
    # (unknown stages have no DEPENDS entry: drawn unconnected, after the rest)
    stages = [Stage(name, LABELS.get(name, name.title()), expected[name]) for name in names]
    by_name = {stage.name: stage for stage in stages}
    for task in tasks:
        by_name[kind(task["name"])].tasks.append(task)
    now = timezone.now()
    for stage in stages:
        summarise(stage, finished, now)
    links = edges(names)
    width, height = layout(stages, links)
    drawn = []
    for parent, child in links:
        a, b = by_name[parent], by_name[child]
        x1, y1 = a.x + NODE_W, a.y + NODE_H / 2
        x2, y2 = b.x, b.y + NODE_H / 2
        bend = (x2 - x1) / 2
        drawn.append(Link(
            f"M{x1:g},{y1:g} C{x1 + bend:g},{y1:g} {x2 - bend:g},{y2:g} {x2:g},{y2:g}",
            "flowing" if b.status == "running" else "done" if a.status == "succeeded" else "",
        ))  # fmt: skip
    return {"stages": stages, "links": drawn, "width": width, "height": height,
            "node_w": NODE_W, "node_h": NODE_H}  # fmt: skip
