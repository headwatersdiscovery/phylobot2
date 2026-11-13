"""Running a claimed job's pipeline. The worker calls `LocalDockerExecutor.run`.

An executor stages the job's inputs, runs Nextflow, reports progress through
`phylobot.jobs.heartbeat`, stores results, and returns the outcome; the worker then
marks the job. A later AWS Batch executor implements the same `run(job) -> Outcome`.

Job inputs live in object storage under jobs/<id>/inputs/ with fixed names
(sequences.fasta, alignment1.fasta..., tree1.nwk...); results and logs go to
jobs/<id>/results/. Only options on PARAM_ALLOWLIST pass from Job.params to Nextflow,
so nothing a user stores can choose a path or an unexpected parameter.
"""

import csv
import json
import logging
import os
import re
import shutil
import subprocess
import threading
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

from django.conf import settings
from django.db import transaction
from django.db.models import Q

from phylobot import (
    alignments,
    catalog,
    jobs,
    machine,
    manifest,
    site_settings,
    storage,
    work_archive,
)
from phylobot.durations import seconds
from phylobot.models import Job, JobAttempt, limits_for

logger = logging.getLogger(__name__)

PARAM_ALLOWLIST = {
    "aligners",
    "tree_builders",
    "modelfinder",
    "fixed_models",
    "models",
    "candidates",
    "bootstrap",
    "alrt",
    "seed",
    "alignment_labels",
    "alignment_keys",
    "skip_combinations",
    "modelfinder_reused",
    "tree_labels",
    "tree_branch_lengths",  # per provided tree: fit or keep (checked in main.nf)
    "gap_methods",
    "properties",
    "property_draws",
    "outgroup",
    "stop_after",
    "seqtype",
    "genetic_code",
    "codon_models",
    "codon_modelfinder",
    "codon_candidates",
    "codeml_models",
    "codeml_trees",
    "codon_trees",
    "codon_amino_trees",
    "branch_site",
    "properties_only",
    "uncertainty_only",
    "esr",
    "polytomies",
}
#: An extension's reused ModelFinder picks (phylobot.extend): <aligner>.model.txt and
#: <aligner>.modelfinder.tsv from the job it extends.
PICKS_DIR = "picks"
#: A properties-only run's staged inputs: the parent job's published
#: combinations, and the TSV naming exactly which ancestor of each to measure.
PROPERTY_DIR = "combinations"
PROPERTY_TARGETS = "property_targets.tsv"
#: An uncertainty run's staged inputs: the parent job's published
#: files, and the run's settings.
UNCERTAINTY_DIR = "published"
UNCERTAINTY_SETTINGS = "uncertainty_settings.json"
#: An ESR job's staged inputs: alignment.fasta, tree.nwk, model.txt and
#: settings.json (levels, the sequences to hold out, seed).
ESR_DIR = "esr"
SIMULATION_INPUTS = {
    "spec.json": "simulate",
    "sim_tree.nwk": "sim_tree",
    "sim_root.fasta": "sim_root",
}
# Reconciliation (phylobot.species): the species table and species tree, staged as files.
SPECIES_INPUTS = {
    "species.tsv": "species_map",
    "species_tree.nwk": "species_tree",
}
GENERIC_FAILURE = "The pipeline failed. The log below may help; if not, please report it."
INTERRUPTED = "The job was interrupted. You can resume it."  # the worker was stopped
DONE = {"COMPLETED", "CACHED"}


@dataclass(frozen=True)
class Outcome:
    status: str  # a terminal Job.Status
    error: str = ""


def as_list(value: Any) -> list[str]:
    """A parameter as a list: accepts a list or a comma-separated string."""
    if value in (None, "", True, False):
        return []
    items = value if isinstance(value, list) else str(value).split(",")
    return [str(item).strip() for item in items if str(item).strip()]


def staged_inputs(inputs: Path) -> dict[str, Any]:
    """Pipeline file parameters from the fixed input names present in `inputs`."""
    params: dict[str, Any] = {}
    if (inputs / "sequences.fasta").is_file():
        params["input"] = str(inputs / "sequences.fasta")
    alignments = sorted(inputs.glob("alignment*.fasta"), key=lambda p: int(p.stem[9:] or 0))
    trees = sorted(inputs.glob("tree*.nwk"), key=lambda p: int(p.stem[4:] or 0))
    if alignments:
        params["alignment"] = [str(p) for p in alignments]
    if trees:
        params["tree"] = [str(p) for p in trees]
    # A simulation (phylobot.simulations): its spec, and the tree and root it names.
    for filename, param in (SIMULATION_INPUTS | SPECIES_INPUTS).items():
        if (inputs / filename).is_file():
            params[param] = str(inputs / filename)
    if (inputs / PICKS_DIR).is_dir():
        params["modelfinder_picks"] = str(inputs / PICKS_DIR)
    # A properties-only run: a directory of staged combinations, and its plan.
    if (inputs / PROPERTY_DIR).is_dir():
        params["properties_dir"] = str(inputs / PROPERTY_DIR)
    if (inputs / PROPERTY_TARGETS).is_file():
        params["property_targets"] = str(inputs / PROPERTY_TARGETS)
    # An uncertainty run: the parent's published files, and the run's settings.
    if (inputs / UNCERTAINTY_DIR).is_dir():
        params["uncertainty_dir"] = str(inputs / UNCERTAINTY_DIR)
    if (inputs / UNCERTAINTY_SETTINGS).is_file():
        params["uncertainty_settings"] = str(inputs / UNCERTAINTY_SETTINGS)
    if (inputs / ESR_DIR).is_dir():
        params["esr_dir"] = str(inputs / ESR_DIR)
    return params


def build_params(job: Job, inputs: Path, outdir: Path) -> dict[str, Any]:
    params = {key: value for key, value in job.params.items() if key in PARAM_ALLOWLIST}
    limits = limits_for(job.owner)
    params |= staged_inputs(inputs)
    params |= {
        "outdir": str(outdir),
        "max_sequences": limits.max_sequences,
        "max_length": limits.max_length,
    }
    return params


def model_count(params: dict[str, Any]) -> int:
    fixed = len(as_list(params.get("models"))) if params.get("fixed_models") else 0
    return int(params.get("modelfinder", True)) + fixed


def codon_model_count(params: dict[str, Any]) -> int:
    """Codon model keys (codon ModelFinder and each codon model): each a combination on
    the codon alignment (coding DNA only)."""
    if not coding(params):
        return 0
    return bool(params.get("codon_modelfinder")) + len(as_list(params.get("codon_models")))


def builders(params: dict[str, Any]) -> list[str]:
    return [t for t in as_list(params.get("tree_builders", "iqtree")) if t != "none"]


def codon_trees(params: dict[str, Any]) -> list[str]:
    """How codon combinations get trees from IQ-TREE and RAxML-NG: 'search' (IQ-TREE
    under the codon model) and/or 'amino' (fitted on amino-acid trees)."""
    value = params.get("codon_trees")
    return [t for t in as_list("search" if value is None else value) if t != "none"]


def codon_search(params: dict[str, Any]) -> bool:
    """Whether codon models get their own IQ-TREE search."""
    return "search" in codon_trees(params) and "iqtree" in builders(params)


def codon_amino_trees(params: dict[str, Any]) -> int:
    """Amino-acid trees each codon model is fitted on, per alignment (FIT_CODON_AA): per
    tree source (IQ-TREE, RAxML-NG) the best by BIC, or every amino-acid model's."""
    if "amino" not in codon_trees(params):
        return 0
    sources = len([t for t in builders(params) if t in ("iqtree", "raxml")])
    return sources * (model_count(params) if params.get("codon_amino_trees") == "all" else 1)


def codon_tree_count(params: dict[str, Any]) -> int:
    """Trees per codon model and alignment: its own search, fits on amino-acid trees, and
    a fit to each provided tree."""
    return codon_search(params) + codon_amino_trees(params) + provided(params, "tree")


def codon_modelfinder(params: dict[str, Any]) -> bool:
    return coding(params) and bool(params.get("codon_modelfinder"))


def codeml_combinations(params: dict[str, Any]) -> int:
    """Codon combinations codeml runs on: each alignment's best, or every one."""
    if not as_list(params.get("codeml_models")) or not codon_model_count(params):
        return 0
    if params.get("stop_after", "ancestors") != "ancestors":
        return 0
    if params.get("codeml_trees", "best") == "all":
        return alignment_count(params) * codon_tree_count(params) * codon_model_count(params)
    return alignment_count(params)


def codeml_runs(params: dict[str, Any]) -> list[str]:
    """One entry per CODEML task of one combination: each model once per starting omega."""
    from phylobot import codeml

    return [
        key
        for key in as_list(params.get("codeml_models"))
        if key in codeml.MODELS
        for _ in codeml.MODELS[key].starts
    ]


def candidate_models(params: dict[str, Any]) -> list[str]:
    """ModelFinder's candidates: `candidates`, or when unset `models` without the mixtures
    (jobs submitted before the two had separate grids; the pipeline falls back the same
    way)."""
    value = params.get("candidates")
    return as_list(
        catalog.candidate_models(as_list(params.get("models"))) if value is None else value
    )


def provided(params: dict[str, Any], kind: str) -> int:
    """How many of the user's own alignments or trees: staged files (params from
    build_params) or their labels (Job.params)."""
    return len(params.get(kind) or params.get(f"{kind}_labels") or [])


def coding(params: dict[str, Any]) -> bool:
    """A coding-DNA job (--seqtype codon)."""
    return params.get("seqtype") == "codon"


def alignment_count(params: dict[str, Any]) -> int:
    aligners = len([a for a in as_list(params.get("aligners", "mafft")) if a != "none"])
    return aligners + provided(params, "alignment")


def alignment_keys(params: dict[str, Any]) -> list[str]:
    """Every alignment's key, in the pipeline's order: the given ones (an inherited
    aligner's key from `alignment_keys`, else provided1... counted among the user's
    own), then the aligners'."""
    keys = as_list(params.get("alignment_keys"))
    out, own = [], 0
    for index in range(provided(params, "alignment")):
        key = keys[index] if index < len(keys) else "provided"
        if key == "provided":
            own += 1
            key = f"provided{own}"
        out.append(key)
    return out + [a for a in as_list(params.get("aligners", "mafft")) if a != "none"]


def planned_combinations(params: dict[str, Any]) -> list[tuple[str, str, str]] | None:
    """An extension that adds models or tree sources (phylobot.extend; amino-acid models
    only): every (alignment, tree source, model) but the extended job's, which it reads
    rather than reruns (`skip_combinations`, as main.nf skips them). None for any other
    job, whose combinations are the whole cross product."""
    skipped = set(as_list(params.get("skip_combinations")))
    if not skipped:
        return None
    sources = [t for t in builders(params) if t in ("iqtree", "raxml")]
    sources += [f"provided{n}" for n in range(1, provided(params, "tree") + 1)]
    models = ["modelfinder"] if params.get("modelfinder", True) else []
    models += as_list(params.get("models")) if params.get("fixed_models") else []
    return [
        (a, t, m)
        for a in alignment_keys(params)
        for t in sources
        for m in models
        if f"{a}/{t}/{m}" not in skipped
    ]


def modelfinder_runs(params: dict[str, Any], planned: list[tuple[str, str, str]]) -> int:
    """MODELFINDER runs of an extension: alignments with a ModelFinder combination still
    to run whose pick is not reused from the extended job."""
    reused = set(as_list(params.get("modelfinder_reused")))
    return len({a for a, _t, m in planned if m == "modelfinder"} - reused)


def combination_count(params: dict[str, Any]) -> int:
    """Alignments x tree sources x models (amino-acid models, and codon models with
    their own tree sources: codon_tree_count); 0 for a job that stops after alignments.
    An extension's are those it runs (planned_combinations)."""
    if params.get("stop_after") == "alignments":
        return 0
    if (planned := planned_combinations(params)) is not None:
        return len(planned)
    trees = len(builders(params)) + provided(params, "tree")
    codon = codon_tree_count(params) * codon_model_count(params)
    return alignment_count(params) * (trees * model_count(params) + codon)


def reconciles(params: dict[str, Any]) -> bool:
    """A job reconciled with a species tree: SPECIES_TREE once, then RECONCILE and
    RECONCILE_EXPORT per combination. Its Job.params say so (`reconcile`); the params
    built for Nextflow (build_params) carry the staged species_map instead."""
    staged = params.get("reconcile") or params.get("species_map")
    return bool(staged) and combination_count(params) > 0


def branch_site(params: dict[str, Any]) -> bool:
    """A branch-site test job (a follow-up from the tree viewer; phylobot.branch_site):
    VALIDATE, CODEML per model and start, BRANCH_SITE_EXPORT, nothing else."""
    return bool(params.get("branch_site"))


def branch_site_runs() -> list[str]:
    """One entry per CODEML task of a branch-site job: each model once per start."""
    from phylobot import branch_site as test

    return [key for key, starts in test.STARTS.items() for _ in starts]


def properties_only(params: dict[str, Any]) -> bool:
    """A properties-only job (a follow-up from the ancestor viewer):
    one PROPERTIES task per combination, nothing else."""
    return bool(params.get("properties_only"))


def uncertainty_only(params: dict[str, Any]) -> bool:
    """An uncertainty run (a follow-up from a job's Uncertainty page):
    one UNCERTAINTY task, nothing else."""
    return bool(params.get("uncertainty_only"))


def property_targets(params: dict[str, Any]) -> int:
    """How many combinations such a job measures — its task count.

    The plan is written twice, because the two callers hold different shapes of params:
    the job page and the ETA have the `Job`'s own `properties_of`, while the worker has
    only the filtered pipeline params and the staged file they point at. Either way the
    answer is the number of combinations, and an unreadable plan counts as one rather
    than as none, so progress can never read "done 1 of 0".
    """
    plan = params.get("properties_of")
    if isinstance(plan, dict) and isinstance(plan.get("combinations"), dict):
        return max(len(plan["combinations"]), 1)
    path = params.get("property_targets")
    if isinstance(path, str):
        try:
            rows = Path(path).read_text().splitlines()[1:]
        except OSError:
            return 1
        return max(len({row.split("\t")[0] for row in rows if row.strip()}), 1)
    return 1


def esr(params: dict[str, Any]) -> bool:
    """An ESR job (phylobot.esr_views): ESR_FIT, ESR_TIP per held-out sequence,
    ESR_COLLECT, nothing else."""
    return bool(params.get("esr"))


def esr_tips(params: dict[str, Any]) -> int:
    """How many sequences an ESR job holds out: from the Job's own `esr_of`, or (the
    worker's filtered params) the staged settings; at least one, as `property_targets`."""
    plan = params.get("esr_of")
    if isinstance(plan, dict) and isinstance(plan.get("tips"), int):
        return max(int(plan["tips"]), 1)
    folder = params.get("esr_dir")
    if isinstance(folder, str):
        try:
            tips = json.loads((Path(folder) / "settings.json").read_text()).get("tips")
        except OSError, ValueError, AttributeError:
            return 1
        return max(len(set(tips)), 1) if isinstance(tips, list) else 1
    return 1


def simulation(params: dict[str, Any]) -> bool:
    """A simulation job (phylobot.simulations): SIM_PREPARE, SIM_ALISIM and SIM_COLLECT
    per replicate, nothing else."""
    return bool(params.get("simulate"))


def simulation_seeds(params: dict[str, Any]) -> list[int]:
    seeds = ((params.get("spec") or {}).get("simulation") or {}).get("seeds") or []
    return [s for s in seeds if isinstance(s, int)]


def expected_tasks(params: dict[str, Any]) -> int:
    """Tasks the pipeline will run: VALIDATE; each aligner (and for coding DNA its
    CODON_ALIGNMENT); an alignment export per alignment; with trees and ModelFinder,
    one MODELFINDER per alignment; per combination a tree (or fit; a RAxML-NG search
    and then a fit) and an export, plus with ancestors a reconstruction, ML/PRANK/GRASP
    gap placement if chosen, and rooting with an outgroup; with a species tree,
    SPECIES_TREE and per combination RECONCILE and RECONCILE_EXPORT, and with GeneRax as
    a tree source its search, fit and topology test per alignment x model; then SUMMARY
    and, with combinations, ANCESTOR_MAP. A
    branch-site job: VALIDATE, its codeml runs and BRANCH_SITE_EXPORT. A properties-only
    job: one PROPERTIES per combination. An uncertainty run: one
    UNCERTAINTY. An ESR job: ESR_FIT, an ESR_TIP per held-out sequence, ESR_COLLECT."""
    if branch_site(params):
        return 1 + len(branch_site_runs()) + 1
    if properties_only(params):
        return property_targets(params)
    if uncertainty_only(params):
        return 1
    if simulation(params):
        return 3 * len(simulation_seeds(params))
    if esr(params):
        return 2 + esr_tips(params)
    if planned_combinations(params) is not None:  # counted step by step, as main.nf skips
        from phylobot.stages import expected_kinds  # (stages imports this module)

        return sum(expected_kinds(params).values())
    aligners = len([a for a in as_list(params.get("aligners", "mafft")) if a != "none"])
    per_combination = 2
    if params.get("stop_after", "ancestors") == "ancestors":
        # parsimony and fitch run inside EXPORT; ml, prank and grasp add a step each
        extra = [m for m in as_list(params.get("gap_methods")) if m in ("ml", "prank", "grasp")]
        rooting = 1 if as_list(params.get("outgroup")) else 0  # ROOT_TREE
        per_combination = 3 + len(extra) + rooting
    exports = alignment_count(params)
    combinations = combination_count(params)
    modelfinder = alignment_count(params) if combinations and params.get("modelfinder", True) else 0
    if combinations and codon_modelfinder(params):
        modelfinder += alignment_count(params)  # MODELFINDER_CODON
    raxml = "raxml" in as_list(params.get("tree_builders", "iqtree"))
    searches = alignment_count(params) * model_count(params) if combinations and raxml else 0
    # GeneRax's tree: TREE_GENERAX and TOPOLOGY_TEST per alignment x model (its fit,
    # FIT_GENERAX, is the combination's tree step)
    if combinations and "generax" in builders(params):
        searches += 2 * alignment_count(params) * model_count(params)
    codon_steps = aligners if coding(params) else 0  # CODON_ALIGNMENT
    codeml = codeml_combinations(params) * (len(codeml_runs(params)) + 1)  # + CODEML_EXPORT
    reconcile = 1 + 2 * combinations if reconciles(params) else 0
    return (
        1 + aligners + codon_steps + exports + modelfinder + searches
        + per_combination * combinations + codeml + reconcile + 1 + (1 if combinations else 0)
    )  # fmt: skip


WORKFLOW_PREFIX = re.compile(r"^(?:[A-Z][A-Z0-9_]*:)+")


def task_name(name: str) -> str:
    """A task's name without its workflow: main.nf runs PHYLOBOT or BRANCH_SITE, so
    Nextflow names tasks "PHYLOBOT:ALIGN_MAFFT (mafft)"; everything here reads
    "ALIGN_MAFFT (mafft)"."""
    return WORKFLOW_PREFIX.sub("", name)


def read_trace(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if "name" in row:
            row["name"] = task_name(row["name"])
    return rows


def cpu_seconds(rows: list[dict[str, str]]) -> float:
    """Compute used: each step's runtime x its CPUs (steps reused from an earlier attempt
    are not counted: that attempt ran them)."""
    return sum(step_cpu(row) for row in rows if row.get("status") != "CACHED")


def step_cpu(row: dict[str, str]) -> float:
    runtime = seconds(row.get("realtime", ""))
    return 0.0 if runtime is None else runtime * int(row.get("cpus") or 1)


@dataclass(frozen=True)
class Usage:
    """An attempt's compute from its trace, split as billing needs it."""

    completed: float  # CPU seconds of steps that ran to completion
    failed: float  # of steps that failed or were stopped
    run: int  # steps run (completed, failed or stopped)
    reused: int  # steps reused from an earlier attempt (CACHED; charged there)


def usage(rows: list[dict[str, str]], ended: datetime | None = None) -> Usage:
    """`ended`: when the run stopped (the trace's last write). A step stopped part-way has
    a start but no runtime; its CPU time is estimated as `ended` minus its start."""
    completed = [row for row in rows if row.get("status") == "COMPLETED"]
    failed = [row for row in rows if row.get("status") in ("FAILED", "ABORTED")]
    return Usage(
        completed=sum(map(step_cpu, completed)),
        failed=sum(step_cpu(row) or stopped_cpu(row, ended) for row in failed),
        run=len(completed) + len(failed),
        reused=sum(row.get("status") == "CACHED" for row in rows),
    )


def stopped_cpu(row: dict[str, str], ended: datetime | None) -> float:
    start = trace_time(row.get("start", ""))
    if not start or ended is None:
        return 0.0
    elapsed = (ended - datetime.fromisoformat(start)).total_seconds()
    return max(elapsed, 0.0) * int(row.get("cpus") or 1)


SUBMITTED = re.compile(r"^\[([0-9a-f]{2}/[0-9a-f]{6})\] Submitted process > (.+)$")


def trace_time(text: str) -> str:
    """A trace timestamp ('2026-09-29 14:32:05.123'; UTC, as the runner sets TZ) as ISO
    8601 with its zone; '' if absent."""
    try:
        return datetime.strptime(text, "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=UTC).isoformat()
    except ValueError:
        return ""


def submitted_tasks(log: Path) -> dict[str, str]:
    """Tasks Nextflow has started, from its console log: {hash: task name}. The trace only
    lists tasks once they finish, so these are how running steps are seen."""
    if not log.is_file():
        return {}
    found = (SUBMITTED.match(line) for line in log.read_text(errors="replace").splitlines())
    return {match[1]: task_name(match[2]) for match in found if match}


def progress_from_trace(
    rows: list[dict[str, str]], total: int, running: dict[str, tuple[str, str]] | None = None
) -> dict[str, Any]:
    """Progress for the job page. `running`: {hash: (task name, first seen, ISO)} for
    started tasks; those not yet in the trace are listed with status RUNNING."""
    done = [row for row in rows if row.get("status") in DONE]
    traced = {row.get("hash") for row in rows}
    return {
        "done": len(done),
        "total": max(total, len(done)),
        "last": done[-1]["name"] if done else "",
        "tasks": [
            {
                "name": row["name"],
                "status": row["status"],
                "exit": row.get("exit", ""),
                "realtime": row.get("realtime", ""),
                "cpus": row.get("cpus", ""),
                # what it really used, against the 2 GB it reserved (for sizing steps)
                "peak_rss": row.get("peak_rss", ""),
                "start": trace_time(row.get("start", "")),
                "complete": trace_time(row.get("complete", "")),
            }
            for row in rows
        ]
        + [
            {"name": name, "status": "RUNNING", "start": seen}
            for key, (name, seen) in (running or {}).items()
            if key not in traced
        ],
    }


def user_error(rows: list[dict[str, str]]) -> str:
    """The message to show the user: VALIDATE's input error if that is what failed."""
    for row in rows:
        if row.get("name", "").startswith("VALIDATE") and row.get("status") == "FAILED":
            err = Path(row.get("workdir", "")) / ".command.err"
            if err.is_file():
                lines = [ln for ln in err.read_text().splitlines() if ln.startswith("Input error:")]
                if lines:
                    return lines[-1]
    return GENERIC_FAILURE


def log_name(task: str) -> str:
    """A task's stored log file name, without .log/.err: "MODELFINDER_mafft"."""
    return re.sub(r"[^A-Za-z0-9.+-]+", "_", task or "task").strip("_")


def collect_task_logs(rows: list[dict[str, str]], dest: Path) -> None:
    """Copy each task's stdout/stderr next to the results, named after the task."""
    dest.mkdir(parents=True, exist_ok=True)
    for row in rows:
        workdir = Path(row.get("workdir", ""))
        name = log_name(row.get("name", "task"))
        for stream in (".command.log", ".command.err"):
            source = workdir / stream
            if source.is_file() and source.stat().st_size:
                shutil.copyfile(source, dest / f"{name}{stream.replace('.command', '')}")


class Cancelled(Exception):
    """The job stopped being ours to run (cancelled or deleted)."""


RunCommand = Callable[[list[str], Path, Callable[[], bool]], int]


def run_command(command: list[str], log: Path, tick: Callable[[], bool]) -> int:
    """Run `command`, sending output to `log`; call `tick` every heartbeat interval.
    Stops the process if `tick` returns False (raises Cancelled) or if it exceeds the
    job time limit (raises TimeoutError). The caller stops any container it started."""
    limit = site_settings.job_time_limit().total_seconds()
    interval = settings.PHYLOBOT_HEARTBEAT_SECONDS
    elapsed = 0.0
    with open(log, "w") as out:
        process = subprocess.Popen(command, stdout=out, stderr=subprocess.STDOUT)
        while True:
            try:
                return process.wait(timeout=interval)
            except subprocess.TimeoutExpired:
                elapsed += interval
                keep_going = tick()
                if keep_going and elapsed <= limit:
                    continue
                process.kill()
                process.wait()
                if not keep_going:
                    raise Cancelled from None
                raise TimeoutError from None


class LocalDockerExecutor:
    """Runs Nextflow in the pinned runner image via the host Docker socket.

    The project is mounted at the same absolute path (PROJECT_ROOT) in the worker, the
    runner and every task container, because Nextflow hands those paths to Docker.
    """

    def __init__(self, run: RunCommand = run_command) -> None:
        self.root = Path(settings.PHYLOBOT_PROJECT_ROOT)
        self.run_command = run

    @staticmethod
    def container(job: Job) -> str:
        """Per attempt, so a superseded attempt's container never blocks or is stopped
        in place of the current one."""
        return f"phylobot-job-{job.id}-{job.attempt}"

    def stop(self, job: Job) -> None:
        """Stop the job's Nextflow container. Killing the `docker run` client alone leaves
        the container (and its pipeline) running."""
        subprocess.run(["docker", "stop", self.container(job)], capture_output=True)

    def stop_earlier(self, job: Job) -> None:
        """Stop any container an earlier attempt left running (one whose worker slept
        through being declared stopped): two Nextflows must never share a work directory."""
        names = [f"phylobot-job-{job.id}"]  # before containers were named per attempt
        names += [f"phylobot-job-{job.id}-{n}" for n in range(1, job.attempt)]
        # Task containers outlive a Nextflow killed outright; they run in the job's folder.
        listed = subprocess.run(
            ["docker", "ps", "-q", "--filter", "name=nxf-"], capture_output=True, text=True
        )
        if ids := listed.stdout.split():
            inspected = subprocess.run(
                ["docker", "inspect", "--format", "{{.Id}} {{.Config.WorkingDir}}", *ids],
                capture_output=True, text=True,
            )  # fmt: skip
            names += orphans(inspected.stdout, self.root / "jobs" / str(job.id))
        subprocess.run(["docker", "stop", *names], capture_output=True)

    def nextflow_command(self, job: Job, workdir: Path, resume: bool = False) -> list[str]:
        """Launched from the job's own folder, so Nextflow's resume cache (.nextflow/) is
        the job's alone; `resume` (a later attempt) reuses the steps that finished."""
        root = str(self.root)
        command = [
            "docker", "run", "--rm",
            "--name", self.container(job),
            "--platform", "linux/amd64",
            "-v", f"{root}:{root}",
            "-v", "/var/run/docker.sock:/var/run/docker.sock",
            "-w", str(workdir),
            "-e", f"NXF_HOME={root}/.nextflow-home",
            "-e", "NXF_ANSI_LOG=false",
            "-e", "TZ=UTC",  # trace timestamps are local time: make them UTC
            # the image's Docker CLI (19.03) takes this but not --platform (nextflow.config)
            "-e", "DOCKER_DEFAULT_PLATFORM=linux/amd64",
            # the CPUs task containers may use (nextflow.config's docker profile)
            "-e", f"PHYLOBOT_DOCKER_CPUSET={settings.PHYLOBOT_DOCKER_CPUSET}",
            settings.PHYLOBOT_NEXTFLOW_IMAGE,
            "nextflow", "-log", str(workdir / "nextflow.log"),
            "run", f"{root}/pipeline/main.nf", "-profile", "docker",
            "-params-file", str(workdir / "params.json"),
            "-work-dir", str(workdir / "work"),
            # this attempt's share of the machine (phylobot.machine), written by run()
            "-c", str(workdir / "resources.config"),
        ]  # fmt: skip
        if resume:
            command.append("-resume")
        return command

    def run(self, job: Job) -> Outcome:
        workdir = self.root / "jobs" / str(job.id)
        inputs, results = workdir / "inputs", workdir / "results"
        trace = results / "pipeline_info" / "trace.tsv"
        log = workdir / "pipeline.log"
        attempt = job.attempt
        if attempt > 1 and not (workdir / "params.json").is_file():
            work_archive.restore(self.root, job)  # another machine's, or a lost disk
        if attempt > 1 and (workdir / "params.json").is_file():
            # A resume: everything stays where Nextflow left it. Inputs are not fetched
            # again and params.json is not rebuilt, since a new file time or a changed
            # limit would change every step's cache key and nothing would be reused.
            self.stop_earlier(job)
            settle_previous(job, workdir, trace)
        else:
            if workdir.exists():
                shutil.rmtree(workdir)
            inputs.mkdir(parents=True)
            storage.download_prefix(f"{job.storage_prefix}inputs/", inputs)
            (workdir / "params.json").write_text(
                json.dumps(build_params(job, inputs, results), indent=2)
            )
        params = json.loads((workdir / "params.json").read_text())
        total = expected_tasks(params)
        # Only a folder Nextflow was launched from holds a resume cache (older jobs ran
        # from the project root): without one, -resume would reuse nothing anyway.
        resume = attempt > 1 and (workdir / ".nextflow").is_dir()
        Job.objects.filter(id=job.id).update(resume_cache=True)
        # Worked out at every launch, a resume's too: it may run on another machine, or
        # beside other jobs. Outside params.json, so no cache key changes.
        job_share = machine.for_launch()
        (workdir / "resources.config").write_text(machine.nextflow_config(job_share))
        JobAttempt.objects.filter(job_id=job.id, number=attempt).update(
            cpus=job_share.cpus, memory_bytes=job_share.memory_bytes
        )
        if job_share.warning:
            logger.warning("job %s: %s", job.id, job_share.warning)

        started: dict[str, tuple[str, str]] = {}  # hash -> (task name, first seen)

        def progress(rows: list[dict[str, str]]) -> dict[str, Any]:
            now = datetime.now(UTC).isoformat()
            for key, name in submitted_tasks(log).items():
                started.setdefault(key, (name, now))
            return progress_from_trace(rows, total, started)

        def tick() -> bool:
            """Heartbeat; False once the job is no longer this attempt's to run (cancelled,
            deleted, or declared stopped and resumed elsewhere)."""
            if not jobs.heartbeat(job.id, progress(read_trace(trace)), attempt):
                return False
            if log.is_file():  # live log for the job page; replaced by the final copy
                try:
                    storage.put_file(f"{job.storage_prefix}results/logs/pipeline.log", log)
                except Exception:
                    logger.warning("job %s: could not upload live log", job.id, exc_info=True)
            return True

        timed_out = False
        try:
            if not tick():
                raise Cancelled
            command = self.nextflow_command(job, workdir, resume)
            code = self.run_command(command, log, tick)
        except Cancelled:
            self.stop(job)
            # Only a deleted job's folder goes: a superseded attempt's folder is the
            # current attempt's too (the worker's maintenance removes deleted jobs' later).
            now = Job.objects.filter(id=job.id).values("deleted_at", "attempt").first()
            if now and now["deleted_at"] is not None and now["attempt"] == attempt:
                shutil.rmtree(workdir, ignore_errors=True)
            return Outcome(Job.Status.CANCELLED)
        except TimeoutError:
            self.stop(job)
            code, timed_out = -1, True

        rows = read_trace(trace)
        jobs.heartbeat(job.id, progress(rows), attempt)
        record(job.id, attempt, rows, trace)
        results.mkdir(parents=True, exist_ok=True)
        collect_task_logs(rows, results / "logs" / "tasks")
        for name in ("pipeline.log", "nextflow.log", "params.json"):
            if (workdir / name).is_file():
                shutil.copyfile(workdir / name, results / "logs" / name)
        manifest.unindex(job.id)  # until the upload is whole (maintenance lists it otherwise)
        try:  # heartbeats while it uploads; stops if the job stops being this attempt's
            storage.upload_tree(
                results,
                f"{job.storage_prefix}results/",
                alive=lambda: jobs.heartbeat(job.id, None, attempt),
                every=settings.PHYLOBOT_HEARTBEAT_SECONDS,
            )
        except storage.UploadStopped:
            logger.info("job %s: upload stopped: no longer attempt %d's", job.id, attempt)
            return Outcome(Job.Status.CANCELLED)
        manifest.indexed(job.id, manifest.directory_sizes(results))

        if timed_out:
            limit = site_settings.job_time_limit()
            return Outcome(Job.Status.FAILED, f"The job exceeded the time limit ({limit}).")
        if code != 0:
            return Outcome(Job.Status.FAILED, user_error(rows))
        shutil.rmtree(workdir, ignore_errors=True)  # all in storage now; kept on failure
        return Outcome(Job.Status.SUCCEEDED)


def orphans(inspected: str, folder: Path) -> list[str]:
    """Container ids from `docker inspect` lines ("<id> <working dir>") working inside
    `folder`, a job's own: only that job's tasks, never another job's or a golden run's."""
    found = []
    for line in inspected.splitlines():
        container, _, workdir = line.partition(" ")
        if workdir and Path(workdir).is_relative_to(folder):
            found.append(container)
    return found


def record(job_id: Any, attempt: int, rows: list[dict[str, str]], trace: Path) -> None:
    ended = datetime.fromtimestamp(trace.stat().st_mtime, UTC) if trace.is_file() else None
    found = usage(rows, ended)
    jobs.record_usage(job_id, attempt, found.completed, found.failed, found.run, found.reused)


def record_unrecorded(root: Path, within: timedelta = timedelta(days=7)) -> int:
    """Record the compute of recently ended attempts whose worker never could (the stale
    sweep failed them: the machine went down, or the worker was killed), from the trace
    their kept folder still holds. Worker maintenance."""
    since = datetime.now(UTC) - within
    pending = JobAttempt.objects.filter(
        completed_cpu_seconds__isnull=True, finished_at__gte=since, status__in=jobs.TERMINAL
    ).select_related("job")[:50]
    count = 0
    for attempt in pending:
        if attempt.number != attempt.job.attempt:
            continue  # its folder now holds a later attempt's trace (settle_previous's case)
        trace = root / "jobs" / str(attempt.job_id) / "results" / "pipeline_info" / "trace.tsv"
        if trace.is_file():
            record(attempt.job_id, attempt.number, read_trace(trace), trace)
            count += 1
    return count


def folder_bytes(folder: Path) -> int:
    """Bytes a folder's files take, without following links (Nextflow stages inputs as
    symbolic links to other steps' files, which would be counted twice)."""
    total = 0
    for directory, _, files in os.walk(folder):
        for name in files:
            try:
                total += (Path(directory) / name).lstat().st_size
            except OSError:
                pass  # removed meanwhile
    return total


def kept_folders(root: Path, running: set[str]) -> dict[str, Path]:
    """Failed jobs' local folders, by job id, that no attempt is running in."""
    folders = {p.name: p for p in (root / "jobs").glob("*") if p.is_dir()}
    return {name: path for name, path in folders.items() if name not in running}


def measure_work_folders(root: Path, running: set[str], limit: int = 5) -> int:
    """Record the size of failed jobs' kept folders not measured yet (`Job.work_bytes`;
    the Settings page shows the total). Worker maintenance."""
    folders = kept_folders(root, running)
    pending = Job.objects.filter(
        status=jobs.Status.FAILED, deleted_at__isnull=True, work_bytes__isnull=True,
        work_expired_at__isnull=True, id__in=[n for n in folders if jobs.is_job_id(n)],
    ).values_list("id", "attempt")[:limit]  # fmt: skip
    count = 0
    for job_id, attempt in pending:
        size = folder_bytes(folders[str(job_id)])
        # Not if it was resumed meanwhile: the folder is the new attempt's now.
        count += Job.objects.filter(id=job_id, attempt=attempt, status=jobs.Status.FAILED).update(
            work_bytes=size
        )
    return count


def expire_work_folders(root: Path, running: set[str], now: datetime | None = None) -> int:
    """Remove failed jobs' kept folders, and their archives, older than the Settings
    page's retention period (site_settings.keep_failed_work), so a later resume starts
    from the beginning.
    Worker maintenance; nothing is removed while the period is unset."""
    keep = site_settings.keep_failed_work()
    if keep is None:
        return 0
    folders = kept_folders(root, running)
    cutoff = (now or datetime.now(UTC)) - keep
    expired = (
        Job.objects.filter(
            status=jobs.Status.FAILED, deleted_at__isnull=True, finished_at__lt=cutoff
        )
        .filter(
            Q(id__in=[n for n in folders if jobs.is_job_id(n)])
            | Q(work_archive_attempt__isnull=False)
        )
        .values_list("id", flat=True)
    )
    count = 0
    for job_id in expired:
        with transaction.atomic():
            # Locked, so a resume clicked meanwhile waits and then finds no cache.
            job = Job.objects.select_for_update().filter(id=job_id).first()
            if job is None or job.status != jobs.Status.FAILED or job.finished_at is None:
                continue
            if job.finished_at >= cutoff:
                continue
            if str(job_id) in folders:
                shutil.rmtree(folders[str(job_id)], ignore_errors=True)
            if job.work_archive_attempt is not None:
                work_archive.remove(job_id)
            job.resume_cache = False
            job.work_bytes = job.work_archive_attempt = job.work_archive_bytes = None
            job.work_expired_at = now or datetime.now(UTC)
            job.save(
                update_fields=[
                    "resume_cache", "work_bytes", "work_archive_attempt", "work_archive_bytes",
                    "work_expired_at",
                ]
            )  # fmt: skip
        logger.info("job %s: removed its work folder, kept for %s", job_id, keep)
        count += 1
    return count


def settle_previous(job: Job, workdir: Path, trace: Path) -> None:
    """Before a resume overwrites them: record the previous attempt's compute if its
    worker never could (the machine went down mid-run), and keep its logs and trace
    under results/logs/attempts/<n>/, uploaded with the rest."""
    previous = job.attempt - 1
    attempt = job.attempts.filter(number=previous).first()
    if attempt is not None and attempt.completed_cpu_seconds is None:
        record(job.id, previous, read_trace(trace), trace)
    keep = workdir / "results" / "logs" / "attempts" / str(previous)
    keep.mkdir(parents=True, exist_ok=True)
    logs = workdir / "results" / "logs"
    for item in list(logs.iterdir()):
        if item.name != "attempts":
            shutil.move(item, keep / item.name)
    for source in (workdir / "pipeline.log", workdir / "nextflow.log", trace):
        if source.is_file() and not (keep / source.name).exists():
            shutil.copyfile(source, keep / source.name)


def run_job(job: Job, executor: LocalDockerExecutor, stop: threading.Event) -> None:
    """Run a claimed job to completion and record the outcome (worker thread entry).
    If the job was deleted meanwhile, remove whatever files it left behind."""
    attempt = job.attempt  # the outcome is this attempt's, whatever has happened since
    try:
        outcome = executor.run(job)
    except Exception:
        logger.exception("job %s: executor error", job.id)
        outcome = Outcome(Job.Status.FAILED, GENERIC_FAILURE)
    if stop.is_set() and outcome.status != Job.Status.SUCCEEDED:
        outcome = Outcome(Job.Status.FAILED, INTERRUPTED)
    if not jobs.mark(job.id, outcome.status, outcome.error, attempt):
        job.refresh_from_db()
        if job.attempt != attempt:  # resumed meanwhile: the rest is the new attempt's
            logger.info("job %s: attempt %d superseded by %d", job.id, attempt, job.attempt)
            return
    logger.info("job %s: %s %s", job.id, outcome.status, outcome.error)
    job.refresh_from_db()
    if job.status == Job.Status.SUCCEEDED and job.deleted_at is None:
        try:  # into My alignments; never a reason to fail a finished job
            alignments.register_job(job)
        except Exception:
            logger.warning("job %s: could not register its alignments", job.id, exc_info=True)
    if job.deleted_at is not None:
        jobs.purge(job)
    else:
        jobs.measure_storage(job)
