"""Users' simulations: the Simulations page's jobs.

A simulation is a Job (it uses compute, so it is queued, limited and accounted like any
other) whose params hold `simulate: true` and the spec (`phylobot.simulate.spec`, the
same engine as the golden datasets); the pipeline runs it with --simulate. Its inputs
are spec.json, and when it is built on one of the user's jobs, that job's tree
(sim_tree.nwk) and an ancestor as the root sequence (sim_root.fasta).

Simulated data must never pass as real: every sequence is
named sim001, sim002, ... A tree from a job is renamed here before it is stored, with
the real names kept beside it (inputs/names.tsv, shown only on the simulation's own
page); a random tree's leaves are renamed by the engine. Each replicate's true alignment
joins My alignments flagged `simulated`, a flag no copy loses.
"""

import json
from dataclasses import dataclass
from typing import Any

from django.db import IntegrityError, transaction
from django.db.models import QuerySet

from phylobot import ancestors, catalog, job_types, manifest, storage
from phylobot.models import Alignment, Job
from phylobot.simulate import parse
from phylobot.simulate.engine import neutral_names
from phylobot.simulate.spec import SHAPES, RandomTree
from phylobot.trees import load_tree

# The pinned IQ-TREE image (pipeline/nextflow.config); recorded in the spec so the
# provenance says what simulated it.
IQTREE_IMAGE = "quay.io/biocontainers/iqtree:3.1.3--h8471819_0"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
MATRICES = [m.name for m in catalog.MATRICES]
RATES = {
    "": "equal rates",
    "G4": "+ gamma (4 categories)",
    "I+G4": "+ invariant sites + gamma",
}
SHAPE_HELP = {
    "yh": "random joins of lineages (Yule-Harding): the usual choice",
    "u": "every topology equally likely",
    "bd": "a birth-death process, with its rates",
    "bal": "as balanced as possible",
    "cat": "a ladder: each split leaves one sequence",
}
assert set(SHAPE_HELP) == set(SHAPES)


def is_simulation(job: Job) -> bool:
    return bool(job.params.get("simulate"))


def model_string(matrix: str, rate: str, alpha: float, invariant: float) -> str:
    """An AliSim model with its parameters, e.g. LG+I{0.1}+G4{0.5}."""
    if matrix not in MATRICES or rate not in RATES:
        raise ValueError(f"unknown model {matrix}+{rate}")
    model = matrix
    if "I" in rate:
        model += f"+I{{{invariant:g}}}"
    if "G4" in rate:
        model += f"+G4{{{alpha:g}}}"
    return model


@dataclass(frozen=True)
class TreeChoice:
    key: str  # the combination
    label: str
    path: str  # the tree file under the job's results/: the rooted one if there is one


def tree_choices(job: Job) -> list[TreeChoice]:
    """The trees a finished job offers to simulate on: one per combination."""
    prefix = f"{job.storage_prefix}results/"
    keys = manifest.paths(job)
    if "summary.json" not in keys:
        return []
    summary = json.loads(storage.read_bytes(prefix + "summary.json"))
    choices = []
    for combo in summary.get("combinations", []):
        key = combo["key"]
        for name in ("tree.rooted.nwk", "tree.nwk"):
            if f"{key}/{name}" in keys:
                choices.append(TreeChoice(key, combo.get("label", key), f"{key}/{name}"))
                break
    return choices


def jobs_with_trees(user: Any, text: str = "") -> QuerySet[Job]:
    """The user's finished analyses that may have trees, newest first; `text` searches
    their names and IDs (the job picker)."""
    found = Job.objects.for_user(user).filter(
        job_types.type_filter("trees"), status=Job.Status.SUCCEEDED
    )
    if text:
        found = found.filter(job_types.search(text))
    return found.order_by("-created_at")


def ancestor_versions(job: Job, combination: str) -> dict[str, dict[str, str]]:
    """{gap version: {node: aligned sequence}} for a combination's ancestors."""
    prefix = f"{job.storage_prefix}results/{combination}/"
    present = {k[len(combination) + 1 :] for k in manifest.paths(job, f"{combination}/")}
    out = {}
    for version, (_, filename) in ancestors.VERSIONS.items():
        if filename in present:
            text = storage.read_bytes(prefix + filename).decode(errors="replace")
            out[version] = ancestors.parse_fasta(text)
    return out


def default_version(job: Job, versions: dict[str, Any]) -> str:
    """The job's own default gap placement if it ran, else the first masked version, else
    IQ-TREE's unmasked reconstruction (an ancestor with every insertion)."""
    for version in [*(job.params.get("gap_methods") or []), "fitch", "parsimony", "ml"]:
        if version in versions:
            return str(version)
    return next(iter(versions), ancestors.ORIGINAL)


def rename_tree(newick: str) -> tuple[str, dict[str, str]]:
    """The tree with its leaves renamed sim001, ... (preorder), and {simulated: original}.
    Internal labels (Node1, ...) are kept: the true ancestors carry the job's names."""
    tree = load_tree(newick)
    leaves = [leaf.taxon.label for leaf in tree.leaf_node_iter()]
    renamed = neutral_names(leaves)
    for leaf in tree.leaf_node_iter():
        leaf.taxon.label = renamed[leaf.taxon.label]
    text = tree.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True)
    return text.strip() + "\n", {new: old for old, new in renamed.items()}


def leaf_count(newick: str) -> int:
    return len(load_tree(newick).leaf_nodes())


def clean_root(sequence: str) -> str:
    """An ungapped protein root sequence; ValueError if it holds anything else."""
    residues = "".join(sequence.split()).replace("-", "").upper()
    if not residues:
        raise ValueError("the root sequence is empty")
    if bad := sorted(set(residues) - AMINO_ACIDS):
        raise ValueError(f"the root sequence may hold the 20 amino acids only, not {''.join(bad)}")
    return residues


def build_spec(
    *,
    description: str,
    tree: dict[str, Any],
    root: dict[str, Any],
    model: str,
    indels: tuple[float, float, float] | None,
    seed: int,
    replicates: int,
) -> dict[str, Any]:
    """A checked spec (raises phylobot.simulate.SpecError). `tree` is the spec's [tree]
    section ({"random": ...} or {} for the staged file), `root` likewise."""
    spec: dict[str, Any] = {
        "description": description,
        "tree": tree,
        "root": root or {"name": "root"},
        "model": {"substitution": model},
        "simulation": {
            "seeds": list(range(seed, seed + replicates)),
            "image": IQTREE_IMAGE,
            "names": "neutral" if "random" in tree else "keep",
        },
    }
    if indels is not None:
        insertion, deletion, mean = indels
        spec["indels"] = {
            "rates": f"{insertion:g},{deletion:g}",
            "sizes": f"GEO{{{mean:g}}},GEO{{{mean:g}}}",
        }
    parse(spec)
    return spec


def residues_written(taxa: int, length: int, replicates: int) -> int:
    """What a simulation stores: every node's row (2n - 1 of them) per replicate."""
    return replicates * (2 * taxa - 1) * length


def describe(params: dict[str, Any]) -> list[str]:
    """The simulation in words, for the job page's plan line."""
    spec = parse(params["spec"])
    source = params.get("simulation_of") or {}
    if isinstance(spec.tree, RandomTree):
        tree = f"a random {SHAPES[spec.tree.shape]} tree of {spec.tree.taxa} sequences"
    else:
        tree = f"the tree of {source.get('label', 'a job')}"
    n = len(spec.seeds)
    steps = [f"simulate {n} replicate{'s' if n != 1 else ''} on {tree}"]
    if source.get("node"):
        steps.append(f"starting from its ancestor {source['node']}")
    model = spec.model.substitution or ""
    steps.append(f"under {model}" + (" with indels" if spec.indels else ", no indels"))
    return steps


@dataclass
class Replicate:
    seed: int
    sequences: int
    columns: int
    gaps: float  # fraction of the true alignment's cells
    files: dict[str, str]  # truth bundle file -> result path
    alignment: Alignment | None


TRUTH_FILES = ("input.fasta", "true_alignment.fasta", "true_tree.nwk", "true_ancestors.fasta")


def replicates(job: Job, keys: set[str]) -> list[Replicate]:
    """The finished replicates, in seed order, with their registered alignments."""
    records = {a.job_key: a for a in Alignment.objects.filter(job=job)}
    out = []
    for seed in (job.params.get("spec") or {}).get("simulation", {}).get("seeds", []):
        folder = f"replicates/{seed}"
        path = f"{folder}/true_alignment.fasta"
        if path not in keys:
            continue
        text = storage.read_bytes(f"{job.storage_prefix}results/{path}").decode(errors="replace")
        rows = list(ancestors.parse_fasta(text).values())
        cells = sum(len(r) for r in rows) or 1
        out.append(
            Replicate(
                seed=seed,
                sequences=len(rows),
                columns=len(rows[0]) if rows else 0,
                gaps=sum(r.count("-") for r in rows) / cells,
                files={f: f"{folder}/{f}" for f in TRUTH_FILES if f"{folder}/{f}" in keys},
                alignment=records.get(f"replicate{seed}"),
            )
        )
    return out


def register(job: Job) -> list[Alignment]:
    """Each replicate's true alignment in My alignments, flagged simulated. Idempotent."""
    if job.status != Job.Status.SUCCEEDED or not is_simulation(job):
        return []
    keys = manifest.paths(job, "replicates/")
    found = []
    for replicate in replicates(job, keys):
        if replicate.alignment is not None:
            found.append(replicate.alignment)
            continue
        key = f"replicate{replicate.seed}"
        try:
            with transaction.atomic():
                record = Alignment.objects.create(
                    owner_id=job.owner_id, name=f"Simulated replicate {replicate.seed}",
                    origin=Alignment.Origin.JOB, seqtype="protein", sequences=replicate.sequences,
                    columns=replicate.columns, job=job, job_key=key,
                    path=replicate.files["true_alignment.fasta"], simulated=True,
                )  # fmt: skip
        except IntegrityError:  # registered meanwhile
            record = Alignment.objects.get(job=job, job_key=key)
        found.append(record)
    return found


def original_names(job: Job) -> dict[str, str]:
    """{simulated name: the job's sequence name} for a simulation on a job's tree."""
    key = f"{job.storage_prefix}inputs/names.tsv"
    if key not in set(storage.list_keys(f"{job.storage_prefix}inputs/")):
        return {}
    lines = storage.read_bytes(key).decode().splitlines()[1:]
    return dict(line.split("\t", 1) for line in lines if "\t" in line)


def replicate_input(job: Job, seed: int) -> str:
    """A replicate's unaligned sequences (input.fasta) as text."""
    path = f"{job.storage_prefix}results/replicates/{int(seed)}/input.fasta"
    return storage.read_bytes(path).decode()


def root_sequence(text: str, node: str) -> str:
    """One ancestor's sequence from a FASTA text, ungapped (KeyError if absent)."""
    return clean_root(ancestors.parse_fasta(text)[node])
