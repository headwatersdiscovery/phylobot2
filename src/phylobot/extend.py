"""Extend a job: a new job that builds on a finished one's results and only adds to them.

A finished job's results never change: a paper cites it, its ZIP records every file's
SHA-256, libraries rebuild their members from its files. So "Extend" never reruns a job
in place. It opens New Job (`?extend=<job>`) with what the parent fixed locked, and
Submit queues a new Job (`Job.extends`). With the parent's seed, it gives what one job
with every setting from the start would have. Two modes, by how far the parent ran:

- **alignments** (B1; the parent stopped after its alignments): the parent's published
  alignments become the new job's inputs, each keeping its aligner's key
  (`alignment_keys`: VALIDATE checks it), so the new job names them as a job that ran
  those aligners itself would. The user adds aligners and goes on to trees and ancestors.
- **results** (B2; the parent has trees, maybe ancestors): the user adds aligners only,
  and every later setting is the parent's (`SubmitForm.inherit`). The new job runs just
  the new alignments; the parent's are read in place (`phylobot.joined`, params
  `inherits`), never copied, so its pages, viewers and ZIP show one job with every
  alignment. Deleting the parent takes those from it, which the delete dialog says.
  A reconciled parent is refused: its species tree is resolved from every tree, so new
  alignments would change its results.
"""

from dataclasses import dataclass, field
from typing import Any

from phylobot import catalog, executors, joined, storage
from phylobot.models import Job

PROVIDED = "provided"  # alignment_keys' entry for the user's own alignment
# Params an extension of a job with trees does not take from it (extended_params).
NOT_INHERITED = {
    "aligners", "alignment_labels", "alignment_keys", "tree_labels", "inherits", "extends",
    "skip_combinations", "modelfinder_reused",
}  # fmt: skip
ALIGNMENTS, RESULTS = "alignments", "results"  # the two modes


@dataclass(frozen=True)
class Alignment:
    key: str  # an aligner's key, or PROVIDED (renumbered provided1... in the new job)
    label: str
    content: bytes  # (results mode: only what the sequences may come from)
    folder: str = ""  # its folder in the parent's results: the key, or provided<n>


@dataclass
class Inheritance:
    """What an extension takes from its parent: the parent's own input sequences (or,
    if it had none, its first alignment's rows), its alignments, aligners' first (in the
    order the parent ran them, so its first aligner stays the primary one), then the
    user's own; in results mode, also its own trees (every tree source runs on every
    alignment, the new ones too)."""

    parent: Job
    fasta: bytes
    alignments: list[Alignment] = field(default_factory=list)
    mode: str = ALIGNMENTS
    trees: list[tuple[str, bytes]] = field(default_factory=list)  # (label, content)

    @property
    def aligners(self) -> list[str]:
        return [a.key for a in self.alignments if a.key != PROVIDED]

    @property
    def folders(self) -> list[str]:
        """Every inherited alignment's folder in the parent's results."""
        return [a.folder for a in self.alignments]


def mode(job: Job) -> str:
    return ALIGNMENTS if job.params.get("stop_after") == "alignments" else RESULTS


def refusal(job: Job) -> str:
    """Why `job` cannot be extended, or "" if it can."""
    params = job.params
    follow_up = params.get("properties_only") or params.get("uncertainty_only")
    if params.get("simulate") or params.get("branch_site") or params.get("esr") or follow_up:
        return "Only an analysis job can be extended."
    if job.imported_at is not None:
        return "An imported job's results were not computed here, so it can't be extended."
    if job.deleted_at is not None or job.purged_at is not None:
        return "This job's files are deleted."
    if job.status != Job.Status.SUCCEEDED:
        return "Only a job that finished can be extended."
    if mode(job) == RESULTS:
        if params.get("reconcile") or "generax" in (params.get("tree_builders") or []):
            return (
                "A reconciled job can't be extended: its species tree is resolved from "
                "every tree, so new alignments would change its results."
            )
        if joined.cached(job).lost:
            return "The job this one extends was deleted, so its results are incomplete."
    return ""


def inherit(job: Job) -> Inheritance:
    """What an extension of `job` inherits. Raises LookupError if files are missing."""
    found = joined.cached(job)  # an extension's parent may be an extension itself
    if not found.summary:
        raise LookupError("the job has no summary")
    codon = job.params.get("seqtype") == "codon"
    order, _ = counts(job)
    entries = sorted(
        found.summary.get("alignments") or [],
        key=lambda e: (
            str(e["key"]).startswith(PROVIDED),
            order.index(e["key"]) if e["key"] in order else len(order),
            int(str(e["key"])[len(PROVIDED) :] or 0) if str(e["key"]).startswith(PROVIDED) else 0,
        ),
    )
    sequences = f"{job.storage_prefix}inputs/sequences.fasta"
    fasta = b""
    inputs = set(storage.list_keys(f"{job.storage_prefix}inputs/"))
    if sequences in inputs:
        fasta = storage.read_bytes(sequences)
    results = mode(job) == RESULTS
    alignments: list[Alignment] = []
    for entry in entries:
        key = str(entry["key"])
        path = f"{key}/alignment.cds.fasta" if codon else f"{key}/alignment.fasta"
        if path not in found.keys:
            raise LookupError(f"the job has no {path}")
        needed = not results or (not fasta and not alignments)
        alignments.append(
            Alignment(
                key=PROVIDED if key.startswith(PROVIDED) else key,
                label=str(entry.get("label") or key),
                content=found.read_bytes(path) if needed else b"",
                folder=key,
            )
        )
    if not alignments:
        raise LookupError("the job has no alignments")
    out = Inheritance(parent=job, fasta=fasta, alignments=alignments, mode=mode(job))
    if results:
        if not fasta:  # the sequences come from the first alignment, as VALIDATE takes them
            out.fasta = ungapped(alignments[0].content)
        labels = job.params.get("tree_labels") or []
        for index, label in enumerate(labels, 1):
            key = f"{job.storage_prefix}inputs/tree{index}.nwk"
            if key not in inputs:
                raise LookupError(f"the job has no tree{index}.nwk")
            out.trees.append((label, storage.read_bytes(key)))
    return out


def ungapped(fasta: bytes) -> bytes:
    """An alignment's rows without gaps."""
    lines = []
    for line in fasta.decode(errors="replace").splitlines():
        lines.append(line if line.startswith(">") else line.replace("-", "").replace(".", ""))
    return ("\n".join(lines) + "\n").encode()


def counts(job: Job) -> tuple[list[str], int]:
    """The parent's aligners (those it inherited first, if it extended a job itself, then
    those it ran) and how many of its alignments were the user's own, from its params
    alone (the New Job summary asks on every change, so no storage reads)."""
    keys = job.params.get("alignment_keys") or []
    labels = job.params.get("alignment_labels") or []
    read = [str(k) for k in (job.params.get("inherits") or {}).get("alignments") or []]
    inherited = [k for k in keys if k != PROVIDED] + [k for k in read if not k.startswith(PROVIDED)]
    ran = [a for a in job.params.get("aligners") or [] if a not in ("none", *inherited)]
    own = sum(1 for k in keys if k == PROVIDED) if keys else len(labels)
    own += sum(1 for k in read if k.startswith(PROVIDED))
    return inherited + ran, own


@dataclass(frozen=True)
class Additions:
    """What an extension of a job with trees adds to every alignment, the inherited ones
    too: tree sources, fixed models, ModelFinder. Every other setting is the parent's."""

    builders: tuple[str, ...] = ()
    models: tuple[str, ...] = ()
    modelfinder: bool = False

    def __bool__(self) -> bool:
        return bool(self.builders or self.models or self.modelfinder)


def additions(parent: Job, builders: list[str], models: list[str], modelfinder: bool) -> Additions:
    """The posted choices that the parent lacks (never a coding-DNA parent's: codon
    combinations are not skipped yet, so they would all run again)."""
    if parent.params.get("seqtype") == "codon":
        return Additions()
    params = parent.params
    had = set(params.get("tree_builders") or [])
    fixed = set(params.get("models") or []) if params.get("fixed_models") else set()
    return Additions(
        builders=tuple(t for t in dict.fromkeys(builders) if t in TREE_SOURCES and t not in had),
        models=tuple(m for m in dict.fromkeys(models) if catalog.is_model(m) and m not in fixed),
        modelfinder=modelfinder and not params.get("modelfinder"),
    )


TREE_SOURCES = ("iqtree", "raxml")  # what an extension can add (not GeneRax: reconciliation)


def extended_params(
    inheritance: Inheritance, aligners: list[str], adds: Additions, trees: list[str]
) -> dict[str, Any]:
    """An extension (results mode)'s params: the parent's, but for the aligners (only the
    new ones run), what it reads from the parent (`inherits`), and the additions: with
    any, the parent's alignments come again as inputs under their keys, and
    `skip_combinations` and `modelfinder_reused` keep what the parent has from running
    again (main.nf). `trees` are the labels of the parent's own trees, which run too."""
    parent = inheritance.parent
    params = {k: v for k, v in parent.params.items() if k not in NOT_INHERITED}
    params["aligners"] = aligners or ["none"]
    if trees:
        params["tree_labels"] = trees
    params["inherits"] = {"job": str(parent.id), "alignments": inheritance.folders}
    params["extends"] = {  # display only
        "job": str(parent.id),
        "name": parent.name,
        "alignments": len(inheritance.alignments),
        "results": True,
    }
    if not adds:
        return params
    had = [t for t in parent.params.get("tree_builders") or [] if t != "none"]
    fixed = list(parent.params.get("models") or []) if parent.params.get("fixed_models") else []
    params["tree_builders"] = had + list(adds.builders)
    params["models"] = fixed + list(adds.models) or None
    params["fixed_models"] = bool(params["models"])
    params["modelfinder"] = bool(parent.params.get("modelfinder") or adds.modelfinder)
    if adds.modelfinder:
        params["candidates"] = None  # IQ-TREE's default set
    if "raxml" in adds.builders and params.get("candidates"):
        # As a job with RAxML-NG from the start: ModelFinder (new alignments only; the
        # inherited picks are reused) never picks a matrix RAxML-NG lacks.
        params["candidates"] = [m for m in params["candidates"] if not catalog.raxml_missing([m])]
    params["alignment_labels"] = [a.label for a in inheritance.alignments]
    params["alignment_keys"] = [a.key for a in inheritance.alignments]
    combos = joined.cached(parent).summary.get("combinations", [])
    params["skip_combinations"] = sorted(c["key"] for c in combos)
    planned = executors.planned_combinations(params) or []
    needing = {a for a, _t, m in planned if m == catalog.MODELFINDER}
    have = {c["aligner"] for c in combos if c["model"] == catalog.MODELFINDER}
    params["modelfinder_reused"] = sorted(needing & have)
    return params


def picks(inheritance: Inheritance, aligners: list[str]) -> dict[str, bytes]:
    """The parent's ModelFinder picks an extension reuses, as its input files
    (executors.PICKS_DIR): each aligner's model.txt and modelfinder.tsv, from one of its
    ModelFinder combinations (every tree source has the same pick)."""
    found = joined.cached(inheritance.parent)
    out = {}
    for aligner in aligners:
        combo = next(
            c["key"] for c in found.summary.get("combinations", [])
            if c["aligner"] == aligner and c["model"] == catalog.MODELFINDER
            and {f"{c['key']}/model.txt", f"{c['key']}/modelfinder.tsv"} <= found.keys
        )  # fmt: skip
        for name in ("model.txt", "modelfinder.tsv"):
            out[f"{executors.PICKS_DIR}/{aligner}.{name}"] = found.read_bytes(f"{combo}/{name}")
    return out


def unusable_picks(inheritance: Inheritance, aligners: list[str]) -> list[tuple[str, str]]:
    """(alignment, model) for each reused ModelFinder pick RAxML-NG can't run (a matrix
    it lacks): adding RAxML-NG would fail on it."""
    out = []
    for name, content in picks(inheritance, aligners).items():
        if name.endswith(".model.txt"):
            model = content.decode(errors="replace").strip()
            if catalog.raxml_missing([model]):
                out.append((name.split("/")[-1].removesuffix(".model.txt"), model))
    return out


def alignment_contents(inheritance: Inheritance) -> list[tuple[str, bytes]]:
    """(label, content) of every inherited alignment, as inputs (an extension with
    additions runs new combinations on them)."""
    found = joined.cached(inheritance.parent)
    codon = inheritance.parent.params.get("seqtype") == "codon"
    name = "alignment.cds.fasta" if codon else "alignment.fasta"
    return [(a.label, found.read_bytes(f"{a.folder}/{name}")) for a in inheritance.alignments]
