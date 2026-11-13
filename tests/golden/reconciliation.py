"""Running and scoring a reconciled job on a gene-family replicate.

Shared by the golden test (test_reconcile.py) and a validation harness, so the pass/fail
test runs exactly what the validation measured: MAFFT, LG+G4 (the simulation model),
IQ-TREE's and GeneRax's trees, both reconciled with the species tree, scored by
`phylobot.evaluation.reconcile`.
"""

from collections.abc import Callable
from pathlib import Path
from typing import Any

from phylobot import reconcile
from phylobot.evaluation import reconcile as R
from phylobot.seqio import read_fasta
from phylobot.species.table import OK, Row, write

from .pipeline import ROOT, remove, run_pipeline

MODEL = "LG+G4"
SOURCES = ("iqtree", "generax")


def species_files(work: Path, genes: dict[str, str], newick: str) -> tuple[Path, Path]:
    """species.tsv (as the web tier writes it) and the species tree, into `work`."""
    rows = [Row(sequence=g, query=s, species=s, source="simulation", status=OK)
            for g, s in genes.items()]  # fmt: skip
    (work / "species.tsv").write_text(write(rows))
    (work / "species_tree.nwk").write_text(newick.strip() + "\n")
    return work / "species.tsv", work / "species_tree.nwk"


def genes(truth: Path) -> dict[str, str]:
    """Gene -> species, from a family replicate's species.tsv."""
    lines = (truth / "species.tsv").read_text().splitlines()[1:]
    return dict(line.split("\t") for line in lines)


def run(work: Path, fasta: Path, genes: dict[str, str], newick: str, seqtype: str,
        seed: int, reuse: bool = False, alrt: int = 0) -> Path:  # fmt: skip
    """The pipeline on one replicate into work/out/run. A finished run is reused only if
    `reuse` (the report's harness); otherwise it is run again from scratch. `alrt`:
    SH-aLRT replicates (0, none, as in the report)."""
    results = work / "out" / "run"
    if (results / "summary.json").exists():
        if reuse:
            return results
        remove(results)
    work.mkdir(parents=True, exist_ok=True)
    species_map, species_tree = species_files(work, genes, newick)

    def rel(path: Path) -> str:
        return str(path.resolve().relative_to(ROOT))

    params: dict[str, Any] = {
        "input": rel(fasta),
        "aligners": ["mafft"],
        "tree_builders": ["iqtree", "generax"],
        "modelfinder": False,
        "fixed_models": True,
        "models": [MODEL],
        "alrt": alrt,
        "bootstrap": 0,
        "gap_methods": [],
        "seed": seed,
        "seqtype": seqtype,
        "species_map": rel(species_map),
        "species_tree": rel(species_tree),
    }
    if seqtype == "codon":
        params["codon_modelfinder"] = False
    return run_pipeline(work / "out", None, params)


def score(truth_dir: Path, results: Path, sequences: bool = True) -> dict[str, dict[str, Any]]:
    """Per tree source: the tree, its duplications, chimeric ancestors and the
    paralog-mixing flags. `results` is the run's folder; `sequences` False skips
    ancestor identity (codon replicates' true ancestors are codons)."""
    folder = results / "mafft"
    events = truth_dir / "true_events.tsv"
    true = R.truth(
        (truth_dir / "true_tree.nwk").read_text(),
        events.read_text() if events.exists() else "event\tgene_node\n",
        read_fasta(truth_dir / "true_ancestors.fasta") if sequences else {},
    )
    au = {}
    au_file = folder / "topology_tests" / f"{MODEL}.tsv"
    if au_file.exists():
        for line in au_file.read_text().splitlines()[1:]:
            cells = line.split("\t")
            au[cells[0]] = cells[-1]

    def files(source: str) -> tuple[str, str]:
        combination = folder / source / MODEL
        return ((combination / "reconciliation.tsv").read_text(),
                (combination / "reconciliation.nwk").read_text())  # fmt: skip

    clades = {source: set(reconcile.rooted_clades(files(source)[1])[0].items())
              for source in SOURCES}  # fmt: skip
    scores = {}
    for source in SOURCES:
        combination = folder / source / MODEL
        ancestors = read_fasta(combination / "ancestors.fasta") if true.ancestors else None
        s = R.score(true, (combination / "tree.nwk").read_text(), files(source), ancestors)
        # the paralog-mixing flag against the *other* tree source, as the ancestor
        # viewer's disagreement note does it (GeneRax's for IQ-TREE's and back)
        other = "iqtree" if source == "generax" else "generax"
        flag = R.flags(true, files(source), files(other))
        # chimeric ancestors the other tree has too: the flag compares the two, so it
        # cannot fire there (the other source's max_chimeric catches them)
        theirs = {below for _, below in clades[other]}
        shared = [n for n, below in clades[source] if n in s.chimeric and below in theirs]
        scores[source] = {
            "rf": s.rf, "ancestors": s.ancestors, "dups_true": s.duplications_true,
            "dups_inferred": s.duplications_inferred, "dups_found": s.duplications_found,
            "precision": s.precision, "chimeric": len(s.chimeric),
            "chimeric_nodes": s.chimeric, "clades_matched": s.clades_matched,
            "identity": s.identity, "losses": s.losses, "flagged": len(flag.flagged),
            "flagged_nodes": flag.flagged, "flag_tp": flag.true_positives,
            "chimeric_shared": sorted(shared),
            "au_p": au.get(source, ""),
        }  # fmt: skip
    return scores


def _missed(m: dict[str, Any]) -> list[str]:
    """Chimeric ancestors the flag should have caught: not flagged, and not shared with
    the other tree (where the flag, a comparison of the two, cannot fire)."""
    shared = set(m.get("chimeric_shared", ()))
    return sorted(set(m["chimeric_nodes"]) - set(m["flagged_nodes"]) - shared)


# expected_metrics.json's keys: (passes, failure message), given a source's scores and
# the threshold
type Rule = tuple[Callable[[dict[str, Any], Any], bool], Callable[[dict[str, Any], Any], str]]
RULES: dict[str, Rule] = {
    "max_rf": (lambda m, x: m["rf"] <= x, lambda m, x: f"RF {m['rf']} > {x}"),
    "max_chimeric": (
        lambda m, x: m["chimeric"] <= x,
        lambda m, x: f"chimeric ancestors {m['chimeric_nodes']}",
    ),
    "duplications_found": (
        lambda m, x: not x or m["dups_found"] == m["dups_true"],
        lambda m, x: f"found {m['dups_found']} of {m['dups_true']} true duplications",
    ),
    "min_duplication_precision": (
        lambda m, x: (m["precision"] or 0.0) >= x,
        lambda m, x: f"duplication precision {m['precision']} < {x}",
    ),
    "chimeric_flagged": (
        lambda m, x: not x or not _missed(m),
        lambda m, x: f"chimeric ancestors not flagged: {_missed(m)}",
    ),
}


def check(scores: dict[str, dict[str, Any]], expected: dict[str, Any]) -> list[str]:
    """The owner's thresholds (expected_metrics.json: per tree source, `RULES`' keys)."""
    failures = []
    for source, thresholds in expected.items():
        if source.startswith("_"):
            continue
        for key, want in thresholds.items():
            passes, message = RULES[key]
            if not passes(scores[source], want):
                failures.append(f"{source}: {message(scores[source], want)}")
    return failures


def summary(scores: dict[str, dict[str, Any]]) -> str:
    lines = []
    for source, m in scores.items():
        identity = "-" if m["identity"] is None else f"{m['identity']:.3f}"
        lines.append(
            f"  {source:<8} RF {m['rf']:>2}  duplications {m['dups_found']}/{m['dups_true']} "
            f"found, {m['dups_inferred']} inferred  chimeric {m['chimeric']} "
            f"(flagged {m['flag_tp']}; {m['flagged']} flags)  identity {identity}  "
            f"AU p {m['au_p'] or '-'}"
        )
    return "\n".join(lines)
