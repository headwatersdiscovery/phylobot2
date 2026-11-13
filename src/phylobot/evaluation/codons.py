"""Score codon reconstructions against a codon golden scenario's truth.

The software check, not the research comparison: on a scenario simulated under the kind
of model IQ-TREE fits (lysin/codon-easy), the codon pipeline's ancestors, IQ-TREE's and
codeml's, must come close to the true ones. Nodes are matched by the subtrees meeting at
each node (`gaps.match_nodes`), which ignores rooting: the truth is rooted at its
midpoint, the inferred trees are unrooted, so the true root has no counterpart and is
left out. Identity is position by position: codon-easy has no indels.

    uv run python -m phylobot.evaluation.codons TRUTH_DIR RESULTS_DIR [FOLDER]
"""

import json
import operator
import sys
from collections.abc import Callable
from pathlib import Path
from statistics import fmean
from typing import Any

from phylobot import codons, gaps
from phylobot.evaluation.sequences import identity
from phylobot.seqio import read_fasta


def codon_identity(inferred: str, truth: str) -> float:
    """Fraction of codons identical between two ungapped coding sequences."""
    a, b = inferred.replace("-", ""), truth.replace("-", "")
    if len(a) != len(b) or len(a) % 3:
        raise ValueError(f"length mismatch ({len(a)} vs {len(b)} nucleotides)")
    pairs = [(a[i : i + 3], b[i : i + 3]) for i in range(0, len(a), 3)]
    return sum(x == y for x, y in pairs) / len(pairs)


def score_codons(
    truth: Path, results: Path, code: codons.GeneticCode = codons.STANDARD, folder: str = ""
) -> dict[str, Any]:
    """Compare results/<folder>{ancestors.fasta, ancestors.cds.fasta} (results/tree.nwk
    for the node names) with truth/{true_tree.nwk, true_ancestors.fasta}. `folder` is ""
    for IQ-TREE's reconstruction, "codeml/<model>" for codeml's."""
    names = gaps.match_nodes(
        (truth / "true_tree.nwk").read_text(), (results / "tree.nwk").read_text()
    )
    true_ancestors = read_fasta(truth / "true_ancestors.fasta")
    proteins = read_fasta(results / folder / "ancestors.fasta")
    coding = read_fasta(results / folder / "ancestors.cds.fasta")
    nodes: list[dict[str, Any]] = []
    protein_ids: list[float] = []
    codon_ids: list[float] = []
    for true_node, inferred_node in sorted(names.items()):
        truth_codons = true_ancestors[true_node]
        truth_protein = codons.translate_row(truth_codons, code)
        protein_ids.append(identity(proteins[inferred_node], truth_protein))
        codon_ids.append(codon_identity(coding[inferred_node], truth_codons))
        nodes.append({
            "true_node": true_node, "inferred_node": inferred_node,
            "identity": protein_ids[-1], "codon_identity": codon_ids[-1],
        })  # fmt: skip
    if not nodes:
        raise ValueError("no inferred node matches a true one")
    return {
        "nodes_true": len(true_ancestors),
        "nodes_compared": len(nodes),
        "mean_identity": fmean(protein_ids),
        "min_identity": min(protein_ids),
        "mean_codon_identity": fmean(codon_ids),
        "min_codon_identity": min(codon_ids),
        "nodes": nodes,
    }


def check_codons(metrics: dict[str, Any], expected: dict[str, Any]) -> list[str]:
    """Every threshold in `expected` (keys as in expected_metrics.json) that is missed."""
    rules: dict[str, tuple[float, Callable[[float, float], bool]]] = {
        "min_nodes_compared": (metrics["nodes_compared"], operator.ge),
        "min_mean_identity": (metrics["mean_identity"], operator.ge),
        "min_identity": (metrics["min_identity"], operator.ge),
        "min_mean_codon_identity": (metrics["mean_codon_identity"], operator.ge),
        "min_codon_identity": (metrics["min_codon_identity"], operator.ge),
    }
    failures = []
    for key, threshold in expected.items():
        if key.startswith("_"):
            continue
        if key not in rules:
            failures.append(f"unknown threshold {key!r}")
        elif not rules[key][1](rules[key][0], threshold):
            failures.append(f"{key}: got {rules[key][0]:.4g}, threshold {threshold}")
    return failures


def summary(label: str, metrics: dict[str, Any]) -> str:
    return (
        f"  {label:<34} {metrics['nodes_compared']}/{metrics['nodes_true']} nodes · "
        f"protein identity {metrics['mean_identity']:.3f} (min {metrics['min_identity']:.3f}) · "
        f"codon identity {metrics['mean_codon_identity']:.3f} "
        f"(min {metrics['min_codon_identity']:.3f})"
    )


def main() -> None:
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    folder = sys.argv[3] if len(sys.argv) == 4 else ""
    metrics = score_codons(Path(sys.argv[1]), Path(sys.argv[2]), folder=folder)
    print(json.dumps({k: v for k, v in metrics.items() if k != "nodes"}, indent=2))


if __name__ == "__main__":
    main()
