"""Score one pipeline result against a simulated truth.

uv run python -m phylobot.evaluation.score TRUTH_DIR RESULTS_DIR --outgroup NAME
"""

import argparse
import operator
from collections.abc import Callable
from dataclasses import asdict
from pathlib import Path
from typing import Any

from phylobot.evaluation.gaps import GapCalls, aligned_identity, column_map, gap_calls
from phylobot.evaluation.sequences import calibration, identity
from phylobot.seqio import read_fasta
from phylobot.trees import clades, load_tree, robinson_foulds


def read_tsv(path: Path) -> list[dict[str, str]]:
    header, *rows = path.read_text().splitlines()
    keys = header.split("\t")
    return [dict(zip(keys, row.split("\t"), strict=True)) for row in rows]


def score(
    truth: Path, results: Path, outgroup: str, ancestors_file: str = "ancestors.fasta"
) -> dict[str, Any]:
    """Compare results/{tree.nwk, <ancestors_file>, alignment.fasta, ...} with
    truth/{true_tree.nwk, true_ancestors.fasta, true_alignment.fasta}.

    Identity is position-by-position when the true ancestors have no gaps (no-indel
    scenarios; unchanged from before indel scoring), otherwise from a global pairwise
    alignment: with indels, equal lengths do not mean the residues line up. Gap calls,
    length errors and residue calibration use a column map between the inferred and
    true alignments (see evaluation.gaps).
    """
    true_newick = (truth / "true_tree.nwk").read_text()
    inferred_newick = (results / "tree.nwk").read_text()
    true_clades = clades(load_tree(true_newick), outgroup)
    inferred_clades = clades(load_tree(inferred_newick), outgroup)
    true_ancestors = read_fasta(truth / "true_ancestors.fasta")
    inferred_ancestors = read_fasta(results / ancestors_file)
    mapping = column_map(
        read_fasta(results / "alignment.fasta"), read_fasta(truth / "true_alignment.fasta")
    )
    truth_has_gaps = any("-" in row for row in true_ancestors.values())
    supports = {row["node"]: row["ufboot"] for row in read_tsv(results / "supports.tsv")}

    site_pp: dict[str, list[tuple[str, float]]] = {}
    for row in read_tsv(results / "ancestors_pp.tsv"):
        site_pp.setdefault(row["node"], []).append((row["state"], float(row["pp"])))

    ancestors: list[dict[str, Any]] = []
    calls: list[tuple[float, bool]] = []
    all_gap_calls = GapCalls()
    length_errors: list[int] = []
    for clade, true_node in sorted(true_clades.items(), key=lambda item: len(item[0])):
        entry: dict[str, Any] = {"clade": sorted(clade), "true_node": true_node}
        inferred_node = inferred_clades.get(clade)
        if inferred_node is None:
            entry["recovered"] = False
            ancestors.append(entry)
            continue
        truth_row = true_ancestors[true_node]
        inferred_row = inferred_ancestors[inferred_node]
        inferred_residues = inferred_row.replace("-", "")
        truth_residues = truth_row.replace("-", "")
        positional = not truth_has_gaps and len(inferred_residues) == len(truth_residues)
        node_calls = gap_calls(inferred_row, truth_row, mapping)
        all_gap_calls += node_calls
        length_errors.append(len(inferred_residues) - len(truth_residues))
        entry |= {
            "recovered": True,
            "inferred_node": inferred_node,
            "support": supports.get(inferred_node) or None,
            "identity": (
                identity(inferred_row, truth_row)
                if positional
                else aligned_identity(inferred_row, truth_row)
            ),
            "length": len(inferred_residues),
            "true_length": len(truth_residues),
            "gap_recall": node_calls.gap_recall,
            "gap_precision": node_calls.gap_precision,
        }
        # Residue calls at inferred columns that map to a true residue (a true gap there
        # is a gap-call error, scored above, not a residue error).
        for column, (state, pp) in enumerate(site_pp[inferred_node]):
            true_column = mapping.get(column)
            if true_column is None or truth_row[true_column] == "-":
                continue
            if inferred_row[column] == "-":
                continue
            calls.append((pp, state.upper() == truth_row[true_column].upper()))
        ancestors.append(entry)

    identities = [a["identity"] for a in ancestors if a["recovered"]]
    return {
        "rf_distance": robinson_foulds(true_newick, inferred_newick),
        "clades_true": len(true_clades),
        "clades_recovered": len(identities),
        "mean_identity": sum(identities) / len(identities) if identities else 0.0,
        "min_identity": min(identities, default=0.0),
        "sites_scored": len(calls),
        # correct residue calls where the true ancestor has a residue: substitution
        # accuracy, separate from gap/length errors (which identity also penalises)
        "residue_accuracy": sum(ok for _, ok in calls) / len(calls) if calls else 0.0,
        "calibration": [asdict(b) for b in calibration(calls)],
        "truth_has_gaps": truth_has_gaps,
        "mean_length_error": sum(length_errors) / len(length_errors) if length_errors else 0.0,
        "mean_abs_length_error": (
            sum(abs(e) for e in length_errors) / len(length_errors) if length_errors else 0.0
        ),
        "gap_columns_scored": all_gap_calls.columns,
        "gap_accuracy": all_gap_calls.accuracy,
        "gap_recall": all_gap_calls.gap_recall,
        "gap_precision": all_gap_calls.gap_precision,
        "lost_residues": all_gap_calls.lost_residues,
        "ancestors": ancestors,
    }


def check(metrics: dict[str, Any], expected: dict[str, Any]) -> list[str]:
    """Return a description of every threshold in `expected` that `metrics` misses."""
    failures = []
    recovered = metrics["clades_recovered"] / metrics["clades_true"]
    rules: dict[str, tuple[float, Callable[[float, float], bool]]] = {
        "max_rf_distance": (metrics["rf_distance"], operator.le),
        "min_clades_recovered_fraction": (recovered, operator.ge),
        "min_mean_identity": (metrics["mean_identity"], operator.ge),
        "min_identity": (metrics["min_identity"], operator.ge),
        "max_mean_abs_length_error": (metrics.get("mean_abs_length_error", 0.0), operator.le),
        "min_gap_recall": (metrics.get("gap_recall") or 0.0, operator.ge),
        "min_gap_precision": (metrics.get("gap_precision") or 0.0, operator.ge),
        "min_gap_accuracy": (metrics.get("gap_accuracy") or 0.0, operator.ge),
    }
    for key, threshold in expected.items():
        if key.startswith("_"):
            continue  # comments / provenance
        if key not in rules:
            failures.append(f"unknown threshold {key!r}")
            continue
        value, ok = rules[key]
        if not ok(value, threshold):
            failures.append(f"{key}: got {value:.4g}, threshold {threshold}")
    return failures


def summary(metrics: dict[str, Any]) -> str:
    lines = [
        f"RF distance {metrics['rf_distance']}; "
        f"clades recovered {metrics['clades_recovered']}/{metrics['clades_true']}; "
        f"identity mean {metrics['mean_identity']:.3f}, min {metrics['min_identity']:.3f}",
    ]
    for a in metrics["ancestors"]:
        name = f"true {a['true_node']:>3} ({len(a['clade'])} taxa)"
        if a["recovered"]:
            lines.append(
                f"  {name} -> {a['inferred_node']:<6} "
                f"support {a['support'] or '-':>4}  identity {a['identity']:.3f}"
            )
        else:
            lines.append(f"  {name} -> not recovered")
    if metrics.get("truth_has_gaps"):

        def pct(value: float | None) -> str:
            return "n/a" if value is None else f"{value:.3f}"

        lines.append(
            f"  gaps: length error mean {metrics['mean_length_error']:+.1f} residues "
            f"(abs {metrics['mean_abs_length_error']:.1f}); column accuracy "
            f"{pct(metrics['gap_accuracy'])}, gap recall {pct(metrics['gap_recall'])}, "
            f"gap precision {pct(metrics['gap_precision'])} over "
            f"{metrics['gap_columns_scored']} columns; {metrics['lost_residues']} true "
            "residues in columns the alignment lacks"
        )
    lines.append(f"  calibration over {metrics['sites_scored']} site calls:")
    for b in metrics["calibration"]:
        if b["count"]:
            lines.append(
                f"    pp {b['low']:.2f}-{b['high']:.2f}: n={b['count']:<5} "
                f"mean pp {b['mean_pp']:.3f}  accuracy {b['accuracy']:.3f}"
            )
    return "\n".join(lines)


def attribution(by_role: dict[str, dict[str, Any]]) -> str:
    """Table comparing the full pipeline with oracle runs (true alignment and/or tree).

    Expects roles "full pipeline", "true alignment", "true tree", "true alignment + tree".
    Identity lost with both oracles is inherent uncertainty plus the reconstruction
    method's error; the extra loss in the full pipeline comes from alignment and tree
    inference.
    """
    gaps = any(m.get("truth_has_gaps") for m in by_role.values())
    heading = "mean ancestor identity; residue accuracy; RF distance; clades recovered"
    if gaps:
        heading += "; mean length error"
    lines = [f"  error attribution ({heading}):"]
    for role, m in by_role.items():
        line = (
            f"    {role:<22} identity {m['mean_identity']:.3f}  "
            f"residues {m.get('residue_accuracy', 0.0):.3f}  RF {m['rf_distance']}  "
            f"clades {m['clades_recovered']}/{m['clades_true']}"
        )
        if gaps:
            line += f"  length {m.get('mean_length_error', 0.0):+.1f}"
        lines.append(line)
    # With indels, identity also penalises ancestral length, so the attribution uses
    # residue accuracy and reports gap placement separately.
    key = "residue_accuracy" if gaps else "mean_identity"
    measure = "residue accuracy" if gaps else "identity"
    both = by_role["true alignment + tree"][key]
    full = by_role["full pipeline"][key]
    lines.append(
        f"    -> by {measure}: with true alignment + tree, ancestors still miss "
        f"{1 - both:.3f} (inherent uncertainty + reconstruction method); "
        f"alignment + tree inference add {both - full:.3f} "
        f"(alignment alone {by_role['true alignment'][key] - full:+.3f}, "
        f"tree alone {by_role['true tree'][key] - full:+.3f} when made true)"
    )
    if gaps:
        lengths = by_role["true alignment + tree"].get("mean_length_error", 0.0)
        lines.append(
            f"    -> gap placement: even with true alignment + tree, ancestors are "
            f"{lengths:+.1f} residues off in length (gap recall "
            f"{by_role['true alignment + tree'].get('gap_recall') or 0.0:.3f})"
        )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("truth", type=Path)
    parser.add_argument("results", type=Path)
    parser.add_argument("--outgroup", required=True)
    args = parser.parse_args()
    metrics = score(args.truth, args.results, args.outgroup)
    print(summary(metrics))


if __name__ == "__main__":
    main()
