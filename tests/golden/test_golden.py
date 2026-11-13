"""Run the real pipeline on each golden replicate and score it against the truth.

Slow and needs Docker, so excluded from `make test`; run with `make golden`.

Each run also feeds the true alignment and true tree back in as provided inputs
("oracles"), giving four combinations. Comparing them attributes the error: what
remains with both oracles is the ancestral reconstruction's own error; the gaps
to the full pipeline are what alignment and tree inference add.

Every run also writes both gap-masked ancestor versions (parsimony, ml); for scenarios
whose true ancestors have gaps, each version is scored and compared.

Thresholds live in each scenario's expected_metrics.json, set by the project owner,
and apply to the full pipeline (mafft/iqtree/modelfinder), scoring ancestors.fasta or,
if the file names one ("gap_method": "ml"), that masked version. Without that file the
test records metrics and skips, so a first run can be reviewed.
"""

import json
import statistics
import tomllib
from pathlib import Path
from typing import Any

import pytest

from phylobot import codons
from phylobot.evaluation import codons as codon_scoring
from phylobot.evaluation.score import attribution, check, score, summary
from phylobot.gaps import METHODS as GAP_METHODS

from .pipeline import run_pipeline

GOLDEN = Path(__file__).resolve().parent
ROOT = GOLDEN.parent.parent


def is_codon(spec: dict[str, Any]) -> bool:
    return str(spec["model"].get("seqtype", "")).startswith("CODON")


def is_family(spec: dict[str, Any]) -> bool:
    """A gene family grown on a species tree (reconcile/): scored by test_reconcile.py."""
    return "family" in spec


CODON_SKIP = pytest.mark.skip(
    reason="codon scenarios are scored by test_codon_pipeline_recovers_truth"
)


def scenario_replicates() -> list[Any]:
    """Every replicate; all but each scenario's first seed are in the full tier only."""
    params = []
    for spec_path in sorted(GOLDEN.glob("*/*/spec.toml")):
        scenario = spec_path.parent
        spec = tomllib.loads(spec_path.read_text())
        for i, seed in enumerate(spec["simulation"]["seeds"]):
            replicate = scenario / "replicates" / str(seed)
            if not replicate.is_dir():
                continue
            marks = [pytest.mark.golden_full] if i > 0 else []
            if is_codon(spec) or is_family(spec):
                continue  # scored by test_codon_pipeline_recovers_truth, test_reconcile.py
            name = str(replicate.relative_to(GOLDEN))
            params.append(pytest.param(replicate, marks=marks, id=name))
    return params


def is_quick(replicate: Path) -> bool:
    spec = tomllib.loads((replicate.parent.parent / "spec.toml").read_text())
    return str(spec["simulation"]["seeds"][0]) == replicate.name


PIPELINE = "mafft/iqtree/modelfinder"
ORACLES = {
    PIPELINE: "full pipeline",
    "provided1/iqtree/modelfinder": "true alignment",
    "mafft/provided1/modelfinder": "true tree",
    "provided1/provided1/modelfinder": "true alignment + tree",
}


@pytest.mark.golden
@pytest.mark.parametrize("replicate", scenario_replicates())
def test_pipeline_recovers_truth(replicate: Path) -> None:
    scenario = replicate.parent.parent
    spec = tomllib.loads((scenario / "spec.toml").read_text())
    out = ROOT / "results" / "golden" / replicate.relative_to(GOLDEN)
    params = {
        "alignment": [str((replicate / "true_alignment.fasta").relative_to(ROOT))],
        "alignment_labels": ["True alignment"],
        "tree": [str((replicate / "true_tree.nwk").relative_to(ROOT))],
        "tree_labels": ["True tree"],
        "gap_methods": list(GAP_METHODS),
        "outgroup": [spec["tree"]["outgroup"]],
    }
    run = run_pipeline(out, replicate / "input.fasta", params)
    outgroup = spec["tree"]["outgroup"]
    by_role = {role: score(replicate, run / key, outgroup) for key, role in ORACLES.items()}
    (out / "metrics.json").write_text(json.dumps(by_role, indent=2) + "\n")
    print(f"\n{replicate.relative_to(GOLDEN)}\n{summary(by_role['full pipeline'])}")
    print(attribution(by_role))

    if by_role["full pipeline"]["truth_has_gaps"]:
        gap_table = {
            method: {
                role: score(replicate, run / key, outgroup, ancestors_file=f)
                for key, role in ORACLES.items()
            }
            for method, f in ancestor_versions().items()
        }
        (out / "gap_metrics.json").write_text(json.dumps(gap_table, indent=2) + "\n")
        print(gap_comparison(gap_table))

    expected_path = scenario / "expected_metrics.json"
    if not expected_path.exists():
        pytest.skip(f"no thresholds yet; metrics written to {out / 'metrics.json'}")
    if not is_quick(replicate):
        return  # full-tier replicates are judged on their average (test_replicate_averages)
    expected = json.loads(expected_path.read_text())
    metrics = by_role["full pipeline"]
    if method := expected.get("gap_method"):
        version = ancestor_versions()[method]
        metrics = score(replicate, run / PIPELINE, outgroup, ancestors_file=version)
    failures = check(metrics, {k: v for k, v in expected.items() if k != "gap_method"})
    assert not failures, "\n".join(failures)


def ancestor_versions() -> dict[str, str]:
    return {"none": "ancestors.fasta"} | {m: f"ancestors.gaps-{m}.fasta" for m in GAP_METHODS}


def gap_comparison(table: dict[str, dict[str, dict[str, object]]]) -> str:
    """Gap placement methods side by side, for the full pipeline and with oracles."""
    lines = ["  gap placement (gap recall / precision; mean length error; residue identity):"]
    for role in ORACLES.values():
        cells = []
        for method, by_role in table.items():
            m = by_role[role]
            recall = m["gap_recall"] or 0.0
            precision = m["gap_precision"]
            cells.append(
                f"{method} {recall:.2f}/{'-' if precision is None else f'{precision:.2f}'} "
                f"{m['mean_length_error']:+.0f} {m['mean_identity']:.3f}"
            )
        lines.append(f"    {role:<22} " + " | ".join(cells))
    return "\n".join(lines)


MULTI = [
    spec.parent
    for spec in sorted(GOLDEN.glob("*/*/spec.toml"))
    if len(tomllib.loads(spec.read_text())["simulation"]["seeds"]) > 1
    and not is_codon(tomllib.loads(spec.read_text()))
    and not is_family(tomllib.loads(spec.read_text()))
]


@pytest.mark.golden
@pytest.mark.golden_full
@pytest.mark.parametrize("scenario", MULTI, ids=[str(s.relative_to(GOLDEN)) for s in MULTI])
def test_replicate_averages(scenario: Path) -> None:
    """Across all replicates of a scenario: each gap method's mean and spread, how often
    each beats the others, and the owner's thresholds applied to the averages.
    Runs after (and reads the metrics written by) test_pipeline_recovers_truth."""
    seeds = tomllib.loads((scenario / "spec.toml").read_text())["simulation"]["seeds"]
    out = ROOT / "results" / "golden" / scenario.relative_to(GOLDEN) / "replicates"
    missing = [s for s in seeds if not (out / str(s) / "gap_metrics.json").exists()]
    if missing:
        pytest.fail(f"no results for seeds {missing}; run the full tier (make golden-full)")
    runs = [json.loads((out / str(s) / "gap_metrics.json").read_text()) for s in seeds]
    print(f"\n{scenario.relative_to(GOLDEN)}: {len(runs)} replicates")
    print(replicate_report(runs))
    expected_path = scenario / "expected_metrics.json"
    if not expected_path.exists():
        pytest.skip("no thresholds yet")
    expected = json.loads(expected_path.read_text())
    method = expected.get("gap_method", "none")
    rows = [run[method]["full pipeline"] for run in runs]
    means = {key: statistics.fmean(float(r[key] or 0.0) for r in rows) for key in AVERAGED}
    means |= {"rf_distance": statistics.fmean(r["rf_distance"] for r in rows)}
    means |= {"clades_recovered": statistics.fmean(r["clades_recovered"] for r in rows)}
    means |= {"clades_true": rows[0]["clades_true"], "min_identity": means["mean_identity"]}
    failures = check(means, {k: v for k, v in expected.items() if k != "gap_method"})
    assert not failures, "averages across replicates: " + "; ".join(failures)


AVERAGED = ("gap_recall", "gap_precision", "gap_accuracy", "mean_abs_length_error",
            "mean_length_error", "mean_identity", "residue_accuracy")  # fmt: skip


def stat(rows: list[dict[str, Any]], key: str, fmt: str = ".3f") -> str:
    """'mean ± sd' of one metric over replicates (missing values count as 0)."""
    values = [float(r[key] or 0.0) for r in rows]
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    return f"{statistics.fmean(values):{fmt}} ± {sd:{fmt}}"


def replicate_report(runs: list[dict[str, Any]]) -> str:
    """Mean ± sd per gap method (full pipeline and with true alignment), and paired wins."""
    methods = list(runs[0])
    lines = []
    for role in ("full pipeline", "true alignment + tree"):
        lines.append(f"  {role}: mean ± sd over replicates")
        lines.append(f"    {'method':<10} {'gap recall':>13} {'precision':>13} "
                     f"{'length error':>14} {'identity':>13}")  # fmt: skip
        for method in methods:
            rows = [run[method][role] for run in runs]
            recall, precision = stat(rows, "gap_recall"), stat(rows, "gap_precision")
            length, identity = stat(rows, "mean_length_error", "+.1f"), stat(rows, "mean_identity")
            lines.append(
                f"    {method:<10} {recall:>13} {precision:>13} {length:>14} {identity:>13}"
            )
    wins = {m: 0 for m in methods}
    for run in runs:
        best = max(methods, key=lambda m: run[m]["full pipeline"]["mean_identity"])
        wins[best] += 1
    lines.append("  best identity (full pipeline), replicates won: "
                 + ", ".join(f"{m} {n}" for m, n in wins.items()))  # fmt: skip
    return "\n".join(lines)


# Codon scenarios (the software check): codon-easy only. The others, codon-selection
# (codon vs amino-acid models, a research question) and codon-indels (gaps on codon
# columns), wait.
VALIDATED_CODON = ("lysin/codon-easy",)
CODON_MODEL = "GY+F3X4+G4"  # the family the scenario was simulated under
CODON_ROLES = {
    f"mafft/provided1/{CODON_MODEL}": "MAFFT alignment, true tree",
    f"provided1/provided1/{CODON_MODEL}": "true alignment + tree",
}


def codon_replicates() -> list[Any]:
    params = []
    for spec_path in sorted(GOLDEN.glob("*/*/spec.toml")):
        scenario = spec_path.parent
        spec = tomllib.loads(spec_path.read_text())
        if not is_codon(spec):
            continue
        name = str(scenario.relative_to(GOLDEN))
        replicate = scenario / "replicates" / str(spec["simulation"]["seeds"][0])
        marks = [] if name in VALIDATED_CODON else [CODON_SKIP]
        if replicate.is_dir():
            params.append(pytest.param(replicate, marks=marks, id=f"{name}/{replicate.name}"))
    return params


@pytest.mark.golden
@pytest.mark.parametrize("replicate", codon_replicates())
def test_codon_pipeline_recovers_truth(replicate: Path) -> None:
    """The codon path end to end on codon-easy: codon checks and alignment
    (back-translation), IQ-TREE's codon fit and ancestors, and codeml M0 with its node
    mapping, against the true ancestors. The true tree is given (codon tree searches
    are slow, and IQ-TREE's own); MAFFT's alignment and the true one are compared."""
    scenario = replicate.parent.parent
    out = ROOT / "results" / "golden" / replicate.relative_to(GOLDEN)
    params = {
        "seqtype": "codon",
        "aligners": ["mafft"],
        "tree_builders": ["none"],
        "modelfinder": False,
        "codon_models": [CODON_MODEL],
        "alignment": [str((replicate / "true_alignment.fasta").relative_to(ROOT))],
        "alignment_labels": ["True alignment"],
        "tree": [str((replicate / "true_tree.nwk").relative_to(ROOT))],
        "tree_labels": ["True tree"],
        "gap_methods": [],
        "codeml_models": ["M0"],
        "codeml_trees": "all",
    }
    run = run_pipeline(out, replicate / "input.fasta", params)
    metrics = {}
    for key, role in CODON_ROLES.items():
        for program, folder in (("iqtree", ""), ("codeml_M0", "codeml/M0")):
            found = codon_scoring.score_codons(replicate, run / key, codons.STANDARD, folder)
            metrics[f"{role} · {program}"] = found
            print(codon_scoring.summary(f"{role} · {program}", found))
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    expected_path = scenario / "expected_metrics.json"
    if not expected_path.exists():
        pytest.skip(f"no thresholds yet; metrics written to {out / 'metrics.json'}")
    expected = json.loads(expected_path.read_text())
    failures = [
        f"{label}: {failure}"
        for label, thresholds in expected.items()
        if not label.startswith("_")
        for failure in codon_scoring.check_codons(metrics[label], thresholds)
    ]
    assert not failures, "\n".join(failures)
