"""Write summary.json and model_comparison.tsv for a job's combinations.

    phylobot-summary SOURCES.tsv FITS.tsv PRIMARY ALIGNMENTS STOP_AFTER

FITS.tsv has one row per combination: aligner, tree, model (keys), then the fit.tsv
fields of that combination (export.FIT_FIELDS). Results live at <outdir>/<key>/, where
key = aligner/tree/model. The web app shows the primary combination first; if PRIMARY
did not run, the first combination in sorted order is used. ALIGNMENTS lists the
job's alignment keys (comma-separated), each published at <outdir>/<key>/alignment.fasta;
STOP_AFTER is the last stage run (alignments, trees or ancestors). A job that stops
after alignments has no combinations.

model_comparison.tsv ranks combinations by BIC *within each alignment*: likelihoods and
information criteria are only comparable on the same data (BIC also depends on the
number of sites). Weights are computed over distinct fits, so when ModelFinder picked a
model that also ran as a fixed model, both rows show the same weight.
"""

import json
import math
import sys
from pathlib import Path
from typing import Any

from phylobot import catalog
from phylobot.pipeline.export import FIT_FIELDS

# Display labels for tool keys; provided inputs take theirs from sources.tsv.
TOOL_LABELS = {
    "mafft": "MAFFT",
    "muscle": "MUSCLE",
    "clustalo": "Clustal Omega",
    "prank": "PRANK",
    "prank_codon": "PRANK (codon model)",
    "iqtree": "IQ-TREE",
    "raxml": "RAxML-NG",
    "generax": "GeneRax (reconciled)",
}
COMPARISON_FIELDS = (
    "aligner", "tree", "model_key", "model", "log_likelihood", "free_parameters",
    "aic", "aicc", "bic", "delta_bic", "w_bic", "w_aic", "w_aicc", "rank_bic", "data",
)  # fmt: skip


def data_type(model_key: str) -> str:
    """'codon' for a codon model's combination, else 'amino'. Likelihoods of the codon
    alignment and of its translation are not comparable, so each type is ranked apart."""
    codon = model_key == catalog.CODON_MODELFINDER or catalog.is_codon_model(model_key)
    return "codon" if codon else "amino"


def tree_label(key: str, sources: dict[tuple[str, str], str], best: str = "") -> str:
    """A tree source's label. A codon model fitted on an amino-acid tree has key
    <source>_aa (that source's amino-acid tree of lowest BIC: `best` names its model key)
    or <source>_aa-<model key> (the tree of that amino-acid model)."""
    if ("tree", key) in sources:
        return sources[("tree", key)]
    source, aa, model = key.partition("_aa")
    if aa and source in TOOL_LABELS and (not model or model.startswith("-")):
        if model:
            which = catalog.label(model[1:])
        else:
            which = f"best: {catalog.label(best)}" if best else "best BIC"
        return f"{TOOL_LABELS[source]} amino-acid tree ({which})"
    return TOOL_LABELS.get(key, key)


def best_amino(fits: list[dict[str, str]], aligner: str, tree: str) -> str:
    """The amino-acid model key whose tree a <tree>_aa codon combination was fitted on:
    lowest BIC among that alignment's and tree source's amino-acid combinations, ties to
    the first model key (as main.nf chooses)."""
    group = [
        (float(f["bic"]), f["model_key"])
        for f in fits
        if f["aligner"] == aligner and f["tree"] == tree and data_type(f["model_key"]) == "amino"
    ]
    return min(group)[1] if group else ""


def read_sources(path: Path) -> dict[tuple[str, str], str]:
    """{(kind, key): label} from sources.tsv."""
    rows = [line.split("\t") for line in path.read_text().splitlines()[1:]]
    return {(kind, key): label for key, kind, label, _ in rows}


def read_fits(path: Path) -> list[dict[str, str]]:
    fields = ("aligner", "tree", "model_key", *FIT_FIELDS)
    rows = [line.split("\t") for line in path.read_text().splitlines() if line.strip()]
    return [dict(zip(fields, row, strict=True)) for row in rows]


def weights(scores: dict[Any, float]) -> dict[Any, float]:
    """Akaike-style weights: exp(-delta/2), normalised."""
    best = min(scores.values())
    raw = {key: math.exp(-(score - best) / 2) for key, score in scores.items()}
    total = sum(raw.values())
    return {key: value / total for key, value in raw.items()}


def compare(fits: list[dict[str, str]]) -> list[dict[str, str]]:
    """Rank combinations by BIC within each alignment and data type (amino-acid models,
    then codon models), with criterion weights."""
    rows: list[dict[str, str]] = []
    groups = sorted({(fit["aligner"], data_type(fit["model_key"])) for fit in fits})
    for aligner, data in groups:
        group = [
            fit for fit in fits
            if fit["aligner"] == aligner and data_type(fit["model_key"]) == data
        ]  # fmt: skip
        distinct = {(fit["tree"], fit["model"]): fit for fit in group}
        scores = {
            criterion: weights({k: float(f[criterion]) for k, f in distinct.items()})
            for criterion in ("bic", "aic", "aicc")
        }
        best_bic = min(float(fit["bic"]) for fit in group)
        ranking = sorted(distinct, key=lambda k: float(distinct[k]["bic"]))
        for fit in sorted(group, key=lambda f: (float(f["bic"]), f["model_key"])):
            key = (fit["tree"], fit["model"])
            rows.append(
                {
                    **{f: fit[f] for f in COMPARISON_FIELDS if f in fit},
                    "delta_bic": f"{float(fit['bic']) - best_bic:.4f}",
                    "w_bic": f"{scores['bic'][key]:.4g}",
                    "w_aic": f"{scores['aic'][key]:.4g}",
                    "w_aicc": f"{scores['aicc'][key]:.4g}",
                    "rank_bic": str(ranking.index(key) + 1),
                    "data": data,
                }
            )
    return rows


def summarize(
    sources_path: Path,
    fits_path: Path,
    primary: str,
    alignment_keys: list[str] | None = None,
    stop_after: str = "ancestors",
) -> None:
    sources = read_sources(sources_path)
    fits = {f"{f['aligner']}/{f['tree']}/{f['model_key']}": f for f in read_fits(fits_path)}

    def labels(key: str) -> dict[str, str]:
        aligner, tree, model = key.split("/")
        best = best_amino(list(fits.values()), aligner, tree[:-3]) if tree.endswith("_aa") else ""
        return {
            "aligner_label": sources.get(("alignment", aligner), TOOL_LABELS.get(aligner, aligner)),
            "tree_label": tree_label(tree, sources, best),
            "model_label": catalog.label(model),
        }

    def label(key: str) -> str:
        return " × ".join(labels(key).values())

    comparison = compare(list(fits.values()))
    by_key = {f"{r['aligner']}/{r['tree']}/{r['model_key']}": r for r in comparison}
    keys = sorted(fits)
    chosen = primary if primary in keys else (keys[0] if keys else "")
    aligned = sorted(alignment_keys or {k.split("/")[0] for k in keys})
    summary = {
        "stop_after": stop_after,
        "alignments": [
            {
                "key": key,
                "label": sources.get(("alignment", key), TOOL_LABELS.get(key, key)),
                "path": f"{key}/alignment.fasta",
            }
            for key in aligned
        ],
        "primary": chosen,
        "combinations": [
            {
                "key": key,
                "aligner": key.split("/")[0],
                "tree": key.split("/")[1],
                "model": key.split("/")[2],
                "label": label(key),
                **labels(key),
                "path": key,
                "fit": {f: fits[key][f] for f in FIT_FIELDS},
                "rank_bic": int(by_key[key]["rank_bic"]),
                "data": data_type(key.split("/")[2]),
                "delta_bic": float(by_key[key]["delta_bic"]),
            }
            for key in keys
        ],
    }
    with open("summary.json", "w") as out:
        json.dump(summary, out, indent=2, ensure_ascii=False)
        out.write("\n")
    lines = ["\t".join(COMPARISON_FIELDS)]
    lines += ["\t".join(row[f] for f in COMPARISON_FIELDS) for row in comparison]
    Path("model_comparison.tsv").write_text("\n".join(lines) + "\n")


def main() -> None:
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    keys = [k for k in sys.argv[4].split(",") if k]
    summarize(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], keys, sys.argv[5])


if __name__ == "__main__":
    main()
