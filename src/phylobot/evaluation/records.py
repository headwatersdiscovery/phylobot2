"""Per-site records for calibration.

One row per node x alignment column: the ML state, its PP, the full posterior, and the
covariates `phylobot.calibration` applies at view time, computed by the same functions
from the same published files (alignment.fasta, tree.nwk, site_rates.tsv). A formula
fitted on these rows is therefore the formula a user's job gets. Truth, where there is
any, fills `true` and `correct`:

- **asr**: every ancestor of a combination (ancestors_posteriors.tsv). On simulated
  data `truth_records` joins the true ancestors the way the golden scorer does: nodes by
  clade on outgroup-rooted trees, columns through extant residues (`column_map`, or
  with `columns="clade"` a map per node from its own descendants: a clade the aligner
  shifted against the rest is then compared with its own homologous residues, so
  placement error and residue error are not one number; PPs speak only to the latter).
- **esr-e0 / esr-e1**: held-out tips from `phylobot.esr` (esr_sites.tsv), truth the
  tip's own residue. A tip's covariates are its own branch and the path to the nearest
  *other* sequence (`calibration.tip_place`), and the column statistics leave the tip
  out: while it is reconstructed, its sequence is unknown.

`present` is 1 where the chosen gap version of the ancestor has a residue, 0 where it
has a gap, empty when no version was asked for, so residue calls can be scored alone.
"""

import csv
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Literal

from phylobot import calibration
from phylobot.alignment_map import majority_columns
from phylobot.evaluation.gaps import column_map
from phylobot.seqio import read_fasta
from phylobot.trees import clades, load_tree

AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"
FIELDS = [
    "source", "node", "site", "ml", "pp_ml", "true", "correct", "present",
    *calibration.KEYS, *AMINO_ACIDS,
]  # fmt: skip


class Combination:
    """A combination's published files, read once."""

    def __init__(self, folder: Path):
        self.folder = folder
        self.alignment = read_fasta(folder / "alignment.fasta")
        self.tree = (folder / "tree.nwk").read_text()
        rates = folder / "site_rates.tsv"
        self.rates = rates.read_text() if rates.exists() else None

    def posteriors(self) -> dict[str, dict[int, dict[str, float]]]:
        """ancestors_posteriors.tsv as {node: {site: {state: p}}}."""
        out: dict[str, dict[int, dict[str, float]]] = {}
        with (self.folder / "ancestors_posteriors.tsv").open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                probs = {a: float(row[a]) for a in AMINO_ACIDS}
                out.setdefault(row["node"], {})[int(row["site"])] = probs
        return out

    def gap_version(self, method: str) -> dict[str, str]:
        return read_fasta(self.folder / f"ancestors.gaps-{method}.fasta")


def _row(
    source: str,
    node: str,
    site: int,
    probs: dict[str, float],
    covariates: dict[str, float],
    true: str = "",
    present: str = "",
) -> dict[str, object]:
    ml = max(AMINO_ACIDS, key=lambda a: probs[a])
    return {
        "source": source, "node": node, "site": site, "ml": ml, "pp_ml": probs[ml],
        "true": true, "correct": "" if not true else int(true == ml), "present": present,
        **covariates, **probs,
    }  # fmt: skip


def ancestor_records(folder: Path, gap_method: str | None = None) -> Iterator[dict[str, object]]:
    """Every ancestor x column of one combination, covariates as the viewer computes them."""
    combination = Combination(folder)
    extants = list(combination.alignment.values())
    masked = combination.gap_version(gap_method) if gap_method else {}
    for node, sites in combination.posteriors().items():
        rows = calibration.covariates(extants, combination.tree, node, combination.rates)
        if rows is None:
            raise ValueError(f"{node} is not a node of {folder}/tree.nwk")
        version = masked.get(node)
        if masked and version is None:
            raise ValueError(f"{node} has no gap-{gap_method} version")
        for site in sorted(sites):
            present = "" if version is None else str(int(version[site - 1] != "-"))
            yield _row("asr", node, site, sites[site], rows[site - 1], present=present)


def truth_records(
    folder: Path, truth: Path, outgroup: str, gap_method: str | None = None,
    columns: Literal["global", "clade"] = "global",
) -> Iterator[dict[str, object]]:  # fmt: skip
    """`ancestor_records` with the simulation's truth (truth/{true_tree.nwk,
    true_ancestors.fasta, true_alignment.fasta}): only the calls whose node is recovered
    and whose column maps to a true residue. A true gap there is a gap error, scored by
    `evaluation.gaps`, not a residue call. With `columns="clade"` a node's descendants
    alone choose its true columns; a column none of them occupies keeps the whole
    alignment's."""
    true_ancestors = read_fasta(truth / "true_ancestors.fasta")
    true_clades = clades(load_tree((truth / "true_tree.nwk").read_text()), outgroup)
    inferred = clades(load_tree((folder / "tree.nwk").read_text()), outgroup)
    match = {inferred[c]: node for c, node in true_clades.items() if c in inferred}
    aligned = read_fasta(folder / "alignment.fasta")
    true_aligned = read_fasta(truth / "true_alignment.fasta")
    mapping = column_map(aligned, true_aligned)
    per_node: dict[str, dict[int, int]] = {}
    if columns == "clade":
        for clade, node in inferred.items():
            if clade in true_clades:
                per_node[node] = majority_columns({n: aligned[n] for n in clade},
                                                  {n: true_aligned[n] for n in clade})  # fmt: skip
    for row in ancestor_records(folder, gap_method):
        true_node = match.get(str(row["node"]))
        site = int(str(row["site"])) - 1
        column = per_node.get(str(row["node"]), {}).get(site, mapping.get(site))
        if true_node is None or column is None:
            continue
        true = true_ancestors[true_node][column].upper()
        if true in AMINO_ACIDS:
            yield row | {"true": true, "correct": int(true == row["ml"])}


def esr_records(folder: Path, esr_sites: Path) -> Iterator[dict[str, object]]:
    """Held-out tips' scored columns (esr_sites.tsv from `phylobot.esr`) with covariates."""
    combination = Combination(folder)
    cache: dict[str, list[dict[str, float]]] = {}
    with esr_sites.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            tip = row["tip"]
            if tip not in cache:
                place = calibration.tip_place(combination.tree, tip)
                if place is None:
                    raise ValueError(f"{tip} is not a tip of {folder}/tree.nwk")
                others = [s for name, s in combination.alignment.items() if name != tip]
                cache[tip] = calibration.covariate_rows(place, others, combination.rates)
            site = int(row["site"])
            probs = {a: float(row[a]) for a in AMINO_ACIDS}
            yield _row(f"esr-{row['level']}", tip, site, probs, cache[tip][site - 1], row["true"])


def write(path: Path, rows: Iterable[dict[str, object]]) -> int:
    """Write records as TSV; returns how many."""
    n = 0
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            n += 1
    return n
