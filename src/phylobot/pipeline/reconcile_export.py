"""A combination's reconciliation, on its published tree (RECONCILE_EXPORT in main.nf).

    phylobot-reconcile-export NAMES.tsv SPECIES_NODES.tsv TREE RECONCILIATION.xml \
        GENERAX_STATS.txt GENERAX_MODEL.txt OUTDIR [--name-nodes]

TREE is the combination's tree with the node names the job publishes (the
reconstruction's Node1...; with --name-nodes, a tree whose nodes `export.name_nodes`
names, as a trees-only job's export does). RECONCILIATION.xml is GeneRax's RecPhyloXML
for that tree (`phylobot.reconcile`).

Writes to OUTDIR, with the sequences' and species' own names:
- reconciliation.tsv: one row per node of tree.nwk (internal nodes and leaves, preorder):
  node, event (speciation, duplication, transfer, leaf), species (the species-tree node
  it maps to), losses (copies lost on the branch above it, in GeneRax's rooting) and
  lost_in (the species that lost them); and last GeneRax's root, named "(root)" (it sits
  on the branch between root_branch's two nodes), with its event and the losses above it;
- reconciliation.json: totals, the root's branch (the nodes either side), likelihoods,
  rates, the reconciliation model and the substitution model GeneRax was given;
- reconciliation.nwk: the tree rooted where GeneRax roots it.
"""

import argparse
import json
import re
from pathlib import Path
from typing import Any

from phylobot import reconcile
from phylobot.pipeline.export import name_nodes, read_names, restorer
from phylobot.pipeline.species_tree import read_names as read_species
from phylobot.species.table import MAX_CELL

EVENT_NAMES = {"S": "speciation", "D": "duplication", "T": "transfer", "leaf": "leaf"}
REC_MODEL = "UndatedDL"
NUMBER = r"-?[0-9.]+(?:e[-+]?[0-9]+)?"


def parse_stats(text: str) -> dict[str, Any]:
    """GeneRax's per-family stats.txt: 'libpll_ll rec_ll' then 'Reconciliation rates = d l'."""
    stats: dict[str, Any] = {}
    lines = text.splitlines()
    if lines and (match := re.fullmatch(rf"\s*({NUMBER})\s+({NUMBER})\s*", lines[0])):
        stats["sequence_log_likelihood"] = float(match[1])
        stats["reconciliation_log_likelihood"] = float(match[2])
    for line in lines[1:]:
        if line.startswith("Reconciliation rates ="):
            values = [float(v) for v in re.findall(NUMBER, line.split("=", 1)[1])]
            stats["rates"] = dict(zip(("duplication", "loss", "transfer"), values, strict=False))
    return stats


def export(
    names_path: Path,
    species_path: Path,
    tree: Path,
    xml: Path,
    stats_path: Path,
    model_path: Path,
    outdir: Path,
    name_tree_nodes: bool = False,
) -> None:
    names = read_names(names_path)
    species = {i: name for i, (name, _) in read_species(species_path.read_text()).items()}
    newick = tree.read_text().strip()
    if name_tree_nodes:
        newick = name_nodes(newick)
    rec = reconcile.parse(xml.read_text())
    rows, ends, rooted = reconcile.annotate(newick, rec)
    rows.append(reconcile.root_row(rec))
    restore = restorer(names)

    def shown(node: str) -> str:
        return names.get(node, node)

    def where(sp: str) -> str:
        return species.get(sp, sp)[:MAX_CELL]

    outdir.mkdir(parents=True, exist_ok=True)
    lines = ["node\tevent\tspecies\tlosses\tlost_in"]
    for row in rows:
        lost = "; ".join(where(s) for s in row.losses)
        lines.append(
            f"{shown(row.node)}\t{EVENT_NAMES[row.event]}\t{where(row.species)}\t"
            f"{len(row.losses)}\t{lost}"
        )
    (outdir / "reconciliation.tsv").write_text("\n".join(lines) + "\n")
    (outdir / "reconciliation.nwk").write_text(restore(rooted) + "\n")
    summary = {
        **reconcile.totals(rows),
        "root_branch": [shown(end) for end in ends],
        "root_event": EVENT_NAMES[rows[-1].event],
        "reconciliation_model": REC_MODEL,
        "generax_model": model_path.read_text().strip(),
        **parse_stats(stats_path.read_text()),
    }
    (outdir / "reconciliation.json").write_text(json.dumps(summary, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(usage=__doc__)
    for name in ("names", "species", "tree", "xml", "stats", "model", "outdir"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--name-nodes", action="store_true")
    args = parser.parse_args()
    export(
        args.names, args.species, args.tree, args.xml, args.stats, args.model, args.outdir,
        args.name_nodes,
    )  # fmt: skip


if __name__ == "__main__":
    main()
