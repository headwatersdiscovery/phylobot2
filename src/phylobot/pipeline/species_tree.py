"""The job's species tree, made ready for reconciliation (SPECIES_TREE in main.nf).

    phylobot-species-tree SPECIES.safe.tsv SPECIES_TREE.safe.nwk SPECIES_NAMES.tsv \
        OUTGROUP.txt OUTDIR --seed N --method guided|arbitrary --gene-trees TREE...

Runs once per job, after every tree search: it resolves the species tree's polytomies
(`phylobot.species.polytomies`; guided by all the combinations' gene trees unless the
user chose arbitrary) and numbers the nodes it adds spr001.... Inputs are VALIDATE's
safe-ID files (sequences s0001..., species sp001..., species-tree nodes spn001...).

Writes to OUTDIR:
- species_tree.resolved.nwk: safe IDs, every internal node labelled (GeneRax keeps the
  labels, so events name our nodes);
- mapping.link: GeneRax's species -> genes map (PHYLDOG format, `sp001:s0001;s0004`);
- species_nodes.tsv: every species-tree node: id, name (the species, Open Tree's taxon
  name, or "ancestor of A and B"), kind (species, taxon, node, or resolved), how a
  resolved node was made (gene trees or arbitrary) and its gene-tree votes;
- species_tree.nwk: the resolved tree with display names, for the job page;
- polytomies.json: the resolution report, with species names.
"""

import argparse
import csv
import io
import json
from pathlib import Path

from phylobot.species import polytomies, speciestree, table
from phylobot.species.speciestree import Node


def read_names(text: str) -> dict[str, tuple[str, str]]:
    """species_names.tsv: {id: (name, kind)}."""
    reader = csv.DictReader(io.StringIO(text), delimiter="\t", quoting=csv.QUOTE_NONE)
    return {row["id"]: (row["name"], row["kind"]) for row in reader}


def first_leaf(node: Node) -> str:
    return min(node.leaves())


def write(
    mapping_path: Path,
    tree_path: Path,
    names_path: Path,
    outgroup_path: Path,
    gene_trees: list[Path],
    outdir: Path,
    seed: int,
    method: str,
) -> None:
    rows = table.read(mapping_path.read_text())
    species_of = {row.sequence: row.species for row in rows}
    names = read_names(names_path.read_text())
    tree = speciestree.parse(tree_path.read_text())
    outgroup = {line.strip() for line in outgroup_path.read_text().splitlines() if line.strip()}
    resolved, report = polytomies.resolve(
        tree, [p.read_text() for p in gene_trees], species_of, seed, method, outgroup
    )

    added = {frozenset(a["species"]): a for p in report["polytomies"] for a in p["added"]}
    nodes: list[tuple[str, str, str, str, str]] = []  # id, name, kind, by, votes
    count = 0

    def display(node: Node) -> str:
        if node.is_leaf():
            return names[node.label][0]
        name, _ = names.get(node.label, ("", ""))
        if name:
            return name
        a, b = sorted(names[first_leaf(c)][0] for c in node.children[:2])
        return f"ancestor of {a} and {b}"

    for node in resolved.preorder():
        if node.is_leaf():
            nodes.append((node.label, names[node.label][0], "species", "", ""))
            continue
        if not node.label:
            count += 1
            node.label = f"spr{count:03d}"
            how = added[frozenset(node.leaves())]
            nodes.append((node.label, display(node), "resolved", how["by"], str(how["votes"])))
        else:
            name, kind = names[node.label]
            nodes.append((node.label, display(node), kind, "", ""))

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "species_tree.resolved.nwk").write_text(resolved.newick() + "\n")
    genes: dict[str, list[str]] = {}
    for row in rows:
        genes.setdefault(row.species, []).append(row.sequence)
    (outdir / "mapping.link").write_text(
        "".join(f"{sp}:{';'.join(g)}\n" for sp, g in sorted(genes.items()))
    )
    (outdir / "species_nodes.tsv").write_text(
        "id\tname\tkind\tresolved_by\tvotes\n" + "".join("\t".join(n) + "\n" for n in nodes)
    )
    shown = {i: name for i, name, *_ in nodes}

    def named(node: Node) -> Node:
        kind = next(n[2] for n in nodes if n[0] == node.label)
        label = shown[node.label] if node.is_leaf() or kind == "taxon" else ""
        return Node(label, [named(c) for c in node.children])

    (outdir / "species_tree.nwk").write_text(named(resolved).newick() + "\n")
    for p in report["polytomies"]:
        p["label"] = shown.get(p["label"], p["label"])
        p["children"] = [sorted(shown[s] for s in group) for group in p["children"]]
        for a in p["added"]:
            a["species"] = sorted(shown[s] for s in a["species"])
    report["outgroup"] = len(report["outgroup"])  # sequence IDs mean nothing on the page
    (outdir / "polytomies.json").write_text(json.dumps(report, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(usage=__doc__)
    for name in ("mapping", "tree", "names", "outgroup", "outdir"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--method", choices=polytomies.METHODS, required=True)
    parser.add_argument("--gene-trees", type=Path, nargs="*", default=[])
    args = parser.parse_args()
    write(
        args.mapping, args.tree, args.names, args.outgroup, args.gene_trees, args.outdir,
        args.seed, args.method,
    )  # fmt: skip


if __name__ == "__main__":
    main()
