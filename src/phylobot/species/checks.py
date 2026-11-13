"""A job's species files, checked and renamed to safe IDs (called by VALIDATE).

The job brings species.tsv and species_tree.nwk (made by `phylobot-species resolve`,
by the web tier, or by hand). Both are user input. Every sequence must have a species
and every species must be a leaf of the tree; species the job does not use are pruned.
Species become sp001... in the files tools read, as sequences become s0001..., and the
tree's internal nodes spn001... (preorder); their names (a species, or a taxon name
such as Open Tree's "Saccharomycetaceae", or none) stay in species_names.tsv.
"""

from dataclasses import dataclass

from phylobot.species import speciestree, table
from phylobot.species.speciestree import Node, SpeciesTreeError, key
from phylobot.species.table import TableError

MAX_LISTED = 10


class SpeciesError(Exception):
    pass


@dataclass(frozen=True)
class Prepared:
    mapping: str  # species.safe.tsv: a species table, columns sequence and species (IDs)
    tree: str  # species_tree.safe.nwk: leaves sp001..., internal nodes spn001...
    names: str  # species_names.tsv: id, name, kind (species, taxon or node)
    report: dict[str, object]


def listing(names: set[str] | list[str]) -> str:
    shown = sorted(names)[:MAX_LISTED]
    more = len(names) - len(shown)
    return ", ".join(repr(n) for n in shown) + (f" and {more} more" if more > 0 else "")


def relabel(node: Node, ids: dict[str, str], internal: list[tuple[str, str]]) -> Node:
    """The tree with leaves renamed to their IDs and internal nodes numbered spn001...
    in preorder; `internal` collects (id, the node's own name or "")."""
    if node.is_leaf():
        return Node(ids[node.label])
    label = f"spn{len(internal) + 1:03d}"
    internal.append((label, node.label))
    return Node(label, [relabel(c, ids, internal) for c in node.children])


def prepare(table_text: str, tree_text: str, safe: dict[str, str]) -> Prepared:
    """`safe` maps the job's sequence names to their safe IDs."""
    try:
        rows = table.read(table_text)
        tree = speciestree.parse(tree_text)
    except (TableError, SpeciesTreeError) as exc:
        raise SpeciesError(str(exc)) from None
    listed = {row.sequence for row in rows}
    missing, extra = set(safe) - listed, listed - set(safe)
    if missing or extra:
        parts = [f"no species for {listing(missing)}"] if missing else []
        parts += [f"not among the sequences: {listing(extra)}"] if extra else []
        raise SpeciesError("species table: " + "; ".join(parts))
    lacking = {row.sequence for row in rows if not row.species}
    if lacking:
        raise SpeciesError(f"species table: sequences without a species: {listing(lacking)}")

    leaves = tree.leaves()
    if "" in leaves:
        raise SpeciesError("species tree: a leaf has no name")
    by_key: dict[str, str] = {}
    for leaf in leaves:
        if key(leaf) in by_key:
            raise SpeciesError(f"species tree: species {leaf!r} appears twice")
        by_key[key(leaf)] = leaf
    absent = {row.species for row in rows if key(row.species) not in by_key}
    if absent:
        raise SpeciesError(f"species not in the species tree: {listing(absent)}")
    used = {by_key[key(row.species)] for row in rows}
    if len(used) < 2:
        raise SpeciesError("reconciliation needs sequences from at least two species")
    pruned = speciestree.prune(tree, used)
    assert pruned is not None
    order = pruned.leaves()
    ids = {name: f"sp{index:03d}" for index, name in enumerate(order, 1)}
    mapping = "sequence\tspecies\n" + "".join(
        f"{safe[row.sequence]}\t{ids[by_key[key(row.species)]]}\n" for row in rows
    )
    internal: list[tuple[str, str]] = []
    safe_tree = relabel(pruned, ids, internal)
    names = "id\tname\tkind\n" + "".join(f"{ids[name]}\t{name}\tspecies\n" for name in order)
    names += "".join(f"{i}\t{n}\t{'taxon' if n else 'node'}\n" for i, n in internal)
    copies: dict[str, int] = {}
    for row in rows:
        copies[by_key[key(row.species)]] = copies.get(by_key[key(row.species)], 0) + 1
    report: dict[str, object] = {
        "species": len(order),
        "pruned_species": len(leaves) - len(order),
        "polytomies": len(pruned.polytomies()),
        "species_with_several_sequences": sum(n > 1 for n in copies.values()),
    }
    return Prepared(mapping, safe_tree.newick() + "\n", names, report)
