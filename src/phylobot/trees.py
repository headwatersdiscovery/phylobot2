"""Tree utilities: loading, splits, clade-based node matching, Robinson-Foulds distance.

Different trees name their internal nodes differently, so a node is identified by
the leaves it separates: its split (unrooted) or, once the tree is rooted on the
outgroup, its clade.
"""

import re
from collections.abc import Iterable

import dendropy
from dendropy.calculate import treecompare

Clade = frozenset[str]


def newick_problem(text: str) -> str | None:
    """What is wrong with malformed Newick that parsers tend to repair silently, or None.

    Scans outside quoted labels and [comments]: no empty elements (",," "(," ",)"
    "()"), no element with two branch lengths, and exactly one terminating ';'.
    """
    stripped: list[str] = []
    quote = comment = False
    for c in text:
        if quote:
            quote = c != "'"
            stripped.append("x")  # a quoted label counts as a name
        elif comment:
            comment = c != "]"
        elif c == "'":
            quote = True
            stripped.append("x")
        elif c == "[":
            comment = True
        elif not c.isspace():
            stripped.append(c)
    if quote or comment:
        return "unterminated quote or comment"
    body = "".join(stripped)
    if not body.endswith(";") or body.count(";") != 1:
        return "must be a single tree ending with ';'"
    for bad in (",,", "(,", ",)", "()"):
        if bad in body:
            return f"empty element ({bad!r})"
    for element in re.split(r"[(),;]", body):
        if element.count(":") > 1:
            return f"more than one branch length in {element!r}"
    return None


def load_tree(
    newick: str,
    taxa: dendropy.TaxonNamespace | None = None,
    prune: Iterable[str] = (),
) -> dendropy.Tree:
    tree = dendropy.Tree.get(
        data=newick,
        schema="newick",
        preserve_underscores=True,
        taxon_namespace=dendropy.TaxonNamespace() if taxa is None else taxa,
    )
    if prune := list(prune):
        tree.prune_taxa_with_labels(prune)
        tree.suppress_unifurcations()
    return tree


def clades(tree: dendropy.Tree, outgroup: str) -> dict[Clade, str]:
    """Map each ingroup clade to the label of the internal node that is its ancestor.

    The tree is re-rooted on the outgroup's branch (a copy; the input is untouched),
    so clades are defined the same way regardless of how the file was rooted.
    """
    tree = tree.clone(depth=1)
    leaf = tree.find_node_with_taxon_label(outgroup)
    if leaf is None:
        raise ValueError(f"outgroup {outgroup!r} not in tree")
    tree.reroot_at_edge(leaf.edge, update_bipartitions=True)
    tree.suppress_unifurcations()
    result: dict[Clade, str] = {}
    for node in tree.postorder_internal_node_iter():
        members = frozenset(n.taxon.label for n in node.leaf_iter())
        if outgroup in members:
            continue  # the root of an outgroup-rooted tree is not reconstructable
        result[members] = node.label or ""
    return result


def robinson_foulds(tree_a: str, tree_b: str, prune: Iterable[str] = ()) -> int:
    """Unrooted Robinson-Foulds distance (number of splits in one tree but not the other)."""
    taxa = dendropy.TaxonNamespace()
    a = load_tree(tree_a, taxa, prune)
    b = load_tree(tree_b, taxa, prune)
    for tree in (a, b):
        tree.is_rooted = False
        tree.encode_bipartitions()
    distance: int = treecompare.symmetric_difference(a, b)
    return distance


def internal_splits(tree: dendropy.Tree) -> list[tuple[str, Clade]]:
    """(label, leaf set) for each internal node, in postorder, as the file roots it."""
    return [
        (node.label or "", frozenset(leaf.taxon.label for leaf in node.leaf_iter()))
        for node in tree.postorder_internal_node_iter()
    ]


def canonical_split(leaves: Clade, everyone: Clade) -> Clade:
    """Unrooted identity of a split: the side without the alphabetically first taxon."""
    return everyone - leaves if min(everyone) in leaves else leaves
