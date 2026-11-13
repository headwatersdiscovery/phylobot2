"""Species trees: from the Open Tree of Life's answer or an uploaded Newick file, to a
rooted tree whose leaves are exactly the job's species.

A species tree here is a topology only (`Node`): the undated reconciliation models
(GeneRax UndatedDL / UndatedDTL) ignore branch lengths, so uploaded lengths are dropped.
Leaves carry species names; internal nodes keep a taxon name where Open Tree gives one,
for display. Unary nodes are always collapsed. Polytomies are kept: they are resolved
later, from the gene trees (`polytomies`).
"""

import re
from collections.abc import Iterator
from dataclasses import dataclass, field

import dendropy

from phylobot.species.opentree import Induced
from phylobot.trees import newick_problem

MAX_TREE_BYTES = 5_000_000
MAX_DEPTH = 500  # Open Tree's deepest chains (Gallus to the root) are about 70
# "Saccharomyces_cerevisiae_ott356221", "'Candida (genus ...) ott1034970'", "mrcaott1ott2"
OPEN_TREE_LABEL = re.compile(r"^(?:(?P<name>.*?)[ _])?(?P<id>ott\d+|mrcaott\d+ott\d+)$")
PLAIN = re.compile(r"^[A-Za-z0-9_.\-]+$")


class SpeciesTreeError(Exception):
    pass


@dataclass
class Node:
    label: str = ""
    children: list[Node] = field(default_factory=list)

    def is_leaf(self) -> bool:
        return not self.children

    def preorder(self) -> Iterator[Node]:
        yield self
        for child in self.children:
            yield from child.preorder()

    def leaves(self) -> list[str]:
        return [n.label for n in self.preorder() if n.is_leaf()]

    def polytomies(self) -> list[Node]:
        return [n for n in self.preorder() if len(n.children) > 2]

    def newick(self) -> str:
        return write(self) + ";"


def key(name: str) -> str:
    """Species names compare with underscores as spaces and runs of spaces as one."""
    return " ".join(name.replace("_", " ").split())


def quote(label: str) -> str:
    text = label.replace(" ", "_")
    if not label or PLAIN.fullmatch(text):
        return text
    return "'" + label.replace("'", "''") + "'"


def write(node: Node) -> str:
    if node.is_leaf():
        return quote(node.label)
    inner = ",".join(write(child) for child in node.children)
    return f"({inner}){quote(node.label) if node.label else ''}"


def parse(text: str, what: str = "species tree") -> Node:
    """Any Newick tree as `Node`s. As Newick defines, an unquoted underscore is a space
    ("Homo_sapiens" is "Homo sapiens"); quoted labels are kept as written."""
    if len(text.encode()) > MAX_TREE_BYTES:
        raise SpeciesTreeError(f"{what} is larger than {MAX_TREE_BYTES // 1_000_000} MB")
    if problem := newick_problem(text):
        raise SpeciesTreeError(f"{what}: {problem}")
    depth = deepest = 0
    for c in text:  # bound the recursion below (`convert`, `write`)
        depth += (c == "(") - (c == ")")
        deepest = max(deepest, depth)
    if deepest > MAX_DEPTH:
        raise SpeciesTreeError(f"{what} is nested more than {MAX_DEPTH} levels deep")
    try:
        tree = dendropy.Tree.get(data=text, schema="newick", suppress_leaf_node_taxa=True)
    except Exception as exc:
        raise SpeciesTreeError(f"{what} is not valid Newick: {exc}") from exc

    def convert(source: dendropy.Node) -> Node:
        return Node(source.label or "", [convert(c) for c in source.child_nodes()])

    return convert(tree.seed_node)


def collapse(node: Node) -> Node:
    """Remove unary nodes; a merged node keeps the lower node's label if it has one."""
    while len(node.children) == 1:
        child = node.children[0]
        if not child.label:
            child.label = node.label
        node = child
    node.children = [collapse(child) for child in node.children]
    return node


def prune(node: Node, keep: set[str]) -> Node | None:
    """The tree with only the leaves in `keep` (unary nodes collapsed); None if empty."""
    if node.is_leaf():
        return node if node.label in keep else None
    children = [c for c in (prune(child, keep) for child in node.children) if c is not None]
    if not children:
        return None
    return collapse(Node(node.label, children))


def from_upload(text: str, species: set[str]) -> Node:
    """An uploaded species tree pruned to the leaves naming one of `species` (compared
    with `key`); leaves keep the tree's spelling. Species not in it are simply absent."""
    tree = parse(text, "the uploaded species tree")
    seen: set[str] = set()
    for node in tree.preorder():
        if node.is_leaf():
            if not node.label:
                raise SpeciesTreeError("the uploaded species tree has a leaf without a name")
            if key(node.label) in seen:
                raise SpeciesTreeError(f"species {node.label!r} appears twice in the uploaded tree")
            seen.add(key(node.label))
    wanted = {key(name) for name in species}
    keep = {leaf for leaf in tree.leaves() if key(leaf) in wanted}
    pruned = prune(tree, keep)
    if pruned is None or len(pruned.leaves()) < 2:
        raise SpeciesTreeError("the uploaded species tree has fewer than two of the species")
    return pruned


def from_open_tree(
    induced: Induced, species_of: dict[int, str]
) -> tuple[Node, dict[int, list[int]]]:
    """The induced subtree with each requested taxon (ott id -> species name) as a leaf
    named for its species, and the requested taxa left out because they contain others
    (e.g. a genus asked for beside one of its species), with the taxa inside each."""
    tree = parse(induced.newick, "the Open Tree of Life's tree")
    by_id: dict[str, Node] = {}
    for node in tree.preorder():
        match = OPEN_TREE_LABEL.match(node.label)
        if match:
            by_id[match["id"]] = node
            name = "" if match["id"].startswith("mrca") else match["name"] or ""
            node.label = name.replace("_", " ")
        else:
            node.label = ""
    targets: dict[int, Node] = {}
    for ott_id in species_of:
        if ott_id in induced.pruned:
            continue
        target = by_id.get(induced.broken.get(ott_id, f"ott{ott_id}"))
        if target is None:
            raise SpeciesTreeError(f"the Open Tree of Life's tree lacks ott{ott_id}")
        targets[ott_id] = target
    nested: dict[int, list[int]] = {}
    for ott_id, node in targets.items():
        below = {id(n) for n in node.preorder()}
        inside = sorted(
            other for other, target in targets.items() if other != ott_id and id(target) in below
        )
        if inside:
            nested[ott_id] = inside
    keep: set[str] = set()
    for ott_id, node in targets.items():
        if ott_id not in nested:
            node.label, node.children = species_of[ott_id], []
            keep.add(species_of[ott_id])
    pruned = prune(tree, keep)
    if pruned is None or sorted(pruned.leaves()) != sorted(keep):
        raise SpeciesTreeError("the Open Tree of Life's tree does not contain every species")
    return pruned, nested
