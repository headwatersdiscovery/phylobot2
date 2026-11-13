"""Resolving the species tree's polytomies, which GeneRax refuses.

By default (owner decision) a polytomy is resolved by the job's own gene trees
("guided"); the user can override with a seeded arbitrary resolution ("arbitrary"), or
upload a species tree that is already resolved.

Guided, at a polytomy whose children hold species C1..Ck: every gene-tree clade whose
species all lie under the polytomy and touch at least two of its children, but not
all, votes for grouping those children (once per gene tree). Groupings are accepted
greedily, most votes first, if compatible with those already accepted (nested or
disjoint); ties go to the job seed. What the votes leave unresolved is joined at
random, also from the seed. Votes come only from where the gene trees speak, so no
duplication or loss is introduced that they do not already imply; where they are
silent the resolution is arbitrary and recorded as such, so events inferred there can
be labelled "species tree unresolved here".

Gene trees are unrooted unless an outgroup roots them. Without one, a split votes for
both of its sides, which is enough to resolve polytomies below the root but cannot
place the root itself: the report says how the root was placed ("root").
"""

import random
from dataclasses import asdict, dataclass, field
from typing import Any

from phylobot.rooting import root_on_outgroup
from phylobot.species.speciestree import Node
from phylobot.trees import load_tree

GUIDED, ARBITRARY = "guided", "arbitrary"
METHODS = (GUIDED, ARBITRARY)

Grouping = frozenset[int]  # indices of a polytomy's children


class PolytomyError(Exception):
    pass


@dataclass(frozen=True)
class Added:
    """A node the resolution added."""

    species: list[str]  # the species below it
    by: str  # "gene trees" or "arbitrary"
    votes: int  # gene trees with a clade supporting it (0 if arbitrary)


@dataclass
class Polytomy:
    label: str  # the taxon name Open Tree gave the node, if any
    root: bool
    children: list[list[str]]  # species below each child, as it was
    informative: int  # gene trees that voted at this polytomy
    added: list[Added] = field(default_factory=list)


def gene_tree_species_sets(
    newick: str, species_of: dict[str, str], outgroup: set[str]
) -> set[frozenset[str]]:
    """The species sets of a gene tree's clades (outgroup-rooted) or of both sides of
    each of its splits (unrooted)."""
    if outgroup:
        newick, _ = root_on_outgroup(newick, outgroup)
    tree = load_tree(newick)
    genes = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    unknown = genes - set(species_of)
    if unknown:
        shown = ", ".join(sorted(unknown)[:10])
        raise PolytomyError(f"gene tree leaves without a species: {shown}")
    every = frozenset(species_of[g] for g in genes)
    sets: set[frozenset[str]] = set()
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            continue
        below = [leaf.taxon.label for leaf in node.leaf_iter()]
        sets.add(frozenset(species_of[g] for g in below))
        if not outgroup and len(below) > 1:  # a leaf's other side is no split at all
            rest = genes - set(below)
            if len(rest) > 1:
                sets.add(frozenset(species_of[g] for g in rest))
    sets.discard(every)
    return sets


def compatible(a: Grouping, b: Grouping) -> bool:
    return not (a & b) or a <= b or b <= a


def resolve(
    tree: Node,
    gene_trees: list[str],
    species_of: dict[str, str],
    seed: int,
    method: str = GUIDED,
    outgroup: set[str] | None = None,
) -> tuple[Node, dict[str, Any]]:
    """The species tree with every polytomy resolved (the input is modified), and a
    report of what was added and why. `species_of` maps gene-tree leaves to species."""
    if method not in METHODS:
        raise PolytomyError(f"unknown polytomy method {method!r} (choose {' or '.join(METHODS)})")
    outgroup = outgroup or set()
    leaves = set(tree.leaves())
    strays = set(species_of.values()) - leaves
    if strays:
        raise PolytomyError(f"species not in the species tree: {', '.join(sorted(strays)[:10])}")
    per_tree = (
        [gene_tree_species_sets(t, species_of, outgroup) for t in gene_trees]
        if method == GUIDED
        else []
    )
    rng = random.Random(seed)
    polytomies = []
    for node in tree.polytomies():  # preorder, so the order (and the seed's use) is fixed
        polytomies.append(resolve_one(node, node is tree, per_tree, rng))
    if not any(p.root for p in polytomies):
        root = "species tree"
    elif method == GUIDED and outgroup:
        root = "gene trees rooted on the outgroup"
    else:
        root = "arbitrary"
        if method == GUIDED:
            root += " (unrooted gene trees cannot place it)"
    report = {
        "method": method,
        "seed": seed,
        "gene_trees": len(gene_trees) if method == GUIDED else 0,
        "outgroup": sorted(outgroup),
        "root": root,
        "polytomies": [asdict(p) for p in polytomies],
    }
    return tree, report


def resolve_one(
    node: Node, root: bool, per_tree: list[set[frozenset[str]]], rng: random.Random
) -> Polytomy:
    children = sorted(node.children, key=lambda c: sorted(c.leaves()))
    species = [frozenset(c.leaves()) for c in children]
    under = frozenset().union(*species)
    k = len(children)
    votes: dict[Grouping, int] = {}
    informative = 0
    for sets in per_tree:
        here: set[Grouping] = set()
        for found in sets:
            if found <= under:
                grouping = frozenset(i for i, s in enumerate(species) if found & s)
                if 2 <= len(grouping) < k:
                    here.add(grouping)
        informative += bool(here)
        for grouping in here:
            votes[grouping] = votes.get(grouping, 0) + 1
    candidates = sorted(votes, key=sorted)
    rng.shuffle(candidates)
    candidates.sort(key=lambda g: -votes[g])  # stable: equal votes keep the seed's order
    accepted: list[Grouping] = []
    for grouping in candidates:
        if all(compatible(grouping, other) for other in accepted):
            accepted.append(grouping)

    polytomy = Polytomy(
        label=node.label,
        root=root,
        children=[sorted(s) for s in species],
        informative=informative,
    )

    def add(new: Node, by: str, count: int) -> Node:
        polytomy.added.append(Added(sorted(new.leaves()), by, count))
        return new

    def binary(items: list[Node]) -> list[Node]:
        items = list(items)
        while len(items) > 2:
            i, j = sorted(rng.sample(range(len(items)), 2))
            joined = add(Node("", [items[i], items[j]]), "arbitrary", 0)
            items = [x for n, x in enumerate(items) if n not in (i, j)] + [joined]
        return items

    def build(group: Grouping) -> list[Node]:
        """The two children of the node holding `group`."""
        inside = [g for g in accepted if g < group]
        tops = [g for g in inside if not any(g < other for other in inside)]
        covered = frozenset().union(*tops) if tops else frozenset()
        items = [children[i] for i in sorted(group - covered)]
        for top in sorted(tops, key=sorted):
            items.append(add(Node("", build(top)), "gene trees", votes[top]))
        return binary(items)

    node.children = build(frozenset(range(k)))
    return polytomy
