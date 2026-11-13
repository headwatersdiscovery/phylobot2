"""Reconciliation events on our trees (pure Python).

GeneRax reconciles a combination's gene tree with the species tree (RECONCILE in the
pipeline, `--strategy EVAL`: the topology is kept). It roots the gene tree and labels
every node with an event: a speciation (S), a duplication (D) or a transfer (T), and the
species-tree node where it happened; losses sit on branches. This module reads its
RecPhyloXML output and puts those events on *our* nodes, matched by the subtrees that
meet at each node (as `gaps.node_partitions` does), so it does not depend on where
either tree is rooted or how its nodes are named.

Losses: GeneRax writes a lost copy as a one-child "speciation" step on the branch above
a gene node, with a `loss` child naming the species that lost it. They are counted on
the branch above the node they lead to, in GeneRax's rooting.

GeneRax's root sits on one of our branches, so it is no node of ours: its event (a
duplication, when a family's two paralog groups meet there) and the losses above it are
a row of its own, named ROOT.

The paralog check (`paralogy`) asks of one ancestor what another tree of the same
alignment and model says about it (GeneRax's reconciled tree for the sequence-only
trees, and the sequence-only tree for GeneRax's): is its clade there too, or do its
sequences fall on both sides of a duplication there? On simulated families either
tree can be the wrong one, so the page reports a disagreement, not a verdict.

Weak duplications (`overlaps`): a duplication whose two copies share few species
(species overlap below WEAK_OVERLAP) is usually a gene-tree error explained as a
duplication plus losses, not a real duplication (Huerta-Cepas et al. 2007); users see
it as "doubtful". On simulated families, most false duplications fall below 0.1 and few
real ones do.
"""

import xml.etree.ElementTree as ElementTree
from dataclasses import dataclass, field

import dendropy

from phylobot.trees import canonical_split, load_tree

Leaves = frozenset[str]

# RecPhyloXML event elements -> our event codes. A transfer leaves its donor with a
# "branchingOut" (or "bifurcationOut" outside the species tree); "transferBack" marks
# the receiving lineage and precedes that node's own event.
EVENTS = {
    "speciation": "S",
    "duplication": "D",
    "branchingOut": "T",
    "bifurcationOut": "T",
    "leaf": "leaf",
}
ROOT = "(root)"  # GeneRax's root in reconciliation.tsv (parentheses: never a node name)
MAX_XML_BYTES = 50_000_000
MAX_DEPTH = 10_000


class ReconcileError(Exception):
    pass


@dataclass
class GeneNode:
    below: Leaves  # leaves under it, in GeneRax's rooting
    event: str  # S, D, T or leaf
    species: str  # the species-tree node it maps to
    children: tuple[Leaves, ...] = ()
    losses: list[str] = field(default_factory=list)  # species that lost a copy above it


@dataclass
class Reconciliation:
    nodes: dict[Leaves, GeneNode]  # by the leaves below each node
    top: Leaves  # every leaf: the root's key in `nodes`

    @property
    def root(self) -> tuple[Leaves, Leaves]:
        """The two sides of GeneRax's root."""
        a, b = self.nodes[self.top].children
        return (a, b) if min(a) < min(b) else (b, a)


def local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def child(element: ElementTree.Element, tag: str) -> ElementTree.Element | None:
    return next((c for c in element if local(c.tag) == tag), None)


def parse(xml: str) -> Reconciliation:
    """GeneRax's `<family>_reconciliated.xml` (RecPhyloXML) for one gene tree."""
    if len(xml.encode()) > MAX_XML_BYTES:
        raise ReconcileError("the reconciliation file is too large")
    try:
        document = ElementTree.fromstring(xml)
    except ElementTree.ParseError as exc:
        raise ReconcileError(f"the reconciliation file is not valid XML: {exc}") from None
    gene_tree = next((e for e in document.iter() if local(e.tag) == "recGeneTree"), None)
    phylogeny = child(gene_tree, "phylogeny") if gene_tree is not None else None
    top = child(phylogeny, "clade") if phylogeny is not None else None
    if top is None:
        raise ReconcileError("the reconciliation file has no gene tree")
    nodes: dict[Leaves, GeneNode] = {}

    def walk(clade: ElementTree.Element, depth: int) -> Leaves | str:
        """The leaves below `clade`, or the species of a loss."""
        if depth > MAX_DEPTH:
            raise ReconcileError("the reconciled gene tree is nested too deeply")
        name_element = child(clade, "name")
        name = (name_element.text or "").strip() if name_element is not None else ""
        rec = child(clade, "eventsRec")
        events = list(rec) if rec is not None else []
        if not events:
            raise ReconcileError(f"gene tree node {name!r} has no event")
        last = events[-1]
        kind, species = local(last.tag), last.get("speciesLocation", "")
        if kind == "loss":
            return species
        results = [walk(c, depth + 1) for c in clade if local(c.tag) == "clade"]
        lost = [r for r in results if isinstance(r, str)]
        kept = [r for r in results if not isinstance(r, str)]
        if kind == "leaf":
            if results or not name:
                raise ReconcileError(f"leaf {name!r} is malformed")
            below: Leaves = frozenset({name})
            if below in nodes:
                raise ReconcileError(f"leaf {name!r} appears twice")
            nodes[below] = GeneNode(below, "leaf", species)
            return below
        if len(kept) == 1:  # a speciation followed by a loss: no node of its own
            nodes[kept[0]].losses[:0] = lost
            return kept[0]
        if len(kept) != 2 or lost:
            raise ReconcileError(f"a gene tree node has {len(kept)} children")
        if kind not in EVENTS:
            raise ReconcileError(f"unknown reconciliation event {kind!r}")
        below = kept[0] | kept[1]
        nodes[below] = GeneNode(below, EVENTS[kind], species, (kept[0], kept[1]))
        return below

    result = walk(top, 0)
    if isinstance(result, str) or len(nodes[result].children) != 2:
        raise ReconcileError("the reconciled gene tree has fewer than two leaves")
    return Reconciliation(nodes, result)


@dataclass(frozen=True)
class Row:
    node: str  # our label (internal node) or the leaf's name
    leaf: bool
    event: str  # S, D, T or leaf
    species: str  # species-tree node
    # Species that lost a copy on the branch above, in GeneRax's rooting (for the two
    # nodes either side of the root, on their halves of the root branch).
    losses: tuple[str, ...]


def annotate(newick: str, rec: Reconciliation) -> tuple[list[Row], list[str], str]:
    """Our tree's nodes with their events: (rows in our tree's preorder, our labels on
    either side of GeneRax's root, our tree rooted there as Newick). Refuses a tree whose
    topology differs from the reconciled one."""
    tree = load_tree(newick)
    everyone = frozenset(leaf.taxon.label for leaf in tree.leaf_node_iter())
    if everyone != rec.top:
        raise ReconcileError("the reconciled tree and the tree have different leaves")
    # A node of an unrooted tree is the set of subtrees meeting at it. The reconciled
    # root meets only two, so it matches none of ours: it sits on one of our branches.
    by_parts: dict[frozenset[Leaves], GeneNode] = {}
    for below, gene in rec.nodes.items():
        if gene.children and below != rec.top:
            by_parts[frozenset({*gene.children, everyone - below})] = gene

    def gene_of(node: dendropy.Node) -> GeneNode | None:
        below = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        if node.is_leaf():
            return rec.nodes.get(below)
        parts = {frozenset(leaf.taxon.label for leaf in c.leaf_iter()) for c in node.child_nodes()}
        if below != everyone:
            parts.add(everyone - below)
        return by_parts.get(frozenset(parts))

    rows: list[Row] = []
    ends: dict[Leaves, str] = {}  # each side of the root -> our node there
    for node in tree.preorder_node_iter():
        label = node.taxon.label if node.is_leaf() else node.label or ""
        match = gene_of(node)
        if match is None:
            raise ReconcileError(f"node {label!r} is not in the reconciled tree")
        rows.append(Row(label, node.is_leaf(), match.event, match.species, tuple(match.losses)))
        if match.below in rec.root:
            ends[match.below] = label
    return rows, [ends[side] for side in rec.root], reroot(tree, rec.root[0])


def reroot(tree: dendropy.Tree, side: Leaves) -> str:
    """The tree rooted at the middle of the branch separating `side` from the rest."""
    tree = tree.clone(depth=1)
    everyone = frozenset(leaf.taxon.label for leaf in tree.leaf_node_iter())
    target = None
    for node in tree.postorder_node_iter():
        below = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        if node is not tree.seed_node and below in (side, everyone - side):
            target = node
            break
    if target is None:
        raise ReconcileError("the reconciled root is not on a branch of the tree")
    length = target.edge.length or 0.0
    half = length / 2
    tree.reroot_at_edge(target.edge, length1=half, length2=half, update_bipartitions=True)
    tree.suppress_unifurcations()
    rooted: str = tree.as_string(
        schema="newick", suppress_rooting=True, unquoted_underscores=True
    ).strip()
    return rooted


def root_row(rec: Reconciliation) -> Row:
    """GeneRax's root: its event, species and the losses above it."""
    top = rec.nodes[rec.top]
    return Row(ROOT, False, top.event, top.species, tuple(top.losses))


def totals(rows: list[Row]) -> dict[str, int]:
    return {
        "duplications": sum(r.event == "D" for r in rows),
        "transfers": sum(r.event == "T" for r in rows),
        "speciations": sum(r.event == "S" for r in rows),
        "losses": sum(len(r.losses) for r in rows),
    }


WEAK_OVERLAP = 0.1  # species overlap below which a duplication is shown as weak (R5)


@dataclass(frozen=True)
class Overlap:
    """How many species a duplication's two copies share, of all the species under it."""

    shared: int
    total: int

    @property
    def ratio(self) -> float:
        return self.shared / self.total if self.total else 0.0

    @property
    def weak(self) -> bool:
        return self.ratio < WEAK_OVERLAP

    def describe(self) -> str:
        return f"its copies share {self.shared} of {self.total} species"


def overlaps(tsv: str, rooted_newick: str) -> dict[str, Overlap]:
    """{node (or ROOT): species overlap} for every duplication in a combination's
    reconciliation.tsv, its copies taken in GeneRax's rooting (reconciliation.nwk)."""
    species: dict[str, str] = {}
    events: dict[str, str] = {}
    for line in tsv.splitlines()[1:]:
        cells = line.split("\t")
        if len(cells) == 5:
            events[cells[0]] = cells[1]
            if cells[1] == "leaf":
                species[cells[0]] = cells[2]
    _, nodes = rooted_clades(rooted_newick)
    found = {}
    for _below, label, children in nodes:
        if events.get(label) != "duplication" or len(children) < 2:
            continue
        sets = [{species.get(gene, gene) for gene in child} for child in children]
        shared = set.intersection(*sets) if len(sets) == 2 else {
            s for i, a in enumerate(sets) for b in sets[i + 1 :] for s in a & b
        }  # fmt: skip
        found[label] = Overlap(len(shared), len(set.union(*sets)))
    return found


def leaf_species(tsv: str) -> dict[str, str]:
    """{sequence: its species' name} from a combination's reconciliation.tsv."""
    found = {}
    for line in tsv.splitlines()[1:]:
        cells = line.split("\t")
        if len(cells) == 5 and cells[1] == "leaf":
            found[cells[0]] = cells[2].replace("_", " ")
    return found


def weak_duplications(tsv: str, rooted_newick: str) -> int:
    return sum(o.weak for o in overlaps(tsv, rooted_newick).values())


Marks = tuple[dict[str, tuple[str, str, Overlap | None]], dict[Leaves, list[str]]]


def marks(tsv: str, rooted_newick: str) -> Marks:
    """What the tree viewer draws from a combination's published reconciliation.tsv and
    reconciliation.nwk: {node name: (event, species, species overlap)} for duplications
    (with their overlap) and transfers (None), and
    {unrooted split (`trees.canonical_split`): species that lost a copy on that branch}.
    Losses are recorded on the branch above a node in GeneRax's rooting; keyed by split,
    they find their branch whatever rooting the viewer draws (the two halves of the
    root branch become one)."""
    tree = load_tree(rooted_newick)
    everyone = frozenset(leaf.taxon.label for leaf in tree.leaf_node_iter())
    below: dict[str, Leaves] = {}
    for node in tree.preorder_node_iter():
        name = node.taxon.label if node.is_leaf() else node.label or ""
        if name:
            below[name] = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
    events: dict[str, tuple[str, str, Overlap | None]] = {}
    losses: dict[Leaves, list[str]] = {}
    overlap = overlaps(tsv, rooted_newick)
    for line in tsv.splitlines()[1:]:
        cells = line.split("\t")
        if len(cells) != 5 or cells[0] == ROOT:  # the root is on a branch: not drawn here
            continue
        node, event, species, _count, lost = cells
        if event in ("duplication", "transfer"):
            events[node] = (event, species, overlap.get(node))
        if lost and node in below:
            split = canonical_split(below[node], everyone)
            losses.setdefault(split, []).extend(s for s in lost.split("; ") if s)
    return events, losses


@dataclass(frozen=True)
class Paralogy:
    """What reconciliation says about one ancestor (`paralogy`)."""

    duplication: bool  # a duplication in its own tree: the gene before the copies split
    reconciled_node: str = ""  # the same clade in the reconciled tree, if it is there
    reconciled_event: str = ""  # and its event there (speciation, duplication, ...)
    mixes: str = ""  # else, the reconciled tree's duplication whose paralogs it mixes
    groups: tuple[tuple[str, ...], ...] = ()  # its sequences in each of those paralogs
    # each group's ancestor in the reconciled tree ("" for a single sequence)
    group_nodes: tuple[str, ...] = ()


def rooted_clades(newick: str) -> tuple[dict[str, Leaves], list[tuple[Leaves, str, list[Leaves]]]]:
    """A reconciliation.nwk's {label: leaves below} and every internal node as (leaves,
    label, children's leaves); the root is unlabelled (GeneRax's, on a branch) and is
    given ROOT's name."""
    tree = load_tree(newick)
    labels: dict[str, Leaves] = {}
    nodes = []
    for node in tree.preorder_internal_node_iter():
        below = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        label = ROOT if node is tree.seed_node else node.label or ""
        children = [frozenset(x.taxon.label for x in c.leaf_iter()) for c in node.child_nodes()]
        if label:
            labels[label] = below
        nodes.append((below, label, children))
    return labels, nodes


def event_names(tsv: str) -> dict[str, str]:
    """{node (or ROOT): event} from a reconciliation.tsv."""
    found = {}
    for line in tsv.splitlines()[1:]:
        cells = line.split("\t")
        if len(cells) == 5:
            found[cells[0]] = cells[1]
    return found


def paralogy(
    node: str, own: tuple[str, str], reconciled: tuple[str, str] | None
) -> Paralogy | None:
    """One ancestor of a combination, from that combination's (reconciliation.tsv,
    reconciliation.nwk) and, if the job has one, the reconciled tree's of the same
    alignment and model. Clades are taken in GeneRax's rooting of each tree. None if
    the node is not in its own reconciliation.

    The ancestor mixes paralogs when its clade is absent from the reconciled tree and
    the smallest reconciled clade containing it is a duplication: its sequences then
    fall in two (or more) paralog groups there, and the reconstruction averages across
    a duplication. A clade that is absent under a speciation is only a different
    topology, not a paralogy problem."""
    own_labels, _ = rooted_clades(own[1])
    if node not in own_labels:
        return None
    duplication = event_names(own[0]).get(node) == "duplication"
    if reconciled is None:
        return Paralogy(duplication)
    clade = own_labels[node]
    _, nodes = rooted_clades(reconciled[1])
    events = event_names(reconciled[0])
    containing = [(below, label, children) for below, label, children in nodes if clade <= below]
    if not containing:
        return Paralogy(duplication)
    below, label, children = min(containing, key=lambda n: len(n[0]))
    if below == clade:
        return Paralogy(duplication, reconciled_node=label, reconciled_event=events.get(label, ""))
    if events.get(label) != "duplication":
        return Paralogy(duplication)
    groups = tuple(sorted(tuple(sorted(clade & c)) for c in children if clade & c))

    def ancestor(group: tuple[str, ...]) -> str:
        if len(group) < 2:
            return ""
        holding = [(b, name) for b, name, _ in nodes if set(group) <= b and name != ROOT]
        return min(holding, key=lambda n: len(n[0]))[1] if holding else ""

    return Paralogy(
        duplication, mixes=label, groups=groups, group_nodes=tuple(ancestor(g) for g in groups)
    )
