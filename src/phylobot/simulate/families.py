"""Gene families: duplication and loss along a species tree.

AliSim cannot duplicate or lose genes, so a family spec (`[family]`) is grown here, in
pure Python, before AliSim runs: a seeded birth-death process of gene copies along the
rooted species tree (the spec's `[tree]`). Every copy duplicates at rate `duplication`
and is lost at rate `loss`, per unit of species-tree branch length; at each species node
every copy speciates into one copy per child. Branch lengths carry over unchanged, so
the gene tree's lengths are in the species tree's units (expected substitutions per
site) and AliSim evolves sequences along it as along any tree.

The process starts with one copy at the top of an optional stem (`origin`, a branch
above the species root), so a duplication can predate the first speciation. Chosen
events build a scenario deliberately ("hidden paralogy"): `duplications` force a
duplication on a species branch, named by the species below it, at a fraction `at` of
its length, and `losses` remove a copy made by one of them (`copy = "1b"`: copy b of the
first forced duplication, with everything descended from it) on a branch at or below it.
A forced duplication applies to every copy on its branch at that point.

The family is redrawn, continuing the same random stream, until at least `min_genes`
copies survive (conditioning on survival, as birth-death simulators do), up to
`ATTEMPTS` times. Lost copies are pruned and the nodes left with one child are
removed, so a duplication or speciation can be real but invisible in the gene tree;
true_events.tsv lists every event with the gene node it became, or "-".

Files written next to AliSim's inputs and copied into the truth bundle by `collect`:

    tree.nwk           the gene tree AliSim runs on (internal nodes Node1, Node2, ...,
                       in preorder; leaves <species>_<n>); AliSim keeps these labels,
                       so true_ancestors.fasta is keyed by them
    species_tree.nwk   the rooted species tree, internal nodes S1, S2, ... (preorder)
    species.tsv        gene -> species
    true_events.tsv    every duplication, loss and speciation, in preorder
    family.json        counts, attempts and the stem length
"""

import json
import math
import random
from dataclasses import dataclass, field
from pathlib import Path

import dendropy
from dendropy.utility.error import DataParseError

from phylobot.simulate.spec import NAME_PATTERN, Family

ATTEMPTS = 1000
MAX_NODES = 50_000  # gene-tree nodes per attempt, lost copies included
SPECIES_TREE = "species_tree.nwk"
SPECIES = "species.tsv"
EVENTS = "true_events.tsv"
SUMMARY = "family.json"
FILES = (SPECIES_TREE, SPECIES, EVENTS, SUMMARY)


class FamilyError(ValueError):
    """A family that cannot be grown on this species tree, with the reason."""


@dataclass
class Species:
    """A species-tree node and the branch above it (for the root, the stem)."""

    name: str
    top: float  # depth of the branch's top, from the start of the stem
    length: float
    children: list[Species] = field(default_factory=list)
    leaves: frozenset[str] = frozenset()


@dataclass
class Forced:
    kind: str  # "duplication" | "loss"
    position: float  # along its branch
    index: int  # 1-based, in the spec's list
    copy: str = ""  # losses: e.g. "1b"


@dataclass
class Gene:
    event: str  # "origin" | "duplication" | "loss" | "speciation" | "leaf"
    branch: str  # the species branch it lies on (for a speciation, the one ending there)
    depth: float
    forced: bool = False
    copies: tuple[str, ...] = ()  # forced-duplication copies it descends from
    children: list[Gene] = field(default_factory=list)
    name: str = ""  # in the pruned gene tree; "" when pruned or collapsed


@dataclass
class Grown:
    origin: Gene
    root: Gene  # the gene tree's root: the survivors' common ancestor
    species: Species
    attempts: int
    family: Family


def species_tree(newick: str) -> Species:
    """The rooted species tree from prepare_tree's text, internal nodes named S1, S2, ...
    in preorder; the root's branch is the stem, added by `grow`."""
    try:
        tree = dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)
    except DataParseError as error:
        raise FamilyError(f"the species tree could not be read: {error}") from error
    made: dict[int, Species] = {}
    seen: set[str] = set()
    internal = 0
    for node in tree.preorder_node_iter():
        if node.taxon is not None:
            name = node.taxon.label
            if not NAME_PATTERN.fullmatch(name) or name in seen:
                raise FamilyError(f"species {name!r}: unique names of letters, digits, _ . | -")
            seen.add(name)
        else:
            internal += 1
            name = f"S{internal}"
        parent = made.get(id(node.parent_node)) if node.parent_node else None
        length = 0.0
        if parent is not None:
            if node.edge.length is None or node.edge.length < 0:
                raise FamilyError(f"species branch {name}: needs a length of 0 or more")
            length = float(node.edge.length)
        here = Species(name, parent.top + parent.length if parent else 0.0, length)
        if parent is not None:
            parent.children.append(here)
        made[id(node)] = here
    root = made[id(tree.seed_node)]
    if any(f"S{i}" in seen for i in range(1, internal + 1)):
        raise FamilyError("species names S1, S2, ... are taken by the tree's internal nodes")
    if len(seen) < 2:
        raise FamilyError("a species tree needs at least two species")
    for node in postorder(root):
        node.leaves = (
            frozenset().union(*(c.leaves for c in node.children))
            if node.children
            else frozenset({node.name})
        )
    return root


def postorder[T: (Species, Gene)](root: T) -> list[T]:
    """Children before parents, without recursion (trees can be thousands deep)."""
    order, stack = [], [root]
    while stack:
        node = stack.pop()
        order.append(node)
        stack.extend(node.children)
    return order[::-1]


def preorder[T: (Species, Gene)](root: T) -> list[T]:
    order, stack = [], [root]
    while stack:
        node = stack.pop()
        order.append(node)
        stack.extend(reversed(node.children))
    return order


def forced_events(family: Family, root: Species) -> dict[str, list[Forced]]:
    """The spec's chosen events by species branch, in order along it; checks they exist."""
    by_clade = {node.leaves: node for node in preorder(root)}
    branch_of: dict[int, Species] = {}
    events: dict[str, list[Forced]] = {}

    def place(clade: tuple[str, ...], at: float, what: str) -> Species:
        node = by_clade.get(frozenset(clade))
        if node is None:
            raise FamilyError(f"{what}: no species branch has exactly {', '.join(clade)} below it")
        if node.length <= 0:
            raise FamilyError(f"{what}: its branch has no length (the stem needs origin > 0)")
        return node

    for i, d in enumerate(family.duplications, 1):
        node = place(d.clade, d.at, f"forced duplication {i}")
        branch_of[i] = node
        events.setdefault(node.name, []).append(Forced("duplication", d.at * node.length, i))
    for i, loss in enumerate(family.losses, 1):
        node = place(loss.clade, loss.at, f"forced loss {i}")
        k = int(loss.copy[:-1])
        if k not in branch_of:
            raise FamilyError(f"forced loss {i}: there is no forced duplication {k}")
        source = branch_of[k]
        position = loss.at * node.length
        dup_position = next(e.position for e in events[source.name] if e.index == k)
        if not (node.leaves < source.leaves or (node is source and position > dup_position)):
            raise FamilyError(
                f"forced loss {i}: must come after duplication {k}, on its branch or below it"
            )
        events.setdefault(node.name, []).append(Forced("loss", position, i, loss.copy))
    for listed in events.values():
        listed.sort(key=lambda e: (e.position, e.kind == "loss", e.index))
    return events


def _attempt(
    family: Family, root: Species, forced: dict[str, list[Forced]], rng: random.Random
) -> Gene:
    """One draw of the process; the full gene tree, lost copies included."""
    origin = Gene("origin", root.name, 0.0)
    rate = family.duplication + family.loss
    count = 0
    # (parent gene, species branch, position along it, copies, next forced event's index)
    stack: list[tuple[Gene, Species, float, tuple[str, ...], int]] = [(origin, root, 0.0, (), 0)]
    while stack:
        parent, species, position, copies, i = stack.pop()
        events = forced.get(species.name, [])
        count += 1
        if count > MAX_NODES:
            raise FamilyError(
                f"the family grew past {MAX_NODES} gene copies: lower the duplication rate"
            )
        while True:
            wait = rng.expovariate(rate) if rate > 0 else math.inf
            t = position + wait
            if i < len(events) and events[i].position <= min(t, species.length):
                event = events[i]
                i += 1
                depth = species.top + event.position
                if event.kind == "duplication":
                    node = Gene("duplication", species.name, depth, True, copies)
                    parent.children.append(node)
                    for side in "ba":  # popped a first
                        tag = (*copies, f"{event.index}{side}")
                        stack.append((node, species, event.position, tag, i))
                    break
                if event.copy in copies:
                    parent.children.append(Gene("loss", species.name, depth, True, copies))
                    break
                position = event.position  # not this copy's: carry on (memoryless)
                continue
            if t < species.length:
                depth = species.top + t
                if rng.random() < family.duplication / rate:
                    node = Gene("duplication", species.name, depth, False, copies)
                    parent.children.append(node)
                    stack.append((node, species, t, copies, i))
                    stack.append((node, species, t, copies, i))
                else:
                    parent.children.append(Gene("loss", species.name, depth, False, copies))
                break
            end = species.top + species.length
            if not species.children:
                parent.children.append(Gene("leaf", species.name, end, False, copies))
            else:
                node = Gene("speciation", species.name, end, False, copies)
                parent.children.append(node)
                for child in reversed(species.children):
                    stack.append((node, child, 0.0, copies, 0))
            break
    return origin


def _survivors(origin: Gene) -> dict[int, int]:
    """Surviving leaves below each gene node, by id."""
    alive: dict[int, int] = {}
    for node in postorder(origin):
        if node.event == "leaf":
            alive[id(node)] = 1
        else:
            alive[id(node)] = sum(alive[id(c)] for c in node.children)
    return alive


def grow(family: Family, newick: str, seed: int) -> Grown:
    """Grow a family on a species tree (prepare_tree's text), seeded; names the pruned
    gene tree's nodes."""
    root = species_tree(newick)
    root.length = family.origin
    for child in root.children:
        _shift(child, family.origin)
    forced = forced_events(family, root)
    rng = random.Random(seed)
    attempts = 0
    while True:
        attempts += 1
        origin = _attempt(family, root, forced, rng)
        alive = _survivors(origin)
        if alive[id(origin)] >= family.min_genes:
            break
        if attempts == ATTEMPTS:
            raise FamilyError(
                f"fewer than {family.min_genes} genes survived in {ATTEMPTS} attempts: "
                "lower the loss rate"
            )
    top = origin
    while sum(1 for c in top.children if alive[id(c)]) < 2:
        top = next(c for c in top.children if alive[id(c)])
    _name(top, alive)
    return Grown(origin, top, root, attempts, family)


def _shift(node: Species, by: float) -> None:
    for n in preorder(node):
        n.top += by


def _kept(node: Gene, alive: dict[int, int]) -> list[Gene]:
    """The node's children in the pruned tree: each surviving child, or the node its
    single-child chain leads to."""
    kept = []
    for child in node.children:
        if not alive[id(child)]:
            continue
        while child.event != "leaf" and sum(1 for c in child.children if alive[id(c)]) == 1:
            child = next(c for c in child.children if alive[id(c)])
        kept.append(child)
    return kept


def _name(root: Gene, alive: dict[int, int]) -> None:
    internal = 0
    per_species: dict[str, int] = {}
    stack = [root]
    while stack:
        node = stack.pop()
        if node.event == "leaf":
            per_species[node.branch] = per_species.get(node.branch, 0) + 1
            node.name = f"{node.branch}_{per_species[node.branch]}"
            continue
        internal += 1
        node.name = f"Node{internal}"
        stack.extend(reversed(_kept(node, alive)))


def gene_newick(grown: Grown) -> str:
    alive = _survivors(grown.origin)
    text: dict[int, str] = {}
    for node in postorder(grown.root):
        if not node.name:
            continue
        if node.event == "leaf":
            text[id(node)] = node.name
            continue
        parts = [f"{text[id(c)]}:{_length(c.depth - node.depth)}" for c in _kept(node, alive)]
        text[id(node)] = f"({','.join(parts)}){node.name}"
    return text[id(grown.root)] + ";"


def species_newick(root: Species) -> str:
    text: dict[int, str] = {}
    for node in postorder(root):
        label = node.name
        if node.children:
            parts = [f"{text[id(c)]}:{_length(c.length)}" for c in node.children]
            label = f"({','.join(parts)}){node.name}"
        text[id(node)] = label
    return text[id(root)] + ";"


def _length(x: float) -> str:
    return f"{x:.10g}"


def events_table(grown: Grown) -> str:
    rows = ["event\tgene_node\tspecies_branch\tdepth\tforced\tcopies"]
    for node in preorder(grown.origin):
        if node.event in ("origin", "leaf"):
            continue
        rows.append(
            "\t".join(
                [
                    node.event,
                    node.name or "-",
                    node.branch,
                    _length(node.depth),
                    "yes" if node.forced else "no",
                    ",".join(node.copies) or "-",
                ]
            )
        )
    return "\n".join(rows) + "\n"


def summary(grown: Grown) -> dict[str, object]:
    nodes = preorder(grown.origin)
    leaves = [n for n in preorder(grown.root) if n.event == "leaf" and n.name]

    def count(event: str, visible: bool | None = None) -> int:
        return sum(
            1 for n in nodes if n.event == event and (visible is None or bool(n.name) == visible)
        )

    per_species: dict[str, int] = {}
    for leaf in leaves:
        per_species[leaf.branch] = per_species.get(leaf.branch, 0) + 1
    return {
        "attempts": grown.attempts,
        "origin": grown.family.origin,
        "genes": len(leaves),
        "species": len(grown.species.leaves),
        "species_without_genes": len(grown.species.leaves - per_species.keys()),
        "duplications": count("duplication"),
        "duplications_visible": count("duplication", visible=True),
        "losses": count("loss"),
        "speciations_visible": count("speciation", visible=True),
    }


def write(grown: Grown, work: Path) -> None:
    """The gene tree for AliSim and the family's truth files, into `work`."""
    (work / "tree.nwk").write_text(gene_newick(grown) + "\n")
    (work / SPECIES_TREE).write_text(species_newick(grown.species) + "\n")
    leaves = [n for n in preorder(grown.root) if n.event == "leaf" and n.name]
    (work / SPECIES).write_text(
        "gene\tspecies\n" + "".join(f"{n.name}\t{n.branch}\n" for n in leaves)
    )
    (work / EVENTS).write_text(events_table(grown))
    (work / SUMMARY).write_text(json.dumps(summary(grown), indent=2) + "\n")
