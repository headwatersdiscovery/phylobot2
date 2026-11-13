"""The same ancestor under the job's other alignments, trees and models.

Each combination reconstructs its own tree and names its nodes in its own way, so "Node7"
in one tree says nothing about "Node7" in another. A node is identified instead by the
extant sequences around it, which every tree of a job shares:

- With an outgroup, by its **clade**: the sequences below it on the outgroup-rooted tree
  (`tree.rooted.nwk`). This is how ancestral reconstruction names ancestors ("the
  ancestor of the fungal Ime2s"), whatever the branching inside the clade. A match whose
  clade agrees but whose subtrees differ is flagged (`Analog.alike` False): the same
  ancestor by name, reached through a different history.
- Without one, the trees are unrooted and "below" means nothing, so by the subtrees that
  meet at it (`gaps.node_partitions`): stricter, the only definition that does not depend
  on where each file happens to put its root.

Pure functions over Newick text. The pipeline's ANCESTOR_MAP step turns every
combination's tree into rows of content-derived IDs (`map_rows`, `ancestor_map.tsv`), and
the viewers match on those (`placed`) rather than parse every tree per request (2.7 s
on a 285-combination job). The IDs hash the sequence names, never positional labels, so
they hold across jobs ("Duplicate and run again" with the input reordered).
Identity between matched ancestors is measured on the comparison the ancestor viewer
draws (`identity`), so what is counted is what is shown.
"""

from __future__ import annotations

import hashlib
from collections.abc import Hashable, Iterable, Mapping
from dataclasses import dataclass
from typing import Protocol

from phylobot.trees import load_tree

GAP = "-"


@dataclass(frozen=True)
class Place:
    """Where an internal node sits, in terms every tree of the job shares."""

    below: frozenset[str]  # the extant sequences below it, as the file roots the tree
    parts: frozenset[frozenset[str]]  # the subtrees (leaf sets) that meet at it


def places(newick: str) -> dict[str, Place]:
    """{label: place} for every labelled internal node. Unlabelled ones (the root that
    outgroup rooting adds by splitting a branch) have no reconstruction to match."""
    tree = load_tree(newick)
    everyone = frozenset(leaf.taxon.label for leaf in tree.leaf_node_iter())
    found: dict[str, Place] = {}
    for node in tree.postorder_internal_node_iter():
        if not node.label:
            continue
        below = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        parts = {frozenset(leaf.taxon.label for leaf in c.leaf_iter()) for c in node.child_nodes()}
        if below != everyone:
            parts.add(everyone - below)
        found[node.label] = Place(below, frozenset(parts))
    return found


ID_LENGTH = 16  # hex digits: 64 bits, so no two clades of a job will share one
MAP = "ancestor_map.tsv"
MAP_FIELDS = ("combination", "node", "clade", "parts", "leaves")


def clade_id(labels: Iterable[str]) -> str:
    """A set of sequences, by their names."""
    text = "\n".join(sorted(labels))
    return hashlib.sha256(text.encode()).hexdigest()[:ID_LENGTH]


def parts_id(parts: Iterable[Iterable[str]]) -> str:
    """The subtrees that meet at a node, each by its `clade_id`."""
    text = "\n".join(sorted(clade_id(part) for part in parts))
    return hashlib.sha256(text.encode()).hexdigest()[:ID_LENGTH]


@dataclass(frozen=True)
class Located:
    """A node of `ancestor_map.tsv`: a Place by its IDs (`index` compares them alike)."""

    below: str  # clade_id of the sequences below it; "" without an outgroup-rooted tree
    parts: str  # parts_id
    leaves: int  # how many sequences are below it


class Where(Protocol):
    @property
    def below(self) -> Hashable: ...

    @property
    def parts(self) -> Hashable: ...


@dataclass(frozen=True)
class MapRow:
    combination: str
    node: str
    place: Located


def map_rows(combination: str, newick: str, rooted: bool) -> list[MapRow]:
    """A combination's nodes for `ancestor_map.tsv`: from its outgroup-rooted tree when
    `rooted` (then with clades), else from its tree (subtrees only: "below" would follow
    wherever the file happens to put its root)."""
    rows = []
    for node, place in places(newick).items():
        clade = clade_id(place.below) if rooted else ""
        located = Located(clade, parts_id(place.parts), len(place.below))
        rows.append(MapRow(combination, node, located))
    return rows


def write_map(rows: Iterable[MapRow]) -> str:
    lines = ["\t".join(MAP_FIELDS)]
    for r in sorted(rows, key=lambda r: (r.combination, r.place.leaves, r.node)):
        lines.append(
            f"{r.combination}\t{r.node}\t{r.place.below}\t{r.place.parts}\t{r.place.leaves}"
        )
    return "\n".join(lines) + "\n"


def read_map(text: str) -> list[MapRow]:
    lines = [line for line in text.splitlines() if line]
    if not lines or tuple(lines[0].split("\t")) != MAP_FIELDS:
        raise ValueError(f"not an {MAP}")
    rows = []
    for line in lines[1:]:
        combination, node, clade, parts, leaves = line.split("\t")
        rows.append(MapRow(combination, node, Located(clade, parts, int(leaves))))
    return rows


def placed(rows: Iterable[MapRow], by_clade: bool) -> dict[str, dict[str, Located]]:
    """{combination: {node: its IDs}}, as `index` takes them. Matching by clade leaves out
    combinations without an outgroup-rooted tree, as matching on the files did."""
    out: dict[str, dict[str, Located]] = {}
    for r in rows:
        if by_clade and not r.place.below:
            continue
        out.setdefault(r.combination, {})[r.node] = r.place
    return out


@dataclass(frozen=True)
class Analog:
    """A node of another combination's tree that is the same ancestor."""

    combination: str
    node: str
    alike: bool  # the same subtrees meet there too, not only the same clade


def _key(place: Where, by_clade: bool) -> Hashable:
    return place.below if by_clade else place.parts


def index(
    trees: Mapping[str, Mapping[str, Where]], combination: str, by_clade: bool
) -> dict[str, list[Analog]]:
    """{node of `combination`: its analogs in every other combination, in `trees` order}.
    `trees` is {combination: places}; a combination without the node's clade (or, without
    an outgroup, its partition) simply contributes nothing for it."""
    lookup = {
        other: {_key(place, by_clade): (label, place) for label, place in mine.items()}
        for other, mine in trees.items()
        if other != combination
    }
    result: dict[str, list[Analog]] = {}
    for node, place in trees.get(combination, {}).items():
        found = []
        for other, table in lookup.items():
            hit = table.get(_key(place, by_clade))
            if hit is not None:
                found.append(Analog(other, hit[0], hit[1].parts == place.parts))
        result[node] = found
    return result


def find(
    trees: Mapping[str, Mapping[str, Where]], combination: str, node: str, by_clade: bool
) -> list[Analog]:
    """One node's analogs (`index`, for a single node)."""
    return index(trees, combination, by_clade).get(node, [])


def missing(
    trees: Mapping[str, Mapping[str, Where]], combination: str, found: Iterable[Analog]
) -> list[str]:
    """The other combinations whose trees lack the node, in `trees` order."""
    have = {a.combination for a in found}
    return [c for c in trees if c != combination and c not in have]


@dataclass(frozen=True)
class Identity:
    """Two ancestors compared on one frame of columns."""

    matches: int  # columns where both have the same residue
    length: int  # the longer ancestor's residue count

    @property
    def fraction(self) -> float:
        return self.matches / self.length if self.length else 0.0

    @property
    def percent(self) -> float:
        return round(100 * self.fraction, 1)


def identity(a: str, b: str, a_length: int, b_length: int) -> Identity:
    """Identity of two ancestors drawn on the same columns (`a` and `b`, gaps as "-"):
    matching residues over the longer ancestor's length, so residues one has and the
    other lacks lower it, as `evaluation.gaps.aligned_identity` does. The lengths are the
    ancestors' own, since a frame can leave some residues of a mapped ancestor undrawn
    (two of them meeting in one column), and those must still count against it."""
    matches = sum(1 for x, y in zip(a, b, strict=True) if x != GAP and x == y)
    return Identity(matches, max(a_length, b_length))
