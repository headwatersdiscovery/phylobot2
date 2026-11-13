"""A mutational path compared across the job's combinations.

A path read from one alignment can be wrong where alignments differ, and the changes
several alignments agree about are the ones to trust. So the path traced in one
combination (`phylobot.paths`) is followed through the others and every change is asked:
do they make it too?

- **The same path elsewhere.** Each combination has its own tree and node labels, so the
  path's ancestors are found by what they are ancestors of (`phylobot.analogs`: the clade
  with an outgroup, else the subtrees meeting at the node), and an extant sequence is
  itself. Between the first and last of them the other tree is walked as it is
  (`paths.between`), so it may pass ancestors the traced tree does not have.
- **The same branch.** A branch of the traced path is compared where the other
  combination has both of its ends; what is compared is the net change between those two
  ancestors there, whether or not its tree puts other nodes between them.
- **The same column.** The other alignment's columns are laid on the traced one's by the
  extant residues they share (`alignment_map.column_map`, a majority vote), as
  `evaluation.disagreement` does for the paper. A change is **the same** when the other
  combination makes exactly that call (before and after) at that column, **different**
  when it changes that column otherwise, and **absent** when it does not change it.

Agreement is not accuracy: every alignment can make the same mistake. What the counts
give is how far a change depends on the alignment and model one happened to choose.

Pure functions over sequences and Newick text; columns are 0-based.
"""

from __future__ import annotations

from collections.abc import Hashable, Mapping
from dataclasses import dataclass, field
from itertools import pairwise
from typing import Any

from phylobot import alignment_map, paths

GAP = "-"
Call = tuple[str, str]  # (before, after), each a residue or "-"

SAME, DIFFERENT, ABSENT = "same", "different", "absent"
# How widely a change is shared, strongest first; NA where nothing could be compared.
ALL, MOST, FEW, NONE, NA = "all", "most", "few", "none", "na"
LEVELS = (ALL, MOST, FEW, NONE, NA)
LEVEL_LABELS = {
    ALL: "in every one",
    MOST: "in most",
    FEW: "in few",
    NONE: "only here",
    NA: "not compared",
}


@dataclass(frozen=True)
class Side:
    """One combination's reconstruction, as far as following a path through it needs."""

    key: str
    aligner: str
    alignment: dict[str, str]  # the extant rows
    ancestors: dict[str, str]  # node -> its row, in the chosen gap version
    newick: str
    places: Mapping[str, Any] = field(default_factory=dict)  # node -> analogs.Where

    def row(self, name: str) -> str | None:
        found = self.ancestors.get(name)
        if found is None:
            found = self.alignment.get(name)
        return None if found is None else found.upper()


def _key(place: Any, by_clade: bool) -> Hashable:
    found: Hashable = place.below if by_clade else place.parts
    return found


def counterparts(main: Side, other: Side, route: list[str], by_clade: bool) -> dict[str, str]:
    """{name on the traced path: the same ancestor (or sequence) in `other`}, for those
    it has: an ancestor by its place in the job's ancestor map, a sequence by its name."""
    lookup = {_key(place, by_clade): node for node, place in other.places.items()}
    found = {}
    for name in route:
        if name in main.alignment:
            if name in other.alignment:
                found[name] = name
            continue
        place = main.places.get(name)
        twin = lookup.get(_key(place, by_clade)) if place is not None else None
        if twin is not None and twin in other.ancestors:
            found[name] = twin
    return found


def route_through(other: Side, matched: list[str]) -> list[str]:
    """`other`'s own path from the first matched node to the last, through whatever its
    tree puts between them; nodes without a reconstruction (the root an outgroup adds)
    are left out. With a single match, that node alone."""
    if not matched:
        return []
    if len(matched) == 1:
        return list(matched)
    walked = paths.between(other.newick, matched[0], matched[-1])
    return [n for n in walked if other.row(n) is not None] or list(matched)


def differences(before: str, after: str) -> dict[int, Call]:
    """{column: (before, after)} where two rows of one alignment differ."""
    return {c: (a, b) for c, (a, b) in enumerate(zip(before, after, strict=True)) if a != b}


def projection(other: dict[str, str], main: dict[str, str]) -> dict[int, int] | None:
    """Each column of `other` -> the column of `main` holding most of its residues, through
    the sequences both alignments hold with the same residues. None where they share too
    few to place anything (another sequence set: nothing there can be compared)."""
    if other == main:
        width = len(next(iter(main.values()), ""))
        return {c: c for c in range(width)}
    shared = [
        n for n in other if n in main and other[n].replace(GAP, "") == main[n].replace(GAP, "")
    ]
    if len(shared) < 2 or 2 * len(shared) < len(main):
        return None
    return alignment_map.column_map({n: other[n] for n in shared}, {n: main[n] for n in shared})


def sync_columns(columns: dict[int, int], width: int, main_width: int) -> list[int]:
    """For scrolling two alignments together: per column of the other alignment, the
    traced alignment's column to hold beside it. `projection`'s votes can step backwards
    and skip columns without residues, so this is the running highest, never past the end."""
    out = []
    high = 0
    for column in range(width):
        high = max(high, columns.get(column, high))
        out.append(min(high, max(main_width - 1, 0)))
    return out


def sync_on(other_row: str, main_row: str) -> list[int] | None:
    """`sync_columns` pinned to one extant sequence instead: each column of the other
    alignment beside the traced alignment's column that holds the same residue of that
    sequence (the next one, where the column is a gap in it). None unless both rows are
    the same sequence."""
    if other_row.replace(GAP, "").upper() != main_row.replace(GAP, "").upper():
        return None
    held = alignment_map.residue_columns(main_row)
    if not held:
        return None
    out = []
    seen = 0
    for character in other_row:
        out.append(held[min(seen, len(held) - 1)])
        if character != GAP:
            seen += 1
    return out


@dataclass(frozen=True)
class Verdict:
    """What one other combination says about a change of the traced path."""

    side: str
    state: str  # SAME, DIFFERENT or ABSENT
    call: Call | None = None  # its own change at that column, when DIFFERENT


@dataclass(frozen=True)
class Support:
    """A change's verdicts, one per combination that could be compared on its branch."""

    verdicts: tuple[Verdict, ...] = ()

    @property
    def compared(self) -> int:
        return len(self.verdicts)

    @property
    def agreeing(self) -> int:
        return sum(1 for v in self.verdicts if v.state == SAME)

    @property
    def fraction(self) -> float | None:
        return self.agreeing / self.compared if self.compared else None

    @property
    def level(self) -> str:
        n, k = self.compared, self.agreeing
        if not n:
            return NA
        if k == n:
            return ALL
        if k == 0:
            return NONE
        return MOST if 2 * k >= n else FEW

    def sides(self, state: str) -> list[str]:
        return [v.side for v in self.verdicts if v.state == state]


@dataclass(frozen=True)
class BranchMatch:
    """One branch of the traced path in one other combination."""

    side: str
    parent: str = ""  # the branch's ends there ("" where its tree has no such ancestor)
    child: str = ""
    comparable: bool = False  # both ends there, and columns that can be laid on ours
    direct: bool = False  # its tree has nothing between the two
    same: int = 0  # of the traced branch's changes
    different: int = 0
    absent: int = 0
    extra: int = 0  # changes it makes on this branch at columns the traced path leaves alone

    @property
    def total(self) -> int:
        return self.same + self.different + self.absent


@dataclass
class Followed:
    """The traced path as one other combination has it."""

    side: Side
    twins: dict[str, str]  # traced name -> its name here
    route: list[str]  # its own path, oldest first
    columns: dict[int, int] | None  # its columns -> the traced alignment's (None: no map)
    # Its changed cells that are the traced path's change too:
    # {its node: {its column: (the traced branch, the traced column)}}.
    shared: dict[str, dict[int, tuple[int, int]]] = field(default_factory=dict)
    missing: list[str] = field(default_factory=list)  # traced ancestors its tree lacks


@dataclass
class Comparison:
    changes: list[dict[int, Call]]  # per traced branch: its changes by column
    support: list[dict[int, Support]]  # per traced branch: each change's support
    branches: list[list[BranchMatch]]  # per traced branch: one per other combination
    followed: list[Followed]

    def levels(self, branch: int) -> dict[str, int]:
        """How many of a branch's changes sit at each level, in LEVELS order."""
        counts = dict.fromkeys(LEVELS, 0)
        for found in self.support[branch].values():
            counts[found.level] += 1
        return counts


def compare(main: Side, route: list[str], others: list[Side], by_clade: bool) -> Comparison:
    """Follow the traced path (`route`, in `main`) through each of `others` and judge
    every change of it. A combination whose alignment cannot be laid on the traced one,
    or whose tree lacks an end of a branch, gives no verdict there: a change is only
    counted against combinations that could have made it."""
    rows = [main.row(name) or "" for name in route]
    changes = [differences(a, b) for a, b in pairwise(rows)]
    verdicts: list[dict[int, list[Verdict]]] = [{c: [] for c in found} for found in changes]
    branches: list[list[BranchMatch]] = [[] for _ in changes]
    followed = []
    for other in others:
        columns = projection(other.alignment, main.alignment)
        twins = counterparts(main, other, route, by_clade)
        theirs = route_through(other, [twins[n] for n in route if n in twins])
        position = {name: i for i, name in enumerate(theirs)}
        track = Followed(
            other, twins, theirs, columns,
            missing=[n for n in route if n not in twins],
        )  # fmt: skip
        for i, (parent, child) in enumerate(pairwise(route)):
            a, b = twins.get(parent, ""), twins.get(child, "")
            if columns is None or not a or not b:
                branches[i].append(BranchMatch(other.key, a, b))
                continue
            own = differences(other.row(a) or "", other.row(b) or "")
            placed: dict[int, tuple[Call, int]] = {}  # traced column -> its call, its column
            for column, call in own.items():
                if column in columns:
                    placed[columns[column]] = (call, column)
            counts = dict.fromkeys((SAME, DIFFERENT, ABSENT), 0)
            agreed: dict[int, int] = {}  # its column -> the traced column
            for column, call in changes[i].items():
                theirs_here = placed.get(column)
                if theirs_here is None:
                    state, shown = ABSENT, None
                elif theirs_here[0] == call:
                    state, shown = SAME, None
                    agreed[theirs_here[1]] = column
                else:
                    state, shown = DIFFERENT, theirs_here[0]
                counts[state] += 1
                verdicts[i][column].append(Verdict(other.key, state, shown))
            # Its cells that carry an agreed change: every node of its own path from just
            # below the branch's top to its bottom that differs there from the row above.
            lo, hi = position.get(a, -1), position.get(b, -1)
            if agreed and 0 <= lo < hi:
                for upper, lower in pairwise(theirs[lo : hi + 1]):
                    up, low = other.row(upper) or "", other.row(lower) or ""
                    for c, column in agreed.items():
                        if up[c] != low[c]:
                            track.shared.setdefault(lower, {})[c] = (i, column)
            branches[i].append(
                BranchMatch(
                    other.key, a, b, True, abs(hi - lo) == 1,
                    counts[SAME], counts[DIFFERENT], counts[ABSENT],
                    sum(1 for column in placed if column not in changes[i]),
                )
            )  # fmt: skip
        followed.append(track)
    support = [
        {column: Support(tuple(found)) for column, found in branch.items()} for branch in verdicts
    ]
    return Comparison(changes, support, branches, followed)


def merged_rows(route: list[str], followed: list[Followed]) -> list[list[str]]:
    """Rows for showing the paths side by side, each ancestor beside itself: one list per
    row, the traced path's name first, then each followed combination's ("" where it has
    no such ancestor). A node only another combination's tree has gets a row of its own,
    between the rows of its matched neighbours; several combinations' share those rows."""
    rows: list[list[str]] = [[name] + [""] * len(followed) for name in route]
    index = {name: i for i, name in enumerate(route)}  # the traced rows, kept up to date
    for slot, track in enumerate(followed, start=1):
        back = {twin: name for name, twin in track.twins.items()}
        at = -1  # the row of the last node placed
        for node in track.route:
            name = back.get(node)
            if name is not None and index[name] > at:
                at = index[name]
                rows[at][slot] = node
                continue
            # Unmatched here: the next row below, if it is an extra one still free in this
            # slot; else a new row there.
            below = at + 1
            if below < len(rows) and not rows[below][0] and not rows[below][slot]:
                at = below
            else:
                rows.insert(below, [""] * (len(followed) + 1))
                index = {n: (i if i < below else i + 1) for n, i in index.items()}
                at = below
            rows[at][slot] = node
    return rows
