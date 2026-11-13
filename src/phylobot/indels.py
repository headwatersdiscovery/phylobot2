"""Indel lengths in an alignment, with a tree or without one.

**With a tree** (`on_tree`): Fitch presence at every internal node (`gaps.fitch_nodes`);
an event is a run of adjacent columns gained or lost on one branch, skipping columns absent
at both ends of it. On simulated data it reads MAFFT alignments with little bias, but
Fitch's tie rule turns deletions into insertions, so the two are not told apart here.

**End gaps are counted apart** (`Event.terminal`): a run reaching either end of what the
branch's two sequences hold. Most are truncated sequences or extensions, not indels, and
they are long, so folding them in would inflate the mean indel length several-fold.

**Without a tree** (`in_alignment`): each distinct interior gap (start, end) once. Only for
alignments with no tree yet: overlapping events read as one long gap, so these lengths are
far too long (PINK1 mean 29.6 against 3.9 on the tree). Never compare the two.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from phylobot.gaps import GAP, fitch_nodes
from phylobot.trees import load_tree

GAPS = re.compile(f"{GAP}+")


@dataclass(frozen=True)
class Event:
    length: int
    insertion: bool  # by Fitch's states, whose ties favour insertions: not to be trusted
    terminal: bool  # reaches an end of what the branch's two sequences hold


@dataclass(frozen=True)
class Lengths:
    """The interior indels' lengths and how many end gaps were left out."""

    counts: dict[int, int]  # length -> interior events
    end_gaps: int
    exposure: float = 0.0  # sum of branch length x the parent's residues; 0 without a tree

    @property
    def on_tree(self) -> bool:
        return self.exposure > 0

    @property
    def events(self) -> int:
        return sum(self.counts.values())

    @property
    def mean(self) -> float:
        return sum(k * n for k, n in self.counts.items()) / self.events if self.counts else 0.0

    @property
    def longest(self) -> int:
        return max(self.counts, default=0)

    @property
    def median(self) -> int:
        seen = 0
        for length in sorted(self.counts):
            seen += self.counts[length]
            if 2 * seen >= self.events:
                return length
        return 0

    @property
    def rate(self) -> float:
        """Interior events per substitution; 0 without a tree."""
        return self.events / self.exposure if self.exposure else 0.0

    def binned(self, last: int) -> list[int]:
        """Events of length 1..`last`, then one bin for everything longer."""
        bins = [self.counts.get(k, 0) for k in range(1, last + 1)]
        return [*bins, sum(n for k, n in self.counts.items() if k > last)]


def branch_events(parent: Sequence[bool], child: Sequence[bool]) -> list[Event]:
    """The runs of residues gained or lost between two rows of presence."""
    columns = [(p, c) for p, c in zip(parent, child, strict=True) if p or c]
    out: list[Event] = []
    start = 0
    for i, (p, c) in enumerate(columns):
        same_run = i > 0 and p != c and columns[i - 1] == (p, c)
        if not same_run:
            start = i
        ends_run = i + 1 == len(columns) or columns[i + 1] != (p, c)
        if p != c and ends_run:
            terminal = start == 0 or i + 1 == len(columns)
            out.append(Event(i - start + 1, insertion=c, terminal=terminal))
    return out


def tree_events(newick: str, alignment: Mapping[str, str]) -> tuple[list[Event], float]:
    """Every event on the tree, and the exposure (sum of branch length x the parent's
    residues). KeyError if a tip of the tree is not in the alignment."""
    tree = load_tree(newick)
    states = fitch_nodes(tree, dict(alignment))
    events: list[Event] = []
    exposure = 0.0
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        parent = states[id(node.parent_node)]
        child = (
            [ch != GAP for ch in alignment[node.taxon.label]]
            if node.is_leaf()
            else states[id(node)]
        )
        events += branch_events(parent, child)
        exposure += (node.edge.length or 0.0) * sum(parent)
    return events, exposure


def on_tree(newick: str, alignment: Mapping[str, str]) -> Lengths:
    events, exposure = tree_events(newick, alignment)
    counts = Counter(e.length for e in events if not e.terminal)
    return Lengths(dict(counts), sum(e.terminal for e in events), exposure)


def in_alignment(alignment: Mapping[str, str]) -> Lengths:
    interior: set[tuple[int, int]] = set()
    ends: set[tuple[int, int]] = set()
    for row in alignment.values():
        first = len(row) - len(row.lstrip(GAP))
        last = len(row.rstrip(GAP))
        for gap in GAPS.finditer(row):
            (interior if first <= gap.start() and gap.end() <= last else ends).add(gap.span())
    return Lengths(dict(Counter(end - start for start, end in interior)), len(ends))


def geometric(mean: float, last: int) -> list[float]:
    """A geometric distribution on 1, 2, ... with this mean, as `Lengths.binned` shares."""
    if mean <= 1:
        return [1.0, *[0.0] * last]
    p = 1 / mean
    bins = [p * (1 - p) ** (k - 1) for k in range(1, last + 1)]
    return [*bins, (1 - p) ** last]
