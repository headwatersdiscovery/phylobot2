"""Where a library's ancestors disagree across combinations.

The same ancestor — a clade — reconstructed under different aligners, trees or models can
carry a different ML (MAP) residue at a site. Such sites can be wrong even where the
posterior is confident, so disagreement is a signal independent of the posterior (disputed
is not the same as low-pp). The library spans combinations, so it can measure this:
for each clade it holds in two or more combinations, the sites where the combinations' ML
residues differ.

The combinations are laid on one reference's columns through the extant residues they share
(`alignment_map`, the same primitive `evaluation.disagreement` uses). Pure.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Any

from phylobot import alignment_map

GAP = "-"


@dataclass(frozen=True)
class State:
    """One combination's reconstruction of a clade: its combination key, the ML (gapped)
    sequence and the combination's extant alignment (the template's columns)."""

    key: str
    template: str
    alignment: dict[str, str]

    @property
    def usable(self) -> bool:
        return bool(self.alignment) and all(
            len(row) == len(self.template) for row in self.alignment.values()
        )


def _place_on(source: State, reference: dict[str, str]) -> dict[int, str]:
    """{reference column: `source`'s ML residue placed there}, by the extant residues the two
    alignments share (`alignment_map`); on a shared reference column the best-agreeing one
    wins, as the library's comparison viewer places them. `reference` is an extant alignment
    (the frame everything is laid on); it need not be one of the combinations."""
    mapping = alignment_map.map_columns(source.alignment, reference)
    best: dict[int, int] = {}
    for column, residue in enumerate(source.template):
        if residue == GAP:
            continue
        ref = mapping.target.get(column)
        if ref is None:  # an insertion the reference lacks: nothing to compare it to
            continue
        if ref not in best or mapping.agreement.get(column, 0.0) > mapping.agreement.get(
            best[ref], 0.0
        ):
            best[ref] = column
    return {ref: source.template[column] for ref, column in best.items()}


def _distinct(states: list[State]) -> list[State]:
    """Usable states, one per combination (gap versions share an alignment, so they are not
    aligner disagreement)."""
    out: list[State] = []
    seen: set[str] = set()
    for state in states:
        if state.usable and state.key not in seen:
            seen.add(state.key)
            out.append(state)
    return out


@dataclass(frozen=True)
class Dispute:
    """One clade's disagreement across its combinations."""

    combinations: int  # distinct combinations compared
    disputed: int  # reference columns where the ML residues differ
    per_key: dict[str, int]  # combination key -> sites it differs from the clade's majority


def compare(states: list[State]) -> Dispute | None:
    """The disputed sites among one clade's reconstructions; None for fewer than two distinct
    combinations. One state per combination is used (gap versions of the same combination
    share an alignment, so they are not aligner disagreement)."""
    usable = _distinct(states)
    if len(usable) < 2:
        return None
    reference = usable[0]
    width = len(reference.template)
    own = {c: reference.template[c] for c in range(width) if reference.template[c] != GAP}
    placed: dict[str, dict[int, str]] = {reference.key: own}
    for state in usable[1:]:
        placed[state.key] = _place_on(state, reference.alignment)
    disputed = 0
    per_key = {state.key: 0 for state in usable}
    for ref in range(width):
        residues = {key: columns.get(ref) for key, columns in placed.items()}
        present = [residue for residue in residues.values() if residue]
        if len(set(present)) > 1:
            disputed += 1
            majority = Counter(present).most_common(1)[0][0]
            for key, residue in residues.items():
                if residue is not None and residue != majority:
                    per_key[key] += 1
    return Dispute(len(usable), disputed, per_key)


def disputed_columns(states: list[State], reference: dict[str, str]) -> set[int]:
    """0-based columns of the `reference` alignment (the viewer's frame) where the clade's
    combinations place different ML residues — the in-viewer disagreement track. All
    combinations are laid on `reference`, so the columns line up with the frame's."""
    usable = _distinct(states)
    if len(usable) < 2 or not reference:
        return set()
    width = len(next(iter(reference.values())))
    placed = [_place_on(state, reference) for state in usable]
    disputed: set[int] = set()
    for column in range(width):
        present = {columns.get(column) for columns in placed}
        present.discard(None)
        if len(present) > 1:
            disputed.add(column)
    return disputed


@dataclass(frozen=True)
class Overview:
    """The library's cross-combination picture, for the summary line."""

    clades: int  # clades reconstructed in two or more combinations
    disputed_sites: int  # disputed reference columns, summed over those clades (disjoint)


def analyse(groups: dict[Any, list[State]]) -> tuple[dict[Any, Dispute], Overview]:
    """{clade key: its Dispute} for every clade with two or more combinations, and the
    library overview. `groups` maps a clade identity (the target's partition) to its states."""
    disputes: dict[Any, Dispute] = {}
    total = 0
    for clade, states in groups.items():
        dispute = compare(states)
        if dispute is not None:
            disputes[clade] = dispute
            total += dispute.disputed
    return disputes, Overview(len(disputes), total)
