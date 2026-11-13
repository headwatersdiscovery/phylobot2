"""Mapping columns between alignments of the same sequences.

Every alignment of a job aligns the same extant sequences, so each extant residue — say
residue 212 of one sequence — is a fixed point every alignment places in some column.
That is enough to map columns between alignments without aligning anything again: a
column of one alignment goes to the column of the other where most of its residues sit.
How cleanly it maps is itself a measurement: where two aligners put the same residues in
a column they agree about homology there, and where they split them, ancestors that
differ at that site differ because of the alignment, not the reconstruction.

Used by the gap scoring of the golden tests (`evaluation.gaps`, inferred against true
alignment) and by the library page's overview, which lays ancestors from several
alignments out on one frame. Pure Python; columns are 0-based throughout.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

GAP = "-"

Residue = tuple[str, int]  # (sequence name, 0-based index in the ungapped sequence)


def residue_columns(row: str) -> list[int]:
    """The column of each residue of an aligned row, in order."""
    return [i for i, c in enumerate(row) if c != GAP]


def column_residues(alignment: dict[str, str]) -> dict[int, frozenset[Residue]]:
    """{column: the residues in it} for every column holding at least one."""
    found: dict[int, set[Residue]] = {}
    for name, row in alignment.items():
        for index, column in enumerate(residue_columns(row)):
            found.setdefault(column, set()).add((name, index))
    return {column: frozenset(residues) for column, residues in found.items()}


def same_sequences(a: dict[str, str], b: dict[str, str]) -> bool:
    """Whether two alignments hold the same sequences (names and residues), the condition
    for mapping between them."""
    return a.keys() == b.keys() and all(
        a[name].replace(GAP, "") == b[name].replace(GAP, "") for name in a
    )


@dataclass(frozen=True)
class ColumnMap:
    """Columns of a source alignment mapped onto a reference alignment."""

    target: dict[int, int]  # source column -> reference column (majority of its residues)
    agreement: dict[int, float]  # source column -> Jaccard of its residues and its target's
    reference_agreement: dict[int, float] = field(default_factory=dict)  # best per ref column

    def exact(self, column: int) -> bool:
        """Both alignments put exactly these residues in one column."""
        return self.agreement.get(column, 0.0) >= 1.0


def majority_columns(source: dict[str, str], reference: dict[str, str]) -> dict[int, int]:
    """Each source column holding a residue -> the reference column where most of its
    residues sit (ties go to the column met first, scanning sequences in order)."""
    votes: dict[int, Counter[int]] = {}
    for name, source_row in source.items():
        source_cols = residue_columns(source_row)
        reference_cols = residue_columns(reference[name])
        if len(source_cols) != len(reference_cols):
            raise ValueError(f"{name}: different residue counts in the two alignments")
        for source_col, reference_col in zip(source_cols, reference_cols, strict=True):
            votes.setdefault(source_col, Counter())[reference_col] += 1
    return {column: counts.most_common(1)[0][0] for column, counts in votes.items()}


def map_columns(source: dict[str, str], reference: dict[str, str]) -> ColumnMap:
    """Map each source column holding a residue to the reference column where most of its
    residues sit (ties go to the column met first, scanning sequences in order), with the
    Jaccard index of the two columns' residue sets: 1 exactly when the aligners agree on
    that column. Both alignments must hold the same sequences (same names, same residues)."""
    target = majority_columns(source, reference)
    source_residues = column_residues(source)
    reference_residues = column_residues(reference)
    agreement: dict[int, float] = {}
    best: dict[int, float] = {}
    for column, mapped in target.items():
        mine, theirs = source_residues[column], reference_residues[mapped]
        score = len(mine & theirs) / len(mine | theirs)
        agreement[column] = score
        best[mapped] = max(best.get(mapped, 0.0), score)
    return ColumnMap(target, agreement, best)


def column_map(inferred: dict[str, str], true: dict[str, str]) -> dict[int, int]:
    """Map inferred-alignment columns to true-alignment columns by majority vote over the
    extant residues they contain. Columns holding no extant residue are left unmapped.
    Both alignments must contain the same sequences (same names, same residues)."""
    return map_columns(inferred, true).target


@dataclass(frozen=True)
class SequenceFrame:
    """An alignment's columns placed on one extant sequence's residues: a column holding a
    residue of it sits at that residue; any other column is an insertion relative to it,
    between residue `slot` and the next."""

    length: int  # residues in the sequence
    residue: list[int | None]  # per column: 1-based residue number, or None
    slot: list[int]  # per column: residues of the sequence before it


def sequence_frame(row: str) -> SequenceFrame:
    """`row` is the chosen extant sequence's aligned row in the alignment being placed."""
    residue: list[int | None] = []
    slot: list[int] = []
    count = 0
    for character in row:
        slot.append(count)
        if character == GAP:
            residue.append(None)
        else:
            count += 1
            residue.append(count)
    return SequenceFrame(count, residue, slot)


def residue_agreement(alignments: list[dict[str, str]], name: str) -> list[float]:
    """For each residue of sequence `name`, how far the alignments agree on what is
    homologous to it: the lowest Jaccard index, over every pair of alignments, between the
    residue sets of the columns that hold it. 1 where every aligner builds the same column."""
    columns = [column_residues(a) for a in alignments]
    positions = [residue_columns(a[name]) for a in alignments]
    length = len(positions[0]) if positions else 0
    scores: list[float] = []
    for index in range(length):
        sets = [columns[k][positions[k][index]] for k in range(len(alignments))]
        lowest = 1.0
        for i in range(len(sets)):
            for j in range(i + 1, len(sets)):
                lowest = min(lowest, len(sets[i] & sets[j]) / len(sets[i] | sets[j]))
        scores.append(lowest)
    return scores
