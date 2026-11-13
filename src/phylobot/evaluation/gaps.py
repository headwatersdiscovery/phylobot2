"""Scoring ancestral sequences that may contain gaps (indel scenarios).

An inferred ancestor is a row of the *inferred* alignment; the true ancestor is a row
of the *true* alignment. The two alignments' columns are linked through the extant
sequences, whose residues occur in both: each inferred column maps to the true column
where most of its residues sit (`column_map`). With that map, an ancestor's
gap/residue call at every inferred column can be checked against the truth
(`gap_calls`). Identity between ancestors of different lengths uses a global pairwise
alignment (`aligned_identity`). The column mapping is `phylobot.alignment_map`'s, shared
with the library page, which maps a job's alignments onto each other the same way.
"""

from dataclasses import dataclass

from phylobot.alignment_map import column_map, residue_columns

__all__ = ["GapCalls", "aligned_identity", "column_map", "gap_calls", "residue_columns"]

GAP = "-"


@dataclass
class GapCalls:
    """Per-column gap/residue calls of an inferred ancestor against its true ancestor,
    over the inferred columns that map to a true column. "Positive" means gap."""

    true_gap: int = 0  # both gap
    false_gap: int = 0  # called gap, truly a residue
    missed_gap: int = 0  # called residue, truly a gap
    true_residue: int = 0  # both residue
    lost_residues: int = 0  # true residues in true columns no inferred column maps to

    def __iadd__(self, other: GapCalls) -> GapCalls:
        for field in self.__dataclass_fields__:
            setattr(self, field, getattr(self, field) + getattr(other, field))
        return self

    @property
    def columns(self) -> int:
        return self.true_gap + self.false_gap + self.missed_gap + self.true_residue

    @property
    def accuracy(self) -> float | None:
        return (self.true_gap + self.true_residue) / self.columns if self.columns else None

    @property
    def gap_recall(self) -> float | None:
        """Share of the truly-gap columns that were called gaps."""
        truly_gap = self.true_gap + self.missed_gap
        return self.true_gap / truly_gap if truly_gap else None

    @property
    def gap_precision(self) -> float | None:
        """Share of the called gaps that truly are gaps (None if no gaps were called)."""
        called = self.true_gap + self.false_gap
        return self.true_gap / called if called else None


def gap_calls(inferred_row: str, true_row: str, mapping: dict[int, int]) -> GapCalls:
    calls = GapCalls()
    for inferred_col, true_col in mapping.items():
        called_gap = inferred_row[inferred_col] == GAP
        truly_gap = true_row[true_col] == GAP
        if called_gap and truly_gap:
            calls.true_gap += 1
        elif called_gap:
            calls.false_gap += 1
        elif truly_gap:
            calls.missed_gap += 1
        else:
            calls.true_residue += 1
    mapped_true = set(mapping.values())
    calls.lost_residues = sum(
        1 for col, char in enumerate(true_row) if char != GAP and col not in mapped_true
    )
    return calls


def aligned_identity(inferred: str, truth: str) -> float:
    """Identity of two ungapped sequences after a global alignment (Needleman-Wunsch,
    linear gaps): matches / length of the longer sequence, so extra or missing
    residues lower it. Case-insensitive."""
    a, b = inferred.replace(GAP, "").upper(), truth.replace(GAP, "").upper()
    if not a or not b:
        return 0.0
    match, mismatch, gap = 2, -1, -2
    # score[j] for the current row, with matches counted along the best path
    previous = [(j * gap, 0) for j in range(len(b) + 1)]
    for i in range(1, len(a) + 1):
        current = [(i * gap, 0)]
        for j in range(1, len(b) + 1):
            same = a[i - 1] == b[j - 1]
            diagonal = (previous[j - 1][0] + (match if same else mismatch),
                        previous[j - 1][1] + same)  # fmt: skip
            up = (previous[j][0] + gap, previous[j][1])
            left = (current[j - 1][0] + gap, current[j - 1][1])
            current.append(max(diagonal, up, left))
        previous = current
    return previous[-1][1] / max(len(a), len(b))
