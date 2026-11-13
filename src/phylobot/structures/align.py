"""Local pairwise alignment of two protein sequences (Smith-Waterman, affine gaps,
BLOSUM62), in numpy: a row at a time, so a 1000 x 1000 alignment takes milliseconds in the
web tier.

Local, because a structure is seldom the same stretch of protein as a job's sequences: a
crystallised domain against full-length sequences, or the reverse.
"""

from dataclasses import dataclass

import numpy as np

ORDER = "ARNDCQEGHILKMFPSTWYV"
_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
OTHER = len(ORDER)  # any other letter (X, B, Z, U...): -1 against everything
BLOSUM62 = np.full((OTHER + 1, OTHER + 1), -1, dtype=np.int32)
BLOSUM62[:OTHER, :OTHER] = np.array(
    [[int(v) for v in row.split()] for row in _ROWS.strip().splitlines()], dtype=np.int32
)
_INDEX = {letter: i for i, letter in enumerate(ORDER)}
GAP_OPEN = 11  # the first residue of a gap
GAP_EXTEND = 1  # each further one
MAX_LENGTH = 6000  # per sequence: 36M cells, the most this will fill


@dataclass(frozen=True)
class Local:
    score: int
    pairs: tuple[tuple[int, int], ...]  # (index in a, index in b) of aligned residues
    identical: int

    @property
    def identity(self) -> float:
        """Identical residues among the aligned pairs."""
        return self.identical / len(self.pairs) if self.pairs else 0.0


def encode(sequence: str) -> np.ndarray:
    return np.array([_INDEX.get(c, OTHER) for c in sequence.upper()], dtype=np.intp)


def local(a: str, b: str) -> Local:
    """The best local alignment of `a` and `b`."""
    n, m = len(a), len(b)
    if not n or not m or n > MAX_LENGTH or m > MAX_LENGTH:
        return Local(0, (), 0)
    ea, eb = encode(a), encode(b)
    steps = np.arange(m, dtype=np.int32) * GAP_EXTEND
    h = np.zeros((n + 1, m + 1), dtype=np.int32)
    below = np.full(m, -(10**6), dtype=np.int32)  # a gap in b, ending in the row above
    for i in range(1, n + 1):
        previous = h[i - 1]
        # A gap in b (a's residue against nothing): opened from above, or extended.
        below = np.maximum(previous[1:] - GAP_OPEN, below - GAP_EXTEND)
        best = np.maximum(np.maximum(previous[:-1] + BLOSUM62[ea[i - 1], eb], below), 0)
        # A gap in a runs along the row being filled. It is never worth opening one from a
        # cell that itself ends such a gap (extending is cheaper), so it can be taken from
        # `best` alone: the largest best[k] - open - extend * (j - k - 1) over k < j.
        reach = np.maximum.accumulate(best + steps)
        along = np.empty(m, dtype=np.int32)
        along[0] = -(10**6)
        along[1:] = reach[:-1] - GAP_OPEN - steps[1:] + GAP_EXTEND
        h[i, 1:] = np.maximum(best, along)
    score = int(h.max())
    if score <= 0:
        return Local(0, (), 0)
    i, j = (int(v) for v in np.unravel_index(int(h.argmax()), h.shape))
    pairs: list[tuple[int, int]] = []
    while i > 0 and j > 0 and h[i, j] > 0:
        here = int(h[i, j])
        if here == int(h[i - 1, j - 1]) + int(BLOSUM62[ea[i - 1], eb[j - 1]]):
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
            continue
        # Which gap ends here: the cell it was opened from is `length` back.
        up = h[i - 1 :: -1, j][:i] - GAP_OPEN - np.arange(i, dtype=np.int32) * GAP_EXTEND
        left = h[i, j - 1 :: -1][:j] - GAP_OPEN - np.arange(j, dtype=np.int32) * GAP_EXTEND
        from_up = np.flatnonzero(up == here)
        from_left = np.flatnonzero(left == here)
        if len(from_up):
            i -= int(from_up[0]) + 1
        elif len(from_left):
            j -= int(from_left[0]) + 1
        else:  # cannot happen; never loop for ever on it
            break
    pairs.reverse()
    identical = sum(a[x].upper() == b[y].upper() for x, y in pairs)
    return Local(score, tuple(pairs), identical)
