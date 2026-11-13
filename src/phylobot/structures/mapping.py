"""Which alignment column each residue of a structure belongs to.

A structure's chain is aligned (`align.local`) with the job's extant sequence most like
it, its *reference*; a residue aligned with a residue of the reference takes that residue's
column, in whichever of the job's alignments is on the page. So the homology statement is
the job's own alignment, and one pairwise alignment serves every aligner and model.

A chain like none of the job's sequences is left unmapped rather than forced: a ligand
peptide, an antibody, or simply the wrong protein.
"""

from collections.abc import Sequence
from dataclasses import dataclass

from phylobot.structures import align, pdb

GAP = "-"
KMER = 3
CANDIDATES = 3  # sequences aligned in full, the best of them by shared k-mers
MIN_PAIRS = 10
MIN_IDENTITY = 0.25
WEAK_IDENTITY = 0.4  # below this the page says the mapping is a guess


@dataclass(frozen=True)
class ChainMap:
    chain: str
    reference: str  # the extant sequence the chain was aligned with
    identity: float  # identical residues among the aligned ones
    residues: tuple[int | None, ...]  # per chain residue: its index in the reference

    @property
    def mapped(self) -> int:
        return sum(r is not None for r in self.residues)

    @property
    def weak(self) -> bool:
        return self.identity < WEAK_IDENTITY


def ungapped(alignment: dict[str, str]) -> dict[str, str]:
    return {name: row.replace(GAP, "").upper() for name, row in alignment.items()}


def _kmers(sequence: str) -> set[str]:
    return {sequence[i : i + KMER] for i in range(len(sequence) - KMER + 1)}


def reference(sequence: str, sequences: dict[str, str]) -> tuple[str, align.Local] | None:
    """The extant sequence most like `sequence`, and their alignment; None if none is
    alike enough."""
    mine = _kmers(sequence)
    shared = {name: len(mine & _kmers(other)) for name, other in sequences.items()}
    ranked = sorted(sequences, key=lambda name: -shared[name])[:CANDIDATES]
    best: tuple[str, align.Local] | None = None
    for name in ranked:
        found = align.local(sequence, sequences[name])
        if best is None or found.score > best[1].score:
            best = (name, found)
    if best is None or len(best[1].pairs) < MIN_PAIRS or best[1].identity < MIN_IDENTITY:
        return None
    return best


def chain_maps(parsed: pdb.Parsed, sequences: dict[str, str]) -> dict[str, ChainMap]:
    """A map for every chain that is like one of `sequences` (ungapped extant sequences)."""
    by_sequence: dict[str, tuple[str, align.Local] | None] = {}
    maps = {}
    for chain in parsed.chains:
        sequence = chain.sequence
        if sequence not in by_sequence:  # the copies of a homo-oligomer are aligned once
            by_sequence[sequence] = reference(sequence, sequences)
        found = by_sequence[sequence]
        if found is None:
            continue
        name, local = found
        residues: list[int | None] = [None] * len(chain.residues)
        for mine, theirs in local.pairs:
            residues[mine] = theirs
        maps[chain.id] = ChainMap(chain.id, name, local.identity, tuple(residues))
    return maps


def columns(chain_map: ChainMap, row: str) -> list[int | None]:
    """Each chain residue's column (from 0) in the alignment whose row for the chain's
    reference is `row`."""
    of_residue = [i for i, c in enumerate(row) if c != GAP]
    return [
        of_residue[r] if r is not None and r < len(of_residue) else None for r in chain_map.residues
    ]


@dataclass(frozen=True)
class Site:
    """One residue of the structure with what the page's sequences have at its column."""

    key: str  # pdb.key: "A:295"
    chain: str
    residue: pdb.Residue
    column: int | None = None  # from 0; None where the residue has no column
    state: str = ""  # the ancestor's state there ("-" for a gap, "" without a column)
    pp: float | None = None
    other: str = ""  # the second sequence's state, when two are compared
    other_pp: float | None = None

    @property
    def mapped(self) -> bool:
        return self.column is not None

    @property
    def change(self) -> str:
        """Between the two sequences: "" (the same, or nothing to compare), "sub",
        "lost" (a residue, then none) or "gained"."""
        if not self.mapped or not self.other or self.state == self.other:
            return ""
        if self.other == GAP:
            return "lost"
        return "gained" if self.state == GAP else "sub"

    @property
    def differs(self) -> str:
        """The ancestor against the structure's own residue: "" (the same, or no column),
        "sub" or "absent" (the ancestor has no residue here)."""
        if not self.mapped or self.state == self.residue.aa:
            return ""
        return "absent" if self.state == GAP else "sub"


def sites(
    chain: pdb.Chain,
    cols: Sequence[int | None],
    row: str,
    pp: Sequence[float] | None = None,
    other: str = "",
    other_pp: Sequence[float] | None = None,
) -> list[Site]:
    """The chain's residues against the ancestor's aligned `row` (and a second sequence's,
    `other`). `pp`: the posterior probability of the ancestor's state per column."""
    out = []
    for residue, column in zip(chain.residues, cols, strict=True):
        key = pdb.key(chain.id, residue)
        if column is None or column >= len(row):
            out.append(Site(key, chain.id, residue))
            continue
        state = row[column].upper()
        out.append(
            Site(
                key,
                chain.id,
                residue,
                column,
                state,
                pp[column] if pp and state != GAP else None,
                other[column].upper() if other else "",
                other_pp[column] if other and other_pp and other[column] != GAP else None,
            )
        )
    return out


@dataclass(frozen=True)
class Between:
    """Residues of the page's sequences in columns the structure has no residue for,
    between two residues that follow each other in the chain: an insertion relative to the
    structure, or a stretch the structure does not resolve."""

    after: str  # key of the residue before
    before: str  # key of the residue after
    first: int  # the columns, from 0, inclusive
    last: int
    count: int  # the ancestor's residues there
    other: int  # the second sequence's


def between(found: Sequence[Site], row: str, other: str = "") -> list[Between]:
    out = []
    mapped = [s for s in found if s.column is not None]
    for a, b in zip(mapped, mapped[1:], strict=False):
        assert a.column is not None and b.column is not None
        if b.column - a.column < 2:
            continue
        span = slice(a.column + 1, b.column)
        mine = sum(c != GAP for c in row[span])
        theirs = sum(c != GAP for c in other[span]) if other else 0
        if mine or theirs:
            out.append(Between(a.key, b.key, a.column + 1, b.column - 1, mine, theirs))
    return out
