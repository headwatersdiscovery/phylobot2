"""The New Job page's input files, before the pipeline's validation sees them (pure).

Sequences can come in several FASTA files (`merge_sequences`): they are combined in
the order given, as long as they hold one type of sequence. A name in more than one
file is either the same sequence again (one copy is kept) or a different one, which
the user resolves: rename the later copies (`name.2`, `name.3`, ...) or keep only the
first file's. Several alignments can be uploaded too; when their sequences differ, the
job can use the sequences they share (`shared_names`, `restrict`), and leaves out the
rest. Everything here is advice and bookkeeping: `phylobot.pipeline.validate` still
checks every file the job runs on.
"""

import re
from dataclasses import dataclass, field

from phylobot.pipeline.validate import InputError, detect_seqtype, parse_fasta

RENAME = "rename"  # resolve a name clash by renaming the later copies
FIRST = "first"  # ... or by keeping only the first file's copy
DUPLICATE_CHOICES = (RENAME, FIRST)
MAX_SEQUENCE_FILES = 20
PASTED = "your pasted sequences"  # the source name of step 2's text
GAP = re.compile(r"[-.]")


def residues(raw: str) -> str:
    """A sequence as compared across files: uppercase, gaps removed."""
    return GAP.sub("", raw).upper()


@dataclass
class FileStats:
    """One input file, as the New Job page lists it."""

    name: str
    sequences: int = 0
    min_length: int = 0  # residues, or nucleotides for DNA (repeated names included)
    max_length: int = 0
    columns: int = 0  # alignments: the aligned length (0 if rows differ in length)
    seqtype: str = ""  # detected: "protein" or "codon"
    error: str = ""
    names: list[str] = field(default_factory=list)

    @property
    def kind(self) -> str:
        return {"codon": "coding DNA", "protein": "protein"}.get(self.seqtype, "")

    @property
    def summary(self) -> str:
        """E.g. "109 coding DNA sequences, 204–238 codons" or "12 sequences × 852 columns"."""
        dna = self.seqtype == "codon"
        kind = "coding DNA sequence" if dna else "protein sequence"
        count = f"{self.sequences} {kind}{'s' * (self.sequences != 1)}"
        if self.columns:
            return f"{count} × {self.columns} columns"
        low, high = (self.min_length // 3, self.max_length // 3) if dna else (
            self.min_length, self.max_length)  # fmt: skip
        span = f"{low}" if low == high else f"{low}–{high}"
        return f"{count}, {span} {'codons' if dna else 'residues'}"


@dataclass(frozen=True)
class Clash:
    """A name in more than one file with different sequences."""

    name: str
    first: str  # the file that has it first
    later: str  # a later file with another sequence under that name
    renamed: str  # its name if the later copies are renamed


@dataclass
class Merged:
    files: list[FileStats]
    fasta: bytes  # the combined FASTA (clashes resolved as chosen, else first copy kept)
    identical: list[tuple[str, str]]  # (name, later file): the same sequence again, dropped
    clashes: list[Clash]
    error: str = ""  # a problem no choice resolves (mixed types, an unreadable file)
    resolution: str = ""  # RENAME, FIRST or "" (not chosen yet)
    incompatible: bool = False  # the error is a mix of coding DNA and protein files

    @property
    def unresolved(self) -> bool:
        return bool(self.clashes) and not self.resolution


def read_stats(
    name: str, content: bytes, aligned: bool = False
) -> tuple[FileStats, list[tuple[str, str, str]]]:
    """(the file's stats, its (name, description, sequence) records)."""
    stats = FileStats(name)
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError:
        stats.error = "not UTF-8 text"
        return stats, []
    if any(not line.startswith(">") and line.strip() for line in text.splitlines()):
        stats.seqtype = detect_seqtype(text)  # known even if the file has another problem
    try:
        records = parse_fasta(text, name, unique=False)  # repeated names: `resolve`
    except InputError as exc:
        stats.error = str(exc).removeprefix(f"{name}: ").removeprefix(f"{name}, ")
        return stats, []
    if not records:
        stats.error = "no sequences"
        return stats, []
    lengths = [len(residues(seq)) for _, _, seq in records]
    stats.sequences = len(records)
    stats.min_length, stats.max_length = min(lengths), max(lengths)
    stats.names = [n for n, _, _ in records]
    if aligned:
        widths = {len(seq) for _, _, seq in records}
        stats.columns = widths.pop() if len(widths) == 1 else 0
        if not stats.columns:
            stats.error = "rows have different lengths: not an alignment"
    return stats, records


def renamed(name: str, taken: set[str]) -> str:
    """`name.2`, or the first of `name.3`, `name.4`, ... not taken."""
    n = 2
    while f"{name}.{n}" in taken:
        n += 1
    return f"{name}.{n}"


def merge_sequences(files: list[tuple[str, bytes]], resolution: str = "") -> Merged:
    """Several FASTA files as one, in order. `resolution` decides name clashes (RENAME
    or FIRST); until it is chosen, the first copy is kept so the rest can be checked."""
    stats, everything = [], []
    for name, content in files:
        file_stats, records = read_stats(name, content)
        stats.append(file_stats)
        everything.append((name, records))
    chosen = resolution if resolution in DUPLICATE_CHOICES else ""
    merged = Merged(stats, b"", [], [], resolution=chosen)
    # Incompatible types first: whatever else is wrong, these files can't go together.
    if len({s.seqtype for s in stats if s.seqtype}) > 1:
        dna = [s.name for s in stats if s.seqtype == "codon"]
        protein = [s.name for s in stats if s.seqtype == "protein"]
        merged.error = (
            f"Incompatible files: {', '.join(dna)} "
            f"{'is' if len(dna) == 1 else 'are'} coding DNA but {', '.join(protein)} "
            f"{'is' if len(protein) == 1 else 'are'} protein. A job takes one type of "
            "sequence: untick or remove one"
        )
        merged.incompatible = True
        return merged
    if bad := [s for s in stats if s.error]:
        merged.error = f"{bad[0].name}: {bad[0].error}"
        return merged
    resolve(merged, everything)
    return merged


def resolve(merged: Merged, sources: list[tuple[str, list[tuple[str, str, str]]]]) -> None:
    """Fill `merged` from records in order: a name seen before (in any file, or earlier
    in the same one) with the same residues is dropped (`identical`); with others it is
    a clash, renamed or dropped as `merged.resolution` says (dropped until chosen)."""
    seen: dict[str, tuple[str, str]] = {}  # name -> (file, residues)
    taken = {n for _, records in sources for n, _, _ in records}
    lines = []
    for filename, records in sources:
        for name, description, raw in records:
            header = name
            if name in seen:
                if residues(raw) == seen[name][1]:
                    merged.identical.append((name, filename))
                    continue
                new = renamed(name, taken)
                taken.add(new)
                merged.clashes.append(Clash(name, seen[name][0], filename, new))
                if merged.resolution != RENAME:
                    continue
                header = new
            else:
                seen[name] = (filename, residues(raw))
            lines.append(f">{header} {description}".rstrip() + f"\n{raw}\n")
    merged.fasta = "".join(lines).encode()


def dedupe_alignment(name: str, content: bytes, resolution: str = "") -> Merged:
    """One alignment with its repeated names resolved as `resolve` does; rows are kept
    as aligned. `merged.error` if it can't be read."""
    stats, records = read_stats(name, content, aligned=True)
    chosen = resolution if resolution in DUPLICATE_CHOICES else ""
    merged = Merged([stats], content, [], [], resolution=chosen)
    if stats.error:
        merged.error = f"{name}: {stats.error}"
        return merged
    resolve(merged, [(name, records)])
    return merged


def shared_names(name_sets: list[set[str]]) -> set[str]:
    """The names every source has."""
    return set.intersection(*name_sets) if name_sets else set()


def restrict(content: bytes, keep: set[str], aligned: bool) -> bytes:
    """A FASTA with only the sequences in `keep`; for an alignment, columns left all
    gaps are dropped too."""
    records = [r for r in parse_fasta(content.decode("utf-8"), "file") if r[0] in keep]
    if aligned and records:
        width = len(records[0][2])
        used = [i for i in range(width) if any(seq[i : i + 1] not in "-." for _, _, seq in records)]
        records = [(n, d, "".join(seq[i] for i in used)) for n, d, seq in records]
    return "".join(f">{n} {d}".rstrip() + f"\n{s}\n" for n, d, s in records).encode()
