"""Minimal FASTA reading and writing."""

from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA file into {name: sequence}, preserving file order."""
    records: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            if name in records:
                raise ValueError(f"duplicate sequence name {name!r} in {path}")
            records[name] = ""
        elif name is not None:
            records[name] += line.strip()
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))
