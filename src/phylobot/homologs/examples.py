"""Saved searches to try the homolog tools on, and the tests' pinned fixtures (remote
results change): S. cerevisiae IME2 (645 aa), blastp, E <= 1e-3, 1,000 hits asked for,
2026-09-29."""

import gzip
from dataclasses import dataclass
from functools import cache
from importlib.resources import files

from phylobot.homologs.hits import SearchResult, parse_ebi_xml, parse_ncbi_xml


@dataclass(frozen=True)
class Example:
    label: str
    source: str
    file: str
    database: str = ""
    searched: str = "2026-09-29"
    took: str = ""
    query_name: str = "Saccharomyces.cerevisiae.IME2"  # in examples/ime2.fasta

    @property
    def sequences_file(self) -> str:
        return self.file.replace(".xml.gz", ".sequences.fasta.gz")


EXAMPLES = {
    "ime2-ebi": Example(
        "Ime2 at EBI", "ebi", "ime2_uniref90.xml.gz", database="uniref90", took="16 min"
    ),
    "ime2-ncbi": Example(
        "Ime2 at NCBI", "ncbi", "ime2_clustered_nr.xml.gz", database="clustered_nr", took="39 s"
    ),
}

# Which examples the "Try a saved search" section on the homologs page lists. The Ime2
# searches above stay in EXAMPLES — they are the tests' pinned fixtures and remain
# reachable at /homologs/examples/<key>/ — but are no longer featured, pending better
# replacements (e.g. MT-CO1). Add a key here to show it; when a real search becomes an
# example, gzip its result + sequences into examples/searches/, add an Example above, and
# list its key here.
FEATURED: tuple[str, ...] = ()


def raw(key: str) -> bytes:
    """The example's search XML as the service returned it."""
    path = files("phylobot").joinpath("examples", "searches", EXAMPLES[key].file)
    return gzip.decompress(path.read_bytes())


@cache
def load(key: str) -> SearchResult:
    """The parsed example (KeyError for an unknown key)."""
    example = EXAMPLES[key]
    if example.source == "ncbi":
        return parse_ncbi_xml(raw(key), 1000)
    return parse_ebi_xml(raw(key))


@cache
def sequences(key: str) -> dict[str, str]:
    """The example's hits' full-length sequences, {hit ID: sequence} (fetched once, as
    the runner fetches them for a real search)."""
    from phylobot.homologs.sequences import parse_fasta

    path = files("phylobot").joinpath("examples", "searches", EXAMPLES[key].sequences_file)
    return parse_fasta(gzip.decompress(path.read_bytes()).decode())


def query(key: str) -> str:
    """The example's query residues (from the bundled Ime2 FASTA)."""
    from phylobot.homologs.sequences import parse_fasta

    text = files("phylobot").joinpath("examples", "ime2.fasta").read_text()
    return parse_fasta(text)[EXAMPLES[key].query_name]
