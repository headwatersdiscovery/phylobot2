"""Hits as FASTA, for a job or a download.

Names are readable and unique, organism then accession ("Yarrowia_lipolytica_Q6C4M5"),
because users read them when choosing an outgroup and on the tree. The description
keeps the protein name, organism and database ID. Everything comes from remote data,
so names are reduced to letters, digits and '_', and descriptions to printable text.
"""

import re

from phylobot.homologs.hits import Hit

WIDTH = 60


def accession(hit_id: str) -> str:
    """'UniRef90_P32581' -> 'P32581'; NCBI accessions stay as they are."""
    return hit_id.split("_", 1)[1] if hit_id.startswith("UniRef") else hit_id


def organism(taxon: str) -> str:
    """The first two words of a taxon name, e.g. genus and species."""
    words = re.sub(r"\(.*?\)", " ", taxon).split()
    return "_".join(words[:2])


def clean(text: str, keep: str) -> str:
    return re.sub(keep, "_", text).strip("_")


def names(hits: list[Hit]) -> dict[str, str]:
    """{hit ID: FASTA name}, unique."""
    out: dict[str, str] = {}
    seen: set[str] = set()
    for hit in hits:
        base = clean(f"{organism(hit.taxon)}_{accession(hit.id)}", r"[^A-Za-z0-9]+")
        name, n = base or "hit", 2
        while name in seen:
            name, n = f"{base}_{n}", n + 1
        seen.add(name)
        out[hit.id] = name
    return out


def description(hit: Hit) -> str:
    text = f"{hit.name} [{hit.taxon}] {hit.id}" if hit.taxon else f"{hit.name} {hit.id}"
    return re.sub(r"[^\x20-\x7e]+", " ", text).replace(">", " ").strip()


def to_fasta(hits: list[Hit], sequences: dict[str, str]) -> tuple[str, list[Hit]]:
    """FASTA of the hits that have a full sequence, and the hits that have none."""
    named = names(hits)
    lines, missing = [], []
    for hit in hits:
        seq = sequences.get(hit.id)
        if not seq:
            missing.append(hit)
            continue
        lines.append(f">{named[hit.id]} {description(hit)}")
        lines.extend(seq[i : i + WIDTH] for i in range(0, len(seq), WIDTH))
    return ("\n".join(lines) + "\n" if lines else ""), missing


def taxids(hits: list[Hit]) -> str:
    """ "name<TAB>taxid" lines for the hits with a taxid, named as in to_fasta."""
    named = names(hits)
    return "".join(f"{named[hit.id]}\t{hit.taxid}\n" for hit in hits if hit.taxid is not None)
