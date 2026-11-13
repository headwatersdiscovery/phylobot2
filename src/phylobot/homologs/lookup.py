"""Find a protein's amino-acid sequence by name and species, upstream of homolog search.

The common starting point is a name, not a sequence: "I want IME2 in Saccharomyces
cerevisiae". This looks that up at UniProt's REST API (https://rest.uniprot.org), which
returns the sequence inline, so the user can pick a hit and drop its residues into the
homolog search below.

Only the name and species the user typed leave the machine, to UniProt. They are user
text, treated as hostile: cleaned to a safe character set here before they reach the
query, so nothing the user types can break out of the UniProt query it is quoted into.
"""

import json
import re
import urllib.parse
from dataclasses import dataclass
from functools import lru_cache
from typing import Any

from phylobot.homologs import remote
from phylobot.homologs.remote import SearchError

BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
TAXONOMY_URL = "https://rest.uniprot.org/taxonomy/search"
# The fields UniProt returns (format=json keeps the nested shape; this only trims it).
FIELDS = "accession,id,protein_name,gene_names,organism_name,length,reviewed,sequence"
MAX_RESULTS = 25
MAX_TERM = 80
# Keep only characters a gene/protein/species name needs; everything else (quotes,
# parentheses, Lucene operators) becomes a space, so the term is safe to quote into the
# UniProt query.
_KEEP = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 ._/+-")

# A short list of common model organisms, offered as a datalist beside the species field
# (scientific name, with a familiar label). Not a limit: any species can be typed.
MODEL_SPECIES = [
    ("Homo sapiens", "human"),
    ("Mus musculus", "mouse"),
    ("Rattus norvegicus", "rat"),
    ("Danio rerio", "zebrafish"),
    ("Xenopus laevis", "frog"),
    ("Gallus gallus", "chicken"),
    ("Bos taurus", "cow"),
    ("Drosophila melanogaster", "fruit fly"),
    ("Caenorhabditis elegans", "nematode"),
    ("Saccharomyces cerevisiae", "baker's yeast"),
    ("Schizosaccharomyces pombe", "fission yeast"),
    ("Candida albicans", "yeast pathogen"),
    ("Aspergillus nidulans", "mould"),
    ("Dictyostelium discoideum", "slime mould"),
    ("Arabidopsis thaliana", "thale cress"),
    ("Escherichia coli", "bacterium"),
]


@dataclass(frozen=True)
class Found:
    """One UniProt entry that matched the name and species."""

    accession: str  # primary accession, e.g. P32581
    entry_id: str  # UniProtKB ID, e.g. IME2_YEAST
    protein: str  # recommended (or submitted) protein name
    gene: str  # first gene name
    organism: str  # scientific name
    taxid: int | None
    length: int
    reviewed: bool  # Swiss-Prot (reviewed) vs TrEMBL
    sequence: str  # the amino-acid residues

    @property
    def url(self) -> str:
        return f"https://www.uniprot.org/uniprotkb/{self.accession}/entry"

    @property
    def job_name(self) -> str:
        """A sensible default search name when the user picks this hit."""
        label = self.gene or self.protein or self.accession
        return f"{label} ({self.organism})"[:100]


def clean_term(text: str, what: str) -> str:
    """A name or species reduced to a safe, quotable term; SearchError if it is empty."""
    kept = "".join(c if c in _KEEP else " " for c in text)
    cleaned = " ".join(kept.split())[:MAX_TERM]
    if not cleaned:
        raise SearchError(f"enter a {what}")
    return cleaned


def build_query(name: str, species: str, reviewed_only: bool) -> str:
    """The UniProt query string: the name as a gene or protein name, in the species if one
    is given (species is optional — without it, UniProt returns matches across organisms)."""
    name = clean_term(name, "protein or gene name")
    query = f'(gene:"{name}" OR protein_name:"{name}")'
    species = _soft_term(species)
    if species:
        query += f' AND organism_name:"{species}"'
    if reviewed_only:
        query += " AND reviewed:true"
    return query


def _protein_name(description: dict[str, Any]) -> str:
    if not isinstance(description, dict):
        return ""
    recommended = (description.get("recommendedName") or {}).get("fullName") or {}
    if value := recommended.get("value"):
        return str(value)
    for submitted in description.get("submissionNames") or []:
        if value := (submitted.get("fullName") or {}).get("value"):
            return str(value)
    return ""


def _gene_name(genes: list[Any]) -> str:
    for gene in genes or []:
        if value := (gene.get("geneName") or {}).get("value"):
            return str(value)
    return ""


def parse(raw: bytes) -> list[Found]:
    """UniProt's JSON into Found rows, reviewed (Swiss-Prot) entries first."""
    try:
        payload = json.loads(raw)
    except (json.JSONDecodeError, ValueError) as exc:
        raise SearchError("UniProt returned an answer that could not be read") from exc
    found: list[Found] = []
    for entry in payload.get("results") or []:
        if not isinstance(entry, dict):
            continue
        sequence = entry.get("sequence") or {}
        residues = str(sequence.get("value") or "")
        if not residues:
            continue  # no sequence to use
        organism = entry.get("organism") or {}
        taxid = organism.get("taxonId")
        found.append(
            Found(
                accession=str(entry.get("primaryAccession") or ""),
                entry_id=str(entry.get("uniProtkbId") or ""),
                protein=_protein_name(entry.get("proteinDescription") or {}),
                gene=_gene_name(entry.get("genes") or []),
                organism=str(organism.get("scientificName") or ""),
                taxid=int(taxid) if isinstance(taxid, int) else None,
                length=int(sequence.get("length") or len(residues)),
                reviewed="unreviewed" not in str(entry.get("entryType") or "").lower(),
                sequence=residues,
            )
        )
    found.sort(key=lambda f: not f.reviewed)  # stable: reviewed first, else UniProt's order
    return found


def search(name: str, species: str, reviewed_only: bool = False) -> list[Found]:
    """Look up `name` in `species` at UniProt; SearchError if the request fails."""
    params = {
        "query": build_query(name, species, reviewed_only),
        "format": "json",
        "fields": FIELDS,
        "size": str(MAX_RESULTS),
    }
    raw = remote.request(f"{BASE_URL}?{urllib.parse.urlencode(params)}")
    return parse(raw)


# ── type-ahead suggestions (the datalists under the two fields) ──────────────────────────
#
# Both fields suggest as the user types, from UniProt, so a species responds to its Latin
# name (Homo sapiens) and its common name (human) alike, and a gene to the symbols that
# exist in the chosen species. Answers are cached per process (lru_cache): the same few
# prefixes recur, and nothing here changes within a session. Suggestions never raise — an
# unreachable or odd answer just means no suggestions, so the field keeps working.

MIN_PREFIX = 2  # characters before suggesting (one letter matches everything)
MAX_SUGGESTIONS = 12
_GENE_TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9.@/_-]*")


@dataclass(frozen=True)
class Suggestion:
    value: str  # dropped into the field when chosen
    label: str  # shown in the datalist option


def _soft_term(text: str) -> str:
    """Like clean_term, but returns "" instead of raising (for a partial, as-you-type value)."""
    kept = "".join(c if c in _KEEP else " " for c in text)
    return " ".join(kept.split())[:MAX_TERM]


@lru_cache(maxsize=512)
def suggest_species(prefix: str) -> tuple[Suggestion, ...]:
    """Species whose scientific or common name matches the typed prefix, at UniProt's
    taxonomy. The scientific name is what goes in the field; the common name rides in the
    label, so "human" and "Homo" both find Homo sapiens."""
    prefix = _soft_term(prefix)
    if len(prefix) < MIN_PREFIX:
        return ()
    params = {
        "query": prefix,
        "fields": "scientific_name,common_name",
        "size": "25",
        "format": "json",
    }
    try:
        payload = json.loads(remote.request(f"{TAXONOMY_URL}?{urllib.parse.urlencode(params)}"))
    except SearchError, json.JSONDecodeError, ValueError:
        return ()
    out: list[Suggestion] = []
    seen: set[str] = set()
    for taxon in payload.get("results") or []:
        if not isinstance(taxon, dict):
            continue
        sci = str(taxon.get("scientificName") or "").strip()
        common = str(taxon.get("commonName") or "").strip()
        if sci and sci.lower() not in seen:
            seen.add(sci.lower())
            out.append(Suggestion(sci, f"{sci} ({common})" if common else sci))
        if len(out) >= MAX_SUGGESTIONS:
            break
    return tuple(out)


@lru_cache(maxsize=512)
def suggest_genes(prefix: str, species: str) -> tuple[Suggestion, ...]:
    """Gene symbols starting with the typed prefix, in the chosen species if one is given
    (genes are species-specific, so scoping both sharpens and speeds the list)."""
    match = _GENE_TOKEN.match(prefix.strip())
    token = match.group(0)[:MAX_TERM] if match else ""
    if len(token) < MIN_PREFIX:
        return ()
    query = f"gene:{token}*"
    species = _soft_term(species)
    if species:
        query += f' AND organism_name:"{species}"'
    params = {"query": query, "fields": "gene_names", "size": "50", "format": "json"}
    try:
        payload = json.loads(remote.request(f"{BASE_URL}?{urllib.parse.urlencode(params)}"))
    except SearchError, json.JSONDecodeError, ValueError:
        return ()
    lowered = token.lower()
    out: list[Suggestion] = []
    seen: set[str] = set()
    for entry in payload.get("results") or []:
        for gene in (entry or {}).get("genes") or []:
            name = str((gene.get("geneName") or {}).get("value") or "")
            if name.lower().startswith(lowered) and name.lower() not in seen:
                seen.add(name.lower())
                out.append(Suggestion(name, name))
    out.sort(key=lambda s: s.value.lower())
    return tuple(out[:MAX_SUGGESTIONS])
