"""Structures from public databases: a PDB entry from RCSB, a predicted model from
AlphaFold DB, and UniProt's list of both for a protein found by name.

What leaves the machine is an identifier (a PDB ID or a UniProt accession) or the protein
and species names the user typed, never a sequence of the job's. Identifiers are checked
against their fixed shapes before they are put into a URL, and a file is only ever
fetched from the two hosts named here: the address AlphaFold's API answers with is
checked too, since it is somebody else's text.
"""

import json
import re
import urllib.parse
from dataclasses import dataclass, field
from typing import Any

from phylobot.homologs import lookup, remote
from phylobot.homologs.remote import SearchError

RCSB, ALPHAFOLD, UPLOAD = "rcsb", "alphafold", "upload"
SOURCES = {RCSB: "PDB", ALPHAFOLD: "AlphaFold DB", UPLOAD: "Uploaded"}
PDB_ID = re.compile(r"[0-9][A-Za-z0-9]{3}")
# UniProt's accession format (https://www.uniprot.org/help/accession_numbers).
ACCESSION = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")
_IN_NAME = re.compile(rf"(?<![A-Za-z0-9])({ACCESSION.pattern})(?![A-Za-z0-9])")
RCSB_FILE = "https://files.rcsb.org/download/{}.{}"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{}"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{}"
ALPHAFOLD_FILES = "https://alphafold.ebi.ac.uk/files/"
FIELDS = (
    "accession,id,protein_name,gene_names,organism_name,length,reviewed,xref_pdb,xref_alphafolddb"
)
MAX_RESULTS = 10
MAX_ENTRIES = 40  # PDB entries listed per protein (ABL1 has 85)


@dataclass(frozen=True)
class Fetched:
    source: str
    accession: str
    title: str
    data: bytes


def clean_id(text: str) -> tuple[str, str]:
    """(source, identifier) for what the user typed: a PDB ID or a UniProt accession,
    which means its AlphaFold model. SearchError for anything else."""
    text = text.strip()
    if text.upper().startswith("AF-"):  # AF-P00519-F1, as AlphaFold DB writes it
        text = text[3:].split("-")[0]
    if PDB_ID.fullmatch(text):
        return RCSB, text.upper()
    if ACCESSION.fullmatch(text.upper()):
        return ALPHAFOLD, text.upper()
    raise SearchError("enter a PDB ID (4 characters, like 2SRC) or a UniProt accession (P12931)")


def clean_ids(text: str, most: int) -> list[tuple[str, str]]:
    """Several identifiers typed together (commas, spaces or new lines between them), each
    checked as `clean_id` does; duplicates once. SearchError names the first bad one."""
    found: list[tuple[str, str]] = []
    for part in re.split(r"[\s,;]+", text.strip()):
        if not part:
            continue
        try:
            one = clean_id(part)
        except SearchError as exc:
            raise SearchError(f"“{part[:20]}” is not a PDB ID or a UniProt accession") from exc
        if one not in found:
            found.append(one)
    if not found:
        raise SearchError("enter a PDB ID (like 2SRC) or a UniProt accession (P12931)")
    if len(found) > most:
        raise SearchError(f"at most {most} structures at a time")
    return found


def accessions_in(names: list[str]) -> dict[str, str]:
    """{sequence name: UniProt accession} for the names that carry one
    (Homo_sapiens_P00519): their AlphaFold models can be attached with one click."""
    found = {}
    for name in names:
        match = _IN_NAME.search(name.replace("_", " ").replace("|", " "))
        if match:
            found[name] = match.group(1)
    return found


def _json(url: str) -> Any:
    try:
        return json.loads(remote.request(url))
    except (json.JSONDecodeError, ValueError) as exc:
        raise SearchError("the structure database's answer could not be read") from exc


def fetch_rcsb(pdb_id: str) -> Fetched:
    """The entry's PDB-format file, or its mmCIF where RCSB has no PDB-format one (the
    entries too large for the format)."""
    if not PDB_ID.fullmatch(pdb_id):
        raise SearchError("not a PDB ID")
    pdb_id = pdb_id.upper()
    try:
        data = remote.request(RCSB_FILE.format(pdb_id, "pdb"))
    except SearchError as exc:
        if "404" not in str(exc):
            raise
        try:
            data = remote.request(RCSB_FILE.format(pdb_id, "cif"))
        except SearchError as again:
            if "404" in str(again):
                raise SearchError(f"RCSB has no entry {pdb_id}") from again
            raise
    title = ""
    try:  # the title is a nicety: the file's own title does when this fails
        entry = _json(RCSB_ENTRY.format(pdb_id))
        title = str((entry.get("struct") or {}).get("title") or "")[:300]
    except SearchError:
        pass
    return Fetched(RCSB, pdb_id, title, data)


def fetch_alphafold(accession: str) -> Fetched:
    if not ACCESSION.fullmatch(accession):
        raise SearchError("not a UniProt accession")
    try:
        answer = _json(ALPHAFOLD_API.format(accession))
    except SearchError as exc:
        if "404" in str(exc):
            raise SearchError(f"AlphaFold DB has no model for {accession}") from exc
        raise
    first = answer[0] if isinstance(answer, list) and answer else {}
    url = str(first.get("pdbUrl") or "") if isinstance(first, dict) else ""
    if not url.startswith(ALPHAFOLD_FILES) or not re.fullmatch(r"[A-Za-z0-9._/-]+", url[8:]):
        raise SearchError(f"AlphaFold DB has no model for {accession}")
    title = f"AlphaFold model of {first.get('uniprotDescription') or accession}"
    if organism := first.get("organismScientificName"):
        title += f" ({organism})"
    return Fetched(ALPHAFOLD, accession, title[:300], remote.request(url))


def fetch(source: str, accession: str) -> Fetched:
    if source == RCSB:
        return fetch_rcsb(accession)
    if source == ALPHAFOLD:
        return fetch_alphafold(accession)
    raise SearchError("unknown structure source")


@dataclass(frozen=True)
class Entry:
    """One PDB entry of a protein, as UniProt lists it."""

    pdb_id: str
    method: str
    resolution: str  # "1.50 A", or "" (NMR, or not given)
    chains: str  # "A/B=83-533"

    @property
    def url(self) -> str:
        return f"https://www.rcsb.org/structure/{self.pdb_id}"


@dataclass(frozen=True)
class Protein:
    accession: str
    entry_id: str
    protein: str
    gene: str
    organism: str
    length: int
    reviewed: bool
    alphafold: bool
    entries: list[Entry] = field(default_factory=list)
    more: int = 0  # PDB entries beyond those listed

    @property
    def url(self) -> str:
        return f"https://www.uniprot.org/uniprotkb/{self.accession}/entry"


def _resolution(entry: Entry) -> float:
    try:
        return float(entry.resolution.split()[0])
    except ValueError, IndexError:
        return 99.0


def parse(raw: bytes) -> list[Protein]:
    """UniProt's answer as proteins with their structures; those with any come first,
    and a protein's entries are best resolution first."""
    try:
        payload = json.loads(raw)
    except (json.JSONDecodeError, ValueError) as exc:
        raise SearchError("UniProt returned an answer that could not be read") from exc
    found = []
    for item in payload.get("results") or []:
        if not isinstance(item, dict):
            continue
        accession = str(item.get("primaryAccession") or "")
        if not ACCESSION.fullmatch(accession):
            continue
        entries, alphafold = [], False
        for xref in item.get("uniProtKBCrossReferences") or []:
            if not isinstance(xref, dict):
                continue
            database, ident = xref.get("database"), str(xref.get("id") or "")
            if database == "AlphaFoldDB":
                alphafold = True
            elif database == "PDB" and PDB_ID.fullmatch(ident):
                props = {
                    str(p.get("key")): str(p.get("value") or "")
                    for p in xref.get("properties") or []
                    if isinstance(p, dict)
                }
                resolution = props.get("Resolution", "")
                entries.append(
                    Entry(
                        ident.upper(),
                        props.get("Method", "")[:40],
                        "" if resolution in ("", "-") else resolution[:20],
                        props.get("Chains", "")[:80],
                    )
                )
        entries.sort(key=_resolution)
        organism = item.get("organism") or {}
        sequence = item.get("sequence") or {}
        found.append(
            Protein(
                accession=accession,
                entry_id=str(item.get("uniProtkbId") or ""),
                protein=lookup._protein_name(item.get("proteinDescription") or {}),
                gene=lookup._gene_name(item.get("genes") or []),
                organism=str(organism.get("scientificName") or ""),
                length=int(sequence.get("length") or 0),
                reviewed="unreviewed" not in str(item.get("entryType") or "").lower(),
                alphafold=alphafold,
                entries=entries[:MAX_ENTRIES],
                more=max(0, len(entries) - MAX_ENTRIES),
            )
        )
    found.sort(key=lambda p: (not p.entries, not p.reviewed))
    return found


def search(name: str, species: str) -> list[Protein]:
    """Proteins named `name` (in `species`, if given) with their PDB entries and
    AlphaFold models. The terms are cleaned as the sequence lookup's are."""
    params = {
        "query": lookup.build_query(name, species, reviewed_only=False),
        "format": "json",
        "fields": FIELDS,
        "size": str(MAX_RESULTS),
    }
    url = f"{lookup.BASE_URL}?{urllib.parse.urlencode(params)}"
    return parse(remote.request(url))
