"""Species from sequence accessions, looked up at NCBI and UniProt.

Most FASTA files come from NCBI or UniProt with the accession kept as the sequence name
(`XP_012345.1`, `sp|P32581|IME2_YEAST`, `UniRef90_P32581`). `find` recognises one;
`Client.look_up` asks the database it belongs to for its organism and NCBI taxid:

- NCBI E-utilities `esummary` (db protein or nuccore), which answers each accession's
  `organism` and `taxid`; NCBI asks for a tool name and contact email, and at most
  three requests a second without an API key;
- UniProt's REST `uniprotkb/accessions`, which answers organism name and id (and
  silently leaves out accessions it does not know);
- Ensembl's REST `lookup/id`, for Ensembl ids whose species prefix PhyloBot does not
  know offline (`headers.ENSEMBL_SPECIES`).

RefSeq's WP_ records are non-redundant proteins: one record for an identical protein in
many genomes, whose organism NCBI gives as their common taxon ("Enterobacteriaceae").
For those, NCBI's Identical Protein Groups list every genome carrying the protein, and
`Found.choices` ranks the species by genomes, for the user to pick one: identical
proteins across species often mean horizontal transfer, which reconciliation under
duplication and loss cannot model, so the choice is theirs.

Only the accessions are sent, never sequences, and only after the user ticks the New
Job page's own consent box for it (owner decision 2026-10-17: separate from Open Tree's).
Answers are bounded in size and time and cached on disk like Open Tree's, so editing
the table does not ask again. Every call goes through a `Transport`, so tests replay
recorded answers.
"""

import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

from phylobot.species.headers import ENSEMBL, clean_taxon
from phylobot.species.known import binomial
from phylobot.species.opentree import OpenTreeError, read_bounded

NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/accessions"
ENSEMBL_URL = "https://rest.ensembl.org/lookup/id"
IPG_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TIMEOUT = 60
MIN_GAP_SECONDS = 0.4  # NCBI: at most 3 requests a second without an API key
NCBI_BATCH = 200
UNIPROT_BATCH = 100
ENSEMBL_BATCH = 500
IPG_BATCH = 10  # groups can hold thousands of genomes each
MAX_CHOICES = 8
MAX_ACCESSIONS = 2000  # per lookup; more are left for the user
CACHE_DAYS = 30

# UniProtKB accession (UniProt's documented format), alone or inside sp|..|, tr|..| or
# UniRef90_...; NCBI RefSeq (XP_, NP_, WP_, YP_, AP_; XM_, NM_ for coding DNA) and
# GenBank protein (three letters and five or seven digits) or versioned nucleotide ids.
UNIPROT = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")
REFSEQ_PROTEIN = re.compile(r"[ANWXY]P_\d{6,9}(?:\.\d{1,3})?")
REFSEQ_NUCLEOTIDE = re.compile(r"[NX][MR]_\d{6,9}(?:\.\d{1,3})?")
GENBANK_PROTEIN = re.compile(r"[A-Z]{3}(?:\d{5}|\d{7})(?:\.\d{1,3})?")
GENBANK_NUCLEOTIDE = re.compile(r"[A-Z]{1,2}\d{5,8}\.\d{1,3}")
UNIPROT_WRAPPED = re.compile(r"^(?:sp|tr)\|([A-Z0-9]{6,10})\|")
UNIREF = re.compile(r"^UniRef(?:50|90|100)_([A-Z0-9]{6,10})$")


class AccessionError(Exception):
    pass


@dataclass(frozen=True)
class Accession:
    database: str  # "uniprot", "protein" or "nuccore" (NCBI), "ensembl"
    id: str


@dataclass(frozen=True)
class Found:
    accession: Accession
    organism: str
    taxid: int | None
    # for a protein found in several species: (species, genomes), most genomes first
    choices: tuple[tuple[str, int], ...] = ()
    genomes: int = 0


def find(name: str) -> Accession | None:
    """The accession a sequence name carries, if it is one we can look up."""
    if match := ENSEMBL.search(name):
        return Accession("ensembl", match[0].split(".")[0])
    if match := UNIPROT_WRAPPED.match(name) or UNIREF.match(name):
        return Accession("uniprot", match[1]) if UNIPROT.fullmatch(match[1]) else None
    # NCBI's "gi|123|ref|XP_1.1|" and similar: try each field
    for part in [name, *name.split("|")]:
        part = part.strip()
        if REFSEQ_PROTEIN.fullmatch(part) or GENBANK_PROTEIN.fullmatch(part):
            return Accession("protein", part)
        if REFSEQ_NUCLEOTIDE.fullmatch(part) or GENBANK_NUCLEOTIDE.fullmatch(part):
            return Accession("nuccore", part)
        if UNIPROT.fullmatch(part):
            return Accession("uniprot", part)
    return None


# (url, form data, or None for GET, or a JSON body as str) -> response bytes; raises
# AccessionError
Transport = Callable[[str, dict[str, str] | str | None], bytes]


class Http:
    """The real transport, paced like NCBI asks."""

    def __init__(self, min_gap: float = MIN_GAP_SECONDS):
        self.min_gap = min_gap
        self.last = 0.0

    def __call__(self, url: str, data: dict[str, str] | str | None) -> bytes:
        wait = self.last + self.min_gap - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        headers = {"User-Agent": "PhyloBot"}
        if isinstance(data, str):  # a JSON body
            body: bytes | None = data.encode()
            headers |= {"Content-Type": "application/json", "Accept": "application/json"}
        else:
            body = urllib.parse.urlencode(data).encode() if data is not None else None
        host = urllib.parse.urlsplit(url).netloc
        req = urllib.request.Request(url, data=body, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as response:
                return read_bounded(response)
        except urllib.error.HTTPError as exc:
            raise AccessionError(f"{host} answered {exc.code}") from exc
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            raise AccessionError(f"could not reach {host}: {exc}") from exc
        except OpenTreeError as exc:  # read_bounded's size and time limits
            raise AccessionError(str(exc).replace("Open Tree of Life", "service")) from exc
        finally:
            self.last = time.monotonic()


class Cache:
    """One accession's answer per file, kept `days` days (misses are not cached)."""

    def __init__(self, directory: Path, days: float = CACHE_DAYS):
        self.directory = directory
        self.seconds = days * 86400

    def path(self, accession: Accession) -> Path:
        key = f"{accession.database}:{accession.id}".encode()
        return self.directory / f"acc-{hashlib.sha256(key).hexdigest()}.json"

    def get(self, accession: Accession) -> Found | None:
        try:
            path = self.path(accession)
            if time.time() - path.stat().st_mtime > self.seconds:
                return None
            data = json.loads(path.read_text())
            taxid = None if data["taxid"] is None else int(data["taxid"])
            choices = tuple((str(name), int(n)) for name, n in data.get("choices") or [])
            genomes = int(data.get("genomes") or 0)
            return Found(accession, str(data["organism"]), taxid, choices, genomes)
        except OSError, ValueError, KeyError, TypeError:
            return None

    def put(self, found: Found) -> None:
        self.directory.mkdir(parents=True, exist_ok=True)
        path = self.path(found.accession)
        partial = path.with_suffix(".part")
        partial.write_text(json.dumps({
            "organism": found.organism, "taxid": found.taxid,
            "choices": [list(c) for c in found.choices], "genomes": found.genomes,
        }))  # fmt: skip
        partial.replace(path)


def chunks[T](items: list[T], size: int) -> Iterable[list[T]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def species_of(organism: str) -> str | None:
    """'Lactococcus lactis subsp. lactis' -> 'Lactococcus lactis'; None for 'Lactococcus
    sp.' and anything else that names no species."""
    name = binomial(organism)
    return name if name and not name.endswith((" sp", " sp.")) else None


def bare(accession: str) -> str:
    return accession.split(".")[0]


class Client:
    def __init__(self, email: str, transport: Transport | None = None, cache: Cache | None = None):
        self.email = email
        self.transport = transport or Http()
        self.cache = cache

    def look_up(self, accessions: Iterable[Accession]) -> dict[Accession, Found]:
        """Each accession's organism and taxid; those the databases do not know are
        absent. AccessionError if a database cannot be asked."""
        wanted = list(dict.fromkeys(accessions))[:MAX_ACCESSIONS]
        found: dict[Accession, Found] = {}
        missing = []
        for accession in wanted:
            cached = self.cache.get(accession) if self.cache else None
            if cached is not None:
                found[accession] = cached
            else:
                missing.append(accession)
        for database in ("protein", "nuccore"):
            batch = [a for a in missing if a.database == database]
            for part in chunks(batch, NCBI_BATCH):
                found |= self.ncbi(database, part)
        uniprot = [a for a in missing if a.database == "uniprot"]
        for part in chunks(uniprot, UNIPROT_BATCH):
            found |= self.uniprot(part)
        ensembl = [a for a in missing if a.database == "ensembl"]
        for part in chunks(ensembl, ENSEMBL_BATCH):
            found |= self.ensembl(part)
        # proteins NCBI places above species: which species carry them
        shared = [f for a, f in found.items() if a in missing and a.database == "protein"
                  and binomial(f.organism) is None]  # fmt: skip
        for group in chunks(shared, IPG_BATCH):
            found |= self.groups(group)
        if self.cache:
            for accession in missing:
                if accession in found:
                    self.cache.put(found[accession])
        return found

    def ncbi(self, database: str, accessions: list[Accession]) -> dict[Accession, Found]:
        data = {
            "db": database,
            "id": ",".join(a.id for a in accessions),
            "retmode": "json",
            "tool": "phylobot",
            "email": self.email,
        }
        try:
            answer = json.loads(self.transport(NCBI_URL, data))
        except ValueError as exc:
            raise AccessionError("NCBI did not answer in JSON") from exc
        result = answer.get("result") if isinstance(answer, dict) else None
        if not isinstance(result, dict):
            return {}  # every id invalid: NCBI answers only an error
        by_id = {bare(a.id): a for a in accessions}
        found = {}
        for uid in result.get("uids") or []:
            entry = result.get(str(uid))
            if not isinstance(entry, dict):
                continue
            accession = by_id.get(bare(str(entry.get("accessionversion") or entry.get("caption"))))
            taxid, organism = entry.get("taxid"), clean_taxon(str(entry.get("organism") or ""))
            if accession is not None and isinstance(taxid, int) and taxid > 0 and organism:
                found[accession] = Found(accession, organism, taxid)
        return found

    def ensembl(self, accessions: list[Accession]) -> dict[Accession, Found]:
        try:
            body = json.dumps({"ids": [a.id for a in accessions]})
            answer = json.loads(self.transport(ENSEMBL_URL, body))
        except ValueError as exc:
            raise AccessionError("Ensembl did not answer in JSON") from exc
        found = {}
        for accession in accessions:
            entry = answer.get(accession.id) if isinstance(answer, dict) else None
            species = entry.get("species") if isinstance(entry, dict) else None
            if isinstance(species, str) and re.fullmatch(r"[a-z]+_[a-z_]+", species):
                organism = clean_taxon(species.replace("_", " ").capitalize())
                found[accession] = Found(accession, organism, None)
        return found

    def groups(self, shared: list[Found]) -> dict[Accession, Found]:
        """NCBI's Identical Protein Groups for proteins found in several genomes: the
        species carrying each, by number of genomes."""
        data = {
            "db": "ipg", "id": ",".join(f.accession.id for f in shared), "rettype": "ipg",
            "retmode": "text", "tool": "phylobot", "email": self.email,
        }  # fmt: skip
        text = self.transport(IPG_URL, data).decode(errors="replace")
        rows = [line.split("\t") for line in text.splitlines()]
        rows = [r for r in rows if len(r) >= 9 and r[0] != "Id"]
        group_of: dict[str, str] = {}  # our accession -> its group id
        wanted = {bare(f.accession.id): f for f in shared}
        for r in rows:
            if bare(r[6]) in wanted:
                group_of[bare(r[6])] = r[0]
        out = {}
        for key, f in wanted.items():
            counts: dict[str, int] = {}
            genomes = 0
            for r in rows:
                if r[0] == group_of.get(key):
                    genomes += 1
                    if species := species_of(clean_taxon(r[8])):
                        counts[species] = counts.get(species, 0) + 1
            ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
            if len(ranked) == 1:  # one species after all (strains and subspecies)
                out[f.accession] = Found(f.accession, ranked[0][0], None, genomes=genomes)
            elif ranked:
                out[f.accession] = Found(f.accession, f.organism, f.taxid,
                                         tuple(ranked[:MAX_CHOICES]), genomes)  # fmt: skip
        return out

    def uniprot(self, accessions: list[Accession]) -> dict[Accession, Found]:
        query = urllib.parse.urlencode({
            "accessions": ",".join(a.id for a in accessions),
            "fields": "accession,organism_name,organism_id",
            "format": "tsv",
        })  # fmt: skip
        text = self.transport(f"{UNIPROT_URL}?{query}", None).decode(errors="replace")
        by_id = {a.id: a for a in accessions}
        found = {}
        for line in text.splitlines()[1:]:
            cells = line.split("\t")
            if len(cells) != 3 or cells[0] not in by_id or not cells[2].isdigit():
                continue
            organism = clean_taxon(cells[1])
            if organism:
                found[by_id[cells[0]]] = Found(by_id[cells[0]], organism, int(cells[2]))
        return found


def describe(found: Found) -> str:
    where = "UniProt" if found.accession.database == "uniprot" else "NCBI"
    return f"from {where} ({found.accession.id})"
