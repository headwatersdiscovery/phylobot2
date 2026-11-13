"""Full-length sequences of search hits.

BLAST reports only each hit's aligned stretch (for Ime2 often just the kinase domain),
but a job or a FASTA download needs whole proteins. After a search, the runner fetches
them by ID in batches:
- EBI hits (UniRef clusters): each cluster's representative sequence from UniProt's
  REST API. The representative's accession is the part of the UniRef ID after
  "UniRef90_": UniProtKB entries come from `/uniprotkb/accessions`, UniParc ones (UPI...)
  from `/uniparc/upis` (~1 s per 100 each). Any still missing are asked for by UniRef ID
  (`/uniref/stream`: same result, but its query can take a minute).
- EBI hits from Swiss-Prot: EBI names them by UniProtKB entry name (VATL1_YEAST), so
  they are asked for by name (`/uniprotkb/stream`, `id:` queries).
- NCBI hits: NCBI E-utilities efetch (~0.5 s per 200; at most 3 requests a second
  without an API key, so requests are spaced).
Only public database IDs are sent. IDs come from the search result (remote data), so
each is checked before use; sequences are checked like any FASTA.
"""

import re
import time
import urllib.parse

from phylobot.homologs import remote
from phylobot.homologs.remote import SearchError

UNIPROT = "https://rest.uniprot.org"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNIREF_BATCH = 100
NCBI_BATCH = 200
NCBI_GAP_SECONDS = 0.4  # E-utilities: at most 3 requests a second without a key
UNIREF_ID = re.compile(r"^UniRef(100|90|50)_[A-Za-z0-9]{1,40}(-[0-9]{1,4})?$")  # -2: isoform
ENTRY_NAME = re.compile(r"^[A-Z0-9]{1,16}_[A-Z0-9]{1,5}$")  # UniProtKB entry name: IME2_YEAST
NCBI_ID = re.compile(r"^[A-Za-z0-9_]{1,40}(\.[0-9]{1,5})?$")
MAX_RESPONSE = 50_000_000
ATTEMPTS = 3  # per batch; the services time out now and then
RESIDUES = set("ABCDEFGHIKLMNPQRSTUVWXYZ")


def parse_fasta(text: str) -> dict[str, str]:
    """{first word of each header: sequence}; sequences upper-cased, '*' dropped, records
    with anything but residue letters left out."""
    records: dict[str, list[str]] = {}
    name = ""
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split(" ", 1)[0]
            records[name] = []
        elif name and line:
            records[name].append(line)
    out = {}
    for key, parts in records.items():
        seq = "".join(parts).upper().rstrip("*")
        if seq and set(seq) <= RESIDUES:
            out[key] = seq
    return out


def _batches(ids: list[str], size: int) -> list[list[str]]:
    return [ids[i : i + size] for i in range(0, len(ids), size)]


def _get(url: str, params: dict[str, str], post: bool = False) -> str:
    """One request, retried; SearchError after ATTEMPTS failures."""
    for attempt in range(1, ATTEMPTS + 1):
        try:
            if post:
                content = remote.request(url, params)
            else:
                content = remote.request(url + "?" + urllib.parse.urlencode(params))
            break
        except SearchError:
            if attempt == ATTEMPTS:
                raise
            time.sleep(2 * attempt)
    if len(content) > MAX_RESPONSE:
        raise SearchError("the sequence service's answer is too large")
    return content.decode(errors="replace")


def _accession(header_word: str) -> str:
    """'sp|P32581|IME2_YEAST' -> 'P32581'; 'UPI001AF2DC2D' stays as it is."""
    parts = header_word.split("|")
    return parts[1] if len(parts) >= 3 else header_word


def fetch_uniref(ids: list[str], email: str = "") -> dict[str, str]:
    """{UniRef ID: representative's sequence}. IDs that are not UniRef IDs are skipped;
    IDs UniProt does not return are simply missing."""
    wanted = [i for i in dict.fromkeys(ids) if UNIREF_ID.fullmatch(i)]
    representative = {i: i.split("_", 1)[1] for i in wanted}
    by_accession: dict[str, str] = {}
    kb = [a for a in representative.values() if not a.startswith("UPI")]
    upi = [a for a in representative.values() if a.startswith("UPI")]
    for batch in _batches(kb, UNIREF_BATCH):
        text = _get(
            f"{UNIPROT}/uniprotkb/accessions", {"accessions": ",".join(batch), "format": "fasta"}
        )
        for word, seq in parse_fasta(text).items():
            by_accession.setdefault(_accession(word), seq)
    for batch in _batches(upi, UNIREF_BATCH):
        text = _get(f"{UNIPROT}/uniparc/upis", {"upis": ",".join(batch), "format": "fasta"})
        by_accession |= parse_fasta(text)
    found = {i: by_accession[a] for i, a in representative.items() if a in by_accession}
    missing = [i for i in wanted if i not in found]
    for batch in _batches(missing, UNIREF_BATCH):
        query = " OR ".join(f"id:{i}" for i in batch)
        found |= parse_fasta(
            _get(f"{UNIPROT}/uniref/stream", {"format": "fasta", "query": f"({query})"})
        )
    return {i: found[i] for i in wanted if i in found}


def fetch_entry_names(ids: list[str], email: str = "") -> dict[str, str]:
    """{UniProtKB entry name: sequence}, for Swiss-Prot hits. Other IDs are skipped."""
    wanted = [i for i in dict.fromkeys(ids) if ENTRY_NAME.fullmatch(i)]
    found: dict[str, str] = {}
    for batch in _batches(wanted, UNIREF_BATCH):
        query = " OR ".join(f"id:{i}" for i in batch)
        text = _get(f"{UNIPROT}/uniprotkb/stream", {"format": "fasta", "query": f"({query})"})
        for word, seq in parse_fasta(text).items():
            parts = word.split("|")  # sp|P25515|VATL1_YEAST
            found.setdefault(parts[2] if len(parts) >= 3 else word, seq)
    return {i: found[i] for i in wanted if i in found}


def fetch_uniprot(ids: list[str], email: str = "") -> dict[str, str]:
    """EBI hits, whichever UniProt database they came from: UniRef IDs or entry names."""
    return fetch_uniref(ids, email) | fetch_entry_names(ids, email)


def fetch_ncbi(ids: list[str], email: str = "") -> dict[str, str]:
    """{NCBI accession as in the search result: sequence}. NCBI answers with versioned
    accessions (NP_012429.1); they are matched back to the unversioned ones."""
    wanted = [i for i in dict.fromkeys(ids) if NCBI_ID.fullmatch(i)]
    found: dict[str, str] = {}
    for n, batch in enumerate(_batches(wanted, NCBI_BATCH)):
        if n:
            time.sleep(NCBI_GAP_SECONDS)
        params = {"db": "protein", "rettype": "fasta", "retmode": "text", "tool": "phylobot"}
        if email:
            params["email"] = email
        params["id"] = ",".join(batch)
        for key, seq in parse_fasta(_get(EFETCH_URL, params, post=True)).items():
            found.setdefault(key.split(".", 1)[0], seq)
            found.setdefault(key, seq)
    return {i: found[i] for i in wanted if i in found}


def to_fasta(sequences: dict[str, str]) -> str:
    """For storage: one record per ID, no descriptions."""
    return "".join(f">{key}\n{seq}\n" for key, seq in sequences.items())
