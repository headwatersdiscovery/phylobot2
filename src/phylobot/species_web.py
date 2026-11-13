"""Species for reconciliation on the New Job page.

The page's species table and Submit both come here (`gather`, then `look_up`). Each
sequence's species comes from, each overriding the one before:

1. its header (`species.resolve.detect`, including UniProt species codes, offline);
2. its accession, looked up at NCBI or UniProt once the user ticks that consent box
   (`accession_species`; only the accessions are sent);
3. a taxid a homolog search carried (`hints`);
4. the user's list, pasted or uploaded (`bulk`);
5. what the user typed in the table (`edits`).

Nothing goes to the Open Tree of Life until its own consent box is ticked (owner
decisions 2026-10-17: the two boxes are separate); an uploaded species tree needs no
consent, since it is matched here. Only species names, taxids and accessions are ever
sent, never sequences.

Reconciliation never blocks a job: `look_up` says what is still missing, and Submit
then runs the job without reconciliation and records why (`Lookup.skipped`).
"""

import csv
import re
import threading
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from tempfile import gettempdir
from typing import Any

from phylobot.pipeline.validate import InputError, parse_fasta
from phylobot.species import accessions, known, resolve
from phylobot.species.accessions import Accession, AccessionError, Found
from phylobot.species.opentree import Cache, Client, Http, OpenTreeError, Transport
from phylobot.species.speciestree import Node, SpeciesTreeError
from phylobot.species.table import NOT_IN_TREE, OK, Row, TableError, clean, number

CACHE_DIR = Path(gettempdir()) / "phylobot-opentree"
SYNTHESIS_SECONDS = 3600  # how long one client (and its synthesis release) is reused
MAX_HINTS = 100_000
# Statuses of the whole lookup, for the page.
NO_SEQUENCES = "no sequences"
NEEDS_CONSENT = "needs consent"
INCOMPLETE = "incomplete"
FAILED = "failed"
READY = "ready"

MAX_LIST_BYTES = 2_000_000  # a pasted or uploaded species list
# Where each row's species came from, for the table's "From" column.
SOURCES = {
    "uniprot": "UniProt header", "uniref": "UniRef header", "ncbi": "NCBI header",
    "genbank": "GenBank header", "uniprot code": "UniProt code", "name": "name",
    "taxid prefix": "id's taxid", "ensembl id": "Ensembl id", "model id": "gene id",
    "guess": "guess", "search": "homolog search", "list": "your list", "you": "you",
    "lookup ncbi": "NCBI lookup", "lookup uniprot": "UniProt lookup",
    "lookup ensembl": "Ensembl lookup",
}  # fmt: skip

# The transports for new clients: the paced HTTP ones, or recordings in tests.
make_transport: Callable[[], Transport] = Http
make_accession_transport: Callable[[], accessions.Transport] = accessions.Http
_client: tuple[float, Client] | None = None
_lock = threading.Lock()


def client() -> Client:
    """One Open Tree client per process, renewed hourly, its answers cached on disk
    (keyed by the synthesis release), so looking the same species up again as the user
    edits the table costs one request."""
    global _client
    with _lock:
        now = time.monotonic()
        if _client is None or now - _client[0] > SYNTHESIS_SECONDS:
            _client = (now, Client(make_transport(), cache=Cache(CACHE_DIR)))
        return _client[1]


def reset_client() -> None:
    global _client
    with _lock:
        _client = None


@dataclass
class Lookup:
    state: str
    rows: list[Row] = field(default_factory=list)
    tree: Node | None = None
    provenance: dict[str, Any] = field(default_factory=dict)
    error: str = ""
    # rows the synthetic tree leaves out -> their parent taxon (name, rank), to offer
    parents: dict[str, tuple[str, str]] = field(default_factory=dict)

    @property
    def pending(self) -> list[Row]:
        return [row for row in self.rows if row.status != OK]

    @property
    def ready(self) -> bool:
        return self.state == READY

    @property
    def species(self) -> int:
        return len(self.tree.leaves()) if self.tree is not None else 0

    @property
    def polytomies(self) -> int:
        return len(self.tree.polytomies()) if self.tree is not None else 0

    @property
    def skipped(self) -> str:
        """Why a job with these species runs without reconciliation ("" if it doesn't)."""
        if self.state == NEEDS_CONSENT:
            return ("the species names were not sent to the Open Tree of Life (its box was "
                    "unticked) and no species tree was uploaded")  # fmt: skip
        if self.state == FAILED:
            return self.error
        if self.state == INCOMPLETE:
            if self.tree is None:
                return "fewer than two of the species are in the species tree"
            names = [row.sequence for row in self.pending]
            shown = ", ".join(names[:3]) + (f" and {len(names) - 3} more" if len(names) > 3 else "")
            return f"{len(names)} sequence(s) had no species in the species tree ({shown})"
        if self.state == NO_SEQUENCES:
            return "there were no sequences to read species from"
        return ""


def shown(row: Row) -> str:
    """What the table showed in a row's cells, unless the user had typed them (""): a
    cell still holding it is not the user's, so a new list or lookup may replace it."""
    return "" if row.source == "you" else f"{row.query}\t{row.taxid or ''}"


def edits(data: Mapping[str, Any]) -> dict[str, tuple[str, int | None]]:
    """The cells the user typed in the table: {sequence: (query, taxid)}. Each row posts
    its sequence name, query, taxid and what was shown (`species_for`, `species_query`,
    `species_taxid`, `species_shown`); a cell unchanged from what was shown is not an
    edit. Everything is cleaned as hostile, and a malformed cell is dropped."""
    get = getattr(data, "getlist", None)
    if get is None:
        return {}
    names, queries, taxids = get("species_for"), get("species_query"), get("species_taxid")
    if not (len(names) == len(queries) == len(taxids)):
        return {}
    was = get("species_shown")
    if len(was) != len(names):
        was = [""] * len(names)  # no record of what was shown: every cell counts
    out: dict[str, tuple[str, int | None]] = {}
    for name, query, taxid, before in zip(names, queries, taxids, was, strict=True):
        try:
            cell = (clean(query, "species"), number(taxid, "taxid"))
            sequence = clean(name, "sequence")
        except TableError:
            continue
        if before and f"{cell[0]}\t{cell[1] or ''}" == before:
            continue
        out[sequence] = cell
    return out


def bulk(text: str) -> dict[str, tuple[str, int | None]]:
    """A species list, pasted or uploaded: one sequence per line, its species (a name,
    or a taxid) and optionally a taxid, separated by tabs, commas or semicolons (each
    line its own). A first line naming the columns (sequence, species or query or
    organism, taxid) is read as a header, so the table's own download comes back in.
    Bad lines are skipped."""
    lines = [line for line in text.splitlines()[:MAX_HINTS] if line.strip()]
    if not lines:
        return {}

    def split(line: str) -> list[str]:  # each line by its own separator (pastes mix them)
        delimiter = next((d for d in ("\t", ",", ";") if d in line), "\t")
        return next(csv.reader([line], delimiter=delimiter))

    rows = [split(line) for line in lines]
    head = [cell.strip().lower() for cell in rows[0]]
    columns = (0, 1, 2)
    if "sequence" in head:
        name_at = head.index("sequence")
        names = ("species", "query", "organism")
        species_at = next((head.index(c) for c in names if c in head), -1)
        taxid_at = head.index("taxid") if "taxid" in head else -1
        columns, rows = (name_at, species_at, taxid_at), rows[1:]
    out: dict[str, tuple[str, int | None]] = {}
    for cells in rows:

        def cell(at: int, cells: list[str] = cells) -> str:
            return cells[at] if 0 <= at < len(cells) else ""

        try:
            name = clean(cell(columns[0]), "sequence")
            species = clean(cell(columns[1]), "species")
            taxid = number(cell(columns[2]), "taxid")
        except TableError:
            continue
        if species.isdigit() and taxid is None:
            species, taxid = "", int(species) if len(species) <= 10 else None
        if name and (species or taxid is not None):
            out[name] = (species, taxid)
    return out


def accession_client() -> accessions.Client:
    from phylobot import site_settings  # the web tier's contact email, for NCBI

    return accessions.Client(
        site_settings.search_email(), make_accession_transport(), accessions.Cache(CACHE_DIR)
    )


@dataclass
class Accessions:
    """The sequences' accessions and what NCBI and UniProt said of them."""

    pending: dict[str, Accession] = field(default_factory=dict)  # sequence -> to look up
    found: dict[str, Found] = field(default_factory=dict)  # sequence -> organism, taxid
    asked: bool = False
    error: str = ""


def accession_species(
    fasta: bytes, consent: bool, get_client: Callable[[], accessions.Client] = accession_client
) -> Accessions:
    """Accessions of the sequences whose headers give no taxid, looked up only if the
    user consented."""
    result = Accessions()
    for row in rows_for(fasta):
        if row.taxid is None and (accession := accessions.find(row.sequence)):
            result.pending[row.sequence] = accession
    if not consent or not result.pending:
        return result
    result.asked = True
    try:
        answers = get_client().look_up(result.pending.values())
    except AccessionError as exc:
        result.error = f"the accessions could not be looked up: {exc}"
        return result
    result.found = {name: answers[a] for name, a in result.pending.items() if a in answers}
    return result


def records(fasta: bytes) -> list[tuple[str, str]]:
    try:
        parsed = parse_fasta(fasta.decode("utf-8", errors="replace"), "sequences")
    except InputError:
        return []
    return [(name, description) for name, description, _ in parsed]


def hints(text: str) -> dict[str, int]:
    """Taxids carried from a homolog search: "name<TAB>taxid" lines (malformed lines
    are ignored)."""
    out: dict[str, int] = {}
    for line in text.splitlines()[:MAX_HINTS]:
        name, _, taxid = line.partition("\t")
        try:
            value = number(taxid, "taxid")
            name = clean(name, "sequence")
        except TableError:
            continue
        if name and value is not None:
            out[name] = value
    return out


def rows_for(
    fasta: bytes,
    typed: Mapping[str, tuple[str, int | None]] | None = None,
    carried: Mapping[str, int] | None = None,
    listed: Mapping[str, tuple[str, int | None]] | None = None,
    found: Mapping[str, Found] | None = None,
) -> list[Row]:
    """One row per sequence: the header's species, its accession's (`found`), a carried
    taxid, the user's list, then the user's own cells, each overriding the one before.
    [] if the FASTA can't be read."""
    typed, carried, listed, found = typed or {}, carried or {}, listed or {}, found or {}
    rows = []
    for row in resolve.detect(records(fasta)):
        if row.sequence in found and row.taxid is None:
            answer = found[row.sequence]
            where = answer.accession.database
            source = f"lookup {'ncbi' if where in ('protein', 'nuccore') else where}"
            note = ""
            if answer.choices:
                more = "+" if len(answer.choices) == accessions.MAX_CHOICES else ""
                note = (
                    f"the identical protein is in {answer.genomes} genomes of "
                    f"{len(answer.choices)}{more} species: choose one (identical proteins "
                    "across species can mean horizontal transfer)"
                )
            row = replace(row, query=answer.organism, taxid=answer.taxid, source=source, note=note)
        if row.sequence in carried and row.taxid is None:
            row = replace(row, taxid=carried[row.sequence], source=row.source or "search")
        if row.sequence in listed:
            query, taxid = listed[row.sequence]
            row = replace(row, query=query, taxid=taxid, source="list", note="")
        if row.sequence in typed:
            query, taxid = typed[row.sequence]
            if (query, taxid) != (row.query, row.taxid):
                row = replace(row, query=query, taxid=taxid, source="you", note="")
        rows.append(row)
    return rows


def enrich(row: Row) -> Row:
    """What the bundled list can add offline: a common name becomes the scientific one
    ("human" -> Homo sapiens), and a bare taxid gets its name, so that Open Tree matches
    it by name in one batch instead of one request per taxid."""
    if row.query and (found := known.common(row.query)) and not known.species(row.query):
        note = f"{row.query!r} is {found.name}"
        return replace(row, query=found.name, taxid=row.taxid or found.taxid, note=note)
    if not row.query and row.taxid is not None and (name := known.by_taxid(row.taxid)):
        return replace(row, query=resolve_name(name))
    return row


def resolve_name(name: str) -> str:
    return known.binomial(name) or name


@dataclass
class Gathered:
    rows: list[Row]
    accessions: Accessions


def gather(
    fasta: bytes,
    data: Mapping[str, Any],
    listed_text: str = "",
    get_accession_client: Callable[[], accessions.Client] = accession_client,
) -> Gathered:
    """Every source of species for a form's sequences, before resolution (`look_up`):
    the posted form `data` gives the table's cells, the carried taxids, the pasted list
    (plus `listed_text`, an uploaded one) and the accession consent."""
    found = accession_species(
        fasta, data.get("accession_consent") in ("on", True), get_accession_client
    )
    listed = bulk(str(data.get("species_paste") or "")[:MAX_LIST_BYTES])
    listed |= bulk(listed_text[:MAX_LIST_BYTES])
    rows = rows_for(
        fasta, edits(data), hints(str(data.get("species_hints") or "")), listed, found.found
    )
    return Gathered([enrich(row) for row in rows], found)


PREFIX = re.compile(r"^([A-Za-z][A-Za-z0-9]{1,19})[._|-]")


def prefix(name: str) -> str | None:
    """A name's lab prefix ("Scer" in "Scer_IME2"); None for accessions, whose prefixes
    (WP_, XP_, sp|, ENS...) are the database's, not the species'."""
    if accessions.find(name) is not None:
        return None
    match = PREFIX.match(name)
    return match[1] if match else None


@dataclass(frozen=True)
class Offer:
    """A species to offer every row named like one the user filled in."""

    prefix: str
    query: str
    taxid: int | None
    source: str  # the row it comes from
    rows: tuple[str, ...]  # the rows it would fill


def offers(rows: Sequence[Row], resolved: bool) -> list[Offer]:
    """Species to offer rows that lack one (or have only a guess), from a row the user
    filled in whose name has the same prefix ("Scer_IME2" -> "Scer_KSS1", "Scer_FUS3"),
    one offer per prefix.
    Offered, never applied: the page shows a button. `resolved`: statuses are real."""
    donors: dict[str, Row] = {}
    for row in rows:
        filled = row.source in ("you", "list") and (row.query or row.taxid is not None)
        if filled and (p := prefix(row.sequence)):
            donors.setdefault(p, row)
    takers: dict[str, list[str]] = {}
    for row in rows:
        # rows without a species, or only a guessed one: the user's own choice beats a guess
        lacking = (not row.query and row.taxid is None) or row.source == "guess"
        if resolved:
            lacking = lacking or row.status != OK
        p = prefix(row.sequence)
        if lacking and p and p in donors and donors[p].sequence != row.sequence:
            donor = donors[p]
            if (donor.query, donor.taxid) != (row.query, row.taxid):
                takers.setdefault(p, []).append(row.sequence)
    return [
        Offer(prefix, donors[prefix].query, donors[prefix].taxid, donors[prefix].sequence,
              tuple(names))
        for prefix, names in takers.items()
    ]  # fmt: skip


def look_up(
    rows: Sequence[Row], consent: bool, upload: str = "", get_client: Callable[[], Client] = client
) -> Lookup:
    """Resolve the rows: on the uploaded species tree if there is one (nothing is sent),
    else at the Open Tree of Life if the user consented, else not at all."""
    rows = list(rows)
    if not rows:
        return Lookup(NO_SEQUENCES)
    try:
        if upload.strip():
            result = resolve.with_upload(rows, upload)
        elif consent:
            result = resolve.with_open_tree(rows, get_client())
        else:
            return Lookup(NEEDS_CONSENT, rows)
    except SpeciesTreeError as exc:
        return Lookup(FAILED, rows, error=f"the species tree could not be used: {exc}")
    except OpenTreeError as exc:
        return Lookup(FAILED, rows, error=f"the Open Tree of Life could not be asked: {exc}")
    state = READY if result.complete else INCOMPLETE
    found = Lookup(state, result.rows, result.tree, result.provenance)
    if consent and not upload.strip():
        found.parents = parents(result.rows, get_client)
    return found


MAX_PARENTS = 50  # rows left out of the synthetic tree whose parent is looked up


def parents(rows: Sequence[Row], get_client: Callable[[], Client]) -> dict[str, tuple[str, str]]:
    """For each row whose taxon the synthetic tree leaves out, its parent taxon: the
    page offers it, as the nearest place in the tree (one cached request per taxon)."""
    out: dict[str, tuple[str, str]] = {}
    asked: dict[int, tuple[str, str] | None] = {}
    for row in rows:
        if row.status != NOT_IN_TREE or row.ott_id is None:
            continue
        if row.ott_id not in asked:
            if len(asked) >= MAX_PARENTS:
                break
            try:
                parent = get_client().parent(row.ott_id)
            except OpenTreeError:
                break
            asked[row.ott_id] = (parent.name, parent.rank) if parent else None
        if (found := asked[row.ott_id]) is not None:
            out[row.sequence] = found
    return out
