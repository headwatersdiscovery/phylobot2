"""species.tsv: the species of each sequence, one row per sequence.

    sequence  query  taxid  ott_id  species  source  status  note

`query` is the name looked up (from the header, or typed by the user), `taxid` the NCBI
taxid if known, `ott_id` the Open Tree taxon it resolved to, and `species` the name the
species tree uses for it: Open Tree's current name (so "Kluyveromyces waltii" becomes
"Lachancea waltii", and synonyms of one taxon become one species), or the leaf of an
uploaded species tree. `status` is "ok" or what still needs the user. Users edit this
file (`query`, `taxid`, `species`) and resolve again, so reading it treats every cell
as hostile.
"""

import csv
import io
from dataclasses import dataclass, fields

COLUMNS = ("sequence", "query", "taxid", "ott_id", "species", "source", "status", "note")
MAX_CELL = 200
MAX_ROWS = 100_000
OK = "ok"
# What a row can still need (status values other than "ok").
UNRESOLVED = "unresolved"  # no taxon found for the query
AMBIGUOUS = "ambiguous"  # several taxa match; the user must choose (by taxid or name)
NOT_IN_TREE = "not in tree"  # a taxon Open Tree's synthesis (or the upload) leaves out


class TableError(Exception):
    pass


@dataclass(frozen=True)
class Row:
    sequence: str
    query: str = ""
    taxid: int | None = None
    ott_id: int | None = None
    species: str = ""
    source: str = ""
    status: str = UNRESOLVED
    note: str = ""


def clean(text: str, what: str) -> str:
    text = " ".join(str(text).split())
    if any(ord(c) < 32 or ord(c) == 127 for c in text):
        raise TableError(f"{what} contains control characters")
    if len(text) > MAX_CELL:
        raise TableError(f"{what} is longer than {MAX_CELL} characters")
    return text


def number(text: str, what: str) -> int | None:
    text = text.strip()
    if not text:
        return None
    if not text.isdigit() or len(text) > 10:
        raise TableError(f"{what} must be a whole number, not {text[:20]!r}")
    return int(text)


def read(text: str) -> list[Row]:
    """species.tsv as rows; TableError on anything malformed. Only `sequence` is
    required; missing columns read as empty."""
    lines = text.splitlines()
    if not lines:
        raise TableError("species table is empty")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t", quoting=csv.QUOTE_NONE)
    if not reader.fieldnames or "sequence" not in reader.fieldnames:
        raise TableError("species table needs a 'sequence' column")
    unknown = set(reader.fieldnames) - set(COLUMNS)
    if unknown:
        raise TableError(f"species table has unknown columns: {', '.join(sorted(unknown))}")
    rows: list[Row] = []
    seen: set[str] = set()
    for line, record in enumerate(reader, 2):
        if line - 1 > MAX_ROWS:
            raise TableError(f"species table has more than {MAX_ROWS} rows")
        if None in record:
            raise TableError(f"species table, line {line}: more cells than columns")
        cells = {k: (v or "") for k, v in record.items()}
        where = f"species table, line {line}"
        sequence = clean(cells["sequence"], f"{where}: sequence")
        if not sequence:
            raise TableError(f"{where}: no sequence name")
        if sequence in seen:
            raise TableError(f"{where}: sequence {sequence!r} is listed twice")
        seen.add(sequence)
        rows.append(
            Row(
                sequence=sequence,
                query=clean(cells.get("query", ""), f"{where}: query"),
                taxid=number(cells.get("taxid", ""), f"{where}: taxid"),
                ott_id=number(cells.get("ott_id", ""), f"{where}: ott_id"),
                species=clean(cells.get("species", ""), f"{where}: species"),
                source=clean(cells.get("source", ""), f"{where}: source"),
                status=clean(cells.get("status", ""), f"{where}: status") or UNRESOLVED,
                note=clean(cells.get("note", ""), f"{where}: note"),
            )
        )
    return rows


def write(rows: list[Row]) -> str:
    out = ["\t".join(COLUMNS)]
    for row in rows:
        cells = []
        for column in fields(Row):
            value = getattr(row, column.name)
            cells.append("" if value is None else " ".join(str(value).split()))
        out.append("\t".join(cells))
    return "\n".join(out) + "\n"
