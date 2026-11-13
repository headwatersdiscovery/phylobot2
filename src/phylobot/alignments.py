"""The user's alignments ("My alignments") and the alignment editor's edits.

Records are `models.Alignment`: uploads, alignments saved from the editor, and every
alignment a job publishes (registered when the job succeeds, `register_job`, or the
first time the user's alignments are listed, `ensure_registered`, for older jobs).

The editor shows an alignment through `Edits`, a small state kept in the page's URL
(columns cut and sequences removed, both as positions in the stored alignment, and the
row order), so a view is shareable, survives reloads and is applied here, on the server,
from the stored rows: the browser sends positions, never sequences. "Save a copy" stores
the edited rows as a new alignment whose `edits` record how it came from its parent
(`provenance`), so residue numbers can be mapped back; "Save (overwrite)" replaces a
saved alignment's rows (never a job's), composing its record (`compose`) and bumping its
`version`, which every cache and copy names.

Uploads are hostile input: they go through the pipeline's own FASTA parser and the same
character and size rules as the New Job page (phylobot.pipeline.validate).
"""

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from django.contrib.auth.models import AbstractBaseUser, AnonymousUser
from django.db import IntegrityError, transaction
from django.db.models import Q, QuerySet
from django.db.models.functions import Coalesce
from django.utils import timezone

from phylobot import codons, storage
from phylobot.models import Alignment, Job, Limits
from phylobot.pipeline import validate

GAP = "-"
SEQTYPES = {"protein": "Protein", "codon": "Coding DNA (codons)"}
MIN_SEQUENCES = 2
MAX_NAME = 100
# Saved alignments per user (uploads and edits; a job's own are not counted). A
# placeholder until the owner decides storage limits.
MAX_SAVED = 200
MAX_RANGE_TOKENS = 10_000  # in one ranges parameter


class AlignmentError(ValueError):
    """A problem with an upload, worded for the user."""


@dataclass(frozen=True)
class Row:
    name: str
    description: str
    row: str  # aligned, upper case, gaps as '-'


def unit(seqtype: str) -> int:
    """Characters per site: 3 for a codon alignment."""
    return 3 if seqtype == "codon" else 1


def site_count(rows: Sequence[Row], seqtype: str) -> int:
    return len(rows[0].row) // unit(seqtype) if rows else 0


def clean_name(text: str) -> str:
    return " ".join(str(text).split())[:MAX_NAME]


def detect_seqtype(text: str) -> str:
    return validate.detect_seqtype(text)


def parse(text: str, seqtype: str, code: codons.GeneticCode, limits: Limits) -> list[Row]:
    """An uploaded alignment, checked: FASTA with unique names, rows of one length, the
    pipeline's characters (amino acids, or A/C/G/T/N for codons) and gaps ("-", or "."
    and "?", stored as "-"), whole-codon gaps for codons, and the user's limits on
    sequences and length."""
    what = "the alignment"
    try:
        records = validate.parse_fasta(text, what)
    except validate.InputError as exc:
        raise AlignmentError(str(exc)) from None
    if len(records) < MIN_SEQUENCES:
        raise AlignmentError(f"need at least {MIN_SEQUENCES} sequences, found {len(records)}")
    if len(records) > limits.max_sequences:
        raise AlignmentError(
            f"at most {limits.max_sequences} sequences allowed, found {len(records)}"
        )
    allowed = (codons.DNA if seqtype == "codon" else validate.AMINO | validate.AMBIGUOUS) | {GAP}
    rows = []
    for name, description, raw in records:
        row = raw.upper().replace(".", GAP).replace("?", GAP)  # "?": missing (PHYLIP)
        try:
            validate.check_residues(name, row, what, allowed)
        except validate.InputError as exc:
            hint = " (is it a protein alignment? choose Protein)" if seqtype == "codon" else ""
            raise AlignmentError(f"{exc}{hint}") from None
        if not row.strip(GAP):
            raise AlignmentError(f"sequence {name!r} is all gaps")
        rows.append(Row(name, description, row))
    lengths = sorted({len(r.row) for r in rows})
    if len(lengths) > 1:
        raise AlignmentError(
            f"the rows have different lengths ({', '.join(map(str, lengths[:5]))}"
            f"{'…' if len(lengths) > 5 else ''}): this is not an alignment"
        )
    if seqtype == "codon":
        for r in rows:
            try:
                codons.translate_row(r.row, code)
            except codons.CodonError as exc:
                raise AlignmentError(f"sequence {r.name!r}: {exc}") from None
    for r in rows:
        length = (len(r.row) - r.row.count(GAP)) // unit(seqtype)
        if length > limits.max_length:
            kind = "codons" if seqtype == "codon" else "residues"
            raise AlignmentError(
                f"sequence {r.name!r} is {length} {kind}; your limit is {limits.max_length}"
            )
    return rows


def read_rows(text: str) -> list[Row]:
    """Stored rows (checked when they were stored, or written by the pipeline)."""
    return [
        Row(name, description, raw.upper().replace(".", GAP))
        for name, description, raw in validate.parse_fasta(text, "alignment")
    ]


def to_fasta(rows: Iterable[Row]) -> str:
    return "".join(
        f">{r.name}{' ' + r.description if r.description else ''}\n{r.row}\n" for r in rows
    )


def normalise_upload(content: bytes) -> bytes:
    """An uploaded alignment as the pipeline reads it: PHYLIP becomes FASTA and "?"
    (missing data) becomes a gap, as IQ-TREE treats it. Anything else is left for the
    usual checks to report; a PHYLIP header the rows do not fit is an AlignmentError."""
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError:
        return content
    text = phylip_to_fasta(text) or text
    if "?" not in text:
        return text.encode()
    lines = text.splitlines(keepends=True)
    return "".join(
        line if line.lstrip().startswith(">") else line.replace("?", GAP) for line in lines
    ).encode()


PHYLIP_HEADER = re.compile(r"\s*(\d+)\s+(\d+)\s*$")
PHYLIP_NAME = 10  # strict PHYLIP: the name is the first ten characters, padded


def phylip_to_fasta(text: str) -> str | None:
    """A PHYLIP alignment as FASTA, or None if `text` does not start with a PHYLIP
    header ("<sequences> <columns>"). Sequential and interleaved, relaxed (the name
    ends at the first space) and strict (the first ten characters) are all read: every
    reading that gives each sequence exactly <columns> characters is tried, and if two
    of them disagree the text is refused as ambiguous. Characters are left to `parse`."""
    lines = [line for line in text.splitlines() if line.strip()]
    header = PHYLIP_HEADER.match(lines[0]) if lines else None
    if not header:
        return None
    count, columns = int(header.group(1)), int(header.group(2))
    shape = f"PHYLIP header says {count} sequences × {columns} columns"
    if count < 1 or columns < 1:
        raise AlignmentError(f"the {shape}")
    body = lines[1:]
    readings = {
        tuple(rows)
        for strict in (False, True)
        for interleaved in (True, False)
        if (rows := _read_phylip(body, count, columns, strict, interleaved)) is not None
    }
    if not readings:
        raise AlignmentError(
            f"the {shape}, but the sequences that follow do not fit it "
            "(in either sequential or interleaved layout)"
        )
    if len(readings) > 1:
        raise AlignmentError(
            "this PHYLIP alignment can be read more than one way (sequential or "
            "interleaved, names of up to 10 characters or longer): paste it as FASTA"
        )
    return "".join(f">{name}\n{row}\n" for name, row in readings.pop())


def _read_phylip(
    lines: list[str], count: int, columns: int, strict: bool, interleaved: bool
) -> list[tuple[str, str]] | None:
    """One reading of a PHYLIP body: (name, row) pairs, or None if it does not fit."""

    def named(line: str) -> tuple[str, str] | None:
        if strict:
            name, rest = line[:PHYLIP_NAME], line[PHYLIP_NAME:]
        else:
            name, rest = (line.split(None, 1) + [""])[:2]
        name = "_".join(name.split())  # a strict name may hold spaces
        return (name, "".join(rest.split())) if name else None

    if len(lines) < count:
        return None
    names: list[str] = []
    parts: list[list[str]] = []
    position = 0
    for _ in range(count):
        if position >= len(lines):
            return None
        first = named(lines[position])
        if first is None:
            return None
        names.append(first[0])
        parts.append([first[1]])
        position += 1
        if not interleaved:  # the row goes on over the next lines until it is full
            while sum(map(len, parts[-1])) < columns and position < len(lines):
                parts[-1].append("".join(lines[position].split()))
                position += 1
    if interleaved:  # later blocks: one line per sequence, in the same order, no names
        rest = lines[position:]
        if len(rest) % count:
            return None
        for i, line in enumerate(rest):
            parts[i % count].append("".join(line.split()))
    elif position != len(lines):
        return None
    rows = ["".join(p) for p in parts]
    if any(len(row) != columns for row in rows):
        return None
    return list(zip(names, rows, strict=True))


def amino_rows(rows: Sequence[Row], seqtype: str, code: codons.GeneticCode) -> list[str]:
    """Each row as amino acids (codons translated; '---' -> '-')."""
    if seqtype != "codon":
        return [r.row for r in rows]
    return [codons.translate_row(r.row, code) for r in rows]


# Ranges: "5-30,77" (1-based, inclusive), as the editor's URL carries positions.


def parse_ranges(text: str, limit: int) -> set[int]:
    """0-based positions below `limit` from "5-30,77"; malformed parts are ignored."""
    out: set[int] = set()
    for token in text.split(",")[:MAX_RANGE_TOKENS]:
        match = re.fullmatch(r"\s*(\d{1,7})\s*(?:-\s*(\d{1,7})\s*)?", token)
        if not match:
            continue
        start = int(match.group(1))
        end = int(match.group(2) or start)
        if start > end:
            start, end = end, start
        out.update(range(max(start, 1) - 1, min(end, limit)))
    return out


def format_ranges(positions: Iterable[int]) -> str:
    """0-based positions as "5-30,77" (1-based)."""
    spans: list[list[int]] = []
    for p in sorted(set(positions)):
        if spans and p == spans[-1][1] + 1:
            spans[-1][1] = p
        else:
            spans.append([p, p])
    return ",".join(f"{a + 1}" if a == b else f"{a + 1}-{b + 1}" for a, b in spans)


def natural_key(name: str) -> list[Any]:
    """Sort key: case-insensitive, with numbers in order (seq2 before seq10)."""
    return [(0, int(part)) if part.isdigit() else (1, part.casefold())
            for part in re.split(r"(\d+)", name) if part]  # fmt: skip


SORTS = {"": "As in the file", "name": "Name, A to Z", "-name": "Name, Z to A"}


@dataclass(frozen=True)
class Edits:
    """What the editor shows of a stored alignment: its columns (sites) cut and rows
    removed, as 0-based positions in the stored alignment, and the row order."""

    cut: frozenset[int] = field(default_factory=frozenset)
    drop: frozenset[int] = field(default_factory=frozenset)
    sort: str = ""

    @classmethod
    def from_query(cls, query: Mapping[str, str], rows: int, sites: int) -> Edits:
        cut = parse_ranges(query.get("cut", ""), sites)
        drop = parse_ranges(query.get("drop", ""), rows)
        if len(cut) >= sites:  # never nothing left
            cut = set()
        if len(drop) >= rows:
            drop = set()
        sort = query.get("sort", "")
        return cls(frozenset(cut), frozenset(drop), sort if sort in SORTS else "")

    def query(self) -> dict[str, str]:
        """The URL parameters (only those in effect)."""
        out = {"cut": format_ranges(self.cut), "drop": format_ranges(self.drop), "sort": self.sort}
        return {k: v for k, v in out.items() if v}

    @property
    def key(self) -> str:
        return "|".join(f"{k}={v}" for k, v in sorted(self.query().items()))

    @property
    def changed(self) -> bool:
        return bool(self.cut or self.drop or self.sort)


def row_order(rows: Sequence[Row], edits: Edits) -> list[int]:
    """Stored row positions in the order shown, without the removed ones."""
    shown = [i for i in range(len(rows)) if i not in edits.drop]
    if edits.sort:
        shown.sort(key=lambda i: natural_key(rows[i].name), reverse=edits.sort == "-name")
    return shown


def kept_sites(sites: int, edits: Edits) -> list[int]:
    return [s for s in range(sites) if s not in edits.cut]


def apply(rows: Sequence[Row], edits: Edits, seqtype: str) -> list[Row]:
    """The rows as the editor shows them."""
    k = unit(seqtype)
    keep = kept_sites(site_count(rows, seqtype), edits)
    out = []
    for i in row_order(rows, edits):
        r = rows[i]
        text = r.row if not edits.cut else "".join(r.row[s * k : s * k + k] for s in keep)
        out.append(Row(r.name, r.description, text))
    return out


def gap_sites(rows: Sequence[Row], seqtype: str) -> list[int]:
    """Sites (of these rows) where every row has a gap."""
    k = unit(seqtype)
    return [
        s for s in range(site_count(rows, seqtype))
        if all(r.row[s * k : s * k + k].strip(GAP) == "" for r in rows)
    ]  # fmt: skip


# One-off edits: parameters naming positions in the view as shown (current columns and
# rows, 1-based), turned into stored positions here; the page then shows the canonical
# URL. An edit that would leave no column or no sequence is refused.
ONE_OFF = ("cut_columns", "keep_columns", "drop_rows", "keep_rows", "cut_gaps", "drop_empty",
           "sort_toggle", "sort_reset")  # fmt: skip


def empty_rows(rows: Sequence[Row]) -> list[int]:
    """Rows (of these) with no residue at all, e.g. after cutting down to one domain."""
    return [i for i, r in enumerate(rows) if not r.row.strip(GAP)]


def edit(rows: Sequence[Row], edits: Edits, query: Mapping[str, str], seqtype: str
         ) -> tuple[Edits, str]:  # fmt: skip
    """`edits` after the query's one-off edits, and a note if one was refused."""
    sites = site_count(rows, seqtype)
    view_sites = kept_sites(sites, edits)
    view_rows = row_order(rows, edits)
    cut, drop, sort = set(edits.cut), set(edits.drop), edits.sort
    note = ""

    def stored(positions: set[int], mapping: list[int]) -> set[int]:
        return {mapping[p] for p in positions}

    if text := query.get("cut_columns", ""):
        cut |= stored(parse_ranges(text, len(view_sites)), view_sites)
    if text := query.get("keep_columns", ""):
        keep = stored(parse_ranges(text, len(view_sites)), view_sites)
        if keep:
            cut |= set(view_sites) - keep
    if text := query.get("drop_rows", ""):
        drop |= stored(parse_ranges(text, len(view_rows)), view_rows)
    if text := query.get("keep_rows", ""):
        keep = stored(parse_ranges(text, len(view_rows)), view_rows)
        if keep:
            drop |= set(view_rows) - keep
    if query.get("cut_gaps"):
        shown = apply(rows, Edits(frozenset(cut), frozenset(drop)), seqtype)
        remaining = kept_sites(sites, Edits(frozenset(cut)))
        cut |= {remaining[s] for s in gap_sites(shown, seqtype)}
    if query.get("drop_empty"):
        shown = apply(rows, Edits(frozenset(cut), frozenset(drop), sort), seqtype)
        order = row_order(rows, Edits(drop=frozenset(drop), sort=sort))
        drop |= {order[i] for i in empty_rows(shown)}
    if query.get("sort_toggle"):
        sort = "-name" if sort == "name" else "name"
    if query.get("sort_reset"):
        sort = ""
    if len(cut) >= sites:
        cut, note = set(edits.cut), "That would remove every column: nothing was cut."
    if len(drop) >= len(rows):
        drop, note = set(edits.drop), "That would remove every sequence: nothing was removed."
    return Edits(frozenset(cut), frozenset(drop), sort), note


def provenance(
    rows: Sequence[Row], edits: Edits, seqtype: str, parent_version: int = 1
) -> dict[str, Any]:
    """How an edited alignment came from its parent: the parent's sites kept (1-based
    ranges), the sequences removed, the order, and which version of the parent. Stored
    with the new alignment."""
    sites = site_count(rows, seqtype)
    return {
        "parent_sites": sites,
        "kept_sites": format_ranges(kept_sites(sites, edits)),
        "removed_sequences": [rows[i].name for i in sorted(edits.drop)],
        "sort": edits.sort,
        "parent_version": parent_version,
    }


def compose(previous: Mapping[str, Any], rows: Sequence[Row], edits: Edits, seqtype: str
            ) -> dict[str, Any]:  # fmt: skip
    """An overwritten alignment's record: `edits` (on `rows`, its current content) after
    the edits `previous` records, so its columns still map to its parent's, or to the
    uploaded file's for an upload."""
    sites = site_count(rows, seqtype)
    if previous.get("kept_sites") and previous.get("parent_sites"):
        base = sorted(parse_ranges(previous["kept_sites"], int(previous["parent_sites"])))
        parent_sites = int(previous["parent_sites"])
    else:
        base, parent_sites = list(range(sites)), sites
    if len(base) != sites:  # a record that does not describe these rows: start afresh
        base, parent_sites, previous = list(range(sites)), sites, {}
    kept = [base[s] for s in kept_sites(sites, edits)]
    removed = list(previous.get("removed_sequences") or [])
    removed += [rows[i].name for i in sorted(edits.drop)]
    out = {
        "parent_sites": parent_sites,
        "kept_sites": format_ranges(kept),
        "removed_sequences": removed,
        "sort": edits.sort or previous.get("sort", ""),
    }
    if "parent_version" in previous:
        out["parent_version"] = previous["parent_version"]
    return out


# Storage and records.


def read(record: Alignment) -> list[Row]:
    return read_rows(storage.read_bytes(record.storage_key).decode("utf-8", errors="replace"))


def genetic_code(record: Alignment) -> codons.GeneticCode:
    return codons.CODES.get(record.genetic_code) or codons.STANDARD


def touch(*records: Alignment) -> None:
    """Note that these were used now (opened, downloaded, used in a job)."""
    Alignment.objects.filter(id__in=[r.id for r in records]).update(last_used_at=timezone.now())


def recent(user: AbstractBaseUser | AnonymousUser, search: str = "") -> QuerySet[Alignment]:
    """The user's alignments, the most recently used first (else by when they were made:
    a job's when the job finished), optionally those whose name or job's name contains
    every word of `search`."""
    found = Alignment.objects.for_user(user).select_related("job")
    for word in search.split()[:8]:
        found = found.filter(Q(name__icontains=word) | Q(job__name__icontains=word))
    used = Coalesce("last_used_at", "job__finished_at", "created_at")
    return found.annotate(used=used).order_by("-used", "job_key")


def saved_count(user: AbstractBaseUser | AnonymousUser) -> int:
    return Alignment.objects.for_user(user).filter(job__isnull=True).count()


def store(owner: Any, name: str, rows: Sequence[Row], seqtype: str, code: int, **fields: Any
          ) -> Alignment:  # fmt: skip
    """A new saved alignment (an upload or an edit) with its file. A copy of simulated
    data is simulated too."""
    parent = fields.get("derived_from")
    if parent is not None and parent.simulated:
        fields["simulated"] = True
    record = Alignment(
        owner=owner, name=clean_name(name) or "Alignment", seqtype=seqtype, genetic_code=code,
        sequences=len(rows), columns=site_count(rows, seqtype), **fields,
    )  # fmt: skip
    storage.put_bytes(record.storage_key, to_fasta(rows).encode())
    record.save()
    return record


def overwrite(record: Alignment, rows: Sequence[Row], edits: Edits) -> None:
    """Replace a saved alignment's rows with `rows` as edited (never a job's: the job
    depends on them). Jobs that used it keep their own copy of their inputs."""
    if record.from_job:
        raise ValueError("a job's alignment can't be overwritten")
    shown = apply(rows, edits, record.seqtype)
    storage.put_bytes(record.storage_key, to_fasta(shown).encode())
    record.edits = compose(record.edits, rows, edits, record.seqtype)
    record.sequences = len(shown)
    record.columns = site_count(shown, record.seqtype)
    record.version += 1
    record.last_used_at = timezone.now()
    record.save()


def delete(record: Alignment) -> None:
    """Delete a saved alignment: its file goes, the record stays. A job's alignments go
    only with the job (phylobot.jobs.delete)."""
    if record.from_job:
        raise ValueError("a job's alignment is deleted with its job")
    record.deleted_at = timezone.now()
    record.save(update_fields=["deleted_at"])
    storage.delete_prefix(record.storage_prefix)


def register_job(job: Job) -> list[Alignment]:
    """Record a finished job's alignments in its owner's workspace: one per alignment
    the job published (an aligner's or the user's own), the codon alignment where there
    is one. Idempotent. A branch-site test publishes none of its own."""
    from phylobot.alignment_views import job_keys, sources  # it imports more than we need
    from phylobot.simulations import is_simulation, register

    if is_simulation(job):  # each replicate's true alignment, flagged simulated
        return register(job)
    if job.status != Job.Status.SUCCEEDED or job.params.get("branch_site"):
        return []
    keys, summary = job_keys(job, inherited=False)  # an extension's parent registers its own
    found = []
    code = int(job.params.get("genetic_code") or 1)
    for source in sources(job, keys, summary):
        existing = Alignment.objects.filter(job=job, job_key=source.key).first()
        if existing is not None:
            found.append(existing)
            continue
        path = source.codon or source.amino
        seqtype = "codon" if source.codon else "protein"
        text = storage.read_bytes(f"{job.storage_prefix}results/{path}").decode(errors="replace")
        try:
            rows = read_rows(text)
        except validate.InputError:
            continue
        try:
            with transaction.atomic():
                record = Alignment.objects.create(
                    owner_id=job.owner_id, name=clean_name(source.label) or source.key,
                    origin=Alignment.Origin.JOB, seqtype=seqtype, genetic_code=code,
                    sequences=len(rows), columns=site_count(rows, seqtype), job=job,
                    job_key=source.key, path=path,
                )  # fmt: skip
        except IntegrityError:  # registered meanwhile (the worker and a page at once)
            record = Alignment.objects.get(job=job, job_key=source.key)
        found.append(record)
    return found


def ensure_registered(user: AbstractBaseUser | AnonymousUser) -> None:
    """Register the alignments of the user's finished jobs that have none yet (jobs from
    before this registry, or whose worker had not registered them)."""
    # (Branch-site tests stay pending, as register_job skips them before reading anything.)
    pending = Job.objects.for_user(user).filter(
        status=Job.Status.SUCCEEDED, alignments__isnull=True
    )
    for job in pending:
        try:
            register_job(job)
        except Exception:  # a job with unreadable results is left out, not an error page
            continue
