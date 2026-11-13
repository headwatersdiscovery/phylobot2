"""Import a job from a PhyloBot ZIP (phylobot.bundle: the job page's "Download ZIP").

The ZIP is hostile input: it may come from anyone, and every byte of it can be edited.
So nothing in it is extracted to disk by name, run, or trusted beyond its type:

- Only the manifest's files are accepted, each under `inputs/`, `results/` or a library
  export's folder, by a path checked part by part (no `..`, absolute paths, backslashes
  or control characters). They go to object storage under a new job id, never to a
  path built from the ZIP.
- Sizes are capped (the upload, the files' total, their number, the manifest, libraries
  and exports; the Settings page sets them, `site_settings.import_limits`), and every
  file is counted as it is read, so a forged header cannot hide a bomb. Its SHA-256 must
  match the manifest's.
- The record's fields are checked by type and bounds, and only finished jobs are taken.
  The job gets a new id; links to jobs this server does not have are dropped, and kept
  aside so importing the linked job later restores them.
- An imported job is never resumed or run: its results were computed elsewhere, so it
  carries `imported_at` (a badge wherever it appears), and the CPU it used there is
  shown as provenance but never counted as ours (`cpu_seconds` stays null). Its files
  count as storage here from the moment of import.

Import is offered in local mode only (owner decision 2026-10-19): cloud mode has no
storage quotas yet.
"""

from __future__ import annotations

import hashlib
import json
import re
import tempfile
import uuid
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any

from django.db import transaction
from django.utils import timezone

from phylobot import analogs, ancestor_names, bundle, clade_collapse, jobs, site_settings, storage
from phylobot.manifest import indexed
from phylobot.models import (
    TAG_PATTERN,
    AncestorName,
    AncestorNote,
    CladeCollapse,
    Job,
    Library,
    LibraryExport,
)

# The size and count limits are on the Settings page (site_settings.import_limits).
# These two stay fixed: a path is part of a storage key (S3 allows 1024 bytes in all),
# and a long error message is cut, never refused.
MAX_PATH = 400
MAX_ERROR = 10_000
# Params naming another job (by `job`).
LINKS = ("simulated_from", "branch_site_of", "extends", "esr_of")
UNLINKED = "unlinked_job"  # the other job's id at the source, until it is imported
STAGED = "imports/"  # uploads waiting for a name (views.import_job), per user
SHA256 = re.compile(r"[0-9a-f]{64}")
PART = re.compile(r"[^\x00-\x1f\x7f/\\]{1,255}")


class Refused(ValueError):
    """Why a ZIP cannot be imported, in words for its user."""


@dataclass
class Archive:
    """A checked ZIP: its manifest and its files (path in the job -> ZIP member)."""

    path: Path
    manifest: dict[str, Any]
    members: dict[str, zipfile.ZipInfo]

    @property
    def job(self) -> dict[str, Any]:
        record: dict[str, Any] = self.manifest["job"]
        return record

    @property
    def name(self) -> str:
        return str(self.job["name"])


def clean_name(name: str) -> str:
    """As jobs.rename: whitespace collapsed, at most Job.name's length."""
    return " ".join(str(name).split())[: jobs.MAX_NAME]


def clashes(owner: Any, name: str) -> bool:
    """Whether one of the owner's jobs already has this name (blank names never clash)."""
    return bool(name) and Job.objects.for_user(owner).filter(name__iexact=name).exists()


def suggest_name(owner: Any, name: str) -> str:
    """'Ime2' -> 'Ime2 (copy)', or 'Ime2 (copy 2)'... if that is taken too."""
    for n in range(1, 1000):
        suffix = " (copy)" if n == 1 else f" (copy {n})"
        candidate = name[: jobs.MAX_NAME - len(suffix)] + suffix
        if not clashes(owner, candidate):
            return candidate
    return name


# --- reading -------------------------------------------------------------------------


def _path(text: str) -> str:
    """A file's path inside the job, checked part by part."""
    if len(text) > MAX_PATH or text.startswith("/"):
        raise Refused(f"The ZIP holds a file with an unacceptable name: {text[:80]!r}.")
    parts = text.split("/")
    if any(p in {".", ".."} or not PART.fullmatch(p) for p in parts):
        raise Refused(f"The ZIP holds a file with an unacceptable name: {text[:80]!r}.")
    return text


def _when(value: Any, field: str) -> datetime | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise Refused(f"The ZIP's record is unreadable ({field}).")
    try:
        moment = datetime.fromisoformat(value)
    except ValueError:
        raise Refused(f"The ZIP's record is unreadable ({field}).") from None
    if timezone.is_naive(moment):
        raise Refused(f"The ZIP's record is unreadable ({field}).")
    return moment


def _number(value: Any, field: str, integer: bool = False) -> Any:
    kinds = (int,) if integer else (int, float)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, kinds) or value < 0 or value > 1e15:
        raise Refused(f"The ZIP's record is unreadable ({field}).")
    return value


def _json(value: Any, field: str) -> dict[str, Any]:
    if not isinstance(value, dict):  # its size is bounded by the record's
        raise Refused(f"The ZIP's record is unreadable ({field}).")
    return value


def _text(value: Any, field: str, limit: int) -> str:
    if not isinstance(value, str):
        raise Refused(f"The ZIP's record is unreadable ({field}).")
    return value[:limit]


def _check_job(record: Any) -> None:
    if not isinstance(record, dict):
        raise Refused("The ZIP's record of the job is missing.")
    _text(record.get("name"), "name", 10_000)
    if record.get("status") not in {s.value for s in jobs.TERMINAL}:
        raise Refused(
            "The job had not finished when it was downloaded: only finished jobs "
            "(succeeded, failed or cancelled) can be imported."
        )
    _json(record.get("params"), "settings")
    _json(record.get("progress") or {}, "progress")
    _text(record.get("error") or "", "error", MAX_ERROR)
    for field in ("created_at", "started_at", "finished_at"):
        _when(record.get(field), field.replace("_", " "))
    for field in ("work_units", "cpu_seconds", "billable_cpu_seconds", "storage_bytes"):
        _number(record.get(field), field.replace("_", " "))
    if _number(record.get("attempt", 1), "attempt", integer=True) is None:
        raise Refused("The ZIP's record is unreadable (attempt).")
    tags = record.get("tags") or []
    if not isinstance(tags, list) or len(tags) > 200:
        raise Refused("The ZIP's record is unreadable (tags).")
    if record.get("imported_from") is not None:
        _json(record["imported_from"], "provenance")
    _check_ancestor_names(record.get("ancestor_names") or [])
    _check_ancestor_notes(record.get("ancestor_notes") or [])
    _check_clade_collapses(record.get("clade_collapses") or [])


CLADE_ID = re.compile(rf"[0-9a-f]{{{analogs.ID_LENGTH}}}")


def _check_ancestor_notes(entries: Any) -> None:
    """Each description on a well-formed clade, text of a bounded size, no clade twice
    (absent from ZIPs written before ancestors had descriptions)."""
    bad = "The ZIP's record is unreadable (ancestor descriptions)."
    if not isinstance(entries, list) or len(entries) > 10_000:
        raise Refused(bad)
    clades: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise Refused(bad)
        clade, parts, text = entry.get("clade"), entry.get("parts"), entry.get("text")
        if not isinstance(clade, str) or not isinstance(parts, str) or not isinstance(text, str):
            raise Refused(bad)
        if not CLADE_ID.fullmatch(clade) or not CLADE_ID.fullmatch(parts) or clade in clades:
            raise Refused(bad)
        if len(text) > AncestorNote.LIMIT:
            raise Refused(bad)
        clades.add(clade)


def _check_ancestor_names(entries: Any) -> None:
    """Each name as the ancestor viewer would have accepted it, on a well-formed clade,
    with no clade or name twice (the database would refuse both later)."""
    bad = "The ZIP's record is unreadable (ancestor names)."
    if not isinstance(entries, list) or len(entries) > 10_000:
        raise Refused(bad)
    clades: set[str] = set()
    names: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise Refused(bad)
        clade, parts, name = entry.get("clade"), entry.get("parts"), entry.get("name")
        if not isinstance(clade, str) or not isinstance(parts, str) or not isinstance(name, str):
            raise Refused(bad)
        if not CLADE_ID.fullmatch(clade) or not CLADE_ID.fullmatch(parts):
            raise Refused(bad)
        if ancestor_names.problem(name, set()) or clade in clades or name.casefold() in names:
            raise Refused(bad)
        clades.add(clade)
        names.add(name.casefold())


def _check_clade_collapses(entries: Any) -> None:
    """Each folded clade on a well-formed clade ID, with no clade twice; the colours and
    label are cleaned on import, so only the structure is checked here (phylobot.clade_collapse)."""
    bad = "The ZIP's record is unreadable (collapsed clades)."
    if not isinstance(entries, list) or len(entries) > 10_000:
        raise Refused(bad)
    clades: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise Refused(bad)
        clade, parts = entry.get("clade"), entry.get("parts")
        if not isinstance(clade, str) or not isinstance(parts, str):
            raise Refused(bad)
        if not CLADE_ID.fullmatch(clade) or not CLADE_ID.fullmatch(parts) or clade in clades:
            raise Refused(bad)
        clades.add(clade)


def _check_libraries(shelf: Any, limits: site_settings.ImportLimits) -> dict[int, set[int]]:
    """Library number (1, 2, ...) -> its export numbers."""
    if not isinstance(shelf, list):
        raise Refused("The ZIP's list of libraries is unreadable.")
    if len(shelf) > limits.libraries:
        raise Refused(f"The job has more than {limits.libraries:,} libraries.")
    numbers = {}
    for n, record in enumerate(shelf, start=1):
        if not isinstance(record, dict):
            raise Refused("The ZIP's list of libraries is unreadable.")
        _text(record.get("name"), "library name", 1000)
        _text(record.get("code"), "library code", 1000)
        if "notes" in record:  # absent from ZIPs written before libraries had notes
            _text(record["notes"], "library notes", 1_000_000)
        _json(record.get("spec"), "library design")
        _number(record.get("spec_version"), "library version", integer=True)
        for field in ("created_at", "updated_at"):
            _when(record.get(field), "library date")
        exports = record.get("exports") or []
        if not isinstance(exports, list):
            raise Refused("The ZIP's list of library exports is unreadable.")
        if len(exports) > limits.exports:
            raise Refused(f"A library in the ZIP has more than {limits.exports:,} exports.")
        seen: set[int] = set()
        for export in exports:
            if not isinstance(export, dict):
                raise Refused("The ZIP's list of library exports is unreadable.")
            number = _number(export.get("number"), "export number", integer=True)
            if not number or number in seen:
                raise Refused("The ZIP's list of library exports is unreadable.")
            seen.add(number)
            _json(export.get("spec"), "export design")
            _number(export.get("spec_version"), "export version", integer=True)
            _number(export.get("member_count"), "export size", integer=True)
            _text(export.get("phylobot_version") or "", "export version", 1000)
            _when(export.get("created_at"), "export date")
        numbers[n] = seen
    return numbers


def _allowed(path: str, exports: dict[int, set[int]]) -> bool:
    parts = path.split("/")
    if parts[0] in bundle.AREAS:
        return len(parts) > 1
    return (
        len(parts) == 5
        and parts[0] == "libraries"
        and parts[2] == "exports"
        and parts[1].isdigit()
        and parts[3].isdigit()
        and int(parts[3]) in exports.get(int(parts[1]), set())
    )


def size(n: int) -> str:
    """1_500_000_000 -> '1.5 GB' (the Settings page's decimal units)."""
    return f"{n / 1e9:g} GB" if n >= 1e9 else f"{n / 1e6:g} MB"


def read(path: Path, limits: site_settings.ImportLimits | None = None) -> Archive:
    """Check a ZIP without unpacking it: everything but the files' checksums, which
    `import_job` verifies as it copies them."""
    limits = limits or site_settings.import_limits()
    if path.stat().st_size > limits.upload_bytes:
        raise Refused(f"The ZIP is larger than {size(limits.upload_bytes)}.")
    try:
        archive = zipfile.ZipFile(path)
    except zipfile.BadZipFile, OSError:
        raise Refused("That file is not a ZIP.") from None
    with archive:
        infos = [i for i in archive.infolist() if not i.is_dir()]
        if len(infos) > limits.files + 1:
            raise Refused(f"The ZIP holds more than {limits.files:,} files.")
        names = [i.filename for i in infos]
        if len(set(names)) != len(names):
            raise Refused("The ZIP holds two files with the same name.")
        found = [n for n in names if PurePosixPath(n).name == bundle.MANIFEST]
        roots = {n[: -len(bundle.MANIFEST)] for n in found if n.count("/") <= 1}
        if len(roots) != 1:
            raise Refused(
                "This ZIP has no PhyloBot record of its job (phylobot-job.json), so it "
                "cannot be imported. ZIPs downloaded before import existed lack one: "
                "download the job again."
            )
        root = roots.pop()
        info = archive.getinfo(root + bundle.MANIFEST)
        if info.file_size > limits.record_bytes:
            raise Refused(f"The ZIP's record is larger than {size(limits.record_bytes)}.")
        try:
            with archive.open(info) as handle:
                manifest = json.loads(handle.read(limits.record_bytes + 1))
        except ValueError, zipfile.BadZipFile, RuntimeError, OSError:
            raise Refused("The ZIP's record is unreadable.") from None
        if not isinstance(manifest, dict) or manifest.get("format") != bundle.FORMAT:
            raise Refused("The ZIP's record is not a PhyloBot job.")
        version = manifest.get("format_version")
        if not isinstance(version, int) or isinstance(version, bool) or version < 1:
            raise Refused("The ZIP's record is unreadable.")
        if version > bundle.FORMAT_VERSION:
            raise Refused("This ZIP was made by a newer PhyloBot: update this one first.")
        if version < bundle.OLDEST_FORMAT_VERSION:
            raise Refused(
                "This ZIP was made by an older PhyloBot whose results this one can no longer "
                "read: run the job again here instead."
            )
        _check_job(manifest.get("job"))
        exports = _check_libraries(manifest.get("libraries") or [], limits)
        files = manifest.get("files")
        if not isinstance(files, dict):
            raise Refused("The ZIP's list of files is unreadable.")
        members: dict[str, zipfile.ZipInfo] = {}
        total = 0
        for item in infos:
            if item is info:
                continue
            if not item.filename.startswith(root):
                raise Refused("The ZIP holds files outside the job's folder.")
            relative = _path(item.filename[len(root) :])
            if not _allowed(relative, exports):
                raise Refused(f"The ZIP holds a file a job cannot have: {relative[:80]!r}.")
            entry = files.get(relative)
            if (
                not isinstance(entry, dict)
                or not isinstance(entry.get("sha256"), str)
                or not SHA256.fullmatch(entry["sha256"])
                or entry.get("size") != item.file_size
            ):
                raise Refused(f"The ZIP's record does not match its file {relative[:80]!r}.")
            if item.flag_bits & 0x1:
                raise Refused("The ZIP is encrypted.")
            total += item.file_size
            if total > limits.unpacked_bytes:
                limit = size(limits.unpacked_bytes)
                raise Refused(f"The ZIP's files add up to more than {limit}.")
            members[relative] = item
        if set(files) != set(members):
            raise Refused("The ZIP is missing files its record lists.")
    return Archive(path=path, manifest=manifest, members=members)


# --- importing -----------------------------------------------------------------------


def _relink(owner: Any, params: dict[str, Any]) -> dict[str, Any]:
    """Point links to other jobs at this server's copies, or keep them aside."""
    params = dict(params)
    for field in LINKS:
        link = params.get(field)
        if not isinstance(link, dict) or "job" not in link:
            continue
        link = dict(link)
        ref = str(link.pop("job"))
        here = (
            Job.objects.for_user(owner)
            .filter(imported_from__id=ref)
            .order_by("-imported_at")
            .first()
        )
        if here is None and jobs.is_job_id(ref):
            here = Job.objects.for_user(owner).filter(id=ref, imported_at__isnull=True).first()
        if here is not None:
            link["job"] = str(here.id)
        else:
            link[UNLINKED] = ref
        params[field] = link
    return params


def _relink_others(job: Job) -> None:
    """Jobs imported before this one that name it get their link back."""
    ref = job.imported_from.get("id")
    for field in LINKS:
        waiting = Job.objects.for_user(job.owner).filter(**{f"params__{field}__{UNLINKED}": ref})
        for other in waiting:
            link = dict(other.params[field])
            link.pop(UNLINKED, None)
            link["job"] = str(job.id)
            other.params = {**other.params, field: link}
            other.save(update_fields=["params"])


def _copy(archive: zipfile.ZipFile, item: zipfile.ZipInfo, digest: str, dest: Path) -> None:
    """Unpack one file to `dest`, counting and hashing what is really there."""
    sha, size = hashlib.sha256(), 0
    try:
        with archive.open(item) as source, dest.open("wb") as out:
            while chunk := source.read(1024 * 1024):
                size += len(chunk)
                if size > item.file_size:
                    raise Refused("The ZIP's files are larger than it says.")
                sha.update(chunk)
                out.write(chunk)
    except zipfile.BadZipFile, EOFError, RuntimeError, OSError:
        raise Refused("The ZIP is damaged.") from None
    if size != item.file_size or sha.hexdigest() != digest:
        raise Refused(f"A file in the ZIP has been changed or damaged: {item.filename[:80]!r}.")


def _provenance(archive: Archive) -> dict[str, Any]:
    record = archive.job
    earlier = record.get("imported_from") or {}
    return {
        "id": str(record.get("id", ""))[:64],
        "phylobot_version": str(archive.manifest.get("phylobot_version", ""))[:40],
        "exported_at": str(archive.manifest.get("exported_at", ""))[:40],
        # CPU used where it ran (a job imported twice keeps the first server's figures).
        "cpu_seconds": record.get("cpu_seconds") or earlier.get("cpu_seconds"),
        "billable_cpu_seconds": (
            record.get("billable_cpu_seconds") or earlier.get("billable_cpu_seconds")
        ),
        "attempts": record.get("attempt", 1),
        **({"earlier": earlier} if earlier else {}),
    }


def import_job(owner: Any, archive: Archive, name: str) -> Job:
    """A new job of `owner`'s from a checked ZIP (`read`), named `name`."""
    manifest, record = archive.manifest, archive.job
    shelf = manifest.get("libraries") or []
    job_id = uuid.uuid4()
    library_ids = [uuid.uuid4() for _ in shelf]
    prefixes = [f"jobs/{job_id}/"] + [f"libraries/{i}/" for i in library_ids]

    def key_of(relative: str) -> str:
        parts = relative.split("/")
        if parts[0] == "libraries":
            n, number, filename = int(parts[1]), parts[3], parts[4]
            return f"libraries/{library_ids[n - 1]}/exports/{number}/{filename}"
        return f"jobs/{job_id}/{relative}"

    try:
        with zipfile.ZipFile(archive.path) as zipped, tempfile.TemporaryDirectory() as tmp:
            scratch = Path(tmp) / "file"
            for relative, item in archive.members.items():
                _copy(zipped, item, manifest["files"][relative]["sha256"], scratch)
                storage.put_file(key_of(relative), scratch)
                scratch.unlink()
        with transaction.atomic():
            job = Job.objects.create(
                id=job_id,
                owner=owner,
                name=clean_name(name),
                status=record["status"],
                params=_relink(owner, record["params"]),
                started_at=_when(record.get("started_at"), "start"),
                finished_at=_when(record.get("finished_at"), "end"),
                progress=record.get("progress") or {},
                error=(record.get("error") or "")[:MAX_ERROR],
                work_units=record.get("work_units"),
                imported_at=timezone.now(),
                imported_from=_provenance(archive),
            )
            if created := _when(record.get("created_at"), "date"):
                Job.objects.filter(id=job.id).update(created_at=created)
            tags = record.get("tags") or []
            job.add_tags([t for t in tags if isinstance(t, str) and TAG_PATTERN.fullmatch(t)])
            AncestorName.objects.bulk_create(
                AncestorName(job=job, clade=e["clade"], parts=e["parts"], name=e["name"])
                for e in record.get("ancestor_names") or []
            )  # checked by _check_ancestor_names
            AncestorNote.objects.bulk_create(
                AncestorNote(job=job, clade=e["clade"], parts=e["parts"], text=e["text"])
                for e in record.get("ancestor_notes") or []
            )  # checked by _check_ancestor_notes
            CladeCollapse.objects.bulk_create(
                CladeCollapse(
                    job=job, clade=e["clade"], parts=e["parts"], **clade_collapse.cleaned_entry(e)
                )
                for e in record.get("clade_collapses") or []
            )  # structure checked by _check_clade_collapses; values cleaned by cleaned_entry
            for library_id, entry in zip(library_ids, shelf, strict=True):
                _import_library(owner, job, library_id, entry)
            _relink_others(job)
            results = {
                relative.removeprefix("results/"): item.file_size
                for relative, item in archive.members.items()
                if relative.startswith("results/")
            }
            indexed(job.id, results)
    except BaseException:
        for prefix in prefixes:
            storage.delete_prefix(prefix)
        raise
    jobs.measure_storage(job)
    job.refresh_from_db()
    return job


def _import_library(owner: Any, job: Job, library_id: uuid.UUID, entry: dict[str, Any]) -> None:
    record = Library.objects.create(
        id=library_id,
        owner=owner,
        job=job,
        name=" ".join(entry["name"].split())[:100] or "Library",
        code=entry["code"][:6],
        notes=(entry.get("notes") or "")[: Library.NOTES_LIMIT],
        spec=entry["spec"],
        spec_version=entry.get("spec_version") or 1,
    )
    dates = {
        field: moment
        for field in ("created_at", "updated_at")
        if (moment := _when(entry.get(field), "library date"))
    }
    if dates:
        Library.objects.filter(id=record.id).update(**dates)
    for export in entry.get("exports") or []:
        receipt = LibraryExport.objects.create(
            library=record,
            number=export["number"],
            spec=export["spec"],
            spec_version=export.get("spec_version") or 1,
            member_count=export.get("member_count") or 0,
            phylobot_version=(export.get("phylobot_version") or "")[:40],
        )
        if moment := _when(export.get("created_at"), "export date"):
            LibraryExport.objects.filter(id=receipt.id).update(created_at=moment)


# --- uploads waiting for a name ------------------------------------------------------


def staged_key(owner: Any, token: str) -> str:
    return f"{STAGED}{owner.pk}/{token}.zip"


def stage(owner: Any, upload: Path) -> str:
    """Keep an upload while its owner picks a name; one per user (the latest)."""
    discard(owner)
    token = uuid.uuid4().hex
    storage.put_file(staged_key(owner, token), upload)
    return token


def fetch(owner: Any, token: str, dest: Path) -> bool:
    """The staged upload, to `dest`; False if there is none."""
    if not re.fullmatch(r"[0-9a-f]{32}", token):
        return False
    key = staged_key(owner, token)
    if key not in storage.list_keys(f"{STAGED}{owner.pk}/"):
        return False
    storage.get_file(key, dest)
    return True


def discard(owner: Any) -> None:
    storage.delete_prefix(f"{STAGED}{owner.pk}/")


def save_upload(upload: Any, dest: Path) -> None:
    """An uploaded file (Django's UploadedFile) to `dest`, refusing it when too large."""
    limit = site_settings.import_limits().upload_bytes
    if upload.size > limit:
        raise Refused(f"The ZIP is larger than {size(limit)}.")
    with dest.open("wb") as out:
        for chunk in upload.chunks():
            out.write(chunk)
