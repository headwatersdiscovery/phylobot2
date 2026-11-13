"""A job's files as one ZIP download (My Jobs and the job page): the inputs as submitted,
every result file and the job's libraries' exports, under one folder named after the job,
with `phylobot-job.json` beside them: the job's record (name, settings, status, times,
tags, libraries) and every file's SHA-256, so another PhyloBot can import it
(phylobot.job_import). Without the manifest it is still an ordinary folder of files.

Streamed: files are read from object storage and compressed chunk by chunk, so the web
tier never holds a whole job in memory or on disk (zipfile writes a data descriptor after
each file when the output is not seekable). The manifest is written last, once every
file's checksum is known.
"""

import hashlib
import io
import json
import re
import zipfile
from collections.abc import Iterator
from datetime import datetime
from importlib import metadata
from typing import Any

from django.utils import timezone

from phylobot import ancestor_names, joined, storage
from phylobot.models import AncestorName, AncestorNote, CladeCollapse, Job, Library

AREAS = ("inputs", "results")
MANIFEST = "phylobot-job.json"
FORMAT = "phylobot-job"
# 2 (2026-10-30): packed properties, posterior byte indexes, the ancestor
# map. Older ZIPs are refused: nothing reads their layout any more.
FORMAT_VERSION = 2
OLDEST_FORMAT_VERSION = 2


def folder_name(job: Job) -> str:
    """'Ime2 fungi' -> 'Ime2_fungi-1a2b3c4d': safe in any file system and unique per job."""
    name = re.sub(r"[^A-Za-z0-9._+-]+", "_", job.name).strip("._")[:60]
    short = str(job.id)[:8]
    return f"{name}-{short}" if name else f"job-{short}"


def phylobot_version() -> str:
    try:
        return metadata.version("phylobot")[:40]
    except metadata.PackageNotFoundError:
        return "unknown"


class _Sink(io.RawIOBase):
    """A write-only, unseekable buffer that hands over what has been written so far."""

    def __init__(self) -> None:
        self.chunks: list[bytes] = []

    def writable(self) -> bool:
        return True

    def write(self, data: bytes) -> int:  # type: ignore[override]
        self.chunks.append(bytes(data))
        return len(data)

    def take(self) -> bytes:
        data = b"".join(self.chunks)
        self.chunks.clear()
        return data


def libraries(job: Job) -> list[Library]:
    """The job's live libraries, oldest first: their place in this list is their number
    in the ZIP (libraries/<n>/), since their ids mean nothing elsewhere."""
    return list(
        job.libraries.filter(deleted_at__isnull=True)
        .order_by("created_at", "id")
        .prefetch_related("exports")
    )


JOINED = "joined:"  # bundle.keys' marker for a table joined from an extension's parent


def keys(job: Job, shelf: list[Library] | None = None) -> list[tuple[str, str]]:
    """(storage key, path in the ZIP) for every file of the job."""
    folder = folder_name(job)
    found = []
    for area in AREAS:
        prefix = f"{job.storage_prefix}{area}/"
        if area == "results":
            # An extension's results include those it reads from its parent
            # (phylobot.joined), so its ZIP stands alone; the tables both share are
            # written joined (JOINED: no single stored file holds them).
            files = joined.results(job)
            for relative in sorted(files.keys):
                shared = files.inherited and relative in (joined.SUMMARY, *joined.JOINED_TABLES)
                key = JOINED + relative if shared else files.key(relative)
                found.append((key, f"{folder}/{area}/{relative}"))
            continue
        for key in storage.list_keys(prefix):
            found.append((key, f"{folder}/{area}/{key[len(prefix) :]}"))
    for n, record in enumerate(libraries(job) if shelf is None else shelf, start=1):
        for export in record.exports.all():
            prefix = export.storage_prefix
            for key in storage.list_keys(prefix):
                found.append(
                    (key, f"{folder}/libraries/{n}/exports/{export.number}/{key[len(prefix) :]}")
                )
    return found


def _when(moment: datetime | None) -> str | None:
    return moment.isoformat() if moment else None


def manifest(job: Job, shelf: list[Library], files: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """What an import needs besides the files. Ids, owners and links to other jobs are
    left out (they mean nothing on another server), except the job's own id, kept as
    provenance and to re-link jobs that name it."""
    return {
        "format": FORMAT,
        "format_version": FORMAT_VERSION,
        "phylobot_version": phylobot_version(),
        "exported_at": timezone.now().isoformat(),
        "job": {
            "id": str(job.id),
            "name": job.name,
            "status": job.status,
            "params": job.params,
            "created_at": _when(job.created_at),
            "started_at": _when(job.started_at),
            "finished_at": _when(job.finished_at),
            "error": job.error,
            "progress": job.progress,
            "attempt": job.attempt,
            "work_units": job.work_units,
            "cpu_seconds": job.cpu_seconds,
            "billable_cpu_seconds": job.billable_cpu_seconds,
            "storage_bytes": job.storage_bytes,
            "tags": sorted(job.tags.values_list("name", flat=True)),
            # Names the user gave ancestors, by clade (phylobot.ancestor_names); an
            # extension's are its parent's, and travel with it since its ZIP stands alone.
            "ancestor_names": [
                {"clade": n.clade, "parts": n.parts, "name": n.name}
                for n in AncestorName.objects.filter(job=ancestor_names.holder(job)).order_by(
                    "name"
                )
            ],
            # Descriptions the user gave ancestors, by clade like the names.
            "ancestor_notes": [
                {"clade": n.clade, "parts": n.parts, "text": n.text}
                for n in AncestorNote.objects.filter(job=ancestor_names.holder(job)).order_by(
                    "clade"
                )
            ],
            # Clades the user folded into triangles in the tree viewer (phylobot.clade_collapse),
            # also by clade; like the names, an extension's are the parent's and travel along.
            "clade_collapses": [
                {
                    "clade": c.clade,
                    "parts": c.parts,
                    "label": c.label,
                    "fill": c.fill,
                    "border_colour": c.border_colour,
                    "border_width": c.border_width,
                }
                for c in CladeCollapse.objects.filter(job=ancestor_names.holder(job)).order_by(
                    "clade"
                )
            ],  # fmt: skip
            # An imported job keeps where it first came from.
            "imported_from": job.imported_from or None,
        },
        "libraries": [
            {
                "name": record.name,
                "code": record.code,
                "notes": record.notes,
                "spec": record.spec,
                "spec_version": record.spec_version,
                "created_at": _when(record.created_at),
                "updated_at": _when(record.updated_at),
                "exports": [
                    {
                        "number": export.number,
                        "created_at": _when(export.created_at),
                        "spec": export.spec,
                        "spec_version": export.spec_version,
                        "member_count": export.member_count,
                        "phylobot_version": export.phylobot_version,
                    }
                    for export in record.exports.all()
                ],
            }
            for record in shelf
        ],
        "files": files,
    }


def results_zip(
    job: Job, under: str, paths: list[str], folder: str, extra: dict[str, str] | None = None
) -> Iterator[bytes]:
    """Some of a job's result files below `under` (`paths`, relative to results/,
    wherever they live) as a streamed ZIP, in `folder` with their paths below `under`;
    `extra` adds small text files beside them. One row's ancestors on the job page."""
    sink = _Sink()
    when = timezone.localtime(job.finished_at or timezone.now()).timetuple()[:6]
    files = joined.results(job)
    with zipfile.ZipFile(sink, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in paths:
            info = zipfile.ZipInfo(f"{folder}/{path.removeprefix(under)}", date_time=when)
            info.compress_type = zipfile.ZIP_DEFLATED
            with archive.open(info, "w") as member:
                for chunk in storage.stream(files.key(path)):
                    member.write(chunk)
                    if data := sink.take():
                        yield data
            yield sink.take()
        for name, text in (extra or {}).items():
            info = zipfile.ZipInfo(f"{folder}/{name}", date_time=when)
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, text)
    yield sink.take()


def zip_stream(job: Job) -> Iterator[bytes]:
    sink = _Sink()
    # Files are dated when the job finished (object storage keeps no useful times).
    when = timezone.localtime(job.finished_at or timezone.now()).timetuple()[:6]
    folder = folder_name(job)
    shelf = libraries(job)
    files: dict[str, dict[str, Any]] = {}
    with zipfile.ZipFile(sink, "w", zipfile.ZIP_DEFLATED) as archive:
        for key, name in keys(job, shelf):
            info = zipfile.ZipInfo(name, date_time=when)
            info.compress_type = zipfile.ZIP_DEFLATED
            digest, size = hashlib.sha256(), 0
            with archive.open(info, "w") as member:
                chunks = (
                    [joined.results(job).read_bytes(key[len(JOINED) :])]
                    if key.startswith(JOINED)
                    else storage.stream(key)
                )
                for chunk in chunks:
                    member.write(chunk)
                    digest.update(chunk)
                    size += len(chunk)
                    if data := sink.take():
                        yield data
            files[name[len(folder) + 1 :]] = {"sha256": digest.hexdigest(), "size": size}
            yield sink.take()
        info = zipfile.ZipInfo(f"{folder}/{MANIFEST}", date_time=when)
        info.compress_type = zipfile.ZIP_DEFLATED
        text = json.dumps(manifest(job, shelf, files), indent=1, sort_keys=True)
        archive.writestr(info, text + "\n")
    yield sink.take()
