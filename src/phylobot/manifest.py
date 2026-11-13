"""Each job's result files (path under results/, size) in Postgres, so pages never list
object storage: listing the 285-combination GFP job's 524,069 objects took 32 s, past
gunicorn's timeout.

Written where results are stored: by the worker after each attempt's upload (from the
local folder it uploaded), by `job_import`, and, for jobs from before the manifest or
whose upload broke off, by worker maintenance listing storage once (`index_unindexed`).
`Job.results_indexed_at` says the rows are whole; it is cleared before an upload starts,
so an upload that fails partway leaves the job to maintenance. Rows are only added or
resized, never removed while the files exist (an attempt's upload never deletes either);
purging the files removes them.

Only results/ is indexed: inputs/ is written once by the web tier at submission and is
small, so its readers still list it.
"""

from pathlib import Path
from uuid import UUID

from django.db import transaction
from django.utils import timezone

from phylobot import storage
from phylobot.models import Job, ResultFile

BATCH = 5000
MAINTENANCE_LIMIT = 20  # jobs indexed per maintenance pass (GFP alone lists for ~32 s)


def results_prefix(job: Job) -> str:
    return f"{job.storage_prefix}results/"


def add(job_id: UUID | str, sizes: dict[str, int]) -> None:
    """Record files (relative to results/) as stored: new paths added, sizes updated."""
    rows = [ResultFile(job_id=job_id, path=path, size=size) for path, size in sizes.items()]
    ResultFile.objects.bulk_create(
        rows, batch_size=BATCH, update_conflicts=True, unique_fields=["job", "path"],
        update_fields=["size"],
    )  # fmt: skip


def directory_sizes(directory: Path) -> dict[str, int]:
    """{path relative to `directory`: bytes} for every file under it, as upload_tree
    names them."""
    return {
        path.relative_to(directory).as_posix(): path.stat().st_size
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }


def unindex(job_id: UUID | str) -> None:
    """Mark a job's rows incomplete (before an upload that may break off)."""
    Job.objects.filter(id=job_id).update(results_indexed_at=None)


def indexed(job_id: UUID | str, sizes: dict[str, int]) -> None:
    """Record a whole upload's files and mark the job's rows complete."""
    with transaction.atomic():
        add(job_id, sizes)
        Job.objects.filter(id=job_id).update(results_indexed_at=timezone.now())


def index_from_storage(job: Job) -> int:
    """List the job's results/ once and record what is there (older jobs, broken
    uploads). Returns how many files."""
    prefix = results_prefix(job)
    sizes = {key[len(prefix) :]: size for key, size in storage.list_sizes(prefix).items()}
    indexed(job.id, sizes)
    return len(sizes)


def index_unindexed(limit: int = MAINTENANCE_LIMIT) -> int:
    """Index finished jobs that have no complete manifest (worker maintenance)."""
    from phylobot.jobs import TERMINAL

    pending = list(
        Job.objects.filter(
            results_indexed_at__isnull=True,
            status__in=TERMINAL,
            deleted_at__isnull=True,
            purged_at__isnull=True,
        ).order_by("-created_at")[:limit]
    )
    for job in pending:
        index_from_storage(job)
    return len(pending)


def forget(job_id: UUID | str) -> None:
    """The files are gone (purged): so are their rows."""
    ResultFile.objects.filter(job_id=job_id).delete()


def sizes(job: Job, under: str = "") -> dict[str, int]:
    """{path relative to results/: bytes} of the job's own result files, optionally only
    those under a folder (`under`, ending in "/")."""
    found = ResultFile.objects.filter(job_id=job.id)
    if under:
        found = found.filter(path__startswith=under)
    return dict(found.values_list("path", "size"))


def paths(job: Job, under: str = "") -> set[str]:
    found = ResultFile.objects.filter(job_id=job.id)
    if under:
        found = found.filter(path__startswith=under)
    return set(found.values_list("path", flat=True))


def has(job: Job, path: str) -> bool:
    return ResultFile.objects.filter(job_id=job.id, path=path).exists()
