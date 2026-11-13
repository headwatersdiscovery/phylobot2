"""My Jobs lists a user's work: pipeline Jobs and homolog searches together.

A `HomologSearch` is a different model from a `Job` — run by `run_searches`, not the
pipeline, with no local compute — so rather than merge the models this adapts both to a
common `WorkRow` and interleaves them by the table's chosen order. This first cut keeps the
rich job filters job-only: a search is shown unless a job-only filter (coding/protein,
origin, tag, or a Type filter that leaves out homolog searches) rules it out. Each row
links to its own detail page; search rows carry no compute figures and are not bulk-
selectable here (delete them from Find homologs).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any

from django.contrib.auth.models import AbstractBaseUser, AnonymousUser
from django.db.models import Q
from django.urls import reverse

from phylobot import job_types
from phylobot.models import HomologSearch, Job

HOMOLOG = "homolog"  # the pseudo-type the Type filter offers for homolog searches
SEARCH_LABEL = "Homolog search"
SEARCH_CAP = 500  # most recent searches considered for the page

_TYPE_RANK = {"simulation": 6, "esr": 5, "dnds": 4, "ancestors": 3, "trees": 2}  # as the pills read
_SEARCH_RANK = 7  # searches sort after simulations, as their own kind
_HEX = "0123456789abcdef-"


@dataclass
class WorkRow:
    """One row of the My Jobs table: a Job or a HomologSearch, by its common fields."""

    kind: str  # "job" | "search"
    obj: Any
    id: str
    name: str
    status: str
    created_at: datetime
    finished_at: datetime | None
    url: str  # its detail page
    cpu_seconds: float | None = None
    storage_bytes: int | None = None
    # Job rows carry these (the template reads them); search rows leave them empty.
    types: list[Any] = field(default_factory=list)
    eta: str = ""
    made: Any = None
    readers: list[str] = field(default_factory=list)
    # Grouped view (`nest`): the row a follow-up job is listed under, and on that row how
    # many follow-ups it holds and what to say of them.
    under: str = ""
    follow_ups: int = 0
    follow_up_note: str = ""

    @property
    def is_search(self) -> bool:
        return self.kind == "search"

    @property
    def short_id(self) -> str:
        return self.id[:8]


def job_row(job: Job, types: list[Any], eta: str, made: Any, readers: list[str]) -> WorkRow:
    return WorkRow(
        "job", job, str(job.id), job.name or "", job.status, job.created_at, job.finished_at,
        reverse("job_detail", args=[job.id]), job.cpu_seconds, job.storage_bytes,
        types, eta, made, readers,
    )  # fmt: skip


def search_row(search: HomologSearch) -> WorkRow:
    return WorkRow(
        "search", search, str(search.id), search.name or "", search.status, search.created_at,
        search.finished_at, reverse("homolog_search", args=[search.id]),
    )  # fmt: skip


# Where a follow-up job's params name the job it was run on (the tree viewer's
# branch-site test, held-out reconstruction, "Measure properties", Uncertainty).
_PARENT_KEYS = (
    ("uncertainty_of", "parent"),
    ("properties_of", "parent"),
    ("branch_site_of", "job"),
    ("esr_of", "job"),
)
NEST_COOKIE = "phylobot_jobs_grouped"  # the My Jobs toggle, remembered per browser


def parent_id(job: Job) -> str:
    """The job a follow-up job was run on, or that an extension extends ("" for a job of
    its own). A duplicate (`rerun_of`) is a job of its own."""
    params = job.params or {}
    for key, name in _PARENT_KEYS:
        of = params.get(key)
        if isinstance(of, dict) and of.get(name):
            return str(of[name])
    return str(job.extends_id) if job.extends_id else ""


def nest(rows: list[WorkRow]) -> list[WorkRow]:
    """The rows with each follow-up job moved under the job it was run on, where that job
    is in the list too (otherwise it stays where it was). A follow-up of a follow-up goes
    under the first job of the chain. Groups keep the table's order, and so do the
    follow-ups within one."""
    listed = {row.id: row for row in rows if not row.is_search}
    parents = {row.id: parent_id(row.obj) for row in listed.values()}

    def root(row_id: str) -> str:
        seen = {row_id}
        while (up := parents.get(row_id, "")) and up in listed and up not in seen:
            seen.add(up)
            row_id = up
        return row_id

    under: dict[str, list[WorkRow]] = {}
    top: list[WorkRow] = []
    for row in rows:
        first = row.id if row.is_search else root(row.id)
        if first == row.id:
            top.append(row)
        else:
            row.under = first
            under.setdefault(first, []).append(row)
    out: list[WorkRow] = []
    for row in top:
        children = under.get(row.id, [])
        if children:
            row.follow_ups = len(children)
            unfinished = [
                f"{sum(c.status == status for c in children)} {status}"
                for status in job_types.STATUS_ORDER
                if status != Job.Status.SUCCEEDED and any(c.status == status for c in children)
            ]
            row.follow_up_note = " · ".join(unfinished)
        out.append(row)
        out.extend(children)
    return out


def _rank(row: WorkRow) -> int:
    if row.is_search:
        return _SEARCH_RANK
    return max((_TYPE_RANK.get(t.key, 1) for t in row.types), default=1)


def _value(row: WorkRow, key: str) -> Any:
    """The row's value for the active sort column; None sorts last whatever the direction."""
    match key:
        case "finished":
            return row.finished_at
        case "name":
            return row.name.lower() or None
        case "status":
            order = job_types.STATUS_ORDER
            return order.index(row.status) if row.status in order else len(order)
        case "cpu":
            return row.cpu_seconds
        case "disk":
            return row.storage_bytes
        case "type":
            return _rank(row)
        case "tags":
            if row.is_search:
                return None
            names = sorted(t.name for t in row.obj.tags.all())
            return names[0] if names else None
    return row.created_at  # "submitted", and the default


def order(rows: list[WorkRow], sort: job_types.Sort) -> list[WorkRow]:
    """The merged rows in the table's order (present values first, nulls always last)."""
    present = [r for r in rows if _value(r, sort.key) is not None]
    absent = [r for r in rows if _value(r, sort.key) is None]
    present.sort(key=lambda r: _value(r, sort.key), reverse=sort.descending)
    return present + absent


def want_searches(filters: job_types.Filters, requested_types: list[str]) -> bool:
    """Whether searches belong in the list: no job-only filter excludes them, and the Type
    filter does not restrict to job types only (homolog selected, or no real type filter —
    an unknown value, dropped for jobs, must not hide searches either)."""
    if filters.data or filters.origin or filters.tag:
        return False
    has_job_type = any(t in job_types.BY_KEY for t in requested_types)
    return HOMOLOG in requested_types or not has_job_type


def searches_for(
    user: AbstractBaseUser | AnonymousUser, filters: job_types.Filters
) -> list[WorkRow]:
    """The user's homolog searches that pass the filters searches can honour (name/id, status,
    date windows), newest first, as WorkRows."""
    rows = HomologSearch.objects.for_user(user)
    if filters.q:
        match = Q(name__icontains=filters.q)
        if filters.q and all(c in _HEX for c in filters.q.lower()):
            match |= Q(id__startswith=filters.q.lower())
        rows = rows.filter(match)
    if filters.statuses:
        rows = rows.filter(status__in=filters.statuses)
    for when, start, end in (
        ("created_at", filters.submitted_from, filters.submitted_to),
        ("finished_at", filters.finished_from, filters.finished_to),
    ):
        if start:
            rows = rows.filter(**{f"{when}__date__gte": start})
        if end:
            rows = rows.filter(**{f"{when}__date__lte": end})
    return [search_row(s) for s in rows[:SEARCH_CAP]]
