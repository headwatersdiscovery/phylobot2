"""The home page of a running PhyloBot: where to start, and where you were.

Whoever sees this page has PhyloBot running already, so it says nothing about installing
it (that is phylobot.com, built from ``website/``). It offers the ways in, the reader's
latest jobs and searches with their state (refreshed while any is unfinished), and
"Jump Back In": the result pages they last had open.

That last list is kept by the browser, not the server (``recent.js``, in localStorage):
every viewer is a GET form whose address is its whole state, so the address is all there is
to remember, and a page view costs no write here. Each page worth returning to names its
kind on ``<body>`` (``skins.PAGES``); the home page asks `alive` which of the
remembered records still exist, so a deleted job's pages drop off the list.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Any

from django.http import HttpRequest, HttpResponse, JsonResponse
from django.shortcuts import render
from django.urls import reverse
from django.utils import timezone
from django.views.decorators.http import require_GET

from phylobot import eta, job_types
from phylobot.models import Alignment, HomologSearch, Job, Library

RECENT_WORK = 6  # jobs and searches listed on the home page
MAX_IDS = 60  # records `alive` looks up in one request (the browser remembers far fewer)

UNFINISHED = (Job.Status.QUEUED, Job.Status.RUNNING)


def ago(when: datetime, now: datetime | None = None) -> str:
    """How long ago, in a few characters: the list is read at a glance."""
    now = now or timezone.now()
    elapsed = max(now - when, timedelta(0))
    minutes = int(elapsed.total_seconds() // 60)
    if minutes < 1:
        return "just now"
    if minutes < 60:
        return f"{minutes} min ago"
    if minutes < 24 * 60:
        return f"{minutes // 60} h ago"
    days = elapsed.days
    if days == 1:
        return "yesterday"
    if days < 14:
        return f"{days} days ago"
    return timezone.localtime(when).strftime("%Y-%m-%d")


@dataclass
class Work:
    """One row of the home page's list: a job or a homolog search."""

    name: str
    url: str
    status: str
    note: str  # "ETA: 5 minutes", "finished 2 h ago"
    active: datetime  # when it last changed state, for the order
    types: list[Any] = field(default_factory=list)  # a job's type pills
    search: bool = False
    simulated: bool = False

    @property
    def unfinished(self) -> bool:
        return self.status in UNFINISHED


def _note(status: str, started: datetime, finished: datetime | None, now: datetime) -> str:
    if status == Job.Status.QUEUED:
        return f"submitted {ago(started, now)}"
    if finished is None:
        return ""
    verb = {"succeeded": "finished", "cancelled": "stopped"}.get(status, status)
    return f"{verb} {ago(finished, now)}"


def recent_work(user: Any, limit: int = RECENT_WORK) -> list[Work]:
    """The user's latest jobs and homolog searches: unfinished ones first, then by when
    each last changed state."""
    now = timezone.now()
    jobs = list(Job.objects.for_user(user).order_by("-created_at")[:limit])
    searches = list(HomologSearch.objects.for_user(user).order_by("-created_at")[:limit])
    typical = eta.typical_durations() if any(j.status == Job.Status.RUNNING for j in jobs) else None
    rows: list[Work] = []
    for job in jobs:
        note = _note(job.status, job.created_at, job.finished_at, now)
        if job.status == Job.Status.RUNNING:
            note = "ETA: " + eta.describe(eta.remaining_seconds(job, typical), short=True)
        rows.append(
            Work(
                name=job.name or str(job.id)[:8],
                url=reverse("job_detail", args=[job.id]),
                status=job.status,
                note=note,
                active=job.finished_at or job.started_at or job.created_at,
                types=job_types.types_of(job.params),
                simulated=bool(job.params.get("simulated")),
            )
        )
    for search in searches:
        rows.append(
            Work(
                name=search.name or str(search.id)[:8],
                url=reverse("homolog_search", args=[search.id]),
                status=search.status,
                note=_note(search.status, search.created_at, search.finished_at, now),
                active=search.finished_at or search.created_at,
                search=True,
            )
        )
    rows.sort(key=lambda row: (row.unfinished, row.active), reverse=True)
    return rows[:limit]


def _work_context(request: HttpRequest) -> dict[str, Any]:
    rows = recent_work(request.user)
    return {"work": rows, "work_unfinished": any(row.unfinished for row in rows)}


def home(request: HttpRequest) -> HttpResponse:
    """The home page: ways in, the latest work, and the pages last open."""
    # The PhyloBot '95 skin's hit counter: the reader's own jobs, deleted ones included.
    jobs = Job.objects
    job_count = jobs.for_user(request.user).count() + jobs.deleted_for(request.user).count()
    context = {**_work_context(request), "job_count": job_count}
    return render(request, "phylobot/home.html", context)


@require_GET
def home_work(request: HttpRequest) -> HttpResponse:
    """The home page's list of jobs and searches alone: it asks again while any of them is
    unfinished."""
    return render(request, "phylobot/_home_work.html", _work_context(request))


def _uuids(values: list[str]) -> list[str]:
    from uuid import UUID

    found: list[str] = []
    for value in values[:MAX_IDS]:
        try:
            found.append(str(UUID(value)))
        except ValueError:
            continue
    return found


@require_GET
def alive(request: HttpRequest) -> JsonResponse:
    """Which of these records (jobs, homolog searches, libraries, alignments: the IDs in the
    addresses "Jump Back In" remembers) are still the reader's and not deleted."""
    ids = _uuids(request.GET.getlist("id"))
    live: set[str] = set()
    if ids:
        for model in (Job, HomologSearch, Library, Alignment):
            found = model.objects.for_user(request.user).filter(id__in=ids)
            live.update(str(pk) for pk in found.values_list("id", flat=True))
    return JsonResponse({"live": sorted(live)})
