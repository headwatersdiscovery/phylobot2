"""Mutational paths the user pinned, to return to without tracing them again.

A saved path is the path page's own address under a name: the combination, the two ends
and the view (which alignments are shown beside it, colours, layout). Opening one goes
through the page like any other address, so everything in it is checked against the job
again; nothing stored is trusted. They belong to the job and are listed on its tree
viewer and on the path page.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any
from urllib.parse import urlencode

from django.urls import reverse

from phylobot.models import Job, SavedPath

# The tree viewer's own display settings: what the path page carries (as "back") so its
# way back leads to the tree as it was drawn. The path and the combination are the path
# page's to say, so they are not among them.
TREE_VIEW_KEYS = (
    "lengths", "node_labels", "spacing", "species", "colour", "width", "labels", "marks",
    "branch_palette", "dnds",
)  # fmt: skip
NAME_LIMIT = 80
MAX_SAVED = 200  # per job: a list to choose from, not an archive


def clean_name(text: str) -> str:
    """A name as stored: one line, single spaces, at most NAME_LIMIT characters."""
    return " ".join(str(text).split())[:NAME_LIMIT]


@dataclass(frozen=True)
class Entry:
    id: int
    name: str
    combination: str
    start: str
    end: str
    url: str


def url_of(job: Job, view: dict[str, Any]) -> str:
    return f"{reverse('path_viewer', args=[job.id])}?{urlencode(view, doseq=True)}"


def listing(job: Job) -> list[Entry]:
    return [
        Entry(p.id, p.name, p.combination, p.start, p.end, url_of(job, p.view))
        for p in SavedPath.objects.filter(job=job)
    ]


def find(job: Job, combination: str, start: str, end: str) -> SavedPath | None:
    return SavedPath.objects.filter(job=job, combination=combination, start=start, end=end).first()


def save(job: Job, name: str, combination: str, start: str, end: str, view: dict[str, Any]) -> str:
    """Pin a path (or rename and re-save one already pinned); "" or why not. The caller
    has checked the combination and both ends against the job."""
    name = clean_name(name)
    if not name:
        return "Give the path a name."
    existing = find(job, combination, start, end)
    if existing is None and SavedPath.objects.filter(job=job).count() >= MAX_SAVED:
        return f"This job already has {MAX_SAVED} saved paths: remove one first."
    SavedPath.objects.update_or_create(
        job=job, combination=combination, start=start, end=end,
        defaults={"name": name, "view": view},
    )  # fmt: skip
    return ""


def remove(job: Job, saved_id: int) -> None:
    SavedPath.objects.filter(job=job, id=saved_id).delete()
