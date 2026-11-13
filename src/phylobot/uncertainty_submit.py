"""Queueing an uncertainty run on a finished job.

The shape of "Measure properties" (`properties_submit`): a new `Job` with its own
accounting and results prefix, whose inputs are the parent's published files, copied
server-side, and whose params say what was asked for (`uncertainty_of`) with a plan key,
so asking twice for the same settings leads to the one run.

Nothing from the request reaches a file except as a choice among the parent job's own
values: the reference must be one of its combinations, a multiplier names a level by its
position among the job's own levels, and a gap version must be one of `VERSIONS`.
`phylobot-uncertainty` checks the settings again in the pipeline.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Any

from django.db import transaction
from django.http import Http404, HttpRequest, HttpResponse, QueryDict
from django.shortcuts import get_object_or_404, redirect
from django.urls import reverse
from django.views.decorators.http import require_POST

from phylobot import ancestors, jobs, joined, results, storage, workload
from phylobot import uncertainty as u
from phylobot import uncertainty_run as run_mod
from phylobot.executors import UNCERTAINTY_DIR, UNCERTAINTY_SETTINGS
from phylobot.gaps import MASKING
from phylobot.models import Job

# What UNCERTAINTY reads of each combination (uncertainty_run.Combination, .prepare).
COMBINATION_FILES = (
    "alignment.fasta",
    "tree.nwk",
    "tree.rooted.nwk",
    "ancestors.fasta",
    "gaps-ml.tsv",
    *(ancestors.VERSIONS[v][1] for v in run_mod.VERSIONS if v != ancestors.ORIGINAL),
)
WHOLE = "ancestors_posteriors.tsv"
INDEX = ancestors.index_name(WHOLE)


@dataclass(frozen=True)
class Level:
    key: str
    label: str


@dataclass(frozen=True)
class Choices:
    """What a run of this job can be asked for: its combinations with ancestors (the
    reference menu), each factor's levels, and the gap versions some combination has."""

    combinations: list[dict[str, Any]]
    reference: str
    levels: dict[str, list[Level]]
    versions: list[str]


def published(job: Job) -> joined.Results:
    found = joined.results(job)
    if "summary.json" not in found.keys or not found.summary:
        raise Http404
    return found


def reader(files: joined.Results) -> run_mod.Read:
    def read(relative: str) -> str | None:
        return files.read(relative) if relative in files.keys else None

    return read


def choices(files: joined.Results) -> Choices:
    combos = [
        c for c in files.summary.get("combinations", [])
        if f"{c['key']}/ancestors.fasta" in files.keys
        and f"{c['key']}/{WHOLE}" in files.keys
    ]  # fmt: skip
    combos = sorted(combos, key=results.combination_order(files.summary.get("primary", "")))
    levels: dict[str, list[Level]] = {f: [] for f in u.FACTORS}
    for c in combos:
        for factor, key, label in (
            ("alignment", c["aligner"], c.get("aligner_label") or c["aligner"]),
            ("tree", c["tree"], c.get("tree_label") or c["tree"]),
            ("model", c["model"], c.get("model_label") or c["model"]),
        ):
            if key not in [lv.key for lv in levels[factor]]:
                levels[factor].append(Level(key, str(label)))
    versions = [
        v for v in run_mod.VERSIONS
        if any(f"{c['key']}/{ancestors.VERSIONS[v][1]}" in files.keys for c in combos)
    ]  # fmt: skip
    levels["gaps"] = [Level(v, version_label(v)) for v in versions]
    return Choices(combos, run_mod.default_reference(combos), levels, versions)


def version_label(version: str) -> str:
    """A gap version in a few words; the unmasked one said plainly, since "None" beside a
    list of gap methods reads as "no gap versions"."""
    if version == ancestors.ORIGINAL:
        return "Unmasked (a residue at every column)"
    return ancestors.short_version_label(version)


def settings_from(post: QueryDict, offered: Choices) -> run_mod.Settings:
    """The run's settings from the form, every value one of the job's own (else 404)."""
    weighting = post.get("weighting", "flat")
    if weighting not in run_mod.WEIGHTINGS:
        raise Http404
    multipliers: dict[str, dict[str, float]] = {}
    for factor, levels in offered.levels.items():
        for i, level in enumerate(levels):
            raw = (post.get(f"weight-{factor}-{i}") or "").strip()
            if not raw:
                continue
            try:
                value = float(raw)
            except ValueError:
                raise Http404 from None
            if not 0 <= value <= 1000:
                raise Http404
            if value != 1:
                multipliers.setdefault(factor, {})[level.key] = value
    versions = post.getlist("versions")
    if any(v not in offered.versions for v in versions):
        raise Http404
    reference = post.get("reference", "") or offered.reference
    if reference not in [c["key"] for c in offered.combinations]:
        raise Http404
    default_versions = [v for v in offered.versions if v in MASKING] or [ancestors.ORIGINAL]
    chosen = tuple(v for v in run_mod.VERSIONS if v in versions)
    return run_mod.Settings(
        weighting,
        multipliers,
        () if list(chosen) == default_versions else chosen,
        "" if reference == offered.reference else reference,
    )


def plan_key(parent: Job, settings: run_mod.Settings) -> str:
    body = json.dumps({"parent": str(parent.id), "settings": settings.to_json()}, sort_keys=True)
    return hashlib.sha256(body.encode()).hexdigest()[:32]


def runs_of(job: Job) -> list[Job]:
    """The uncertainty runs asked from this job, newest first."""
    return list(
        Job.objects.for_user(job.owner)
        .filter(params__uncertainty_of__parent=str(job.id))
        .order_by("-created_at")
    )


def existing_run(owner: Any, key: str) -> Job | None:
    """A run of exactly this plan that succeeded or is still to finish (failed and
    cancelled ones can be asked for again)."""
    return (
        Job.objects.for_user(owner)
        .filter(params__uncertainty_of__key=key)
        .exclude(status__in=[Job.Status.FAILED, Job.Status.CANCELLED])
        .order_by("created_at")
        .first()
    )


def staged_names(keys: set[str], combination: str) -> list[str]:
    """The files of one combination UNCERTAINTY reads: the posteriors table with its byte
    index (it reads one node at a time)."""
    names = [
        f"{combination}/{name}"
        for name in (*COMBINATION_FILES, INDEX)
        if f"{combination}/{name}" in keys
    ]
    return [*names, f"{combination}/{WHOLE}"]


def stage(job: Job, files: joined.Results, offered: Choices, settings: run_mod.Settings) -> None:
    dest = f"{job.storage_prefix}inputs/{UNCERTAINTY_DIR}/"
    # The joined summary: an extension's run covers the combinations it inherits.
    storage.put_bytes(f"{dest}summary.json", files.read_bytes("summary.json"))
    for combo in offered.combinations:
        for name in staged_names(files.keys, combo["key"]):
            if name in files.keys:
                storage.copy_key(files.key(name), f"{dest}{name}")
    body = json.dumps(settings.to_json(), indent=1) + "\n"
    storage.put_bytes(f"{job.storage_prefix}inputs/{UNCERTAINTY_SETTINGS}", body.encode())


@require_POST
def submit(request: HttpRequest, job_id: str) -> HttpResponse:
    """Queue an uncertainty run on one of the user's finished jobs."""
    parent = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if parent.status != Job.Status.SUCCEEDED or parent.deleted_at is not None:
        raise Http404
    files = published(parent)
    offered = choices(files)
    if not offered.combinations:
        raise Http404
    settings = settings_from(request.POST, offered)
    key = plan_key(parent, settings)
    of = {
        "parent": str(parent.id),
        "settings": settings.to_json(),
        "combinations": len(offered.combinations),
        "key": key,
    }
    params: dict[str, Any] = {
        "uncertainty_only": True,
        "uncertainty_of": of,
        "aligners": ["none"],
        "tree_builders": ["none"],
        "seed": parent.params.get("seed", 1),
    }
    name = f"Uncertainty: {parent.name or str(parent.id)[:8]}"
    size = workload.IME2_SIZE
    page = reverse("uncertainty", args=[parent.id])
    with transaction.atomic():
        Job.objects.select_for_update().filter(id=parent.id).first()  # one submission at a time
        already = existing_run(parent.owner, key)
        if already is not None:
            return redirect(f"{page}?run={already.id}")
        job = Job.objects.create(
            owner=parent.owner,
            name=name[:100],
            params=params,
            work_units=workload.estimate(params, size).units,
            input_size=workload.size_record(size),
        )
        job.tags.set(parent.tags.all())
        stage(job, files, offered, settings)
        transaction.on_commit(lambda: jobs.enqueue(job))
    return redirect(f"{page}?run={job.id}")


def context(job: Job) -> dict[str, Any]:
    """The job page's part for an uncertainty run: the way back to the job it measured."""
    of = job.params.get("uncertainty_of") or {}
    if not of:
        return {}
    parent_id = str(of.get("parent") or "")
    parent = Job.objects.for_user(job.owner).filter(id=parent_id).first() if parent_id else None
    url = f"{reverse('uncertainty', args=[parent.id])}?run={job.id}" if parent else ""
    return {"of": of, "parent": parent, "url": url}
