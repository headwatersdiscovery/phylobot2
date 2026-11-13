"""Measuring properties on ancestors a finished job already published, as a follow-up
job (`phylobot.properties_run` is the pure side).

The ancestor viewer's property panel lists what the registry can measure and this
ancestor lacks; its button posts here. The branch-site test is the precedent and the
shape is the same: a new `Job` with its own accounting and **its own results prefix**
(§18.3 — a job's files are its own, and deleting either job must not strip the other),
whose inputs are the parent's published ancestors, staged here, and whose params say what
was asked for (`properties_of`) with a dedup key, so asking twice leads to the one run.

Node matching is the web tier's (`phylobot.analogs` needs rooted trees and a job-wide
view); the pipeline is handed an explicit list of combinations and nodes and never
repeats it.

Nothing from the request reaches a file except as a choice among the parent job's own
values: every combination, node, gap version and property key is checked against what
that job published and against the registry.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
from typing import Any
from urllib.parse import urlencode

from django.db import transaction
from django.http import Http404, HttpRequest, HttpResponse
from django.shortcuts import get_object_or_404, redirect
from django.urls import reverse
from django.views.decorators.http import require_POST

from phylobot import (
    analog_views,
    analogs,
    ancestors,
    jobs,
    joined,
    properties,
    properties_run,
    property_views,
    storage,
    workload,
)
from phylobot.executors import PROPERTY_DIR, PROPERTY_TARGETS
from phylobot.gaps import MASKING
from phylobot.models import Job

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Published:
    """Enough of a finished job for `analog_views` to match nodes across its trees. It
    takes `Any` deliberately (the viewers' `AncestorSelection` is the other caller), so
    this does not drag the whole viewer in to read four fields."""

    job: Job
    summary: dict[str, Any]
    keys: set[str]
    results: joined.Results  # where each file lives (an extension's parent, for some)


def published(job: Job) -> Published:
    found = joined.results(job)
    if "summary.json" not in found.keys or not found.summary:
        raise Http404
    return Published(job, found.summary, found.keys, found)


def plan_key(
    parent: Job, targets: list[properties_run.Target], keys: list[str], version: str, draws: int
) -> str:
    """What makes two runs the same work: the same ancestors, properties, gap version
    and draw count, on the same parent. Asking twice must not queue it twice."""
    body = json.dumps(
        {
            "parent": str(parent.id),
            "targets": sorted((t.combination, t.node) for t in targets),
            "properties": sorted(keys),
            "gaps": version,
            "draws": draws,
        },
        sort_keys=True,
    )
    return hashlib.sha256(body.encode()).hexdigest()[:32]


def attempt_key(base: str, attempt: int) -> str:
    """The same work, drawn again. A recompute asks for the *same* plan with another
    seed, so the plan alone can no longer identify a run: the attempt number joins it,
    and attempt 0 is the first."""
    return base if attempt == 0 else f"{base}.{attempt}"


def attempts_of(owner: Any, base: str) -> int:
    """How many runs of this plan there have been, failed ones included: a failed draw
    still used its seed, and reusing it would reproduce the failure."""
    return Job.objects.for_user(owner).filter(params__properties_of__base=base).count()


def under_way(owner: Any, base: str) -> Job | None:
    """A run of this plan that has not finished. Both buttons check it, so a second
    click while the first run is queued leads to that run rather than to a second one —
    the one case where "recompute" plainly does not mean "again"."""
    return (
        Job.objects.for_user(owner)
        .filter(params__properties_of__base=base)
        .filter(status__in=[Job.Status.QUEUED, Job.Status.RUNNING])
        .order_by("created_at")
        .first()
    )


def measurements_of(runs: list[Job], combination: str, node: str, version: str) -> int:
    """How many of these runs measured this ancestor and gap version."""
    found = 0
    for run in runs:
        of = run.params.get("properties_of") or {}
        if of.get("gaps") != version:
            continue
        if node in (of.get("combinations") or {}).get(combination, []):
            found += 1
    return found


def next_seed(parent: Job, measured_already: int) -> int:
    """The seed the next run of this ancestor draws with: the parent's, moved on by
    every measurement this ancestor already has.

    Not the attempt number of the *plan*, which is the subtle one. PROPERTIES seeds each
    node as `<seed>:<node>:<version>`, so a run at the parent's seed reproduces the
    parent's draw exactly — down to the last value. A recompute asked for one ancestor
    after a four-ancestor run shares no plan with it, but would share its seed and hand
    back the same numbers under the name of a new draw. Counting the measurements this
    ancestor already has is what makes "draw again" true.
    """
    try:
        base = int(parent.params.get("seed") or 1)
    except TypeError, ValueError:
        base = 1
    return base + measured_already


def existing_run(owner: Any, key: str) -> Job | None:
    """The user's run of exactly this plan, if one has succeeded or is queued or running:
    the oldest, so every click leads to the same job. Failed and cancelled runs do not
    count, so a run that failed can be asked for again."""
    return (
        Job.objects.for_user(owner)
        .filter(params__properties_of__key=key)
        .exclude(status__in=[Job.Status.FAILED, Job.Status.CANCELLED])
        .order_by("created_at")
        .first()
    )


def runs_for(job: Job) -> list[Job]:
    """Property runs measured on this job's ancestors, newest first (§18.3: their results
    live under their own job, and the viewer reads both). With extensions
    (`phylobot.joined`), runs asked from a relative too: one measured the same files only
    where `measured_in` matches, which the reader checks per combination."""
    family = {str(j.id) for j in relatives(job)}
    return list(
        Job.objects.for_user(job.owner)
        .filter(params__properties_of__parent__in=family)
        .order_by("-created_at")
    )


def relatives(job: Job) -> list[Job]:
    """The job, the jobs it inherits results from, and those that inherit from it."""
    out, up = [job], job
    while (parent := joined.inheritance(up)[0]) is not None:
        out.append(parent)
        up = parent
    return out + joined.dependents(job)


def measured_in(run: Job, combination: str) -> str:
    """The id of the job whose files a properties run measured for `combination`: the
    one it recorded (runs asked from an extension), else the job it was asked from."""
    of = run.params.get("properties_of") or {}
    return str((of.get("owners") or {}).get(combination) or of.get("parent") or "")


def measurable(version: str) -> bool:
    """A gap version a posterior can be drawn from: IQ-TREE's reconstruction and its
    masked copies. PRANK's and GRASP's carry no posteriors."""
    return version in (ancestors.ORIGINAL, *MASKING)


def stage(
    parent: Job,
    job: Job,
    targets: list[properties_run.Target],
    version: str,
    files: joined.Results | None = None,
) -> None:
    """Copy what the pipeline needs into the new job's inputs, server-side.

    Per combination: the chosen version's ancestors, and the posteriors of the nodes being
    measured only, read through the parent's byte index (one ranged read per node), so
    staging one ancestor copies one ancestor's rows rather than the whole table.
    """
    dest = f"{job.storage_prefix}inputs/{PROPERTY_DIR}/"
    wanted: dict[str, list[str]] = {}
    for target in targets:
        wanted.setdefault(target.combination, []).append(target.node)
    rows = ["combination\tnode"]
    whole = "ancestors_posteriors.tsv"
    for combination, nodes in wanted.items():
        # An extension's inherited combinations are read from the job that holds them.
        source = files.prefix(combination) if files else f"{parent.storage_prefix}results/"
        fasta = ancestors.VERSIONS[version][1]
        storage.copy_key(f"{source}{combination}/{fasta}", f"{dest}{combination}/{fasta}")
        table = f"{source}{combination}/{whole}"
        try:
            index_text = storage.read_bytes(f"{source}{combination}/{ancestors.index_name(whole)}")
            index = ancestors.read_index(index_text.decode())
        except Exception:  # no index: the whole table, as it stands
            storage.copy_key(table, f"{dest}{combination}/{whole}")
        else:
            runs = [(0, index.header)]
            runs += [run for node in sorted(set(nodes)) for run in index.nodes.get(node, [])]
            kept = b"".join(storage.read_range(table, start, n) for start, n in runs)
            storage.put_bytes(f"{dest}{combination}/{whole}", kept)
        rows += [f"{combination}\t{node}" for node in sorted(nodes)]
    storage.put_bytes(
        f"{job.storage_prefix}inputs/{PROPERTY_TARGETS}", ("\n".join(rows) + "\n").encode()
    )


def context(job: Job) -> dict[str, Any]:
    """The job page's part for a properties run: what it measured, and the way back to
    the ancestor it measured. Empty for every other kind of job."""
    of = job.params.get("properties_of") or {}
    if not of:
        return {}
    parent_id = str(of.get("parent") or "")  # written by `submit`, so a valid UUID if set
    parent = Job.objects.for_user(job.owner).filter(id=parent_id).first() if parent_id else None
    viewer = ""
    if parent is not None:
        query = {
            "combination": of.get("combination", ""),
            "node": of.get("node", ""),
            "gaps": of.get("gaps", ancestors.ORIGINAL),
        }
        viewer = f"{reverse('ancestor_viewer', args=[parent.id])}?{urlencode(query)}"
    labels = [
        properties.REGISTRY[key].label
        for key in of.get("properties", [])
        if key in properties.REGISTRY
    ]
    return {"of": of, "parent": parent, "viewer_url": viewer, "labels": labels}


@require_POST
def submit(request: HttpRequest, job_id: str) -> HttpResponse:
    """Queue a properties run on one of the user's finished jobs."""
    parent = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if parent.status != Job.Status.SUCCEEDED:
        raise Http404
    at = published(parent)
    combos = {c["key"]: c for c in at.summary.get("combinations", [])}
    combination = combos.get(request.POST.get("combination", ""))
    if combination is None:
        raise Http404
    version = request.POST.get("gaps", ancestors.ORIGINAL)
    if not measurable(version):
        raise Http404
    fasta = ancestors.VERSIONS[version][1]
    if f"{combination['key']}/{fasta}" not in at.keys:
        raise Http404
    placed = analog_views.trees(at)
    key = combination["key"]
    node = request.POST.get("node", "")
    if key not in placed or node not in placed[key]:
        raise Http404

    chosen = [k for k in request.POST.getlist(properties_run.FIELD) if k in properties.REGISTRY]
    if not chosen:
        raise Http404
    scope = properties_run.scope_of(request.POST.get(properties_run.SCOPE_FIELD))
    found = analogs.find(placed, key, node, analog_views.by_clade(at))
    nodes = {combo: sorted(places) for combo, places in placed.items()}
    targets = properties_run.expand(scope, key, node, {a.combination: a.node for a in found}, nodes)
    # Only combinations that published this version: a target whose files are not there
    # would stage nothing and measure nothing (§18.4 — the offer must not promise it).
    targets = [t for t in targets if f"{t.combination}/{fasta}" in at.keys]
    if not targets:
        raise Http404

    # Raw request text: the viewer's own number field, so a silent clamp beats a 404.
    draws = properties.draws_in_range(request.POST.get("draws"))
    base = plan_key(parent, targets, chosen, version, draws)
    # "Recompute" means this plan again, so it needs a key of its own; anything else is
    # the first attempt, and clicking it twice leads to the one run.
    attempt = attempts_of(parent.owner, base) if request.POST.get("redraw") else 0
    plan = attempt_key(base, attempt)
    runs = runs_for(parent)
    measured_already = measurements_of(runs, key, node, version) + (
        1 if property_views.read(at.results.prefix(key), key, node, version) else 0
    )
    seed = next_seed(parent, measured_already)
    of = {
        "parent": str(parent.id),
        "node": node,
        "combination": key,
        "gaps": version,
        "scope": scope,
        "properties": chosen,
        "draws": draws,
        "ancestors": len(targets),
        "attempt": attempt,
        "seed": seed,
        # One list per combination, which is what the staged plan and the pipeline use;
        # it is also what the task count is read from on the job page.
        "combinations": {
            combo: sorted(t.node for t in targets if t.combination == combo)
            for combo in sorted({t.combination for t in targets})
        },
        "base": base,
        "key": plan,
        # Whose files each combination's ancestors are (an extension's inherited ones are
        # its parent's), so every job that reads those files shows the numbers.
        "owners": {
            combo: str(at.results.owner(combo).id)
            for combo in sorted({t.combination for t in targets})
        },
    }
    params: dict[str, Any] = {
        "properties_only": True,
        "properties": chosen,
        "property_draws": draws,
        # The pipeline measures the masked versions it is given; "original" is the
        # unmasked reconstruction and needs no gap method.
        "gap_methods": [] if version == ancestors.ORIGINAL else [version],
        "properties_of": of,
        "aligners": ["none"],
        "tree_builders": ["none"],
        "seed": seed,
    }
    label = (
        properties.REGISTRY[chosen[0]].label if len(chosen) == 1 else f"{len(chosen)} properties"
    )
    name = f"Ancestral properties: {label} on {node} in {parent.name or str(parent.id)[:8]}"
    size = workload.IME2_SIZE
    with transaction.atomic():
        # One submission at a time per parent, so two clicks queue one run.
        Job.objects.select_for_update().filter(id=parent.id).first()
        already = under_way(parent.owner, base) or existing_run(parent.owner, plan)
        if already is not None:
            return redirect("job_detail", job_id=already.id)
        job = Job.objects.create(
            owner=parent.owner,
            name=name[:100],
            params=params,
            work_units=workload.estimate(params, size).units,
            input_size=workload.size_record(size),
        )
        job.tags.set(parent.tags.all())
        stage(parent, job, targets, version, at.results)
        transaction.on_commit(lambda: jobs.enqueue(job))
    return redirect("job_detail", job_id=job.id)
