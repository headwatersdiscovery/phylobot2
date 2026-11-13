"""Ancestral libraries on the web (phylobot.library): the `Library` records, "Add to
library" from the tree and ancestor viewers, the top-bar chip, and the library page's targets.

Targets are collected across pages: the viewers post to the user's *current* library for
that job — the one they changed last, created on the first add — and redirect back. Every
posted value is checked against the job's own results (`views.select_ancestor`), so
nothing from the request is stored except as a choice among them.
"""

import dataclasses
import json
import re
import secrets
import uuid
from dataclasses import dataclass
from importlib import metadata
from typing import Any
from urllib.parse import urlencode

from django.contrib import messages
from django.db import models, transaction
from django.http import Http404, HttpRequest, HttpResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.utils import timezone
from django.utils.dateparse import parse_datetime
from django.views.decorators.http import require_POST

from phylobot import (
    about,
    alignment_views,
    ancestors,
    gaps,
    joined,
    library,
    library_disagreement,
    library_export,
    library_glyphs,
    library_index,
    results,
    storage,
    trees,
)
from phylobot.models import Job, Library, LibraryExport, Tag

# The axes the New Job page's defaults give more than one of (forms.DEFAULT_ALIGNERS...).
PRESET_AXES = ("aligners", "models", "gaps")


def phylobot_version() -> str:
    """What produced an export, recorded on its receipt."""
    try:
        return metadata.version("phylobot")[:40]
    except metadata.PackageNotFoundError:
        return "unknown"


def libraries_of(user: Any, job: Job) -> list[Library]:
    """The user's libraries on this job, the current one first."""
    return list(Library.objects.for_user(user).filter(job=job))


def current_library(user: Any, job: Job) -> Library | None:
    return Library.objects.for_user(user).filter(job=job).first()


def create_library(user: Any, job: Job) -> Library:
    number = Library.objects.filter(owner_id=user.pk, job=job).count() + 1
    return Library.objects.create(
        owner_id=user.pk,
        job=job,
        name=f"Library {number}",
        code=library.default_code(job.name, number),
        spec=library.new_spec(secrets.randbelow(2**31)),
        spec_version=library.SPEC_VERSION,
    )


def node_identity(selection: Any, node: str) -> tuple[frozenset[frozenset[str]], int]:
    """The node's partition (`gaps.node_partitions`, on the combination's tree as IQ-TREE
    wrote it) and how many sequences are below it on the tree as the viewers draw it."""
    by_label = {
        label: parts for parts, label in gaps.node_partitions(selection.read("tree.nwk")).items()
    }
    if node not in by_label:
        raise Http404
    drawn = trees.load_tree(selection.tree())
    below = next(
        (len(n.leaf_nodes()) for n in drawn.preorder_internal_node_iter() if n.label == node),
        0,
    )
    return by_label[node], below


def back_to(request: HttpRequest, job: Job) -> str:
    """Where to return after an action: the page it came from if that is one of this
    job's pages (never anywhere else), else the job page."""
    target = request.POST.get("next", "")
    if target.startswith(f"/jobs/{job.id}/") and "//" not in target and "\\" not in target:
        return target
    return reverse("job_detail", args=[job.id])


def owned_job(request: HttpRequest, job_id: str) -> Job:
    return get_object_or_404(Job.objects.for_user(request.user), id=job_id)


def owned_library(request: HttpRequest, job: Job, library_id: Any) -> Library:
    """One of the user's libraries on this job; 404 for anyone else's, a deleted one, or
    an id that is not one at all."""
    try:
        wanted = uuid.UUID(str(library_id))
    except ValueError:
        raise Http404 from None
    return get_object_or_404(Library.objects.for_user(request.user), id=wanted, job=job)


NEW = "new"


def destination(request: HttpRequest, job: Job, chosen: str) -> Library:
    """The library an action writes to: one of the user's on this job by id, a new one
    ("new"), or with nothing chosen the current one (created if there is none). The
    caller holds the job's row lock."""
    if chosen == NEW:
        return create_library(request.user, job)
    if chosen:
        return owned_library(request, job, chosen)
    found = current_library(request.user, job)
    return found if found is not None else create_library(request.user, job)


@require_POST
def add(request: HttpRequest, job_id: str) -> HttpResponse:
    """Add one or more ancestors of one reconstruction (POST combination, gaps, node...)
    to a library: the posted `library` (an id, or "new"), else the current one."""
    from phylobot.views import select_ancestor  # views imports this module

    job = owned_job(request, job_id)
    selection = select_ancestor(job, request.POST)
    if selection is None:
        raise Http404
    # select_ancestor falls back to defaults; here a value that is not the job's own is
    # an error, never a quiet substitution.
    if selection.combination["key"] != request.POST.get(
        "combination"
    ) or selection.version != request.POST.get("gaps"):
        raise Http404
    nodes = list(dict.fromkeys(request.POST.getlist("node")))
    known = set(selection.node_names) & set(selection.ancestors)
    if not nodes or any(node not in known for node in nodes):
        raise Http404
    entries = [
        library.target_entry(selection.combination["key"], node, selection.version, *identity)
        for node in nodes
        for identity in [node_identity(selection, node)]
    ]
    named = nodes[0] if len(nodes) == 1 else f"{len(nodes)} ancestors"
    return save_targets(request, job, entries, named)


def save_targets(
    request: HttpRequest, job: Job, entries: list[dict[str, object]], named: str
) -> HttpResponse:
    """Add `entries` to the posted library (an id, or "new"), else the current one, say
    so, and go back to the page the request came from."""
    with transaction.atomic():
        # One change at a time per job, so two quick clicks neither lose a target nor
        # create two libraries.
        Job.objects.select_for_update().filter(id=job.id).first()
        record = destination(request, job, request.POST.get("library", ""))
        record.spec, added = library.add_targets(record.spec, entries)
        record.save()
    if added == len(entries):
        messages.success(request, f"Added {named} to {record.name}.")
    elif added:
        messages.success(
            request, f"Added {added} of {named} to {record.name}; the rest were in it."
        )
    else:
        messages.info(request, f"{named[:1].upper()}{named[1:]} already in {record.name}.")
    return redirect(back_to(request, job))


def summary_of(job: Job) -> dict[str, Any]:
    """The job's summary.json, or {} if it has none (an unfinished or purged job)."""
    try:
        return dict(joined.cached(job).summary)  # an extension's: joined with its parent's
    except ValueError:
        return {}


@dataclass(frozen=True)
class Loaded:
    """A target's reconstruction as stored in the job: its aligned row and, where the
    version has them, its per-column distributions."""

    template: str
    dists: list[list[tuple[str, float]]] | None
    alignment: dict[str, str]  # the combination's extant rows, the template's columns


def load_targets(job: Job, targets: list[dict[str, Any]], base: Any = None) -> list[Loaded | None]:
    """Each target's reconstruction, or None where the job no longer has it. The job's
    files are listed once (`base`, an `AncestorSelection` of the job, if the caller has
    one), then each combination and version is read from that listing
    (`AncestorSelection.other`); only an exact match counts, never a fallback, which
    would show another ancestor."""
    from phylobot.views import select_ancestor  # views imports this module

    if base is None and targets:
        base = select_ancestor(job, {})
    selections: dict[tuple[str, str], Any] = {}
    alignments: dict[str, dict[str, str]] = {}
    loaded: list[Loaded | None] = []
    for target in targets:
        key, node, version = library.target_key(target)
        if (key, version) not in selections:
            selections[key, version] = None if base is None else base.other(key, version)
        selection = selections[key, version]
        if selection is None or node not in selection.ancestors:
            loaded.append(None)
            continue
        template = selection.ancestors[node]
        dists = None
        table = "ancestors_posteriors.tsv"
        if ancestors.has_posteriors(version) and (
            f"{key}/{ancestors.version_folder(version)}{table}" in selection.keys
        ):
            dists = ancestors.distributions(selection.read_node(table, node), node)
            if len(dists) != len(template):  # never seen; drawn uncoloured rather than wrong
                dists = None
        if key not in alignments:
            present = f"{key}/alignment.fasta" in selection.keys
            alignments[key] = (
                ancestors.parse_fasta(selection.read("alignment.fasta")) if present else {}
            )
        loaded.append(Loaded(template, dists, alignments[key]))
    return loaded


def paralogy_flags(job: Job, targets: list[dict[str, Any]]) -> list[dict[str, str] | None]:
    """Each target's reconciliation flag, as the ancestor viewer's note says it
    (`views.paralogy_context`): the job's two trees disagree about it, or it sits at a
    (not weak) duplication. Only reconciled jobs are read (one more listing of the job)."""
    from phylobot.executors import reconciles
    from phylobot.views import paralogy_context, select_ancestor  # views imports this module

    if not targets or not reconciles(job.params or {}):
        return [None] * len(targets)
    base = select_ancestor(job, {})
    selections: dict[str, Any] = {}
    flags: list[dict[str, str] | None] = []
    for target in targets:
        key, node, _ = library.target_key(target)
        if key not in selections:
            selections[key] = None if base is None else base.other(key, ancestors.ORIGINAL)
        selection = selections[key]
        found = paralogy_context(selection, node) if selection is not None else None
        flag = None
        if found and found["disagree"]:
            flag = {"label": "trees disagree", "title": (
                f"{found['own_label'].capitalize()} and {found['other_label']} group this "
                "ancestor's sequences differently: treat it as uncertain, and compare it "
                "with the other tree's ancestors.")}  # fmt: skip
        elif (
            found
            and found["found"].duplication
            and not (found["overlap"] and found["overlap"].weak)
            and found["found"].reconciled_event in ("", "duplication")
        ):
            flag = {"label": "pre-duplication", "title": (
                "A gene duplication at this node: the gene before its copies split.")}  # fmt: skip
        flags.append(flag)
    return flags


def extant_lengths(loaded: list[Loaded | None]) -> list[int]:
    """The extant sequences' lengths, which both tables and the length chart compare with."""
    return library_glyphs.extant_lengths([data.alignment for data in loaded if data is not None])


@dataclass
class Viewers:
    """The library's alignment viewers and the layout menus' state."""

    frames: list[tuple[library_glyphs.Frame, dict[str, Any]]]  # each with its viewer data
    layout: library_glyphs.Layout
    alignments: list[tuple[str, str]]
    extant: list[str]
    # Each drawable target's entry, by its place in the library: what its members are
    # drawn on (`member_viewers`).
    entries: dict[int, library_glyphs.Entry] = dataclasses.field(default_factory=dict)
    # Cross-combination disagreement across the targets (phylobot.library_disagreement).
    disagreement: library_disagreement.Overview | None = None


def target_rows(
    job: Job,
    record: Library,
    threshold: float,
    palette: str,
    alignment: str = "",
    extant: list[str] | None = None,
    cell: float | str = "fit",
    scheme: str = library_glyphs.DEFAULT_SCHEME,
    loaded: list[Loaded | None] | None = None,
) -> tuple[list[dict[str, Any]], Viewers]:
    """The targets as the library page lists them — links back to each viewer and the
    sequence as a site map — and the alignment viewers across them all, laid out as asked
    (`library_glyphs.resolve_layout` replaces anything unknown)."""
    combos = {c["key"]: c for c in summary_of(job).get("combinations", [])}
    ancestor_url = reverse("ancestor_viewer", args=[job.id])
    tree_url = reverse("tree_viewer", args=[job.id])
    targets = record.targets
    if loaded is None:
        loaded = load_targets(job, targets)
    in_range = library_glyphs.extant_range(extant_lengths(loaded))
    flags = paralogy_flags(job, targets)
    rows = []
    entries: list[library_glyphs.Entry] = []
    by_index: dict[int, library_glyphs.Entry] = {}
    # The targets grouped by clade (the stored partition), for cross-combination disagreement.
    clade_states: dict[Any, list[library_disagreement.State]] = {}
    row_clade: dict[int, Any] = {}
    for index, (target, data) in enumerate(zip(targets, loaded, strict=True)):
        key, node, version = library.target_key(target)
        combo = combos.get(key, {})
        here = {"aligner": combo.get("aligner", ""), "combination": key, "node": node}
        anchor = f"target-{index + 1}"
        row: dict[str, Any] = {
            "anchor": anchor,
            "key_json": json.dumps([key, node, version]),
            "combination": key,
            "combination_label": combo.get("label", key),
            "node": node,
            "version": version,
            "version_label": ancestors.version_label(version),
            "version_short": ancestors.short_version_label(version),
            "below": target.get("below") or "",
            "missing": data is None,
            "ancestor_url": f"{ancestor_url}?{urlencode(here | {'gaps': version})}",
            "tree_url": f"{tree_url}?{urlencode(here)}",
            "paralogy": flags[index],
        }
        if data is not None:
            site_map = library_glyphs.site_map(data.template, data.dists, threshold)
            sequence_url = reverse("library_sequence", args=[job.id, record.id])
            row |= {
                "map": site_map,
                "length_warning": in_range.check(site_map.length) if in_range else None,
                "fasta": site_map.fasta(f"{node} {combo.get('label', key)} {version}"),
                "map_url": f"{sequence_url}?{urlencode(here | {'version': version})}",
            }
            label = node if version == ancestors.ORIGINAL else f"{node} · {version}"
            row["short_label"] = label  # as the viewers' row labels name it
            title = f"{node} · {row['combination_label']} · {row['version_label']}"
            aligner = combo.get("aligner", key)
            entries.append(
                library_glyphs.Entry(
                    label, anchor, title, data.template, data.dists, aligner,
                    combo.get("aligner_label", aligner), data.alignment,
                )
            )  # fmt: skip
            by_index[index] = entries[-1]
            partition = target.get("partition")
            if partition is not None:
                clade = tuple(tuple(part) for part in partition)
                clade_states.setdefault(clade, []).append(
                    library_disagreement.State(key, data.template, data.alignment)
                )
                row_clade[index] = clade
        rows.append(row)
    disputes, overview = library_disagreement.analyse(clade_states)
    for index, row in enumerate(rows):
        dispute = disputes.get(row_clade.get(index))
        row["disputed"] = dispute.per_key.get(row["combination"]) if dispute else None
    layout = library_glyphs.resolve_layout(entries, alignment, extant or [])
    least = max(library_glyphs.threshold_percent(threshold), library_glyphs.FLOOR)
    built = library_glyphs.frames(entries, layout) if entries else []
    # The in-viewer disagreement track: on a shared-columns layout, mark the reference frame's
    # columns where an ancestor's combinations give different ML residues (library_disagreement).
    track = None
    if layout.alignment not in (library_glyphs.OWN, library_glyphs.ON_SEQUENCE):
        reference = next((e for e in entries if e.aligner == layout.alignment and e.fits), None)
        if reference is not None and clade_states:
            width = len(reference.template)
            disputed: set[int] = set()
            for states in clade_states.values():
                disputed |= library_disagreement.disputed_columns(states, reference.alignment)
            if disputed:
                values = [1.0 if c in disputed else 0.0 for c in range(width)]
                track = (width, library_glyphs.disagreement_track(values, "column"))
    frames = []
    for i, frame in enumerate(built):
        frame = dataclasses.replace(frame, index=i)
        if track is not None and frame.units == track[0]:
            frame = dataclasses.replace(frame, tracks=[*frame.tracks, track[1]])
        frames.append(frame)
    viewers = Viewers(
        [(f, f.data(palette, least, cell, scheme)) for f in frames],
        layout,
        library_glyphs.alignment_options(entries) if entries else [],
        library_glyphs.extant_options(entries),
        by_index,
        overview,
    )
    return rows, viewers


MEMBERS_SHOWN = {"selected": "The selected set", "all": "Every member"}
DEFAULT_MEMBERS_SHOWN = "selected"
MEMBER_ROWS = 24  # rows in view before a members' viewer scrolls (a plate's row is 12)


def parse_members_shown(text: object) -> str:
    key = str(text)
    return key if key in MEMBERS_SHOWN else DEFAULT_MEMBERS_SHOWN


def member_viewers(
    members: list[dict[str, Any]],
    viewers: Viewers,
    shown: str,
    palette: str,
    threshold_pct: int,
    cell: float | str = "fit",
    scheme: str = library_glyphs.DEFAULT_SCHEME,
    alignment: str = "",
) -> tuple[list[tuple[library_glyphs.Frame, dict[str, Any]]], str]:
    """The members as rows of alignment viewers (L5c), laid out as the ancestors' viewers
    are, so a column means the same thing in both. Each member sits beneath the members
    of its own ancestor (its representative source), the selected ones first, coloured by
    that ancestor's posterior for the state it carries and outlined where it departs from
    the ML sequence. A member can only be drawn on its ancestor's columns, which is why it
    is placed through its ancestor's entry and never onto another target's frame.

    `members` are `design_context`'s rows, already in the table's order. `alignment` is the
    ancestors' layout ("Map ancestors onto"), passed so one control maps both viewers onto
    the same columns; returns the drawn frames and the resolved alignment."""
    groups: dict[int, list[library_glyphs.Entry]] = {}
    for row in members:
        if shown == "selected" and not row["picked"]:
            continue
        member = row["member"]
        index = member.source.combination - 1  # `design_targets` numbers targets from 1
        target = viewers.entries.get(index)
        if target is None or len(member.sequence) != len(target.template):
            continue
        pairs = zip(target.template, member.sequence, strict=True)
        marked = frozenset(i for i, (was, now) in enumerate(pairs) if was != now)
        groups.setdefault(index, []).append(
            dataclasses.replace(
                target, label=row["name"], anchor=row["anchor"], template=member.sequence,
                title=row["name"], member=True, marked=marked, cut=not row["picked"],
            )
        )  # fmt: skip
    entries = [e for index in sorted(groups) for e in groups[index]]
    if not entries:
        return [], library_glyphs.OWN
    layout = library_glyphs.resolve_layout(entries, alignment, list(viewers.layout.extant))
    out = []
    for i, frame in enumerate(library_glyphs.frames(entries, layout)):
        frame = dataclasses.replace(frame, index=i, group="m")
        data = frame.data(palette, threshold_pct, cell, scheme)
        data["height"] = min(len(frame.rows), MEMBER_ROWS)
        out.append((frame, data))
    return out, layout.alignment


# Above this many members the chart stops looking past the chosen budget: spanning costs
# a distance per member per pick, and a user with hundreds of single-site alternates would
# otherwise wait seconds to see a curve they did not ask for.
CHART_POOL_LIMIT = 250
STRATEGY_LABELS = {
    "supported": "Most supported",
    "spanning": "Spanning",
    "per_node": "Per ancestor",
}


def combination_weights(job: Job) -> dict[str, float]:
    """Each combination's BIC weight from `model_comparison.tsv`: how much of the fit the
    models give it, which is what makes a member supported by well-fitting choices outrank
    one supported by poor ones. One small file per job; {} for a job run before it existed.
    Weights are comparable within an alignment only, which is how `model_comparison.tsv`
    itself ranks them."""
    files = joined.cached(job)
    if "model_comparison.tsv" not in files.keys:
        return {}
    weights: dict[str, float] = {}
    for row in results.read_tsv(files.read("model_comparison.tsv")):
        try:
            weights[f"{row['aligner']}/{row['tree']}/{row['model_key']}"] = float(row["w_bic"])
        except KeyError, TypeError, ValueError:
            continue
    return weights


def design_targets(
    job: Job, record: Library, threshold: float, loaded: list[Loaded | None] | None = None
) -> tuple[list[library.Target], dict[tuple[int, str, str], str]]:
    """The library's targets as the generator wants them. Only reconstructions with
    posteriors can be drawn from: PRANK's and GRASP's have none, so they are targets the
    page draws but not ones a variant can be made of."""
    entries = record.targets
    weights = combination_weights(job)
    ancestor_url = reverse("ancestor_viewer", args=[job.id])
    out: list[library.Target] = []
    urls: dict[tuple[int, str, str], str] = {}
    if loaded is None:
        loaded = load_targets(job, entries)
    for index, (entry, data) in enumerate(zip(entries, loaded, strict=True), start=1):
        if data is None or data.dists is None:
            continue
        key, node, version = library.target_key(entry)
        gap = "" if version == ancestors.ORIGINAL else version
        source = library.Source(
            combination=index, node=node, gap=gap, weight=weights.get(key, 0.0), label=key
        )
        out.append(library.Target(source=source, template=data.template, dists=data.dists))
        here = {"combination": key, "node": node, "gaps": version}
        urls[source.key] = f"{ancestor_url}?{urlencode(here)}"
    return out, urls


def posted_design(request: HttpRequest, record: Library, threshold: float) -> library.Design:
    """The design to show: the saved one, with anything in the address on top of it. The
    URL carries the controls until Save writes them, as every other viewer's do."""
    design = library.design_from(record.spec, threshold)
    axes = design.axes
    get = request.GET
    if "samples" in get:
        axes = dataclasses.replace(axes, samples=library.clamp_int(get.get("samples"), 0, 10_000))
    if "subs" in get:
        axes = dataclasses.replace(
            axes, max_substitutions=library.clamp_int(get.get("subs"), 0, library.MAX_SUBSTITUTIONS)
        )
    # "Include AltAll" is affirmative and on by default; a hidden "off" field accompanies the
    # box, so an absent/unticked box reads "off". The threshold param (`altall`) rides along
    # on every design submit, so its presence is the signal the design form was used.
    if "altall" in get or "altall_on" in get:
        axes = dataclasses.replace(axes, alt_all=get.get("altall_on") != "off")
    budget = design.budget
    if "budget" in get:
        budget = library.clamp_int(get.get("budget"), 1, library.MAX_BUDGET)
    strategy = get.get("strategy", design.strategy)
    return library.Design(
        axes=axes,
        budget=budget,
        strategy=strategy if strategy in library.STRATEGIES else design.strategy,
    )


BUDGET_SLIDER_BASE = 192  # the slider's default reach


def budget_slider_max(budget: int) -> int:
    """The slider's upper bound: 192 by default, or twice the budget when it exceeds 192, so
    the thumb never pins at the end (library.js mirrors this live as the user types)."""
    return budget * 2 if budget > BUDGET_SLIDER_BASE else BUDGET_SLIDER_BASE


def design_context(
    job: Job,
    record: Library,
    design: library.Design,
    threshold: float,
    loaded: list[Loaded | None] | None = None,
) -> dict[str, Any]:
    """The Design section: the live count, the members the budget takes, and coverage.

    `count_candidates` runs first and refuses over `CANDIDATE_LIMIT` before anything is
    generated, so a wide library says so instead of hanging (the New Job page's pattern).
    The targets' posteriors are already read for the page's glyphs, so generating adds
    arithmetic, not I/O.
    """
    if loaded is None:
        loaded = load_targets(job, record.targets)
    targets, urls = design_targets(job, record, threshold, loaded)
    plain_strategies, advanced_strategies, advanced_open = library.strategy_choices(design.strategy)
    out: dict[str, Any] = {
        "design": design,
        "axes": design.axes,
        # The fill strategies, split into the plain one and the two behind "Advanced"
        # (library.strategy_choices; the disclosure opens when a saved spec uses one).
        "strategies": plain_strategies,
        # The "?" popup beside "Fill them with": built from the About page's own
        # citation objects (about.strategy_help), so the two can never cite differently.
        "strategy_help": about.strategy_help(),
        "advanced_strategies": advanced_strategies,
        "advanced_open": advanced_open,
        "sample_presets": library.SAMPLE_PRESETS,
        "max_substitutions": library.MAX_SUBSTITUTIONS,
        "max_budget": library.MAX_BUDGET,
        # The slider reaches 192 by default (twice the budget once it is set higher); the
        # number field keeps the full range for anyone ordering more than that.
        "budget_slider_max": budget_slider_max(design.budget),
        "candidate_limit": library.CANDIDATE_LIMIT,
        "designable": len(targets),
        "ambiguous": sum(len(library.ambiguous_swaps(t, design.axes.threshold)) for t in targets),
        "members": [],
        "too_many": 0,
    }
    if not targets:
        return out
    candidates = library.count_candidates(targets, design.axes)
    out["candidates"] = candidates
    if candidates > library.CANDIDATE_LIMIT:
        out["too_many"] = candidates
        return out
    seed = record.spec.get("seed", 0)
    members = library.reduce(library.generate(targets, design.axes, seed=seed))
    # The job's default gap method leads the controls, so a small library holds the
    # masked ancestor the user would have ordered, not the unmasked one.
    preferred = gaps.default_methods(bool(job.params.get("outgroup")))
    picked = library.select(
        members, budget=design.budget, strategy=design.strategy, preferred_gaps=preferred
    )
    chosen = {id(m) for m in picked.members}
    coverage = library.coverage(picked.members, targets, threshold=design.axes.threshold)
    by_source = {t.source.key: t for t in targets}
    in_range = library_glyphs.extant_range(extant_lengths(loaded))
    rows = []
    for number, member in enumerate(members, start=1):
        target = by_source.get(member.source.key)
        template = target.template if target else None
        comparable = template is not None and len(template) == len(member.sequence)
        rows.append({
            "member": member,
            "name": member.name(record.code),
            "variant_label": library.variant_label(member.variant, member.kind),
            "fits": library.name_fits(member.name(record.code)),
            "picked": id(member) in chosen,
            "differences": len(member.differences(template)) if template and comparable else None,
            "mean_pp": member.mean_posterior(target) if target else None,
            "length_warning": in_range.check(member.length) if in_range else None,
            "source_url": urls.get(member.source.key, ""),
            "anchor": f"member-{number}",  # its row in the table and in the viewer
        })  # fmt: skip
    rows.sort(key=lambda r: (not r["picked"],))  # the selected set first, order kept
    # What the budget and the strategy buy, drawn over the same axes. Every strategy's
    # selection at budget k is the first k of its picks, so one run each gives the whole
    # curve. Spanning is the costly one (a distance per member per pick), so a wide pool
    # charts only as far as it is worth looking.
    chart_max = min(len(members), max(design.budget * 2, 48))
    if len(members) > CHART_POOL_LIMIT:
        chart_max = min(chart_max, design.budget)
    curves = {}
    for key in library.STRATEGIES:
        order = library.select(members, budget=chart_max, strategy=key, preferred_gaps=preferred)
        curves[key] = (STRATEGY_LABELS[key], library.coverage_curve(
            order.members, targets, threshold=design.axes.threshold
        ))  # fmt: skip
    chart = library_glyphs.coverage_chart(curves, coverage.total, design.budget)
    out |= {
        "budget_slider_max": budget_slider_max(design.budget),
        "chart": chart,
        "chart_capped": len(members) > CHART_POOL_LIMIT,
        "members": rows,
        "distinct": len(members),
        "selection": picked,
        "coverage": coverage,
        "dropped": picked.dropped,
    }
    return out


@require_POST
def design(request: HttpRequest, job_id: str, library_id: str) -> HttpResponse:
    """Save the design into the spec. Everything it holds is a number or one of a fixed
    set of keys, checked by `library.design_from` on the way back out."""
    job = owned_job(request, job_id)
    with transaction.atomic():
        Job.objects.select_for_update().filter(id=job.id).first()
        record = owned_library(request, job, library_id)
        saved = design_from_post(request)
        view = library.View(
            palette=request.POST.get("palette", ""),
            scheme=request.POST.get("scheme", ""),
            zoom=request.POST.get("zoom", ""),
            alignment=request.POST.get("aln", ""),
            extant=tuple(request.POST.getlist("extant")),
            members=request.POST.get("members", ""),
        )
        record.spec = library.with_design(
            record.spec, saved, view, timezone.now().isoformat(timespec="seconds")
        )
        record.save()
    messages.success(request, "Saved the design and the view.")
    return redirect(back_to(request, job))


def design_from_post(request: HttpRequest) -> library.Design:
    """The design a form sent. Save and Export read it the same way, so an export freezes
    what the page was showing rather than whatever was saved last."""
    threshold = library.parse_threshold(request.POST.get("altall", ""))
    strategy = request.POST.get("strategy", "")
    return library.Design(
        axes=library.Axes(
            alt_all=request.POST.get("alt_all") == "1",
            threshold=threshold,
            samples=library.clamp_int(request.POST.get("samples"), 0, 10_000),
            max_substitutions=library.clamp_int(
                request.POST.get("subs"), 0, library.MAX_SUBSTITUTIONS
            ),
        ),
        budget=library.clamp_int(request.POST.get("budget"), 1, library.MAX_BUDGET),
        strategy=strategy if strategy in library.STRATEGIES else library.DEFAULT_STRATEGY,
    )


def export_rows(
    job: Job, record: Library, design: library.Design, threshold: float
) -> tuple[list[library_export.Row], dict[str, Any]]:
    """The members an export writes, and what the methods paragraph needs. Built from the
    same generator the page shows, so an export is what was on screen."""
    context = design_context(job, record, design, threshold)
    combos = {c["key"]: c for c in summary_of(job).get("combinations", [])}
    entries = record.targets
    info: dict[tuple[int, str, str], library_export.Source] = {}
    templates: dict[tuple[int, str, str], str] = {}
    weights = combination_weights(job)
    index = 0
    for entry, data in zip(entries, load_targets(job, entries), strict=True):
        index += 1
        if data is None or data.dists is None:
            continue
        key, node, version = library.target_key(entry)
        gap = "" if version == ancestors.ORIGINAL else version
        combo = combos.get(key, {})
        info[index, node, gap] = library_export.Source(
            combination=key,
            label=combo.get("label", key),
            node=node,
            version=version,
            version_label=ancestors.version_label(version),
            short_label=ancestors.short_version_label(version),
            weight=weights.get(key, 0.0),
        )
        templates[index, node, gap] = data.template
    rows = []
    for row in context.get("members", []):
        member = row["member"]
        source = info.get(member.source.key)
        if source is None:
            continue
        rows.append(
            library_export.Row(
                name=row["name"],
                member=member,
                source=source,
                template=templates[member.source.key],
                selected=bool(row["picked"]),
                mean_pp=row["mean_pp"],
            )
        )
    return rows, context


@require_POST
def export(request: HttpRequest, job_id: str, library_id: str) -> HttpResponse:
    """Freeze the current design: write the four documents and a copy of the spec under
    the library's own prefix, and leave a `LibraryExport` receipt. Nothing is overwritten
    — each export is numbered, so what was ordered can always be read back."""
    job = owned_job(request, job_id)
    record = owned_library(request, job, library_id)
    design = design_from_post(request)
    threshold = design.axes.threshold
    rows, context = export_rows(job, record, design, threshold)
    if not rows:
        messages.error(request, "There is nothing to export yet: this library has no members.")
        return redirect(back_to(request, job))

    coverage = context.get("coverage")
    methods = library_export.Methods(
        library_name=record.name,
        code=record.code,
        job_name=job.name or "Job",
        design=design,
        sources=sorted({row.source for row in rows}, key=lambda s: (s.combination, s.node)),
        members=context.get("distinct", len(rows)),
        selected=sum(1 for row in rows if row.selected),
        ambiguous=context.get("ambiguous", 0),
        coverage=coverage.sentence if coverage else "",
        seed=int(record.spec.get("seed", 0) or 0),
        version=phylobot_version(),
    )
    with transaction.atomic():
        Job.objects.select_for_update().filter(id=job.id).first()
        number = (record.exports.aggregate(models.Max("number"))["number__max"] or 0) + 1
        receipt = LibraryExport.objects.create(
            library=record,
            number=number,
            spec=library.with_design(record.spec, design, library.view_from(record.spec)),
            spec_version=record.spec_version,
            member_count=methods.selected,
            phylobot_version=methods.version,
        )
        documents = library_export.documents(
            rows, methods, json.dumps(receipt.spec, indent=1, sort_keys=True)
        )
        for name, text in documents.items():
            storage.put_bytes(f"{receipt.storage_prefix}{name}", text.encode())
    messages.success(
        request, f"Export {number}: {methods.selected} sequence(s) written, with a receipt."
    )
    return redirect(back_to(request, job))


def export_file(
    request: HttpRequest, job_id: str, library_id: str, number: int, name: str
) -> HttpResponse:
    """One file of one export, as it was written."""
    if name not in library_export.FILES:
        raise Http404
    job = owned_job(request, job_id)
    record = owned_library(request, job, library_id)
    receipt = get_object_or_404(LibraryExport, library=record, number=number)
    key = f"{receipt.storage_prefix}{name}"
    if key not in storage.list_keys(key):
        raise Http404
    stem = library_export.safe_name(record.name)
    response = HttpResponse(
        storage.read_bytes(key),
        content_type="text/plain; charset=utf-8",
    )
    response["Content-Disposition"] = f'attachment; filename="{stem}-{number}-{name}"'
    return response


@require_POST
def restore(request: HttpRequest, job_id: str, library_id: str, number: int) -> HttpResponse:
    """Put an export's design back into the editable library. Its targets are restored
    too: a design means nothing without the ancestors it was drawn from."""
    job = owned_job(request, job_id)
    with transaction.atomic():
        Job.objects.select_for_update().filter(id=job.id).first()
        record = owned_library(request, job, library_id)
        receipt = get_object_or_404(LibraryExport, library=record, number=number)
        kept_seed = record.spec.get("seed")
        record.spec = dict(receipt.spec) | ({"seed": kept_seed} if kept_seed else {})
        record.save()
    messages.success(request, f"Restored the settings of export {number}.")
    return redirect("library_page", job_id=job.id, library_id=record.id)


@require_POST
def duplicate(request: HttpRequest, job_id: str, library_id: str) -> HttpResponse:
    """Another library on this job with the same ancestors and the same design: the way to
    keep several configurations side by side (owner decision 2026-10-08), rather than a
    second kind of saved thing inside a library. The seed is copied too, so the draws are
    the same and the only difference is what the copy is changed to."""
    job = owned_job(request, job_id)
    with transaction.atomic():
        Job.objects.select_for_update().filter(id=job.id).first()
        record = owned_library(request, job, library_id)
        copy = Library.objects.create(
            owner_id=record.owner_id,  # the library's own owner, already checked above
            job=job,
            name=f"{record.name} (copy)"[:100],
            code=record.code,
            spec=dict(record.spec),
            spec_version=record.spec_version,
        )
    messages.success(request, f"Copied {record.name} to {copy.name}.")
    return redirect("library_page", job_id=job.id, library_id=copy.id)


def library_page(request: HttpRequest, job_id: str, library_id: str) -> HttpResponse:
    """A library's targets: the overview glyph, and per ancestor its site map. The
    threshold (AltAll's) and the palette are GET parameters, as in the viewers.
    The design controls, members and preview join it in L5."""
    job = owned_job(request, job_id)
    record = owned_library(request, job, library_id)
    # A saved view is where the page starts; anything in the address still wins, so a
    # shared link shows what its sender saw.
    saved_view = library.view_from(record.spec)
    threshold = library.parse_threshold(
        request.GET.get("altall", ""), default=library.design_from(record.spec).axes.threshold
    )
    palette = request.GET.get("palette") or saved_view.palette
    if palette not in ancestors.PALETTES:
        palette = ancestors.DEFAULT_PALETTE
    zoom = library_glyphs.parse_zoom(request.GET.get("zoom") or saved_view.zoom or "fit")
    scheme = library_glyphs.parse_scheme(request.GET.get("scheme") or saved_view.scheme)
    shown = parse_members_shown(request.GET.get("members") or saved_view.members)
    # Both halves of the page draw on the same reconstructions: read them once.
    loaded = load_targets(job, record.targets)
    rows, viewers = target_rows(
        job, record, threshold, palette, request.GET.get("aln") or saved_view.alignment,
        request.GET.getlist("extant") or list(saved_view.extant), zoom, scheme, loaded,
    )  # fmt: skip
    design = posted_design(request, record, threshold)
    # "I can order" stays with the library as soon as it is changed, not only on Save: it
    # used to live in the address alone, and a reload from the Controls panel (which sends
    # the view, not the design) put the saved number back.
    if "budget" in request.GET and design.budget != library.design_from(record.spec).budget:
        with transaction.atomic():
            Job.objects.select_for_update().filter(id=job.id).first()
            record = owned_library(request, job, library_id)
            record.spec = dict(record.spec) | {"budget": design.budget}
            record.save(update_fields=["spec"])
    designed = design_context(job, record, design, threshold, loaded)
    # A library that spans several aligners offers the Design section its own members-layout
    # menu ("Show separate" or on one aligner's columns, where the agree/disagree tracks appear);
    # the ancestors use the Controls panel's "Compare ancestors" (aln) control.
    layout_aligners = [
        (key, label)
        for key, label in viewers.alignments
        if key not in (library_glyphs.OWN, library_glyphs.ON_SEQUENCE)
    ]
    threshold_pct = max(library_glyphs.threshold_percent(threshold), library_glyphs.FLOOR)
    # The members follow the ancestors' layout ("Map ancestors onto" in the Controls panel):
    # one control maps both the ancestors' and the members' viewers onto the chosen alignment.
    member_frames, _ = member_viewers(
        designed["members"], viewers, shown, palette, threshold_pct, zoom, scheme,
        viewers.layout.alignment,
    )  # fmt: skip
    # Name each member's source combination readably (MAFFT, PRANK…) for the table, so rows
    # from different aligners are told apart; the viewer names them through `_named`.
    for row in designed.get("members", []):
        entry = viewers.entries.get(row["member"].source.combination - 1)
        row["aligner_label"] = entry.aligner_label if entry else ""
    context = {
        "job": job,
        "library": record,
        "multi_aligner": len(layout_aligners) > 1,
        # Axes with one choice (or none): the page suggests a rerun with the defaults (L6).
        "thin_axes": [
            axis
            for axis in library_index.axes(job.params or {})
            if axis.thin and axis.key in PRESET_AXES
        ],
        "rows": rows,
        "length_chart": library_glyphs.length_chart(
            extant_lengths(loaded),
            [
                (
                    r["short_label"],
                    f"{r['node']} ({r['combination_label']}, {r['version_label']})",
                    r["map"].length,
                )
                for r in rows
                if not r["missing"]
            ],
        ),
        "viewers": viewers,
        "palette_stops": {key: library_glyphs.palette_stops(key) for key in ancestors.PALETTES},
        "zoom": zoom,
        "threshold": threshold,
        "threshold_pct": threshold_pct,
        "colour_classes": library_glyphs.colour_classes(palette),
        "no_posteriors": library_glyphs.NO_POSTERIORS,
        "threshold_min": library.MIN_THRESHOLD,
        "threshold_max": library.MAX_THRESHOLD,
        "palette": palette,
        "palettes": [(k, v[0]) for k, v in ancestors.PALETTES.items()],
        "scheme": scheme,
        "saved_at": library.saved_at(record.spec),
        "saved_at_display": parse_datetime(library.saved_at(record.spec)),
        # A design saved before the page recorded the moment still counts as saved.
        "has_design": library_index.designed(record.spec),
        "exports": library_export.listing(list(record.exports.all())),
        "schemes": list(library_glyphs.SCHEMES.items()),
        "by_posterior": scheme == library_glyphs.BY_POSTERIOR,
        "scheme_legend": alignment_views.legend(scheme),
        # Sequences open before a palette change (the page's script sends them).
        "open_rows": {a for a in request.GET.getlist("open") if a in {r["anchor"] for r in rows}},
        "legend": [(pp, ancestors.colour(pp, palette)) for pp in (0.0, 0.25, 0.5, 0.75, 1.0)],
        **designed,
        "member_frames": member_frames,
        "members_shown": shown,
        "members_shown_options": list(MEMBERS_SHOWN.items()),
        "others": [lib for lib in libraries_of(request.user, job) if lib.id != record.id],
        # Libraries on the user's other jobs cannot be destinations: a target names a
        # combination, node and version of *this* job's results, so copying one across
        # would store a target the other job cannot resolve. The dialog says so rather
        # than leaving the list looking short.
        "elsewhere": Library.objects.for_user(request.user)
        .filter(job__deleted_at__isnull=True)
        .exclude(job=job)
        .count(),
        "code_limit": library.CODE_LIMIT,
        "notes_limit": Library.NOTES_LIMIT,
        "excel_url": f"{request.path}?{urlencode(excel_query(request), doseq=True)}",
    }
    if request.GET.get("export") == "xlsx":
        return members_workbook(context)
    return render(request, "phylobot/library.html", context)


def excel_query(request: HttpRequest) -> dict[str, list[str]]:
    """The page's own address with `export=xlsx`: the download is what the page shows.
    (library.js rebuilds it from the address on click, since palette, zoom and the design
    change the address without a reload.)"""
    query = {k: request.GET.getlist(k) for k in request.GET if k != "export"}
    return query | {"export": ["xlsx"]}


def members_workbook(context: dict[str, Any]) -> HttpResponse:
    """The members viewer of "Design the set" as an Excel workbook (phylobot.sites_excel),
    as the tree viewer's "Download Excel" does for a mutational path."""
    from phylobot import sites_excel  # openpyxl: only needed here

    job, record = context["job"], context["library"]
    design = context["design"]
    about = [
        ("Library", f"{record.name} ({record.code})"),
        ("Job", f"{job.name or '(unnamed)'} ({job.id})"),
        ("Uncertainty threshold", f"{context['threshold']:.2f}"),
        ("Budget", str(design.budget)),
        ("Strategy", STRATEGY_LABELS.get(design.strategy, design.strategy)),
        ("Drawn", MEMBERS_SHOWN[context["members_shown"]]),
        ("Palette", dict(context["palettes"])[context["palette"]]),
        ("Exported", timezone.now().strftime("%Y-%m-%d %H:%M UTC")),
    ]
    frames = [frame for frame, _ in context["member_frames"]]
    response = HttpResponse(
        sites_excel.library_workbook(frames, context["members"], context["palette"], about),
        content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
    name = re.sub(r"[^A-Za-z0-9._+-]", "_", f"{record.code}.members")
    response["Content-Disposition"] = f'attachment; filename="{name}.xlsx"'
    return response


def library_sequence(request: HttpRequest, job_id: str, library_id: str) -> HttpResponse:
    """One target's site map, fetched when its row is opened (the page would otherwise
    carry every residue of every ancestor). The target is named by its key, never its
    position, so a page that is out of date cannot fetch another ancestor."""
    job = owned_job(request, job_id)
    record = owned_library(request, job, library_id)
    key = tuple(request.GET.get(k, "") for k in ("combination", "node", "version"))
    target = next((t for t in record.targets if library.target_key(t) == key), None)
    if target is None:
        raise Http404
    (data,) = load_targets(job, [target])
    if data is None:
        raise Http404
    threshold = library.parse_threshold(request.GET.get("altall", ""))
    context = {
        "m": library_glyphs.site_map(data.template, data.dists, threshold),
        "threshold_pct": max(library_glyphs.threshold_percent(threshold), library_glyphs.FLOOR),
    }
    return render(request, "phylobot/_library_site_map.html", context)


DELETED_SHOWN = 100  # deleted libraries listed beneath the live ones (the admin has all)


def my_libraries(request: HttpRequest) -> HttpResponse:
    """Every library of the user's, the last changed first, one card each. The library is
    the unit here, not the job (owner decision 2026-10-07): a card lists the job or jobs it
    draws on as a sub-item. Libraries on deleted jobs are left out: their results, and so
    their ancestors, are gone.

    Nothing on this page generates members (owner decision 2026-10-07): every number comes
    from the spec, the job's params or an export receipt, so thirty libraries cost one
    query and no storage reads.
    """
    from phylobot.views import is_uuid  # views imports this module

    records = (
        Library.objects.for_user(request.user)
        .filter(job__deleted_at__isnull=True)
        .select_related("job")
        .prefetch_related("exports", "job__tags")
    )
    state = request.GET.get("state", "")
    if state not in library_index.STATES:
        state = ""
    tag = request.GET.get("tag", "").lstrip("#").lower()
    if tag:  # libraries have no tags of their own: they filter by their job's
        records = records.filter(job__tags__name=tag)
    job_id = request.GET.get("job", "")
    if is_uuid(job_id):
        records = records.filter(job_id=job_id)

    cards = [library_index.card(record, list(record.exports.all())) for record in records]
    if state:
        cards = [c for c in cards if c.state == state]

    deleted = (
        Library.objects.deleted_for(request.user)
        .select_related("job")
        .order_by("-deleted_at")[:DELETED_SHOWN]
    )
    context = {
        "cards": cards,
        "total": len(cards),
        "deleted": deleted,
        "state": state,
        "state_choices": [(key, library_index.STATE_LABELS[key]) for key in library_index.STATES],
        "tag": tag,
        "job_id": job_id if is_uuid(job_id) else "",
        "all_tags": Tag.objects.filter(owner_id=request.user.pk).values_list("name", flat=True),
        "filtered": bool(state or tag or is_uuid(job_id)),
    }
    return render(request, "phylobot/my_libraries.html", context)


@require_POST
def new(request: HttpRequest, job_id: str) -> HttpResponse:
    """Start another library on this job; it becomes the current one."""
    job = owned_job(request, job_id)
    if job.status != Job.Status.SUCCEEDED:
        raise Http404
    with transaction.atomic():
        Job.objects.select_for_update().filter(id=job.id).first()
        record = create_library(request.user, job)
    return redirect("library_page", job_id=job.id, library_id=record.id)


@require_POST
def rename(request: HttpRequest, job_id: str, library_id: str) -> HttpResponse:
    """Change a library's name, code or notes (each only when posted: the title's pencil
    posts the name alone). The code is cleaned, never refused: member names must stay
    well formed (`library.short_code`)."""
    job = owned_job(request, job_id)
    record = owned_library(request, job, library_id)
    name = " ".join(request.POST.get("name", "").split())[:100]
    if name:
        record.name = name
    if "code" in request.POST:
        record.code = library.short_code(request.POST.get("code", ""))
    if "notes" in request.POST:
        record.notes = request.POST["notes"].replace("\r\n", "\n")[: Library.NOTES_LIMIT]
    record.save()
    return redirect("library_page", job_id=job.id, library_id=record.id)


@require_POST
def remove(request: HttpRequest, job_id: str, library_id: str) -> HttpResponse:
    """Remove one target, named by its key (combination, node, version)."""
    job = owned_job(request, job_id)
    key = tuple(request.POST.get(k, "") for k in ("combination", "node", "version"))
    with transaction.atomic():
        Job.objects.select_for_update().filter(id=job.id).first()
        record = owned_library(request, job, library_id)
        record.spec = library.remove_target(record.spec, key)  # type: ignore[arg-type]
        record.save()
    return redirect(back_to(request, job))


ACTIONS = ("remove", "copy", "move")


def posted_keys(request: HttpRequest) -> set[tuple[str, str, str]]:
    """The ticked targets: each checkbox's value is its key as JSON, never a position."""
    keys: set[tuple[str, str, str]] = set()
    for raw in request.POST.getlist("target"):
        try:
            value = json.loads(raw)
        except ValueError:
            raise Http404 from None
        if not (
            isinstance(value, list) and len(value) == 3 and all(isinstance(v, str) for v in value)
        ):
            raise Http404
        keys.add((value[0], value[1], value[2]))
    return keys


@require_POST
def targets(request: HttpRequest, job_id: str, library_id: str) -> HttpResponse:
    """Remove the ticked targets, or copy or move them to other libraries of this job
    (POST action, target..., to: library ids or "new"; several only for a copy). Targets
    the library no longer holds (a page out of date) are skipped; the rest are acted on."""
    job = owned_job(request, job_id)
    action = request.POST.get("action", "")
    if action not in ACTIONS:
        raise Http404
    wanted = posted_keys(request)
    if not wanted:
        messages.info(request, "Tick the ancestors to act on first.")
        return redirect(back_to(request, job))
    destinations = [to for to in request.POST.getlist("to") if to]
    if action == "copy" and not destinations:
        raise Http404
    if action == "move" and len(destinations) != 1:
        raise Http404  # moving to several libraries is not one action
    with transaction.atomic():
        Job.objects.select_for_update().filter(id=job.id).first()
        record = owned_library(request, job, library_id)
        chosen = [t for t in record.targets if library.target_key(t) in wanted]
        copies: list[tuple[str, int]] = []  # (library name, how many were new)
        for to in dict.fromkeys(destinations):
            dest = destination(request, job, to)
            if dest.id == record.id:
                messages.info(request, "Choose another library to copy or move to.")
                return redirect(back_to(request, job))
            dest.spec, added = library.add_targets(dest.spec, chosen)
            dest.save()
            copies.append((dest.name, added))
        if action in ("remove", "move"):
            for target in chosen:
                record.spec = library.remove_target(record.spec, library.target_key(target))
            record.save()
    count = len(chosen)
    what = f"{count} ancestor{'s' if count != 1 else ''}"
    if action == "remove":
        messages.success(request, f"Removed {what} from {record.name}.")
    else:
        verb = "Copied" if action == "copy" else "Moved"
        where = ", ".join(
            f"{name} ({count - added} already there)" if count - added else name
            for name, added in copies
        )
        messages.success(request, f"{verb} {what} to {where}.")
    return redirect(back_to(request, job))


@require_POST
def delete(request: HttpRequest, job_id: str, library_id: str) -> HttpResponse:
    """Delete a library: it leaves the user's lists, and the record stays (as jobs)."""
    job = owned_job(request, job_id)
    record = owned_library(request, job, library_id)
    record.deleted_at = timezone.now()
    record.save(update_fields=["deleted_at"])
    messages.info(request, f"Deleted {record.name}.")
    return redirect("job_detail", job_id=job.id)


def context(request: HttpRequest) -> dict[str, Any]:
    """The top-bar chip on a job's pages: the library shown on its own page, else the
    job's current one, which is also where the viewers' "Add to library" buttons add.
    One query, and only on URLs that name a job."""
    match = getattr(request, "resolver_match", None)
    user = getattr(request, "user", None)
    if match is None or user is None or not user.is_authenticated:
        return {}
    job_id = match.kwargs.get("job_id")
    if job_id is None:
        return {}
    found = Library.objects.for_user(user).filter(job_id=job_id, job__deleted_at__isnull=True)
    library_id = match.kwargs.get("library_id")
    record = (found.filter(id=library_id).first() if library_id else None) or found.first()
    if record is None:
        return {"library_choices": []}  # the viewers' buttons then say "a new library"
    # The viewers' "Add to" menus: the job's libraries, the current one first.
    return {"library_chip": record, "library_choices": list(found)}
