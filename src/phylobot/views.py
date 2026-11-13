import csv
import io
import json
import logging
import re
import tempfile
import uuid
from collections import Counter
from collections.abc import Collection
from dataclasses import dataclass, replace
from pathlib import Path, PurePosixPath
from typing import Any
from urllib.parse import urlencode

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import login_not_required
from django.contrib.auth.models import User
from django.db import connection, transaction
from django.db.models import Count, Sum
from django.http import (
    Http404,
    HttpRequest,
    HttpResponse,
    JsonResponse,
    QueryDict,
    StreamingHttpResponse,
)
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.utils import timezone
from django.views.decorators.http import require_GET, require_http_methods, require_POST

from phylobot import (
    about as about_page,
)
from phylobot import (
    alignment_views,
    alignments,
    analog_views,
    analogs,
    ancestor_names,
    ancestors,
    branch_site,
    branch_site_views,
    bundle,
    calibration,
    catalog,
    clade_collapse,
    codeml,
    codons,
    diagnosis,
    esr_views,
    eta,
    example_sets,
    extend,
    filelist,
    gaps,
    homolog_views,
    indel_views,
    input_views,
    inputs,
    job_import,
    job_types,
    jobs,
    joined,
    library,
    library_views,
    machine,
    manifest,
    paths,
    properties,
    properties_run,
    properties_submit,
    property_compare,
    property_views,
    reconcile,
    results,
    saved_paths,
    searches,
    simulation_views,
    site_settings,
    species_web,
    stages,
    storage,
    trees,
    treeview,
    uncertainty_submit,
    work_records,
    worklist,
    workload,
)
from phylobot.executors import alignment_count, combination_count, log_name
from phylobot.forms import (
    DEFAULT_MODELS,
    RECONCILED_TREE,
    CalibrationForm,
    SiteSettingsForm,
    SubmitForm,
    analysis_flags,
    check_input,
    codon_grid,
    codon_tree_params,
    describe,
    grid,
    initial_from_job,
    section_defaults,
    sequence_names,
    ungapped,
    used_files,
)
from phylobot.models import (
    Alignment,
    AncestorName,
    HomologSearch,
    Job,
    SiteSettings,
    Structure,
    Tag,
    UserCalibration,
    limits_for,
    parse_tags,
)
from phylobot.pipeline.validate import detect_seqtype
from phylobot.species import resolve as species_resolve
from phylobot.species.table import OK as SPECIES_OK
from phylobot.species.table import Row as SpeciesRow
from phylobot.templatetags.phylobot_extras import duration

logger = logging.getLogger(__name__)

MAX_OUTGROUP_LIST = 2000  # sequences listed in the outgroup checklist
# Shown inline as plain text; everything else (notably HTML reports, which contain
# script and user-supplied names) is download-only.
TEXT_SUFFIXES = {".fasta", ".tsv", ".nwk", ".txt", ".json", ".log", ".err", ".csv"}
LOG_TAIL_LINES = 40
LIVE_LOG = "logs/pipeline.log"  # under results/; uploaded while the job runs


def current_user(request: HttpRequest) -> User:
    """The signed-in user (local mode: the local user); middleware guarantees one."""
    assert isinstance(request.user, User)
    return request.user


def about(request: HttpRequest) -> HttpResponse:
    """What PhyloBot is, the 2016 paper, and every method it uses, with citations."""
    return render(
        request,
        "phylobot/about.html",
        {
            "reference": about_page.PHYLOBOT_2016,
            "methods": about_page.METHODS,
            "matrices": about_page.MATRICES,
            "variants": about_page.VARIANTS,
            "raxml_missing": about_page.RAXML_MISSING,
            "comparison_tools": about_page.COMPARISON_TOOLS,
            "comparison": about_page.COMPARISON,
            "comparison_legend": about_page.COMPARISON_LEGEND,
        },
    )


def copied_job(request: HttpRequest, job_id: str) -> Job | None:
    """The user's own job that "Duplicate and run again" copies (None if not given)."""
    if not job_id:
        return None
    if not is_uuid(job_id):
        raise Http404
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if job.params.get("simulate") or job.params.get("esr"):  # not an analysis to copy
        raise Http404
    return job


def plan_after_alignments(extension: extend.Inheritance | None) -> list[str]:
    """What a parent with trees did after its alignments, which an extension's new
    alignments repeat (the locked steps 4 and 5)."""
    if extension is None or extension.mode != extend.RESULTS:
        return []
    params = {**extension.parent.params, "aligners": [], "left_out_count": 0, "extends": None}
    steps = describe(params, 0, len(extension.trees))
    return [step for step in steps if not step.startswith("check and translate")]


def extended_job(request: HttpRequest, job_id: str) -> Job | None:
    """The user's own finished job that "Extend" builds on (None if not given)."""
    if not job_id:
        return None
    if not is_uuid(job_id):
        raise Http404
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if extend.refusal(job):
        raise Http404
    return job


def copied_inputs(job: Job) -> tuple[str, dict[str, list[tuple[str, bytes]]]]:
    """A job's stored inputs: its sequences (text) and uploaded alignments and trees
    ((label, content), in order, labelled as when submitted)."""
    prefix = f"{job.storage_prefix}inputs/"
    keys = set(storage.list_keys(prefix))
    fasta = ""
    if f"{prefix}sequences.fasta" in keys:
        fasta = storage.read_bytes(f"{prefix}sequences.fasta").decode(errors="replace")
    reuse: dict[str, list[tuple[str, bytes]]] = {}
    for kind, ext in (("alignment", "fasta"), ("tree", "nwk")):
        labels = job.params.get(f"{kind}_labels") or []
        for index, label in enumerate(labels, 1):
            key = f"{prefix}{kind}{index}.{ext}"
            if key in keys:
                reuse.setdefault(kind, []).append((label, storage.read_bytes(key)))
    return fasta, reuse


def submit(request: HttpRequest) -> HttpResponse:
    limits = limits_for(request.user)
    from_homologs: dict[str, Any] | None = None  # step 2 filled from a homolog review
    from_simulation: dict[str, Any] | None = None  # step 2 filled from a simulation
    copy_selected: dict[str, set[str]] | None = None  # model grids of a copied job
    source = copied_job(request, request.POST.get("from_job") or request.GET.get("from") or "")
    extending = request.POST.get("extend_job") or request.GET.get("extend") or ""
    parent = extended_job(request, extending)
    extension = None
    if parent is not None:
        try:
            extension = extend.inherit(parent)
        except LookupError:
            raise Http404 from None
    reuse: dict[str, list[tuple[str, bytes]]] = {}
    owner = request.user
    if request.method == "POST":
        if source is not None:
            _, reuse = copied_inputs(source)
        form = SubmitForm(
            request.POST, request.FILES, limits=limits, reuse=reuse, owner=owner,
            extension=extension,
        )  # fmt: skip
        if form.is_valid():
            with transaction.atomic():
                params = form.params()
                job = Job.objects.create(
                    owner=current_user(request),
                    name=form.cleaned_data["name"],
                    params=params,
                    work_units=workload.estimate(params, form.size()).units,
                    input_size=workload.size_record(form.size()),
                    rerun_of=source,
                    extends=parent,
                )
                if (origin := source or parent) is not None:
                    job.tags.set(origin.tags.all())
                for filename, content in form.input_files().items():
                    storage.put_bytes(f"{job.storage_prefix}inputs/{filename}", content)
                alignments.touch(*(form.saved[k] for k in form.cleaned_data["saved_alignments"]))
                transaction.on_commit(lambda: jobs.enqueue(job))
            return redirect("job_detail", job_id=job.id)
    elif extension is not None:
        assert parent is not None
        initial = {
            "name": parent.name,
            "sequences_text": extension.fasta.decode(errors="replace"),
            "seqtype": parent.params.get("seqtype") or "auto",
            "genetic_code": parent.params.get("genetic_code") or 1,
            "aligners": extension.aligners,
            "outgroup": parent.params.get("outgroup") or [],
            "seed": parent.params.get("seed", 1),
            "extend_job": parent.id,
        }
        form = SubmitForm(initial=initial, limits=limits, owner=owner, extension=extension)
    elif source is not None:
        fasta, reuse = copied_inputs(source)
        initial, copy_selected = initial_from_job(source.params, source.name, fasta)
        form = SubmitForm(initial=initial | {"from_job": source.id}, limits=limits, owner=owner)
    else:
        from_homologs = homolog_views.job_prefill(request)
        initial = {}
        if from_homologs:
            initial = {
                "sequences_text": from_homologs["fasta"],
                "species_hints": from_homologs["taxids"],
            }
        # "Use in a new job" on a simulation's replicate: its unaligned sequences.
        from_simulation = simulation_views.job_prefill(request)
        if from_simulation:
            initial = {
                "sequences_text": from_simulation["fasta"],
                "simulated_from": f"{from_simulation['job'].id}:{from_simulation['seed']}",
            }
        # "Use in a new job" on My alignments: that alignment chosen (form: the user's only).
        if chosen := request.GET.getlist("alignment"):
            initial["saved_alignments"] = chosen
        form = SubmitForm(initial=initial, limits=limits, owner=owner)
    # Unbound: every candidate (IQ-TREE's default ModelFinder set), the default full searches.
    if form.is_bound:
        grids = ("candidates", "models", "codon_candidates", "codon_models")
        selected = {g: set(request.POST.getlist(g)) for g in grids}
    elif copy_selected is not None:
        selected = copy_selected
    else:
        selected = {
            "candidates": set(catalog.all_models()),
            "models": set(DEFAULT_MODELS),
            "codon_candidates": set(catalog.all_codon_models()),
            "codon_models": set(),
        }
    names = sequence_names(submitted_fasta(request)) if form.is_bound else []
    if not form.is_bound and from_homologs:
        names = sequence_names(from_homologs["fasta"].encode())
    if not form.is_bound and from_simulation:
        names = sequence_names(from_simulation["fasta"].encode())
    outgroup_checked = set(request.POST.getlist("outgroup"))
    if not form.is_bound and source is not None:
        fasta_text = str(form.initial.get("sequences_text") or "")
        first = reuse.get("alignment", [("", b"")])[0][1]
        names = sequence_names(fasta_text.encode() if fasta_text else ungapped(first))
        outgroup_checked = set(source.params.get("outgroup") or [])
    if not form.is_bound and extension is not None:
        first = extension.fasta or ungapped(extension.alignments[0].content)
        names = sequence_names(first)
        outgroup_checked = set(extension.parent.params.get("outgroup") or [])
    chosen_alignments = list(form.saved.values())  # the user's own, as posted or linked
    if not form.is_bound and chosen_alignments and not names:
        loaded = chosen_alignments[0]
        names = sequence_names(ungapped(storage.read_bytes(loaded.storage_key)))
    context = {
        "outgroup_names": names[:MAX_OUTGROUP_LIST],
        "outgroup_width": name_width(names[:MAX_OUTGROUP_LIST]),
        "outgroup_checked": outgroup_checked,
        "form": form,
        "copied": source,
        "extension": extension,
        "extension_plan": plan_after_alignments(extension),
        "extension_sources": [("iqtree", "IQ-TREE 3"), ("raxml", "RAxML-NG 2")],
        "extension_models": (
            set(extension.parent.params.get("models") or [])
            if extension is not None and extension.parent.params.get("fixed_models")
            else set()
        ),
        # The parent's sequences, for the page's checks and summary (SubmitForm ignores them).
        "extension_fasta": (
            (extension.fasta or ungapped(extension.alignments[0].content)).decode(errors="replace")
            if extension is not None
            else ""
        ),
        "copied_files": {kind: [label for label, _ in files] for kind, files in reuse.items()},
        "grid": grid(),
        "variants": catalog.VARIANTS,
        "fixed_grid": grid(mixtures=True),
        "fixed_variants": catalog.VARIANTS + catalog.MIXTURE_VARIANTS,
        "codon_grid": codon_grid(),
        "codon_variants": catalog.CODON_VARIANTS,
        "paml_codes": ",".join(
            str(c.ncbi) for c in codons.GENETIC_CODES if c.paml_icode is not None
        ),
        "raxml_missing": catalog.RAXML_MISSING,
        "section_defaults": section_defaults(),
        "selected": selected,
        "limits": limits,
        "from_homologs": from_homologs,
        "from_simulation": from_simulation,
        "chosen_alignments": chosen_alignments,
        "examples": example_sets.EXAMPLES,
    }
    if form.is_bound:  # redrawn after an error: the species table as the user left it
        context |= species_context(request)
    return render(request, "phylobot/submit.html", context)


def can_edit_site_settings(request: HttpRequest) -> bool:
    """Local mode: the one (local) user. Cloud mode: staff only."""
    if settings.PHYLOBOT_MODE == "local":
        return True
    return bool(request.user.is_authenticated and request.user.is_staff)


WORK_ESTIMATE_JOBS = 25  # the viewer's recent jobs compared with their estimates
OUTLIER = 2.0  # flag jobs that used more than 2x, or under half, of their estimate


@require_GET
def work_estimate(request: HttpRequest) -> HttpResponse:
    """How the New Job page's work estimate is made: this server's calibration, example
    jobs, the reference costs, and the viewer's own recent jobs, estimated beside used.
    Open to every user in both modes; the job table never shows other users' jobs."""
    spu, learned_from = workload.calibration()
    finished = (
        Job.objects.for_user(request.user)
        .filter(status=Job.Status.SUCCEEDED, work_units__gt=0, cpu_seconds__gt=0)
        .order_by("-finished_at")[:WORK_ESTIMATE_JOBS]
    )
    rows = []
    for job in finished:
        assert job.work_units and job.cpu_seconds  # filtered above
        expected = job.work_units * (spu or workload.SECONDS_PER_UNIT)
        ratio = job.cpu_seconds / expected
        rows.append({
            "job": job,
            "expected": expected,
            "ratio": ratio,
            "outlier": ratio > OUTLIER or ratio < 1 / OUTLIER,
        })  # fmt: skip
    examples = []
    for size in workload.EXAMPLE_SIZES:
        estimate = workload.estimate(workload.EXAMPLE_PARAMS, size, spu)
        examples.append({
            "size": size,
            "units": estimate.units,
            "seconds": estimate.units * (spu or workload.SECONDS_PER_UNIT),
            "largest": max(estimate.groups, key=lambda g: g[2])[1],
        })  # fmt: skip
    context = {
        "seconds_per_unit": spu,
        "learned_from": learned_from,
        "min_history": workload.MIN_HISTORY,
        "history": workload.HISTORY,
        "reference": workload.REFERENCE_SIZE,
        "costs": workload.reference_costs(),
        "examples": examples,
        "rows": rows,
        "outlier": OUTLIER,
    }
    return render(request, "phylobot/work_estimate.html", context)


@require_GET
def work_records_page(request: HttpRequest) -> HttpResponse:
    """Every user's finished jobs, stage by stage: the estimate beside the CPU time used
    (phylobot.work_records). The admin's: others get a 404. ?limit= jobs (newest first);
    ?format=csv for one line per job and stage."""
    if not can_edit_site_settings(request):
        raise Http404
    try:
        limit = min(max(int(request.GET.get("limit", "")), 1), work_records.MAX_LIMIT)
    except ValueError:
        limit = work_records.DEFAULT_LIMIT
    rows = work_records.records(limit)
    if request.GET.get("format") == "csv":
        response = HttpResponse(work_records.to_csv(rows), content_type="text/csv")
        response["Content-Disposition"] = 'attachment; filename="work-records.csv"'
        return response
    spu, learned_from = workload.calibration()
    per_unit = spu or workload.SECONDS_PER_UNIT
    return render(
        request,
        "phylobot/work_records.html",
        {
            "rows": [
                {
                    "record": r,
                    "ratio": r.ratio(per_unit),
                    "stages": [(s, s.ratio(per_unit)) for s in r.stages],
                }
                for r in rows
            ],
            "summary": work_records.summary(rows, per_unit),
            "seconds_per_unit": spu,
            "per_unit": per_unit,
            "learned_from": learned_from,
            "history": workload.HISTORY,
            "limit": limit,
            "local_mode": settings.PHYLOBOT_MODE == "local",
        },
    )


def initial(setting: site_settings.Setting, stored: Any) -> Any:
    """A stored value as the Settings form shows it (an on/off switch posts "1"/"0")."""
    if setting.switch and stored is not None:
        return "1" if stored else "0"
    return stored


@require_http_methods(["GET", "POST"])
def site_settings_page(request: HttpRequest) -> HttpResponse:
    """Site-wide settings (phylobot.site_settings), overriding .env without a restart.
    Others get a 404, as for any page they may not see."""
    if not can_edit_site_settings(request):
        raise Http404
    row = site_settings.stored()
    if request.method == "POST":
        form = SiteSettingsForm(request.POST)
        if form.is_valid():
            row = row or SiteSettings(pk=1)
            for s in site_settings.SETTINGS:
                chosen = form.cleaned_data[s.field]
                if s.field == "search_email":
                    chosen = chosen or ""  # a text field: blank, not NULL
                setattr(row, s.field, chosen)
            row.updated_by = current_user(request)
            row.save()
            logger.info("site settings changed by %s", request.user)
            messages.success(request, "Settings saved. They apply from now on; no restart needed.")
            return redirect("site_settings")
    else:
        form = SiteSettingsForm(
            initial={
                s.field: initial(s, getattr(row, s.field, None)) for s in site_settings.SETTINGS
            }
        )
    failed = Job.objects.filter(
        status=Job.Status.FAILED, deleted_at__isnull=True, work_expired_at__isnull=True
    )
    kept = failed.filter(work_bytes__isnull=False).aggregate(
        jobs=Count("id"), size=Sum("work_bytes")
    )
    archived = failed.filter(work_archive_attempt__isnull=False).aggregate(
        jobs=Count("id"), size=Sum("work_archive_bytes")
    )
    local_mode = settings.PHYLOBOT_MODE == "local"
    # Local mode has no accounts, so the per-user framing is misleading: Billing is hidden
    # and a few labels drop "per user" / "for every user" (they are just this machine's
    # limits across jobs). The group key stays the same (stable section ids); only the
    # heading and these labels change.
    local_group_titles = {site_settings.USERS: "Limits Across Jobs"}
    local_labels = {
        "searches_per_day": "Searches per day",
        "searches_in_flight_per_user": "Searches waiting at a service",
    }
    grouped: dict[str, list[dict[str, Any]]] = {}
    for s in site_settings.SETTINGS:
        if local_mode and s.group == site_settings.BILLING:
            continue
        grouped.setdefault(s.group, []).append(
            {
                "field": form[s.field],
                "setting": s,
                "label": local_labels.get(s.field, s.label) if local_mode else s.label,
                "value": site_settings.value(s.field, row),
                "source": site_settings.source(s.field, row),
                "kept": kept if s.field == "keep_failed_work_days" else None,
                "archived": archived if s.field == "archive_failed_work" else None,
            }
        )
    sections = [
        {
            "key": group,
            "title": local_group_titles.get(group, group) if local_mode else group,
            "rows": rows,
        }
        for group, rows in grouped.items()
    ]
    return render(
        request,
        "phylobot/site_settings.html",
        {
            "form": form,
            "sections": sections,
            "row": row,
            "local_mode": local_mode,
            "calibration_group": site_settings.CALIBRATION,
            "jobs_group": site_settings.JOBS,
            "machine": machine_panel(),
        },
    )


def machine_panel() -> dict[str, Any]:
    """The Settings page's "This machine": what the worker found, the totals in force
    (detected, then .env, then this page), and the share each job gets at launch."""
    found = machine.latest()
    detected = (
        machine.Detected(found.cpus, found.isolated, found.memory_bytes, found.note)
        if found is not None
        else None
    )
    totals = machine.totals(detected)
    slots = machine.slots()
    share = machine.share(totals, slots) if totals.cpus and totals.memory_bytes else None
    return {
        "found": found,
        "totals": totals,
        "slots": slots,
        "share": share,
        "task_memory": totals.task_bytes,
        "head_memory": machine.HEAD_MEMORY,
    }


# Settings and My calibration offer the same two resets: back to the layer below, or every
# term 0 (ln T = 0, so T = 1 at every site: calibrated equals raw).
CALIBRATION_RESETS = {"default", "none"}


@require_POST
def calibration_defaults(request: HttpRequest) -> HttpResponse:
    """The Settings page's resets of the site's calibration defaults: to the built-in
    coefficients (blank), or to no effect (0). Touches no other setting."""
    if not can_edit_site_settings(request):
        raise Http404
    to = request.POST.get("to", "")
    if to not in CALIBRATION_RESETS:
        raise Http404
    row = site_settings.stored() or SiteSettings(pk=1)
    for term in site_settings.CALIBRATION_TERMS:
        setattr(row, f"calibration_{term}", None if to == "default" else 0.0)
    row.updated_by = current_user(request)
    row.save()
    messages.success(
        request,
        "Calibration defaults reset to the built-in coefficients."
        if to == "default"
        else "Calibration defaults set to no effect: T = 1, calibrated equals raw.",
    )
    return redirect(reverse("site_settings") + "#calibrated-posterior-probabilities")


@require_http_methods(["GET", "POST"])
def my_calibration(request: HttpRequest) -> HttpResponse:
    """Every user's own calibration coefficients, over the site's defaults: a blank field
    uses the site's value. Save, reset to the site's defaults, or set to no effect."""
    user = current_user(request)
    own = site_settings.own_calibration(user)
    terms = site_settings.CALIBRATION_TERMS
    if request.method == "POST":
        action = request.POST.get("action", "save")
        if action in CALIBRATION_RESETS:
            own = own or UserCalibration(user=user)
            for term in terms:
                setattr(own, term, None if action == "default" else 0.0)
            own.save()
            messages.success(
                request,
                "Your coefficients now follow this site's defaults."
                if action == "default"
                else "Your coefficients are set to no effect: T = 1, calibrated equals raw.",
            )
            return redirect("my_calibration")
        form = CalibrationForm(request.POST)
        if form.is_valid():
            own = own or UserCalibration(user=user)
            for term in terms:
                setattr(own, term, form.cleaned_data[term])
            own.save()
            messages.success(request, "Your coefficients are saved. Viewers use them at once.")
            return redirect("my_calibration")
    else:
        form = CalibrationForm(initial={t: getattr(own, t, None) for t in terms})
    site = site_settings.site_calibration()
    in_effect = site_settings.calibration_coefficients(user)
    covariates = {c.key: c for c in calibration.COVARIATES}
    rows = []
    for term in terms:
        rows.append({
            "field": form[term],
            "label": "Intercept (ln T at the reference site)" if term == "intercept"
            else covariates[term].label,
            "meaning": "" if term == "intercept"
            else f"{covariates[term].meaning} 0 at {covariates[term].reference}.",
            "built_in": site_settings.built_in_term(term),
            "site": site[term],
            "in_effect": in_effect.intercept if term == "intercept" else in_effect.betas[term],
        })  # fmt: skip
    return render(
        request,
        "phylobot/my_calibration.html",
        {
            "rows": rows,
            "form": form,
            "own": own,
            "coefficients": in_effect,
            "can_edit_site": can_edit_site_settings(request),
            "t_bounds": (calibration.T_MIN, calibration.T_MAX),
        },
    )


@require_GET
def example(request: HttpRequest, name: str) -> HttpResponse:
    entry = example_sets.BY_KEY.get(name)
    if entry is None:
        raise Http404
    # ?variant=subset loads the quick-start subset; anything else loads the full set.
    variant = "subset" if request.GET.get("variant") == "subset" else "full"
    from importlib.resources import files

    text = files("phylobot").joinpath(example_sets.resource(name, variant)).read_text()
    form = SubmitForm(initial={"sequences_text": text}, limits=limits_for(request.user))
    context = {
        "field": form["sequences_text"],
        # the outgroup checklist and the sequence check are refreshed too (out-of-band
        # swaps); the example's outgroup is ticked
        "outgroup_names": sequence_names(text.encode()),
        "outgroup_width": name_width(sequence_names(text.encode())),
        "outgroup_checked": set(entry.outgroup),
        "check": check_input(text.encode(), "auto", 1, limits_for(request.user)),
        "oob": True,
    }
    return render(request, "phylobot/_sequences_field.html", context)


def submitted_fasta(request: HttpRequest) -> bytes:
    """The FASTA in a submit-form POST: the pasted text, the ticked sequence files
    combined, or else the first ticked alignment with its gaps removed (all bounded)."""
    limit = limits_for(request.user).max_upload_bytes
    merged = input_views.merged_sequences(request)
    text = request.POST.get("sequences_text", "")
    if merged is not None:  # a problem: the pasted text as is, so the check names it
        return text.encode() if merged.error and text.strip() else merged.fasta
    uploads = used_files(request.FILES.getlist("alignment_files"), request.POST, "alignment")
    if not text.strip() and uploads:
        first = inputs.dedupe_alignment(
            uploads[0].name or "", input_views.read_upload(uploads[0], limit),
            request.POST.get("alignment_duplicates", ""),
        )  # fmt: skip
        lines = first.fasta.decode(errors="replace").splitlines()
        text = "\n".join(ln if ln.startswith(">") else re.sub(r"[-.]", "", ln) for ln in lines)
    chosen = [key for key in request.POST.getlist("saved_alignments") if is_uuid(key)]
    if not text.strip() and chosen:  # the first alignment chosen from My alignments
        record = Alignment.objects.for_user(request.user).filter(id=chosen[0]).first()
        if record is not None:
            return ungapped(storage.read_bytes(record.storage_key)[:limit])
    return text.encode()


def name_width(names: list[str]) -> int:
    """The outgroup checklist's column width, in characters: its longest name, so short
    names pack into more columns (and long ones never overlap)."""
    return min(max((len(n) for n in names), default=12), 60)


def species_context(request: HttpRequest) -> dict[str, Any]:
    """The reconciliation block's species table for the sequences in a submit-form POST
    (species_web): from the headers, accessions (with consent), lists and the user's
    cells, resolved on an uploaded species tree or at the Open Tree of Life (with its
    own consent). Rows that need the user come first; the consent boxes live in the
    table, shown only when they apply."""
    if request.POST.get("reconcile") != "on":
        return {"species": None}
    limit = limits_for(request.user).max_upload_bytes

    def uploaded(field: str) -> str:
        upload = request.FILES.get(field)
        if upload is None:
            return ""
        upload.seek(0)  # the form may have read it already
        return str(upload.read(limit).decode(errors="replace"))

    text = uploaded("species_tree_file")
    opentree = request.POST.get("species_consent") == "on"
    gathered = species_web.gather(
        submitted_fasta(request), request.POST, uploaded("species_table_file")
    )
    found = species_web.look_up(gathered.rows, opentree, text)
    drawing = None
    if found.tree is not None:
        drawing = treeview.layout(found.tree.newick(), "cladogram")
        # Newick writes spaces as underscores; species read better with spaces.
        drawing.leaves = [
            replace(leaf, name=leaf.name.replace("_", " ")) for leaf in drawing.leaves
        ]
    resolved = found.state not in (species_web.NEEDS_CONSENT, species_web.NO_SEQUENCES)
    # a protein shared by several species, until the user picks one of them
    choices = {
        row.sequence: answer.choices
        for row in found.rows
        if (answer := gathered.accessions.found.get(row.sequence)) is not None
        and answer.choices
        and row.query == answer.organism
    }

    def needs(row: SpeciesRow) -> bool:
        if row.sequence in choices:
            return True
        if not resolved:
            return not row.query and row.taxid is None
        return row.status != SPECIES_OK

    rows = [
        {
            "row": row,
            "shown": species_web.shown(row),
            "source": species_web.SOURCES.get(row.source, row.source),
            "needs": needs(row),
            "guess": row.source == "guess",
            "choices": choices.get(row.sequence, ()),
            # a taxon the synthetic tree leaves out can stand in at its parent
            "parent": found.parents.get(row.sequence),
        }
        for row in found.rows
    ]
    rows.sort(key=lambda item: not item["needs"])  # stable: the rest keep their order
    pending = sum(1 for row in found.rows if needs(row))
    return {
        "species": found,
        "species_rows": rows,
        "species_pending": pending,
        "species_named": sum(1 for r in found.rows if r.query or r.taxid is not None),
        "species_guesses": sum(1 for row in found.rows if row.source == "guess"),
        "species_offers": species_web.offers(found.rows, resolved),
        "species_accessions": gathered.accessions,
        "species_consents": {
            "opentree": opentree,
            "accessions": request.POST.get("accession_consent") == "on",
        },
        "species_drawing": drawing,
        "species_upload": bool(text.strip()),
        "species_citations": species_resolve.CITATIONS,
    }


@require_POST
def species_table(request: HttpRequest) -> HttpResponse:
    """The species table as it stands, as a TSV download: edit it in a spreadsheet and
    upload it back as a species list (species_web.bulk reads the header)."""
    context = species_context(request)
    found = context["species"]
    lines = ["sequence\tspecies\ttaxid\tin_species_tree\tstatus\tnote"]
    for row in found.rows if found is not None else []:
        cells = [row.sequence, row.query, str(row.taxid or ""), row.species, row.status, row.note]
        lines.append("\t".join(" ".join(c.split()) for c in cells))
    response = HttpResponse("\n".join(lines) + "\n", content_type="text/tab-separated-values")
    response["Content-Disposition"] = 'attachment; filename="species.tsv"'
    return response


@require_POST
def species_check(request: HttpRequest) -> HttpResponse:
    """The New Job page's species table (HTMX)."""
    return render(request, "phylobot/_species.html", species_context(request))


@require_POST
def outgroup_options(request: HttpRequest) -> HttpResponse:
    """The outgroup checklist for the sequences currently in the form (HTMX)."""
    names = sequence_names(submitted_fasta(request))
    context = {
        "outgroup_names": names[:MAX_OUTGROUP_LIST],
        "outgroup_width": name_width(names[:MAX_OUTGROUP_LIST]),
        "outgroup_checked": set(request.POST.getlist("outgroup")) & set(names),
    }
    return render(request, "phylobot/_outgroup.html", context)


def genetic_code(request: HttpRequest) -> int:
    """The form's genetic code (1 if missing; VALIDATE refuses codes not offered)."""
    value = request.POST.get("genetic_code", "1")
    return int(value) if value.isdigit() and len(value) <= 2 else 1


def seqtype(request: HttpRequest, fasta: bytes) -> str:
    """The sequence type in effect: chosen, or detected from the sequences."""
    chosen = request.POST.get("seqtype", "auto")
    if chosen in ("protein", "codon"):
        return chosen
    return detect_seqtype(fasta.decode("utf-8", errors="replace"))


@require_POST
def sequence_check(request: HttpRequest) -> HttpResponse:
    """The sequences' type and size, or what the pipeline would refuse them for (HTMX)."""
    fasta = submitted_fasta(request)
    check = None
    if fasta.strip():
        check = check_input(
            fasta, request.POST.get("seqtype", "auto"), genetic_code(request),
            limits_for(request.user),
        )  # fmt: skip
    context = {"check": check, "files": input_views.sequence_files_context(request)}
    return render(request, "phylobot/_sequence_check.html", context)


def replicates(request: HttpRequest, box: str, field: str) -> int:
    """A branch-support replicate count from the form: 0 if unticked, else the count
    (bounded as the form bounds it; 1000 if missing or invalid)."""
    if request.POST.get(box) != "on":
        return 0
    try:
        return min(max(int(request.POST.get(field, "")), 1000), 10000)
    except ValueError:
        return 1000


def preview_reconciliation(request: HttpRequest, stop_after: str, amino: bool) -> dict[str, Any]:
    """Reconciliation as the summary describes and estimates it: on if ticked with
    somewhere to get the species tree from (the species table says whether every
    sequence has a species; SubmitForm decides at Submit)."""
    if stop_after == "alignments" or request.POST.get("reconcile") != "on":
        return {}
    if request.POST.get("species_consent") != "on" and not request.FILES.get("species_tree_file"):
        return {"reconcile_skipped": "the species names are not sent to the Open Tree of "
                "Life until you tick its box (or upload a species tree)"}  # fmt: skip
    params: dict[str, Any] = {
        "reconcile": True,
        "polytomies": "arbitrary" if request.POST.get("polytomies") == "arbitrary" else "guided",
    }
    builders = [t for t in request.POST.getlist("tree_builders") if t in ("iqtree", "raxml")]
    if request.POST.get("reconciled_tree") == "on" and amino:
        params["tree_builders"] = [*builders, RECONCILED_TREE]
    return params


@require_POST
def preview(request: HttpRequest) -> HttpResponse:
    """What the current form choices will run: combinations and ModelFinder candidates."""
    picks = {
        grid: [m for m in request.POST.getlist(grid) if catalog.is_model(m)]
        for grid in ("models", "candidates")
    }
    modelfinder = request.POST.get("modelfinder") == "on"
    fixed_models = request.POST.get("fixed_models") == "on"
    uploads = request.FILES.getlist("alignment_files")
    alignment_files = len(used_files(uploads, request.POST, "alignment"))
    alignment_files += len(request.POST.getlist("saved_alignments"))
    tree_files = len(request.FILES.getlist("tree_files"))
    # A duplicated job reuses its own uploads unless new ones are chosen (SubmitForm).
    source_id = request.POST.get("from_job", "")
    if source_id and is_uuid(source_id):
        source = Job.objects.for_user(request.user).filter(id=source_id).first()
        if source is not None:
            alignment_files = alignment_files or len(source.params.get("alignment_labels") or [])
            tree_files = tree_files or len(source.params.get("tree_labels") or [])
    # An extension: the parent's alignments are this job's first ones (SubmitForm.inherit).
    parent_id = request.POST.get("extend_job", "")
    parent = None
    if parent_id and is_uuid(parent_id):
        parent = Job.objects.for_user(request.user).filter(id=parent_id).first()
    if parent is not None and extend.refusal(parent):
        parent = None
    inherited, own = extend.counts(parent) if parent is not None else ([], 0)
    alignment_files += len(inherited) + own
    aligners = [a for a in request.POST.getlist("aligners") if a not in inherited]
    stop_after = SubmitForm.stop_after({
        "run_trees": request.POST.get("run_trees") == "on",
        "run_ancestors": request.POST.get("run_ancestors") == "on",
    })  # fmt: skip
    params = {
        "stop_after": stop_after,
        "aligners": aligners,
        "tree_builders": request.POST.getlist("tree_builders"),
        "modelfinder": modelfinder,
        "fixed_models": fixed_models,
        "models": picks["models"] if fixed_models else [],
        "candidates": picks["candidates"] if modelfinder else [],
        # Replicate counts (0 when unticked): describe() needs which measures run,
        # the work estimate how many replicates.
        "alrt": replicates(request, "sh_alrt", "alrt"),
        "bootstrap": replicates(request, "ufboot", "bootstrap"),
        "alignment_labels": ["x"] * alignment_files,
        "tree_labels": ["x"] * tree_files,
        "outgroup": request.POST.getlist("outgroup"),
    }
    if parent is not None:
        params["extends"] = {"name": parent.name, "alignments": len(inherited) + own}
    fasta = submitted_fasta(request)
    sequence_type = seqtype(request, fasta)
    params["seqtype"] = sequence_type
    if sequence_type == "codon":
        params["genetic_code"] = genetic_code(request)
    else:  # PRANK's codon model is for coding DNA only (the page hides it)
        params["aligners"] = [a for a in aligners if a != "prank_codon"]
    # As SubmitForm.clean: codon models only, amino-acid models only, or both.
    amino, codon = analysis_flags(sequence_type, request.POST.get("analyses", "both"))
    codon_ticked: list[str] = []  # codon ModelFinder's candidates
    if not amino:
        params |= {"modelfinder": False, "fixed_models": False, "models": [], "candidates": []}
    if codon:
        grids = {
            grid: [m for m in request.POST.getlist(grid) if catalog.is_codon_model(m)]
            for grid in ("codon_models", "codon_candidates")
        }
        if request.POST.get("codon_modelfinder") == "on":
            params["codon_modelfinder"] = True
            params["codon_candidates"] = codon_ticked = grids["codon_candidates"]
        if request.POST.get("codon_fixed") == "on":
            params["codon_models"] = grids["codon_models"]
        params |= codon_tree_params(
            request.POST.get("codon_search") == "on",
            request.POST.get("codon_on_amino") == "on",
            request.POST.get("codon_amino_trees"),
        )
        if request.POST.get("codeml") == "on" and stop_after == "ancestors":
            params["codeml_models"] = [
                m for m in request.POST.getlist("codeml_models") if m in codeml.MODELS
            ]
            trees = request.POST.get("codeml_trees")
            params["codeml_trees"] = trees if trees in codeml.WHERE else "best"
    raxml = "raxml" in request.POST.getlist("tree_builders")
    ticked = picks["candidates"]
    default_set = catalog.is_default_set(ticked, raxml=raxml)
    candidates = None if default_set else catalog.modelfinder_candidates(ticked)
    params["gap_methods"] = request.POST.getlist("gap_methods")
    params |= preview_reconciliation(request, stop_after, amino)
    if parent is not None and extend.mode(parent) == extend.RESULTS:
        # As Submit would queue it (SubmitForm.inherit_settings): the new aligners, and
        # what is added to every alignment; the rest is the parent's.
        adds = extend.additions(
            parent, request.POST.getlist("tree_builders"), request.POST.getlist("models"),
            request.POST.get("modelfinder") == "on",
        )  # fmt: skip
        if sequence_type != "codon":
            aligners = [a for a in aligners if a != "prank_codon"]
        try:
            params = extend.extended_params(
                extend.inherit(parent),
                aligners,
                adds,
                list(parent.params.get("tree_labels") or []),
            )  # fmt: skip
        except LookupError:
            raise Http404 from None
        alignment_files = 0  # the inherited are not this job's own
        ticked = list(parent.params.get("candidates") or catalog.all_models())
        raxml = "raxml" in (parent.params.get("tree_builders") or [])
        default_set = catalog.is_default_set(ticked, raxml=raxml)
        candidates = None if default_set else catalog.modelfinder_candidates(ticked)
    context = {
        "size": (size := workload.size_of(fasta, sequence_type == "codon")),
        "estimate": (
            estimate := workload.estimate(
                params, size or workload.IME2_SIZE, workload.seconds_per_unit()
            )
        ),
        "compute_time": None if estimate.minutes is None else eta.describe(estimate.minutes * 60),
        "steps": describe(params, alignment_files, tree_files),
        "alignment_count": alignment_count(params),
        "combinations": combination_count(params),
        "params": params,
        "candidates": candidates,
        "widened": candidates is not None and len(candidates.models()) > len(set(ticked)),
        "codon_candidates": (
            None
            if catalog.is_default_codon_set(codon_ticked)
            else catalog.codon_candidates(codon_ticked)
        ),  # fmt: skip
    }
    return render(request, "phylobot/_preview.html", context)


DELETED_SHOWN = 500  # most recent deleted jobs listed (the admin has them all)
JOBS_SHOWN = 200  # most recent jobs listed that match the filters
USAGE = ("cpu", "disk")  # My Jobs: the sort keys that share the Usage column
# My Jobs' columns, each with the sort keys whose links share its header: tags sit under
# a job's name, and the finish time (or ETA) under its status, so the table stays narrow.
JOB_COLUMNS = (
    ("job", ("name", "tags")), ("type", ("type",)), ("submitted", ("submitted",)),
    ("status", ("status", "finished")), ("usage", USAGE),
)  # fmt: skip


BACKFILL_PER_PAGE = 25  # older jobs whose sequences are counted per My Jobs visit


def count_sequences(shown: list[Job]) -> None:
    """Jobs from before `sequence_count` was recorded at Submit: count their input once
    (the sequences, else the first uploaded alignment) and keep it in their params."""
    lacking = [
        job for job in shown
        if "sequence_count" not in job.params and not job.params.get("simulate")
        and not job.params.get("branch_site") and not job.params.get("esr")
        and job.purged_at is None
    ][:BACKFILL_PER_PAGE]  # fmt: skip
    for job in lacking:
        count = 0
        for name in ("sequences.fasta", "alignment1.fasta"):
            try:
                text = storage.read_bytes(f"{job.storage_prefix}inputs/{name}")
            except Exception:  # missing (or unreadable) input: try the next
                continue
            count = sum(1 for line in text.splitlines() if line.startswith(b">"))
            break
        job.params["sequence_count"] = count or None
        Job.objects.filter(id=job.id).update(params=job.params)


def job_label(job: Job) -> str:
    return job.name or str(job.id)[:8]


def my_jobs(request: HttpRequest) -> HttpResponse:
    # "Group follow-up jobs": a display preference like the skin, so the toggle's link
    # (?group=on|off) is remembered in a cookie and the address stays the list's own.
    if (asked := request.GET.get("group")) in ("on", "off"):
        query = request.GET.copy()
        query.pop("group")
        here = reverse("my_jobs") + (f"?{query.urlencode()}" if query else "")
        response = redirect(here)
        response.set_cookie(
            worklist.NEST_COOKIE, "1" if asked == "on" else "0",
            max_age=365 * 24 * 3600, samesite="Lax",
        )  # fmt: skip
        return response
    grouped = request.COOKIES.get(worklist.NEST_COOKIE) == "1"
    everything = Job.objects.for_user(request.user)
    filters = job_types.Filters.parse(request.GET)
    sort = job_types.Sort.parse(request.GET)
    deleted = Job.objects.deleted_for(request.user).prefetch_related("tags", "reruns")
    # The Type filter may ask only for homolog searches (a pseudo-type that Filters drops):
    # then no pipeline jobs are listed.
    requested_types = request.GET.getlist("type")
    # Only the explicit "homolog" pseudo-type lists searches alone; an unknown value is just
    # dropped (Filters already ignores it), so jobs still show.
    only_searches = worklist.HOMOLOG in requested_types and not any(
        t in job_types.BY_KEY for t in requested_types
    )
    jobs_qs = everything.none() if only_searches else filters.apply(everything)
    shown = list(sort.apply(jobs_qs).prefetch_related("tags")[:JOBS_SHOWN])
    etas: dict[Any, str] = {}
    if any(job.status == Job.Status.RUNNING for job in shown):
        typical = eta.typical_durations()  # once per page, not per job
        etas = {
            job.id: eta.describe(eta.remaining_seconds(job, typical), short=True)
            for job in shown
            if job.status == Job.Status.RUNNING
        }
    count_sequences(shown)
    # Extensions that read results from a listed job (phylobot.joined), for the delete
    # dialog: one query for the page.
    readers: dict[str, list[str]] = {}
    for child in everything.filter(
        deleted_at__isnull=True, params__inherits__job__in=[str(job.id) for job in shown]
    ):
        readers.setdefault(child.params["inherits"]["job"], []).append(job_label(child))
    # Homolog searches interleave with the jobs (phylobot.worklist, Approach A).
    job_rows = [
        worklist.job_row(
            job,
            job_types.types_of(job.params),
            etas.get(job.id, ""),
            job_types.composition(job.params),
            readers.get(str(job.id), []),
        )  # fmt: skip
        for job in shown
    ]
    search_rows = (
        worklist.searches_for(request.user, filters)
        if worklist.want_searches(filters, requested_types)
        else []
    )
    rows = worklist.order(job_rows + search_rows, sort)[:JOBS_SHOWN]
    if grouped:
        rows = worklist.nest(rows)
    headers = sort.headers(request.GET)
    group_query = request.GET.copy()
    group_query["group"] = "off" if grouped else "on"
    context = {
        "jobs": shown,
        "rows": rows,
        "grouped": grouped,
        "group_query": group_query.urlencode(),
        "has_searches": bool(search_rows),
        "deleted": [
            (job, job_types.types_of(job.params))
            for job in deleted.order_by("-deleted_at")[:DELETED_SHOWN]
        ],
        "filters": filters,
        "homolog_type": worklist.HOMOLOG,
        "homolog_label": worklist.SEARCH_LABEL,
        "homolog_selected": worklist.HOMOLOG in requested_types,
        "sort": sort,
        "columns": [
            {
                "key": column,
                "links": [h for h in headers if h["key"] in keys],
                "sort": next((h["sort"] for h in headers if h["key"] in keys and h["sort"]), ""),
            }
            for column, keys in JOB_COLUMNS
        ],
        # Counts span both kinds: total is everything the user has, matching is what passes
        # the filters (a job-only filter leaves every search out of matching).
        "total": everything.count() + HomologSearch.objects.for_user(request.user).count(),
        "matching": (0 if only_searches else jobs_qs.count()) + len(search_rows),
        "job_types": job_types.TYPES,
        "statuses": Job.Status.values,
        "data_choices": job_types.DATA,
        "origin_choices": job_types.ORIGIN,
        "all_tags": Tag.objects.filter(owner_id=request.user.pk).values_list("name", flat=True),
        "can_import": settings.PHYLOBOT_MODE == "local",
        "import_limit": job_import.size(site_settings.import_limits().upload_bytes),
    }
    return render(request, "phylobot/my_jobs.html", context)


def safe_next(request: HttpRequest, default: str) -> str:
    """Where to go after a POST: a My Jobs URL from the form, never an outside one."""
    target = request.POST.get("next", "")
    if target.startswith("/jobs/") and "//" not in target and "\\" not in target:
        return target
    return default


@require_POST
def bulk_tags(request: HttpRequest) -> HttpResponse:
    """Add or remove tags on the selected jobs (My Jobs toolbar)."""
    back = safe_next(request, reverse("my_jobs"))
    ids = [value for value in request.POST.getlist("job") if is_uuid(value)]
    selected = list(Job.objects.for_user(request.user).filter(id__in=ids))
    try:
        names = parse_tags(request.POST.get("tags", ""))
    except ValueError as exc:
        messages.error(request, str(exc))
        return redirect(back)
    if not selected or not names:
        messages.error(request, "Select one or more jobs and enter a tag.")
        return redirect(back)
    adding = request.POST.get("action") == "add"
    for job in selected:
        job.add_tags(names) if adding else job.remove_tags(names)
    labels = " ".join(f"#{n}" for n in names)
    plural = "s" if len(selected) != 1 else ""
    verb = f"Added {labels} to" if adding else f"Removed {labels} from"
    messages.success(request, f"{verb} {len(selected)} job{plural}.")
    return redirect(back)


@require_POST
def job_tags(request: HttpRequest, job_id: str) -> HttpResponse:
    """Add or remove tags on one job (job page); returns the tags fragment for HTMX."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    error = ""
    try:
        names = parse_tags(request.POST.get("tags", ""))
    except ValueError as exc:
        names, error = [], str(exc)
    if request.POST.get("action") == "remove":
        job.remove_tags(names)
    elif names:
        job.add_tags(names)
    if request.headers.get("HX-Request"):
        return render(request, "phylobot/_job_tags.html", {"job": job, "tag_error": error})
    if error:
        messages.error(request, error)
    return redirect("job_detail", job_id=job.id)


@require_POST
def delete_jobs(request: HttpRequest) -> HttpResponse:
    """Delete the selected jobs and homolog searches (both are listed on My Jobs,
    phylobot.worklist): their files go, their records stay (see jobs.delete / searches.delete).
    IDs that are not the user's own, live items are ignored."""
    ids = [value for value in request.POST.getlist("job") if is_uuid(value)]
    selected = list(Job.objects.for_user(request.user).filter(id__in=ids))
    count = jobs.delete(selected)
    search_ids = [value for value in request.POST.getlist("search") if is_uuid(value)]
    chosen = list(HomologSearch.objects.for_user(request.user).filter(id__in=search_ids))
    for search in chosen:
        searches.delete(search)
    parts = []
    if count:
        parts.append(f"{count} job{'s' if count != 1 else ''}")
    if chosen:
        parts.append(f"{len(chosen)} search{'es' if len(chosen) != 1 else ''}")
    if parts:
        messages.success(request, f"Deleted {' and '.join(parts)}.")
    return redirect(safe_next(request, reverse("my_jobs")))


@require_POST
def rename_job(request: HttpRequest, job_id: str) -> HttpResponse:
    """Rename from the job page (HTMX: the title part comes back, done editing)."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    jobs.rename(job, request.POST.get("name", ""))
    if request.headers.get("HX-Request"):
        return render(request, "phylobot/_job_title.html", {"job": job})
    return redirect("job_detail", job_id=job.id)


@require_POST
def stop_job(request: HttpRequest, job_id: str) -> HttpResponse:
    """Stop a queued or running job (jobs.stop): its runner stops the containers at its
    next heartbeat. Refused quietly if it had already finished. A stopped job is
    cancelled, so it cannot be resumed; its files are kept."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if jobs.stop(job):
        messages.success(request, "Stopping the job; its steps stop within a few seconds.")
    return redirect(safe_next(request, reverse("job_detail", args=[job.id])))


@require_POST
def resume_job(request: HttpRequest, job_id: str) -> HttpResponse:
    """Resume a failed job as its next attempt (jobs.resume): its finished steps are
    reused. Refused quietly if it cannot be resumed (or someone already did)."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if jobs.resume(job):
        messages.success(request, f"Resumed as attempt {job.attempt}.")
    return redirect("job_detail", job_id=job.id)


@require_GET
def job_download(request: HttpRequest, job_id: str) -> StreamingHttpResponse:
    """The job's inputs and results as one ZIP."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    response = StreamingHttpResponse(bundle.zip_stream(job), content_type="application/zip")
    response["Content-Disposition"] = f'attachment; filename="{bundle.folder_name(job)}.zip"'
    return response


def ancestors_zip(request: HttpRequest, job_id: str) -> StreamingHttpResponse:
    """One combination's ancestor files as a ZIP (the job page's "download" per row):
    ?combination=<key>. The key must be one of the job's combinations with ancestors."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    key = request.GET.get("combination", "")
    files = joined.results(job)
    if not key or f"{key}/ancestors.fasta" not in files.keys:
        raise Http404
    paths = results.ancestor_files(key, files.keys)
    named = ancestor_names_tsv(job, files, key)
    folder = re.sub(r"[^A-Za-z0-9._+-]", "_", f"{bundle.folder_name(job)}_{key}_ancestors")
    response = StreamingHttpResponse(
        bundle.results_zip(job, f"{key}/", paths, folder, named), content_type="application/zip"
    )
    response["Content-Disposition"] = f'attachment; filename="{folder}.zip"'
    return response


def results_export(request: HttpRequest, job_id: str, anchor: str, fmt: str) -> HttpResponse:
    """One of the job page's results tables (an alignment and data type, by its anchor)
    as CSV or Excel, with unrounded statistics and each row's combination."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if fmt not in ("csv", "xlsx") or job.status != Job.Status.SUCCEEDED:
        raise Http404
    files = joined.results(job)
    if joined.SUMMARY not in files.keys:
        raise Http404
    found = [r for r in alignment_results(job, files.summary, files.keys, files)
             if r.anchor == anchor and r.table.export_rows]  # fmt: skip
    if not found:
        raise Http404
    table = found[0].table
    stem = re.sub(
        r"[^A-Za-z0-9._+-]",
        "_",
        f"{bundle.folder_name(job)}_{anchor.removeprefix('alignment-')}_models",
    )
    if fmt == "csv":
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(table.export_columns)
        writer.writerows(table.export_rows)
        response = HttpResponse(out.getvalue(), content_type="text/csv; charset=utf-8")
    else:
        from phylobot import sites_excel  # openpyxl: only needed here

        about = [("Job", job.name or str(job.id)), ("Job ID", str(job.id)),
                 ("Table", table.title), ("Notes", table.subtitle)]  # fmt: skip
        response = HttpResponse(
            sites_excel.results_workbook(table.export_columns, table.export_rows, about),
            content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )
    response["Content-Disposition"] = f'attachment; filename="{stem}.{fmt}"'
    return response


def ancestor_names_tsv(job: Job, files: joined.Results, key: str) -> dict[str, str]:
    """{"ancestor_names.tsv": text} for a combination whose nodes have names, else {}:
    the names live in PhyloBot, not in the pipeline's files."""
    if analogs.MAP not in files.keys or not ancestor_names.exist(job):
        return {}
    by_clade = bool(job.params.get("outgroup"))
    nodes = analogs.placed(analogs.read_map(files.read(analogs.MAP)), by_clade).get(key, {})
    names = ancestor_names.for_nodes(job, nodes, by_clade)
    if not names:
        return {}
    lines = ["node\tname\tsame_branching"]
    for node, named in sorted(names.items()):
        lines.append(f"{node}\t{named.name}\t{'yes' if named.alike else 'no'}")
    return {"ancestor_names.tsv": "\n".join(lines) + "\n"}


def local_only() -> None:
    """Import is for local mode only for now: cloud mode has no storage quotas yet."""
    if settings.PHYLOBOT_MODE != "local":
        raise Http404


@require_POST
def import_job(request: HttpRequest) -> HttpResponse:
    """Import a job from a ZIP downloaded from any PhyloBot (phylobot.job_import). If one
    of the user's jobs has its name, the upload waits in storage while they choose
    another (`import_job_name`)."""
    local_only()
    upload = request.FILES.get("archive")
    if upload is None:
        messages.error(request, "Choose a ZIP to import.")
        return redirect("my_jobs")
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "job.zip"
        try:
            job_import.save_upload(upload, path)
            archive = job_import.read(path)
            name = job_import.clean_name(archive.name)
            if job_import.clashes(request.user, name):
                token = job_import.stage(request.user, path)
                return redirect("import_job_name", token=token)
            job = job_import.import_job(request.user, archive, name)
        except job_import.Refused as refusal:
            messages.error(request, f"Could not import that ZIP. {refusal}")
            return redirect("my_jobs")
    messages.success(request, f"Imported “{job.name or 'Job'}”.")
    return redirect("job_detail", job_id=job.id)


@require_http_methods(["GET", "POST"])
def import_job_name(request: HttpRequest, token: str) -> HttpResponse:
    """Name an imported job whose name one of the user's jobs already has."""
    local_only()
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "job.zip"
        if not job_import.fetch(request.user, token, path):
            messages.error(request, "That import is no longer waiting: choose the ZIP again.")
            return redirect("my_jobs")
        try:
            archive = job_import.read(path)
        except job_import.Refused as refusal:
            job_import.discard(request.user)
            messages.error(request, f"Could not import that ZIP. {refusal}")
            return redirect("my_jobs")
        original = job_import.clean_name(archive.name)
        name, error = job_import.suggest_name(request.user, original), ""
        if request.method == "POST":
            name = job_import.clean_name(request.POST.get("name", ""))
            if not name:
                error = "Give the imported job a name."
            elif job_import.clashes(request.user, name):
                error = f"You already have a job named “{name}”."
            else:
                try:
                    job = job_import.import_job(request.user, archive, name)
                except job_import.Refused as refusal:
                    messages.error(request, f"Could not import that ZIP. {refusal}")
                    return redirect("my_jobs")
                finally:
                    job_import.discard(request.user)
                messages.success(request, f"Imported “{job.name}”.")
                return redirect("job_detail", job_id=job.id)
    return render(
        request,
        "phylobot/import_job.html",
        {"token": token, "original": original, "name": name, "error": error,
         "record": archive.job},
    )  # fmt: skip


@require_POST
def import_job_cancel(request: HttpRequest, token: str) -> HttpResponse:
    local_only()
    job_import.discard(request.user)
    return redirect("my_jobs")


def is_uuid(text: str) -> bool:
    try:
        uuid.UUID(text)
    except ValueError:
        return False
    return True


INPUT_PREVIEW_BYTES = 200_000  # larger inputs are offered as a download only


def job_detail(request: HttpRequest, job_id: str) -> HttpResponse:
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    context = job_context(job) | input_context(job)
    # Ancestral libraries need reconstructed ancestors: a finished job that ran that far.
    context["libraries_possible"] = (
        job.status == Job.Status.SUCCEEDED
        and not job.params.get("branch_site")
        and not job.params.get("simulate")
        and not job.params.get("esr")
        and job.params.get("stop_after", "ancestors") == "ancestors"
    )
    context["libraries"] = library_views.libraries_of(request.user, job)
    context["structures"] = list(Structure.objects.filter(job=ancestor_names.holder(job)))
    context["dependents"] = [job_label(j) for j in joined.dependents(job)]
    context["extensions"] = list(job.extensions.filter(deleted_at__isnull=True))
    return render(request, "phylobot/job_detail.html", context)


def input_context(job: Job) -> dict[str, Any]:
    """The job's uploaded inputs, as submitted (not in the polled status fragment)."""
    prefix = f"{job.storage_prefix}inputs/"
    try:
        keys = [k[len(prefix) :] for k in storage.list_keys(prefix)]
    except Exception:
        logger.exception("job %s: could not list inputs", job.id)
        keys = []
    context: dict[str, Any] = {"inputs": keys, "input_text": "", "input_stats": None}
    if "sequences.fasta" in keys:
        data = storage.read_bytes(prefix + "sequences.fasta")
        context["input_size"] = len(data)
        if len(data) <= INPUT_PREVIEW_BYTES:
            context["input_text"] = data.decode(errors="replace")
    if manifest.has(job, "validation.json"):
        report_key = f"{job.storage_prefix}results/validation.json"
        context["input_stats"] = json.loads(storage.read_bytes(report_key))
    return context


@require_GET
def job_folder(request: HttpRequest, job_id: str, folder: str) -> HttpResponse:
    """One folder of a job's results, for jobs too large to list on their page."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    relative = PurePosixPath(folder)
    if relative.is_absolute() or ".." in relative.parts or not relative.parts:
        raise Http404
    under = relative.as_posix() + "/"
    files = joined.cached(job)
    sizes = files.sizes(under)
    if not sizes:
        raise Http404
    context = {
        "job": job,
        "folder": under,
        "crumbs": filelist.breadcrumbs(under),
        "listing": filelist.folder_listing(sizes, under, files.summary),
    }
    return render(request, "phylobot/job_folder.html", context)


def job_status(request: HttpRequest, job_id: str) -> HttpResponse:
    """The live part of the job page, polled by HTMX until the job finishes."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    # refresh: the page's "On this page" nav is replaced too (out of band), as the
    # results it lists appear.
    return render(request, "phylobot/_job_status.html", job_context(job) | {"refresh": True})


def extend_label(job: Job) -> str:
    """What extending the job would add, for the Extend button below its steps graph
    ("" if it can't be extended)."""
    if extend.refusal(job):
        return ""
    if extend.mode(job) == extend.ALIGNMENTS:
        return "more aligners, trees and ancestors"
    return "more aligners, tree sources and models"


def job_context(job: Job) -> dict[str, Any]:
    results_prefix = f"{job.storage_prefix}results/"
    finished = job.status in jobs.TERMINAL
    # From the manifest (phylobot.manifest), never a storage listing. A running job's
    # files go up when it ends; only its live log is there before.
    own = manifest.sizes(job)
    if not finished and LIVE_LOG not in own:
        try:
            if storage.exists(results_prefix + LIVE_LOG):
                own[LIVE_LOG] = 0
        except Exception:
            logger.exception("job %s: could not look for its log", job.id)
    keys = sorted(own)
    progress = job.progress or {}
    done, total = progress.get("done", 0), progress.get("total", 0)
    end = job.finished_at or timezone.now()
    context: dict[str, Any] = {
        "job": job,
        "finished": finished,
        "done": done,
        "total": total,
        "percent": round(100 * done / total) if total else 0,
        "last_task": progress.get("last", ""),
        "elapsed": duration((end - job.started_at).total_seconds()) if job.started_at else "",
        "all_tags": Tag.objects.filter(owner_id=job.owner_id).values_list("name", flat=True),
        "eta": None if finished else eta.describe(eta.remaining_seconds(job)),
        "combinations": combination_count(job.params),
        "extendable": not extend.refusal(job),
        "plan": describe(
            job.params,
            len(job.params.get("alignment_labels") or []),
            len(job.params.get("tree_labels") or []),
        ),
        "log_tail": "",
        "files": None,  # below, with the summary's labels if it has one
        "file_count": len(keys),
        "indexing": finished and job.results_indexed_at is None and job.purged_at is None,
        "summary": None,
        "alignment_results": [],
        "outgroup": job.params.get("outgroup") or [],
        "stages": stages.graph(job.params, progress, finished),
        "extend_label": extend_label(job),
        "failure": failure(job, keys) if job.status == Job.Status.FAILED else None,
        "resumable": jobs.resumable(job),
        "attempts": attempt_rows(job),
    }
    if LIVE_LOG in keys:
        text = storage.read_bytes(results_prefix + LIVE_LOG).decode(errors="replace")
        context["log_tail"] = "\n".join(text.splitlines()[-LOG_TAIL_LINES:])
    if job.params.get("simulate"):
        context |= simulation_views.context(job, set(keys))
    if job.params.get("esr"):
        context |= esr_views.context(job, set(keys))
    if job.params.get("properties_only"):
        context["properties_run"] = properties_submit.context(job)
    if job.params.get("uncertainty_only"):
        context["uncertainty_run"] = uncertainty_submit.context(job)
    if job.params.get("branch_site"):
        context |= branch_site_views.context(
            job, keys, lambda relative: storage.read_bytes(results_prefix + relative).decode()
        )
    if job.status == Job.Status.SUCCEEDED and "summary.json" in keys:
        # With the results it inherits (phylobot.joined): the page of an extension is
        # the page of one job that ran every alignment.
        files = joined.results(job, set(keys))
        keys = sorted(files.keys)
        summary = files.summary
        context["summary"] = summary
        context["inherited"] = files
        context["files"] = filelist.folder_listing(
            files.sizes() if files.inherited else own, "", summary
        )
        context["file_count"] = len(keys)
        if files.inherited:  # the joined job's, not only the ones this job ran
            context["combinations"] = len(summary.get("combinations", []))
        context["alignment_results"] = alignment_results(job, summary, files.keys, files)
        context["with_codons"] = any(r.data == "codon" for r in context["alignment_results"])
        context |= indel_context(job, summary, files, context["alignment_results"])
    else:
        context["files"] = filelist.folder_listing(own)
    return context


def indel_context(
    job: Job,
    summary: dict[str, Any],
    files: joined.Results,
    results_list: list[results.AlignmentResult],
) -> dict[str, Any]:
    """The job page's row of indel-length thumbnails. An alignment that cannot be measured
    is left out: the histograms are never a reason for the page to fail."""
    measured = []
    sources = alignment_views.sources(job, files.keys, summary)
    shown = [r.aligner for r in results_list if r.first]  # the table above them, in its order
    for source in sorted(
        sources, key=lambda s: shown.index(s.key) if s.key in shown else len(shown)
    ):
        if not source.amino:
            continue
        try:
            measured.append(indel_views.measure(job, source.key, source.label, source.amino,
                                                 summary, files.keys, files.read))  # fmt: skip
        except Exception:
            logger.exception("job %s: could not measure indels in %s", job.id, source.key)
    return {
        "indel_thumbnails": indel_views.thumbnails(measured),
    }


def attempt_rows(job: Job) -> list[dict[str, Any]]:
    """The job's attempts for its page, once there is more than one."""
    if job.attempt < 2:
        return []
    rows = []
    for attempt in job.attempts.all():
        used = None
        if attempt.completed_cpu_seconds is not None:
            used = attempt.completed_cpu_seconds + (attempt.failed_cpu_seconds or 0)
        end = attempt.finished_at or (timezone.now() if attempt.number == job.attempt else None)
        rows.append(
            {
                "attempt": attempt,
                "status": attempt.status or job.status,
                "elapsed": duration((end - attempt.started_at).total_seconds())
                if end and attempt.started_at
                else "",
                "used": used,
            }
        )
    return rows


def failure(job: Job, keys: list[str]) -> diagnosis.Diagnosis | None:
    """Why a failed job failed (phylobot.diagnosis), from its stored logs."""
    prefix = f"{job.storage_prefix}results/"

    def read(relative: str) -> str:
        if relative not in keys:
            return ""
        return storage.read_bytes(prefix + relative)[-200_000:].decode(errors="replace")

    try:
        logs = diagnosis.Logs(read("logs/pipeline.log"), read("logs/nextflow.log"))
        name, _ = diagnosis.failed_task(job.progress or {}, logs)
        if name:
            stem = f"logs/tasks/{log_name(name)}"
            logs.tasks[name], logs.outputs[name] = read(f"{stem}.err"), read(f"{stem}.log")
        return diagnosis.diagnose(job.error, job.progress or {}, logs)
    except Exception:  # never a reason the job page fails
        logger.exception("job %s: could not diagnose the failure", job.id)
        return None


def alignment_results(
    job: Job, summary: dict[str, Any], keys: Collection[str], files: joined.Results | None = None
) -> list[results.AlignmentResult]:
    """The per-alignment results tables (see phylobot.results)."""
    files = files or joined.Results(job)
    read = files.read

    comparison = read("model_comparison.tsv") if "model_comparison.tsv" in keys else None
    combos = summary.get("combinations", [])
    # One small file per combination: read together (285 of each on the GFP job).
    per_combination = files.read_many(
        [f"{c['key']}/{name}" for c in combos for name in ("modelfinder.tsv", "rooting.json")
         if f"{c['key']}/{name}" in keys]
    )  # fmt: skip
    modelfinder = {
        c["key"]: per_combination[f"{c['key']}/modelfinder.tsv"]
        for c in combos
        if f"{c['key']}/modelfinder.tsv" in per_combination
    }
    # One alignment file per aligner: <aligner>/alignment.fasta (jobs from before it
    # existed: the copy in each of the aligner's combinations).
    alignments: dict[str, str] = {}
    for entry in summary.get("alignments", []):
        if entry["path"] in keys:
            alignments[entry["key"]] = read(entry["path"])
    for c in combos:
        path = f"{c['key']}/alignment.fasta"
        if c["aligner"] not in alignments and path in keys:
            alignments[c["aligner"]] = read(path)
    viewer = reverse("ancestor_viewer", args=[job.id])
    tree = reverse("tree_viewer", args=[job.id])
    zipped = reverse("ancestors_zip", args=[job.id])
    rooting = {
        c["key"]: json.loads(per_combination[f"{c['key']}/rooting.json"])
        for c in combos
        if f"{c['key']}/rooting.json" in per_combination
    }
    reconciliation = {
        c["key"]: json.loads(read(f"{c['key']}/reconciliation.json"))
        for c in combos
        if f"{c['key']}/reconciliation.json" in keys
    }
    for key, totals in reconciliation.items():  # weak duplications (reconcile.overlaps)
        if {f"{key}/reconciliation.tsv", f"{key}/reconciliation.nwk"} <= set(keys):
            totals["weak"] = reconcile.weak_duplications(
                read(f"{key}/reconciliation.tsv"), read(f"{key}/reconciliation.nwk")
            )
    topology: dict[str, dict[str, str]] = {}
    for aligner, model in dict.fromkeys((c["aligner"], c["model"]) for c in combos):
        path = f"{aligner}/topology_tests/{model}.tsv"
        if path in keys:
            topology |= results.topology_rows(aligner, model, read(path))
    out = results.alignment_results(
        summary,
        comparison,
        modelfinder,
        lambda relative: reverse("job_file", args=[job.id, relative]),
        alignments,
        set(keys),
        rooting,
        lambda key: f"{viewer}?{urlencode({'combination': key})}",
        lambda key: f"{tree}?{urlencode({'combination': key})}",
        reconciliation,
        topology,
        lambda key: f"{zipped}?{urlencode({'combination': key})}",
    )
    compare_codon_trees(out, read, set(keys))
    available = set(keys)
    for table in out:
        if table.data != "codon":
            continue
        for combo in combos:
            folder = f"{combo['key']}/codeml"
            if combo["aligner"] != table.aligner or f"{folder}/codeml.tsv" not in available:
                continue
            table.codeml.append(
                results.codeml_result(
                    combo,
                    read(f"{folder}/codeml.tsv"),
                    read(f"{folder}/selection.tsv"),
                    read(f"{folder}/beb.tsv"),
                    lambda relative: reverse("job_file", args=[job.id, relative]),
                    lambda key, model: (
                        f"{viewer}?"
                        + urlencode({"combination": key, "gaps": ancestors.CODEML_PREFIX + model})
                    ),
                    available,
                )
            )
    return out


def compare_codon_trees(
    tables: list[results.AlignmentResult], read: Any, available: set[str]
) -> None:
    """For each alignment with both kinds of model: how far the best codon-model tree is
    from the best amino-acid-model tree (unrooted Robinson-Foulds)."""
    best = {(t.aligner, t.data): t for t in tables}
    for (aligner, data), table in best.items():
        amino = best.get((aligner, "amino"))
        if data != "codon" or amino is None:
            continue
        source = table.best.split("/")[1]
        if source == amino.best.split("/")[1] and source.startswith("provided"):
            table.comparison = (
                "Both best trees are fits of the same tree of yours (topology kept), so they "
                "have the same topology; add IQ-TREE's search to compare the topologies "
                "each kind of model finds."
            )
            continue
        paths = [f"{key}/tree.nwk" for key in (table.best, amino.best)]
        if not all(path in available for path in paths):
            continue
        codon_tree, amino_tree = (read(path) for path in paths)
        leaves = len(trees.load_tree(codon_tree).leaf_nodes())
        distance = trees.robinson_foulds(codon_tree, amino_tree)
        table.comparison = (
            f"Best codon-model tree vs best amino-acid-model tree: Robinson-Foulds distance "
            f"{distance} of {2 * (leaves - 3)} (0 = the same topology)."
        )


def job_file(request: HttpRequest, job_id: str, path: str) -> StreamingHttpResponse:
    """A result file, e.g. /jobs/<id>/files/summary.json."""
    return serve_job_file(request, job_id, "results", path)


def job_input(request: HttpRequest, job_id: str, path: str) -> StreamingHttpResponse:
    """An uploaded input, as submitted, e.g. /jobs/<id>/inputs/sequences.fasta."""
    return serve_job_file(request, job_id, "inputs", path)


def serve_job_file(
    request: HttpRequest, job_id: str, area: str, path: str
) -> StreamingHttpResponse:
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    relative = PurePosixPath(path)
    if relative.is_absolute() or ".." in relative.parts:
        raise Http404
    prefix = f"{job.storage_prefix}{area}/"
    key = prefix + relative.as_posix()
    body: Any = None
    if area == "results":  # wherever the file lives (phylobot.joined)
        files = joined.results(job)
        if relative.as_posix() not in files.keys:
            raise Http404
        key = files.key(relative.as_posix())
        if files.inherited and relative.as_posix() in (joined.SUMMARY, *joined.JOINED_TABLES):
            body = iter([files.read_bytes(relative.as_posix())])  # joined, not one file
    elif key not in set(storage.list_keys(prefix)):
        raise Http404
    filename = re.sub(r"[^A-Za-z0-9._+-]", "_", relative.name)
    response = StreamingHttpResponse(body if body is not None else storage.stream(key))
    if relative.suffix in TEXT_SUFFIXES:
        response["Content-Type"] = "text/plain; charset=utf-8"
        response["Content-Disposition"] = f'inline; filename="{filename}"'
    else:
        response["Content-Type"] = "application/octet-stream"
        response["Content-Disposition"] = f'attachment; filename="{filename}"'
    response["Content-Security-Policy"] = "sandbox"
    return response


@dataclass
class AncestorSelection:
    """A validated choice of combination, gap version, node and extant sequence for the
    ancestor viewer. Every value is one of the job's own (never raw request text)."""

    job: Job
    summary: dict[str, Any]
    keys: set[str]
    combination: dict[str, Any]
    versions: list[str]
    version: str
    ancestors: dict[str, str]  # the chosen version's rows ({} without ancestors)
    node_names: list[str]  # nodes with a reconstruction ([] without ancestors)
    requires: str = "ancestors.fasta"  # the file every offered combination has
    files: joined.Results | None = None  # where each file lives (set by select_ancestor)

    @property
    def results(self) -> joined.Results:
        return self.files or joined.Results(self.job, self.keys, self.summary)

    def read_ancestor(self, name: str) -> str:
        """An ancestor file (ancestors_pp.tsv, ancestors_posteriors.tsv) of the chosen
        version: IQ-TREE's in the combination's folder, codeml's in codeml/<model>/."""
        return self.read(ancestors.version_folder(self.version) + name)

    def read_node(self, name: str, node: str) -> str:
        """One node's rows of a posteriors table (with its header), by ranged reads
        through EXPORT's byte index; the whole table where there is no index. The
        callers parse either the same way."""
        table = ancestors.version_folder(self.version) + name
        index = f"{self.combination['key']}/{ancestors.index_name(table)}"
        if index in self.keys:
            found = ancestors.read_index(self.results.read(index))
            path = f"{self.combination['key']}/{table}"
            text = found.text(node, lambda start, n: self.results.read_range(path, start, n))
            if text is not None:
                return text
        return self.read_ancestor(name)

    def read(self, name: str, combination: str = "") -> str:
        """A result file of this combination (or another of the job's)."""
        return self.results.read(f"{combination or self.combination['key']}/{name}")

    def other(self, combination: str, version: str) -> AncestorSelection | None:
        """The same job's `version` of another combination's ancestors, from the files
        already listed (no further listing of storage, which is the slow part). None if
        the combination or that version is not the job's. Its `node_names` are that
        version's nodes, without reading the original reconstruction as well."""
        combo = next(
            (c for c in self.summary.get("combinations", []) if c["key"] == combination), None
        )
        if combo is None or f"{combination}/{self.requires}" not in self.keys:
            return None
        files = version_files(combination)
        versions = [v for v in files if files[v] in self.keys]
        if version not in versions:
            return None
        rows = ancestors.parse_fasta(self.results.read(files[version]))
        return replace(
            self, combination=combo, versions=versions, version=version, ancestors=rows,
            node_names=list(rows),
        )  # fmt: skip

    def tree(self, combination: str = "") -> str:
        """The combination's tree: outgroup-rooted if there is one (descendants, and so
        the node menu, only mean something on a rooted tree)."""
        combination = combination or self.combination["key"]
        rooted = f"{combination}/tree.rooted.nwk" in self.keys
        return self.read("tree.rooted.nwk" if rooted else "tree.nwk", combination)


def version_files(combination: str) -> dict[str, str]:
    """{version: its ancestors file under results/} for one combination, whether or not
    the job wrote it."""
    files = {v: f"{combination}/{ancestors.VERSIONS[v][1]}" for v in ancestors.VERSIONS}
    for model in codeml.MODELS:  # codeml's reconstructions, where it ran
        files[ancestors.CODEML_PREFIX + model] = f"{combination}/codeml/{model}/ancestors.fasta"
    return files


def select_ancestor(
    job: Job, query: Any, requires: str = "ancestors.fasta"
) -> AncestorSelection | None:
    """Resolve the viewer's query parameters against the job's results. Unknown or
    missing values fall back to defaults. Only combinations with `requires` are offered
    (the tree viewer asks for tree.nwk: trees-only jobs have no ancestors); None if
    there are none."""
    if job.status != Job.Status.SUCCEEDED:
        return None
    found = joined.results(job)  # with the results it inherits
    keys, summary = found.keys, found.summary
    if "summary.json" not in keys or not summary:
        return None
    primary = summary.get("primary", "")
    combos = sorted(summary.get("combinations", []), key=results.combination_order(primary))
    combos = [c for c in combos if f"{c['key']}/{requires}" in keys]
    if not combos:
        return None
    by_key = {c["key"]: c for c in combos}
    requested = by_key.get(query.get("combination", ""))
    aligner = query.get("aligner", "")
    if aligner not in {c["aligner"] for c in combos}:
        aligner = (requested or by_key.get(primary) or combos[0])["aligner"]
    mine = [c for c in combos if c["aligner"] == aligner]
    # Default: the alignment's best combination by BIC (rank 1; ModelFinder's row on a
    # tie, as the results table merges them). Jobs without ranks keep display order.
    fallback = min(
        mine, key=lambda c: (c.get("rank_bic", 0), not catalog.is_modelfinder(c["model"]))
    )
    combination: dict[str, Any] = requested if requested and requested in mine else fallback
    files = version_files(combination["key"])
    versions = [v for v in files if files[v] in keys]
    version = query.get("gaps", "")
    if version not in versions:
        # The job's first chosen gap method, else the unmasked reconstruction.
        chosen = [m for m in job.params.get("gap_methods") or [] if m in versions]
        version = chosen[0] if chosen else ancestors.ORIGINAL

    def fasta(version: str) -> dict[str, str]:
        if files[version] not in keys:
            return {}
        return ancestors.parse_fasta(found.read(files[version]))

    node_names = list(fasta(ancestors.ORIGINAL))
    rows = fasta(version)
    return AncestorSelection(
        job, summary, keys, combination, versions, version, rows, node_names, requires, found
    )


def ancestor_viewer(request: HttpRequest, job_id: str) -> HttpResponse:
    """One ancestral node, site by site, with posterior probabilities (Fig. 5C of the
    2016 paper). The controls re-request this view; HTMX gets only the viewer part."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    selection = select_ancestor(job, request.GET)
    context: dict[str, Any] = {"job": job, "selection": selection}
    if selection is not None:
        context |= ancestor_context(selection, request.GET)
        context["library_nodes"] = [context["node"]]  # for "Add to library"
    htmx = bool(request.headers.get("HX-Request"))
    # The + / × buttons are one-off actions: the address bar gets the resulting view,
    # never the button, so reloading or sharing it does not add or remove again.
    one_off = "add_extant" in request.GET or "remove_extant" in request.GET
    if selection is not None and one_off and not htmx:
        return redirect(context["canonical_url"])
    template = "_ancestor_viewer.html" if htmx else "ancestors.html"
    response = render(request, f"phylobot/{template}", context)
    if selection is not None and htmx:
        response["HX-Push-Url"] = context["canonical_url"]
    return response


MAX_EXTANT_COLUMNS = 10


def extant_columns(query: Any, names: list[str]) -> list[str]:
    """The extant sequences shown beside the ancestor (?extant=...&extant=...), after
    the + (add_extant: the next sequence not yet shown) or × (remove_extant=<index>;
    never the first) buttons. Unknown names are dropped; at least one is shown."""
    chosen = [name for name in query.getlist("extant") if name in names] or names[:1]
    if "add_extant" in query and len(chosen) < MAX_EXTANT_COLUMNS:
        chosen.append(next((n for n in names if n not in chosen), names[0]))
    remove = query.get("remove_extant", "")
    if remove.isdigit() and 0 < int(remove) < len(chosen):
        del chosen[int(remove)]
    return chosen[:MAX_EXTANT_COLUMNS]


def combination_menus(
    selection: AncestorSelection,
) -> tuple[dict[str, str], list[tuple[str, str]]]:
    """The alignment menu ({key: label}) and the chosen alignment's model menu
    ([(combination key, label)]), shared by the ancestor and tree viewers."""
    combos = sorted(
        (c for c in selection.summary["combinations"]
         if f"{c['key']}/{selection.requires}" in selection.keys),
        key=results.combination_order(selection.summary.get("primary", "")),
    )  # fmt: skip
    several_trees = len({c["tree"] for c in combos}) > 1

    def model_label(c: dict[str, Any]) -> str:
        fitted = c.get("fit", {}).get("model", "")
        label = (
            f"ModelFinder → {fitted}"
            if catalog.is_modelfinder(c["model"]) and fitted
            else (c.get("model_label", c["model"]))
        )
        return f"{c.get('tree_label', c['tree'])} · {label}" if several_trees else label

    aligners = {c["aligner"]: c.get("aligner_label", c["aligner"]) for c in combos}
    aligner = selection.combination["aligner"]
    return aligners, [(c["key"], model_label(c)) for c in combos if c["aligner"] == aligner]


# supports.tsv columns, in the order shown ("SH-aLRT/UFBoot", as IQ-TREE writes them).
SUPPORT_MEASURES = treeview.SUPPORT_MEASURES


def support_values(selection: AncestorSelection) -> dict[str, dict[str, float]]:
    """{node: {measure: support %}} for the chosen combination (supports.tsv); measures
    that did not run, and nodes without a value, are left out."""
    if f"{selection.combination['key']}/supports.tsv" not in selection.keys:
        return {}
    values: dict[str, dict[str, float]] = {}
    for row in results.read_tsv(selection.read("supports.tsv")):
        for measure, _ in treeview.SUPPORT_MEASURES:
            try:
                values.setdefault(row["node"], {})[measure] = float(row.get(measure) or "")
            except ValueError:
                continue
    return {node: v for node, v in values.items() if v}


def agreement(selection: AncestorSelection, newick: str) -> dict[str, tuple[int, int]]:
    """{node: (k, n)}: how many of the job's n other alignment x model combinations of
    the same tree source have the branch above each node of `newick` (as drawn). Branches
    are compared as unrooted splits, so neither tree's rooting matters (like bootstrap
    support); {} when there is no other combination to compare with."""
    combo = selection.combination
    others = [
        c["key"]
        for c in selection.summary.get("combinations", [])
        if c["key"] != combo["key"]
        and c["tree"] == combo["tree"]
        and f"{c['key']}/tree.nwk" in selection.keys
    ]
    if not others:
        return {}
    drawn = trees.load_tree(newick)
    everyone = frozenset(leaf.taxon.label for leaf in drawn.leaf_node_iter())
    split_sets = []
    for key in others:
        other = trees.load_tree(selection.read("tree.nwk", key))
        if frozenset(leaf.taxon.label for leaf in other.leaf_node_iter()) != everyone:
            continue  # different sequences (e.g. a provided alignment): not comparable
        split_sets.append(
            {trees.canonical_split(clade, everyone) for _, clade in trees.internal_splits(other)}
        )
    if not split_sets:
        return {}
    found = {}
    for label, clade in trees.internal_splits(drawn):
        # a branch with one sequence (or none) on one side is in every tree: no label
        if label and 1 < len(clade) < len(everyone) - 1:
            split = trees.canonical_split(clade, everyone)
            found[label] = (sum(split in s for s in split_sets), len(split_sets))
    return found


def free_ratio_rates(selection: AncestorSelection) -> dict[frozenset[str], paths.BranchRate]:
    """codeml's free-ratio dN/dS per branch for the chosen combination, if it ran."""
    rates_file = f"codeml/{codeml.FREE_RATIOS}/branches.tsv"
    if f"{selection.combination['key']}/{rates_file}" not in selection.keys:
        return {}
    return paths.branch_rates(selection.read(rates_file))


def tree_display(
    query: Any, supports: dict[str, dict[str, float]], rates: dict[Any, Any]
) -> tuple[treeview.Display, dict[str, Any]]:
    """The tree viewer's branch display from its query (colour, width: "", a support
    measure or "dnds"; branch_palette: ancestors.PALETTES; labels: all, a measure or
    none; marks: numbers or stars; dnds: show or hide), offering only what this tree
    has; and the menus for the template."""
    names = dict(treeview.SUPPORT_MEASURES)
    measures = [m for m, _ in treeview.SUPPORT_MEASURES if any(m in v for v in supports.values())]
    styles = [("", "None")] + [
        (m, "Agreement across combos" if m == treeview.AGREE else f"{names[m]} support")
        for m in measures
    ]
    if rates:
        styles.append((treeview.DNDS, "dN/dS (codeml free ratios)"))
    allowed = {key for key, _ in styles}
    colour = query.get("colour", "") if query.get("colour", "") in allowed else ""
    width = query.get("width", "") if query.get("width", "") in allowed else ""
    label_choices = [("all", "All")] if len(measures) > 1 else []
    label_choices += [
        (m, ("Agreement (k of n)" if m == treeview.AGREE else names[m])
            if len(measures) > 1 else "Show")
        for m in measures
    ]  # fmt: skip
    label_choices += [("none", "None")] if measures else []
    statistics = [m for m in measures if m != treeview.AGREE]
    # By default the tree's own statistics; agreement is shown when chosen (or "All").
    default = "all" if not statistics or len(statistics) == len(measures) else (
        statistics[0] if len(statistics) == 1 else "stats"
    )  # fmt: skip
    if default == "stats":
        label_choices.insert(0, ("stats", " / ".join(names[m] for m in statistics)))
    labels = query.get("labels", default)
    if labels not in {key for key, _ in label_choices}:
        labels = default
    if labels == "all":
        shown = tuple(measures)
    elif labels == "stats":
        shown = tuple(statistics)
    else:
        shown = () if labels == "none" else (labels,)
    stars = query.get("marks") == "stars"
    dnds_labels = bool(rates) and query.get("dnds") == "show"
    palette = query.get("branch_palette", "")
    if palette not in ancestors.PALETTES:
        palette = ancestors.DEFAULT_PALETTE
    display = treeview.Display(colour, width, shown, stars, dnds_labels, palette)
    resolved = {"colour": colour, "width": width, "labels": labels,
                "marks": "stars" if stars else "numbers"}  # fmt: skip
    if colour:
        resolved["branch_palette"] = palette
    if rates:
        resolved["dnds"] = "show" if dnds_labels else "hide"
    menus = {
        "styles": styles, "colour": colour, "width": width,
        "label_choices": label_choices, "labels": labels, "stars": stars,
        "has_rates": bool(rates), "dnds_labels": dnds_labels,
        "palettes": [(k, label) for k, (label, _) in ancestors.PALETTES.items()],
        "palette": palette,
        "query": {k: v for k, v in resolved.items() if v},
    }  # fmt: skip
    return display, menus


def node_support(selection: AncestorSelection) -> tuple[dict[str, str], list[str]]:
    """({node: support text}, names of the measures shown) for the chosen combination
    (supports.tsv). Only measures that ran are shown; with two, the text is 'SH/UF'.
    Jobs from before SH-aLRT have a ufboot column only."""
    if f"{selection.combination['key']}/supports.tsv" not in selection.keys:
        return {}, []
    rows = results.read_tsv(selection.read("supports.tsv"))
    shown = [(c, name) for c, name in SUPPORT_MEASURES if any(row.get(c) for row in rows)]
    text = {
        row["node"]: "/".join(row.get(c) or "-" for c, _ in shown)
        for row in rows
        if any(row.get(c) for c, _ in shown)
    }
    return text, [name for _, name in shown]


@dataclass(frozen=True)
class PropertySource:
    """One file the property panel reads: the job's own PROPERTIES output, or a later
    run's (a properties run owns its results, so the panel is
    assembled from several prefixes rather than written back into this job)."""

    prefix: str
    summary: dict[str, Any]
    run: Job | None  # None for the job's own measurement

    @property
    def note(self) -> str:
        if self.run is None:
            return ""
        when = self.run.finished_at or self.run.created_at
        return f"measured later, on {when:%Y-%m-%d}"


def property_sources(
    selection: AncestorSelection,
    combination: str,
    node: str,
    keys: set[str],
    runs: list[Job],
) -> list[PropertySource]:
    """Every stored measurement of this ancestor and gap version, the job's own first
    and later runs after it (so a property measured twice shows the newer number).

    Read once per request: the panel (what was measured) and the offer (what was not)
    are answers from the same files.
    """
    version = selection.version
    out = []
    if property_views.has(keys, combination):
        prefix = selection.results.prefix(combination)  # an inherited one: the parent's
        if body := property_views.read(prefix, combination, node, version):
            out.append(PropertySource(prefix, body, None))
    holder = str(selection.results.owner(combination).id)
    for run in reversed(runs):
        of = run.params.get("properties_of") or {}
        if run.status != Job.Status.SUCCEEDED or of.get("gaps") != version:
            continue
        if properties_submit.measured_in(run, combination) != holder:
            continue  # a relative's run, on files this job does not read
        if node not in (of.get("combinations") or {}).get(combination, []):
            continue
        prefix = f"{run.storage_prefix}results/"
        if body := property_views.read(prefix, combination, node, version):
            out.append(PropertySource(prefix, body, run))
    return out


def property_data(
    selection: AncestorSelection, combination: str, node: str, sources: list[PropertySource]
) -> tuple[dict[str, Any] | None, dict[str, list[float]], dict[str, str]]:
    """One ancestor's measurements, however many runs they came from: the merged summary,
    each property's draws read from the file **its own number** came from, and a note of
    where each number came from. The panel and the comparison page both go through here,
    so neither can pair a histogram with another run's summary."""
    summary = property_views.merge([source.summary for source in sources])
    if not summary:
        return None, {}, {}
    measured = [k for k in summary.get("properties", {}) if k in properties.REGISTRY]
    histograms: dict[str, list[float]] = {}
    provenance: dict[str, str] = {}
    for source in reversed(sources):  # newest first: it wins in `merge` too
        left = [k for k in measured if k not in histograms]
        histograms |= property_views.read_histograms(
            source.prefix, combination, node, selection.version, left
        )
        for key in source.summary.get("properties", {}):
            provenance.setdefault(key, source.note)
    return summary, histograms, provenance


def property_panel(
    selection: AncestorSelection,
    combination: str,
    node: str,
    sources: list[PropertySource],
    altall_in: tuple[str, list[list[tuple[str, float]]], int] | None,
    threshold: float,
    coding: bool,
) -> property_views.Panel | None:
    """The ancestor viewer's property posteriors for the chosen node and gap version.

    Reads each property's draws (the histogram cannot be rebinned from a summary), and
    computes **f(AltAll)** here rather than reading it: AltAll depends on the viewer's
    threshold slider, which the pipeline never saw, so a stored value would be wrong as
    soon as the reader moved it.
    """
    summary, histograms, provenance = property_data(selection, combination, node, sources)
    if not summary:
        return None
    altall_values: dict[str, float] = {}
    # Not for codon views: there `altall_in` is the codon row, and the property functions
    # take amino acids — handing them codons would read "ATGAAA" as six residues and
    # produce a confident, meaningless number. The draws are amino-acid based either way
    # (EXPORT sums codon posteriors into amino acids), so there is nothing to compare to.
    if altall_in is not None and not coding:
        template, dists, width = altall_in
        try:
            sequence = library.alt_all(template, dists, threshold, width).sequence
        except ValueError:  # a shape mismatch is the viewer's problem, not the panel's
            sequence = ""
        if sequence:
            # The same masked template the draws came from, so f(AltAll) is comparable
            # with the histogram rather than with a differently gapped sequence.
            # Pure-Python properties only: the predictors need torch, which the web
            # tier does not carry, and a model call per page view would cost seconds.
            altall_values = {
                key: properties.REGISTRY[key].fn(sequence)
                for key in summary.get("properties", {})
                if key in properties.REGISTRY and properties.REGISTRY[key].image == properties.TOOLS
            }
    return property_views.of(
        summary,
        histograms,
        altall_values,
        about_page.PROPERTY_METHODS,
        f"{threshold:g}" if not coding else "",
        provenance,
    )


def property_redraw(
    selection: AncestorSelection,
    combination: str,
    node: str,
    sources: list[PropertySource],
    query: Any,
) -> dict[str, Any] | None:
    """The "draw again" control: which properties it would re-measure, and with what seed.

    A property posterior is a sample, so two runs of the same plan give slightly
    different numbers. The seed is shown because that is what makes the difference
    reproducible rather than mysterious, and it is the next attempt's, computed from the
    same plan key `properties_submit` will build — the control and the run must agree
    about which attempt this is.
    """
    summary = property_views.merge([source.summary for source in sources])
    measured = properties_run.measured(summary)
    if not measured or not properties_submit.measurable(selection.version):
        return None
    keys = [key for key in properties.REGISTRY if key in measured]
    return {
        "labels": [properties.REGISTRY[key].label for key in keys],
        "draws": properties.draws_in_range(query.get("draws")),
        "min_draws": properties.MIN_DRAWS,
        "max_draws": properties.MAX_DRAWS,
        # Every measurement this ancestor already has moves the seed on, so "draw again"
        # draws again: `sources` is exactly that list (the job's own and each later run).
        "seed": properties_submit.next_seed(selection.job, len(sources)),
    }


def property_offer(
    selection: AncestorSelection,
    combination: str,
    node: str,
    summary: dict[str, Any] | None,
    query: Any,
    length: int,
    runs: list[Job],
) -> dict[str, Any] | None:
    """What could still be measured for this ancestor.

    None when the registry holds nothing this ancestor lacks — the ordinary state of a
    job run today, and the reason the offer is built from the *difference* between the
    registry and the stored summary rather than from a flag: when B6 adds a property,
    every finished job starts offering it without anything being migrated.

    The scope counts come from the trees, so "in every combination that recovered it" is
    the same clade matching the viewer's "in the other combinations" section uses
    (`phylobot.analogs`, cached per tree file). The pipeline is later handed this explicit
    list and never repeats the matching (§18.5).
    """
    # Only a version a run could actually measure: PRANK's and GRASP's carry no
    # posteriors to draw from, and codeml's live under their own folder, which
    # PROPERTIES_ONLY does not stage. The offer and `properties_submit` ask the same
    # question, or the panel would show a button that refuses.
    if not properties_submit.measurable(selection.version):
        return None
    offered = properties_run.missing(summary)
    # The reader's own number, kept in the URL with the rest of the checklist so the
    # cost below it follows the field.
    draws = properties.draws_in_range(query.get("draws"))
    placed = analog_views.trees(selection)
    clade = analog_views.by_clade(selection)
    found = analogs.find(placed, combination, node, clade) if combination in placed else []
    nodes = {key: sorted(places) for key, places in placed.items()}
    # A combination whose rooted tree is missing is absent from `placed`; the viewer still
    # has its node list, so "every ancestor in this combination" stays answerable.
    nodes.setdefault(combination, sorted(selection.ancestors))
    counts = {
        key: len(
            properties_run.expand(
                key, combination, node, {a.combination: a.node for a in found}, nodes
            )
        )
        for key in properties_run.SCOPES
    }
    scope = properties_run.scope_of(query.get(properties_run.SCOPE_FIELD))
    asked = query.getlist(properties_run.FIELD) if properties_run.FIELD in query else None
    chosen = properties_run.chosen(asked, offered)
    targets = counts[scope]
    return {
        "rows": properties_run.offers(summary, targets, draws, length),
        "min_draws": properties.MIN_DRAWS,
        "max_draws": properties.MAX_DRAWS,
        "chosen": chosen,
        "field": properties_run.FIELD,
        "scope_field": properties_run.SCOPE_FIELD,
        "scope": scope,
        "scopes": [(key, label, counts[key]) for key, label in properties_run.SCOPES.items()],
        "targets": targets,
        "draws": draws,
        "cost": properties_run.cost_text(properties_run.total_cost(chosen, targets, draws, length)),
        "measured": sorted(properties_run.measured(summary)),
        # Runs already under way on this job's ancestors: the answer to "did my click
        # do anything?", which is otherwise invisible until the run finishes.
        "pending": [run for run in runs if run.status in (Job.Status.QUEUED, Job.Status.RUNNING)],
    }


def altall_inputs(
    selection: AncestorSelection,
    node: str,
    coding: bool,
    dists: list[list[tuple[str, float]]] | None = None,
) -> tuple[str, list[list[tuple[str, float]]], int] | None:
    """(aligned row, per-column distributions, column width) for AltAll, or None when the
    chosen reconstruction cannot offer one. The viewer and the downloads both go through
    this, so the count on the page and the file the user gets can never disagree.
    `dists` reuses the viewer's, which it has already read."""
    if not ancestors.has_posteriors(selection.version):
        return None
    if not coding:
        return (
            selection.ancestors[node],
            (
                dists
                if dists is not None
                else ancestors.distributions(
                    selection.read_node("ancestors_posteriors.tsv", node), node
                )
            ),
            1,
        )
    coding_file = ancestors.coding_file(selection.version)
    if f"{selection.combination['key']}/{coding_file}" not in selection.keys:
        return None
    row = ancestors.parse_fasta(selection.read(coding_file)).get(node)
    if row is None:
        return None
    return (
        row,
        (
            dists
            if dists is not None
            else ancestors.codon_distributions(
                selection.read_node(ancestors.CODON_POSTERIORS, node), node
            )
        ),
        3,
    )


@dataclass(frozen=True)
class Calibration:
    """The viewer's calibrated probabilities for one node: per column, its covariates and
    temperature, with the coefficients that made them (`phylobot.calibration`)."""

    coefficients: calibration.Coefficients
    covariates: list[dict[str, float]]
    temperatures: list[float]
    ancestor: str  # the node's row in the version shown: T is summarised over its residues

    @property
    def summary(self) -> calibration.Summary | None:
        return calibration.summary(
            [t for t, char in zip(self.temperatures, self.ancestor, strict=True) if char != "-"]
        )


def calibration_for(
    selection: AncestorSelection, node: str, alignment: dict[str, str]
) -> Calibration | None:
    """The node's temperatures under the coefficients in force (Settings), or None when
    they cannot be computed (the node is not in the combination's tree). Covariates come
    from the combination's own files: the alignment, tree.nwk (unrooted, so they do not
    depend on the outgroup) and site_rates.tsv."""
    keys, combo = selection.keys, selection.combination["key"]
    if f"{combo}/tree.nwk" not in keys:
        return None
    rates = selection.read("site_rates.tsv") if f"{combo}/site_rates.tsv" in keys else None
    rows = calibration.covariates(list(alignment.values()), selection.read("tree.nwk"), node, rates)
    if rows is None:
        return None
    # The viewer's own coefficients: jobs are seen only by their owner (Job.for_user).
    coefficients = site_settings.calibration_coefficients(selection.job.owner)
    temps = calibration.temperatures(coefficients, rows)
    return Calibration(coefficients, rows, temps, selection.ancestors[node])


def calibration_choice(query: Any, selection: AncestorSelection, states: str) -> str:
    """ "calibrated" only where it applies: IQ-TREE's (or codeml's) amino-acid posteriors.
    Codon probabilities stay raw (coefficients are fitted on amino acids)."""
    wanted = query.get("pp") == calibration.CALIBRATED
    applies = ancestors.has_posteriors(selection.version) and states == "amino"
    return calibration.CALIBRATED if wanted and applies else calibration.RAW


def ancestor_context(selection: AncestorSelection, query: Any) -> dict[str, Any]:
    combo, keys = selection.combination, selection.keys
    aligners, models = combination_menus(selection)
    support, _ = node_support(selection)
    called = ancestor_names.for_selection(selection)
    nodes = ancestors.node_choices(selection.tree(), selection.node_names, support, called)
    node = query.get("node", "")
    # Node names are per tree: coming from another combination, keep the same ancestor
    # (`phylobot.analogs`: its clade, or without an outgroup its partition).
    came_from = query.get("shown", "")
    if came_from != combo["key"] and f"{came_from}/ancestors.fasta" in keys:
        placed = analog_views.trees(selection)
        same = [
            a.node
            for a in analogs.find(placed, came_from, node, analog_views.by_clade(selection))
            if a.combination == combo["key"]
        ]
        node = same[0] if same else ""
    if node not in selection.node_names:
        node = nodes[0].name
    alignment = ancestors.parse_fasta(selection.read("alignment.fasta"))
    extants = extant_columns(query, list(alignment))
    palette = query.get("palette", "")
    if palette not in ancestors.PALETTES:
        palette = ancestors.DEFAULT_PALETTE
    lines = query.get("lines", "")
    if lines not in ancestors.LINE_STYLES:
        lines = ancestors.DEFAULT_LINES
    # States: amino acids, or codons for codon reconstructions (IQ-TREE's codon models and
    # codeml; their masked versions share the original's codon posteriors).
    codon_file = ancestors.version_folder(selection.version) + ancestors.CODON_POSTERIORS
    codons_available = (
        ancestors.has_posteriors(selection.version) and f"{combo['key']}/{codon_file}" in keys
    )
    states = "codons" if codons_available and query.get("states") == "codons" else "amino"
    notes = None
    if states == "codons":
        dists = ancestors.codon_distributions(
            selection.read_node(ancestors.CODON_POSTERIORS, node), node
        )
        code = codons.CODES.get(int(selection.job.params.get("genetic_code") or 1))
        notes = (code or codons.STANDARD).table
    else:
        dists = ancestors.distributions(selection.read_node("ancestors_posteriors.tsv", node), node)
    # Calibrated probabilities replace the raw ones for everything below that reads `dists`
    # (tiles, support plots, AltAll), so the page never mixes the two.
    pp = calibration_choice(query, selection, states)
    calibrated = None
    if pp == calibration.CALIBRATED:
        calibrated = calibration_for(selection, node, alignment)
    if calibrated is None:
        pp = calibration.RAW
    else:
        dists = calibration.apply(dists, calibrated.temperatures)
    # Only ml gap placement gives probabilities for gaps (gaps-ml.tsv).
    scored_gaps = selection.version == "ml" and f"{combo['key']}/gaps-ml.tsv" in keys
    p_gap = None
    if scored_gaps:
        p_gap = ancestors.gap_probabilities(selection.read("gaps-ml.tsv"), node)
    rows = ancestors.site_rows(
        selection.version,
        selection.ancestors[node],
        dists,
        [alignment[n] for n in extants],
        palette,
        p_gap,
        notes,
    )
    if calibrated is not None:
        rows = [
            replace(row, temperature=t)
            for row, t in zip(rows, calibrated.temperatures, strict=True)
        ]
    base = {"combination": combo["key"], "node": node, "gaps": selection.version}
    if states == "codons":
        base["states"] = "codons"
    if pp == calibration.CALIBRATED:
        base["pp"] = pp
    download = reverse("ancestor_download", args=[selection.job.id])
    has_posteriors = ancestors.has_posteriors(selection.version)
    altall = library.parse_threshold(query.get("altall"))
    # Every unmasked column's best alternative, threshold-free, so the page can report the
    # count at any threshold without going back to the server as the slider moves.
    altall_in = altall_inputs(selection, node, states == "codons", dists)
    ancestor_length = len(selection.ancestors[node].replace("-", ""))
    # One read for both the panel (what was measured) and the offer (what was not),
    # including any later properties run on this job's ancestors.
    property_runs = properties_submit.runs_for(selection.job)
    property_stored = property_sources(selection, combo["key"], node, keys, property_runs)
    altall_alternatives = (
        library.alternatives(altall_in[0], altall_in[1], altall_in[2]) if altall_in else []
    )

    def extra(what: str) -> dict[str, str]:
        """Only the AltAll downloads carry its threshold; the rest keep the links they had."""
        return {"altall": f"{altall:g}"} if what.startswith("altall") else {}

    view = base | {
        "aligner": combo["aligner"], "extant": extants, "palette": palette, "lines": lines,
        "altall": f"{altall:g}",
    }  # fmt: skip
    viewer = reverse("ancestor_viewer", args=[selection.job.id])
    return {
        "canonical_url": f"{viewer}?{urlencode(view, doseq=True)}",
        "tree_url": reverse("tree_viewer", args=[selection.job.id])
        + "?"
        + urlencode({"aligner": combo["aligner"], "combination": combo["key"], "node": node}),
        "excel_url": f"{download}?{urlencode(view | {'what': 'xlsx'}, doseq=True)}",
        "properties_csv_url": reverse("ancestor_properties_csv", args=[selection.job.id])
        + "?"
        + urlencode(view, doseq=True),
        # Only worth offering where there is something to compare against.
        "compare_url": (
            reverse("job_properties", args=[selection.job.id]) + "?" + urlencode(base)
            if len(analog_views.combinations(selection)) > 1
            else ""
        ),
        "aligners": aligners,
        "aligner": combo["aligner"],
        "aligner_label": aligners[combo["aligner"]],
        "models": models,
        "combination": combo["key"],
        "nodes": nodes,
        "node": node,
        "called": called.get(node),
        # Naming needs the node's clade from the job's ancestor map.
        "nameable": combo["key"] in analog_views.trees(selection),
        "description": description_of(selection, node),
        "name_rules": ancestor_names.RULES,
        "versions": [(v, ancestors.version_label(v)) for v in selection.versions],
        "version": selection.version,
        "version_label": ancestors.version_label(selection.version),
        "has_posteriors": has_posteriors,
        "codons_available": codons_available,
        "states": states,
        "scored_gaps": scored_gaps,
        # Raw (IQ-TREE's) or calibrated probabilities; the menu only where they apply.
        "pp": pp,
        "calibration_available": has_posteriors and states == "amino",
        "calibration": calibrated,
        "covariates": calibration.COVARIATES,
        "t_bounds": (calibration.T_MIN, calibration.T_MAX),
        "plots": ancestors.node_plots(
            selection.version,
            selection.ancestors[node],
            dists,
            palette,
            lines,
            "codon" if states == "codons" else "residue",
        ),  # fmt: skip
        "unit": "codon" if states == "codons" else "residue",
        # Properties were measured per node and gap version by PROPERTIES, so the panel
        # follows the viewer's own menus with one small read. Absent for jobs that did
        # not measure them (the box unticked, or run before the feature).
        "properties": property_panel(
            selection, combo["key"], node, property_stored, altall_in, altall, states == "codons"
        ),
        # Draw again, with a new seed, for the properties already measured (§18). Absent
        # when nothing is measured here, since there would be nothing to redraw.
        "property_redraw": property_redraw(selection, combo["key"], node, property_stored, query),
        # What the registry holds and this ancestor does not: the offer to measure it.
        # Absent when there is nothing left to offer.
        "property_offer": property_offer(
            selection,
            combo["key"],
            node,
            property_views.merge([source.summary for source in property_stored]),
            query,
            ancestor_length,
            property_runs,
        ),
        "line_styles": [(k, label) for k, (label, _) in ancestors.LINE_STYLES.items()],
        "lines": lines,
        "file_stem": re.sub(
            r"[^A-Za-z0-9._+-]",
            "_",
            f"{ancestor_names.fasta_id(node, called)}.{combo['key'].replace('/', '_')}"
            f".{selection.version}",
        ),
        "sequence_names": list(alignment),
        "extants": extants,
        "can_add_extant": len(extants) < MAX_EXTANT_COLUMNS,
        "palettes": [(k, label) for k, (label, _) in ancestors.PALETTES.items()],
        "palette": palette,
        "legend": [(pp, ancestors.colour(pp, palette)) for pp in (0.2, 0.4, 0.6, 0.8, 0.9, 1.0)],
        "rows": rows,
        "tile_columns": range(max((len(r.tiles) for r in rows), default=1)),
        "length": ancestor_length,
        "downloads": [
            (
                CODON_DOWNLOADS.get(what, label) if states == "codons" else label,
                f"{download}?{urlencode(base | extra(what) | {'what': what})}",
                DOWNLOAD_NOTES.get(what, ""),
            )
            for what, label in ANCESTOR_DOWNLOADS.items()
            if (what != "gaps" or scored_gaps)
            and (not what.startswith("altall") or altall_in is not None)
        ],
        "altall_threshold": f"{altall:g}",
        "altall_available": altall_in is not None,
        "altall_count": library.count_at(altall_alternatives, altall),
        "altall_unit": "codon" if states == "codons" else "site",
        # Sorted, so the slider's count is a binary search in the browser rather than a
        # scan of every site on each move.
        "altall_pps": sorted((round(s.pp, 5) for s in altall_alternatives), reverse=True),
        # The alternatives' probabilities as a histogram, so the threshold slider has
        # something to cut: the bars right of the line are the sites AltAll switches.
        "altall_plot": ancestors.altall_plot([s.pp for s in altall_alternatives], palette),
        "altall_min": library.MIN_THRESHOLD,
        "altall_max": library.MAX_THRESHOLD,
        # The same ancestor in the job's other combinations (`phylobot.analogs`).
        "comparison": analog_views.comparison(selection, node, query, palette, altall, extants[0]),
        # Duplications at it, or paralogs it mixes in the reconciled tree (`reconcile`).
        "paralogy": paralogy_context(selection, node),
    }


# Per-node downloads from the viewer: what -> link text.
ANCESTOR_DOWNLOADS = {
    "sequence": "ML sequence (FASTA)",
    "aligned": "aligned ML sequence (FASTA)",
    "altall": "AltAll sequence (FASTA)",
    "altall_aligned": "aligned AltAll sequence (FASTA)",
    "altall_sites": "AltAll's ambiguous sites (TSV)",
    "posteriors": "posterior probabilities of every state (TSV)",
    "states": "most probable state per site (TSV)",
    "gaps": "gap probabilities (TSV)",  # ml gap placement only
}
#: What each download is, for the viewer's table. Keyed as ANCESTOR_DOWNLOADS; a key
#: without an entry simply gets no description, so adding a download never breaks the
#: page (a test keeps them in step).
DOWNLOAD_NOTES = {
    "sequence": "The reconstructed ancestor itself, gaps removed: the sequence to order "
    "or to put in a figure.",
    "aligned": "The same sequence on the alignment's columns, gaps kept, so it lines up "
    "with the extant sequences and with ancestors from this alignment.",
    "altall": "The AltAll control at the threshold set under Support: every ambiguous site "
    "switched to its next-most-probable state. A result that holds for both this and the "
    "ML sequence does not rest on the uncertain sites.",
    "altall_aligned": "AltAll on the alignment's columns, for the same comparisons as "
    "the aligned ML sequence.",
    "altall_sites": "Only the sites AltAll switches, with both states and their "
    "probabilities: the evidence behind the control.",
    "posteriors": "Every state's posterior probability at every site — the full "
    "distribution this page draws, for your own analysis.",
    "states": "One row per site: the most probable state and its probability. The "
    "compact form, when the full distribution is more than you need.",
    "gaps": "Per site, the probability that this ancestor had a gap, from the "
    "gap/residue reconstruction (ML gap placement only).",
}
CODON_DOWNLOADS = {  # the same, with the viewer showing codons
    "sequence": "ML coding sequence (FASTA)",
    "aligned": "aligned ML coding sequence (FASTA)",
    "posteriors": "posterior probabilities of codons, 0.01 or more (TSV)",
    "states": "most probable codon per site (TSV)",
}


def property_candidates(
    selection: AncestorSelection,
    node: str,
    shown: set[str] | None,
    top: int | None = None,
) -> list[property_compare.Candidate]:
    """This ancestor and its analogs in every other combination, with their property
    draws. Matching is `phylobot.analogs` (clade, or partition without an outgroup), the
    same machinery the tree viewer's popovers use, on the same gap version throughout:
    comparing a masked ancestor with an unmasked one would confound gap placement with
    the method choice this page is about.

    A combination that lacks the clade, or that measured no properties, contributes
    nothing — `missing_combinations` reports those separately, because an absent ancestor
    is a result, not a gap. With `top`, only each alignment's
    `top` best by BIC are read, and the combination shown in the viewer.
    """
    here = selection.combination["key"]
    combos = {
        c["key"]: c
        for c in analog_views.combinations(selection)
        if c["key"] == here or property_compare.within_top(c, top)
    }
    placed = analog_views.trees(selection)
    found = analogs.find(placed, here, node, analog_views.by_clade(selection)) if placed else []
    several = len({c["tree"] for c in combos.values()}) > 1
    members = [(here, node)] + [(a.combination, a.node) for a in found]
    # Properties measured after the job finished live under their own run (§18.3), so
    # each candidate is read the same way the viewer's panel reads this one.
    runs = properties_submit.runs_for(selection.job)
    out = []
    for key, name in members:
        if key not in combos:
            continue
        sources = property_sources(selection, key, name, selection.keys, runs)
        summary, values, _where = property_data(selection, key, name, sources)
        if not summary or not values:
            continue
        out.append(
            property_compare.Candidate(
                key=key,
                node=name,
                label=analog_views.label(combos[key], several),
                aligner=combos[key]["aligner"],
                shown=shown is None or key in shown,
                values=values,
                summary=summary,
            )
        )
    return out


@require_GET
def job_properties(request: HttpRequest, job_id: str) -> HttpResponse:
    """Compare one ancestor's properties across every combination that recovered it.
    One GET form, so every view has a shareable URL."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    selection = select_ancestor(job, request.GET)
    context: dict[str, Any] = {"job": job, "selection": selection}
    if selection is not None:
        node = request.GET.get("node", "")
        if node not in selection.node_names:
            node = selection.node_names[0] if selection.node_names else ""
        asked = set(request.GET.getlist("show"))
        top = property_compare.top_of(request.GET.get("top", property_compare.DEFAULT_TOP))
        candidates = property_candidates(selection, node, asked or None, top)
        # Means are drawn unless explicitly switched off, so a URL that never mentions
        # them (a link from elsewhere) still shows them. The form always submits a
        # value — see the hidden field in the template.
        show_means = request.GET.get("means", "on") != "off"
        palette = property_compare.palette_of(request.GET.get("palette"))
        keys = [k for k in properties.REGISTRY if any(k in c.values for c in candidates)]
        panels = [
            p
            for k in keys
            if (
                p := property_compare.panel(
                    k,
                    properties.REGISTRY[k].label,
                    properties.REGISTRY[k].units,
                    candidates,
                    show_means=show_means,
                )
            )
            is not None
        ]
        every = analog_views.combinations(selection)
        here = selection.combination["key"]
        combos = {
            c["key"]: c for c in every if c["key"] == here or property_compare.within_top(c, top)
        }
        several = len({c["tree"] for c in every}) > 1
        have = {c.key for c in candidates}
        per_alignment = Counter((c["aligner"], c.get("data", "amino")) for c in every)
        most = max(per_alignment.values(), default=0)
        aligners, models = combination_menus(selection)
        support, _ = node_support(selection)
        view = {
            "aligner": selection.combination["aligner"],
            "combination": selection.combination["key"],
            "node": node,
            "gaps": selection.version,
        }
        context |= {
            "node": node,
            "show_means": show_means,
            "palette": palette,
            "palettes": [
                (k, label, note) for k, (label, note) in property_compare.PALETTES.items()
            ],
            "nodes": ancestors.node_choices(
                selection.tree(),
                selection.node_names,
                support,
                ancestor_names.for_selection(selection),
            ),  # fmt: skip
            "combination": selection.combination["key"],
            "aligner": selection.combination["aligner"],
            "aligners": aligners,
            "models": models,
            "versions": [(v, ancestors.version_label(v)) for v in selection.versions],
            "version": selection.version,
            "version_label": ancestors.version_label(selection.version),
            "candidates": candidates,
            # What the menu shows: "All" when N covers every alignment's combinations.
            "top": property_compare.ALL if top is None or top >= most else top,
            "top_choices": [n for n in property_compare.TOP_CHOICES if n < most],
            "top_default": property_compare.DEFAULT_TOP,
            "per_alignment": most,
            "alignment_groups": len(per_alignment),
            "compared": len(combos),
            "combination_count": len(every),
            # Each panel with what the template needs beside it, rather than a dict the
            # template would have to index (Django has no such filter, and adding a
            # generic one to reach into a mapping hides where the data came from).
            "panels": [
                {"panel": p, "method": about_page.PROPERTY_METHODS.get(p.key)} for p in panels
            ],
            "missing": [analog_views.label(c, several) for k, c in combos.items() if k not in have],
            "viewer_url": reverse("ancestor_viewer", args=[job.id]) + "?" + urlencode(view),
            "canonical_url": reverse("job_properties", args=[job.id])
            + "?"
            + urlencode(
                view
                | {
                    "show": [c.key for c in candidates if c.shown],
                    "top": property_compare.ALL if top is None else top,
                    "palette": palette,
                    "means": "on" if show_means else "off",
                },
                doseq=True,
            ),
        }
    return render(request, "phylobot/properties.html", context)


@require_POST
def ancestor_name(request: HttpRequest, job_id: str) -> HttpResponse:
    """Name an ancestor (POST combination, gaps, node, ancestor_name; or remove), then
    back to it. The clade comes from the job's own ancestor map, never the request, and
    the name is checked here whatever the browser allowed (phylobot.ancestor_names)."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    selection = select_ancestor(job, request.POST)
    node = request.POST.get("node", "")
    if selection is None or node not in selection.node_names:
        raise Http404
    place = analog_views.trees(selection).get(selection.combination["key"], {}).get(node)
    if place is None:
        raise Http404
    clade = ancestor_names.match_key(place, analog_views.by_clade(selection))
    back = {
        "aligner": selection.combination["aligner"],
        "combination": selection.combination["key"],
        "node": node,
        "gaps": selection.version,
    }
    url = reverse("ancestor_viewer", args=[job.id]) + "?" + urlencode(back)
    # Named from the tree viewer's popover: back to the tree as it was. Only this job's
    # tree viewer, by path and query: never another site or page.
    tree = reverse("tree_viewer", args=[job.id])
    asked = str(request.POST.get("return_to", ""))
    if asked == tree or asked.startswith(tree + "?"):
        url = asked
    if request.POST.get("remove"):
        ancestor_names.remove(job, clade)
        messages.success(request, f"{node}'s name removed.")
        return redirect(url)
    # The ancestor viewer's Name and Description form posts both; the tree popover the name.
    if "description" in request.POST:
        ancestor_names.save_note(job, clade, str(place.parts), str(request.POST["description"]))
    name = str(request.POST.get("ancestor_name", "")).strip()
    current = AncestorName.objects.filter(job=ancestor_names.holder(job), clade=clade).first()
    if name == (current.name if current else ""):
        if "description" in request.POST:
            messages.success(request, f"{node}'s description saved.")
        return redirect(url)
    if not name:  # the name cleared: as "Remove the name"
        ancestor_names.remove(job, clade)
        messages.success(request, f"{node}'s name removed.")
        return redirect(url)
    sequences = set(ancestors.parse_fasta(selection.read("alignment.fasta")))
    error = ancestor_names.problem(name, sequences) or ancestor_names.save(
        job, clade, str(place.parts), name
    )
    if error:
        messages.error(request, f"Not named: {error}")
    else:
        messages.success(request, f"{node} is now {name}, in every combination that has it.")
    return redirect(url)


def description_of(selection: AncestorSelection, node: str) -> str:
    """The user's description of `node` (by its clade, so shared across combinations)."""
    place = analog_views.trees(selection).get(selection.combination["key"], {}).get(node)
    if place is None:
        return ""
    return ancestor_names.note_for(selection.job, place, analog_views.by_clade(selection))


@require_POST
def collapse_clade(request: HttpRequest, job_id: str) -> HttpResponse:
    """Fold a clade into a triangle in the tree viewer, restyle it, or expand it again
    (POST combination, aligner, node, action, and for styling label/fill/border_colour/
    border_width). Like ancestor_name, the clade comes from the job's own ancestor map,
    never the request, so it follows the clade onto every combination (phylobot.clade_collapse)."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    selection = select_ancestor(job, request.POST)
    node = request.POST.get("node", "")
    if selection is None or node not in selection.node_names:
        raise Http404
    place = analog_views.trees(selection).get(selection.combination["key"], {}).get(node)
    if place is None:
        raise Http404
    clade = ancestor_names.match_key(place, analog_views.by_clade(selection))
    # Back to the tree exactly as it was, but only this job's tree viewer (never elsewhere).
    tree = reverse("tree_viewer", args=[job.id])
    asked = str(request.POST.get("return_to", ""))
    url = asked if asked == tree or asked.startswith(tree + "?") else tree
    if request.POST.get("action") == "expand":
        clade_collapse.remove(job, clade)
        messages.success(request, f"{node} expanded.")
    else:
        clade_collapse.save(
            job, clade, str(place.parts),
            request.POST.get("label", ""),
            request.POST.get("fill", ""),
            request.POST.get("border_colour", ""),
            request.POST.get("border_width", ""),
        )  # fmt: skip
        messages.success(request, f"{node} collapsed, in every combination that has it.")
    return redirect(url)


@require_GET
def ancestor_properties_csv(request: HttpRequest, job_id: str) -> HttpResponse:
    """The property panel as it is shown, as CSV: one row per property, the columns the
    page displays plus the provenance needed to read them later (node, combination, gap
    version, draws, the AltAll threshold f(AltAll) was taken at).

    Built from the same `property_views.Panel` the page renders, so the file and the page
    can never disagree — the ancestor viewer's Excel download works the same way.
    """
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    selection = select_ancestor(job, request.GET)
    if selection is None:
        raise Http404
    node = request.GET.get("node", "")
    if node not in selection.node_names:
        raise Http404
    combination = selection.combination["key"]
    states_codons = request.GET.get("states") == "codons"
    altall = library.parse_threshold(request.GET.get("altall"))
    panel = property_panel(
        selection,
        combination,
        node,
        property_sources(
            selection, combination, node, selection.keys, properties_submit.runs_for(selection.job)
        ),
        altall_inputs(selection, node, states_codons),
        altall,
        states_codons,
    )
    if panel is None:
        raise Http404
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow([
        "property", "label", "units", "mean", "sd",
        f"low_{panel.level}pc", f"high_{panel.level}pc",
        "f_ML", "shift_sd", "f_AltAll", "altall_threshold",
        "draws", "node", "combination", "gaps", "sampling",
    ])  # fmt: skip
    for row in panel.rows:
        writer.writerow([
            row.key, row.label, row.units, row.mean, row.sd, row.low, row.high,
            row.ml, row.shift, row.altall, panel.altall_threshold,
            # Per property: runs measured at different draw counts can sit in one panel.
            row.draws, node, combination, selection.version, panel.sampling,
        ])  # fmt: skip
    stem = re.sub(r"[^A-Za-z0-9._+-]", "_", f"{node}.{combination.replace('/', '_')}")
    response = HttpResponse(out.getvalue(), content_type="text/csv; charset=utf-8")
    response["Content-Disposition"] = f'attachment; filename="{stem}.properties.csv"'
    return response


@require_GET
def ancestor_download(request: HttpRequest, job_id: str) -> HttpResponse:
    """One node's sequence or site table, e.g. ?combination=...&node=Node3&what=posteriors."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    what = request.GET.get("what", "")
    selection = select_ancestor(job, request.GET)
    if selection is not None and what == "xlsx":
        return site_table_workbook(selection, request.GET)
    node = request.GET.get("node", "")
    if selection is None or what not in ANCESTOR_DOWNLOADS or node not in selection.node_names:
        raise Http404
    tag = selection.combination["key"].replace("/", "_")
    # A named ancestor is its name in the files (phylobot.ancestor_names), its node after it.
    called = ancestor_names.for_selection(selection)
    head = ancestor_names.fasta_id(node, called)
    record = head if head == node else f"{head} {node}"
    coding = request.GET.get("states") == "codons" and ancestors.has_posteriors(selection.version)
    combo_key = selection.combination["key"]
    # Calibrated downloads are the page's probabilities (`phylobot.calibration`), never
    # IQ-TREE's file relabelled; the sequences themselves are the same either way.
    calibrated = None
    if calibration_choice(request.GET, selection, "codons" if coding else "amino") == "calibrated":
        alignment = ancestors.parse_fasta(selection.read("alignment.fasta"))
        calibrated = calibration_for(selection, node, alignment)
    stamp = f" pp=calibrated ({calibrated.coefficients.source})" if calibrated else ""
    if coding and what in ("sequence", "aligned"):
        path = f"{combo_key}/{ancestors.coding_file(selection.version)}"
        if path not in selection.keys:
            raise Http404
        row = ancestors.parse_fasta(selection.read(ancestors.coding_file(selection.version)))[node]
        sequence = row.replace("-", "") if what == "sequence" else row
        body = f">{record} {combo_key} gaps={selection.version} coding DNA\n{sequence}\n"
        filename = f"{head}.{tag}.{selection.version}.cds.fasta"
    elif what in ("sequence", "aligned"):
        row = selection.ancestors[node]
        sequence = row.replace("-", "") if what == "sequence" else row
        body = f">{record} {selection.combination['key']} gaps={selection.version}\n{sequence}\n"
        filename = f"{head}.{tag}.{selection.version}.fasta"
    elif what.startswith("altall"):
        # AltAll follows the viewer's States menu, as the ML downloads do: amino acids, or
        # codons for a codon reconstruction. Versions without posteriors (PRANK, GRASP)
        # have no alternatives to offer.
        inputs = altall_inputs(selection, node, coding)
        if inputs is None:
            raise Http404
        template, dists, width = inputs
        if calibrated is not None:
            dists = calibration.apply(dists, calibrated.temperatures)
        threshold = library.parse_threshold(request.GET.get("altall"))
        design = library.alt_all(template, dists, threshold, width)
        kind = "cds." if coding else ""
        stem = f"{head}.{tag}.{selection.version}.altall{threshold:g}"
        if calibrated is not None:
            stem += ".calibrated"
        if what == "altall_sites":
            body, filename = design.sites_tsv(), f"{stem}.{kind}sites.tsv"
        else:
            sequence = design.sequence if what == "altall_aligned" else design.unaligned
            body = (
                f">{head}_AltAll{record[len(head) :]} {selection.combination['key']} "
                f"gaps={selection.version} "
                f"altall threshold={threshold:g} switched={design.ambiguous}{stamp}"
                f"{' coding DNA' if coding else ''}\n{sequence}\n"
            )
            filename = f"{stem}.{kind}fasta"
    else:
        source = {"posteriors": "ancestors_posteriors.tsv", "states": "ancestors_pp.tsv",
                  "gaps": "gaps-ml.tsv"}[what]  # fmt: skip
        if coding and what != "gaps":
            source = {"posteriors": ancestors.CODON_POSTERIORS,
                      "states": "ancestors_codon_pp.tsv"}[what]  # fmt: skip
        name = source
        if what != "gaps":  # codeml's versions have their own posteriors
            source = ancestors.version_folder(selection.version) + source
        if f"{selection.combination['key']}/{source}" not in selection.keys:
            raise Http404
        # Posteriors are read by node through their byte index; the small tables whole.
        if calibrated is not None and what in ("posteriors", "states"):
            raw = selection.read_node("ancestors_posteriors.tsv", node)
            body = calibrated_table(raw, node, calibrated, what)
            filename = f"{head}.{tag}.{what}.calibrated.tsv"
        else:
            whole = what != "posteriors"
            text = selection.read(source) if whole else selection.read_node(name, node)
            body = "\n".join(ancestors.posterior_lines(text, node)) + "\n"
            filename = f"{head}.{tag}.{what}.tsv"
    response = HttpResponse(body, content_type="text/plain; charset=utf-8")
    filename = re.sub(r"[^A-Za-z0-9._+-]", "_", filename)
    response["Content-Disposition"] = f'attachment; filename="{filename}"'
    return response


def calibrated_table(raw: str, node: str, calibrated: Calibration, what: str) -> str:
    """ancestors_posteriors.tsv's rows for `node` (what="posteriors": every state, in its
    column order) or ancestors_pp.tsv's (what="states": the most probable state), with
    the calibrated probabilities in place of IQ-TREE's, and a T column."""
    header, *_ = ancestors.posterior_lines(raw, node)
    states = header.split("\t")[2:]
    dists = calibration.apply(ancestors.distributions(raw, node), calibrated.temperatures)
    if what == "states":
        lines = ["node\tsite\tstate\tpp\ttemperature"]
        lines += [
            f"{node}\t{site}\t{dist[0][0]}\t{dist[0][1]:.5f}\t{t:.4f}"
            for site, (dist, t) in enumerate(zip(dists, calibrated.temperatures, strict=True), 1)
        ]
    else:
        lines = ["\t".join(["node", "site", *states, "temperature"])]
        for site, (dist, t) in enumerate(zip(dists, calibrated.temperatures, strict=True), 1):
            p = dict(dist)
            lines.append("\t".join([node, str(site), *(f"{p[x]:.5f}" for x in states), f"{t:.4f}"]))
    return "\n".join(lines) + "\n"


def site_table_workbook(selection: AncestorSelection, query: Any) -> HttpResponse:
    """The viewer's site table as it is shown (same node, sequences and palette), as an
    Excel workbook with the tiles' colours (phylobot.sites_excel)."""
    from phylobot import sites_excel  # openpyxl: only needed here

    context = ancestor_context(selection, query)
    combo = selection.combination
    model = dict(context["models"])[combo["key"]]
    node = next(c.label for c in context["nodes"] if c.name == context["node"])
    plots = context["plots"]
    about = [
        ("Job", f"{selection.job.name or '(unnamed)'} ({selection.job.id})"),
        ("Alignment", context["aligner_label"]),
        ("Model", model),
        ("Ancestral node", node),
        ("Reconstruction / gap masking", context["version_label"]),
        (
            "States",
            "codons (with the amino acid each encodes)"
            if context["states"] == "codons"
            else "amino acids",
        ),
        (
            "Posterior probabilities",
            f"calibrated (T per site; coefficients: {context['calibration'].coefficients.source})"
            if context["calibration"]
            else "raw (IQ-TREE's)",
        ),
        ("Extant sequences", ", ".join(context["extants"])),
        ("Palette", dict(context["palettes"])[context["palette"]]),
        ("States shown", f"probability {ancestors.MIN_SHOWN} or more (the TSV downloads have all)"),
        ("Exported", timezone.now().strftime("%Y-%m-%d %H:%M UTC")),
    ]
    if plots:
        summary = f"mean posterior probability {plots.mean:.3f}, standard deviation "
        summary += f"{plots.sd:.3f}, over {plots.sites} sites"
        about.insert(5, ("ML residues", summary))
    response = HttpResponse(
        sites_excel.workbook(context, about),
        content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
    response["Content-Disposition"] = f'attachment; filename="{context["file_stem"]}.sites.xlsx"'
    return response


def tree_viewer(request: HttpRequest, job_id: str) -> HttpResponse:
    """One combination's tree, drawn left to right (phylobot.treeview); ancestral nodes
    link to the ancestor viewer. Controls re-request this view, as in the ancestor viewer."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    selection = select_ancestor(job, request.GET, requires="tree.nwk")
    context: dict[str, Any] = {"job": job, "selection": selection}
    if selection is not None:
        combo = selection.combination
        aligners, models = combination_menus(selection)
        mode = request.GET.get("lengths", "")
        if mode not in treeview.BRANCH_MODES:
            mode = treeview.DEFAULT_MODE
        # "all" / "named" (only ancestors the user named) / "none": node name labels;
        # support labels keep their own control. Not asked for: every name on a small
        # tree, only the user's own on a crowded one (decided once the tree is drawn).
        node_labels = request.GET.get("node_labels", "")
        labels_asked = node_labels in ("all", "named", "none")
        rooted = f"{combo['key']}/tree.rooted.nwk" in selection.keys
        tree_file = "tree.rooted.nwk" if rooted else "tree.nwk"
        tag = combo["key"].replace("/", "_")
        newick = selection.tree()
        # The path is only highlighted here; its page (path_views) shows what changed.
        traced = path_context(selection, request.GET, newick, with_table=False)
        supports = support_values(selection)
        agree = agreement(selection, newick)
        for node, (k, n) in agree.items():
            supports.setdefault(node, {})[treeview.AGREE] = round(100 * k / n, 1)
        rates = free_ratio_rates(selection)
        display, display_menus = tree_display(request.GET, supports, rates)
        bs_panel = branch_site_views.panel(job, combo, selection.keys, newick, request.GET)
        recon = reconciliation_context(selection)
        show_species = bool(recon) and request.GET.get("species") == "show"
        folded = clade_collapse.for_selection(selection)
        drawing = treeview.layout(
            newick,
            mode,
            ancestors=set(selection.node_names),
            outgroup=set(job.params.get("outgroup") or []),
            path=traced["path"],
            supports=supports,
            rates=rates,
            display=display,
            foreground=bs_panel.foreground,
            events=recon["events"] if recon else None,
            losses=recon["losses"] if recon else None,
            agreement=agree,
            tip_notes=recon["species"] if recon and show_species else None,
            collapsed=folded,
        )
        spacing = treeview.slider(request.GET.get("spacing"), treeview.TIP_SPACING)
        drawing = treeview.stretch(drawing, spacing)
        if not labels_asked:
            crowded = len(drawing.nodes) > treeview.CROWDED_NODES
            node_labels = "named" if crowded else "all"
        names = dict(treeview.SUPPORT_MEASURES)
        support_measures = [
            names[m] for m, _ in treeview.SUPPORT_MEASURES
            if m in display.labels and m != treeview.AGREE
            and any(m in v for v in supports.values())
        ]  # fmt: skip
        agreement_shown = treeview.AGREE in display.labels and bool(agree)
        agreement_total = next(iter(agree.values()))[1] if agree else 0
        highlight = request.GET.get("node", "")
        ancestor_url = reverse("ancestor_viewer", args=[job.id])
        here = {"aligner": combo["aligner"], "combination": combo["key"]}
        called = ancestor_names.for_selection(selection)
        drawing.nodes = [
            replace(
                n,
                href=f"{ancestor_url}?{urlencode(here | {'node': n.name})}",
                called=called[n.name].name if n.name in called else "",
                called_alike=called[n.name].alike if n.name in called else True,
            )
            if n.link
            else n
            for n in drawing.nodes
        ]
        # A folded clade's ancestor is not drawn as a node, so carry its name onto the
        # triangle (shown at the apex dot, where that ancestor sits).
        drawing.collapsed = [
            replace(cl, called=called[cl.name].name, called_alike=called[cl.name].alike)
            if cl.name in called
            else cl
            for cl in drawing.collapsed
        ]
        linked = {n.name for n in drawing.nodes if n.link}
        # Each ancestral node's popover: the same ancestor in the other combinations.
        popovers = [
            replace(p, called=called[p.node].name) if p.node in called else p
            for p in analog_views.popovers(selection, linked)
        ]
        view = {"aligner": combo["aligner"], "combination": combo["key"], "lengths": mode}
        view |= display_menus["query"]
        view |= {"node_labels": node_labels} if labels_asked else {}
        view |= {} if spacing == 1.0 else {"spacing": f"{spacing:g}"}
        view |= {"species": "show"} if show_species else {}
        # What the path page keeps to come back to this tree as it is drawn now.
        back = urlencode({k: v for k, v in view.items() if k in saved_paths.TREE_VIEW_KEYS})
        view |= traced.pop("query")
        view |= bs_panel.query
        tree_url = reverse("tree_viewer", args=[job.id])
        path_page = {"aligner": combo["aligner"], "combination": combo["key"],
                     "gaps": selection.version, "path_from": traced["path_from"],
                     "path_to": traced["path_to"]} | ({"back": back} if back else {})  # fmt: skip
        context |= traced
        context |= {
            "path_page_url": f"{reverse('path_viewer', args=[job.id])}?{urlencode(path_page)}",
            "path_viewer_url": reverse("path_viewer", args=[job.id]),
            "saved_paths": saved_paths.listing(job),
            "aligners": aligners,
            "aligner": combo["aligner"],
            "models": models,
            "combination": combo["key"],
            "name_rules": ancestor_names.RULES,
            "modes": list(treeview.BRANCH_MODES.items()),
            "mode": mode,
            "node_labels": node_labels,
            "labels_auto": not labels_asked,
            "spacing": f"{spacing:g}",
            "spacing_bounds": treeview.TIP_SPACING,
            "drawing": drawing,
            "bs_panel": bs_panel,
            "bs_scopes": list(branch_site.SCOPES.items()),
            "support_measures": support_measures,
            "agreement_shown": agreement_shown,
            "show_species": show_species,
            "agreement_total": agreement_total,
            "support_stars": display.stars,
            "stars": treeview.STARS,
            "display": display_menus,
            "highlight": highlight if highlight in linked else "",
            "popovers": popovers,
            "clade_styles": folded,
            "clade_fills": clade_collapse.FILL_SWATCHES,
            "clade_borders": clade_collapse.BORDER_SWATCHES,
            "clade_border_widths": clade_collapse.BORDER_WIDTHS,
            "analog_note": analog_views.match_note(analog_views.by_clade(selection)),
            "library_highlight_nodes": [highlight] if highlight in linked else [],
            "rooted": rooted,
            "reconciliation": recon,
            "outgroup": job.params.get("outgroup") or [],
            "newick_url": reverse("job_file", args=[job.id, f"{combo['key']}/{tree_file}"]),
            "file_stem": re.sub(r"[^A-Za-z0-9._+-]", "_", f"tree.{tag}.{mode}"),
            "canonical_url": f"{tree_url}?{urlencode(view, doseq=True)}",
        }
    htmx = bool(request.headers.get("HX-Request"))
    one_off = "compare_set" in request.GET or "path_clear" in request.GET
    if selection is not None and one_off and not htmx:
        return redirect(context["canonical_url"])
    template = "_tree_viewer.html" if htmx else "tree.html"
    response = render(request, f"phylobot/{template}", context)
    if selection is not None and htmx:
        response["HX-Push-Url"] = context["canonical_url"]
    return response


def reconciliation_context(selection: AncestorSelection) -> dict[str, Any] | None:
    """The combination's reconciliation with the species tree, if the job ran one
    (phylobot.reconcile): node events and branch losses for the drawing, its totals and
    file links for the page."""
    combo = selection.combination["key"]
    files = ("reconciliation.tsv", "reconciliation.nwk", "reconciliation.json")
    if any(f"{combo}/{name}" not in selection.keys for name in files):
        return None
    events, losses = reconcile.marks(
        selection.read("reconciliation.tsv"), selection.read("reconciliation.nwk")
    )
    job = selection.job
    return {
        "events": {name: treeview.Event(*event) for name, event in events.items()},
        "weak": sum(1 for _, _, overlap in events.values() if overlap and overlap.weak),
        "species": reconcile.leaf_species(selection.read("reconciliation.tsv")),
        "losses": losses,
        "summary": json.loads(selection.read("reconciliation.json")),
        "tsv_url": reverse("job_file", args=[job.id, f"{combo}/reconciliation.tsv"]),
        "nwk_url": reverse("job_file", args=[job.id, f"{combo}/reconciliation.nwk"]),
        "species_tree_url": (
            reverse("job_file", args=[job.id, "species/species_tree.nwk"])
            if "species/species_tree.nwk" in selection.keys
            else ""
        ),
    }


def other_tree(selection: AncestorSelection) -> str | None:
    """The tree to compare this combination's with (same alignment and model): GeneRax's
    reconciled tree for a sequence-only tree, and IQ-TREE's (else RAxML-NG's, else a
    provided one) for GeneRax's. None if the job has no such reconciled combination."""
    combo, keys = selection.combination, selection.keys
    names = ("reconciliation.tsv", "reconciliation.nwk")
    if combo.get("tree") == "generax":
        sources = ["iqtree", "raxml"] + sorted(
            {c["tree"] for c in selection.summary.get("combinations", [])
             if c["tree"].startswith("provided")}
        )  # fmt: skip
    else:
        sources = ["generax"]
    for source in sources:
        key = f"{combo['aligner']}/{source}/{combo['model']}"
        if all(f"{key}/{name}" in keys for name in names):
            return key
    return None


def au_values(selection: AncestorSelection) -> dict[str, float]:
    """{tree source: AU p} from the alignment and model's topology test, if it ran."""
    combo = selection.combination
    path = f"{combo['aligner']}/topology_tests/{combo['model']}.tsv"
    if path not in selection.keys:
        return {}
    found = {}
    text = selection.read(f"{combo['model']}.tsv", f"{combo['aligner']}/topology_tests")
    for row in results.read_tsv(text):
        try:
            found[row["tree"]] = float(row["p_au"])
        except KeyError, ValueError:
            continue
    return found


def paralogy_context(selection: AncestorSelection, node: str) -> dict[str, Any] | None:
    """What reconciliation says about the viewed ancestor (`reconcile.paralogy`), and
    whether the job's other tree of the same alignment and model agrees: the R5
    validation found that either the sequence-only tree or GeneRax's can be the wrong
    one, so a disagreement is reported as such, with the AU test's view, rather than as
    a verdict. None when there is nothing to say."""
    combo, keys = selection.combination, selection.keys
    names = ("reconciliation.tsv", "reconciliation.nwk")
    if any(f"{combo['key']}/{name}" not in keys for name in names):
        return None
    own = (selection.read(names[0]), selection.read(names[1]))
    other = other_tree(selection)
    against = (selection.read(names[0], other), selection.read(names[1], other)) if other else None
    found = reconcile.paralogy(node, own, against)
    if found is None:
        return None
    overlap = reconcile.overlaps(*own).get(node) if found.duplication else None
    disagree = against is not None and not found.reconciled_node
    if not (found.duplication or found.reconciled_node or disagree):
        return None
    job = selection.job
    viewer, tree = reverse("ancestor_viewer", args=[job.id]), reverse("tree_viewer", args=[job.id])
    own_source = combo["tree"]
    other_source = other.split("/")[1] if other else ""
    label = {"generax": "GeneRax's reconciled tree", "iqtree": "IQ-TREE's tree",
             "raxml": "RAxML-NG's tree"}  # fmt: skip

    def ancestor_url(name: str) -> str:
        if not name or other is None:
            return ""
        return f"{viewer}?{urlencode({'combination': other, 'node': name})}"

    au = au_values(selection)
    own_p, other_p = au.get(own_source), au.get(other_source)
    verdict = ""
    if own_p is not None and other_p is not None:
        rejected = [s for s, p in ((own_source, own_p), (other_source, other_p)) if p < 0.05]
        if not rejected:
            verdict = "neither"
        elif len(rejected) == 1:
            verdict = "own" if rejected[0] == own_source else "other"
    shown = 6  # sequences named per group before "and N more"
    return {
        "found": found,
        "overlap": overlap,
        "disagree": disagree,
        "own_label": label.get(own_source, f"{own_source}'s tree"),
        "other_label": label.get(other_source, f"{other_source}'s tree"),
        "own_p": own_p,
        "other_p": other_p,
        "verdict": verdict,
        # the clade in GeneRax's rooting of this tree, as the comparison takes it
        "clade": sorted(reconcile.rooted_clades(own[1])[0].get(node, ())),
        "own_short": label.get(own_source, own_source).split("'")[0],
        "other_short": label.get(other_source, other_source).split("'")[0],
        "root": found.mixes == reconcile.ROOT,
        "reconciled_url": ancestor_url(found.reconciled_node),
        "duplication_url": (
            f"{tree}?{urlencode({'combination': other, 'node': found.mixes})}"
            if found.mixes and found.mixes != reconcile.ROOT
            else f"{tree}?{urlencode({'combination': other or combo['key']})}"
        ),
        "groups": [
            {
                "names": group[:shown],
                "more": max(0, len(group) - shown),
                "url": ancestor_url(n),
                "node": n,
            }
            for group, n in zip(found.groups, found.group_nodes, strict=True)
        ],  # fmt: skip
    }


def path_workbook(context: dict[str, Any]) -> HttpResponse:
    """The traced mutational path as an Excel workbook (phylobot.sites_excel)."""
    from phylobot import sites_excel  # openpyxl: only needed here

    job = context["job"]
    about = [
        ("Job", f"{job.name or '(unnamed)'} ({job.id})"),
        ("Alignment", context["aligners"][context["aligner"]]),
        ("Model", dict(context["models"])[context["combination"]]),
        ("Path", " → ".join(context["path"])),
        ("Gap/indel masking", dict(context["versions"])[context["version"]]),
        ("Columns", "only those that change" if context["only_changes"] else "all"),
        *(
            [
                ("Silent changes", "shown" if context["show_silent"] else "hidden"),
                ("Cells", "codons" if context["cells"] == "codons" else "amino acids"),
            ]
            if context["has_codons"]
            else []
        ),
        ("Palette", dict(context["palettes"])[context["palette"]]),
        *([("Compared with", context["compared_with"])] if context.get("compared_with") else []),
        ("Exported", timezone.now().strftime("%Y-%m-%d %H:%M UTC")),
    ]
    response = HttpResponse(
        sites_excel.path_workbook(
            context["path_table"], context["palette"], about, context.get("path_support")
        ),
        content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
    name = re.sub(r"[^A-Za-z0-9._+-]", "_", f"path.{context['path'][0]}-{context['path'][-1]}")
    response["Content-Disposition"] = f'attachment; filename="{name}.xlsx"'
    return response


def path_context(
    selection: AncestorSelection, query: Any, newick: str, with_table: bool = True
) -> dict[str, Any]:
    """A mutational path (phylobot.paths) checked against the tree: its two ends and,
    once traced, the route and (`with_table`) the path table. The tree viewer asks for
    the route alone, to highlight it; the path page (phylobot.path_views) for the table.
    `query` keys: path_from, path_to, trace, cols (changes|all), compare (extant names;
    compare_set=all|none sets them), palette, gaps; for codon reconstructions also silent
    (show|hide) and cells (amino|codons). The returned "query" holds the resolved values
    for the page's URL. path_clear (the Clear button) drops every path setting."""
    if query.get("path_clear"):
        query = QueryDict(mutable=True)
    called = ancestor_names.for_selection(selection)
    ancestors_ = ancestors.node_choices(newick, selection.node_names, {}, called)
    start = query.get("path_from", "")
    if start not in {n.name for n in ancestors_}:
        start = ""
    below = paths.descendants(newick, start) if start else []
    end = query.get("path_to", "")
    if end not in {n.name for n in below if n.extant or n.name in selection.node_names}:
        end = ""
    alignment = ancestors.parse_fasta(selection.read("alignment.fasta"))
    names = list(alignment)
    compare = [n for n in query.getlist("compare") if n in alignment]
    choice = query.get("compare_set", "")
    compare = names if choice == "all" else [] if choice == "none" else compare
    palette = query.get("palette", "")
    if palette not in ancestors.PALETTES:
        palette = ancestors.DEFAULT_PALETTE
    only_changes = query.get("cols", "changes") != "all"
    # Codon reconstructions: the version's coding ancestors and the extant codon alignment.
    combo_key = selection.combination["key"]
    coding = ancestors.coding_file(selection.version)
    has_codons = {f"{combo_key}/{coding}", f"{combo_key}/alignment.cds.fasta"} <= selection.keys
    show_silent = query.get("silent", "show") != "hide"
    cells = "codons" if has_codons and query.get("cells") == "codons" else "amino"
    traced = bool(query.get("trace")) and bool(start and end)
    route = paths.path(newick, start, end) if traced else []
    table = None
    rates = None  # codeml's free-ratio dN/dS per branch, when it ran
    if route and with_table:
        on_tree = [n for n in route if n not in alignment]
        sequences = {n: selection.ancestors[n] for n in on_tree} | alignment
        probabilities = {}
        if selection.version in (ancestors.ORIGINAL, *gaps.MASKING):
            probabilities = paths.ml_probabilities(
                selection.read_ancestor("ancestors_pp.tsv"), set(on_tree)
            )
        codon_rows = codon_pp = None
        if has_codons:
            codon_rows = ancestors.parse_fasta(selection.read("alignment.cds.fasta"))
            coded = ancestors.parse_fasta(selection.read(coding))
            codon_rows |= {n: coded[n] for n in on_tree}
            pp_file = ancestors.version_folder(selection.version) + "ancestors_codon_pp.tsv"
            if cells == "codons" and f"{combo_key}/{pp_file}" in selection.keys:
                codon_pp = paths.ml_probabilities(selection.read(pp_file), set(on_tree))
        rates = free_ratio_rates(selection) or None
        table = paths.trace(
            route, sequences, probabilities, set(alignment), compare, palette, only_changes,
            codon_rows, codon_pp, show_silent, cells, rates,
        )  # fmt: skip
        if called:  # branches between named ancestors say so ("Anc.GFP (Node4) → Node5")
            table.branches = [
                replace(
                    b,
                    parent=ancestor_names.label(b.parent, called),
                    child=ancestor_names.label(b.child, called),
                )
                for b in table.branches
            ]
    resolved: dict[str, Any] = {"gaps": selection.version, "palette": palette}
    if start:
        resolved["path_from"] = start
    if end:
        resolved["path_to"] = end
    if traced:
        resolved |= {"trace": "1", "cols": "changes" if only_changes else "all"}
        if has_codons:
            resolved |= {"silent": "show" if show_silent else "hide", "cells": cells}
    if compare:
        resolved["compare"] = compare
    return {
        "query": resolved,
        "path_starts": ancestors_,
        "path_from": start,
        "path_ends": [n for n in below if n.extant or n.name in selection.node_names],
        "path_to": end,
        "path": route,
        "path_table": table,
        "traced": traced,
        "only_changes": only_changes,
        "has_codons": has_codons,
        "path_rates": bool(rates),
        "show_silent": show_silent,
        "cells": cells,
        "sequence_names": names,
        "compare": compare,
        "palettes": [(k, label) for k, (label, _) in ancestors.PALETTES.items()],
        "palette": palette,
        "versions": [(v, ancestors.version_label(v)) for v in selection.versions],
        "version": selection.version,
        "legend": [(pp, ancestors.colour(pp, palette)) for pp in (0.2, 0.4, 0.6, 0.8, 0.9, 1.0)],
    }


@login_not_required
def healthz(request: HttpRequest) -> JsonResponse:
    """Report whether the web tier can reach Postgres and object storage."""
    checks: dict[str, str] = {}
    try:
        with connection.cursor() as cursor:
            cursor.execute("SELECT 1")
        checks["postgres"] = "ok"
    except Exception:
        logger.exception("healthz: postgres check failed")
        checks["postgres"] = "error"
    if storage.bucket_reachable():
        checks["storage"] = "ok"
    else:
        logger.error("healthz: object storage bucket not reachable")
        checks["storage"] = "error"
    healthy = all(status == "ok" for status in checks.values())
    return JsonResponse(checks, status=200 if healthy else 503)
