"""My alignments (`/alignments/`) and the alignment editor (`/alignments/<id>/`).

The editor is the alignment viewer's canvas (`alignment.js`) over one stored alignment,
with editing added: its edits live in the page's URL (`alignments.Edits`), so the page
is one GET form like the viewers and every edit is applied here from the stored rows.
Selections are made in the browser (`alignment_editor.js`) and sent as positions in the
view shown (`alignments.ONE_OFF`); the response's URL is the canonical one. "Save" makes
a new alignment; nothing here changes a stored one.
"""

import re
from typing import Any
from urllib.parse import urlencode

from django.contrib import messages
from django.contrib.auth.models import User
from django.core.cache import cache
from django.db.models import Prefetch
from django.http import Http404, HttpRequest, HttpResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.utils.timesince import timesince
from django.views.decorators.http import require_GET, require_POST

from phylobot import alignment_views, alignments, alignview
from phylobot.forms import AlignmentUploadForm
from phylobot.models import Alignment, Job, limits_for

JOBS_SHOWN = 50  # jobs listed under "From your jobs"
CACHE_SECONDS = 3600


def current_user(request: HttpRequest) -> User:
    """The signed-in user (local mode: the local user); middleware guarantees one."""
    assert isinstance(request.user, User)
    return request.user


def owned(request: HttpRequest, alignment_id: Any) -> Alignment:
    return get_object_or_404(
        Alignment.objects.for_user(request.user).select_related("job", "derived_from__job"),
        id=alignment_id,
    )


def file_stem(record: Alignment) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", record.title).strip("_")[:80] or "alignment"


@require_GET
def my_alignments(request: HttpRequest) -> HttpResponse:
    return page(request, AlignmentUploadForm(limits=limits_for(request.user)))


def page(request: HttpRequest, form: AlignmentUploadForm, status: int = 200) -> HttpResponse:
    """Saved alignments (uploads and edits: rename, delete). Each job's alignments
    (read-only: they belong to the job) follow from `jobs_section`, fetched by the page
    once it has loaded, so the saved ones show at once however many jobs there are."""
    saved = Alignment.objects.for_user(request.user).filter(job__isnull=True)
    context = {
        "form": form,
        "saved": saved.select_related("derived_from__job"),
        "limits": limits_for(request.user),
        "max_saved": alignments.MAX_SAVED,
    }
    return render(request, "phylobot/my_alignments.html", context, status=status)


@require_GET
def jobs_section(request: HttpRequest) -> HttpResponse:
    """My alignments' "From your jobs" (HTMX): registers older jobs' alignments first."""
    alignments.ensure_registered(request.user)
    from_jobs = (
        Job.objects.for_user(request.user)
        .filter(alignments__isnull=False)
        .distinct()
        .prefetch_related(
            Prefetch("alignments", queryset=Alignment.objects.order_by("created_at", "job_key")),
            "tags",
        )
        .order_by("-created_at")
    )
    jobs = list(from_jobs[:JOBS_SHOWN])
    context = {"jobs": jobs, "more_jobs": max(0, from_jobs.count() - len(jobs))}
    return render(request, "phylobot/_my_alignments_jobs.html", context)


@require_POST
def upload(request: HttpRequest) -> HttpResponse:
    form = AlignmentUploadForm(request.POST, request.FILES, limits=limits_for(request.user))
    if not form.is_valid():
        return page(request, form, status=400)
    if alignments.saved_count(request.user) >= alignments.MAX_SAVED:
        form.add_error(None, f"You have {alignments.MAX_SAVED} saved alignments, the most "
                             "allowed: delete some first.")  # fmt: skip
        return page(request, form, status=400)
    upload_file = form.cleaned_data.get("file")
    record = alignments.store(
        current_user(request), form.cleaned_data["name"], form.rows, form.sequence_type,
        form.cleaned_data["genetic_code"], origin=Alignment.Origin.UPLOAD,
        filename=alignments.clean_name(upload_file.name if upload_file else "")[:200],
    )  # fmt: skip
    return redirect("alignment_editor", alignment_id=record.id)


@require_POST
def upload_check(request: HttpRequest) -> HttpResponse:
    """The upload form's live check (HTMX), as soon as a file is chosen or text pasted:
    what the alignment holds, or why it would not be saved, and the choice for a name
    used for different rows. Nothing is saved."""
    form = AlignmentUploadForm(request.POST, request.FILES, limits=limits_for(request.user))
    empty = not request.FILES.get("file") and not request.POST.get("text", "").strip()
    valid = not empty and form.is_valid()
    context = {
        "empty": empty,
        "form": form,
        "valid": valid,
        "repeats": form.repeats,
        "rows": len(form.rows) if valid else 0,
        "columns": len(form.rows[0].row) if valid and form.rows else 0,
        "codon": form.sequence_type == "codon",
        "duplicates": request.POST.get("duplicates", ""),
    }
    return render(request, "phylobot/_alignment_upload_check.html", context)


@require_GET
def download(request: HttpRequest, alignment_id: str) -> HttpResponse:
    """The stored alignment as FASTA; `?states=amino` translates a codon alignment."""
    record = owned(request, alignment_id)
    rows = alignments.read(record)
    alignments.touch(record)
    stem = file_stem(record)
    if record.seqtype == "codon" and request.GET.get("states") == "amino":
        amino = alignments.amino_rows(rows, "codon", alignments.genetic_code(record))
        rows = [alignments.Row(r.name, r.description, a) for r, a in zip(rows, amino, strict=True)]
        stem += ".protein"
    return fasta_response(rows, f"{stem}.fasta")


def fasta_response(rows: list[alignments.Row], filename: str) -> HttpResponse:
    response = HttpResponse(alignments.to_fasta(rows), content_type="text/plain; charset=utf-8")
    response["Content-Disposition"] = f'attachment; filename="{filename}"'
    return response


@require_POST
def rename(request: HttpRequest, alignment_id: str) -> HttpResponse:
    record = owned(request, alignment_id)
    if record.from_job:
        raise Http404  # named by its job and aligner
    if name := alignments.clean_name(request.POST.get("name", "")):
        record.name = name
        record.save(update_fields=["name", "updated_at"])
    return redirect(back(request))


@require_POST
def delete(request: HttpRequest, alignment_id: str) -> HttpResponse:
    record = owned(request, alignment_id)
    if record.from_job:
        raise Http404  # deleted with its job
    alignments.delete(record)
    messages.info(request, f"Deleted {record.name}.")
    return redirect("my_alignments")


@require_POST
def duplicate(request: HttpRequest, alignment_id: str) -> HttpResponse:
    """A saved copy of a saved alignment, named "<name> (copy)"."""
    record = owned(request, alignment_id)
    if record.from_job:
        raise Http404  # "Edit a copy" is the job alignments' way
    if alignments.saved_count(request.user) >= alignments.MAX_SAVED:
        messages.error(request, f"You have {alignments.MAX_SAVED} saved alignments, the most "
                                "allowed: delete some first.")  # fmt: skip
        return redirect("my_alignments")
    rows = alignments.read(record)
    suffix = " (copy)"
    name = record.name[: alignments.MAX_NAME - len(suffix)] + suffix
    edits = alignments.provenance(rows, alignments.Edits(), record.seqtype, record.version)
    new = alignments.store(
        current_user(request),
        name,
        rows,
        record.seqtype,
        record.genetic_code,
        origin=Alignment.Origin.EDIT,
        derived_from=record,
        edits={**edits, "duplicate": True},
    )
    messages.success(request, f"Duplicated as “{new.name}”.")
    return redirect("my_alignments")


def back(request: HttpRequest) -> str:
    """Where a form asked to return to: a page of this site, else My alignments."""
    target = request.POST.get("next", "")
    return target if target.startswith("/alignments/") else reverse("my_alignments")


@require_POST
def save(request: HttpRequest, alignment_id: str) -> HttpResponse:
    """The alignment as edited (the posted edit state): over this one ("Save
    (overwrite)", saved alignments only), or as a new alignment ("Save a copy")."""
    record = owned(request, alignment_id)
    rows = alignments.read(record)
    sites = alignments.site_count(rows, record.seqtype)
    edits = alignments.Edits.from_query(request.POST, len(rows), sites)
    editor_url = f"{reverse('alignment_editor', args=[record.id])}?{urlencode(edits.query())}"
    if request.POST.get("mode") == "overwrite":
        if record.from_job:
            raise Http404  # the job depends on it
        if str(record.version) != request.POST.get("version", ""):
            # Positions name the rows the page showed; these have changed since.
            stale = "This alignment was overwritten after the page was opened: nothing was saved."
            messages.error(request, f"{stale} Here it is as it is now.")
            return redirect("alignment_editor", alignment_id=record.id)
        if not edits.changed:
            messages.info(request, "No edits: nothing to save.")
            return redirect(editor_url)
        alignments.overwrite(record, rows, edits)
        messages.success(request, f"Saved “{record.name}” (overwritten).")
        return redirect("alignment_editor", alignment_id=record.id)
    name = alignments.clean_name(request.POST.get("name", "")) or f"{record.title} (edited)"
    if alignments.saved_count(request.user) >= alignments.MAX_SAVED:
        messages.error(request, f"You have {alignments.MAX_SAVED} saved alignments, the most "
                                "allowed: delete some on My alignments first.")  # fmt: skip
        return redirect(editor_url)
    new = alignments.store(
        current_user(request), name, alignments.apply(rows, edits, record.seqtype),
        record.seqtype, record.genetic_code, origin=Alignment.Origin.EDIT,
        derived_from=record,
        edits=alignments.provenance(rows, edits, record.seqtype, record.version),
    )  # fmt: skip
    messages.success(request, f"Saved as “{new.name}”: it is in My alignments.")
    return redirect("alignment_editor", alignment_id=new.id)


PICKER_PAGE = 30


@require_GET
def picker(request: HttpRequest) -> HttpResponse:
    """The New Job page's "Load from My alignments" list (HTMX): the most recently used
    first, filtered as the user types (`q`), a page at a time as they scroll (`page`)."""
    if request.GET.get("page", "0") == "0":
        alignments.ensure_registered(request.user)
    search = request.GET.get("q", "")[:100]
    page = request.GET.get("page", "0")
    number = int(page) if page.isdigit() and len(page) < 6 else 0
    start = number * PICKER_PAGE
    found = list(alignments.recent(request.user, search)[start : start + PICKER_PAGE + 1])
    for a in found:  # "used 5 minutes ago", "job finished 2 days ago"
        when = a.last_used_at or a.used  # type: ignore[attr-defined]
        what = "used" if a.last_used_at else "job finished" if a.job_id else "saved"
        ago = timesince(when)
        a.ago = f"{what} just now" if ago.startswith("0\xa0") else f"{what} {ago} ago"  # type: ignore[attr-defined]
    context = {
        "items": found[:PICKER_PAGE],
        "next_page": number + 1 if len(found) > PICKER_PAGE else None,
        "q": search,
        "first": number == 0,
    }
    return render(request, "phylobot/_alignment_picker_rows.html", context)


@require_GET
def editor(request: HttpRequest, alignment_id: str) -> HttpResponse:
    record = owned(request, alignment_id)
    rows = alignments.read(record)
    sites = alignments.site_count(rows, record.seqtype)
    edits = alignments.Edits.from_query(request.GET, len(rows), sites)
    htmx = bool(request.headers.get("HX-Request"))
    note = ""
    if any(request.GET.get(k) for k in alignments.ONE_OFF):
        edits, note = alignments.edit(rows, edits, request.GET, record.seqtype)
        if note:
            messages.warning(request, note)
    if not htmx:  # opened (an edit's swap is not a new use)
        alignments.touch(record)
    context = editor_context(request, record, rows, edits) | {"htmx": htmx}
    export = request.GET.get("export", "")
    if export in ("fasta", "clustal", "tracks"):
        return exported(context, export)
    if not htmx and (any(k in request.GET for k in alignments.ONE_OFF)
                     or any(k in request.GET for k in alignment_views.ONE_OFF)):  # fmt: skip
        return redirect(context["canonical_url"])
    template = "_alignment_editor.html" if htmx else "alignment_editor.html"
    response = render(request, f"phylobot/{template}", context)
    if htmx:
        response["HX-Push-Url"] = context["canonical_url"]
    return response


def editor_context(
    request: HttpRequest, record: Alignment, rows: list[alignments.Row], edits: alignments.Edits
) -> dict[str, Any]:
    query = request.GET
    seqtype = record.seqtype
    codon = seqtype == "codon"
    code = alignments.genetic_code(record)
    stored_sites = alignments.site_count(rows, seqtype)
    shown = alignments.apply(rows, edits, seqtype)
    names = [r.name for r in shown]
    amino = alignments.amino_rows(shown, seqtype, code)
    sites = alignments.site_count(shown, seqtype)

    scheme_menu = alignment_views.SCHEMES | (alignment_views.CODON_SCHEMES if codon else {})
    scheme = query.get("scheme", "")
    if scheme not in scheme_menu:
        scheme = alignment_views.DEFAULT_SCHEME
    reference = query.get("reference", "")
    if reference not in names:
        reference = ""
    zoom = alignment_views.parse_zoom(query.get("zoom", ""))
    height = query.get("height", "")
    if height not in alignment_views.HEIGHTS:
        height = alignment_views.DEFAULT_HEIGHT

    offered = list(alignview.COMPOSITION)
    shown_tracks = alignment_views.track_list(query, offered)
    exporting = query.get("export") == "tracks"
    key = f"alignedit:v1:{record.id}:{record.version}:{edits.key}"
    computed: dict[str, alignview.Track] | None = cache.get(key)
    if computed is None:
        computed = alignview.composition_tracks(amino)
        cache.set(key, computed, CACHE_SECONDS)
    ref_amino = amino[names.index(reference)] if reference else ""
    profile = computed["hydropathy_profile"]
    profile.series = [alignview.Series(profile.label, alignview.profile(ref_amino or amino[0]))]
    tracks = [computed[k] for k in (offered if exporting else shown_tracks)]

    view: dict[str, Any] = edits.query() | {
        "scheme": scheme, "zoom": zoom, "height": height, "tracks_set": "1",
        "track": shown_tracks,
    }  # fmt: skip
    if reference:
        view["reference"] = reference
    here = reverse("alignment_editor", args=[record.id])
    canonical = f"{here}?{urlencode(view, doseq=True)}"
    kept = alignments.kept_sites(stored_sites, edits)
    gaps = len(alignments.gap_sites(shown, seqtype))
    empty = len(alignments.empty_rows(shown))
    data = {
        "names": names,
        "ancestor": [False] * len(names),
        "rows": [r.row for r in shown],
        "amino": amino if codon else [],
        "codon": codon,
        "extant": len(names),
        "reference": names.index(reference) if reference else -1,
        "scheme": scheme,
        "cell": zoom,
        "rowFromCell": True,
        "height": height,
        "tracks": [t.json() for t in tracks if t.key in shown_tracks],
        "ancestorViewer": "",
        "editable": True,
        "sort": edits.sort,
        "stateKey": f"{record.id}|{record.version}|{edits.key}",
        "origin": [s + 1 for s in kept] if edits.cut else [],
    }
    derived = record.derived_from
    if derived is not None and (derived.deleted_at or derived.job and derived.job.deleted_at):
        derived = None
    return {
        "record": record,
        "derived": derived,
        # The parent has been overwritten since this was made from it.
        "derived_older": derived is not None
        and record.edits.get("parent_version", 1) != derived.version,
        "edits": edits,
        "state": edits.query(),
        "changed": edits.changed,
        "canonical_url": canonical,
        "reset_url": f"{here}?{urlencode({k: v for k, v in view.items() if k not in edits.query()}, doseq=True)}",  # noqa: E501
        "export_url": canonical,
        "schemes": list(scheme_menu.items()),
        "scheme": scheme,
        "legend": alignment_views.legend(scheme),
        "references": names,
        "reference": reference,
        "zoom": zoom,
        "heights": list(alignment_views.HEIGHTS.items()),
        "height": height,
        "shown": [(k, alignview.TRACKS[k]) for k in shown_tracks],
        "addable": [
            (
                group,
                [
                    (k, alignview.TRACKS[k])
                    for k in offered
                    if k not in shown_tracks and alignview.TRACKS[k].group == group
                ],
            )
            for group in dict.fromkeys(
                alignview.TRACKS[k].group for k in offered if k not in shown_tracks
            )
        ],  # fmt: skip
        "can_add": len(shown_tracks) < alignment_views.MAX_TRACKS,
        "data": data,
        "sequences": len(names),
        "columns": sites,
        "stored_sequences": len(rows),
        "stored_columns": stored_sites,
        "removed_sequences": len(rows) - len(names),
        "removed_columns": stored_sites - sites,
        "gap_columns": gaps,
        "empty_rows": empty,
        "sorts": alignments.SORTS,
        "tracks": tracks,
        "download_rows": {r.name: r.row for r in shown},
        "shown_rows": shown,
        "states": "codons" if codon else "amino",
        "file_stem": file_stem(record) + (".edited" if edits.changed else ""),
        "default_name": f"{record.title} (edited)"[: alignments.MAX_NAME],
        "seqtypes": alignments.SEQTYPES,
    }


def exported(context: dict[str, Any], what: str) -> HttpResponse:
    """The alignment as edited (FASTA, with descriptions; or Clustal) or its tracks."""
    stem = context["file_stem"]
    if what == "fasta":
        return fasta_response(context["shown_rows"], f"{stem}.fasta")
    return alignment_views.download(context | {"ready": True}, what)
