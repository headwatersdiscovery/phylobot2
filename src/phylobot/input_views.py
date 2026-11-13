"""The New Job page's lists of uploaded files (HTMX): sequence files in step 2 and
alignment files in step 3, each with a tick box, what it holds and a remove button.

The files themselves stay in the page until Submit (submit.html's script keeps every
file the user chose in the form's file input; `forms.used_files` leaves out the
unticked ones), so these views only describe what was posted: per-file counts, sequence
files' name clashes (`phylobot.inputs.merge_sequences`), and whether the alignments hold
the same sequences, with the choice to use the ones they share.
"""

from typing import Any

from django.core.files.uploadedfile import UploadedFile
from django.http import HttpRequest, HttpResponse
from django.shortcuts import render
from django.views.decorators.http import require_POST

from phylobot import alignments, inputs, storage
from phylobot.forms import MAX_PROVIDED, chosen_alignments, is_uuid, sequence_names, used_files
from phylobot.models import limits_for
from phylobot.pipeline.validate import MIN_SEQUENCES

MAX_NAMES_SHOWN = 8  # missing names listed per alignment


def read_upload(upload: UploadedFile[Any], limit: int) -> bytes:
    """An upload's content (bounded), from the start however often it is read."""
    upload.seek(0)
    data: bytes = upload.read(limit)
    return data


def listed(request: HttpRequest, field: str, kind: str, aligned: bool) -> list[dict[str, Any]]:
    """Every posted file of a file list, ticked or not, with its stats."""
    limit = limits_for(request.user).max_upload_bytes
    used = {f.name for f in used_files(request.FILES.getlist(field), request.POST, kind)}
    rows, seen = [], set()
    for upload in request.FILES.getlist(field):
        if upload.name in seen:
            continue
        seen.add(upload.name)
        too_big = (upload.size or 0) > limit
        content = read_upload(upload, limit)
        problem = ""
        if aligned and not too_big:  # PHYLIP and "?" as the submitted job reads them
            try:
                content = alignments.normalise_upload(content)
            except alignments.AlignmentError as exc:
                problem = str(exc)
        stats, _ = inputs.read_stats(upload.name or "", content, aligned)
        if problem:
            stats.error = problem
        if too_big:
            stats.error = f"larger than your limit of {limit:,} bytes"
        rows.append(
            {"name": upload.name, "used": upload.name in used, "stats": stats, "content": content}
        )
    return rows


def merged_sequences(request: HttpRequest) -> inputs.Merged | None:
    """Step 2's sequences with repeated names resolved: the ticked files combined, or
    the pasted text; None without either."""
    limit = limits_for(request.user).max_upload_bytes
    files = used_files(request.FILES.getlist("sequences_file"), request.POST, "sequences")
    contents = [(f.name or "", read_upload(f, limit)) for f in files[: inputs.MAX_SEQUENCE_FILES]]
    text = request.POST.get("sequences_text", "")
    if not contents and text.strip():
        contents = [(inputs.PASTED, text.encode()[:limit])]
    if not contents:
        return None
    return inputs.merge_sequences(contents, request.POST.get("duplicate_names", ""))


def sequence_files_context(request: HttpRequest) -> dict[str, Any]:
    """`_sequence_files.html`'s context."""
    return {
        "rows": listed(request, "sequences_file", "sequences", aligned=False),
        "merged": merged_sequences(request),
        "duplicate_names": request.POST.get("duplicate_names", ""),
    }


@require_POST
def alignment_files(request: HttpRequest) -> HttpResponse:
    """Step 3's uploaded alignments, and whether every alignment (uploaded or chosen
    from My alignments) holds the same sequences as each other and as step 2's."""
    limit = limits_for(request.user).max_upload_bytes
    rows = listed(request, "alignment_files", "alignment", aligned=True)
    choice = request.POST.get("alignment_duplicates", "")
    sources: list[tuple[str, set[str]]] = []
    clashes, identical = [], []
    for row in rows:
        if row["used"] and not row["stats"].error:  # repeated names resolved as chosen
            deduped = inputs.dedupe_alignment(row["name"], row["content"], choice)
            clashes += deduped.clashes
            identical += deduped.identical
            sources.append((row["name"], set(sequence_names(deduped.fasta))))
    keys = [k for k in request.POST.getlist("saved_alignments") if is_uuid(k)]
    for record in chosen_alignments(request.user, keys)[:MAX_PROVIDED]:
        content = storage.read_bytes(record.storage_key)[:limit]
        sources.append((record.title, set(sequence_names(content))))
    merged = merged_sequences(request)
    if merged is not None and not merged.error:
        sources.append(("your sequences (step 2)", set(sequence_names(merged.fasta))))
    context: dict[str, Any] = {
        "rows": rows,
        "too_many": len(rows) > MAX_PROVIDED,
        "clashes": clashes,
        "identical": identical,
        "alignment_duplicates": choice,
    }
    names = [n for _, n in sources]
    if len(sources) > 1 and len({frozenset(n) for n in names}) > 1:
        everyone = set().union(*names)
        shared = inputs.shared_names(names)
        context["mismatch"] = {
            "sources": [
                {
                    "label": label,
                    "count": len(n),
                    "missing": sorted(everyone - n)[:MAX_NAMES_SHOWN],
                    "missing_count": len(everyone - n),
                }
                for label, n in sources
            ],
            "shared": len(shared),
            "left_out": len(everyone - shared),
            "enough": len(shared) >= MIN_SEQUENCES,
            "minimum": MIN_SEQUENCES,
        }
        context["shared_names"] = request.POST.get("shared_names") == "on"
    return render(request, "phylobot/_alignment_files.html", context)
