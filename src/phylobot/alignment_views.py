"""The alignment viewer (`/jobs/<id>/alignment/`).

One GET form, as the ancestor and tree viewers, so every view has a shareable URL. The
grid and tracks are drawn by `alignment.js` in a canvas from JSON this view embeds;
metrics come from `phylobot.alignview`. Every value taken from the query is checked
against the job's own alignments, combinations, sequence names and the track registry.
"""

import json
import re
from dataclasses import dataclass
from typing import Any
from urllib.parse import urlencode

from django.core.cache import cache
from django.http import Http404, HttpRequest, HttpResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse

from phylobot import (
    alignments,
    alignview,
    ancestors,
    catalog,
    codons,
    indel_views,
    joined,
    library_glyphs,
    manifest,
    results,
    storage,
)
from phylobot.models import Alignment, AlignmentTracks, Job

# Colour schemes: key -> menu label (the colours are in alignment.js; legends below).
SCHEMES = {
    "clustalx": "Clustal X",
    "differences": "Differences",
    "chemistry": "Chemistry",
    "hydropathy": "Hydropathy (Kyte-Doolittle)",
    "identity": "Percentage identity",
    "structure": "Secondary-structure propensity (Chou-Fasman)",
    "none": "None",
}
CODON_SCHEMES = {"synonymous": "Synonymous differences"}
DEFAULT_SCHEME = "clustalx"
# Zoom: a slider (as on the library page) from "fit" to library_glyphs.MAX_CELL pixels
# per residue; alignment.js derives each row's height from the width. The old menu's
# keys still arrive in saved links, so they map to the widths they stood for.
OLD_ZOOMS = {"overview": 3.0, "small": 6.0, "medium": 11.0, "large": 15.0}
DEFAULT_ZOOM = 11.0


def parse_zoom(text: str) -> float | str:
    """The viewer's zoom from request text: "fit", a cell width (clamped to the slider's
    range), or an old menu key; anything else is the default."""
    if text in OLD_ZOOMS:
        return OLD_ZOOMS[text]
    if text == "fit":
        return "fit"
    try:
        value = float(text)
    except ValueError:
        return DEFAULT_ZOOM
    return library_glyphs.parse_zoom(text) if value == value else DEFAULT_ZOOM


# Viewer height: fit the window (its bottom, and so the horizontal scrollbar, on screen
# once scrolled to), or a number of sequences.
HEIGHTS = {"fit": "Fit the window", "10": "10 sequences", "20": "20 sequences",
           "40": "40 sequences"}  # fmt: skip
DEFAULT_HEIGHT = "fit"
MAX_TRACKS = 12
EXPORTS = ("fasta", "clustal", "tracks")


@dataclass(frozen=True)
class Source:
    """One alignment the job has: its amino-acid file and, for coding DNA, its codon
    alignment (storage keys relative to the job's results/ or inputs/)."""

    key: str
    label: str
    amino: str  # results-relative path ("" for a codon-only input)
    codon: str = ""  # results-relative path, or "inputs/..." for a branch-site input


def job_keys(job: Job, inherited: bool = True) -> tuple[set[str], dict[str, Any]]:
    """The job's result files and summary; with `inherited`, joined with the results an
    extension reads from its parent (`phylobot.joined`)."""
    if not inherited:
        keys = joined.listing(job)
        summary = (
            json.loads(storage.read_bytes(f"{job.storage_prefix}results/summary.json"))
            if "summary.json" in keys
            else {}
        )
        return keys, summary
    found = joined.cached(job)
    return found.keys, found.summary


def sources(job: Job, keys: set[str], summary: dict[str, Any]) -> list[Source]:
    """The job's alignments in display order. A branch-site test has one: the codon
    alignment it was given (its input), no amino-acid file of its own."""
    if job.params.get("branch_site"):
        return [Source("tested", "Codon alignment tested", "", "inputs/alignment1.fasta")]
    out: list[Source] = []
    for entry in summary.get("alignments", []):
        if entry["path"] in keys:
            codon = f"{entry['key']}/alignment.cds.fasta"
            out.append(Source(entry["key"], entry.get("label", entry["key"]), entry["path"],
                              codon if codon in keys else ""))  # fmt: skip
    seen = {s.key for s in out}
    for combo in summary.get("combinations", []):  # jobs from before <aligner>/alignment.fasta
        path = f"{combo['key']}/alignment.fasta"
        if combo["aligner"] not in seen and path in keys:
            seen.add(combo["aligner"])
            out.append(Source(combo["aligner"], combo.get("aligner_label", combo["aligner"]), path))
    return out


def read(job: Job, path: str) -> str:
    if path.startswith("inputs/"):
        return storage.read_bytes(f"{job.storage_prefix}{path}").decode()
    return joined.cached(job).read(path)  # an inherited one: the parent's


def genetic_code(job: Job) -> codons.GeneticCode:
    return codons.CODES.get(int(job.params.get("genetic_code") or 1)) or codons.STANDARD


def model_label(combo: dict[str, Any]) -> str:
    fitted = combo.get("fit", {}).get("model", "")
    if catalog.is_modelfinder(combo["model"]) and fitted:
        return f"ModelFinder → {fitted}"
    return str(combo.get("model_label", combo["model"]))


def ancestor_versions(combo_key: str, keys: set[str]) -> list[str]:
    """The gap versions (and codeml reconstructions) this combination has."""
    files = {v: ancestors.VERSIONS[v][1] for v in ancestors.VERSIONS}
    for version in list(files):
        files[version] = f"{combo_key}/{files[version]}"
    for key in sorted(k for k in keys if k.startswith(f"{combo_key}/codeml/")):
        if key.endswith("/ancestors.fasta"):
            model = key.split("/")[-2]
            files[ancestors.CODEML_PREFIX + model] = key
    return [v for v, path in files.items() if path in keys]


def version_path(combo_key: str, version: str, coding: bool) -> str:
    if coding:
        return f"{combo_key}/{ancestors.coding_file(version)}"
    if model := ancestors.codeml_model(version):
        return f"{combo_key}/codeml/{model}/ancestors.fasta"
    return f"{combo_key}/{ancestors.VERSIONS[version][1]}"


def branch_site_series(
    job: Job, owner: Any, source: Source, width: int, keys: set[str]
) -> list[alignview.Series]:
    """P(positive) per column from each branch-site test on this alignment: the job's
    own result for a test, else the user's tests of this job and aligner."""
    if job.params.get("branch_site"):
        tests = [job] if "branch_site/sites.tsv" in keys else []
    else:
        tests = list(
            Job.objects.for_user(owner)
            .filter(params__branch_site_of__job=str(job.id), status=Job.Status.SUCCEEDED,
                    deleted_at__isnull=True)
            .order_by("created_at")
        )  # fmt: skip
        tests = [t for t in tests if t.params["branch_site_of"].get("aligner") == source.key]
    out = []
    for test in tests:
        key = f"{test.storage_prefix}results/branch_site/sites.tsv"
        if not manifest.has(test, "branch_site/sites.tsv"):
            continue
        of = test.params.get("branch_site_of") or {}
        label = f"{of.get('node', '?')} ({of.get('scope', 'branch')})"
        values = alignview.tsv_column(storage.read_bytes(key).decode(), "p_positive", width)
        out.append(alignview.Series(label, values))
    return out


CACHE_SECONDS = 3600


def composition(
    job: Job, source: Source, rows: list[str], reference: str
) -> dict[str, alignview.Track]:
    """The alignment-only tracks, cached per job and alignment (a finished job's results
    never change; 1.4 s at 500 x 5000). Only the hydropathy profile depends on the
    reference row, so it is recomputed each time (cheap)."""
    key = f"alignview:v1:{job.id}:{source.key}"
    tracks: dict[str, alignview.Track] | None = cache.get(key)
    if tracks is None:
        tracks = alignview.composition_tracks(rows)
        cache.set(key, tracks, CACHE_SECONDS)
    profile = tracks["hydropathy_profile"]
    profile.series = [alignview.Series(profile.label, alignview.profile(reference or rows[0]))]
    return tracks


def legend(scheme: str) -> list[tuple[str, str]]:
    """(label, CSS colour) pairs; alignment.js uses the same colours."""
    return LEGENDS.get(scheme, [])


LEGENDS = {
    "clustalx": [("Hydrophobic", "#80a0f0"), ("Positive", "#f01505"), ("Negative", "#c048c0"),
                 ("Polar", "#15c015"), ("Cysteine", "#f08080"), ("Glycine", "#f09048"),
                 ("Proline", "#c0c000"), ("Aromatic", "#15a4a4"),
                 ("Below its column's threshold", "transparent")],
    "chemistry": [("Acidic D E", "#eaa57b"), ("Basic K R H", "#7fb1d4"),
                  ("Aliphatic A V L I M", "#f5ee9a"), ("Aromatic F W Y", "#e3b3cc"),
                  ("Polar S T N Q", "#7fcbb4"), ("Special G P C", "#b0d9f2")],
    "differences": [("Same as the reference (.)", "transparent"),
                    ("Different (coloured by chemistry)", "#7fb1d4"),
                    ("Gap where the reference has a residue, or the reverse", "#c8ccd4")],
    "hydropathy": [("Hydrophilic (−4.5)", "#6f9fcf"), ("0", "#ffffff"),
                   ("Hydrophobic (+4.5)", "#ec9e4a")],
    "identity": [("Matches the consensus, > 80%", "#6464ff"), ("> 60%", "#9999ff"),
                 ("> 40%", "#ccccff")],
    "structure": [("Helix former", "#eaa57b"), ("Strand former", "#7fb1d4"),
                  ("Turn former", "#7fcbb4"), ("Stronger shade: higher propensity", "transparent")],
    "synonymous": [("Same codon as the reference (.)", "transparent"),
                   ("Synonymous (same amino acid)", "#a8d8f2"),
                   ("Non-synonymous", "#f0c56b"), ("Gap difference", "#c8ccd4")],
}  # fmt: skip


def track_list(
    query: Any, available: list[str], remembered: list[str] | None = None
) -> list[str]:
    """The tracks shown (?track=..., in order), after the one-off edits: add_track=<key>,
    remove_track=<i>, track_up=<i>. Unknown and unavailable keys are dropped. A URL that
    names none starts from `remembered` (what the reader last chose for this alignment;
    an empty list is "no tracks"), else from the defaults."""
    if "track" in query or "tracks_set" in query:
        chosen = [k for k in query.getlist("track") if k in available]
    elif remembered is not None:
        chosen = [k for k in remembered if k in available]
    else:
        chosen = [k for k in alignview.DEFAULT_TRACKS if k in available]
    chosen = list(dict.fromkeys(chosen))
    add = query.get("add_track", "")
    if add in available and add not in chosen:
        chosen.append(add)
    remove = query.get("remove_track", "")
    if remove.isdigit() and int(remove) < len(chosen):
        del chosen[int(remove)]
    up = query.get("track_up", "")
    if up.isdigit() and 0 < int(up) < len(chosen):
        i = int(up)
        chosen[i - 1], chosen[i] = chosen[i], chosen[i - 1]
    return chosen[:MAX_TRACKS]


ONE_OFF = ("add_track", "remove_track", "track_up")


def remembered_tracks(user: Any, job: Job, aligner: str) -> list[str] | None:
    """The tracks this user last chose for this alignment (None: never chose)."""
    if not getattr(user, "is_authenticated", False):
        return None
    row = AlignmentTracks.objects.filter(owner=user, job=job, aligner=aligner).first()
    if row is None or not isinstance(row.tracks, list):
        return None
    return [k for k in row.tracks if isinstance(k, str)]


def remember_tracks(user: Any, job: Job, aligner: str, tracks: list[str]) -> None:
    if getattr(user, "is_authenticated", False):
        AlignmentTracks.objects.update_or_create(
            owner=user, job=job, aligner=aligner, defaults={"tracks": tracks}
        )


def alignment_viewer(request: HttpRequest, job_id: str) -> HttpResponse:
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    context: dict[str, Any] = {"job": job, "ready": False}
    keys: set[str] = set()
    summary: dict[str, Any] = {}
    if job.status == Job.Status.SUCCEEDED:
        keys, summary = job_keys(job)
    available = sources(job, keys, summary) if summary else []
    if available:
        context |= viewer_context(request, job, keys, summary, available)
        export = request.GET.get("export", "")
        if export in EXPORTS:
            return download(context, export)
        measured = indel_measurement(job, context["source"], keys, summary)
        if export == "indels" and measured:
            body = indel_views.tsv(measured)
            response = HttpResponse(body, content_type="text/plain; charset=utf-8")
            filename = f"{context['file_stem']}.indel-lengths.tsv"
            response["Content-Disposition"] = f'attachment; filename="{filename}"'
            return response
        if measured:
            context["indel_chart"] = indel_views.large(measured)
    htmx = bool(request.headers.get("HX-Request"))
    if context["ready"] and any(k in request.GET for k in ONE_OFF) and not htmx:
        return redirect(context["canonical_url"])
    context["htmx"] = htmx
    template = "_alignment_viewer.html" if htmx else "alignment.html"
    response = render(request, f"phylobot/{template}", context)
    if context["ready"] and htmx:
        response["HX-Push-Url"] = context["canonical_url"]
    return response


def indel_measurement(
    job: Job, source: Source, keys: set[str], summary: dict[str, Any]
) -> indel_views.Measured | None:
    if not source.amino:  # a branch-site test's codon alignment
        return None
    return indel_views.measure(
        job, source.key, source.label, source.amino, summary, keys, lambda p: read(job, p)
    )


def viewer_context(
    request: HttpRequest, job: Job, keys: set[str], summary: dict[str, Any], available: list[Source]
) -> dict[str, Any]:
    query = request.GET
    by_key = {s.key: s for s in available}
    source = by_key.get(query.get("aligner", ""), available[0])
    states = (
        "codons"
        if source.codon and (query.get("states") == "codons" or not source.amino)
        else "amino"
    )
    code = genetic_code(job)

    extant_amino: dict[str, str]
    extant_codon: dict[str, str] = {}
    if source.codon and states == "codons" or not source.amino:
        extant_codon = alignview.parse_fasta(read(job, source.codon))
    if source.amino:
        extant_amino = alignview.parse_fasta(read(job, source.amino))
    else:
        extant_amino = {n: codons.translate_row(r, code) for n, r in extant_codon.items()}
    names = list(extant_amino)
    width = max((len(r) for r in extant_amino.values()), default=0)

    # The combination (tree source + model) for row order, ancestors and tree tracks.
    primary = summary.get("primary", "")
    combos = sorted(
        (c for c in summary.get("combinations", [])
         if c["aligner"] == source.key and f"{c['key']}/tree.nwk" in keys),
        key=results.combination_order(primary),
    )  # fmt: skip
    tree_sources = list(dict.fromkeys(c["tree"] for c in combos))
    tree_labels = {c["tree"]: c.get("tree_label", c["tree"]) for c in combos}
    combo: dict[str, Any] | None = None
    if combos:
        tree_source = query.get("tree_source", "")
        requested = next((c for c in combos if c["key"] == query.get("combination")), None)
        if tree_source not in tree_sources:
            tree_source = requested["tree"] if requested else combos[0]["tree"]
        mine = [c for c in combos if c["tree"] == tree_source]
        best = min(
            mine, key=lambda c: (c.get("rank_bic", 0), not catalog.is_modelfinder(c["model"]))
        )
        combo = requested if requested in mine else best
    combo_key = combo["key"] if combo else ""
    rooted = combo is not None and f"{combo_key}/tree.rooted.nwk" in keys
    newick = (
        read(job, f"{combo_key}/{'tree.rooted.nwk' if rooted else 'tree.nwk'}") if combo else ""
    )

    # Ancestors as rows (owner decision: a toggle, off by default).
    versions = ancestor_versions(combo_key, keys) if combo else []
    show_ancestors = bool(versions) and query.get("ancestors") == "on"
    version = query.get("gaps", "")
    if version not in versions:
        chosen = [m for m in job.params.get("gap_methods") or [] if m in versions]
        version = chosen[0] if chosen else ancestors.ORIGINAL
    ancestor_amino: dict[str, str] = {}
    ancestor_codon: dict[str, str] = {}
    codon_ancestors_missing = False
    if show_ancestors:
        ancestor_amino = alignview.parse_fasta(read(job, version_path(combo_key, version, False)))
        if states == "codons":
            path = version_path(combo_key, version, True)
            if path in keys:
                ancestor_codon = alignview.parse_fasta(read(job, path))
            else:
                codon_ancestors_missing = True
                ancestor_amino = {}

    order = "tree" if combo and query.get("order") == "tree" else "alignment"
    rows = list(names)
    if order == "tree":
        listed = [n for n in alignview.tree_order(newick, job.params.get("outgroup") or [],
                                                  internal=bool(ancestor_amino))
                  if n in extant_amino or n in ancestor_amino]  # fmt: skip
        rows = listed + [n for n in names if n not in listed]
        rows += [n for n in ancestor_amino if n not in rows]
    elif ancestor_amino:
        preorder = alignview.tree_order(newick, internal=True)
        rows += [n for n in preorder if n in ancestor_amino]
        rows += [n for n in ancestor_amino if n not in rows]

    scheme_menu = SCHEMES | (CODON_SCHEMES if states == "codons" else {})
    scheme = query.get("scheme", "") if query.get("scheme", "") in scheme_menu else DEFAULT_SCHEME
    reference = query.get("reference", "")
    if reference not in rows:
        reference = ""  # the consensus (owner decision)
    zoom = parse_zoom(query.get("zoom", ""))
    height = query.get("height", "") if query.get("height", "") in HEIGHTS else DEFAULT_HEIGHT

    # Tracks offered: composition always; tree tracks where this combination has them.
    offered = list(alignview.COMPOSITION)
    if combo:
        if f"{combo_key}/site_rates.tsv" in keys:
            offered.append("rate")
        if f"{combo_key}/site_lh.tsv" in keys:
            offered.append("site_lh")
        if f"{combo_key}/ancestors.fasta" in keys:
            offered.append("changes")
        if f"{combo_key}/ancestors_pp.tsv" in keys:
            offered.append("ancestor_pp")
        if any(v != ancestors.ORIGINAL and not ancestors.codeml_model(v) for v in versions):
            offered.append("gap_calls")
    bs_series = branch_site_series(job, request.user, source, width, keys) if source.codon else []
    if bs_series:
        offered.append("branch_site")
    exporting = query.get("export") == "tracks"
    shown = track_list(query, offered, remembered_tracks(request.user, job, source.key))
    if any(name in query for name in ONE_OFF):
        # The reader changed the tracks: this alignment opens with them from now on.
        remember_tracks(request.user, job, source.key, shown)
    extant_rows = [extant_amino[n] for n in names]
    ref_amino = extant_amino.get(reference, "") or ancestor_amino.get(reference, "")
    computed = composition(job, source, extant_rows, ref_amino)
    wanted = offered if exporting else shown

    def tree_track(key: str) -> alignview.Track | None:
        kind = alignview.TRACKS[key]
        series: list[alignview.Series] = []
        low: float | None = None
        high: float | None = None
        if key == "rate":
            series = [alignview.Series(kind.label, alignview.tsv_column(
                read(job, f"{combo_key}/site_rates.tsv"), "rate", width))]  # fmt: skip
            low = 0.0
        elif key == "site_lh":
            series = [alignview.Series(kind.label, alignview.tsv_column(
                read(job, f"{combo_key}/site_lh.tsv"), "log_likelihood", width))]  # fmt: skip
        elif key == "changes":
            sequences = extant_amino | alignview.parse_fasta(
                read(job, f"{combo_key}/ancestors.fasta")
            )
            values = alignview.substitutions(read(job, f"{combo_key}/tree.nwk"), sequences, width)
            series, low = [alignview.Series(kind.label, values)], 0.0
        elif key == "ancestor_pp":
            values = alignview.mean_pp(read(job, f"{combo_key}/ancestors_pp.tsv"), width)
            series, low, high = [alignview.Series(kind.label, values)], 0.0, 1.0
        elif key == "gap_calls":
            masked = {
                dict(results.GAP_VERSIONS).get(v, v): alignview.parse_fasta(
                    read(job, version_path(combo_key, v, False)))
                for v in versions if v != ancestors.ORIGINAL and not ancestors.codeml_model(v)
            }  # fmt: skip
            series, low, high = alignview.gap_shares(masked, width), 0.0, 1.0
        elif key == "branch_site":
            series, low, high = bs_series, 0.0, 1.0
        if not series:
            return None
        return alignview.Track(key, kind.label, kind.description, series, low, high)

    tracks = []
    for key in wanted:
        track = computed.get(key) or tree_track(key)
        if track is not None:
            tracks.append(track)

    # The view's canonical query: only resolved values.
    view: dict[str, Any] = {"aligner": source.key}
    if source.codon and source.amino:
        view["states"] = states
    if combo:
        view |= {"tree_source": combo["tree"], "combination": combo_key, "order": order}
        if show_ancestors:
            view |= {"ancestors": "on", "gaps": version}
    view |= {"scheme": scheme, "zoom": zoom, "height": height, "tracks_set": "1", "track": shown}
    if reference:
        view["reference"] = reference
    here = reverse("alignment_viewer", args=[job.id])

    all_amino = extant_amino | ancestor_amino
    all_codon = extant_codon | ancestor_codon
    grid_rows = [all_codon[n] if states == "codons" else all_amino[n] for n in rows]
    ancestor_viewer = ""
    if combo and f"{combo_key}/ancestors.fasta" in keys:
        query_part = {"aligner": source.key, "combination": combo_key}
        if show_ancestors:
            query_part["gaps"] = version
        ancestor_viewer = f"{reverse('ancestor_viewer', args=[job.id])}?{urlencode(query_part)}"
    data = {
        "names": rows,
        "ancestor": [n in ancestor_amino and n not in extant_amino for n in rows],
        "rows": grid_rows,
        "amino": [all_amino[n] for n in rows] if states == "codons" else [],
        "codon": states == "codons",
        "extant": len(names),
        "reference": rows.index(reference) if reference else -1,
        "scheme": scheme,
        "cell": zoom,
        "rowFromCell": True,
        "height": height,
        "tracks": [t.json() for t in tracks if t.key in shown],
        "ancestorViewer": ancestor_viewer,
    }
    stem = re.sub(r"[^A-Za-z0-9._+-]", "_", f"alignment.{source.key}")
    return {
        "ready": True,
        "editor_url": editor_url(
            job,
            source,
            {
                "scheme": scheme,
                "zoom": zoom,
                "height": height,
                "tracks_set": "1",
                "track": [k for k in shown if k in alignview.COMPOSITION],
            },
        ),  # fmt: skip
        "editor_label": source.label,
        "canonical_url": f"{here}?{urlencode(view, doseq=True)}",
        "export_url": f"{here}?{urlencode(view, doseq=True)}",
        "sources": [(s.key, s.label) for s in available],
        "source": source,
        "states": states,
        "tree_sources": [(t, tree_labels[t]) for t in tree_sources],
        "tree_source": combo["tree"] if combo else "",
        "models": [
            (c["key"], model_label(c)) for c in combos if combo and c["tree"] == combo["tree"]
        ],
        "combination": combo_key,
        "order": order,
        "has_ancestors": bool(versions),
        "show_ancestors": show_ancestors,
        "codon_ancestors_missing": codon_ancestors_missing,
        "versions": [(v, ancestors.version_label(v)) for v in versions],
        "version": version,
        "schemes": list(scheme_menu.items()),
        "scheme": scheme,
        "legend": legend(scheme),
        "references": rows,
        "reference": reference,
        "zoom": zoom,
        "heights": list(HEIGHTS.items()),
        "height": height,
        "shown": [(k, alignview.TRACKS[k]) for k in shown],
        "addable": [
            (
                group,
                [
                    (k, alignview.TRACKS[k])
                    for k in offered
                    if k not in shown and alignview.TRACKS[k].group == group
                ],
            )
            for group in dict.fromkeys(alignview.TRACKS[k].group for k in offered if k not in shown)
        ],  # fmt: skip
        "can_add": len(shown) < MAX_TRACKS,
        "data": data,
        "sequences": len(names),
        "columns": width,
        "tracks": tracks,
        "download_rows": dict(zip(rows, grid_rows, strict=True)),
        "file_stem": stem + (".cds" if states == "codons" else ""),
        "tree_url": (
            f"{reverse('tree_viewer', args=[job.id])}?{urlencode({'combination': combo_key})}"
            if combo
            else ""
        ),
    }


def editor_url(job: Job, source: Source, view: dict[str, Any]) -> str:
    """This alignment in the alignment editor (as registered in My alignments), shown as
    here: colours, zoom, height and the tracks the editor has."""
    if job.params.get("branch_site"):
        return ""
    record = Alignment.objects.filter(job=job, job_key=source.key).first()
    if record is None:
        try:
            alignments.register_job(job)
        except Exception:
            return ""
        record = Alignment.objects.filter(job=job, job_key=source.key).first()
    if record is None:
        return ""
    return f"{reverse('alignment_editor', args=[record.id])}?{urlencode(view, doseq=True)}"


def download(context: dict[str, Any], what: str) -> HttpResponse:
    """The alignment as shown (rows, order, ancestors if shown) as FASTA or Clustal,
    or every track as TSV."""
    if not context.get("ready"):
        raise Http404
    stem = context["file_stem"]
    rows = context["download_rows"]
    if what == "fasta":
        body = "".join(f">{name}\n{row}\n" for name, row in rows.items())
        filename = f"{stem}.fasta"
    elif what == "clustal":
        body = alignview.clustal(rows, protein=context["states"] != "codons")
        filename = f"{stem}.aln"
    else:
        body = alignview.tracks_tsv(context["tracks"], context["columns"])
        filename = f"{stem}.tracks.tsv"
    response = HttpResponse(body, content_type="text/plain; charset=utf-8")
    response["Content-Disposition"] = f'attachment; filename="{filename}"'
    return response
