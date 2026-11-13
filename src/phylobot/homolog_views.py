"""Find homologs: the search form, each search's page, and the review of its hits.

`/homologs/` composes a search and lists the user's searches and the saved examples.
Every search, real (`/homologs/<id>/`) or saved example (`/homologs/examples/<key>/`),
has one landing page (homolog_search.html): its status while it runs, the error if it
failed, the review when it finished.

The review works the same for a real search and a saved example (a `Review`): filters
(with sliders), hand-picked includes and excludes, histograms, the selected and
not-selected hits, FASTA downloads, and "use these sequences in a new job", which opens
the New Job page with step 2 filled in. All review state lives in the query string, so
every view (and download) has a shareable URL and nothing is stored until a job is.
"""

import logging
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any
from urllib.parse import quote, urlencode

from django.conf import settings
from django.contrib.auth.models import User
from django.http import Http404, HttpRequest, HttpResponse, JsonResponse, QueryDict
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.views.decorators.http import require_GET, require_http_methods, require_POST

from phylobot import searches
from phylobot.forms import HomologFilterForm, HomologSearchForm, SequenceLookupForm
from phylobot.homologs import examples as homolog_examples
from phylobot.homologs import export
from phylobot.homologs import filters as homolog_filters
from phylobot.homologs import lookup as homolog_lookup
from phylobot.homologs import plots as homolog_plots
from phylobot.homologs import remote as homolog_remote
from phylobot.homologs import sources as homolog_sources
from phylobot.homologs.hits import Hit, SearchResult
from phylobot.models import HomologSearch

logger = logging.getLogger(__name__)

# Where each service's hit IDs lead, for the hit tables.
HIT_LINKS = {
    "ebi": "https://www.uniprot.org/uniref/{id}",
    "ncbi": "https://www.ncbi.nlm.nih.gov/protein/{id}",
}
TABLES = ("hits",)  # the hit table, whose open/closed state the page keeps
# The hit table's sortable columns (`?sort=<key>`, or `-<key>` for descending).
SORTS = ("use", "rank", "protein", "organism", "coverage", "identity", "evalue", "length",
         "id", "note")  # fmt: skip


def current_user(request: HttpRequest) -> User:
    """The signed-in user (local mode: the local user); middleware guarantees one."""
    assert isinstance(request.user, User)
    return request.user


@dataclass
class Review:
    """A finished search to review: a real search or a saved example."""

    result: SearchResult
    load_sequences: Callable[[], dict[str, str]]
    label: str
    service: str
    panel_url: str  # the review page (the filters re-request it)
    download_url: Callable[[str], str]  # "selected", "all" or "raw" -> its URL
    raw: Callable[[], bytes]  # the service's raw result
    job_param: tuple[str, str]  # what the New Job page is told to load: (name, value)


def example_review(key: str) -> Review:
    if key not in homolog_examples.EXAMPLES:
        raise Http404
    example = homolog_examples.EXAMPLES[key]
    return Review(
        result=homolog_examples.load(key),
        load_sequences=lambda: homolog_examples.sequences(key),
        label=example.label,
        service=example.source,
        panel_url=reverse("homolog_example", args=[key]),
        download_url=lambda which: reverse("homolog_example_download", args=[key, which]),
        raw=lambda: homolog_examples.raw(key),
        job_param=("homologs_example", key),
    )


def search_review(search: HomologSearch) -> Review:
    from phylobot import storage

    return Review(
        result=searches.result(search),
        load_sequences=lambda: searches.full_sequences(search),
        label=search.name or "Your search",
        service=search.service,
        panel_url=reverse("homolog_search", args=[search.id]),
        download_url=lambda which: reverse("homolog_search_download", args=[search.id, which]),
        raw=lambda: storage.read_bytes(search.result_key),
        job_param=("homologs", str(search.id)),
    )


@dataclass
class State:
    """The review's settings as the query string gives them."""

    form: HomologFilterForm
    settings: homolog_filters.Settings
    include: list[str]
    exclude: list[str]
    show: set[str]
    base: str = "filters"
    sort: str = ""  # "" = search order


def read_state(query: QueryDict) -> State:
    """Filters missing from the query string take their defaults, so a link that gives
    only some of them (e.g. ?min_coverage=55) still applies those."""
    bound = any(name in query for name in HomologFilterForm.base_fields)
    data = None
    if bound:
        defaults = HomologFilterForm()
        data = {
            name: query.get(name, defaults[name].initial) for name in HomologFilterForm.base_fields
        }
    form = HomologFilterForm(data)
    chosen = homolog_filters.Settings()
    if bound and form.is_valid():
        chosen = form.cleaned_data["settings"]
    return State(
        form=form,
        settings=chosen,
        include=query.getlist("include"),
        exclude=query.getlist("exclude"),
        show={name for name in query.getlist("show") if name in TABLES},
        base=query.get("base", "filters"),  # select() falls back to "filters"
        sort=sort if (sort := query.get("sort", "")).lstrip("-") in SORTS else "",
    )


def state_query(
    review: Review, state: State, selection: homolog_filters.Selection, choices: bool = True
) -> str:
    """The review's state as a query string (only known hit IDs), for links; without the
    hand-picked choices if `choices` is False."""
    s = state.settings
    params: list[tuple[str, str]] = [
        ("max_evalue", f"{s.max_evalue:g}"),
        ("min_coverage", f"{s.min_coverage:g}"),
        ("min_identity", f"{s.min_identity:g}"),
        ("max_identity", f"{s.max_identity:g}"),
    ]
    if choices:
        if selection.base != "filters":
            params.append(("base", selection.base))
        params += [("include", i) for i in selection.include_ids]
        params += [("exclude", i) for i in selection.exclude_ids]
    params += [("show", table) for table in sorted(state.show)]
    if state.sort:
        params.append(("sort", state.sort))
    return urlencode(params)


def select(review: Review, state: State) -> homolog_filters.Selection:
    return homolog_filters.select(
        review.result.hits, state.settings, state.include, state.exclude, state.base,
        always=homolog_filters.query_ids(review.result.hits),
    )  # fmt: skip


def hit_link(service: str, hit_id: str) -> str:
    """The hit's page at its database: EBI's Swiss-Prot hits are UniProtKB entry names,
    its other databases' UniRef clusters."""
    if service == "ebi" and not hit_id.startswith("UniRef"):
        return f"https://www.uniprot.org/uniprotkb/{quote(hit_id, safe='')}"
    return HIT_LINKS[service].format(id=quote(hit_id, safe=""))


def hit_rows(review: Review, state: State, selection: homolog_filters.Selection
             ) -> list[dict[str, Any]]:  # fmt: skip
    """One row per hit, sorted as `state.sort` says (search order by default), with the
    note: the user's own choice, why the filters remove it, a missing sequence."""
    selected_ids = {h.id for h in selection.selected}
    kept_ids = {h.id for h in selection.by_filters}
    included = {h.id for h in selection.included}
    excluded = {h.id for h in selection.excluded}
    have = set(review.load_sequences())

    def row(rank: int, hit: Hit) -> dict[str, Any]:
        # The query is kept by default (homolog_filters.is_query_hit); it shows as the query,
        # not as a hand-pick, unless the user has unticked it (then it is an exclude).
        is_query = hit.id in selection.always and hit.id not in excluded
        picked = "included" if hit.id in included else "excluded" if hit.id in excluded else ""
        choice = "" if is_query else picked
        if is_query:
            reason, notes = "", ["the query sequence (kept despite the identity filter)"]
        else:
            reason = "" if hit.id in kept_ids else homolog_filters.reason(hit, state.settings)
            notes = [f"{choice} by you"] if choice else []
            notes += [reason] if reason else []
        notes += [] if hit.id in have else ["no full-length sequence"]
        return {
            "rank": rank,
            "hit": hit,
            "link": hit_link(review.service, hit.id),
            "selected": hit.id in selected_ids,
            "default": hit.id in selection.default,  # what a tick is measured against
            "choice": choice,
            "is_query": is_query,
            "reason": reason,
            "has_sequence": hit.id in have,
            "note": "; ".join(notes),
        }

    rows = [row(rank, hit) for rank, hit in enumerate(review.result.hits, 1)]
    key = state.sort.lstrip("-")
    if key and key != "rank":
        by: dict[str, Callable[[dict[str, Any]], Any]] = {
            "use": lambda r: not r["selected"],
            "protein": lambda r: r["hit"].name.casefold(),
            "organism": lambda r: r["hit"].taxon.casefold(),
            "coverage": lambda r: r["hit"].coverage,
            "identity": lambda r: r["hit"].identity,
            "evalue": lambda r: r["hit"].evalue,
            "length": lambda r: r["hit"].length,
            "id": lambda r: r["hit"].id.casefold(),
            "note": lambda r: (r["note"] == "", r["note"].casefold()),  # blank notes last
        }
        rows.sort(key=by[key])  # stable: ties stay in search order
    if state.sort.startswith("-"):
        rows.reverse()
    return rows


def review_context(review: Review, state: State) -> dict[str, Any]:
    """The review panel: filters, counts, histograms, the hit tables and the actions."""
    selection = select(review, state)
    have = set(review.load_sequences())
    rows = hit_rows(review, state, selection)
    query = state_query(review, state, selection)
    usable = sum(1 for h in selection.selected if h.id in have)
    job_name, job_value = review.job_param
    return {
        "review": review,
        "service": homolog_sources.SOURCES[review.service].label,
        "search": review.result,
        "panel_url": review.panel_url,
        "form": state.form,
        "selection": selection,
        "stages": selection.stages,
        "histograms": homolog_plots.filter_histograms(review.result.hits, state.settings),
        "rows": rows,
        "not_selected": len(rows) - len(selection.selected),
        "show": state.show,
        "sort": state.sort,
        "sort_columns": SORT_COLUMNS,
        "usable": usable,
        "no_sequence": len(selection.selected) - usable,
        "downloads": {
            which: f"{review.download_url(which)}?{query}"
            for which in ("selected", "all", "raw", "csv")
        },
        "use_url": f"{reverse('submit')}?{job_name}={quote(job_value)}&{query}",
        "reset_url": f"{review.panel_url}?{state_query(review, state, selection, choices=False)}",
    }


# The hit table's columns: (sort key, heading, numeric).
SORT_COLUMNS = [
    ("use", "Use", False), ("rank", "#", True), ("protein", "Protein", False),
    ("organism", "Organism", False), ("coverage", "Coverage", True),
    ("identity", "Identity", True), ("evalue", "E-value", True), ("length", "Length", True),
    ("id", "ID", False), ("note", "Note", False),
]  # fmt: skip


def render_review(
    request: HttpRequest, review: Review, context: dict[str, Any], template: str
) -> HttpResponse:
    context |= review_context(review, read_state(request.GET))
    if request.headers.get("HX-Request"):
        return render(request, "phylobot/_homolog_filters.html", context)
    return render(request, template, context)


# --- pages -------------------------------------------------------------------------


@require_http_methods(["GET", "POST"])
def homologs(request: HttpRequest) -> HttpResponse:
    """The search form, the user's searches, and saved example searches (their review
    opens here). A submitted search is queued for run_searches and followed on its page."""
    user = current_user(request)
    mode = settings.PHYLOBOT_MODE
    form = HomologSearchForm(
        request.POST if request.method == "POST" else None,
        mode=mode,
        typical=searches.typical_seconds(),
    )
    if request.method == "POST" and form.is_valid():
        service, database = form.cleaned_data["database"].split(":", 1)
        try:
            searches.check_can_search(user)
            search = searches.create(
                user,
                homolog_remote.Query(form.cleaned_data["query"]),
                service,
                database,
                name=form.cleaned_data["name"],
            )
        except (searches.LimitError, homolog_remote.SearchError) as exc:
            form.add_error(None, str(exc))
        else:
            return redirect("homolog_search", search_id=search.id)
    key = request.GET.get("example", "")
    if key:  # links from before examples had their own page
        if key not in homolog_examples.EXAMPLES:
            raise Http404
        rest = request.GET.copy()
        del rest["example"]
        url = reverse("homolog_example", args=[key])
        return redirect(f"{url}?{rest.urlencode()}" if rest else url)
    context: dict[str, Any] = {
        "search_form": form,
        "lookup_form": SequenceLookupForm(),
        "model_species": homolog_lookup.MODEL_SPECIES,
        "examples": example_rows(),
        "local_mode": mode == "local",
        "searches": HomologSearch.objects.for_user(user)[:50],
    }
    return render(request, "phylobot/homologs.html", context)


@require_GET
def find_sequence(request: HttpRequest) -> HttpResponse:
    """Upstream of homolog search: look up a protein's amino-acid sequence by name and
    species at UniProt (phylobot.homologs.lookup), so a user who has only a name can find a
    sequence to search with. The chosen hit fills the search box below (homologs.js).

    Only the name and species leave the machine; they are cleaned before the query is built.
    HTMX requests get just the results fragment; without JavaScript the whole page renders
    with the results in place."""
    user = current_user(request)
    form = SequenceLookupForm(request.GET or None)
    results: list[homolog_lookup.Found] | None = None
    error = ""
    if request.GET and form.is_valid():
        try:
            results = homolog_lookup.search(
                form.cleaned_data["protein"],
                form.cleaned_data["species"],
                form.cleaned_data["reviewed_only"],
            )
        except homolog_remote.SearchError as exc:
            error = str(exc)
    elif request.GET:
        # Bound but invalid: only this fragment is swapped in, so surface the reason here
        # (the first field or consent error) rather than silently showing nothing.
        error = str(
            next((e for errors in form.errors.values() for e in errors), "Please fill in the form.")
        )
    lookup_context: dict[str, Any] = {
        "lookup_form": form,
        "lookup_results": results,
        "lookup_error": error,
        "lookup_done": results is not None or bool(error),
        "model_species": homolog_lookup.MODEL_SPECIES,
    }
    if request.headers.get("HX-Request"):
        return render(request, "phylobot/_sequence_lookup_results.html", lookup_context)
    # No JavaScript: render the whole page with the results in the lookup section.
    mode = settings.PHYLOBOT_MODE
    context: dict[str, Any] = {
        "search_form": HomologSearchForm(mode=mode, typical=searches.typical_seconds()),
        "examples": example_rows(),
        "local_mode": mode == "local",
        "searches": HomologSearch.objects.for_user(user)[:50],
        **lookup_context,
    }
    return render(request, "phylobot/homologs.html", context)


@require_GET
def suggest(request: HttpRequest) -> JsonResponse:
    """Type-ahead for the lookup fields (phylobot.homologs.lookup): `?kind=species|gene&
    q=<prefix>` (gene also takes `&species=`). Returns {"items": [{value, label}]} from
    UniProt, cached; failures and short prefixes return an empty list, never an error, so
    the field keeps working. Signed-in users only (as the rest of /homologs/)."""
    current_user(request)
    kind = request.GET.get("kind", "")
    prefix = request.GET.get("q", "")
    if kind == "species":
        items = homolog_lookup.suggest_species(prefix)
    elif kind == "gene":
        items = homolog_lookup.suggest_genes(prefix, request.GET.get("species", ""))
    else:
        items = ()
    return JsonResponse({"items": [{"value": s.value, "label": s.label} for s in items]})


def example_rows() -> list[dict[str, Any]]:
    """The featured saved examples as rows like the user's own searches (homologs page)."""
    rows = []
    for key in homolog_examples.FEATURED:
        example = homolog_examples.EXAMPLES[key]
        rows.append({"key": key, "example": example, "result": homolog_examples.load(key)})
    return rows


@require_GET
def homolog_example(request: HttpRequest, key: str) -> HttpResponse:
    """A saved example's landing page: the same page as a finished search of your own."""
    review = example_review(key)
    context: dict[str, Any] = {
        "saved_example": homolog_examples.EXAMPLES[key],
        "saved_example_key": key,
    }
    return render_review(request, review, context, "phylobot/homolog_search.html")


def search_context(search: HomologSearch) -> dict[str, Any]:
    source = homolog_sources.SOURCES.get(search.service)
    database = source.databases.get(search.database, search.database) if source else search.database
    return {
        "homolog_search": search,
        "service_label": source.label if source else search.service,
        "database_label": database,
        "active": search.status in searches.ACTIVE,
    }


def own_search(request: HttpRequest, search_id: str) -> HomologSearch:
    return get_object_or_404(HomologSearch.objects.for_user(request.user), id=search_id)


@require_GET
def homolog_search(request: HttpRequest, search_id: str) -> HttpResponse:
    """One search: its progress while it waits at the service, then the review."""
    search = own_search(request, search_id)
    context = search_context(search)
    if search.status == HomologSearch.Status.SUCCEEDED:
        try:
            review = search_review(search)
        except Exception:
            logger.exception("search %s: stored result unreadable", search.id)
            context["result_error"] = True
        else:
            return render_review(request, review, context, "phylobot/homolog_search.html")
    return render(request, "phylobot/homolog_search.html", context)


@require_GET
def homolog_search_status(request: HttpRequest, search_id: str) -> HttpResponse:
    """The search's status, polled by HTMX while it waits; reloads the page when done."""
    search = own_search(request, search_id)
    response = render(request, "phylobot/_homolog_search_status.html", search_context(search))
    if search.status not in searches.ACTIVE:
        response["HX-Refresh"] = "true"
    return response


def next_page(request: HttpRequest, default: str) -> str:
    """Where to go after a POST: a Find homologs page named by the form, never elsewhere."""
    target = request.POST.get("next", "")
    if target.startswith("/homologs/") and "//" not in target and "\\" not in target:
        return target
    return default


@require_POST
def delete_homolog_searches(request: HttpRequest) -> HttpResponse:
    """Delete the selected searches (one from its page, or many from the list): files go,
    records stay. IDs that are not the user's own live searches are ignored."""
    from django.contrib import messages

    ids = [value for value in request.POST.getlist("search") if is_uuid(value)]
    chosen = list(HomologSearch.objects.for_user(request.user).filter(id__in=ids))
    for search in chosen:
        searches.delete(search)
    if chosen:
        n = len(chosen)
        messages.success(request, f"Deleted {n} search{'es' if n != 1 else ''}.")
    return redirect(next_page(request, reverse("homologs")))


@require_POST
def rename_homolog_search(request: HttpRequest, search_id: str) -> HttpResponse:
    """Rename from the search's page (HTMX: the title part comes back, done editing)."""
    search = own_search(request, search_id)
    searches.rename(search, request.POST.get("name", ""))
    if request.headers.get("HX-Request"):
        return render(request, "phylobot/_homolog_title.html", {"homolog_search": search})
    return redirect("homolog_search", search_id=search.id)


def after_duplicate(request: HttpRequest, copy: HomologSearch) -> HttpResponse:
    from django.contrib import messages

    messages.success(request, f"Created “{copy.name}”, a copy you can rename or delete.")
    return redirect(next_page(request, reverse("homolog_search", args=[copy.id])))


@require_POST
def duplicate_homolog_search(request: HttpRequest, search_id: str) -> HttpResponse:
    """Copy a finished search as a new one of the user's (nothing is sent again)."""
    search = own_search(request, search_id)
    if search.status != HomologSearch.Status.SUCCEEDED:
        raise Http404
    return after_duplicate(request, searches.duplicate(search, current_user(request)))


@require_POST
def duplicate_homolog_example(request: HttpRequest, key: str) -> HttpResponse:
    """Copy a saved example into the user's own searches."""
    if key not in homolog_examples.EXAMPLES:
        raise Http404
    return after_duplicate(request, searches.duplicate_example(key, current_user(request)))


def is_uuid(text: str) -> bool:
    import uuid

    try:
        uuid.UUID(text)
    except ValueError:
        return False
    return True


# --- downloads ---------------------------------------------------------------------


def download(request: HttpRequest, review: Review, which: str) -> HttpResponse:
    """The selected hits or all hits as FASTA (full-length sequences), or the service's
    raw result."""
    stem = "".join(c if c.isalnum() else "-" for c in review.label.lower()).strip("-")[:40]
    if which == "raw":
        response = HttpResponse(review.raw(), content_type="application/xml")
        filename = f"{stem or 'search'}.blast.xml"
    elif which in ("selected", "all"):
        hits = select(review, read_state(request.GET)).selected
        if which == "all":
            hits = review.result.hits
        fasta, _ = export.to_fasta(hits, review.load_sequences())
        response = HttpResponse(fasta, content_type="text/plain; charset=utf-8")
        filename = f"{stem or 'search'}.{which}.fasta"
    elif which == "csv":
        response = HttpResponse(hits_csv(review, read_state(request.GET)),
                                content_type="text/csv; charset=utf-8")  # fmt: skip
        filename = f"{stem or 'search'}.hits.csv"
    else:
        raise Http404
    response["Content-Disposition"] = f'attachment; filename="{filename}"'
    return response


def hits_csv(review: Review, state: State) -> str:
    """The hit table as shown (same order, every hit), for a spreadsheet."""
    import csv
    import io

    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(["selected", "rank", "id", "protein", "organism", "taxid", "coverage_pct",
                     "identity_pct", "evalue", "bits", "length", "note"])  # fmt: skip
    for r in hit_rows(review, state, select(review, state)):
        h = r["hit"]
        writer.writerow([
            "yes" if r["selected"] else "no", r["rank"], h.id, h.name, h.taxon,
            h.taxid if h.taxid is not None else "", f"{h.coverage:.1f}", f"{h.identity:.1f}",
            f"{h.evalue:.3g}", f"{h.bits:g}", h.length, r["note"],
        ])  # fmt: skip
    return out.getvalue()


@require_GET
def homolog_search_download(request: HttpRequest, search_id: str, which: str) -> HttpResponse:
    search = own_search(request, search_id)
    if search.status != HomologSearch.Status.SUCCEEDED:
        raise Http404
    return download(request, search_review(search), which)


@require_GET
def homolog_example_download(request: HttpRequest, key: str, which: str) -> HttpResponse:
    return download(request, example_review(key), which)


# --- into a job --------------------------------------------------------------------


def job_prefill(request: HttpRequest) -> dict[str, Any] | None:
    """For the New Job page: the selected hits as step 2's sequences, when it is opened
    with ?homologs=<search> or ?homologs_example=<key> (plus the review's state). The
    sequences then go through the normal form and its validation; nothing is submitted."""
    if "homologs" in request.GET:
        search = own_search(request, request.GET["homologs"])
        if search.status != HomologSearch.Status.SUCCEEDED:
            raise Http404
        review = search_review(search)
    elif "homologs_example" in request.GET:
        review = example_review(request.GET["homologs_example"])
    else:
        return None
    selection = select(review, read_state(request.GET))
    fasta, missing = export.to_fasta(selection.selected, review.load_sequences())
    return {
        "fasta": fasta,
        # each hit's taxid, for the reconciliation block's species table (species_web)
        "taxids": export.taxids(selection.selected),
        "count": len(selection.selected) - len(missing),
        "missing": len(missing),
        "label": review.label,
        "included": len(selection.included),
        "excluded": len(selection.excluded),
        "back_url": f"{review.panel_url}?{state_query(review, read_state(request.GET), selection)}",
    }
