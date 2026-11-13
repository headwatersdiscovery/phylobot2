"""My usage: the user's CPU time and storage, over a chosen window (phylobot.usage)."""

import csv
import io

from django.http import Http404, HttpRequest, HttpResponse
from django.shortcuts import render
from django.utils import timezone
from django.views.decorators.http import require_GET

from phylobot import usage
from phylobot.models import Job

# What the page reads: not the inputs or logs, only the stored figures.
FIELDS = ("id", "name", "status", "params", "progress", "created_at", "finished_at",
          "deleted_at", "purged_at", "cpu_seconds", "storage_bytes", "imported_at")  # fmt: skip


@require_GET
def my_usage(request: HttpRequest) -> HttpResponse:
    now = timezone.now()
    today = usage.local_day(now)
    live = Job.objects.for_user(request.user).only(*FIELDS)
    deleted = Job.objects.deleted_for(request.user).only(*FIELDS)
    jobs = [*live, *deleted]
    first = min((usage.local_day(j.created_at) for j in jobs), default=today)
    filters = usage.Filters.parse(request.GET, today, first)
    data = usage.compute(jobs, filters, now)
    if what := request.GET.get("export"):
        return export(request, data, what, today)
    downloads: dict[str, dict[str, str]] = {}
    for key in usage.EXPORTS:
        for fmt in ("csv", "xlsx"):
            query = request.GET.copy()
            query["export"], query["format"] = key, fmt
            downloads.setdefault(key, {})[fmt] = f"{request.path}?{query.urlencode()}"
    context = {
        "downloads": downloads,
        "u": data,
        "presets": usage.PRESETS,
        "kinds": usage.KINDS,
        "kind_help": usage.KIND_HELP,
        "today": today,
        "first_day": first,
        "has_jobs": bool(jobs),
    }
    return render(request, "phylobot/my_usage.html", context)


def export(request: HttpRequest, data: usage.Usage, what: str, today: object) -> HttpResponse:
    """One panel's figures, for the window and filters on the page, as CSV or Excel."""
    fmt = request.GET.get("format", "csv")
    if what not in usage.EXPORTS or fmt not in ("csv", "xlsx"):
        raise Http404
    columns, rows = usage.export_table(data, what)
    window = data.filters.window
    stem = f"phylobot_usage_{what}_{window.start:%Y%m%d}-{window.end:%Y%m%d}"
    if fmt == "csv":
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
        response = HttpResponse(out.getvalue(), content_type="text/csv; charset=utf-8")
    else:
        from phylobot import sites_excel  # openpyxl: only needed here

        about = [("Usage", usage.EXPORTS[what]), ("Window", data.filters.window.label),
                 ("From", f"{data.filters.window.start}"), ("To", f"{data.filters.window.end}"),
                 ("Downloaded", f"{today}")]  # fmt: skip
        response = HttpResponse(
            sites_excel.results_workbook(columns, rows, about, name=usage.EXPORTS[what]),
            content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )
    response["Content-Disposition"] = f'attachment; filename="{stem}.{fmt}"'
    return response
