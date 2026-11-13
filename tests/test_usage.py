"""My usage (phylobot.usage, usage_views): windows, buckets, the figures and the page."""

from datetime import UTC, date, datetime, timedelta
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.http import QueryDict
from django.test import Client

from phylobot import usage
from phylobot.models import Job
from phylobot.modes import local_user

TODAY = date(2026, 10, 16)
NOW = datetime(2026, 10, 16, 18, tzinfo=UTC)


def filters(query: str = "", first: date = date(2026, 9, 1)) -> usage.Filters:
    return usage.Filters.parse(QueryDict(query), TODAY, first)


def test_windows() -> None:
    assert filters().window == usage.Window(date(2026, 7, 19), TODAY, "90d")  # the default
    assert filters("range=7d").window.start == date(2026, 10, 10)
    assert filters("range=all").window == usage.Window(date(2026, 9, 1), TODAY, "all")
    custom = filters("range=custom&from=2026-10-01&to=2026-12-31").window
    assert custom == usage.Window(date(2026, 10, 1), TODAY, "custom")  # never past today
    assert filters("range=custom&from=bad").window.start == date(2026, 9, 1)
    assert filters("range=nonsense").window.preset == "90d"
    assert filters("range=30d").window.previous() == usage.Window(
        date(2026, 8, 18), date(2026, 9, 16), ""
    )
    assert filters("range=all").window.previous() is None
    assert filters("kind=dnds&kind=bogus").kinds == ["dnds"]
    assert filters("deleted=0").include_deleted is False
    assert filters("deleted=0&deleted=1").include_deleted is True  # the box ticked


def test_buckets() -> None:
    week = usage.Window(date(2026, 10, 1), date(2026, 10, 16), "custom")
    assert usage.unit_for(week) == "day" and len(usage.buckets(week, "day")) == 16
    year = usage.Window(date(2025, 10, 17), TODAY, "1y")
    assert usage.unit_for(year) == "week"
    weeks = usage.buckets(year, "week")
    assert weeks[0].start == year.start and weeks[-1].end == TODAY  # clipped to the window
    assert all(b.start.weekday() == 0 for b in weeks[1:])
    long = usage.Window(date(2023, 1, 15), TODAY, "all")
    months = usage.buckets(long, "month")
    assert usage.unit_for(long) == "month" and months[1].start == date(2023, 2, 1)
    assert months[-1] == usage.Bucket(date(2026, 10, 1), TODAY)


def test_formats_and_scales() -> None:
    assert usage.nice_ceiling(17.2) == (20, 5)
    assert usage.nice_ceiling(0.9) == (1.0, 0.25)
    assert [usage.duration(s) for s in (45, 600, 3 * 3600 + 120, 40 * 3600)] == [
        "45 s", "10 min", "3 h 2 min", "40 h",
    ]  # fmt: skip
    assert usage.size(5.4 * 1024 * 1024) == "5.4 MB"
    assert usage.percent(1, 300) == "<1%" and usage.percent(0, 0) == "0%"
    assert usage.stage_of("PHYLOBOT:ALIGN_MAFFT (mafft)") == "align"
    assert usage.stage_of("TREE_IQTREE (mafft-iqtree)") == "trees"


def job(**fields: Any) -> Job:
    """An unsaved Job for the pure computation."""
    fields.setdefault("params", {})
    fields.setdefault("status", Job.Status.SUCCEEDED)
    return Job(**fields)


def at(day: int, hour: int = 12) -> datetime:
    return datetime(2026, 10, day, hour, tzinfo=UTC)


@pytest.mark.django_db
def test_compute() -> None:
    tasks = [
        {"name": "PHYLOBOT:ALIGN_MAFFT (mafft)", "realtime": "10s", "cpus": "2"},
        {"name": "PHYLOBOT:TREE_IQTREE (x)", "realtime": "1m", "cpus": "4"},
    ]
    jobs = [
        job(created_at=at(10), finished_at=at(10, 13), cpu_seconds=260, storage_bytes=1000,
            progress={"tasks": tasks}),
        job(created_at=at(12), finished_at=at(12, 13), cpu_seconds=100, storage_bytes=500,
            deleted_at=at(14), params={"codeml_models": ["M0"], "seqtype": "codon"}),
        job(created_at=at(13), finished_at=at(13, 13), cpu_seconds=40,
            status=Job.Status.FAILED, params={"simulate": True}),
        job(created_at=datetime(2026, 9, 1, tzinfo=UTC), cpu_seconds=999),  # before the window
    ]  # fmt: skip
    result = usage.compute(jobs, filters("range=7d"), NOW)
    tiles = {t.label: t for t in result.tiles}
    assert tiles["CPU time"].value == "6 min" and tiles["CPU time"].change.startswith("none")
    assert tiles["Storage now"].value == "1000 bytes"  # the deleted job's 500 are freed
    assert "500 bytes freed" in tiles["Storage now"].detail
    assert tiles["Jobs submitted"].value == "3"  # 10 to 16 October; not 1 September
    assert tiles["CPU on failed jobs"].value == "40 s"
    # Storage history: held from finish until deletion.
    storage = {p.label: p.value for p in result.storage.points}
    assert storage["Sat 10 Oct 2026"] == "1000 bytes"
    assert storage["Mon 12 Oct 2026"] == "1.5 KB"
    assert storage["Wed 14 Oct 2026"] == "1000 bytes"
    assert [label for label, _, _ in result.kind_table] == ["Ancestors", "dN/dS", "Simulations"]

    whole = usage.compute(jobs, filters("range=30d"), NOW)
    assert [(b.label, b.value) for b in whole.stages] == [
        ("Alignment", "20 s"),
        ("Tree searches and fits", "4 min"),
    ]
    only = usage.compute(jobs, filters("range=30d&kind=simulation"), NOW)
    assert only.jobs_in_window == 1
    live = usage.compute(jobs, filters("range=30d&deleted=0"), NOW)
    assert live.jobs_in_window == 2 and "freed" not in live.tiles[1].detail


@pytest.mark.django_db
def test_my_usage_page(client: Client) -> None:
    assert "No jobs yet" in client.get("/usage/").content.decode()
    owner = local_user()
    mine = Job.objects.create(owner=owner, name="mine", status=Job.Status.SUCCEEDED,
                              cpu_seconds=7200, storage_bytes=2048)  # fmt: skip
    Job.objects.filter(id=mine.id).update(finished_at=mine.created_at + timedelta(minutes=5))
    Job.objects.create(owner=User.objects.create_user("eve"), name="theirs", cpu_seconds=99999)
    page = client.get("/usage/").content.decode()
    assert 'aria-current="page">My Usage</a>' in page  # under Account, in local mode
    assert '<div class="tile-value">2 h</div>' in page and "2.0 KB" in page
    assert ">mine</a>" in page and "theirs" not in page
    assert 'data-tip-label="' in page and "Show as table" in page
    assert '<input type="radio" name="range" value="90d" checked>' in page
    custom = client.get("/usage/?range=custom&from=2020-01-01").content.decode()
    assert 'value="custom" checked' in custom

    # Each panel downloads its figures for the window on the page, as CSV or Excel.
    assert "export=cpu&amp;format=csv" in page and "export=jobs&amp;format=xlsx" in page
    jobs_csv = client.get("/usage/?range=30d&export=jobs&format=csv")
    assert jobs_csv["Content-Type"].startswith("text/csv")
    assert "phylobot_usage_jobs_" in jobs_csv["Content-Disposition"]
    header, row = jobs_csv.content.decode().splitlines()
    assert header.startswith("Job ID,Name,Kind,Status") and "CPU seconds" in header
    assert row.startswith(f"{mine.id},mine,") and ",7200,2 h,2048," in row and row.endswith(",yes")
    cpu_csv = client.get("/usage/?range=30d&export=cpu&format=csv").content.decode().splitlines()
    assert cpu_csv[0] == "Day from,To,CPU seconds,CPU time" and cpu_csv[-1].split(",")[2] == "7200"
    storage = client.get("/usage/?range=30d&export=storage&format=xlsx")
    assert storage.content[:2] == b"PK"  # a workbook
    assert client.get("/usage/?export=nope").status_code == 404


class TestKinds:
    """Which KIND counts a job (`usage.kind_of`).

    The page breaks jobs down by kind with `for key in KINDS: [j for j in inside if
    kind_of(j.params) == key]`, so a kind_of result that is not a KINDS key makes the
    job vanish from the chart and the table with no sign that anything is missing.
    """

    def test_every_job_type_lands_in_a_kind(self) -> None:
        """The guard the dropped jobs needed: whatever job_types grows next, kind_of
        must still answer with a slot this page can draw."""
        cases: list[dict[str, Any]] = [
            {},
            {"stop_after": "alignments"},
            {"stop_after": "trees"},
            {"stop_after": "ancestors"},
            {"codeml_models": ["M0"]},
            {"branch_site": True},
            {"simulate": True},
            {"properties_only": True},
            {"uncertainty_only": True},
            {"esr": True},
            {"reconcile": True},
            {"reconcile": True, "codeml_models": ["M0"]},
        ]
        for params in cases:
            assert usage.kind_of(params) in usage.KINDS, params

    def test_a_properties_follow_up_is_its_own_kind(self) -> None:
        assert usage.kind_of({"properties_only": True}) == "properties"

    def test_reconciling_does_not_make_a_job_a_kind_of_its_own(self) -> None:
        """Reconciliation is something an ancestors job did, not a stage past ancestors.
        It is last in `types_of`, so taking that blindly dropped every reconciled job."""
        assert usage.kind_of({"reconcile": True}) == "ancestors"
        assert usage.kind_of({"reconcile": True, "codeml_models": ["M0"]}) == "dnds"

    def test_each_kind_has_its_own_colour_slot(self) -> None:
        slots = [slot for _, slot in usage.KINDS.values()]
        assert len(set(slots)) == len(slots)


def test_measuring_properties_is_its_own_stage() -> None:
    """PROPERTIES is costed apart from the reconstruction it follows: a neural predictor
    can cost hours per ancestor, and under "ancestors" that would be invisible."""
    from phylobot import workload

    assert usage.stage_of("PHYLOBOT:PROPERTIES (Node1)") == "properties"
    assert usage.stage_of("PROPERTIES_ONLY:PROPERTIES (Node1)") == "properties"
    assert usage.STAGES["properties"] == workload.GROUP_LABELS["properties"]
    assert set(usage.STAGES) == set(workload.GROUP_LABELS)
