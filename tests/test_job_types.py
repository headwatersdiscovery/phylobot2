"""Job types and the My Jobs filters (phylobot.job_types)."""

import re
from datetime import UTC, datetime
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.test import Client

from phylobot.job_types import types_of
from phylobot.models import Job
from phylobot.modes import local_user
from tests.conftest import FakeStorage

pytestmark = pytest.mark.django_db

KINDS: dict[str, dict[str, Any]] = {
    "full": {"stop_after": "ancestors", "seqtype": "protein"},
    "old": {},  # before stop_after existed: a full analysis
    "aligned": {"stop_after": "alignments"},
    "treed": {"stop_after": "trees"},
    "dnds": {"stop_after": "ancestors", "seqtype": "codon", "codeml_models": ["M0"]},
    "branchsite": {"branch_site": True, "seqtype": "codon"},
    "simulation": {"simulate": True, "spec": {}},
    "on-sim": {"stop_after": "trees", "simulated": True},
}


def keys(params: dict[str, Any]) -> list[str]:
    return [t.key for t in types_of(params)]


def test_types_of() -> None:
    assert keys(KINDS["full"]) == keys(KINDS["old"]) == ["alignments", "trees", "ancestors"]
    assert keys(KINDS["aligned"]) == ["alignments"]
    assert keys(KINDS["treed"]) == ["alignments", "trees"]
    assert keys(KINDS["dnds"]) == ["alignments", "trees", "ancestors", "dnds"]
    assert keys(KINDS["branchsite"]) == ["dnds"]
    assert keys(KINDS["simulation"]) == ["simulation"]


def shown(client: Client, query: str) -> set[str]:
    page = client.get(f"/jobs/?{query}").content.decode()
    return {name for name in KINDS if f">{name}</a>" in page}


def test_my_jobs_filters(client: Client) -> None:
    owner = local_user()
    for name, params in KINDS.items():
        Job.objects.create(owner=owner, name=name, params=params)
    Job.objects.create(owner=User.objects.create_user("bob"), name="full", params={})
    finished = datetime(2026, 10, 1, 12, tzinfo=UTC)
    Job.objects.filter(name="dnds").update(status=Job.Status.SUCCEEDED, finished_at=finished)
    Job.objects.get(name="treed").add_tags(["review"])

    assert shown(client, "") == set(KINDS)
    analyses = set(KINDS) - {"branchsite", "simulation"}
    assert shown(client, "type=alignments") == analyses
    assert shown(client, "type=trees") == analyses - {"aligned"}
    assert shown(client, "type=ancestors") == {"full", "old", "dnds"}
    assert shown(client, "type=dnds") == {"dnds", "branchsite"}
    assert shown(client, "type=dnds&type=simulation") == {"dnds", "branchsite", "simulation"}
    assert shown(client, "type=simulation") == {"simulation"}
    assert shown(client, "type=bogus") == set(KINDS)  # unknown values are dropped
    assert shown(client, "status=succeeded") == {"dnds"}
    assert shown(client, "status=queued&type=simulation") == {"simulation"}
    assert shown(client, "data=codon") == {"dnds", "branchsite"}
    assert "dnds" not in shown(client, "data=protein") and "old" in shown(client, "data=protein")
    assert shown(client, "origin=simulated") == {"simulation", "on-sim"}
    assert shown(client, "origin=real") == set(KINDS) - {"simulation", "on-sim"}
    assert shown(client, "tag=review") == {"treed"}
    assert shown(client, "q=TREE") == {"treed"}
    assert shown(client, "finished_from=2026-10-01&finished_to=2026-10-01") == {"dnds"}
    assert shown(client, "finished_from=2026-10-02") == set()
    assert shown(client, "finished_to=not-a-date") == set(KINDS)
    dnds = Job.objects.get(name="dnds", owner=owner)
    assert shown(client, f"q={str(dnds.id)[:8]}") == {"dnds"}

    page = client.get("/jobs/?type=ancestors&status=queued").content.decode()
    assert "Showing 2 of 8 jobs" in page and '<span class="filter-count">2</span>' in page
    assert "No jobs match these filters" in client.get("/jobs/?q=zzz").content.decode()


def test_my_jobs_type_column(client: Client) -> None:
    Job.objects.create(owner=local_user(), name="x", params=KINDS["dnds"])
    page = client.get("/jobs/").content.decode()
    for label in ("Alignment", "Trees", "Ancestors", "dN/dS"):
        assert f">{label}</a>" in page
    assert 'href="/jobs/?type=dnds"' in page


def order(client: Client, query: str) -> list[str]:
    page = client.get(f"/jobs/?{query}").content.decode()
    names = re.findall(r'data-name="([^"]*)"', page)
    return names


def test_my_jobs_sorting(client: Client) -> None:
    owner = local_user()
    jobs = {}
    for name, params in (("beta", KINDS["simulation"]), ("Alpha", KINDS["aligned"]),
                         ("", KINDS["full"]), ("gamma", KINDS["dnds"])):  # fmt: skip
        jobs[name] = Job.objects.create(owner=owner, name=name, params=params)
    when = datetime(2026, 10, 1, 12, tzinfo=UTC)
    Job.objects.filter(id=jobs["beta"].id).update(status=Job.Status.FAILED, finished_at=when)
    Job.objects.filter(id=jobs["gamma"].id).update(
        status=Job.Status.SUCCEEDED, finished_at=when.replace(day=2)
    )
    jobs["gamma"].add_tags(["zeta"])
    jobs["Alpha"].add_tags(["eta"])
    unnamed = "(unnamed)"

    assert order(client, "") == ["gamma", unnamed, "Alpha", "beta"]  # newest first
    assert order(client, "sort=submitted&dir=asc") == ["beta", "Alpha", unnamed, "gamma"]
    assert order(client, "sort=name") == ["Alpha", "beta", "gamma", unnamed]  # unnamed last
    assert order(client, "sort=name&dir=desc") == ["gamma", "beta", "Alpha", unnamed]
    assert order(client, "sort=type") == ["Alpha", unnamed, "gamma", "beta"]
    assert order(client, "sort=tags") == ["Alpha", "gamma", unnamed, "beta"]
    assert order(client, "sort=status") == [unnamed, "Alpha", "gamma", "beta"]
    assert order(client, "sort=finished") == ["gamma", "beta", unnamed, "Alpha"]  # newest first
    assert order(client, "sort=finished&dir=asc") == ["beta", "gamma", unnamed, "Alpha"]
    assert order(client, "sort=bogus") == order(client, "")

    page = client.get("/jobs/?sort=name&type=simulation").content.decode()
    assert 'aria-sort="ascending"' in page and "?type=simulation&amp;sort=name&amp;dir=desc" in page
    assert '<input type="hidden" name="sort" value="name">' in page  # the filters keep it
    assert 'href="/jobs/?sort=name&amp;dir=asc">Clear filters' in page


def test_deleted_jobs_keep_their_types(client: Client) -> None:
    from phylobot import jobs

    job = Job.objects.create(owner=local_user(), name="gone", params=KINDS["branchsite"])
    jobs.delete([job])
    deleted = client.get("/jobs/").content.decode().split("Deleted Jobs")[1]
    assert '<span class="job-type codeml" title=' in deleted and ">dN/dS</span>" in deleted
    assert "?type=" not in deleted  # the filters cover live jobs only: no links


def test_storage_is_measured_once_and_kept_after_purge(
    client: Client, fake_storage: FakeStorage
) -> None:
    from phylobot import jobs

    owner = local_user()
    done = Job.objects.create(
        owner=owner, name="done", status=Job.Status.SUCCEEDED, cpu_seconds=259
    )
    queued = Job.objects.create(owner=owner, name="queued")
    fake_storage.objects[f"{done.storage_prefix}inputs/sequences.fasta"] = b"x" * 1000
    fake_storage.objects[f"{done.storage_prefix}results/a.fasta"] = b"x" * 2048
    fake_storage.objects["jobs/other/results/a.fasta"] = b"x" * 99

    assert jobs.measure_unmeasured() == 1  # finished jobs only
    done.refresh_from_db()
    assert done.storage_bytes == 3048 and queued.storage_bytes is None
    assert jobs.measure_unmeasured() == 0  # once
    row = client.get("/jobs/").content.decode().split(">done</a>")[1].split("</tr>")[0]
    assert "4 min 19 s" in row and "3.0\xa0KB" in row

    jobs.delete([done])
    done.refresh_from_db()
    assert done.storage_bytes == 3048 and done.purged_at is not None  # the record stays
    assert not jobs.measure_storage(done)  # never re-measured once its files are gone


def test_sort_by_usage(client: Client) -> None:
    owner = local_user()
    for name, cpu, size in (("small", 10, 5_000_000), ("big", 999, 1_000), ("none", None, None)):
        Job.objects.create(owner=owner, name=name, cpu_seconds=cpu, storage_bytes=size)
    assert order(client, "sort=cpu") == ["big", "small", "none"]  # largest first, unknown last
    assert order(client, "sort=disk") == ["small", "big", "none"]
    assert order(client, "sort=disk&dir=asc") == ["big", "small", "none"]
    page = client.get("/jobs/?sort=disk").content.decode()
    assert 'class="sortable usage multi" aria-sort="descending"' in page


def test_reconciled_jobs_get_their_own_type_and_filter(client: Client) -> None:
    from phylobot.job_types import composition

    reconciled = {"stop_after": "ancestors", "reconcile": True}
    assert keys(reconciled) == ["alignments", "trees", "ancestors", "reconciled"]
    assert keys({"stop_after": "trees", "reconcile_skipped": "no consent"}) == [
        "alignments",
        "trees",
    ]  # fmt: skip  (asked for, but ran without it: no pill)
    for name, params in (("recon", reconciled), ("plain", KINDS["full"])):
        Job.objects.create(owner=local_user(), name=name, params=params)
    page = client.get("/jobs/?type=reconciled").content.decode()
    assert ">recon</a>" in page and ">plain</a>" not in page
    assert composition({"simulate": True}) is None


def test_my_jobs_says_what_each_job_was_made_of(client: Client, fake_storage: FakeStorage) -> None:
    from phylobot.job_types import composition

    params = {
        "stop_after": "ancestors", "aligners": ["mafft", "muscle"], "alignment_labels": ["mine"],
        "modelfinder": True, "fixed_models": True, "models": ["LG+G4", "WAG+G4"],
        "tree_builders": ["iqtree", "generax"], "sequence_count": 7,
    }  # fmt: skip
    made = composition(params)
    assert made is not None
    assert made.text == "7 sequences · 3 alignments × 3 models · 2 tree sources"
    assert "Alignments: MAFFT, MUSCLE, your mine" in made.details
    assert "Models: ModelFinder, LG+G4, WAG+G4" in made.details
    aligned_only = composition({"stop_after": "alignments", "aligners": ["mafft"]})
    assert aligned_only is not None and aligned_only.text == "1 alignment"
    # an older job, from before the count was kept: counted once from its input, then kept
    old = Job.objects.create(owner=local_user(), name="old", params={"aligners": ["mafft"]})
    fake_storage.objects[f"{old.storage_prefix}inputs/sequences.fasta"] = b">a\nMK\n>b\nMK\n"
    page = client.get("/jobs/").content.decode()
    line = '<span class="phrase">2 sequences</span> · <span class="phrase">1 alignment</span>'
    assert line in page
    old.refresh_from_db()
    assert old.params["sequence_count"] == 2
