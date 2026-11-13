"""Homolog searches interleaved with jobs on My Jobs (phylobot.worklist, Approach A)."""

from datetime import UTC, datetime

import pytest
from django.contrib.auth.models import User
from django.test import Client

from phylobot.models import HomologSearch, Job
from phylobot.modes import local_user

pytestmark = pytest.mark.django_db


def make(owner: User) -> tuple[Job, HomologSearch]:
    """A job and a search for the user, the search the newer of the two."""
    job = Job.objects.create(owner=owner, name="MT-CO1 tree", params={"stop_after": "trees"})
    Job.objects.filter(id=job.id).update(created_at=datetime(2026, 11, 8, 9, tzinfo=UTC))
    search = HomologSearch.objects.create(
        owner=owner, name="MT-CO1 homologs", service="ebi", database="uniref90",
        query_length=500, status=HomologSearch.Status.SUCCEEDED,
    )  # fmt: skip
    HomologSearch.objects.filter(id=search.id).update(
        created_at=datetime(2026, 11, 9, 9, tzinfo=UTC)
    )
    return job, search


def page(client: Client, query: str = "") -> str:
    return client.get(f"/jobs/?{query}").content.decode()


def test_a_search_appears_on_my_jobs(client: Client) -> None:
    job, search = make(local_user())
    html = page(client)
    # Both are listed; the search links to its own page and wears a "Homolog search" chip.
    assert f'href="/jobs/{job.id}/"' in html
    assert f'href="/homologs/{search.id}/"' in html
    assert "MT-CO1 homologs" in html and "Homolog search" in html
    # Newest first by default: the search (newer) comes before the job.
    assert html.index("MT-CO1 homologs") < html.index("MT-CO1 tree")
    # Not bulk-selectable (no job checkbox for it), no compute shown.
    assert f'name="job" value="{search.id}"' not in html


def test_the_homolog_type_filter(client: Client) -> None:
    job, search = make(local_user())
    only = page(client, "type=homolog")  # searches alone
    assert "MT-CO1 homologs" in only and "MT-CO1 tree" not in only
    jobs_only = page(client, "type=trees")  # a real job type hides searches
    assert "MT-CO1 tree" in jobs_only and "MT-CO1 homologs" not in jobs_only
    both = page(client, "type=trees&type=homolog")  # together
    assert "MT-CO1 tree" in both and "MT-CO1 homologs" in both


def test_job_only_filters_hide_searches(client: Client) -> None:
    make(local_user())
    assert "MT-CO1 homologs" not in page(client, "data=protein")  # a job-only filter
    assert "MT-CO1 homologs" not in page(client, "origin=real")
    # A shared filter (status) still applies to searches.
    assert "MT-CO1 homologs" in page(client, "status=succeeded")
    assert "MT-CO1 homologs" not in page(client, "status=running")


def test_searches_count_toward_the_total(client: Client) -> None:
    make(local_user())  # one job + one search
    assert "2 jobs" in page(client)


def test_search_rows_have_select_duplicate_and_delete(client: Client) -> None:
    job, search = make(local_user())
    html = page(client)
    assert f'name="search" value="{search.id}"' in html  # bulk-selectable
    assert f'data-delete="{search.id}"' in html  # its own delete
    assert f'id="dup-{search.id}"' in html and 'form="dup-' in html  # duplicate (succeeded)
    # A running search cannot be duplicated: no duplicate form for it.
    HomologSearch.objects.filter(id=search.id).update(status=HomologSearch.Status.RUNNING)
    assert f'id="dup-{search.id}"' not in page(client)


@pytest.mark.usefixtures("fake_storage")
def test_delete_jobs_endpoint_deletes_searches_too(client: Client) -> None:
    job, search = make(local_user())
    response = client.post("/jobs/delete/", {"job": str(job.id), "search": str(search.id)})
    assert response.status_code == 302
    assert not Job.objects.for_user(local_user()).filter(id=job.id).exists()
    assert not HomologSearch.objects.for_user(local_user()).filter(id=search.id).exists()
    # Records are kept (soft delete), as for jobs.
    assert HomologSearch.objects.get(id=search.id).deleted_at is not None


# --- Grouped view: follow-up jobs under the job they were run on ---------------------


def family(owner: User) -> dict[str, Job]:
    """A job, three kinds of follow-up on it, a follow-up of one of those, and a stranger."""
    parent = Job.objects.create(owner=owner, name="Parent", status=Job.Status.SUCCEEDED)
    pid = str(parent.id)
    esr = Job.objects.create(
        owner=owner, name="ESR on parent", status=Job.Status.RUNNING,
        params={"esr": True, "esr_of": {"kind": "job", "job": pid}},
    )
    unc = Job.objects.create(
        owner=owner, name="Uncertainty of parent", status=Job.Status.SUCCEEDED,
        params={"uncertainty_of": {"parent": pid}},
    )
    extension = Job.objects.create(owner=owner, name="Extension", status=Job.Status.SUCCEEDED,
                                   extends=parent)  # fmt: skip
    grandchild = Job.objects.create(
        owner=owner, name="Properties of the extension", status=Job.Status.FAILED,
        params={"properties_of": {"parent": str(extension.id)}},
    )
    orphan = Job.objects.create(
        owner=owner, name="ESR on a job not listed", status=Job.Status.SUCCEEDED,
        params={"esr_of": {"kind": "job", "job": "00000000-0000-0000-0000-000000000000"}},
    )
    alone = Job.objects.create(owner=owner, name="Alone", status=Job.Status.SUCCEEDED,
                               rerun_of=parent)  # fmt: skip
    return {"parent": parent, "esr": esr, "unc": unc, "extension": extension,
            "grandchild": grandchild, "orphan": orphan, "alone": alone}  # fmt: skip


def test_nest_moves_follow_ups_under_their_job() -> None:
    from phylobot import worklist

    jobs = family(local_user())
    assert worklist.parent_id(jobs["esr"]) == str(jobs["parent"].id)
    assert worklist.parent_id(jobs["extension"]) == str(jobs["parent"].id)
    assert worklist.parent_id(jobs["alone"]) == ""  # a duplicate is a job of its own
    order = ["alone", "orphan", "grandchild", "extension", "unc", "esr", "parent"]
    rows = [worklist.job_row(jobs[k], [], "", None, []) for k in order]
    nested = worklist.nest(rows)
    by_name = {row.name: row for row in nested}
    # Groups keep the table's order, follow-ups theirs; a follow-up of a follow-up goes to
    # the first job of the chain; one whose job is not listed stays where it was.
    assert [row.name for row in nested] == [
        "Alone", "ESR on a job not listed", "Parent", "Properties of the extension",
        "Extension", "Uncertainty of parent", "ESR on parent",
    ]  # fmt: skip
    parent = by_name["Parent"]
    assert parent.follow_ups == 4 and parent.follow_up_note == "1 running · 1 failed"
    under = ("Properties of the extension", "Extension", "Uncertainty of parent", "ESR on parent")
    assert all(by_name[n].under == parent.id for n in under)
    assert by_name["Alone"].under == "" and by_name["ESR on a job not listed"].under == ""


def test_grouping_is_a_toggle_off_by_default(client: Client) -> None:
    family(local_user())
    flat = page(client)
    assert 'role="switch" aria-checked="false"' in flat and "follow-up-row" not in flat
    assert "follow-up job" not in flat.split("<tbody>")[1]
    # The switch's link turns it on, is remembered, and leaves the address as it was.
    response = client.get("/jobs/", {"group": "on", "status": "succeeded"})
    assert response.status_code == 302 and response["Location"] == "/jobs/?status=succeeded"
    grouped = page(client)
    assert 'role="switch" aria-checked="true"' in grouped
    assert grouped.count('class="follow-up-row"') == 4 and "4 follow-up jobs" in grouped
    assert "1 running · 1 failed" in grouped
    # Follow-ups come folded, so a bulk action cannot reach one the reader has not opened.
    folded = grouped.split('<tr class="follow-up-row"')[1:]
    assert all(" hidden>" in row.split("\n")[0] for row in folded)
    assert client.get("/jobs/", {"group": "off"}).status_code == 302
    assert "follow-up-row" not in page(client)
