"""Homolog searches as records (phylobot.searches): the runner's steps against a fake
service, service and user limits, deletion, and the pages."""

import re
from dataclasses import replace
from datetime import timedelta
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.test import Client, override_settings
from django.utils import timezone

from phylobot import searches
from phylobot.homologs import examples, remote, sequences, sources
from phylobot.homologs.remote import Query, SearchError
from phylobot.models import HomologSearch
from phylobot.modes import local_user
from phylobot.settings import middleware_for
from tests.conftest import FakeStorage

pytestmark = [pytest.mark.django_db, pytest.mark.usefixtures("email")]

QUERY = Query("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK")
Status = HomologSearch.Status


@pytest.fixture
def email(settings: Any) -> None:
    settings.PHYLOBOT_SEARCH_EMAIL = "someone@example.org"


class FakeService:
    """A search service that answers from a script: states to report, then the
    pinned Ime2 result."""

    def __init__(self) -> None:
        self.states: list[str | Exception] = []
        self.submitted: list[tuple[str, str, str]] = []
        self.submit_error: Exception | None = None

    def submit(self, query: Query, email: str, database: str) -> tuple[str, int]:
        if self.submit_error:
            raise self.submit_error
        self.submitted.append((query.residues, email, database))
        return f"remote-{len(self.submitted)}", 10

    def state(self, remote_id: str) -> str:
        answer = self.states.pop(0)
        if isinstance(answer, Exception):
            raise answer
        return answer

    def fetch(self, remote_id: str) -> bytes:
        return examples.raw("ime2-ebi")


@pytest.fixture
def service(monkeypatch: pytest.MonkeyPatch) -> FakeService:
    fake = FakeService()
    source = replace(
        sources.SOURCES["ebi"],
        submit=fake.submit,
        state=fake.state,
        fetch=fake.fetch,
        fetch_sequences=lambda ids, email: dict(examples.sequences("ime2-ebi")),
        max_in_flight=2,
    )
    monkeypatch.setitem(sources.SOURCES, "ebi", source)
    return fake


def new(owner: User | None = None) -> HomologSearch:
    return searches.create(owner or local_user(), QUERY, "ebi", "uniref90", name="test")


def run_due() -> None:
    """One pass of the runner, whatever the leases say."""
    HomologSearch.objects.filter(status__in=searches.ACTIVE).update(next_check_at=timezone.now())
    for search in searches.claim_due():
        searches.advance(search)


def fresh(search: HomologSearch) -> HomologSearch:
    return HomologSearch.objects.get(id=search.id)


# --- runner ------------------------------------------------------------------------


def test_search_runs_to_the_end(fake_storage: FakeStorage, service: FakeService) -> None:
    search = new()
    assert fake_storage.objects[search.query_key] == QUERY.fasta.encode()
    assert search.status == Status.QUEUED and search.query_length == len(QUERY.residues)

    run_due()
    search = fresh(search)
    assert service.submitted == [(QUERY.residues, "someone@example.org", "uniref90")]
    assert (search.status, search.remote_id) == (Status.RUNNING, "remote-1")
    assert search.submitted_at is not None

    service.states = [remote.PENDING, remote.FINISHED]
    run_due()
    assert fresh(search).status == Status.RUNNING
    assert fresh(search).last_checked_at is not None  # shown while it waits
    run_due()
    search = fresh(search)
    assert search.status == Status.SUCCEEDED and search.next_check_at is None
    assert (search.hit_count, search.truncated, search.database_sequences) == (
        1000,
        True,
        121_494_192,
    )
    assert search.database_release.startswith("2026-08-27")
    assert fake_storage.objects[search.result_key] == examples.raw("ime2-ebi")
    assert fake_storage.objects[search.sequences_key].count(b">") == 1000  # full length
    assert len(searches.result(search).hits) == 1000


def test_claims_lease_each_search_once(fake_storage: FakeStorage, service: FakeService) -> None:
    new()
    first = searches.claim_due()
    assert len(first) == 1
    assert searches.claim_due() == []  # leased: another runner skips it
    HomologSearch.objects.update(next_check_at=timezone.now() - timedelta(seconds=1))
    assert len(searches.claim_due()) == 1  # lease expired (a runner died): taken again


def test_service_that_does_not_answer_is_retried_then_failed(
    fake_storage: FakeStorage, service: FakeService
) -> None:
    search = new()
    run_due()
    service.states = [SearchError("timeout")] * (remote.MAX_STATUS_FAILURES + 1)
    for _ in range(remote.MAX_STATUS_FAILURES):
        run_due()
        assert fresh(search).status == Status.RUNNING
    run_due()
    search = fresh(search)
    assert search.status == Status.FAILED and "did not answer" in search.error


def test_failed_at_the_service(fake_storage: FakeStorage, service: FakeService) -> None:
    search = new()
    run_due()
    service.states = [remote.FAILED]
    run_due()
    assert fresh(search).status == Status.FAILED
    assert "failed at EBI" in fresh(search).error


def test_submit_errors_are_retried(fake_storage: FakeStorage, service: FakeService) -> None:
    search = new()
    service.submit_error = SearchError("EBI answered 500")
    run_due()
    search = fresh(search)
    assert (search.status, search.check_failures) == (Status.QUEUED, 1)
    service.submit_error = None
    run_due()
    assert fresh(search).status == Status.RUNNING


def test_service_and_user_limits_hold_searches_in_the_queue(
    fake_storage: FakeStorage, service: FakeService, settings: Any
) -> None:
    other = User.objects.create_user("other")
    mine = [new() for _ in range(2)]
    theirs = new(other)
    run_due()
    # The fake service takes 2 at once: the third waits, whoever owns it.
    assert [fresh(s).status for s in mine] == [Status.RUNNING] * 2
    assert fresh(theirs).status == Status.QUEUED
    settings.PHYLOBOT_SEARCHES_IN_FLIGHT_PER_USER = 1
    service.states = [remote.PENDING] * 5  # the running ones are checked too
    HomologSearch.objects.filter(id=mine[1].id).update(status=Status.QUEUED)
    run_due()
    assert fresh(mine[1]).status == Status.QUEUED  # mine[0] is my one slot
    assert fresh(theirs).status == Status.RUNNING


def test_old_searches_are_abandoned(fake_storage: FakeStorage, service: FakeService) -> None:
    search = new()
    HomologSearch.objects.filter(id=search.id).update(created_at=timezone.now() - timedelta(days=1))
    run_due()
    assert "too long" in fresh(search).error


def test_delete_cancels_and_purges_but_keeps_the_record(
    fake_storage: FakeStorage, service: FakeService
) -> None:
    search = new()
    run_due()
    searches.delete(fresh(search))
    search = fresh(search)
    assert search.status == Status.CANCELLED and search.deleted_at and search.purged_at
    assert not fake_storage.list_keys(search.storage_prefix)
    assert searches.claim_due() == [] and not HomologSearch.objects.for_user(local_user())


def test_deleted_while_fetching_leaves_no_files(
    fake_storage: FakeStorage, service: FakeService, monkeypatch: pytest.MonkeyPatch
) -> None:
    search = new()
    run_due()
    service.states = [remote.FINISHED]

    def fetch_then_deleted(remote_id: str) -> bytes:
        searches.delete(fresh(search))
        return examples.raw("ime2-ebi")

    monkeypatch.setattr(service, "fetch", fetch_then_deleted)
    monkeypatch.setitem(
        sources.SOURCES, "ebi", replace(sources.SOURCES["ebi"], fetch=fetch_then_deleted)
    )
    run_due()
    assert fresh(search).status == Status.CANCELLED
    assert not fake_storage.list_keys(search.storage_prefix)


def test_limits_before_creating(settings: Any) -> None:
    settings.PHYLOBOT_SEARCH_EMAIL = ""
    with pytest.raises(searches.LimitError, match="PHYLOBOT_SEARCH_EMAIL"):
        searches.check_can_search(local_user())


def test_daily_limit(fake_storage: FakeStorage, settings: Any) -> None:
    settings.PHYLOBOT_SEARCHES_PER_DAY = 2
    new(), new()
    with pytest.raises(searches.LimitError, match="2 searches"):
        searches.check_can_search(local_user())


# --- pages -------------------------------------------------------------------------

FORM = {"name": "Ime2", "database": "ebi:uniref90", "consent": "on", "query": QUERY.residues}


def test_start_a_search(client: Client, fake_storage: FakeStorage) -> None:
    response = client.post("/homologs/", {**FORM, "query": f">private name\n{QUERY.residues}"})
    search = HomologSearch.objects.get()
    assert response.status_code == 302 and response["Location"] == f"/homologs/{search.id}/"
    assert (search.name, search.service, search.database) == ("Ime2", "ebi", "uniref90")
    assert b"private name" not in fake_storage.objects[search.query_key]
    page = client.get(f"/homologs/{search.id}/").content.decode()
    assert "queued" in page and 'hx-trigger="every 5s"' in page
    assert search.name in client.get("/homologs/").content.decode()  # listed


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        ({"consent": ""}, "Tick this box to search"),
        ({"query": "MKT"}, "minimum"),
        ({"database": "ebi:uniprotkb"}, "Select a valid choice"),
    ],
)
def test_search_form_refuses(
    client: Client, fake_storage: FakeStorage, fields: dict[str, str], message: str
) -> None:
    response = client.post("/homologs/", {**FORM, **fields})
    assert response.status_code == 200 and message in response.content.decode()
    assert not HomologSearch.objects.exists()


def test_search_refused_without_contact_email(
    client: Client, fake_storage: FakeStorage, settings: Any
) -> None:
    settings.PHYLOBOT_SEARCH_EMAIL = ""
    page = client.post("/homologs/", FORM).content.decode()
    assert "PHYLOBOT_SEARCH_EMAIL" in page and not HomologSearch.objects.exists()


@override_settings(PHYLOBOT_MODE="cloud", MIDDLEWARE=middleware_for("cloud"))
def test_cloud_mode_cannot_choose_ncbi(client: Client, fake_storage: FakeStorage) -> None:
    client.force_login(User.objects.create_user("someone"))
    response = client.post("/homologs/", {**FORM, "database": "ncbi:clustered_nr"})
    assert "Select a valid choice" in response.content.decode()


def test_finished_search_opens_the_review(
    client: Client, fake_storage: FakeStorage, service: FakeService
) -> None:
    search = new()
    run_due()
    service.states = [remote.FINISHED]
    run_due()
    page = client.get(f"/homologs/{search.id}/").content.decode()
    assert re.findall(r'class="count">([\d,]+)<', page) == ["1,000", "1,000", "74", "62", "63"]
    assert f'hx-get="/homologs/{search.id}/"' in page
    panel = client.get(f"/homologs/{search.id}/?min_coverage=55", HTTP_HX_REQUEST="true")
    assert panel.content.decode().lstrip().startswith('<section id="homolog-filters"')
    status = client.get(f"/homologs/{search.id}/status/")
    assert status["HX-Refresh"] == "true"  # finished: the polling page reloads


def test_others_searches_are_not_found(client: Client, fake_storage: FakeStorage) -> None:
    theirs = new(User.objects.create_user("other"))
    for url in (f"/homologs/{theirs.id}/", f"/homologs/{theirs.id}/status/"):
        assert client.get(url).status_code == 404
    for url in (f"/homologs/{theirs.id}/rename/", f"/homologs/{theirs.id}/duplicate/"):
        assert client.post(url, {"name": "mine now"}).status_code == 404
    client.post("/homologs/delete/", {"search": [str(theirs.id)]})
    assert fresh(theirs).deleted_at is None and fresh(theirs).name == "test"  # untouched


def test_delete_one_or_many(client: Client, fake_storage: FakeStorage) -> None:
    a, b, c = new(), new(), new()
    response = client.post("/homologs/delete/", {"search": [str(a.id)]})
    assert response.status_code == 302 and response["Location"] == "/homologs/"
    assert fresh(a).deleted_at is not None and fresh(b).deleted_at is None
    assert client.get(f"/homologs/{a.id}/").status_code == 404
    response = client.post(
        "/homologs/delete/",
        {"search": [str(b.id), str(c.id), "not-a-uuid"], "next": "https://evil.example/"},
    )
    assert response["Location"] == "/homologs/"  # never an outside page
    assert all(fresh(s).deleted_at for s in (b, c))
    assert not HomologSearch.objects.for_user(local_user()).exists()
    assert HomologSearch.objects.count() == 3  # records kept


def test_rename(client: Client, fake_storage: FakeStorage) -> None:
    search = new()
    response = client.post(
        f"/homologs/{search.id}/rename/",
        {"name": "  Ime2   fungi <b>x</b> "},
        HTTP_HX_REQUEST="true",
    )
    assert fresh(search).name == "Ime2 fungi <b>x</b>"
    text = response.content.decode()
    assert text.lstrip().startswith('<div class="editable-title"') and "&lt;b&gt;" in text
    client.post(f"/homologs/{search.id}/rename/", {"name": "x" * 500})
    assert len(fresh(search).name) == searches.MAX_NAME


def finish(search: HomologSearch, fake_storage: FakeStorage) -> HomologSearch:
    fake_storage.put_bytes(search.result_key, examples.raw("ime2-ebi"))
    fake_storage.put_bytes(
        search.sequences_key, sequences.to_fasta(examples.sequences("ime2-ebi")).encode()
    )
    now = timezone.now()
    HomologSearch.objects.filter(id=search.id).update(
        status=Status.SUCCEEDED,
        remote_id="r1",
        submitted_at=now - timedelta(minutes=10),
        finished_at=now,
    )
    return fresh(search)


def test_duplicate_copies_a_finished_search(client: Client, fake_storage: FakeStorage) -> None:
    original = finish(new(), fake_storage)
    response = client.post(f"/homologs/{original.id}/duplicate/")
    copy = HomologSearch.objects.exclude(id=original.id).get()
    assert response["Location"] == f"/homologs/{copy.id}/"
    assert (copy.name, copy.status, copy.copied_from_id) == (
        "test (copy)",
        Status.SUCCEEDED,
        original.id,
    )
    assert (copy.hit_count, copy.service, copy.database) == (1000, "ebi", "uniref90")
    for key in ("query_key", "result_key", "sequences_key"):
        assert (
            fake_storage.objects[getattr(copy, key)] == fake_storage.objects[getattr(original, key)]
        )
    assert "A copy of" in client.get(f"/homologs/{copy.id}/").content.decode()
    # A second copy is numbered; copies send nothing, so they don't count as searches.
    client.post(f"/homologs/{original.id}/duplicate/", {"next": "/homologs/"})
    assert HomologSearch.objects.filter(name="test (copy 2)").exists()
    assert searches.typical_seconds() == {"ebi:uniref90": 600.0}  # the original only
    searches.delete(original)
    assert len(searches.result(copy).hits) == 1000  # the copy keeps its own files


def test_unfinished_searches_cannot_be_duplicated(
    client: Client, fake_storage: FakeStorage
) -> None:
    assert client.post(f"/homologs/{new().id}/duplicate/").status_code == 404


def test_duplicate_a_saved_example(client: Client, fake_storage: FakeStorage) -> None:
    response = client.post("/homologs/examples/ime2-ebi/duplicate/")
    copy = HomologSearch.objects.get()
    assert response["Location"] == f"/homologs/{copy.id}/"
    assert (copy.name, copy.copied_from_example, copy.query_length) == (
        "Ime2 at EBI (copy)",
        "ime2-ebi",
        645,
    )
    assert fake_storage.objects[copy.query_key].startswith(b">query\nMVEKRSRQ")
    assert (
        "Duplicate into your searches"
        in client.get("/homologs/examples/ime2-ebi/").content.decode()
    )
    assert client.post("/homologs/examples/nope/duplicate/").status_code == 404


def test_saved_examples_are_read_only(client: Client) -> None:
    page = client.get("/homologs/examples/ime2-ebi/").content.decode()
    assert "read-only" in page and "data-rename" not in page and "Delete search" not in page


def test_searches_list_has_selection_and_row_actions(
    client: Client, fake_storage: FakeStorage
) -> None:
    done = finish(new(), fake_storage)
    new()
    page = client.get("/homologs/").content.decode()
    assert page.count('<input type="checkbox" name="search"') == 2 and 'id="select-all"' in page
    assert page.count("data-delete=") == 2
    assert page.count('formaction="/homologs/') == 1 and f"/homologs/{done.id}/duplicate/" in page


def test_every_database_choice_is_described() -> None:
    from phylobot.forms import search_choices

    assert {value for value, _ in search_choices("local")} == set(sources.ABOUT)
    assert {a.speed for a in sources.ABOUT.values()} <= {"fast", "slow"}


def test_typical_seconds_is_the_median_per_database(fake_storage: FakeStorage) -> None:
    now = timezone.now()
    for minutes in (10, 16, 30):
        search = new()
        HomologSearch.objects.filter(id=search.id).update(
            status=Status.SUCCEEDED,
            submitted_at=now - timedelta(minutes=minutes),
            finished_at=now,
        )
    new()  # still queued: not counted
    assert searches.typical_seconds() == {"ebi:uniref90": 16 * 60}


def test_search_page_shows_speed_and_recent_times(
    client: Client, fake_storage: FakeStorage
) -> None:
    page = client.get("/homologs/").content.decode()
    assert page.count('class="speed fast"') == 2 and page.count('class="speed slow"') == 3
    assert "Ime2: 16 min" in page and "Ime2: 39 s" in page
    assert "recent searches here" not in page
    search = new()
    now = timezone.now()
    HomologSearch.objects.filter(id=search.id).update(
        status=Status.SUCCEEDED, submitted_at=now - timedelta(minutes=12), finished_at=now
    )
    assert "recent searches here: typically 12 min 0 s" in client.get("/homologs/").content.decode()


def test_search_finished_before_sequences_were_fetched(fake_storage: FakeStorage) -> None:
    search = new()
    assert searches.full_sequences(search) == {}  # no sequences.fasta: none, not an error


def test_runner_survives_database_errors(monkeypatch: pytest.MonkeyPatch) -> None:
    """A runner started before a migration must wait, not exit."""
    from phylobot.management.commands import run_searches

    calls = []

    def claim() -> list[HomologSearch]:
        calls.append(1)
        if len(calls) == 1:
            raise RuntimeError("column does not exist")
        raise KeyboardInterrupt  # stop the loop

    slept: list[float] = []
    monkeypatch.setattr(searches, "claim_due", claim)
    monkeypatch.setattr("phylobot.management.commands.run_searches.time.sleep", slept.append)
    monkeypatch.setattr("phylobot.management.commands.run_searches.signal.signal", lambda *a: None)
    with pytest.raises(KeyboardInterrupt):
        run_searches.Command().handle()
    assert len(calls) == 2 and slept == [run_searches.ERROR_SECONDS]
