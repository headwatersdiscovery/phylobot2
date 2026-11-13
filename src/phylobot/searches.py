"""Homolog searches: create them, follow them at the remote service, delete them.

A search is not a pipeline Job (owner decision 2026-09-29): it uses no compute here and
mostly waits on EBI or NCBI, for up to a quarter of an hour. So instead of a thread per
search, one runner (manage.py run_searches) repeatedly claims the searches that are due
and moves each one step on:

    queued --submit--> running --check ... check--> succeeded (result and the hits'
                                                    full sequences stored)
                                                \\-> failed

- Claiming pushes `next_check_at` forward by a lease under SELECT ... FOR UPDATE SKIP
  LOCKED, so several runners never take the same search, and a runner that dies leaves
  its searches to be picked up again when the lease runs out.
- Status only moves forward: each update is conditional on the status it started from,
  so a search cancelled (deleted) meanwhile stays cancelled.
- Service limits: at most `max_in_flight` searches waiting at each service across the
  site (EBI asks for batches of at most 30), and a per-user number (site_settings);
  others wait their turn in the queue.
- A service that does not answer is retried; after remote.MAX_STATUS_FAILURES failures
  in a row the search fails.
- Deleting a search removes its files and keeps the record, as for jobs.
"""

import logging
import statistics
from datetime import timedelta
from typing import Any

from django.contrib.auth.models import AbstractBaseUser
from django.db import transaction
from django.utils import timezone

from phylobot import site_settings, storage
from phylobot.homologs import remote, sources
from phylobot.homologs import sequences as hit_sequences
from phylobot.homologs.hits import ParseError, SearchResult
from phylobot.homologs.remote import Query, SearchError
from phylobot.models import HomologSearch

logger = logging.getLogger(__name__)

Status = HomologSearch.Status
ACTIVE = (Status.QUEUED, Status.RUNNING)
LEASE = timedelta(minutes=5)  # a claimed search is left alone this long
WAIT_FOR_SLOT = timedelta(seconds=30)  # a queued search waiting for a service slot


class LimitError(Exception):
    """The user may not start another search now (shown to them)."""


def check_can_search(owner: AbstractBaseUser) -> None:
    """Raise LimitError if the site cannot search, or this user has searched enough today."""
    if not site_settings.search_email():
        raise LimitError(
            "Searching is not set up: EBI and NCBI require a contact email. Enter one "
            "on the Settings page (or set PHYLOBOT_SEARCH_EMAIL in .env)."
        )
    since = timezone.now() - timedelta(days=1)
    today = (
        HomologSearch.objects.filter(owner_id=owner.pk, created_at__gte=since)
        .filter(copied_from__isnull=True, copied_from_example="")  # copies send nothing
        .count()
    )
    if today >= site_settings.searches_per_day():
        raise LimitError(
            f"You have started {today} searches in the last 24 hours, the most allowed. "
            "Try again later."
        )


def create(
    owner: AbstractBaseUser, query: Query, service: str, database: str, name: str = ""
) -> HomologSearch:
    """Queue a search. The query is stored in object storage, not the database."""
    source = sources.get(service)
    if database not in source.databases:
        raise SearchError(f"unknown database {database!r} for {service}")
    with transaction.atomic():
        search = HomologSearch.objects.create(
            owner_id=owner.pk,
            name=name,
            service=service,
            database=database,
            query_length=len(query.residues),
            next_check_at=timezone.now(),
        )
        storage.put_bytes(search.query_key, query.fasta.encode())
    return search


def typical_seconds(recent: int = 20) -> dict[str, float]:
    """Median time at the service (submitted to finished) of the last `recent` finished
    searches per "service:database", across the site, for the search page."""
    finished = (
        HomologSearch.objects.filter(
            status=Status.SUCCEEDED, submitted_at__isnull=False, finished_at__isnull=False
        )
        .filter(copied_from__isnull=True, copied_from_example="")
        .order_by("-finished_at")
        .values_list("service", "database", "submitted_at", "finished_at")[:500]
    )
    times: dict[str, list[float]] = {}
    for service, database, submitted, done in finished:
        if submitted is None or done is None:  # excluded by the query; for the types
            continue
        key = f"{service}:{database}"
        if len(times.setdefault(key, [])) < recent:
            times[key].append((done - submitted).total_seconds())
    return {key: statistics.median(values) for key, values in times.items()}


def claim_due(limit: int = 20) -> list[HomologSearch]:
    """Take up to `limit` active searches whose next check is due, leasing each one."""
    now = timezone.now()
    with transaction.atomic():
        due = list(
            HomologSearch.objects.select_for_update(skip_locked=True)
            .filter(status__in=ACTIVE, next_check_at__lte=now, deleted_at__isnull=True)
            .order_by("next_check_at")[:limit]
        )
        ids = [search.id for search in due]
        HomologSearch.objects.filter(id__in=ids).update(next_check_at=now + LEASE)
    return due


def _update(search: HomologSearch, expected: str, **fields: Any) -> bool:
    """Apply `fields` if the search is still `expected` (forward-only); False otherwise."""
    changed = HomologSearch.objects.filter(
        id=search.id, status=expected, deleted_at__isnull=True
    ).update(**fields)
    if changed:
        for name, value in fields.items():
            setattr(search, name, value)
    return changed == 1


def _fail(search: HomologSearch, message: str) -> None:
    logger.info("search %s failed: %s", search.id, message)
    _update(
        search,
        search.status,
        status=Status.FAILED,
        error=message,
        finished_at=timezone.now(),
        next_check_at=None,
    )


def _retry_later(search: HomologSearch, exc: Exception, delay: int) -> None:
    """A service that did not answer: try again, unless it has failed too often."""
    failures = search.check_failures + 1
    if failures > remote.MAX_STATUS_FAILURES:
        _fail(search, f"The search service did not answer: {exc}")
        return
    logger.info("search %s: retry %d after: %s", search.id, failures, exc)
    _update(
        search,
        search.status,
        check_failures=failures,
        next_check_at=timezone.now() + timedelta(seconds=delay),
    )


def advance(search: HomologSearch) -> None:
    """Move one claimed search one step on. Never raises for service trouble."""
    source = sources.get(search.service)
    now = timezone.now()
    if now - search.created_at > site_settings.search_time_limit():
        _fail(search, "The search took too long at the service and was abandoned.")
        return
    try:
        if search.status == Status.QUEUED:
            _submit(search, source)
        elif search.status == Status.RUNNING:
            _check(search, source)
    except SearchError as exc:
        _retry_later(search, exc, source.poll_seconds)
    except ParseError as exc:
        _fail(search, f"The service's answer could not be read: {exc}")


def _submit(search: HomologSearch, source: sources.Source) -> None:
    waiting = HomologSearch.objects.filter(status=Status.RUNNING, deleted_at__isnull=True)
    if (
        waiting.filter(service=search.service).count() >= source.max_in_flight
        or waiting.filter(owner_id=search.owner_id).count()
        >= site_settings.searches_in_flight_per_user()
    ):
        _update(search, Status.QUEUED, next_check_at=timezone.now() + WAIT_FOR_SLOT)
        return
    query = remote.clean_query(storage.read_bytes(search.query_key).decode())
    remote_id, first_wait = source.submit(query, site_settings.search_email(), search.database)
    now = timezone.now()
    logger.info("search %s: submitted to %s as %s", search.id, search.service, remote_id)
    _update(
        search,
        Status.QUEUED,
        status=Status.RUNNING,
        remote_id=remote_id,
        submitted_at=now,
        check_failures=0,
        next_check_at=now + timedelta(seconds=first_wait),
    )


def _check(search: HomologSearch, source: sources.Source) -> None:
    state = source.state(search.remote_id)
    logger.info("search %s: %s says %s", search.id, search.service, state)
    if state == remote.PENDING:
        _update(
            search,
            Status.RUNNING,
            check_failures=0,
            last_checked_at=timezone.now(),
            next_check_at=timezone.now() + timedelta(seconds=source.poll_seconds),
        )
    elif state == remote.FAILED:
        _fail(search, f"The search failed at {source.label}.")
    else:
        content = source.fetch(search.remote_id)
        result = source.parse(content)
        # Whole proteins for jobs and downloads; a failure here retries the step.
        full = source.fetch_sequences([h.id for h in result.hits], site_settings.search_email())
        storage.put_bytes(search.result_key, content)
        storage.put_bytes(search.sequences_key, hit_sequences.to_fasta(full).encode())
        stored = _update(
            search,
            Status.RUNNING,
            status=Status.SUCCEEDED,
            finished_at=timezone.now(),
            next_check_at=None,
            check_failures=0,
            program=result.program[:40],
            database_release=result.database_created[:40],
            database_sequences=result.database_sequences,
            hit_count=len(result.hits),
            truncated=result.truncated,
        )
        if not stored:  # deleted while we fetched: its files were already purged
            storage.delete_prefix(search.storage_prefix)
            return
        logger.info("search %s: %d hits", search.id, len(result.hits))


def result(search: HomologSearch) -> SearchResult:
    """The finished search's hits, read from its stored result."""
    return sources.get(search.service).parse(storage.read_bytes(search.result_key))


def full_sequences(search: HomologSearch) -> dict[str, str]:
    """The finished search's hits' full-length sequences, {hit ID: sequence}; empty for a
    search finished before sequences were fetched."""
    if search.sequences_key not in storage.list_keys(search.storage_prefix):
        return {}
    return hit_sequences.parse_fasta(storage.read_bytes(search.sequences_key).decode())


MAX_NAME = 100


def rename(search: HomologSearch, name: str) -> None:
    """Rename a search (its owner's label; blank is allowed and shows as "Search")."""
    search.name = " ".join(name.split())[:MAX_NAME]
    HomologSearch.objects.filter(id=search.id).update(name=search.name)


def copy_name(owner: AbstractBaseUser, name: str) -> str:
    """ "<name> (copy)", or "(copy 2)", "(copy 3)"... if the owner already has that name."""
    base = (name or "Search")[: MAX_NAME - 12]
    taken = set(
        HomologSearch.objects.for_user(owner)
        .filter(name__startswith=f"{base} (copy")
        .values_list("name", flat=True)
    )
    candidate, n = f"{base} (copy)", 2
    while candidate in taken:
        candidate, n = f"{base} (copy {n})", n + 1
    return candidate


def _copy(
    owner: AbstractBaseUser,
    name: str,
    service: str,
    database: str,
    files: dict[str, bytes],
    result: SearchResult,
    **provenance: Any,
) -> HomologSearch:
    """A finished search of `owner`'s holding these files. Nothing is sent anywhere."""
    now = timezone.now()
    with transaction.atomic():
        copy = HomologSearch.objects.create(
            owner_id=owner.pk,
            name=copy_name(owner, name),
            status=Status.SUCCEEDED,
            service=service,
            database=database,
            query_length=result.query_length,
            finished_at=now,
            program=result.program[:40],
            database_release=result.database_created[:40],
            database_sequences=result.database_sequences,
            hit_count=len(result.hits),
            truncated=result.truncated,
            **provenance,
        )
        for attribute, content in files.items():
            storage.put_bytes(getattr(copy, attribute), content)
    return copy


def duplicate(search: HomologSearch, owner: AbstractBaseUser) -> HomologSearch:
    """Copy a finished search (query, hits, sequences) as a new search of `owner`'s."""
    if search.status != Status.SUCCEEDED:
        raise SearchError("only a finished search can be duplicated")
    files = {"query_key": storage.read_bytes(search.query_key),
             "result_key": storage.read_bytes(search.result_key)}  # fmt: skip
    if search.sequences_key in storage.list_keys(search.storage_prefix):
        files["sequences_key"] = storage.read_bytes(search.sequences_key)
    copy = _copy(owner, search.name, search.service, search.database, files,
                 result(search), copied_from=search)  # fmt: skip
    HomologSearch.objects.filter(id=copy.id).update(remote_id=search.remote_id)
    return copy


def duplicate_example(key: str, owner: AbstractBaseUser) -> HomologSearch:
    """Copy a saved example into `owner`'s own searches (so they can rename or delete it)."""
    from phylobot.homologs import examples

    example = examples.EXAMPLES[key]
    files = {
        "query_key": Query(examples.query(key)).fasta.encode(),
        "result_key": examples.raw(key),
        "sequences_key": hit_sequences.to_fasta(examples.sequences(key)).encode(),
    }
    return _copy(owner, example.label, example.source, example.database, files,
                 examples.load(key), copied_from_example=key)  # fmt: skip


def delete(search: HomologSearch) -> None:
    """Delete for its owner: cancel it if unfinished, remove its files, keep the record."""
    now = timezone.now()
    fields: dict[str, Any] = {"deleted_at": now, "next_check_at": None}
    if search.status in ACTIVE:
        fields |= {"status": Status.CANCELLED, "finished_at": now}
    HomologSearch.objects.filter(id=search.id, deleted_at__isnull=True).update(**fields)
    storage.delete_prefix(search.storage_prefix)
    HomologSearch.objects.filter(id=search.id).update(purged_at=timezone.now())
