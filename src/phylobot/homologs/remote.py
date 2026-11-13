"""What every remote homolog search shares: the query checks, HTTP, and waiting.

A remote search sends the query's residues to a public service, so it is opt-in. Only
the residues are sent, under the fixed name "query"; the user's sequence name never is.
Treat the query as hostile: it is checked here before anything leaves the machine.
"""

import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass

MIN_LENGTH = 20
MAX_LENGTH = 5000
RESIDUES = set("ACDEFGHIKLMNPQRSTVWYBZJX")
EMAIL = re.compile(r"^[^@\s]{1,64}@[^@\s]{1,255}\.[^@\s]+$")
TIMEOUT = 60  # seconds without any data before a request fails
TOTAL_TIMEOUT = 300  # seconds for a whole answer, however slowly it arrives
CHUNK = 1 << 16
MAX_STATUS_FAILURES = 5  # consecutive unanswered status checks before giving up

# A service's job state, as each client reports it.
PENDING, FINISHED, FAILED = "pending", "finished", "failed"


class SearchError(Exception):
    pass


@dataclass(frozen=True)
class Query:
    residues: str

    @property
    def fasta(self) -> str:
        return f">query\n{self.residues}\n"


def clean_query(text: str) -> Query:
    """One protein sequence, raw or FASTA, as a Query; raises SearchError otherwise."""
    lines = [line.strip() for line in text.strip().splitlines() if line.strip()]
    if sum(line.startswith(">") for line in lines) > 1:
        raise SearchError("give one query sequence, not several")
    residues = "".join(line for line in lines if not line.startswith(">"))
    residues = re.sub(r"\s+", "", residues).upper().replace("-", "").rstrip("*")
    if not residues:
        raise SearchError("the query sequence is empty")
    bad = sorted(set(residues) - RESIDUES)
    if bad:
        raise SearchError(f"the query contains characters that are not amino acids: {''.join(bad)}")
    if len(residues) < MIN_LENGTH:
        raise SearchError(f"the query is {len(residues)} residues; the minimum is {MIN_LENGTH}")
    if len(residues) > MAX_LENGTH:
        raise SearchError(f"the query is {len(residues)} residues; the maximum is {MAX_LENGTH}")
    return Query(residues)


def check_email(email: str) -> str:
    if not EMAIL.fullmatch(email):
        raise SearchError("the search service requires a contact email (PHYLOBOT_SEARCH_EMAIL)")
    return email


def request(url: str, data: dict[str, str] | None = None) -> bytes:
    """GET, or POST `data` form-encoded; SearchError on any failure."""
    body = urllib.parse.urlencode(data).encode() if data is not None else None
    req = urllib.request.Request(url, data=body, headers={"User-Agent": "PhyloBot"})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as response:
            # TIMEOUT bounds each read, not the whole answer: a server that trickles
            # bytes could hold a request (and the search runner) for an hour. So the
            # answer is read in pieces against an overall deadline.
            deadline = time.monotonic() + TOTAL_TIMEOUT
            parts: list[bytes] = []
            while chunk := response.read(CHUNK):
                parts.append(chunk)
                if time.monotonic() > deadline:
                    raise SearchError(f"the search service took over {TOTAL_TIMEOUT} s to answer")
            return b"".join(parts)
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode(errors="replace")[:500]
        raise SearchError(f"the search service answered {exc.code}: {detail}") from exc
    except (urllib.error.URLError, TimeoutError) as exc:
        raise SearchError(f"could not reach the search service: {exc}") from exc


def wait(
    state: Callable[[], str],
    poll_seconds: float,
    first_seconds: float = 0,
    limit_seconds: float = 3600,
) -> None:
    """Block until `state()` reports FINISHED; SearchError if it FAILED or took too long.

    `state` raises SearchError for an unanswered or unexpected check; services time out
    now and then under load, so a long search tolerates a few of those in a row.
    """
    deadline = time.monotonic() + limit_seconds
    time.sleep(first_seconds)
    failures = 0
    while True:
        try:
            current = state()
            failures = 0
        except SearchError:
            failures += 1
            if failures > MAX_STATUS_FAILURES:
                raise
            current = PENDING
        if current == FINISHED:
            return
        if current == FAILED:
            raise SearchError("the search failed at the search service")
        if time.monotonic() > deadline:
            raise SearchError(f"the search did not finish in {limit_seconds:.0f} s")
        time.sleep(poll_seconds)
