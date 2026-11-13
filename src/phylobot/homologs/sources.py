"""The remote search services a user can choose between.

Each service is driven in steps (submit, check, fetch, parse) so one runner can follow
many searches at once without a thread waiting on each (phylobot.homologs.runner).
`search` runs all the steps in a row, for scripts and one-off use.
"""

from collections.abc import Callable
from dataclasses import dataclass

from phylobot.homologs import ebi, ncbi, remote, sequences
from phylobot.homologs.hits import SearchResult, parse_ebi_xml, parse_ncbi_xml
from phylobot.homologs.remote import Query, SearchError


@dataclass(frozen=True)
class Source:
    label: str
    databases: dict[str, str]  # service's name -> label
    default: str
    note: str
    local_only: bool  # the service's usage limits rule out the hosted site
    # Searches in flight at the service at once, across the whole site: EBI asks for
    # batches of at most 30; NCBI's pacing (one request per 10 s) allows one.
    max_in_flight: int
    poll_seconds: int  # between status checks of one search
    submit: Callable[[Query, str, str], tuple[str, int]]  # -> (service's ID, first wait s)
    state: Callable[[str], str]  # -> remote.PENDING / FINISHED / FAILED
    fetch: Callable[[str], bytes]  # the finished search's XML
    parse: Callable[[bytes], SearchResult]
    # The hits' full-length sequences, {hit ID: sequence}, given their IDs and the email.
    fetch_sequences: Callable[[list[str], str], dict[str, str]]


def _ebi_submit(query: Query, email: str, database: str) -> tuple[str, int]:
    return ebi.submit(query, email, database), ebi.POLL_SECONDS


def _ncbi_submit(query: Query, email: str, database: str) -> tuple[str, int]:
    rid, estimate = ncbi.submit(query, email, database)
    return rid, min(max(estimate, ncbi.MIN_GAP_SECONDS), ncbi.POLL_SECONDS)


SOURCES = {
    "ebi": Source(
        label="EBI (UniProt databases)",
        databases=ebi.DATABASES,
        default="uniref90",
        note="UniProt taxonomy and annotation; slow when busy (Ime2: 16 min).",
        local_only=False,
        max_in_flight=25,
        poll_seconds=ebi.POLL_SECONDS,
        submit=_ebi_submit,
        state=ebi.state,
        fetch=ebi.fetch,
        parse=parse_ebi_xml,
        fetch_sequences=sequences.fetch_uniprot,
    ),
    "ncbi": Source(
        label="NCBI",
        databases=ncbi.DATABASES,
        default="clustered_nr",
        note="Fast (Ime2: under a minute); NCBI limits heavy use, so local mode only.",
        local_only=True,
        max_in_flight=1,
        poll_seconds=ncbi.POLL_SECONDS,
        submit=_ncbi_submit,
        state=ncbi.state,
        fetch=ncbi.fetch,
        parse=lambda content: parse_ncbi_xml(content, ncbi.MAX_HITS),
        fetch_sequences=sequences.fetch_ncbi,
    ),
}
DEFAULT_SOURCE = "ebi"


@dataclass(frozen=True)
class About:
    """What the search page says beside a database choice."""

    speed: str  # "fast" or "slow"
    timing: str  # measured where we have a measurement; "expected" otherwise
    note: str


# Keyed "service:database", as the search form's choices. Keep one entry per database
# (tested). Timings: Ime2 (645 aa), 2026-09-29.
ABOUT = {
    "ebi:uniref90": About(
        "slow",
        "Ime2: 16 min",
        "UniProt clustered at 90% identity (121 million clusters). Best taxonomy and "
        "annotation for reviewing hits. Recommended when you can wait.",
    ),
    "ebi:uniprotkb_swissprot": About(
        "fast",
        "expected under a few minutes (small database; not yet measured)",
        "Reviewed UniProt entries only (about 570,000). Well annotated, but few "
        "sequences per family.",
    ),
    "ebi:uniref50": About(
        "slow",
        "expected several minutes: quicker than UniRef90, still slow (not yet measured)",
        "UniProt clustered at 50% identity. For very large families; each hit stands for "
        "more varied sequences, so its species name is less specific.",
    ),
    "ncbi:clustered_nr": About(
        "fast",
        "Ime2: 39 s and 61 s",
        "NCBI's non-redundant proteins, clustered (the largest choice). No UniProt "
        "annotation. Offered only when PhyloBot runs on your own machine.",
    ),
    "ncbi:swissprot": About(
        "slow",
        "unreliable: 3 of 3 searches still waiting after 25 min to 2 h (2026-09-30), while "
        "ClusteredNR answered in under a minute; use EBI Swiss-Prot instead",
        "NCBI's copy of reviewed UniProt entries. Curated only, few sequences per family.",
    ),
}


def get(source: str) -> Source:
    if source not in SOURCES:
        raise SearchError(f"unknown search service {source!r}")
    return SOURCES[source]


def search(source: str, query: Query, email: str, database: str = "") -> tuple[SearchResult, bytes]:
    """Run a search to the end at the chosen service (its default database if none):
    the parsed result, and the raw XML to keep. Blocks for as long as the search takes."""
    chosen = get(source)
    remote_id, first = chosen.submit(query, email, database or chosen.default)
    remote.wait(lambda: chosen.state(remote_id), chosen.poll_seconds, first_seconds=first)
    content = chosen.fetch(remote_id)
    return chosen.parse(content), content
