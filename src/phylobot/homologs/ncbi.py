"""blastp at NCBI (https://blast.ncbi.nlm.nih.gov), through its URL API.

Fast (Ime2 against ClusteredNR: 39 s, against 16.5 min at EBI), but NCBI's guidelines
limit it to light use: at most one request every 10 s, one status check per search per
minute, and users with more than 100 searches a day are slowed or blocked; "not intended
for projects that involve a large number of BLAST searches". So it suits local mode, not
the hosted site, where every user's searches would come from one address.
https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html
"""

import re
import threading
import time
import urllib.parse

from phylobot.homologs import remote
from phylobot.homologs.remote import Query, SearchError

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
DATABASES = {"clustered_nr": "ClusteredNR", "swissprot": "Swiss-Prot"}
MAX_HITS = 1000  # as at EBI; NCBI may return more (Ime2: 2,000 for 1,000 asked)
EVALUE = "1e-3"
MIN_GAP_SECONDS = 10  # NCBI: no more than one request every 10 s
POLL_SECONDS = 60  # NCBI: no more than one status check per search per minute
TOOL = "phylobot"
RID = re.compile(r"^[A-Z0-9]{6,20}$")
STATES = {
    "WAITING": remote.PENDING,
    "READY": remote.FINISHED,
    "FAILED": remote.FAILED,
    "UNKNOWN": remote.FAILED,  # expired or never existed
}

_lock = threading.Lock()
_last_request = 0.0


def _request(params: dict[str, str], post: bool = False) -> bytes:
    """One request to NCBI, spaced at least MIN_GAP_SECONDS after the previous one from
    this process."""
    global _last_request
    with _lock:
        pause = _last_request + MIN_GAP_SECONDS - time.monotonic()
        if pause > 0:
            time.sleep(pause)
        try:
            if post:
                return remote.request(URL, params)
            return remote.request(URL + "?" + urllib.parse.urlencode(params))
        finally:
            _last_request = time.monotonic()


def check_rid(rid: str) -> str:
    if not RID.fullmatch(rid):
        raise SearchError(f"unexpected NCBI request ID {rid[:100]!r}")
    return rid


def submit(query: Query, email: str, database: str = "clustered_nr") -> tuple[str, int]:
    """Start a search; return NCBI's request ID and its estimated seconds to completion."""
    if database not in DATABASES:
        raise SearchError(f"unknown NCBI database {database!r}")
    params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": database,
        "EXPECT": EVALUE,
        "HITLIST_SIZE": str(MAX_HITS),
        "EMAIL": remote.check_email(email),
        "TOOL": TOOL,
        "QUERY": query.fasta,
    }
    page = _request(params, post=True).decode(errors="replace")
    rid = re.search(r"^\s*RID = (\S+)\s*$", page, re.MULTILINE)
    if not rid:
        raise SearchError("NCBI did not accept the search")
    rtoe = re.search(r"^\s*RTOE = (\d+)\s*$", page, re.MULTILINE)
    return check_rid(rid[1]), int(rtoe[1]) if rtoe else POLL_SECONDS


def state(rid: str) -> str:
    params = {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": check_rid(rid)}
    page = _request(params).decode(errors="replace")
    status = re.search(r"Status=([A-Z]+)", page)
    if not status or status[1] not in STATES:
        raise SearchError("unexpected answer from NCBI")
    return STATES[status[1]]


def fetch(rid: str) -> bytes:
    """The finished search as BLAST XML."""
    return _request({"CMD": "Get", "FORMAT_TYPE": "XML", "RID": check_rid(rid)})
