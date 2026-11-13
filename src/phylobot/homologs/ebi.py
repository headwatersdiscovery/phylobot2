"""blastp at the EBI Job Dispatcher (https://www.ebi.ac.uk/jdispatcher/), REST API.

Owner decision 2026-09-29: UniRef90 (Swiss-Prot and UniRef50 as alternatives). EBI's
phmmer has no UniRef and returns at most 100 hits; blastp returns up to 1,000. Searches
of UniRef90 can take a quarter of an hour (Ime2: 16.5 min).
"""

import re

from phylobot.homologs import remote
from phylobot.homologs.remote import Query, SearchError

BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
DATABASES = {"uniref90": "UniRef90", "uniprotkb_swissprot": "Swiss-Prot", "uniref50": "UniRef50"}
MAX_HITS = 1000  # the largest value EBI's blastp accepts for scores and alignments
EVALUE = "1e-3"  # report threshold (EBI accepts only fixed values); users filter tighter
POLL_SECONDS = 15
JOB_ID = re.compile(r"^ncbiblast-[A-Za-z0-9-]{1,80}$")
# EBI's job states. Any other answer (EBI sometimes serves an HTML error page) is a
# failed check, retried like a timeout.
STATES = {
    "QUEUED": remote.PENDING,
    "RUNNING": remote.PENDING,
    "FINISHED": remote.FINISHED,
    "FAILURE": remote.FAILED,
    "ERROR": remote.FAILED,
    "NOT_FOUND": remote.FAILED,
}


def check_job_id(job_id: str) -> str:
    if not JOB_ID.fullmatch(job_id):
        raise SearchError(f"unexpected EBI job ID {job_id[:100]!r}")
    return job_id


def submit(query: Query, email: str, database: str = "uniref90") -> str:
    """Start a search; return EBI's job ID."""
    if database not in DATABASES:
        raise SearchError(f"unknown EBI database {database!r}")
    params = {
        "email": remote.check_email(email),
        "program": "blastp",
        "stype": "protein",
        "database": database,
        "exp": EVALUE,
        "scores": str(MAX_HITS),
        "alignments": str(MAX_HITS),
        "title": "phylobot",
        "sequence": query.fasta,
    }
    return check_job_id(remote.request(f"{BASE_URL}/run", params).decode().strip())


def state(job_id: str) -> str:
    answer = remote.request(f"{BASE_URL}/status/{check_job_id(job_id)}").decode().strip()
    if answer not in STATES:
        raise SearchError(f"unexpected answer from EBI: {answer[:100]!r}")
    return STATES[answer]


def fetch(job_id: str) -> bytes:
    """The finished search as EBI's XML."""
    return remote.request(f"{BASE_URL}/result/{check_job_id(job_id)}/xml")
