"""End to end: the real stack, driven over HTTP the way a browser drives it.

Needs `make up` (use `make e2e`, which starts the stack first). Submits the Ime2 example
through the web form, follows the job to completion, checks every part of the job page
and its downloads, tags it, then deletes it and checks its files are gone but its record
is kept. Standard library only (urllib with a cookie jar for the session and CSRF).
"""

import html
import http.cookiejar
import json
import os
import re
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
BASE = os.environ.get("E2E_BASE_URL", "http://127.0.0.1:8000")
IME2 = (ROOT / "src" / "phylobot" / "examples" / "ime2.fasta").read_text()
JOB_TIMEOUT = 15 * 60


class Browser:
    """Just enough of a browser: cookies, CSRF tokens, form posts, redirects."""

    def __init__(self) -> None:
        self.cookies = http.cookiejar.CookieJar()
        self.opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(self.cookies))
        self.csrf = ""

    def get(self, path: str) -> tuple[int, str, str]:
        """(status, final URL path, body)."""
        try:
            with self.opener.open(BASE + path, timeout=60) as response:
                body = response.read().decode()
                self.remember_csrf(body)
                return response.status, urllib.parse.urlparse(response.url).path, body
        except urllib.error.HTTPError as error:
            return error.code, path, error.read().decode(errors="replace")

    def post(
        self, path: str, fields: list[tuple[str, str]], referer: str = "/"
    ) -> tuple[int, str, str]:
        data = urllib.parse.urlencode([("csrfmiddlewaretoken", self.csrf), *fields]).encode()
        headers = {"Referer": BASE + referer}
        request = urllib.request.Request(BASE + path, data=data, headers=headers)
        try:
            with self.opener.open(request, timeout=120) as response:
                body = response.read().decode()
                self.remember_csrf(body)
                return response.status, urllib.parse.urlparse(response.url).path, body
        except urllib.error.HTTPError as error:
            return error.code, path, error.read().decode(errors="replace")

    def remember_csrf(self, body: str) -> None:
        if match := re.search(r'name="csrfmiddlewaretoken" value="([^"]+)"', body):
            self.csrf = match.group(1)


def storage_keys(prefix: str) -> list[str]:
    """Object storage contents, checked from inside the web container."""
    code = (
        "import django, os; os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phylobot.settings');"
        "django.setup(); from phylobot import storage;"
        f"print('\\n'.join(storage.list_keys({prefix!r})))"
    )
    out = subprocess.run(
        ["docker", "compose", "exec", "-T", "web", "python", "-c", code],
        cwd=ROOT, capture_output=True, text=True, check=True,
    )  # fmt: skip
    return [line for line in out.stdout.splitlines() if line.strip()]


def text(fragment: str) -> str:
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", fragment))).strip()


@pytest.fixture(scope="module")
def browser() -> Browser:
    b = Browser()
    try:
        status, _, _ = b.get("/healthz")
    except OSError:
        status = 0
    if status != 200:
        pytest.fail(f"The stack is not up at {BASE} (run `make e2e`, which starts it).")
    return b


@pytest.mark.e2e
def test_submit_follow_inspect_tag_and_delete(browser: Browser) -> None:
    # 1. The submit page, and the example loader.
    status, _, page = browser.get("/jobs/new/")
    assert status == 200 and "<h1>New job</h1>" in page
    assert page.count('type="checkbox" name="models"') == 252
    for aligner in ("mafft", "muscle", "clustalo", "prank"):
        assert f'value="{aligner}"' in page
    _, _, example = browser.get("/examples/ime2/")
    assert "Saccharomyces.cerevisiae.IME2" in example

    # 2. Invalid input is rejected on the form, before any job exists.
    _, _, before = browser.get("/jobs/")
    status, path, page = browser.post("/jobs/new/", [
        ("sequences_text", ">a\nMKT\n>b\nMKV\n"), ("aligners", "mafft"),
        ("tree_builders", "iqtree"), ("modelfinder", "on"), ("bootstrap", "1000"), ("seed", "1"),
        ("run_trees", "on"), ("run_ancestors", "on"),
    ])  # fmt: skip
    assert status == 200 and path == "/jobs/new/" and "need at least 4 sequences" in page
    _, _, after = browser.get("/jobs/")
    assert after.count("data-delete=") == before.count("data-delete=")

    # 3. A real job: 2 aligners x (ModelFinder + LG+G4) = 4 combinations.
    name = f"e2e {time.strftime('%Y-%m-%d %H:%M:%S')}"
    status, path, page = browser.post("/jobs/new/", [
        ("name", name), ("sequences_text", IME2),
        ("aligners", "mafft"), ("aligners", "muscle"), ("tree_builders", "iqtree"),
        ("modelfinder", "on"), ("fixed_models", "on"), ("models", "LG+G4"),
        ("gap_methods", "parsimony"), ("outgroup", "Yarrowia.lipolytica.IME2"),
        ("sh_alrt", "on"), ("alrt", "1000"), ("ufboot", "on"), ("bootstrap", "1000"),
        ("seed", "1"), ("run_trees", "on"), ("run_ancestors", "on"),
    ])  # fmt: skip
    match = re.fullmatch(r"/jobs/([0-9a-f-]{36})/", path)
    assert status == 200 and match, f"expected the job page, got {path}"
    job_id = match.group(1)
    assert "4 combinations" in page

    # 4. Follow it to completion through the live status fragment.
    seen_progress, deadline, final = False, time.monotonic() + JOB_TIMEOUT, ""
    while time.monotonic() < deadline:
        _, _, fragment = browser.get(f"/jobs/{job_id}/status/")
        state = re.search(r'class="status (\w+)"', fragment)
        assert state, fragment[:500]
        seen_progress |= bool(re.search(r"\d+ of \d+ steps", fragment))
        if state.group(1) in {"succeeded", "failed", "cancelled"}:
            final = state.group(1)
            break
        assert 'hx-trigger="every 5s"' in fragment  # still polling while unfinished
        time.sleep(5)
    assert final == "succeeded", text(fragment)[:1000]
    assert seen_progress and "hx-trigger" not in fragment  # polling stops when done

    # 5. The finished job page.
    status, _, page = browser.get(f"/jobs/{job_id}/")
    assert status == 200
    assert "Input sequences" in page and "7 sequences, 369–800 residues" in text(page)
    results = page[page.index(">Results</button></h2>") : page.index("<h2>All files</h2>")]
    for label in ("MAFFT", "MUSCLE"):
        assert f">Alignment: {label}</button></h3>" in results
    assert 'id="alignments-overview">Alignments</h3>' in results  # overview first
    assert results.count('<table class="results">') == 2 + 2  # 2 alignments + 2 candidates
    assert results.count("/alignment.fasta") == 2 * 2  # overview + each heading
    # With only LG+G4 ticked ModelFinder picks LG+G4, so each alignment's two combinations
    # share one results row (results.merge_modelfinder_duplicates).
    assert results.count("(also ticked)") == 2
    assert results.count("ancestors.zip?combination=") == 2 and results.count("/tree.nwk") == 2
    assert "ModelFinder\u00a0→" in results and "background-color: rgb(" in results

    # Downloads hold real content.
    _, _, fasta = browser.get(f"/jobs/{job_id}/files/mafft/iqtree/modelfinder/ancestors.fasta")
    records = re.findall(r"^>(\S+)\n([A-Z]+)$", fasta, re.MULTILINE)
    assert len(records) == 5 and all(len(seq) > 300 for _, seq in records)
    assert results.count("shared presence, ") == 2  # each row's ZIP holds a masked version
    assert results.count(">rooted</a>") == 2
    masked_path = "mafft/iqtree/modelfinder/ancestors.gaps-parsimony.fasta"
    _, _, masked = browser.get(f"/jobs/{job_id}/files/{masked_path}")
    unmasked = dict(re.findall(r"^>(\S+)\n([A-Z-]+)$", fasta, re.MULTILINE))
    for node, row in re.findall(r"^>(\S+)\n([A-Z-]+)$", masked, re.MULTILINE):
        assert len(row) == len(unmasked[node]) and "-" in row  # same columns, some masked
        assert all(m in ("-", u) for m, u in zip(row, unmasked[node], strict=True))
    _, _, supports = browser.get(f"/jobs/{job_id}/files/mafft/iqtree/modelfinder/supports.tsv")
    header, *nodes = [line.split("\t") for line in supports.splitlines()]
    assert header == ["node", "ufboot", "sh_alrt"]  # both measures ticked on the form
    supported = [(uf, sh) for node, uf, sh in nodes if node != "Node1"]  # top node: none
    assert supported and all(uf and sh for uf, sh in supported)
    _, _, comparison = browser.get(f"/jobs/{job_id}/files/model_comparison.tsv")
    rows = [line.split("\t") for line in comparison.splitlines()[1:]]
    assert sorted({row[0] for row in rows}) == ["mafft", "muscle"] and len(rows) == 4
    _, _, submitted = browser.get(f"/jobs/{job_id}/inputs/sequences.fasta")
    assert submitted == IME2
    summary = json.loads(browser.get(f"/jobs/{job_id}/files/summary.json")[2])
    assert len(summary["combinations"]) == 4

    # 6. My Jobs: listed, taggable, filterable.
    _, _, listing = browser.get("/jobs/")
    assert name in listing
    status, _, listing = browser.post(
        "/jobs/tags/", [("action", "add"), ("tags", "#e2e"), ("job", job_id), ("next", "/jobs/")],
        referer="/jobs/",
    )  # fmt: skip
    assert "Added #e2e to 1 job." in listing
    _, _, filtered = browser.get("/jobs/?tag=e2e")
    assert name in filtered

    # 7. Delete: files gone, record kept.
    assert storage_keys(f"jobs/{job_id}/")
    status, path, listing = browser.post("/jobs/delete/", [("job", job_id)], referer="/jobs/")
    assert path == "/jobs/" and "Deleted 1 job." in listing
    assert browser.get(f"/jobs/{job_id}/")[0] == 404
    deleted = listing[listing.index("Deleted jobs") :]
    assert name in deleted and "#e2e" in deleted
    assert storage_keys(f"jobs/{job_id}/") == []
