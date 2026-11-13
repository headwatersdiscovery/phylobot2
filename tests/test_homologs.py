"""Homolog search: query checks, the EBI and NCBI clients (mocked), hit parsing,
filters, their histograms and the filter form."""

from typing import Any

import pytest

from phylobot.forms import HomologFilterForm
from phylobot.homologs import ebi, examples, filters, ncbi, plots, remote, sources
from phylobot.homologs.hits import (
    Hit,
    ParseError,
    SearchResult,
    coverage,
    parse_description,
    parse_ebi_xml,
    parse_ncbi_xml,
)
from phylobot.homologs.remote import SearchError

QUERY = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"


def hit(**fields: Any) -> Hit:
    defaults: dict[str, Any] = dict(
        id="UniRef90_X",
        name="p",
        taxon="t",
        taxid=1,
        members=1,
        length=100,
        evalue=1e-50,
        bits=200.0,
        identity=60.0,
        coverage=95.0,
    )
    return Hit(**{**defaults, **fields})


def ebi_ime2() -> SearchResult:
    """The pinned Ime2 search at EBI (see phylobot.homologs.examples)."""
    return parse_ebi_xml(examples.raw("ime2-ebi"))


def ncbi_ime2() -> SearchResult:
    """The pinned Ime2 search at NCBI; it returned 2,000 hits for 1,000 asked."""
    return parse_ncbi_xml(examples.raw("ime2-ncbi"), 1000)


# --- query -------------------------------------------------------------------------


def test_query_from_fasta_drops_the_name() -> None:
    query = remote.clean_query(f">secret project name\n{QUERY[:30]}\n{QUERY[30:].lower()}*\n")
    assert query.residues == QUERY
    assert query.fasta == f">query\n{QUERY}\n"


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("", "empty"),
        (">a\nMKT\n>b\nMKT", "one query"),
        (QUERY + "1", "not amino acids"),
        (QUERY + "; rm -rf /", "not amino acids"),
        ("MKT", "minimum"),
        ("M" * (remote.MAX_LENGTH + 1), "maximum"),
    ],
)
def test_bad_queries_rejected(text: str, message: str) -> None:
    with pytest.raises(SearchError, match=message):
        remote.clean_query(text)


# --- waiting -----------------------------------------------------------------------


def test_wait_survives_slow_status_answers() -> None:
    answers: list[str | Exception] = [
        SearchError("timeout"),
        SearchError("unexpected answer"),
        remote.PENDING,
        remote.FINISHED,
    ]

    def state() -> str:
        answer = answers.pop(0)
        if isinstance(answer, Exception):
            raise answer
        return answer

    remote.wait(state, poll_seconds=0)
    assert answers == []


def test_wait_reports_failed_search() -> None:
    with pytest.raises(SearchError, match="failed"):
        remote.wait(lambda: remote.FAILED, poll_seconds=0)


def test_wait_gives_up_when_the_service_keeps_failing() -> None:
    def state() -> str:
        raise SearchError("unexpected answer")

    with pytest.raises(SearchError, match="unexpected answer"):
        remote.wait(state, poll_seconds=0)


# --- EBI client --------------------------------------------------------------------

Sent = list[tuple[str, dict[str, str] | None]]


@pytest.fixture
def ebi_requests(monkeypatch: pytest.MonkeyPatch) -> Sent:
    sent: Sent = []
    answers = {"run": b"ncbiblast-R20260929-095212-0255-17798067-p1m\n"}

    def fake(url: str, data: dict[str, str] | None = None) -> bytes:
        sent.append((url, data))
        return answers.get(url.rsplit("/", 1)[-1], b"FINISHED")

    monkeypatch.setattr(remote, "request", fake)
    return sent


def test_ebi_submit_sends_only_residues(ebi_requests: Sent) -> None:
    job_id = ebi.submit(remote.clean_query(f">my name\n{QUERY}"), "someone@example.org")
    assert job_id == "ncbiblast-R20260929-095212-0255-17798067-p1m"
    url, data = ebi_requests[0]
    assert url.endswith("/ncbiblast/run")
    assert data is not None
    assert data["sequence"] == f">query\n{QUERY}\n"
    assert data["database"] == "uniref90" and data["program"] == "blastp"
    assert "my name" not in str(data)


def test_ebi_submit_refuses_bad_database_and_email(ebi_requests: Sent) -> None:
    query = remote.clean_query(QUERY)
    with pytest.raises(SearchError, match="database"):
        ebi.submit(query, "someone@example.org", database="uniprotkb")
    with pytest.raises(SearchError, match="email"):
        ebi.submit(query, "")
    assert ebi_requests == []


def test_ebi_states(monkeypatch: pytest.MonkeyPatch) -> None:
    answers = [b"RUNNING", b"<!doctype html><title>Error: 500</title>"]
    monkeypatch.setattr(remote, "request", lambda url, data=None: answers.pop(0))
    assert ebi.state("ncbiblast-R1") == remote.PENDING
    with pytest.raises(SearchError, match="unexpected answer"):
        ebi.state("ncbiblast-R1")


@pytest.mark.parametrize("job_id", ["../../run", "ncbiblast-x/../y", "phmmer-R1", ""])
def test_ebi_job_ids_are_checked_before_use(job_id: str) -> None:
    with pytest.raises(SearchError, match="job ID"):
        ebi.state(job_id)


# --- NCBI client -------------------------------------------------------------------


@pytest.fixture
def ncbi_requests(monkeypatch: pytest.MonkeyPatch) -> Sent:
    sent: Sent = []

    def fake(url: str, data: dict[str, str] | None = None) -> bytes:
        sent.append((url, data))
        if data is not None:
            return b"<!--QBlastInfoBegin\n    RID = BS9DVWPC016\n    RTOE = 30\nQBlastInfoEnd\n-->"
        return b"<!--QBlastInfoBegin\n\tStatus=READY\nQBlastInfoEnd\n-->"

    monkeypatch.setattr(remote, "request", fake)
    monkeypatch.setattr(ncbi, "MIN_GAP_SECONDS", 0)
    return sent


def test_ncbi_submit_sends_only_residues(ncbi_requests: Sent) -> None:
    rid, estimate = ncbi.submit(remote.clean_query(f">my name\n{QUERY}"), "someone@example.org")
    assert (rid, estimate) == ("BS9DVWPC016", 30)
    url, data = ncbi_requests[0]
    assert data is not None
    assert data["QUERY"] == f">query\n{QUERY}\n"
    assert (data["CMD"], data["DATABASE"], data["TOOL"]) == ("Put", "clustered_nr", "phylobot")
    assert "my name" not in str(data)
    assert ncbi.state(rid) == remote.FINISHED
    assert ncbi_requests[1][0].endswith("CMD=Get&FORMAT_OBJECT=SearchInfo&RID=BS9DVWPC016")


def test_ncbi_refuses_bad_input(ncbi_requests: Sent) -> None:
    with pytest.raises(SearchError, match="database"):
        ncbi.submit(remote.clean_query(QUERY), "someone@example.org", database="nr; x")
    with pytest.raises(SearchError, match="request ID"):
        ncbi.state("BS9&CMD=Delete")
    assert ncbi_requests == []


def test_ncbi_spaces_its_requests(monkeypatch: pytest.MonkeyPatch) -> None:
    """NCBI asks for at least 10 s between requests."""
    clock = [1000.0]
    slept: list[float] = []
    monkeypatch.setattr("phylobot.homologs.ncbi.time.monotonic", lambda: clock[0])
    monkeypatch.setattr("phylobot.homologs.ncbi.time.sleep", lambda s: slept.append(s))
    monkeypatch.setattr(remote, "request", lambda url, data=None: b"Status=WAITING")
    monkeypatch.setattr(ncbi, "_last_request", 0.0)
    ncbi.state("BS9DVWPC016")
    clock[0] += 3
    ncbi.state("BS9DVWPC016")
    assert slept == [ncbi.MIN_GAP_SECONDS - 3]


def test_sources_registry() -> None:
    assert set(sources.SOURCES) == {"ebi", "ncbi"}
    assert sources.SOURCES["ncbi"].local_only and not sources.SOURCES["ebi"].local_only
    for source in sources.SOURCES.values():
        assert source.default in source.databases
    with pytest.raises(SearchError, match="unknown search service"):
        sources.search("google", remote.clean_query(QUERY), "someone@example.org")


# --- hits --------------------------------------------------------------------------


def test_uniref_description() -> None:
    assert parse_description(
        "Meiosis-induced protein kinase IME2 n=5 Tax=Saccharomyces TaxID=4930 RepID=IME2_YEAST"
    ) == ("Meiosis-induced protein kinase IME2", "Saccharomyces", 4930, 5)
    assert parse_description("Something else") == ("Something else", "", None, None)


def test_coverage_merges_overlapping_hsps() -> None:
    assert coverage(100, [(1, 50), (41, 60)]) == 60.0
    assert coverage(100, [(60, 41)]) == 20.0
    assert coverage(10, [(5, 30)]) == 60.0  # clipped to the query


def test_parse_ebi_ime2() -> None:
    search = ebi_ime2()
    assert (search.program, search.database) == ("BLASTP 2.16.0+", "uniref90")
    assert search.database_created.startswith("2026-08-27")
    assert search.database_sequences == 121_494_192
    assert search.query_length == 645
    assert len(search.hits) == 1000 and search.truncated
    assert search.hits[0] == Hit(
        id="UniRef90_P32581",
        name="Meiosis induction protein kinase IME2/SME1",
        taxon="Saccharomyces",
        taxid=4930,
        members=10,
        length=645,
        evalue=0.0,
        bits=1310.0,
        identity=100.0,
        coverage=100.0,
    )


def test_parse_ncbi_ime2() -> None:
    search = ncbi_ime2()
    assert (search.program, search.database) == ("BLASTP 2.17.0+", "clustered_nr")
    assert search.database_created == ""  # NCBI's XML does not say
    assert search.database_sequences == 545_824_732
    assert len(search.hits) == 2000 and search.truncated
    top = search.hits[0]
    assert (top.id, top.name, top.taxon) == (
        "NP_012429",
        "protein kinase IME2",
        "Saccharomyces cerevisiae S288C",
    )
    assert (top.identity, top.coverage, top.taxid) == (100.0, 100.0, None)


@pytest.mark.parametrize(
    ("content", "message"),
    [
        (b"not xml", "malformed"),
        (b"<EBIApplicationResult xmlns='http://www.ebi.ac.uk/schema'/>", "lacks"),
    ],
)
def test_bad_ebi_xml_rejected(content: bytes, message: str) -> None:
    with pytest.raises(ParseError, match=message):
        parse_ebi_xml(content)


def test_bad_ncbi_xml_rejected() -> None:
    with pytest.raises(ParseError, match="lacks BlastOutput"):
        parse_ncbi_xml(b"<html/>", 1000)
    with pytest.raises(ParseError, match="lacks"):
        parse_ncbi_xml(
            b"<BlastOutput><BlastOutput_query-len>9</BlastOutput_query-len></BlastOutput>", 1
        )


# --- filters -----------------------------------------------------------------------


def test_filters_count_each_stage() -> None:
    hits = [
        hit(id="keep"),
        hit(id="weak", evalue=1e-3),
        hit(id="fragment", coverage=40.0),
        hit(id="too close", identity=99.0),
        hit(id="too far", identity=20.0),
    ]
    kept, stages = filters.apply(hits, filters.Settings())
    assert [h.id for h in kept] == ["keep"]
    assert [(s.name, s.remaining) for s in stages] == [
        ("search", 5),
        ("evalue", 4),
        ("coverage", 3),
        ("identity", 1),
    ]


def test_the_query_self_hit_is_kept_despite_the_identity_filter() -> None:
    # The query appears as a full-length, fully identical hit; the identity filter would
    # drop it (> 90%), but it is the sequence the reconstruction starts from.
    hits = [hit(id="keep", identity=60.0), hit(id="query", identity=100.0, coverage=100.0)]
    assert filters.query_ids(hits) == frozenset({"query"})
    always = filters.query_ids(hits)
    s = filters.select(hits, filters.Settings(), [], [], always=always)
    assert [h.id for h in s.selected] == ["keep", "query"]  # the query is in, by default
    assert "query" in s.default and s.always == frozenset({"query"})
    assert [h.id for h in s.included] == []  # not reported as a hand-pick
    assert s.stages[-1].description == "with the query kept"
    # The user can still untick it (an exclude, since it sits in the default now).
    dropped = filters.select(hits, filters.Settings(), [], ["query"], always=always)
    assert [h.id for h in dropped.selected] == ["keep"]
    # "Select none" starts from nothing, query included.
    none = filters.select(hits, filters.Settings(), [], [], base="none", always=always)
    assert none.selected == []


def test_presets_and_settings_checks() -> None:
    recent = filters.Settings.preset("recent", min_coverage=80)
    assert (recent.min_identity, recent.min_coverage) == (50.0, 80)
    with pytest.raises(ValueError, match="above"):
        filters.Settings(min_identity=80, max_identity=70)
    with pytest.raises(ValueError, match="between"):
        filters.Settings(min_coverage=120)


@pytest.mark.parametrize(
    ("search", "min_coverage", "counts"),
    [
        # Most hits align only the kinase domain (~55% of Ime2): 70% removes them, and
        # with them true Ime2 orthologs whose C-terminal half has diverged (Yarrowia,
        # Candida); 50% lets nearly every kinase through. Hence a user setting, shown
        # with the coverage histogram.
        (ebi_ime2, 70, [1000, 1000, 74, 62]),
        (ebi_ime2, 50, [1000, 1000, 979, 967]),
        (ncbi_ime2, 70, [2000, 2000, 34, 30]),
    ],
)
def test_ime2_filter_counts(search: Any, min_coverage: float, counts: list[int]) -> None:
    _, stages = filters.apply(search().hits, filters.Settings(min_coverage=min_coverage))
    assert [s.remaining for s in stages] == counts


# --- histograms --------------------------------------------------------------------


def test_histograms_split_each_bin_at_the_thresholds() -> None:
    hits = [
        hit(coverage=52, identity=40),
        hit(coverage=68, identity=40),
        hit(coverage=72, identity=40),
        hit(coverage=100, identity=95),
        hit(coverage=100, identity=20, evalue=1.0),  # removed before coverage
    ]
    cov, ident = plots.filter_histograms(hits, filters.Settings(min_coverage=70))
    assert (cov.total, cov.kept, ident.total, ident.kept) == (4, 2, 2, 1)
    by_label = {bar.label: (bar.kept_count, bar.removed_count) for bar in cov.bars}
    assert by_label["50-55%"] == (0, 1)
    assert by_label["65-70%"] == (0, 1)
    assert by_label["70-75%"] == (1, 0)
    assert by_label["95-100%"] == (1, 0)  # 100% falls in the last bin
    assert [m.label for m in cov.markers] == ["70%"]
    assert [m.label for m in ident.markers] == ["30%", "90%"]


def test_histogram_bar_geometry() -> None:
    """Segments stack from the baseline, removed above kept, with a gap between."""
    three = [hit(), hit(), hit()]
    h = plots.histogram("x", "X", "x", three, lambda h: 50.0, 40, 60)
    bar = next(b for b in h.bars if b.label == "50-55%")
    assert bar.kept is not None and bar.removed is None
    assert bar.kept.y + bar.kept.height == pytest.approx(h.frame.bottom)
    top = float(h.y_ticks[-1][1].replace(",", ""))
    assert bar.kept.height == pytest.approx(3 / top * h.frame.plot_height, abs=0.1)
    two = [hit(coverage=41), hit(coverage=44)]
    mixed = plots.histogram("x", "X", "x", two, lambda h: h.coverage, 42, 100)
    both = next(b for b in mixed.bars if b.label == "40-45%")
    assert both.kept is not None and both.removed is not None
    assert both.removed.y + both.removed.height == pytest.approx(both.kept.y - plots.GAP)


# --- filter form -------------------------------------------------------------------


def test_filter_form_suggests_70_and_accepts_lower() -> None:
    assert HomologFilterForm().fields["min_coverage"].initial == 70
    form = HomologFilterForm(
        {"max_evalue": "1e-05", "min_coverage": "55", "min_identity": "30", "max_identity": "90"}
    )
    assert form.is_valid(), form.errors
    assert form.cleaned_data["settings"] == filters.Settings(min_coverage=55)


@pytest.mark.parametrize(
    "data",
    [
        {"max_evalue": "1e-05", "min_coverage": "150", "min_identity": "30", "max_identity": "90"},
        {"max_evalue": "7", "min_coverage": "70", "min_identity": "30", "max_identity": "90"},
        {"max_evalue": "1e-05", "min_coverage": "70", "min_identity": "80", "max_identity": "40"},
        {"max_evalue": "1e-05", "min_coverage": "nan", "min_identity": "30", "max_identity": "90"},
    ],
)
def test_filter_form_rejects_bad_values(data: dict[str, str]) -> None:
    assert not HomologFilterForm(data).is_valid()


def test_a_trickling_answer_is_cut_off(monkeypatch: pytest.MonkeyPatch) -> None:
    """A server sending a byte now and then must not hold a request indefinitely."""
    clock = [0.0]

    class Trickle:
        def __enter__(self) -> Trickle:
            return self

        def __exit__(self, *args: object) -> None:
            pass

        def read(self, size: int) -> bytes:
            clock[0] += 50  # each read arrives within TIMEOUT, but they never stop
            return b"x"

    monkeypatch.setattr(
        "phylobot.homologs.remote.urllib.request.urlopen", lambda *a, **k: Trickle()
    )
    monkeypatch.setattr("phylobot.homologs.remote.time.monotonic", lambda: clock[0])
    with pytest.raises(SearchError, match="took over"):
        remote.request("https://example.org/")
