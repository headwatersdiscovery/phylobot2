"""Find a Sequence by Name (phylobot.homologs.lookup): the UniProt query is built from
cleaned user text, the JSON is parsed reviewed-first, and the page offers each hit to the
homolog search below. The network is always mocked."""

import json
from typing import Any

import pytest
from django.test import Client

from phylobot.homologs import lookup, remote
from phylobot.homologs.remote import SearchError

pytestmark = pytest.mark.django_db

# Two UniProt entries: one unreviewed (TrEMBL), one reviewed (Swiss-Prot), plus one with no
# sequence (dropped). Trimmed to the fields the parser reads.
UNIPROT_JSON = {
    "results": [
        {
            "entryType": "UniProtKB unreviewed (TrEMBL)",
            "primaryAccession": "A0A1B2C3",
            "uniProtkbId": "A0A1B2C3_YEAST",
            "organism": {"scientificName": "Saccharomyces cerevisiae", "taxonId": 4932},
            "proteinDescription": {"submissionNames": [{"fullName": {"value": "Protein kinase"}}]},
            "genes": [{"geneName": {"value": "IME2"}}],
            "sequence": {"value": "MKLVAKQ", "length": 7},
        },
        {
            "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "primaryAccession": "P32581",
            "uniProtkbId": "IME2_YEAST",
            "organism": {"scientificName": "Saccharomyces cerevisiae", "taxonId": 559292},
            "proteinDescription": {
                "recommendedName": {"fullName": {"value": "Meiosis induction protein kinase IME2"}}
            },
            "genes": [{"geneName": {"value": "IME2"}}],
            "sequence": {"value": "MKLVAKQRST", "length": 10},
        },
        {  # no sequence: dropped
            "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "primaryAccession": "NOSEQ1",
            "organism": {"scientificName": "Saccharomyces cerevisiae"},
            "sequence": {},
        },
    ]
}


@pytest.fixture(autouse=True)
def _clear_suggestion_cache() -> Any:
    """Suggestions are lru_cached per process; clear between tests so a mocked answer in one
    test is not served to the next."""
    lookup.suggest_species.cache_clear()
    lookup.suggest_genes.cache_clear()
    yield
    lookup.suggest_species.cache_clear()
    lookup.suggest_genes.cache_clear()


def answer(monkeypatch: pytest.MonkeyPatch, payload: Any) -> list[str]:
    """Mock remote.request to return `payload` (dict -> JSON, else bytes); record URLs."""
    urls: list[str] = []
    body = json.dumps(payload).encode() if not isinstance(payload, bytes) else payload

    def fake(url: str, data: Any = None) -> bytes:
        urls.append(url)
        return body

    monkeypatch.setattr(remote, "request", fake)
    return urls


@pytest.fixture
def uniprot(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Mock the network: record the URLs requested, answer with the sample JSON."""
    urls: list[str] = []

    def fake(url: str, data: dict[str, str] | None = None) -> bytes:
        urls.append(url)
        return json.dumps(UNIPROT_JSON).encode()

    monkeypatch.setattr(remote, "request", fake)
    return urls


# ── the lookup logic ────────────────────────────────────────────────────────────────────


def test_terms_are_cleaned_to_a_safe_quotable_form() -> None:
    assert lookup.clean_term("  IME2 ", "name") == "IME2"
    assert lookup.clean_term('kinase "x" (y)', "name") == "kinase x y"  # quotes/parens out
    assert lookup.clean_term("x" * 200, "name") == "x" * lookup.MAX_TERM
    with pytest.raises(SearchError, match="enter a species"):
        lookup.clean_term('"()"', "species")


def test_the_query_cannot_be_broken_out_of() -> None:
    q = lookup.build_query('IME2" OR reviewed:true OR "', "Saccharomyces cerevisiae", False)
    assert '(gene:"IME2 OR reviewed' in q  # the quotes are gone, so it stays one term
    assert q.count('"') == 6  # three quoted terms (gene, protein_name, organism), no more
    assert "reviewed:true" not in q  # the ":" was stripped: no injected operator survives
    assert 'organism_name:"Saccharomyces cerevisiae"' in q
    assert lookup.build_query("IME2", "yeast", True).endswith(" AND reviewed:true")


def test_species_is_optional() -> None:
    q = lookup.build_query("IME2", "", True)
    assert "organism_name" not in q  # no species: UniProt matches across organisms
    assert q == '(gene:"IME2" OR protein_name:"IME2") AND reviewed:true'


def test_parse_reads_entries_reviewed_first_and_drops_the_sequenceless() -> None:
    found = lookup.parse(json.dumps(UNIPROT_JSON).encode())
    assert [f.accession for f in found] == ["P32581", "A0A1B2C3"]  # reviewed first
    best = found[0]
    assert best.reviewed and best.gene == "IME2" and best.length == 10
    assert best.protein == "Meiosis induction protein kinase IME2"
    assert best.sequence == "MKLVAKQRST" and best.organism == "Saccharomyces cerevisiae"
    assert best.job_name == "IME2 (Saccharomyces cerevisiae)"
    assert "uniprot.org/uniprotkb/P32581" in best.url
    second = found[1]
    assert not second.reviewed and second.protein == "Protein kinase"  # from submissionNames


def test_search_builds_the_request_and_parses_it(uniprot: list[str]) -> None:
    found = lookup.search("IME2", "Saccharomyces cerevisiae", reviewed_only=True)
    assert len(found) == 2 and found[0].accession == "P32581"
    assert uniprot and uniprot[0].startswith(lookup.BASE_URL)
    assert "reviewed%3Atrue" in uniprot[0] and "format=json" in uniprot[0]


def test_an_unreadable_answer_is_an_error() -> None:
    with pytest.raises(SearchError, match="could not be read"):
        lookup.parse(b"not json")


# ── the page ────────────────────────────────────────────────────────────────────────────


def hx(client: Client, **params: str) -> str:
    response = client.get("/homologs/find-sequence/", params, headers={"HX-Request": "true"})
    return response.content.decode()


def test_the_lookup_section_is_on_the_homologs_page(client: Client) -> None:
    page = client.get("/homologs/").content.decode()
    assert "Find a Sequence by Name" in page
    assert 'hx-get="/homologs/find-sequence/"' in page
    assert 'id="species-list"' in page and "Saccharomyces cerevisiae" in page  # the datalist
    assert "Do you have a protein or gene name" in page  # the updated help text
    assert "Start with an amino acid sequence" in page  # the updated Search help text


def test_the_page_has_a_toc_and_a_collapsible_lookup(client: Client) -> None:
    page = client.get("/homologs/").content.decode()
    # The "On this page" nav links to each section that is present.
    assert 'data-toc aria-label="On this page"' in page
    for anchor in ("#find-sequence", "#search"):
        assert f'href="{anchor}"' in page
    assert 'id="find-sequence"' in page and 'id="search"' in page
    # The lookup section folds (fold.js): a heading button controls its body.
    assert 'class="fold-toggle" aria-expanded="true" aria-controls="find-sequence-body"' in page
    assert 'class="fold-body" id="find-sequence-body"' in page


def test_the_consents_are_named_and_explained_in_a_popup(client: Client) -> None:
    page = client.get("/homologs/").content.decode()
    # The lookup consent: new wording, a distinct id (no collision with the search consent),
    # and a "?" that opens the explanation.
    assert "Consent to send search to UniProt" in page
    assert 'id="lookup_consent"' in page and 'popovertarget="help-lookup-consent"' in page
    assert 'id="help-lookup-consent"' in page and "Your own sequences are never sent" in page
    # The search consent: the service is named (filled in by JS) with its own "?" popup.
    assert "Consent to send my sequence to <span data-search-service>" in page
    assert 'popovertarget="help-search-consent"' in page
    assert "Agree to send the residues of your sequence" in page  # moved into the popup


def test_a_lookup_lists_hits_to_use(client: Client, uniprot: list[str]) -> None:
    html = hx(client, protein="IME2", species="Saccharomyces cerevisiae", consent="on")
    assert "P32581" in html and "Meiosis induction protein kinase IME2" in html
    assert "reviewed" in html  # the Swiss-Prot badge
    assert 'class="link use-sequence"' in html and 'data-sequence="MKLVAKQRST"' in html
    assert 'data-name="IME2 (Saccharomyces cerevisiae)"' in html


def test_a_lookup_without_a_species_still_searches(client: Client, uniprot: list[str]) -> None:
    html = hx(client, protein="IME2", consent="on")  # no species
    assert "P32581" in html  # UniProt still returns matches
    # No organism clause in the query (organism_name without a colon is just a return field).
    assert uniprot and "organism_name%3A" not in uniprot[0]


def test_without_consent_the_lookup_is_refused(client: Client, uniprot: list[str]) -> None:
    html = hx(client, protein="IME2", species="yeast")  # no consent
    assert "Tick this box" in html
    assert not uniprot  # nothing was sent to UniProt


def test_a_network_failure_is_shown(client: Client, monkeypatch: pytest.MonkeyPatch) -> None:
    def boom(url: str, data: Any = None) -> bytes:
        raise SearchError("could not reach the search service")

    monkeypatch.setattr(remote, "request", boom)
    html = hx(client, protein="IME2", species="yeast", consent="on")
    assert "could not reach the search service" in html


def test_no_matches_says_so(client: Client, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(remote, "request", lambda url, data=None: b'{"results": []}')
    html = hx(client, protein="Nonesuch", species="yeast", consent="on")
    assert "No matches at UniProt" in html


# ── type-ahead suggestions ──────────────────────────────────────────────────────────────

SPECIES_JSON = {
    "results": [
        {"scientificName": "Homo sapiens", "commonName": "Human"},
        {"scientificName": "Homo sapiens neanderthalensis", "commonName": ""},
        {"scientificName": "Homo sapiens", "commonName": "Human"},  # duplicate
    ]
}
GENES_JSON = {
    "results": [
        {"genes": [{"geneName": {"value": "IME2"}}]},
        {"genes": [{"geneName": {"value": "IME1"}}, {"geneName": {"value": "IME2"}}]},
        {"genes": [{"geneName": {"value": "RIM11"}}]},  # does not start with the prefix
    ]
}


def test_species_suggestions_match_latin_and_common(monkeypatch: pytest.MonkeyPatch) -> None:
    urls = answer(monkeypatch, SPECIES_JSON)
    items = lookup.suggest_species("human")  # a common name reaches UniProt's taxonomy
    values = [s.value for s in items]
    assert values == ["Homo sapiens", "Homo sapiens neanderthalensis"]  # deduped
    assert items[0].label == "Homo sapiens (Human)"  # the common name rides in the label
    assert urls and "taxonomy/search" in urls[0]
    assert lookup.suggest_species("h") == ()  # too short to ask


def test_gene_suggestions_filter_by_prefix_and_scope_to_species(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    urls = answer(monkeypatch, GENES_JSON)
    items = lookup.suggest_genes("IME", "Saccharomyces cerevisiae")
    assert [s.value for s in items] == ["IME1", "IME2"]  # prefix-filtered, deduped, sorted
    assert "gene%3AIME%2A" in urls[0]  # gene:IME*
    assert "Saccharomyces+cerevisiae" in urls[0] or "Saccharomyces%20cerevisiae" in urls[0]


def test_suggestions_never_raise(monkeypatch: pytest.MonkeyPatch) -> None:
    def boom(url: str, data: Any = None) -> bytes:
        raise SearchError("down")

    monkeypatch.setattr(remote, "request", boom)
    assert lookup.suggest_species("Homo") == () and lookup.suggest_genes("IME", "") == ()


def test_the_suggest_endpoint_returns_json(client: Client, monkeypatch: pytest.MonkeyPatch) -> None:
    answer(monkeypatch, SPECIES_JSON)
    response = client.get("/homologs/suggest/", {"kind": "species", "q": "Homo"})
    assert response.status_code == 200
    data = response.json()
    assert data["items"][0] == {"value": "Homo sapiens", "label": "Homo sapiens (Human)"}
    # Unknown kind and short prefixes are simply empty, never an error.
    assert client.get("/homologs/suggest/", {"kind": "nonsense", "q": "x"}).json() == {"items": []}


def test_the_page_loads_the_type_ahead(client: Client) -> None:
    page = client.get("/homologs/").content.decode()
    assert 'src="/static/phylobot/sequence_lookup.js"' in page
    assert 'data-suggest="species"' in page and 'data-suggest="gene"' in page
    assert 'id="gene-list"' in page and 'data-suggest-url="/homologs/suggest/"' in page
