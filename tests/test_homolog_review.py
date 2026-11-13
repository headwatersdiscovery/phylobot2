"""Reviewing a search: hand-picked includes and excludes, full-length sequences, FASTA
export and downloads, and taking the selection into a new job."""

import csv
import io
import re
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.test import Client

from phylobot import searches
from phylobot.homologs import examples, export, filters, remote, sequences, sources
from phylobot.homologs.hits import Hit
from phylobot.homologs.remote import Query
from phylobot.models import HomologSearch
from phylobot.modes import local_user
from tests.conftest import FakeStorage


def hit(hid: str, **fields: Any) -> Hit:
    defaults: dict[str, Any] = dict(
        id=hid, name="kinase", taxon="Candida albicans (strain SC5314)", taxid=1,
        members=1, length=100, evalue=1e-50, bits=200.0, identity=60.0, coverage=95.0,
    )  # fmt: skip
    return Hit(**{**defaults, **fields})


# --- selection ---------------------------------------------------------------------


def test_choices_override_the_filters() -> None:
    hits = [hit("a"), hit("b", coverage=40), hit("c"), hit("d", identity=10)]
    s = filters.select(hits, filters.Settings(), include=["b", "zzz"], exclude=["c"])
    assert [h.id for h in s.by_filters] == ["a", "c"]
    assert [h.id for h in s.selected] == ["a", "b"]  # search order kept
    assert ([h.id for h in s.included], [h.id for h in s.excluded]) == (["b"], ["c"])
    assert s.stages[-1].name == "choices" and s.stages[-1].remaining == 2


def test_choices_that_agree_with_the_filters_change_nothing() -> None:
    hits = [hit("a"), hit("b", coverage=40)]
    s = filters.select(hits, filters.Settings(), include=["a"], exclude=["b"])
    assert (s.included, s.excluded) == ([], [])
    assert [st.name for st in s.stages][-1] == "identity"  # no "choices" stage


def test_reason_names_the_first_failing_filter() -> None:
    assert filters.reason(hit("x", coverage=40, identity=10), filters.Settings()).startswith(
        "fails covers"
    )
    assert filters.reason(hit("x"), filters.Settings()) == ""


# --- export ------------------------------------------------------------------------


def test_fasta_names_are_readable_unique_and_safe() -> None:
    hits = [
        hit("UniRef90_P32581", taxon="Saccharomyces"),
        hit("UniRef90_Q6C4M5", taxon="Yarrowia lipolytica (strain CLIB 122 / E 150)"),
        hit("XP_1", taxon="Yarrowia lipolytica"),
        hit("XP_1b", taxon="Yarrowia lipolytica", name="evil>\nname\t"),
    ]
    named = export.names(hits)
    assert named["UniRef90_P32581"] == "Saccharomyces_P32581"
    assert named["UniRef90_Q6C4M5"] == "Yarrowia_lipolytica_Q6C4M5"
    assert len(set(named.values())) == 4
    fasta, missing = export.to_fasta(hits, {h.id: "MKT" * 30 for h in hits[:3]})
    assert [h.id for h in missing] == ["XP_1b"]
    headers = [line for line in fasta.splitlines() if line.startswith(">")]
    assert headers[0] == ">Saccharomyces_P32581 kinase [Saccharomyces] UniRef90_P32581"
    assert all(len(line) <= export.WIDTH for line in fasta.splitlines() if line[0] != ">")
    evil, _ = export.to_fasta(hits[3:], {"XP_1b": "MKT"})
    assert evil.count(">") == 1 and "\t" not in evil  # no header injection


# --- fetching full sequences -------------------------------------------------------


def test_fetch_uniref_routes_by_representative(monkeypatch: pytest.MonkeyPatch) -> None:
    asked: list[str] = []

    def fake(url: str, data: dict[str, str] | None = None) -> bytes:
        asked.append(url)
        if "/uniprotkb/accessions" in url:
            return b">sp|P32581|IME2_YEAST x\nMKTA\n>tr|Q9UPZ9-2|X y\nMKTC\n"
        if "/uniparc/upis" in url:
            return b">UPI001AF2DC2D status=active\nMKTB\n"
        return b">UniRef90_A0A1 n=1\nMKTD\n"  # the slow fallback

    monkeypatch.setattr(remote, "request", fake)
    ids = ["UniRef90_P32581", "UniRef90_UPI001AF2DC2D", "UniRef90_Q9UPZ9-2", "UniRef90_A0A1",
           "UniRef90_bad id", "../x"]  # fmt: skip
    found = sequences.fetch_uniref(ids)
    assert found == {
        "UniRef90_P32581": "MKTA",
        "UniRef90_UPI001AF2DC2D": "MKTB",
        "UniRef90_Q9UPZ9-2": "MKTC",
        "UniRef90_A0A1": "MKTD",
    }
    assert [u.split("?")[0].rsplit("/", 2)[-2:] for u in asked] == [
        ["uniprotkb", "accessions"],
        ["uniparc", "upis"],
        ["uniref", "stream"],
    ]


def test_fetch_swissprot_hits_by_entry_name(monkeypatch: pytest.MonkeyPatch) -> None:
    """EBI names Swiss-Prot hits by entry name, not UniRef ID: they used to get no
    sequences at all, so "Use these 0 sequences"."""
    asked: list[str] = []

    def fake(url: str, data: dict[str, str] | None = None) -> bytes:
        asked.append(url)
        return b">sp|P25515|VATL1_YEAST V-type\nMTEL\n>sp|P23380|VATL_DROME V\nMSSE\n"

    monkeypatch.setattr(remote, "request", fake)
    found = sources.SOURCES["ebi"].fetch_sequences(["VATL1_YEAST", "VATL_DROME", "bad id"], "")
    assert found == {"VATL1_YEAST": "MTEL", "VATL_DROME": "MSSE"}
    assert len(asked) == 1 and "/uniprotkb/stream?" in asked[0] and "id%3AVATL1_YEAST" in asked[0]


def test_fetch_ncbi_matches_versions_and_retries(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = []

    def fake(url: str, data: dict[str, str] | None = None) -> bytes:
        calls.append(data)
        if len(calls) == 1:
            raise remote.SearchError("timeout")
        return b">NP_012429.1 protein kinase IME2 [S. cerevisiae]\nMKTA*\n>X_2.1 bad\nMK1\n"

    monkeypatch.setattr(remote, "request", fake)
    monkeypatch.setattr("phylobot.homologs.sequences.time.sleep", lambda s: None)
    assert sequences.fetch_ncbi(["NP_012429", "X_2", "no;pe"]) == {"NP_012429": "MKTA"}
    assert calls[-1] is not None and calls[-1]["id"] == "NP_012429,X_2"


def test_examples_have_full_sequences() -> None:
    ebi = examples.sequences("ime2-ebi")
    assert len(ebi) == 1000 and len(ebi["UniRef90_P32581"]) == 645
    assert len(examples.sequences("ime2-ncbi")) == 1996  # 4 withdrawn at NCBI


# --- pages -------------------------------------------------------------------------

pytestmark = pytest.mark.django_db
BASE = (
    "/homologs/examples/ime2-ebi/?max_evalue=1e-05&min_coverage=70&min_identity=30&max_identity=90"
)


def counts(page: str) -> list[str]:
    return re.findall(r'class="count">([\d,]+)<', page)


def test_review_has_sliders_actions_and_both_tables(client: Client) -> None:
    page = client.get(BASE).content.decode()
    assert page.count('type="range"') == 3
    assert "Use these 63 sequences in a new job" in page
    assert "Search results: 63 selected /\n        937 not selected" in page
    assert page.count('class="pick"') == 1000 and page.count('<table class="hits">') == 1


def test_hand_picked_choices(client: Client) -> None:
    yarrowia = "UniRef90_Q6C329"  # Yarrowia IME2: removed by coverage (57%)
    top_kept = "UniRef90_A0A6C1E9R5"
    page = client.get(f"{BASE}&include={yarrowia}&exclude={top_kept}").content.decode()
    assert counts(page)[-1] == "63"  # 62 + query, then +1 -1
    assert "after your own choices" in page and "1 included" in page and "1 excluded" in page
    assert "included by you" in page and "excluded by you" in page
    assert f'<input type="hidden" name="include" value="{yarrowia}">' in page
    reset = re.search(r'<a class="head-action" href="([^"]+)" data-reset', page)
    assert reset and "include=" not in reset.group(1) and "min_coverage=70" in reset.group(1)
    # The note says why the filters removed the hit the user put back.
    row = page.split(f'data-id="{yarrowia}"')[1].split("</tr>")[0]
    assert "included by you" in row and "fails covers" in row


def ids_in_order(page: str) -> list[str]:
    return re.findall(r'class="pick" data-id="([^"]+)"', page)


def test_sorting(client: Client) -> None:
    search_order = ids_in_order(client.get(BASE).content.decode())
    page = client.get(f"{BASE}&sort=identity").content.decode()
    identity = [
        float(x) for x in re.findall(r'<td class="num">(\d+)%</td><td class="num">\d+\.\d', page)
    ]
    assert identity == sorted(identity) and ids_in_order(page) != search_order
    assert 'aria-sort="ascending"' in page.split('data-sort="identity"')[0].rsplit("<th", 1)[1]
    down = ids_in_order(client.get(f"{BASE}&sort=-rank").content.decode())
    assert down == search_order[::-1]
    use = client.get(f"{BASE}&sort=use").content.decode()
    ticks = re.findall(r'aria-label="Use [^"]+"( checked)?>', use)
    assert ticks[:63] == [" checked"] * 63 and " checked" not in ticks[63:]
    assert ids_in_order(client.get(f"{BASE}&sort=bogus").content.decode()) == search_order
    note = client.get(f"{BASE}&sort=note").content.decode().split("<tbody>")[1]
    rows = note.split("</tr>")
    assert "fails" in rows[0] and "fails" not in rows[-2]  # notes first, blank ones last


def test_select_all_none_and_reset(client: Client) -> None:
    everything = client.get(f"{BASE}&base=all").content.decode()
    assert counts(everything)[-1] == "1,000" and "after your own choices" in everything
    assert '<input type="hidden" name="base" value="all">' in everything
    assert "Use these 1,000 sequences" in everything
    one_out = client.get(f"{BASE}&base=all&exclude=UniRef90_A0A6C1E9R5").content.decode()
    assert counts(one_out)[-1] == "999"
    nothing = client.get(f"{BASE}&base=none&include=UniRef90_Q6C329").content.decode()
    assert counts(nothing)[-1] == "1" and "62 excluded" in nothing and "1 included" in nothing
    reset = re.search(r'<a class="head-action" href="([^"]+)" data-reset', nothing)
    assert reset and "base=" not in reset.group(1) and "include=" not in reset.group(1)


def test_csv(client: Client) -> None:
    got = client.get(
        f"/homologs/examples/ime2-ebi/download/csv/?{BASE.split('?')[1]}&sort=-identity"
    )
    assert got["Content-Disposition"].endswith('.hits.csv"')
    rows = list(csv.reader(io.StringIO(got.content.decode())))
    assert rows[0][:3] == ["selected", "rank", "id"] and len(rows) == 1001
    assert sum(r[0] == "yes" for r in rows[1:]) == 63
    identity = [float(r[7]) for r in rows[1:]]
    assert identity == sorted(identity, reverse=True)
    assert any(r[-1].startswith("fails") for r in rows[1:])


def test_uniprotkb_descriptions() -> None:
    from phylobot.homologs.hits import parse_description

    text = ("V-type proton ATPase subunit c OS=Saccharomyces cerevisiae (strain ATCC 204508 / "
            "S288c) OX=559292 GN=VMA3 PE=1 SV=1")  # fmt: skip
    assert parse_description(text) == (
        "V-type proton ATPase subunit c",
        "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)",
        559292,
        None,
    )
    assert parse_description("something else") == ("something else", "", None, None)


def test_swissprot_hits_link_to_uniprotkb() -> None:
    from phylobot.homolog_views import hit_link

    assert hit_link("ebi", "VATL1_YEAST") == "https://www.uniprot.org/uniprotkb/VATL1_YEAST"
    assert hit_link("ebi", "UniRef90_P32581") == "https://www.uniprot.org/uniref/UniRef90_P32581"


def test_downloads(client: Client) -> None:
    selected = client.get(f"/homologs/examples/ime2-ebi/download/selected/?{BASE.split('?')[1]}")
    assert selected["Content-Disposition"].endswith('.selected.fasta"')
    text = selected.content.decode()
    assert text.count(">") == 63 and ">Saccharomyces_A0A6C1E9R5 " in text
    everything = client.get("/homologs/examples/ime2-ebi/download/all/").content.decode()
    assert everything.count(">") == 1000
    raw = client.get("/homologs/examples/ime2-ebi/download/raw/")
    assert raw.content == examples.raw("ime2-ebi")
    assert client.get("/homologs/examples/ime2-ebi/download/other/").status_code == 404
    assert client.get("/homologs/examples/nope/download/all/").status_code == 404


def test_use_in_a_new_job(client: Client) -> None:
    page = client.get(BASE).content.decode()
    use = re.search(r'<a class="button primary" href="([^"]+)"', page)
    assert use
    url = use.group(1).replace("&amp;", "&")
    job_page = client.get(url).content.decode()
    assert "Step 2 holds 63 sequences from" in job_page
    textarea = re.search(r"<textarea[^>]*>(.*?)</textarea>", job_page, re.S)
    assert textarea and textarea.group(1).count("&gt;") == 63
    assert "Yarrowia" not in textarea.group(1)  # removed by the 70% coverage filter
    assert 'value="Saccharomyces_A0A6C1E9R5"' in job_page  # the outgroup checklist


@pytest.fixture
def finished(fake_storage: FakeStorage) -> HomologSearch:
    search = searches.create(local_user(), Query("M" * 30), "ebi", "uniref90", name="Ime2")
    fake_storage.put_bytes(search.result_key, examples.raw("ime2-ebi"))
    fasta = sequences.to_fasta(examples.sequences("ime2-ebi"))
    fake_storage.put_bytes(search.sequences_key, fasta.encode())
    HomologSearch.objects.filter(id=search.id).update(status=HomologSearch.Status.SUCCEEDED)
    return search


def test_real_search_downloads_and_job(client: Client, finished: HomologSearch) -> None:
    got = client.get(f"/homologs/{finished.id}/download/selected/").content.decode()
    assert got.count(">") == 63
    page = " ".join(
        client.get(f"/jobs/new/?homologs={finished.id}&min_coverage=55").content.decode().split()
    )
    assert "Step 2 holds 627 sequences from Ime2" in page
    assert "more than your limit of 500" in page


def test_others_cannot_download_or_prefill(client: Client, finished: HomologSearch) -> None:
    HomologSearch.objects.filter(id=finished.id).update(owner=User.objects.create_user("other"))
    assert client.get(f"/homologs/{finished.id}/download/all/").status_code == 404
    assert client.get(f"/jobs/new/?homologs={finished.id}").status_code == 404
