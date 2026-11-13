"""The New Job page's reconciliation block (species_web, SubmitForm, the species table):
species from the headers, the user's corrections, the consent gate, and Submit running
the job with or without reconciliation. Open Tree answers are replayed (conftest's
no_open_tree refuses any other request)."""

import re
from pathlib import Path
from typing import Any

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from django.http import QueryDict
from django.test import Client

from phylobot import species_web
from phylobot.homologs import export
from phylobot.homologs.hits import Hit
from phylobot.models import Job
from phylobot.species.opentree import Client as OpenTreeClient
from phylobot.species.table import OK, read
from tests.conftest import FakeStorage
from tests.fixtures.opentree import Replay

ROOT = Path(__file__).resolve().parent.parent
IME2 = (ROOT / "src" / "phylobot" / "examples" / "ime2.fasta").read_text()
IME2_NAMES = re.findall(r"^>(\S+)", IME2, re.M)
# The Ime2 species as an uploaded tree (names with underscores, one extra to prune).
IME2_TREE = (
    "((((Saccharomyces_cerevisiae,Kluyveromyces_waltii),Kluyveromyces_lactis),"
    "((Candida_albicans,Candida_parapsilosis),Candida_lusitaniae)),"
    "(Yarrowia_lipolytica,Homo_sapiens));"
)


def recorded(monkeypatch: pytest.MonkeyPatch, name: str = "ime2") -> Replay:
    replay = Replay(name)
    monkeypatch.setattr(species_web, "make_transport", lambda: replay)
    species_web.reset_client()
    return replay


def test_rows_read_headers_then_carried_taxids_then_the_users_cells() -> None:
    fasta = b">Saccharomyces.cerevisiae.IME2\nMK\n>q2 kinase [Candida albicans]\nMK\n>x\nMK\n"
    rows = species_web.rows_for(fasta)
    assert [(r.query, r.source) for r in rows] == [
        ("Saccharomyces cerevisiae", "name"),
        ("Candida albicans", "ncbi"),
        ("", ""),
    ]
    rows = species_web.rows_for(fasta, {"x": ("Homo sapiens", None)}, {"q2": 5476})
    assert (rows[1].taxid, rows[1].source) == (5476, "ncbi")
    assert (rows[2].query, rows[2].source) == ("Homo sapiens", "you")
    assert species_web.rows_for(fasta, carried={"x": 9606})[2].source == "search"
    # a cell posted back unchanged is not an edit
    echoed = {"Saccharomyces.cerevisiae.IME2": ("Saccharomyces cerevisiae", None)}
    same = species_web.rows_for(fasta, echoed)
    assert same[0].source == "name"
    assert species_web.rows_for(b"not fasta") == []


def test_cells_are_cleaned_as_hostile() -> None:
    data = QueryDict(mutable=True)
    data.setlist("species_for", ["a", "b", "c"])
    data.setlist("species_query", ["Homo  sapiens", "x" * 500, "Mus musculus"])
    data.setlist("species_taxid", ["9606", "1", "ten"])
    assert species_web.edits(data) == {"a": ("Homo sapiens", 9606)}  # b too long, c not a number
    data.setlist("species_taxid", ["9606"])  # lists of different lengths: ignored
    assert species_web.edits(data) == {}
    assert species_web.hints("a\t9606\nb\tnope\n\t5\nc\t4932\n") == {"a": 9606, "c": 4932}


def test_nothing_is_sent_without_consent() -> None:
    rows = species_web.rows_for(IME2.encode())
    found = species_web.look_up(rows, consent=False)  # conftest refuses any request
    assert found.state == species_web.NEEDS_CONSENT and found.tree is None
    assert "unticked" in found.skipped
    assert species_web.look_up([], consent=True).state == species_web.NO_SEQUENCES


def test_open_tree_with_consent(monkeypatch: pytest.MonkeyPatch) -> None:
    recorded(monkeypatch)
    found = species_web.look_up(species_web.rows_for(IME2.encode()), consent=True)
    assert found.ready and found.skipped == "" and found.species == 7
    assert found.provenance["synth_id"] == "opentree16.1"
    assert all(row.status == OK for row in found.rows)


def test_an_uploaded_tree_needs_no_consent_and_is_pruned() -> None:
    rows = species_web.rows_for(IME2.encode())
    found = species_web.look_up(rows, consent=False, upload=IME2_TREE)
    assert found.ready and found.tree is not None and found.species == 7
    assert "Homo_sapiens" not in found.tree.newick()
    bad = species_web.look_up(rows, consent=False, upload="((a,b);")
    assert bad.state == species_web.FAILED and "species tree could not be used" in bad.skipped


def test_unresolved_rows_say_which(monkeypatch: pytest.MonkeyPatch) -> None:
    rows = species_web.rows_for(IME2.encode() + b">mystery\nMKV\n")
    found = species_web.look_up(rows, consent=False, upload=IME2_TREE)
    assert found.state == species_web.INCOMPLETE
    assert [row.sequence for row in found.pending] == ["mystery"]
    assert found.skipped == "1 sequence(s) had no species in the species tree (mystery)"


def test_homolog_hits_carry_their_taxids() -> None:
    hit = Hit("UniRef90_P32581", "IME2", "Saccharomyces cerevisiae", 4932, 3, 645, 0, 900, 99, 100)
    bare = Hit("XP_1", "IME2-like", "", None, None, 600, 0, 500, 40, 90)
    assert export.taxids([hit, bare]) == "Saccharomyces_cerevisiae_P32581\t4932\n"


# The page and Submit.

pytestmark = pytest.mark.django_db


def boxes_of(page: str) -> dict[str, str]:
    """The reconciliation block's checkbox tags by name."""
    names = r'name="(reconcile|species_consent|reconciled_tree)"'
    return {n: tag for tag in re.findall(r"<input[^>]*>", page) for n in re.findall(names, tag)}


def post(client: Client, **fields: Any) -> Any:
    data = {"aligners": ["mafft"], "tree_builders": ["iqtree"], "modelfinder": "on",
            "run_trees": "on", "run_ancestors": "on", "bootstrap": 1000, "seed": 1,
            "sequences_text": IME2, "reconcile": "on", **fields}  # fmt: skip
    return client.post("/jobs/new/", data)


def test_the_block_is_ticked_and_consent_is_not(client: Client) -> None:
    page = client.get("/jobs/new/").content.decode()
    boxes = boxes_of(page)
    assert "checked" in boxes["reconcile"]  # owner decision 2026-10-17
    assert "checked" not in boxes["reconciled_tree"]
    # the consent box comes with the species table, unticked
    table = client.post("/species-check/", {"sequences_text": IME2, "reconcile": "on"})
    text = table.content.decode()
    assert "checked" not in boxes_of(text)["species_consent"]
    assert "Get the species tree from the Open Tree of Life" in text


def test_species_table_without_consent(client: Client) -> None:
    page = client.post("/species-check/", {"sequences_text": IME2, "reconcile": "on"})
    text = page.content.decode()
    assert text.count('name="species_for"') == 7 and "Species found for 7 of 7" in text
    assert 'value="Saccharomyces cerevisiae"' in text and "Get the species tree" in text
    off = client.post("/species-check/", {"sequences_text": IME2}).content.decode()
    assert "species_for" not in off


def test_species_table_with_consent_draws_the_tree(
    client: Client, monkeypatch: pytest.MonkeyPatch
) -> None:
    recorded(monkeypatch)
    text = client.post(
        "/species-check/", {"sequences_text": IME2, "reconcile": "on", "species_consent": "on"}
    ).content.decode()
    assert "All 7 sequences have a species" in text and "opentree16.1" in text
    assert "Lachancea waltii" in text  # Open Tree's current name for Kluyveromyces waltii
    assert 'aria-label="Species tree, 7 species"' in text and "Hinchliff" in text


def test_submit_without_consent_runs_unreconciled_and_says_why(
    client: Client, fake_storage: FakeStorage
) -> None:
    post(client, reconciled_tree="on")
    job = Job.objects.get()
    assert "reconcile" not in job.params and "unticked" in job.params["reconcile_skipped"]
    assert job.params["tree_builders"] == ["iqtree"]  # no GeneRax without a species tree
    assert not any("species" in key for key in fake_storage.objects)
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "skip reconciliation with the species tree: the species names were not sent" in page


def test_submit_with_consent_stores_the_species_and_tree(
    client: Client, fake_storage: FakeStorage, monkeypatch: pytest.MonkeyPatch
) -> None:
    recorded(monkeypatch)
    post(client, species_consent="on", reconciled_tree="on", polytomies="arbitrary")
    job = Job.objects.get()
    assert job.params["reconcile"] is True and job.params["polytomies"] == "arbitrary"
    assert job.params["tree_builders"] == ["iqtree", "generax"]
    assert job.params["species"]["synth_id"] == "opentree16.1"
    assert job.params["species"]["species"] == 7
    inputs = f"jobs/{job.id}/inputs/"
    rows = read(fake_storage.objects[inputs + "species.tsv"].decode())
    assert [row.sequence for row in rows] == IME2_NAMES and all(r.status == OK for r in rows)
    assert fake_storage.objects[inputs + "species_tree.nwk"].endswith(b";\n")
    assert b"opentree16.1" in fake_storage.objects[inputs + "species_provenance.json"]


def test_submit_with_an_uploaded_species_tree_and_a_correction(
    client: Client, fake_storage: FakeStorage
) -> None:
    fasta = IME2 + ">mystery\n" + "MKV" * 40 + "\n"
    post(
        client, sequences_text=fasta,
        species_tree_file=SimpleUploadedFile("species.nwk", IME2_TREE.encode()),
        species_for=["mystery"], species_query=["Homo sapiens"], species_taxid=[""],
    )  # fmt: skip
    job = Job.objects.get()
    assert job.params["reconcile"] is True and job.params["species"]["source"] == "upload"
    stored = fake_storage.objects[f"jobs/{job.id}/inputs/species.tsv"].decode()
    rows = {r.sequence: r for r in read(stored)}
    assert (rows["mystery"].species, rows["mystery"].source) == ("Homo sapiens", "you")


def test_unticked_means_no_reconciliation_and_nothing_said(
    client: Client, fake_storage: FakeStorage
) -> None:
    post(client, reconcile="", species_consent="on", reconciled_tree="on")
    job = Job.objects.get()
    assert "reconcile" not in job.params and "reconcile_skipped" not in job.params
    assert job.params["tree_builders"] == ["iqtree"]


def test_reconciled_tree_needs_amino_acid_models(client: Client) -> None:
    lysin = (ROOT / "src" / "phylobot" / "examples" / "lysin.fasta").read_text()
    response = post(
        client, sequences_text=lysin, seqtype="codon", analyses="codon",
        codon_modelfinder="on", codon_search="on", reconciled_tree="on",
    )  # fmt: skip
    assert "needs amino-acid models" in response.content.decode() and not Job.objects.exists()


def test_preview_describes_reconciliation(client: Client) -> None:
    base = {"sequences_text": IME2, "aligners": ["mafft"], "tree_builders": ["iqtree"],
            "modelfinder": "on", "run_trees": "on", "run_ancestors": "on",
            "reconcile": "on"}  # fmt: skip
    text = client.post("/preview/", base).content.decode()
    assert "Skip reconciliation with the species tree" in text
    response = client.post("/preview/", base | {"species_consent": "on", "reconciled_tree": "on"})
    page = response.content.decode()
    assert "Reconcile every tree with the species tree" in page
    assert "Search a reconciled tree with GeneRax" in page


def test_duplicate_keeps_the_choices_but_asks_consent_again(
    client: Client, fake_storage: FakeStorage, monkeypatch: pytest.MonkeyPatch
) -> None:
    recorded(monkeypatch)
    post(client, species_consent="on", reconciled_tree="on", polytomies="arbitrary")
    job = Job.objects.get()
    page = client.get(f"/jobs/new/?from={job.id}").content.decode()
    boxes = boxes_of(page)
    assert "checked" in boxes["reconcile"] and "checked" in boxes["reconciled_tree"]
    assert "checked" not in boxes.get("species_consent", "")  # asked again, in the table
    assert re.search(r'value="arbitrary"[^>]*checked', page)
    assert 'name="tree_builders" value="iqtree"' in page and 'value="generax"' not in page


def test_homolog_hand_off_fills_the_taxids(client: Client) -> None:
    page = client.get("/jobs/new/?homologs_example=ime2-ebi").content.decode()
    hints = re.search(r'name="species_hints" value="([^"]*)"', page)
    assert hints and re.search(r"^\w+\t\d+$", hints.group(1), re.M)


def test_open_tree_client_is_shared_and_renewed(monkeypatch: pytest.MonkeyPatch) -> None:
    recorded(monkeypatch)
    first = species_web.client()
    assert species_web.client() is first and isinstance(first, OpenTreeClient)
    monkeypatch.setattr(species_web, "SYNTHESIS_SECONDS", -1)
    assert species_web.client() is not first
