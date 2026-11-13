"""Species without typing them: UniProt species codes (offline), accessions looked up at
NCBI and UniProt (with consent), species lists, and offers from a row with the same name
prefix (phylobot.species.codes, .accessions, species_web). The canned answers below are
the services' real answers for public Ime2 accessions (2026-10-17), trimmed."""

import json
import re
from typing import Any

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from django.http import QueryDict
from django.test import Client

from phylobot import species_web
from phylobot.models import Job
from phylobot.species import accessions, headers, known
from phylobot.species.accessions import Accession, AccessionError
from phylobot.species.table import OK, Row, read
from tests.conftest import FakeStorage

NCBI_ANSWER = {
    "header": {"type": "esummary", "version": "0.3"},
    "result": {
        "uids": ["6322355", "2480749293"],
        "6322355": {"caption": "NP_012429", "accessionversion": "NP_012429.1", "taxid": 559292,
                    "organism": "Saccharomyces cerevisiae S288C"},
        "2480749293": {"caption": "CAI6736627", "accessionversion": "CAI6736627.1",
                       "taxid": 4932, "organism": "Saccharomyces cerevisiae"},
    },
}  # fmt: skip
UNIPROT_ANSWER = (
    "Entry\tOrganism\tOrganism (ID)\n"
    "A7TQK5\tVanderwaltozyma polyspora (strain ATCC 22028 / DSM 70294) (Kluyveromyces "
    "polysporus)\t436907\n"
)


class Services:
    """A transport answering like NCBI and UniProt; it records what was sent."""

    def __init__(self) -> None:
        self.sent: list[tuple[str, Any]] = []

    def __call__(self, url: str, data: dict[str, str] | str | None) -> bytes:
        self.sent.append((url, data))
        if url.startswith(accessions.NCBI_URL):
            return json.dumps(NCBI_ANSWER).encode()
        if url.startswith(accessions.UNIPROT_URL):
            return UNIPROT_ANSWER.encode()
        raise AccessionError(f"unexpected {url}")


@pytest.mark.parametrize(
    ("name", "description", "expected"),
    [
        ("sp|P32581|IME2_YEAST", "Meiosis-inducing kinase", ("Saccharomyces cerevisiae", 4932)),
        ("Q6C4M5_YARLI", "", ("Yarrowia lipolytica", 4952)),  # the species, not the strain
        ("tr|A0A6C1E9R5|A0A6C1E9R5_9SACH", "", ("Saccharomycetaceae", 4893)),  # a group
        ("XP_1", "kinase [organism=Candida albicans] [gene=IME2]", ("Candida albicans", None)),
        ("hit7", "IME2 [Candida tropicalis] UniRef90_X", ("Candida tropicalis", None)),
        ("IME2_NOTACODE", "", None),
        # identifiers that imply their species
        ("9606.ENSP00000269305", "", ("Homo sapiens", 9606)),  # STRING / eggNOG
        ("7227_0:001c7b", "", ("Drosophila melanogaster", 7227)),  # OrthoDB
        ("ENSMUSP00000012345.2", "", ("Mus musculus", 10090)),
        ("ENSP00000269305", "", ("Homo sapiens", 9606)),
        ("AT1G01010.1", "", ("Arabidopsis thaliana", 3702)),
        ("YAL001C", "", ("Saccharomyces cerevisiae", 4932)),
        ("FBpp0070000", "", ("Drosophila melanogaster", 7227)),
        ("WBGene00000001", "", ("Caenorhabditis elegans", 6239)),
        ("PF3D7_0100100", "", ("Plasmodium falciparum", 5833)),
        # a species UniProt knows, anywhere, in any case
        ("IME2_Saccharomyces_cerevisiae", "", ("Saccharomyces cerevisiae", None)),
        ("kinase-homo_sapiens-1", "", ("Homo sapiens", None)),
        ("Protein_kinase_7", "kinase from Mus musculus", ("Mus musculus", None)),
        # the old conventions are untouched (Open Tree resolves the synonym)
        ("Kluyveromyces.waltii.IME2", "", ("Kluyveromyces waltii", None)),
        # guesses: only when one species fits
        ("Scer_IME2", "", ("Saccharomyces cerevisiae", 4932)),
        ("S.cerevisiae_IME2", "", ("Saccharomyces cerevisiae", 4932)),
        ("H_sapiens_X", "", ("Homo sapiens", 9606)),
        ("Gene_1", "", None),
        ("seq1", "", None),
    ],
)
def test_headers_read_codes_and_genbank(name: str, description: str, expected: Any) -> None:
    found = headers.detect(name, description)
    assert (None if found is None else (found.name, found.taxid)) == expected


def test_guesses_say_so() -> None:
    found = headers.detect("Scer_IME2")
    assert found is not None and found.source == "guess" and "'Scer'" in found.note


def test_known_names_common_names_and_abbreviations() -> None:
    assert known.common("human") == known.Known("Homo sapiens", 9606)
    assert known.common("fruit fly") is not None  # several species: the model organism
    assert known.common("Fruit Fly").name == "Drosophila melanogaster"  # type: ignore[union-attr]
    assert known.common("yeast") is None  # too many species to choose from
    assert known.abbreviation("Dmel") == known.Known("Drosophila melanogaster", 7227)
    assert known.abbreviation("Abcd") is None
    assert known.by_taxid(9606) == "Homo sapiens"


def test_a_common_name_in_a_cell_becomes_the_species() -> None:
    data = QueryDict(mutable=True)
    data.setlist("species_for", ["a", "b"])
    data.setlist("species_query", ["zebrafish", ""])
    data.setlist("species_taxid", ["", "10090"])
    data.setlist("species_shown", ["", ""])
    a, b = species_web.gather(b">a\nMK\n>b\nMK\n", data).rows
    assert (a.query, a.taxid, a.note) == ("Danio rerio", 7955, "'zebrafish' is Danio rerio")
    assert (b.query, b.taxid) == ("Mus musculus", 10090)  # a bare taxid gets its name


def test_species_codes_are_bundled_with_their_release() -> None:
    assert "CC BY 4.0" in known.release() and "2026_03" in known.release()
    yeast = known.code("YEAST")
    assert yeast is not None and (yeast.taxid, yeast.group) == (559292, False)
    group = known.code("9SACH")
    assert group is not None and group.group


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("XP_012345.1", Accession("protein", "XP_012345.1")),
        ("WP_020694038.1", Accession("protein", "WP_020694038.1")),
        ("ENSACAP00000000001", Accession("ensembl", "ENSACAP00000000001")),
        ("NP_012429", Accession("protein", "NP_012429")),
        ("CAI6736627", Accession("protein", "CAI6736627")),
        ("XM_012345678.2", Accession("nuccore", "XM_012345678.2")),
        ("sp|P32581|IME2_YEAST", Accession("uniprot", "P32581")),
        ("UniRef90_A7TQK5", Accession("uniprot", "A7TQK5")),
        ("gi|6322355|ref|NP_012429.1|", Accession("protein", "NP_012429.1")),
        ("Scer_IME2", None),
        ("Kluyveromyces.waltii.IME2", None),
    ],
)
def test_accessions_are_recognised(name: str, expected: Accession | None) -> None:
    assert accessions.find(name) == expected


def test_client_batches_parses_and_caches(tmp_path: Any) -> None:
    services = Services()
    client = accessions.Client("me@example.org", services, accessions.Cache(tmp_path))
    wanted = [Accession("protein", "NP_012429"), Accession("protein", "CAI6736627.1"),
              Accession("protein", "XP_000000001"), Accession("uniprot", "A7TQK5")]  # fmt: skip
    found = client.look_up(wanted)
    assert {a.id: (f.organism, f.taxid) for a, f in found.items()} == {
        "NP_012429": ("Saccharomyces cerevisiae S288C", 559292),
        "CAI6736627.1": ("Saccharomyces cerevisiae", 4932),
        "A7TQK5": ("Vanderwaltozyma polyspora", 436907),  # strain notes dropped
    }
    url, data = services.sent[0]
    assert data["tool"] == "phylobot" and data["email"] == "me@example.org"
    assert "XP_000000001" in data["id"]  # only accessions are sent
    asked = len(services.sent)
    client.look_up(wanted[:2] + wanted[3:])  # all cached: nothing sent
    assert len(services.sent) == asked


IPG = (
    "Id\tSource\tNucleotide Accession\tStart\tStop\tStrand\tProtein\tProtein Name\tOrganism"
    "\tStrain\tAssembly\n"
    "27962999\tRefSeq\tNZ_1\t1\t9\t+\tWP_000062431.1\tDUF4942\tEscherichia coli\tA\tG1\n"
    "27962999\tRefSeq\tNZ_2\t1\t9\t+\tWP_000062431.1\tDUF4942\tEscherichia coli\tB\tG2\n"
    "27962999\tINSDC\tCP3\t1\t9\t+\tEQY54909.1\tDUF4942\tKlebsiella pneumoniae\tC\tG3\n"
    "27962999\tINSDC\tCP4\t1\t9\t+\tEQY54910.1\tDUF4942\tKlebsiella sp.\tD\tG4\n"
    "41321826\tRefSeq\tNZ_5\t1\t9\t+\tWP_003131952.1\tX\tLactococcus lactis subsp. lactis\t\tG5\n"
    "41321826\tRefSeq\tNZ_6\t1\t9\t+\tWP_003131952.1\tX\tLactococcus lactis\t\tG6\n"
)


def test_shared_proteins_list_their_species_and_ensembl_names_one(tmp_path: Any) -> None:
    sent: list[Any] = []

    def transport(url: str, data: Any) -> bytes:
        sent.append((url, data))
        if url == accessions.NCBI_URL:
            return json.dumps({"result": {
                "uids": ["1", "2"],
                "1": {"accessionversion": "WP_000062431.1", "taxid": 543,
                      "organism": "Enterobacteriaceae"},
                "2": {"accessionversion": "WP_003131952.1", "taxid": 1357,
                      "organism": "Lactococcus"},
            }}).encode()  # fmt: skip
        if url == accessions.IPG_URL:
            return IPG.encode()
        if url == accessions.ENSEMBL_URL:
            assert json.loads(data) == {"ids": ["ENSACAP00000000001"]}  # a JSON body
            return json.dumps({"ENSACAP00000000001": {"species": "anolis_carolinensis"}}).encode()
        raise AccessionError(url)

    client = accessions.Client("", transport, accessions.Cache(tmp_path))
    wanted = [Accession("protein", "WP_000062431.1"), Accession("protein", "WP_003131952.1"),
              Accession("ensembl", "ENSACAP00000000001")]  # fmt: skip
    found = client.look_up(wanted)
    shared, single, lizard = (found[a] for a in wanted)
    assert shared.organism == "Enterobacteriaceae" and shared.genomes == 4
    assert shared.choices == (("Escherichia coli", 2), ("Klebsiella pneumoniae", 1))  # no "sp."
    assert (single.organism, single.taxid, single.choices) == ("Lactococcus lactis", None, ())
    assert lizard.organism == "Anolis carolinensis"
    asked = len(sent)
    assert client.look_up(wanted) == found and len(sent) == asked  # the cache keeps choices


def test_a_shared_protein_row_asks_the_user_to_choose() -> None:
    choices = (("Escherichia coli", 472), ("Salmonella enterica", 133))
    answer = accessions.Found(Accession("protein", "WP_1"), "Enterobacteriaceae", 543, choices, 745)
    [row] = species_web.rows_for(b">WP_1\nMK\n", found={"WP_1": answer})
    assert row.query == "Enterobacteriaceae" and row.source == "lookup ncbi"
    assert "745 genomes of 2 species: choose one" in row.note and "horizontal" in row.note


def test_lists_read_tabs_commas_headers_and_taxids() -> None:
    assert species_web.bulk("a\tHomo sapiens\nc\t4932\n\nd\n") == {
        "a": ("Homo sapiens", None),
        "c": ("", 4932),
    }  # fmt: skip  (d has no species: skipped)
    mixed = species_web.bulk("a\tHomo sapiens\nb,Mus musculus,10090\n")
    assert mixed == {"a": ("Homo sapiens", None), "b": ("Mus musculus", 10090)}
    table = "sequence\tspecies\ttaxid\tin_species_tree\nx\tHomo sapiens\t9606\tHomo sapiens\n"
    assert species_web.bulk(table) == {"x": ("Homo sapiens", 9606)}
    assert species_web.bulk("Name,Species\nx,Danio rerio\n") == {
        "Name": ("Species", None),
        "x": ("Danio rerio", None),
    }  # no "sequence" column: not a header


def test_cells_count_only_when_changed_from_what_was_shown() -> None:
    data = QueryDict(mutable=True)
    data.setlist("species_for", ["a", "b", "c"])
    data.setlist("species_query", ["Homo sapiens", "Mus musculus", "Danio rerio"])
    data.setlist("species_taxid", ["", "", ""])
    data.setlist("species_shown", ["Homo sapiens\t", "Rattus\t", ""])
    # a: as shown (not the user's); b: changed; c: typed by the user earlier
    assert species_web.edits(data) == {"b": ("Mus musculus", None), "c": ("Danio rerio", None)}


def test_a_list_replaces_what_the_table_showed_but_not_what_was_typed() -> None:
    fasta = b">a\nMK\n>b\nMK\n"
    data = QueryDict(mutable=True)
    data.setlist("species_for", ["a", "b"])
    data.setlist("species_query", ["Homo sapiens", "Mus musculus"])
    data.setlist("species_taxid", ["", ""])
    data.setlist("species_shown", ["Homo sapiens\t", ""])  # b was typed
    data["species_paste"] = "a\tDanio rerio\nb\tGallus gallus\n"
    rows = species_web.gather(fasta, data).rows
    assert [(r.query, r.source) for r in rows] == [("Danio rerio", "list"), ("Mus musculus", "you")]


def test_accessions_are_looked_up_only_with_consent() -> None:
    fasta = b">NP_012429\nMK\n>sp|P32581|IME2_YEAST\nMK\n>UniRef90_A7TQK5\nMK\n>odd\nMK\n"
    services = Services()

    def client() -> accessions.Client:
        return accessions.Client("", services)

    without = species_web.gather(fasta, {}, get_accession_client=client)
    assert not services.sent and not without.accessions.asked
    # the UniProt code gave P32581's species offline: only two accessions are pending
    assert sorted(without.accessions.pending) == ["NP_012429", "UniRef90_A7TQK5"]
    found = species_web.gather(fasta, {"accession_consent": "on"}, get_accession_client=client)
    rows = {r.sequence: r for r in found.rows}
    assert (rows["NP_012429"].taxid, rows["NP_012429"].source) == (559292, "lookup ncbi")
    assert rows["UniRef90_A7TQK5"].query == "Vanderwaltozyma polyspora"
    assert rows["sp|P32581|IME2_YEAST"].source == "uniprot code"
    assert {d["id"] for u, d in services.sent if d} == {"NP_012429"}


def test_lookup_failure_is_reported_not_raised() -> None:
    def broken() -> accessions.Client:
        def transport(url: str, data: Any) -> bytes:
            raise AccessionError("could not reach eutils.ncbi.nlm.nih.gov")

        return accessions.Client("", transport)

    got = species_web.gather(
        b">NP_012429\nMK\n", {"accession_consent": "on"}, get_accession_client=broken
    )
    assert "could not be looked up" in got.accessions.error and got.rows[0].taxid is None


def test_offers_from_a_row_with_the_same_prefix() -> None:
    rows = [
        Row("Scer_IME2", query="Saccharomyces cerevisiae", source="you", status=OK),
        Row("Scer_KSS1"),
        Row("Scer_FUS3", query="Saccharomyces cerevisiae", source="name", status=OK),
        Row("Scer_SMK1"),
        Row("Calb_IME2"),
    ]
    rows += [  # accession prefixes are the database's: never offered between accessions
        Row("WP_000062431.1", query="Escherichia coli", source="you", status=OK),
        Row("WP_020694038.1", status="not in tree"),
    ]
    [offer] = species_web.offers(rows, resolved=True)
    assert (offer.prefix, offer.source) == ("Scer", "Scer_IME2")
    assert offer.query == "Saccharomyces cerevisiae"
    assert offer.rows == ("Scer_KSS1", "Scer_SMK1")


# The page.

pytestmark = pytest.mark.django_db


def test_table_offers_the_lookup_and_reads_codes_offline(client: Client) -> None:
    fasta = ">sp|P32581|IME2_YEAST\nMKV\n>NP_012429.1\nMKV\n>Scer_KSS1\nMKV\n"
    page = client.post("/species-check/", {"sequences_text": fasta, "reconcile": "on"})
    text = page.content.decode()
    assert "UniProt code" in text and 'value="Saccharomyces cerevisiae"' in text
    assert 'value="4932"' in text  # species-level taxid, though YEAST is strain S288c
    assert "Look up 1 accession number at NCBI, UniProt or Ensembl" in " ".join(text.split())
    assert 'name="accession_consent">' in text  # unticked


def test_table_with_accession_consent(client: Client, monkeypatch: pytest.MonkeyPatch) -> None:
    services = Services()
    monkeypatch.setattr(species_web, "make_accession_transport", lambda: services)
    fasta = ">NP_012429.1\nMKV\n>CAI6736627\nMKV\n"
    text = client.post(
        "/species-check/", {"sequences_text": fasta, "reconcile": "on", "accession_consent": "on"}
    ).content.decode()
    assert "NCBI lookup" in text and 'value="559292"' in text
    assert "Found 2 of 2." in text and 'name="accession_consent" checked' in text


def test_rows_needing_the_user_come_first_and_offers_apply_to_many(client: Client) -> None:
    fasta = ">Scer_IME2\nMKV\n>Scer_KSS1\nMKV\n>Scer_FUS3\nMKV\n>Homo_sapiens_X\nMKV\n>odd\nMKV\n"
    text = client.post("/species-check/", {
        "sequences_text": fasta, "reconcile": "on",
        "species_for": ["Scer_IME2"], "species_query": ["Saccharomyces paradoxus"],
        "species_taxid": [""], "species_shown": ["Saccharomyces cerevisiae\t4932"],
    }).content.decode()  # fmt: skip
    order = re.findall(r'name="species_for" value="([^"]+)"', text)
    assert order[0] == "odd"  # the only row without a species
    assert re.search(r'data-rows="Scer_KSS1&#10;Scer_FUS3"', text)
    assert "to 2 other sequences named Scer" in " ".join(text.split())


def test_download_and_upload_the_table(client: Client, fake_storage: FakeStorage) -> None:
    fasta = ">Scer_IME2\nMKV\n>Scer_KSS1\nMKV\n"
    tsv = client.post("/species-table/", {"sequences_text": fasta, "reconcile": "on"})
    assert tsv["Content-Disposition"] == 'attachment; filename="species.tsv"'
    lines = tsv.content.decode().splitlines()
    assert lines[0].startswith("sequence\tspecies\ttaxid") and len(lines) == 3
    # filled in and pasted back, with the user's own species tree: the job reconciles
    filled = (
        "sequence\tspecies\ttaxid\n"
        "Scer_IME2\tSaccharomyces cerevisiae\t\nScer_KSS1\tCandida albicans\t\n"
    )
    protein = "MKVLAAGIRSTW" * 4
    response = client.post("/jobs/new/", {
        "aligners": ["mafft"], "tree_builders": ["iqtree"], "modelfinder": "on",
        "run_trees": "on", "run_ancestors": "on", "seed": 1, "reconcile": "on",
        "sequences_text": f">Scer_IME2\n{protein}\n>Scer_KSS1\n{protein}A\n>Calb_X\n{protein}C\n"
        f">Hsap_X\n{protein}D\n",
        "species_paste": filled + "Calb_X\tCandida albicans\nHsap_X,Homo sapiens,9606\n",
        "species_tree_file": SimpleUploadedFile(
            "s.nwk", b"((Saccharomyces_cerevisiae,Candida_albicans),Homo_sapiens);"
        ),
    })  # fmt: skip
    assert response.status_code == 302
    job = Job.objects.get()
    assert job.params.get("reconcile") is True, job.params.get("reconcile_skipped")
    stored = fake_storage.objects[f"jobs/{job.id}/inputs/species.tsv"].decode()
    assert {r.sequence: (r.source, r.species) for r in read(stored)} == {
        "Scer_IME2": ("list", "Saccharomyces cerevisiae"),
        "Scer_KSS1": ("list", "Candida albicans"),
        "Calb_X": ("list", "Candida albicans"),
        "Hsap_X": ("list", "Homo sapiens"),
    }


def test_a_taxon_left_out_of_the_tree_is_offered_its_parent() -> None:
    from phylobot.species.opentree import Taxon

    class Taxonomy:
        asked: list[int] = []

        def parent(self, ott_id: int) -> Taxon:
            self.asked.append(ott_id)
            return Taxon(277145, "Rhodospirillales", "order", (), True, ())

    rows = [Row("WP_1", query="Reyranella", ott_id=4780359, status="not in tree"),
            Row("WP_2", query="Reyranella", ott_id=4780359, status="not in tree"),
            Row("x", query="Homo sapiens", ott_id=770315, status=OK)]  # fmt: skip
    taxonomy = Taxonomy()
    found = species_web.parents(rows, lambda: taxonomy)  # type: ignore[arg-type,return-value]
    assert found == {"WP_1": ("Rhodospirillales", "order"), "WP_2": ("Rhodospirillales", "order")}
    assert taxonomy.asked == [4780359]  # once per taxon
