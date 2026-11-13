"""Species for reconciliation (phylobot.species, reconciliation R1).

The Open Tree of Life is replayed from answers recorded on 2026-10-16
(tests/fixtures/opentree/, synthesis opentree16.1); no test touches the network.
"""

import json
import sys
from pathlib import Path
from typing import Any

import pytest

from phylobot.pipeline import validate
from phylobot.species import checks, headers, polytomies, resolve, speciestree, table
from phylobot.species.opentree import Answer, Cache, Client, OpenTreeError
from phylobot.species.polytomies import PolytomyError
from phylobot.species.speciestree import Node, SpeciesTreeError
from phylobot.species.table import AMBIGUOUS, NOT_IN_TREE, OK, UNRESOLVED, Row, TableError
from tests.fixtures.opentree import Replay
from tests.fixtures.opentree.record import AWKWARD, example_records

FIXTURES = Path(__file__).parent / "fixtures" / "opentree"


def resolved(name: str, rows: list[Row]) -> resolve.Resolution:
    return resolve.with_open_tree(rows, Client(Replay(name)))


def by_sequence(result: resolve.Resolution) -> dict[str, Row]:
    return {row.sequence: row for row in result.rows}


def clades(node: Node) -> set[frozenset[str]]:
    return {frozenset(n.leaves()) for n in node.preorder() if not n.is_leaf()}


# Headers


@pytest.mark.parametrize(
    ("name", "description", "expected"),
    [
        ("Saccharomyces.cerevisiae.IME2", "", ("Saccharomyces cerevisiae", None, "name")),
        ("Yarrowia_lipolytica_Q6C4M5", "", ("Yarrowia lipolytica", None, "name")),
        ("Haliotis_discus_hannai", "", ("Haliotis discus hannai", None, "name")),
        ("Homo_sapiens|P12345", "", ("Homo sapiens", None, "name")),
        ("Candida_sp_ABC", "", ("Candida sp", None, "name")),
        (
            "sp|P32581|IME2_YEAST",
            "Meiosis-inducing kinase OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) "
            "OX=559292 GN=IME2 PE=1 SV=2",
            ("Saccharomyces cerevisiae", 559292, "uniprot"),
        ),
        (
            "UniRef90_P32581",
            "Meiosis-inducing kinase n=12 Tax=Saccharomyces TaxID=4930 RepID=IME2_YEAST",
            ("Saccharomyces", 4930, "uniref"),
        ),
        (
            "XP_500001.1",
            "kinase [Yarrowia lipolytica CLIB122]",
            ("Yarrowia lipolytica CLIB122", None, "ncbi"),
        ),
    ],
)
def test_headers_name_their_species(name: str, description: str, expected: tuple[Any, ...]) -> None:
    found = headers.detect(name, description)
    assert found is not None
    assert (found.name, found.taxid, found.source) == expected


@pytest.mark.parametrize("name", ["seq1", "IME2", "s0001", "ABC_def", "lowercase_genus"])
def test_headers_without_a_species(name: str) -> None:
    assert headers.detect(name, "a protein") is None


# species.tsv


def test_table_round_trip() -> None:
    rows = [
        Row(
            "a",
            query="Homo sapiens",
            taxid=9606,
            ott_id=770315,
            species="Homo sapiens",
            source="uniprot",
            status=OK,
        ),  # fmt: skip
        Row("b", note="no species found in the header"),
    ]
    assert table.read(table.write(rows)) == rows


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("", "empty"),
        ("query\nHomo\n", "'sequence' column"),
        ("sequence\tcolour\na\tred\n", "unknown columns"),
        ("sequence\tquery\na\t\na\t\n", "listed twice"),
        ("sequence\ttaxid\na\t96O6\n", "whole number"),
        ("sequence\tquery\na\t" + "x" * 300 + "\n", "longer than"),
        ("sequence\tquery\na\tb\tc\n", "more cells"),
        ("sequence\tquery\n\tx\n", "no sequence name"),
    ],
)
def test_table_refuses_malformed_input(text: str, message: str) -> None:
    with pytest.raises(TableError, match=message):
        table.read(text)


# Open Tree: the examples


def test_ime2_resolves_with_current_names() -> None:
    result = resolved("ime2", resolve.detect(example_records("ime2")))
    rows = by_sequence(result)
    assert result.complete
    assert rows["Kluyveromyces.waltii.IME2"].species == "Lachancea waltii"
    renamed = rows["Kluyveromyces.waltii.IME2"]
    assert renamed.note == "'Kluyveromyces waltii' is now 'Lachancea waltii'"
    assert rows["Candida.lusitaniae.IME2"].species == "Clavispora lusitaniae"
    assert result.tree is not None
    assert sorted(result.tree.leaves()) == sorted(row.species for row in result.rows)
    # Open Tree leaves the yeasts' families as a four-way polytomy at the root, and
    # Saccharomycetaceae as a three-way one.
    assert [len(p.children) for p in result.tree.polytomies()] == [4, 3]
    assert result.provenance["synth_id"] == "opentree16.1"
    assert len(result.provenance["citations"]) == 4


def test_lysin_synonyms_become_one_species() -> None:
    result = resolved("lysin", resolve.detect(example_records("lysin")))
    rows = by_sequence(result)
    assert result.complete
    synonym, accepted = rows["Haliotis_sieboldii"], rows["Haliotis_gigantea"]
    assert synonym.species == accepted.species == "Haliotis gigantea"
    assert rows["Haliotis_discus_hannai"].species == "Haliotis discus hannai"  # a subspecies
    assert result.tree is not None
    assert len(result.tree.leaves()) == len(result.rows) - 1
    # The abalones are almost one star: only the two H. tuberculata subspecies are grouped.
    assert clades(result.tree) == {
        frozenset(result.tree.leaves()),
        frozenset({"Haliotis tuberculata coccinea", "Haliotis tuberculata tuberculata"}),
    }


def test_awkward_names_say_what_the_user_must_do() -> None:
    result = resolved("awkward", AWKWARD)
    rows = by_sequence(result)
    assert not result.complete
    status = {name: row.status for name, row in rows.items()}
    assert status == {
        "misspelt": UNRESOLVED,
        "hybrid": OK,
        "genus": AMBIGUOUS,
        "inside_genus": OK,
        "nonsense": UNRESOLVED,
        "by_taxid": OK,
        "strain_species": OK,
        "pinned": OK,
        "pinned_pruned": NOT_IN_TREE,
        "gene_word": OK,
        "broken": OK,
        "bird": OK,
        "taxid_only": OK,
    }
    assert "did you mean Saccharomyces cerevisiae" in rows["misspelt"].note
    assert "Candida albicans" in rows["genus"].note
    # A strain (from a UniProt taxid) joins its species when the species is in the job.
    assert rows["by_taxid"].ott_id == 170091
    assert rows["by_taxid"].species == rows["strain_species"].species == "Saccharomyces cerevisiae"
    assert rows["pinned"].species == "Clavispora lusitaniae"
    assert rows["gene_word"].species == "Candida parapsilosis"
    assert rows["gene_word"].note == "matched as 'Candida parapsilosis'"
    assert rows["taxid_only"].species == "Yarrowia lipolytica"
    assert "stand-in" in rows["broken"].note
    assert result.provenance["broken"] == {"ott632176": "ott985906"}
    assert result.tree is not None
    placed = {row.species for row in result.rows if row.species}
    assert sorted(result.tree.leaves()) == sorted(placed)


def test_pruned_ids_are_dropped_and_the_tree_fetched_again() -> None:
    newick = (
        "((Saccharomyces_cerevisiae_ott356221)Saccharomyces_ott908546,"
        "Candida_albicans_ott154960)x_ott1;"
    )
    calls: list[list[int]] = []

    def transport(endpoint: str, payload: dict[str, Any]) -> Answer:
        if endpoint == "tree_of_life/about":
            return Answer(200, {"synth_id": "opentree16.1", "taxonomy_version": "3.7draft3"})
        calls.append(payload["ott_ids"])
        if 6376729 in payload["ott_ids"]:  # Open Tree's real answer, 2026-10-16
            return Answer(400, {
                "message": "[/v3/tree_of_life/induced_subtree] Error: node_id 'ott6376729' "
                "was not found!",
                "unknown": {"ott6376729": "pruned_ott_id"},
            })  # fmt: skip
        return Answer(200, {"newick": newick, "broken": {}, "supporting_studies": []})

    induced = Client(transport).induced_subtree([356221, 6376729, 154960])
    assert calls == [[154960, 356221, 6376729], [154960, 356221]]
    assert induced.pruned == (6376729,)


def test_open_tree_errors_are_reported() -> None:
    def transport(endpoint: str, payload: dict[str, Any]) -> Answer:
        if endpoint == "tree_of_life/about":
            return Answer(200, {"synth_id": "opentree16.1", "taxonomy_version": "3.7draft3"})
        return Answer(500, {"message": "internal error"})

    with pytest.raises(OpenTreeError, match="answered 500: internal error"):
        Client(transport).match_names(["Homo sapiens"])


def test_cache_answers_repeat_requests(tmp_path: Path) -> None:
    rows = resolve.detect(example_records("ime2"))
    first = Replay("ime2")
    resolve.with_open_tree(rows, Client(first, Cache(tmp_path)))
    again = Replay("ime2")
    resolve.with_open_tree(rows, Client(again, Cache(tmp_path)))
    assert len(first.calls) == 3
    assert again.calls == ["tree_of_life/about"]  # the release is always asked


# Species trees


def test_species_tree_newick_round_trips_awkward_labels() -> None:
    tree = Node("Saccharomycetales", [
        Node("Saccharomyces cerevisiae"),
        Node("Candida (genus in Candida/Lodderomyces clade)", [Node("Candida albicans"),
                                                               Node("O'Brien's yeast")]),
    ])  # fmt: skip
    parsed = speciestree.parse(tree.newick())
    assert parsed.newick() == tree.newick()
    expected = ["Candida albicans", "O'Brien's yeast", "Saccharomyces cerevisiae"]
    assert sorted(parsed.leaves()) == expected


def test_upload_is_pruned_to_the_job_species() -> None:
    newick = (
        "((Homo_sapiens:1,Pan_troglodytes:1):1,(Mus_musculus:1,Rattus_norvegicus:1):1,"
        "Gallus_gallus:2);"
    )
    tree = speciestree.from_upload(newick, {"Homo sapiens", "Mus  musculus", "Gallus gallus"})
    assert tree.newick() == "(Homo_sapiens,Mus_musculus,Gallus_gallus);"


@pytest.mark.parametrize(
    ("newick", "message"),
    [
        ("((A,B),(A,C));", "appears twice"),
        ("((A,B),(C,));", "empty element"),
        ("((A,X),Y);", "fewer than two"),
    ],
)
def test_upload_refuses_bad_trees(newick: str, message: str) -> None:
    with pytest.raises(SpeciesTreeError, match=message):
        speciestree.from_upload(newick, {"A", "B", "C"})


def test_resolving_against_an_upload_contacts_nothing() -> None:
    rows = [Row("g1", query="Homo sapiens"), Row("g2", query="Homo_sapiens"),
            Row("g3", query="Mus musculus"), Row("g4", query="Danio rerio"), Row("g5")]  # fmt: skip
    result = resolve.with_upload(rows, "((Homo_sapiens,Pan_troglodytes),Mus_musculus);")
    status = {row.sequence: (row.status, row.species) for row in result.rows}
    assert status == {
        "g1": (OK, "Homo sapiens"),
        "g2": (OK, "Homo sapiens"),
        "g3": (OK, "Mus musculus"),
        "g4": (NOT_IN_TREE, ""),
        "g5": (UNRESOLVED, ""),
    }
    assert result.tree is not None and result.tree.newick() == "(Homo_sapiens,Mus_musculus);"


# Polytomies

IME2_TREE = (
    "(Clavispora_lusitaniae,(Saccharomyces_cerevisiae,Kluyveromyces_lactis,Lachancea_waltii)"
    "Saccharomycetaceae,(Candida_albicans,Candida_parapsilosis)Candida,Yarrowia_lipolytica)"
    "Saccharomycetales;"
)
IME2_SPECIES = {
    "kwal": "Lachancea waltii", "scer": "Saccharomyces cerevisiae",
    "clus": "Clavispora lusitaniae", "ylip": "Yarrowia lipolytica",
    "cpar": "Candida parapsilosis", "calb": "Candida albicans", "klac": "Kluyveromyces lactis",
}  # fmt: skip
IME2_GENE_TREE = (
    "(ylip:0.5,((calb:0.1,cpar:0.1):0.2,clus:0.3):0.1,(scer:0.2,(klac:0.1,kwal:0.1):0.1):0.2);"
)
CTG = frozenset({"Candida albicans", "Candida parapsilosis", "Clavispora lusitaniae"})
KLUYVEROMYCES = frozenset({"Kluyveromyces lactis", "Lachancea waltii"})


def species_tree(newick: str = IME2_TREE) -> Node:
    return speciestree.parse(newick)


def is_binary(tree: Node) -> bool:
    return all(len(n.children) in (0, 2) for n in tree.preorder())


def test_guided_by_outgroup_rooted_gene_trees() -> None:
    tree, report = polytomies.resolve(
        species_tree(), [IME2_GENE_TREE], IME2_SPECIES, seed=1, outgroup={"ylip"}
    )
    assert is_binary(tree)
    assert sorted(tree.leaves()) == sorted(IME2_SPECIES.values())
    root = [frozenset(c.leaves()) for c in tree.children]
    assert frozenset({"Yarrowia lipolytica"}) in root  # the outgroup's side of the root
    assert CTG in clades(tree) and KLUYVEROMYCES in clades(tree)
    assert report["root"] == "gene trees rooted on the outgroup"
    added = [a for p in report["polytomies"] for a in p["added"]]
    assert all(a["by"] == "gene trees" and a["votes"] == 1 for a in added)
    assert [p["label"] for p in report["polytomies"]] == ["Saccharomycetales", "Saccharomycetaceae"]


def test_without_an_outgroup_the_root_is_left_to_chance_and_said_so() -> None:
    tree, report = polytomies.resolve(species_tree(), [IME2_GENE_TREE], IME2_SPECIES, seed=1)
    assert is_binary(tree)
    assert KLUYVEROMYCES in clades(tree)  # below the root, splits are enough
    assert report["root"] == "arbitrary (unrooted gene trees cannot place it)"


def test_majority_of_gene_trees_wins_and_votes_are_counted() -> None:
    other = "(ylip,(calb,cpar,clus),((scer,klac),kwal));"  # S. cerevisiae with K. lactis
    tree, report = polytomies.resolve(
        species_tree(), [IME2_GENE_TREE, IME2_GENE_TREE, other], IME2_SPECIES, seed=3,
        outgroup={"ylip"},
    )  # fmt: skip
    assert KLUYVEROMYCES in clades(tree)
    family = next(p for p in report["polytomies"] if p["label"] == "Saccharomycetaceae")
    assert family["informative"] == 3
    assert family["added"] == [{"species": sorted(KLUYVEROMYCES), "by": "gene trees", "votes": 2}]


def test_silent_gene_trees_leave_an_arbitrary_resolution_marked_as_such() -> None:
    star = "(ylip,calb,cpar,clus,scer,klac,kwal);"
    tree, report = polytomies.resolve(species_tree(), [star], IME2_SPECIES, seed=5)
    assert is_binary(tree)
    assert {a["by"] for p in report["polytomies"] for a in p["added"]} == {"arbitrary"}
    assert all(p["informative"] == 0 for p in report["polytomies"])


def test_paralogs_vote_by_species() -> None:
    # Two copies of each Saccharomycetaceae gene after a duplication; both copies group
    # K. lactis with L. waltii.
    species = {**IME2_SPECIES, "scer2": "Saccharomyces cerevisiae",
               "klac2": "Kluyveromyces lactis", "kwal2": "Lachancea waltii"}  # fmt: skip
    gene_tree = "(ylip,(calb,cpar,clus),((scer,(klac,kwal)),(scer2,(klac2,kwal2))));"
    tree, _ = polytomies.resolve(species_tree(), [gene_tree], species, seed=1, outgroup={"ylip"})
    assert KLUYVEROMYCES in clades(tree)


def test_arbitrary_method_follows_the_seed_and_ignores_gene_trees() -> None:
    def run(seed: int) -> str:
        tree, report = polytomies.resolve(
            species_tree(), [IME2_GENE_TREE], IME2_SPECIES, seed=seed, method="arbitrary"
        )
        assert report["gene_trees"] == 0 and report["root"] == "arbitrary"
        return tree.newick()

    assert run(7) == run(7)
    assert len({run(seed) for seed in range(20)}) > 1


def test_lysin_star_resolved_from_gene_trees() -> None:
    names = ["A", "B", "C", "D", "E", "F"]
    star = speciestree.parse("(" + ",".join(f"Haliotis_{n}" for n in names) + ")Haliotis;")
    species = {n.lower(): f"Haliotis {n}" for n in names}
    gene = "(((a,b),c),((d,e),f));"
    tree, report = polytomies.resolve(star, [gene], species, seed=2)
    assert frozenset({"Haliotis A", "Haliotis B"}) in clades(tree)
    assert frozenset({"Haliotis D", "Haliotis E"}) in clades(tree)
    assert is_binary(tree)
    assert report["root"].startswith("arbitrary")


@pytest.mark.parametrize(
    ("species", "method", "message"),
    [
        ({**IME2_SPECIES, "kwal": "Not in tree"}, "guided", "not in the species tree"),
        ({k: v for k, v in IME2_SPECIES.items() if k != "kwal"}, "guided", "without a species"),
        (IME2_SPECIES, "random", "unknown polytomy method"),
    ],
)
def test_polytomy_inputs_checked(species: dict[str, str], method: str, message: str) -> None:
    with pytest.raises(PolytomyError, match=message):
        polytomies.resolve(species_tree(), [IME2_GENE_TREE], species, seed=1, method=method)


# VALIDATE

SAFE = {name: f"s{i:04d}" for i, name in enumerate(IME2_SPECIES, 1)}
WITHOUT_KWAL = {k: v for k, v in IME2_SPECIES.items() if k != "kwal"}


def species_table(species: dict[str, str]) -> str:
    rows = [Row(seq, query=name, species=name, status=OK) for seq, name in species.items()]
    return table.write(rows)


def test_checks_rename_species_and_prune_the_tree() -> None:
    bigger = IME2_TREE.replace("Yarrowia_lipolytica", "(Yarrowia_lipolytica,Homo_sapiens)")
    prepared = checks.prepare(species_table(IME2_SPECIES), bigger, SAFE)
    rows = [line.split("\t") for line in prepared.names.splitlines()[1:]]
    names = {i: n for i, n, kind in rows if kind == "species"}
    assert "Homo sapiens" not in names.values() and len(names) == 7
    tree = speciestree.parse(prepared.tree)
    assert sorted(tree.leaves()) == sorted(names)
    internal = {i: (n, kind) for i, n, kind in rows if kind != "species"}
    assert [n.label for n in tree.preorder() if not n.is_leaf()] == sorted(internal)
    assert internal["spn001"] == ("Saccharomycetales", "taxon")
    assert len(internal) == 3  # (Yarrowia, Homo) went with Homo sapiens
    mapping = {r.sequence: r.species for r in table.read(prepared.mapping)}
    assert mapping[SAFE["kwal"]] == next(i for i, n in names.items() if n == "Lachancea waltii")
    assert prepared.report == {"species": 7, "pruned_species": 1, "polytomies": 2,
                               "species_with_several_sequences": 0}  # fmt: skip


@pytest.mark.parametrize(
    ("species", "tree", "message"),
    [
        (WITHOUT_KWAL, IME2_TREE, "no species for 'kwal'"),
        ({**IME2_SPECIES, "extra": "Homo sapiens"}, IME2_TREE, "not among the sequences"),
        ({**IME2_SPECIES, "kwal": ""}, IME2_TREE, "without a species"),
        ({**IME2_SPECIES, "kwal": "Homo sapiens"}, IME2_TREE, "not in the species tree"),
        (dict.fromkeys(IME2_SPECIES, "Lachancea waltii"), IME2_TREE, "at least two species"),
        (IME2_SPECIES, "((a,b);", "not valid Newick"),
    ],
)
def test_checks_refuse(species: dict[str, str], tree: str, message: str) -> None:
    with pytest.raises(checks.SpeciesError, match=message):
        checks.prepare(species_table(species), tree, SAFE)


def test_validate_writes_species_files_only_when_asked(tmp_path: Path) -> None:
    fasta = "".join(f">{name}\nMKTAYIAKQRQISFVKSHFSRQ\n" for name in IME2_SPECIES)
    (tmp_path / "input.fasta").write_text(fasta)
    manifest = {"sequences": True, "max_sequences": 50, "max_length": 500, "max_provided": 2}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    validate.validate(tmp_path / "manifest.json", tmp_path)
    assert not (tmp_path / "species.safe.tsv").exists()
    assert "species" not in json.loads((tmp_path / "validation.json").read_text())

    (tmp_path / "species.tsv").write_text(species_table(IME2_SPECIES))
    (tmp_path / "species_tree.nwk").write_text(IME2_TREE)
    (tmp_path / "manifest.json").write_text(json.dumps({**manifest, "species": True}))
    validate.validate(tmp_path / "manifest.json", tmp_path)
    assert (tmp_path / "species.safe.tsv").read_text().startswith("sequence\tspecies\ns0001\tsp")
    assert json.loads((tmp_path / "validation.json").read_text())["species"]["species"] == 7

    (tmp_path / "species_tree.nwk").write_text(IME2_TREE.replace("Yarrowia_lipolytica", "Homo"))
    with pytest.raises(validate.InputError, match="Yarrowia lipolytica"):
        validate.validate(tmp_path / "manifest.json", tmp_path)


# The command line


def test_cli_with_an_uploaded_tree(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from phylobot.species import cli

    fasta = tmp_path / "in.fasta"
    names = ["Homo_sapiens_1", "Mus_musculus_2", "seq3"]
    fasta.write_text("".join(f">{n}\nMKTAYIAK\n" for n in names))
    upload = tmp_path / "species.nwk"
    upload.write_text("((Homo_sapiens,Pan_troglodytes),Mus_musculus);")
    out = tmp_path / "out"

    def run(*argv: str | Path) -> int:
        monkeypatch.setattr(sys, "argv", ["phylobot-species", *map(str, argv)])
        with pytest.raises(SystemExit) as stop:
            cli.main()
        return int(stop.value.code or 0)

    assert run("resolve", fasta, out, "--tree", upload) == 1  # seq3 has no species yet
    rows = table.read((out / "species.tsv").read_text())
    assert [r.status for r in rows] == [OK, OK, UNRESOLVED]
    edited = [r if r.sequence != "seq3" else Row("seq3", query="Homo sapiens") for r in rows]
    (out / "edited.tsv").write_text(table.write(edited))
    assert run("resolve", fasta, out, "--tree", upload, "--table", out / "edited.tsv") == 0
    assert (out / "species_tree.nwk").read_text() == "(Homo_sapiens,Mus_musculus);\n"

    gene = tmp_path / "gene.nwk"
    gene.write_text("((Homo_sapiens_1,seq3),Mus_musculus_2);")
    resolved_tree, report = tmp_path / "resolved.nwk", tmp_path / "report.json"
    species_files = (out / "species_tree.nwk", out / "species.tsv")
    assert run("polytomies", *species_files, resolved_tree, report, "--gene-trees", gene) == 0
    assert json.loads(report.read_text())["polytomies"] == []
