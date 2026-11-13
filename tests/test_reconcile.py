"""Reconciliation events on our trees (phylobot.reconcile) and the pipeline's
reconciliation steps (phylobot.pipeline.species_tree, reconcile_export).

GeneRax's real output is in tests/fixtures/generax/ (see its README)."""

import json
from pathlib import Path

import pytest

from phylobot import reconcile
from phylobot.pipeline import reconcile_export, species_tree
from phylobot.reconcile import ReconcileError
from phylobot.trees import load_tree

FIXTURES = Path(__file__).parent / "fixtures" / "generax"
XML = (FIXTURES / "discordant.xml").read_text()
GENE_TREE = (FIXTURES / "discordant_gene.nwk").read_text()


def by_node(rows: list[reconcile.Row]) -> dict[str, reconcile.Row]:
    return {row.node: row for row in rows}


def test_events_and_losses_land_on_our_nodes() -> None:
    rows, ends, rooted = reconcile.annotate(GENE_TREE, reconcile.parse(XML))
    nodes = by_node(rows)
    assert [r.node for r in rows if r.event == "D"] == ["NodeR", "NodeD", "NodeC"]
    assert nodes["NodeA"].event == nodes["NodeB"].event == "S"
    expected = {"duplications": 3, "transfers": 0, "speciations": 2, "losses": 12}
    assert reconcile.totals(rows) == expected
    # K. lactis's copy was lost in S. cerevisiae and L. waltii on the way down.
    lost = nodes["Kluyveromyces.lactis.IME2"].losses
    assert lost == ("Saccharomyces_cerevisiae", "Lachancea_waltii")
    assert nodes["Yarrowia.lipolytica.IME2"].species == "Yarrowia_lipolytica"
    # GeneRax roots the tree between Yarrowia and the rest.
    assert ends == ["NodeC", "Yarrowia.lipolytica.IME2"]
    top = [c for c in load_tree(rooted).seed_node.child_nodes()]
    assert sorted(len(list(c.leaf_iter())) for c in top) == [1, 6]


def test_rows_follow_our_tree_however_it_is_rooted() -> None:
    rerooted = (
        "(Yarrowia.lipolytica.IME2,Candida.lusitaniae.IME2,(Saccharomyces.cerevisiae.IME2,"
        "((Kluyveromyces.lactis.IME2,Candida.albicans.IME2)NodeA,"
        "(Kluyveromyces.waltii.IME2,Candida.parapsilosis.IME2)NodeB)NodeR)NodeD)NodeC;"
    )
    first, _, _ = reconcile.annotate(GENE_TREE, reconcile.parse(XML))
    second, _, _ = reconcile.annotate(rerooted, reconcile.parse(XML))
    assert {(r.node, r.event) for r in first} == {(r.node, r.event) for r in second}


def test_a_different_topology_is_refused() -> None:
    other = (
        "((Kluyveromyces.lactis.IME2,Kluyveromyces.waltii.IME2)NodeA,"
        "(Candida.albicans.IME2,Candida.parapsilosis.IME2)NodeB,"
        "(Saccharomyces.cerevisiae.IME2,(Candida.lusitaniae.IME2,Yarrowia.lipolytica.IME2)NodeC)NodeD)NodeR;"
    )
    with pytest.raises(ReconcileError, match="not in the reconciled tree"):
        reconcile.annotate(other, reconcile.parse(XML))


@pytest.mark.parametrize(
    ("xml", "message"),
    [
        ("<recPhylo>", "not valid XML"),
        ("<recPhylo><spTree/></recPhylo>", "no gene tree"),
        (XML.replace("<duplication", "<mystery").replace("</duplication", "</mystery"), "unknown"),
    ],
)
def test_malformed_reconciliations_refused(xml: str, message: str) -> None:
    with pytest.raises(ReconcileError, match=message):
        reconcile.parse(xml)


def write_species_inputs(tmp: Path) -> None:
    (tmp / "species.safe.tsv").write_text(
        "sequence\tspecies\n" + "".join(f"s{i:04d}\tsp00{i}\n" for i in range(1, 5))
    )
    (tmp / "species_tree.safe.nwk").write_text("(sp001,sp002,sp003,sp004)spn001;\n")
    (tmp / "species_names.tsv").write_text(
        "id\tname\tkind\nsp001\tYarrowia lipolytica\tspecies\nsp002\tCandida albicans\tspecies\n"
        "sp003\tCandida parapsilosis\tspecies\nsp004\tSaccharomyces cerevisiae\tspecies\n"
        "spn001\tSaccharomycetales\ttaxon\n"
    )
    (tmp / "outgroup.txt").write_text("s0001\n")
    (tmp / "gene.nwk").write_text("(s0001,((s0002,s0003),s0004));\n")


def test_species_tree_step_resolves_and_names_nodes(tmp_path: Path) -> None:
    write_species_inputs(tmp_path)
    out = tmp_path / "out"
    species_tree.write(
        tmp_path / "species.safe.tsv", tmp_path / "species_tree.safe.nwk",
        tmp_path / "species_names.tsv", tmp_path / "outgroup.txt", [tmp_path / "gene.nwk"],
        out, seed=1, method="guided",
    )  # fmt: skip
    assert (out / "species_tree.resolved.nwk").read_text() == (
        "(sp001,(sp004,(sp002,sp003)spr002)spr001)spn001;\n"
    )
    expected = "sp001:s0001\nsp002:s0002\nsp003:s0003\nsp004:s0004\n"
    assert (out / "mapping.link").read_text() == expected
    nodes = [line.split("\t") for line in (out / "species_nodes.tsv").read_text().splitlines()]
    assert nodes[0] == ["id", "name", "kind", "resolved_by", "votes"]
    assert ["spr002", "ancestor of Candida albicans and Candida parapsilosis", "resolved",
            "gene trees", "1"] in nodes  # fmt: skip
    assert (out / "species_tree.nwk").read_text() == (
        "(Yarrowia_lipolytica,(Saccharomyces_cerevisiae,(Candida_albicans,Candida_parapsilosis)))"
        "Saccharomycetales;\n"
    )
    report = json.loads((out / "polytomies.json").read_text())
    assert report["polytomies"][0]["label"] == "Saccharomycetales"
    assert report["root"] == "gene trees rooted on the outgroup"


def test_reconcile_export_names_everything(tmp_path: Path) -> None:
    names = {"Kluyveromyces.lactis.IME2": "K. lactis IME2"}
    (tmp_path / "names.tsv").write_text(
        "id\tname\tdescription\n" + "".join(f"{k}\t{v}\t\n" for k, v in names.items())
    )
    (tmp_path / "species_nodes.tsv").write_text(
        "id\tname\tkind\tresolved_by\tvotes\nspn002\tSaccharomycetales core\ttaxon\t\t\n"
        "Saccharomyces_cerevisiae\tSaccharomyces cerevisiae\tspecies\t\t\n"
    )
    (tmp_path / "model.txt").write_text("LG+G\n")
    out = tmp_path / "out"
    reconcile_export.export(
        tmp_path / "names.tsv", tmp_path / "species_nodes.tsv", FIXTURES / "discordant_gene.nwk",
        FIXTURES / "discordant.xml", FIXTURES / "stats.txt", tmp_path / "model.txt", out,
    )  # fmt: skip
    lines = (out / "reconciliation.tsv").read_text().splitlines()
    assert lines[0] == "node\tevent\tspecies\tlosses\tlost_in"
    assert "NodeR\tduplication\tSaccharomycetales core\t0\t" in lines
    lost = "Saccharomyces cerevisiae; Lachancea_waltii"  # named where species_nodes.tsv names it
    assert f"K. lactis IME2\tleaf\tKluyveromyces_lactis\t2\t{lost}" in lines
    summary = json.loads((out / "reconciliation.json").read_text())
    assert summary["duplications"] == 3 and summary["losses"] == 12
    assert summary["root_branch"] == ["NodeC", "Yarrowia.lipolytica.IME2"]
    assert summary["rates"]["duplication"] > 0.1
    assert summary["generax_model"] == "LG+G"
    rooted = load_tree((out / "reconciliation.nwk").read_text())
    assert [len(list(c.leaf_iter())) for c in rooted.seed_node.child_nodes()] == [1, 6]


def test_parse_stats() -> None:
    stats = reconcile_export.parse_stats("-9497.26 -2.6e-06\nReconciliation rates = 1e-07 0.5 \n")
    assert stats == {
        "sequence_log_likelihood": -9497.26,
        "reconciliation_log_likelihood": -2.6e-06,
        "rates": {"duplication": 1e-07, "loss": 0.5},
    }


def test_marks_for_the_viewer(tmp_path: Path) -> None:
    (tmp_path / "names.tsv").write_text("id\tname\tdescription\n")
    (tmp_path / "species.tsv").write_text("id\tname\tkind\tresolved_by\tvotes\n")
    (tmp_path / "model.txt").write_text("LG+G\n")
    out = tmp_path / "out"
    reconcile_export.export(
        tmp_path / "names.tsv", tmp_path / "species.tsv", FIXTURES / "discordant_gene.nwk",
        FIXTURES / "discordant.xml", FIXTURES / "stats.txt", tmp_path / "model.txt", out,
    )  # fmt: skip
    events, losses = reconcile.marks(
        (out / "reconciliation.tsv").read_text(), (out / "reconciliation.nwk").read_text()
    )
    assert set(events) == {"NodeR", "NodeD", "NodeC"}
    assert events["NodeR"][:2] == ("duplication", "spn002")
    assert events["NodeR"][2] == reconcile.Overlap(0, 4)  # no species with both: weak
    assert events["NodeR"][2].weak
    assert sum(len(v) for v in losses.values()) == 12
    assert losses[frozenset({"Kluyveromyces.lactis.IME2"})] == [
        "Saccharomyces_cerevisiae", "Lachancea_waltii",
    ]  # fmt: skip


def test_job_description_mentions_reconciliation() -> None:
    from phylobot.forms import describe

    params = {"aligners": ["mafft"], "tree_builders": ["iqtree"], "modelfinder": True}
    assert not any("species tree" in step for step in describe(params, 0, 0))
    steps = describe(params | {"reconcile": True, "polytomies": "guided"}, 0, 0)
    assert any("reconcile every tree with the species tree" in step for step in steps)
    assert not any("GeneRax for each" in step for step in steps)
    steps = describe(params | {"reconcile": True, "tree_builders": ["iqtree", "generax"]}, 0, 0)
    assert any("search a reconciled tree with GeneRax" in step for step in steps)


HEADER = "node\tevent\tspecies\tlosses\tlost_in\n"


def tsv(*rows: tuple[str, str]) -> str:
    return HEADER + "".join(f"{node}\t{event}\tsp\t0\t\n" for node, event in rows)


# The reconciled tree: two paralog groups (1 and 2) of species A and B, split at the root.
RECONCILED = (
    tsv(("N1", "speciation"), ("N2", "speciation"), (reconcile.ROOT, "duplication")),
    "((A1,B1)N1,(A2,B2)N2);",
)


def test_an_ancestor_mixing_paralogs_is_flagged() -> None:
    # The ML tree groups by species: its M1 averages A's two paralogs.
    own = (tsv(("M1", "duplication"), ("M2", "duplication")), "((A1,A2)M1,(B1,B2)M2);")
    found = reconcile.paralogy("M1", own, RECONCILED)
    assert found == reconcile.Paralogy(
        duplication=True, mixes=reconcile.ROOT, groups=(("A1",), ("A2",)), group_nodes=("", "")
    )
    bigger = ("((A1,B1,C1)N1,(A2,B2)N2);", "((A1,C1,A2)M1,(B1,B2)M2);")
    found = reconcile.paralogy("M1", (own[0], bigger[1]), (RECONCILED[0], bigger[0]))
    assert found is not None and found.groups == (("A1", "C1"), ("A2",))
    assert found.group_nodes == ("N1", "")


def test_an_ancestor_in_the_reconciled_tree_is_matched() -> None:
    own = (tsv(("M1", "speciation"), ("M2", "speciation")), "((A1,B1)M1,(A2,B2)M2);")
    assert reconcile.paralogy("M2", own, RECONCILED) == reconcile.Paralogy(
        False, "N2", "speciation"
    )
    assert reconcile.paralogy("M2", own, None) == reconcile.Paralogy(False)
    assert reconcile.paralogy("Nowhere", own, RECONCILED) is None


def test_a_clade_absent_under_a_speciation_is_not_paralogy() -> None:
    reconciled = (tsv(("N1", "speciation"), ("N2", "speciation")), "(((A,B)N1,C)N2,D);")
    own = (tsv(("M1", "speciation"), ("M2", "speciation")), "(((A,C)M1,B)M2,D);")
    assert reconcile.paralogy("M1", own, reconciled) == reconcile.Paralogy(False)


def test_the_root_row_carries_generax_root(tmp_path: Path) -> None:
    rec = reconcile.parse(XML)
    row = reconcile.root_row(rec)
    assert row.node == reconcile.ROOT and row.event == "S"
    (tmp_path / "names.tsv").write_text("id\tname\tdescription\n")
    (tmp_path / "species.tsv").write_text("id\tname\tkind\tresolved_by\tvotes\n")
    (tmp_path / "model.txt").write_text("LG+G\n")
    out = tmp_path / "out"
    reconcile_export.export(
        tmp_path / "names.tsv", tmp_path / "species.tsv", FIXTURES / "discordant_gene.nwk",
        FIXTURES / "discordant.xml", FIXTURES / "stats.txt", tmp_path / "model.txt", out,
    )  # fmt: skip
    last = (out / "reconciliation.tsv").read_text().splitlines()[-1]
    assert last.startswith("(root)\tspeciation")
    summary = json.loads((out / "reconciliation.json").read_text())
    assert summary["root_event"] == "speciation"


def test_species_overlap_marks_weak_duplications() -> None:
    header = "node\tevent\tspecies\tlosses\tlost_in\n"
    leaves = {"a1": "A", "a2": "A", "b1": "B", "c1": "C", "c2": "C", "d1": "D"}
    rows = "".join(f"{g}\tleaf\t{sp}\t0\t\n" for g, sp in leaves.items())
    # a duplication within species A (copies a1, a2): real, overlap 1 of 1
    recent = "((a1,a2)N1,b1);"
    found = reconcile.overlaps(header + "N1\tduplication\tA\t0\t\n" + rows, recent)
    assert found["N1"] == reconcile.Overlap(1, 1) and not found["N1"].weak
    # copies (a1, b1) and (c1, d1): no species in both, a typical tree error
    wrong = "((a1,b1)N2,(c1,d1)N3)N4;"
    tsv = header + "(root)\tduplication\tS\t0\t\n" + rows
    assert reconcile.overlaps(tsv, wrong)[reconcile.ROOT] == reconcile.Overlap(0, 4)
    assert reconcile.weak_duplications(tsv, wrong) == 1
