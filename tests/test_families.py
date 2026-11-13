"""Gene families: duplication and loss along a species tree (phylobot.simulate.families).

Pure Python; AliSim on a grown family runs in tests/golden/test_simulate.py.
"""

import json
import math
import random
from pathlib import Path
from typing import Any

import pytest

from phylobot.seqio import read_fasta, write_fasta
from phylobot.simulate import SpecError, parse
from phylobot.simulate import families as F
from phylobot.simulate.engine import collect, prepare

SPECIES = "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);"
IMAGE = "quay.io/biocontainers/iqtree:3.1.3--h8471819_0"


def spec(family: dict[str, Any], **simulation: Any) -> Any:
    return parse(
        {
            "tree": {},
            "root": {"length": 50},
            "model": {"substitution": "LG"},
            "family": family,
            "simulation": {"seeds": [1], **simulation},
        }
    )


def family(**values: Any) -> Any:
    family = spec({"duplication": 0, "loss": 0, **values}).family
    assert family is not None
    return family


HIDDEN_PARALOGY = {
    "duplication": 0,
    "loss": 0,
    "origin": 0.1,
    "duplications": [{"clade": ["A", "B", "C", "D"]}],
    "losses": [{"clade": ["A", "B"], "copy": "1b"}, {"clade": ["C", "D"], "copy": "1a"}],
}


def test_without_events_the_gene_tree_is_the_species_tree() -> None:
    grown = F.grow(family(), SPECIES, 1)
    assert F.gene_newick(grown) == "((A_1:0.1,B_1:0.1)Node2:0.2,(C_1:0.1,D_1:0.1)Node3:0.2)Node1;"
    assert F.species_newick(grown.species) == "((A:0.1,B:0.1)S2:0.2,(C:0.1,D:0.1)S3:0.2)S1;"
    rows = F.events_table(grown).splitlines()[1:]
    assert [r.split("\t")[:3] for r in rows] == [
        ["speciation", "Node1", "S1"],
        ["speciation", "Node2", "S2"],
        ["speciation", "Node3", "S3"],
    ]


def test_hidden_paralogy() -> None:
    """A duplication before the first speciation and reciprocal losses: the gene tree
    has the species tree's shape, but its root is a duplication."""
    grown = F.grow(family(**HIDDEN_PARALOGY), SPECIES, 1)
    assert F.gene_newick(grown) == (
        "((A_1:0.1,B_1:0.1)Node2:0.25,(C_1:0.1,D_1:0.1)Node3:0.25)Node1;"
    )
    assert F.events_table(grown) == (
        "event\tgene_node\tspecies_branch\tdepth\tforced\tcopies\n"
        "duplication\tNode1\tS1\t0.05\tyes\t-\n"
        "speciation\t-\tS1\t0.1\tno\t1a\n"
        "speciation\tNode2\tS2\t0.3\tno\t1a\n"
        "loss\t-\tS3\t0.2\tyes\t1a\n"
        "speciation\t-\tS1\t0.1\tno\t1b\n"
        "loss\t-\tS2\t0.2\tyes\t1b\n"
        "speciation\tNode3\tS3\t0.3\tno\t1b\n"
    )
    summary = F.summary(grown)
    assert summary["duplications_visible"] == 1 and summary["losses"] == 2


def test_the_same_seed_grows_the_same_family() -> None:
    rates = family(duplication=2, loss=1.5)
    assert F.gene_newick(F.grow(rates, SPECIES, 7)) == F.gene_newick(F.grow(rates, SPECIES, 7))
    assert F.gene_newick(F.grow(rates, SPECIES, 7)) != F.gene_newick(F.grow(rates, SPECIES, 8))


@pytest.mark.parametrize("seed", range(40))
def test_events_are_consistent_with_both_trees(seed: int) -> None:
    """A duplication's children stay on its species branch; a speciation's children
    start on the child branches of its species node; depths only increase."""
    grown = F.grow(family(duplication=1.5, loss=1.0, origin=0.05), SPECIES, seed)
    children = {s.name: [c.name for c in s.children] for s in F.preorder(grown.species)}
    for node in F.preorder(grown.origin):
        for child in node.children:
            assert child.depth >= node.depth
            if node.event == "duplication":
                assert child.branch == node.branch
            elif node.event == "speciation":
                assert child.branch in children[node.branch]
    leaves = [n for n in F.preorder(grown.root) if n.event == "leaf" and n.name]
    for leaf in leaves:  # every surviving gene ends at its species' tip
        assert not children[leaf.branch]
        assert math.isclose(leaf.depth, 0.05 + 0.3)
    named = [n for n in F.preorder(grown.origin) if n.name and n.event != "leaf"]
    assert {n.event for n in named} <= {"duplication", "speciation"}
    assert len(named) == len(leaves) - 1  # a binary gene tree


def leaves_per_species(origin: F.Gene) -> dict[str, int]:
    counts: dict[str, int] = {}
    for node in F.preorder(origin):
        if node.event == "leaf":
            counts[node.branch] = counts.get(node.branch, 0) + 1
    return counts


def test_copy_numbers_match_birth_death_theory() -> None:
    """Unconditioned linear birth-death from one copy over time T: the mean number of
    copies is exp((lambda - mu) T); with loss only, a copy survives with exp(-mu T)."""
    root = F.species_tree("(A:1,B:1);")
    rng = random.Random(2026)
    n = 3000
    rates = family(duplication=0.6, loss=0.3)
    mean = sum(leaves_per_species(F._attempt(rates, root, {}, rng)).get("A", 0) for _ in range(n))
    variance = (0.9 / 0.3) * math.exp(0.3) * (math.exp(0.3) - 1)
    assert abs(mean / n - math.exp(0.3)) < 4 * math.sqrt(variance / n)

    losses = family(loss=0.5)
    alive = sum("A" in leaves_per_species(F._attempt(losses, root, {}, rng)) for _ in range(n))
    p = math.exp(-0.5)
    assert abs(alive / n - p) < 4 * math.sqrt(p * (1 - p) / n)


def test_conditioning_on_survival() -> None:
    grown = F.grow(family(loss=2.0, min_genes=3), "(A:1,B:1,C:1,D:1);", 3)
    genes = F.summary(grown)["genes"]
    assert grown.attempts > 1 and isinstance(genes, int) and genes >= 3
    with pytest.raises(F.FamilyError, match="survived in 1000 attempts"):
        F.grow(family(loss=100), SPECIES, 1)


def test_a_runaway_family_is_refused() -> None:
    with pytest.raises(F.FamilyError, match="grew past"):
        F.grow(family(duplication=100), "(A:1,B:1);", 1)


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ({"duplications": [{"clade": ["A", "C"]}]}, "no species branch has exactly A, C"),
        ({"duplications": [{"clade": ["A", "B", "C", "D"]}]}, "the stem needs origin > 0"),
        (
            {"duplications": [{"clade": ["A", "B"]}], "losses": [{"clade": ["C"], "copy": "1a"}]},
            "on its branch or below it",
        ),
        (
            {
                "duplications": [{"clade": ["A", "B"], "at": 0.6}],
                "losses": [{"clade": ["A", "B"], "copy": "1a", "at": 0.4}],
            },
            "must come after duplication 1",
        ),
    ],
)
def test_chosen_events_must_fit_the_species_tree(values: dict[str, Any], message: str) -> None:
    with pytest.raises(F.FamilyError, match=message):
        F.grow(family(**values), SPECIES, 1)


def test_species_tree_names() -> None:
    with pytest.raises(F.FamilyError, match="taken by the tree's internal nodes"):
        F.species_tree("((S1:0.1,B:0.1):0.2,C:0.1);")
    with pytest.raises(F.FamilyError, match="Duplicate taxon labels: A"):
        F.species_tree("((A:0.1,A:0.1):0.2,C:0.1);")
    with pytest.raises(F.FamilyError, match="needs a length"):
        F.species_tree("((A,B):0.2,C:0.1);")


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ({"duplication": -1}, "family.duplication"),
        ({"loss": "0.1"}, "family.loss"),
        ({"origin": 11}, "family.origin"),
        ({"min_genes": 2}, "family.min_genes"),
        ({"speciation": 1}, "unknown speciation"),
        ({"duplications": [{"clade": ["A"], "at": 1}]}, "below 1"),
        ({"duplications": [{"clade": ["A;B"]}]}, "letters, digits"),
        ({"duplications": [{"clade": ["A", "A"]}]}, "each species once"),
        ({"losses": [{"clade": ["A"], "copy": "1a"}]}, "there is no duplication 1"),
        ({"duplications": [{"clade": ["A"]}], "losses": [{"clade": ["A"], "copy": "1c"}]}, "1b"),
    ],
)
def test_bad_family_specs_are_refused(values: dict[str, Any], message: str) -> None:
    with pytest.raises(SpecError, match=message):
        spec({"duplication": 0.1, "loss": 0.1, **values})


def test_a_family_needs_a_species_tree_file() -> None:
    data = {
        "tree": {"random": "yh", "taxa": 10, "lengths": [0.02, 0.1, 0.3]},
        "root": {"length": 50},
        "model": {"substitution": "LG"},
        "family": {"duplication": 0.1, "loss": 0.1},
        "simulation": {"seeds": [1]},
    }
    with pytest.raises(SpecError, match="must be a file"):
        parse(data)


def fake_alisim(work: Path) -> None:
    """AliSim's outputs on the hidden-paralogy gene tree: it keeps our node labels."""
    tree = (work / "tree.nwk").read_text().strip()
    names = ["Node1", "Node2", "A_1", "B_1", "Node3", "C_1", "D_1"]
    write_fasta(work / "sim.fa", {n: "MKQL" for n in names})
    (work / "sim.full.treefile").write_text(tree.replace(";", ":0;") + "\n")


@pytest.mark.parametrize("names", ["keep", "neutral"])
def test_prepare_and_collect_a_family(names: str, tmp_path: Path) -> None:
    work, out = tmp_path / "work", tmp_path / "out"
    work.mkdir()
    (tmp_path / "species.nwk").write_text(SPECIES + "\n")
    data = spec(HIDDEN_PARALOGY, names=names, image=IMAGE)
    plan = prepare(data, 1, work, tmp_path / "species.nwk")
    assert plan.commands[0][plan.commands[0].index("-t") + 1] == "tree.nwk"
    assert "Node1" in (work / "tree.nwk").read_text()
    fake_alisim(work)
    collect(data, plan, work, out, "IQ-TREE version 3.1.3", IMAGE)

    genes = ["A_1", "B_1", "C_1", "D_1"]
    simulated = [f"sim00{i}" for i in range(1, 5)] if names == "neutral" else genes
    assert list(read_fasta(out / "input.fasta")) == simulated
    assert list(read_fasta(out / "true_ancestors.fasta")) == ["Node1", "Node2", "Node3"]
    species = [f"{g}\t{s}" for g, s in zip(simulated, "ABCD", strict=True)]
    assert (out / "species.tsv").read_text() == "gene\tspecies\n" + "".join(
        f"{row}\n" for row in species
    )
    assert (out / "species_tree.nwk").read_text().startswith("((A:0.1,B:0.1)S2")
    assert json.loads((out / "family.json").read_text())["duplications_visible"] == 1
    checked = json.loads((out / "provenance.json").read_text())["sha256"]
    assert set(F.FILES) <= set(checked)
