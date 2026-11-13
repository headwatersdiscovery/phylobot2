"""ESR's pure parts: fixing a model from a report, the held-out inputs, pruning and
reinserting by split (whatever root IQ-TREE writes), reading the held-out node back.
The IQ-TREE side is in tests/golden/test_esr.py."""

from pathlib import Path

import dendropy
import pytest

from phylobot import esr
from phylobot.trees import load_tree

FIXTURES = Path(__file__).parent / "fixtures" / "iqtree"

# The Ime2 golden replicate's IQ-TREE tree: unrooted, written with a trifurcating top node.
TREE = (
    "(Yarrowia.lipolytica.IME2:1.42,(Candida.lusitaniae.IME2:0.92,(Candida.parapsilosis.IME2"
    ":0.49,Candida.albicans.IME2:0.27)Node3:0.25)Node2:0.017,(Saccharomyces.cerevisiae.IME2"
    ":0.40,(Kluyveromyces.lactis.IME2:0.41,Kluyveromyces.waltii.IME2:0.16)Node6:0.10)Node5"
    ":0.37)Node1;"
)
TIPS = [
    "Yarrowia.lipolytica.IME2", "Candida.lusitaniae.IME2", "Candida.parapsilosis.IME2",
    "Candida.albicans.IME2", "Saccharomyces.cerevisiae.IME2", "Kluyveromyces.lactis.IME2",
    "Kluyveromyces.waltii.IME2",
]  # fmt: skip


def splits(newick: str) -> dict[frozenset[str], float]:
    """Each edge as its smaller side (ties broken by sort), with its length."""
    tree = load_tree(newick)
    everyone = frozenset(leaf.taxon.label for leaf in tree.leaf_node_iter())
    out: dict[frozenset[str], float] = {}
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            continue
        side = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        key = min(side, everyone - side, key=lambda s: (len(s), sorted(s)))
        out[key] = out.get(key, 0.0) + float(node.edge.length)  # two root edges are one
    return out


def reroot_at(newick: str, tip: str) -> str:
    """The same unrooted tree written from another tip's parent, as IQ-TREE does."""
    tree = load_tree(newick)
    tree.reroot_at_node(tree.find_node_with_taxon_label(tip).parent_node,
                        suppress_unifurcations=False)  # fmt: skip
    return str(tree.as_string(schema="newick", suppress_rooting=True))


class TestFixedModel:
    def test_gamma(self) -> None:
        report = (FIXTURES / "ime2_lg_g4.iqtree").read_text()
        assert esr.model_name(report) == "LG+G4"
        assert esr.fixed_model(report) == "LG+G4{0.7923}"
        assert esr.log_likelihood(report) == -4412.9327

    def test_every_component(self) -> None:
        report = (FIXTURES / "ime2_lg_f_i_r4.iqtree").read_text()
        fixed = esr.fixed_model(report)
        assert fixed.startswith("LG+F{0.0696,0.0471,")
        assert "+I{0.0527}+R4{0.2491,0.1497,0.2389,0.5279,0.2270,0.8400,0.2323,2.7803}" in fixed
        assert fixed.count(",") == 19 + 7

    def test_unknown_component_is_refused(self) -> None:
        report = (FIXTURES / "ime2_lg_g4.iqtree").read_text()
        with pytest.raises(ValueError, match="cannot fix"):
            esr.fixed_model(
                report.replace("Model of substitution: LG+G4", "Model of substitution: LG+C20+G4")
            )


class TestHeldOut:
    def test_alignment(self) -> None:
        out = esr.held_out_alignment({"a": "AC-D", "b": "ACED", "c": "A-ED"}, "b")
        assert out == {"a": "AC-D?", "c": "A-ED?", **{d: "????L" for d in esr.DUMMIES}}

    def test_reserved_names(self) -> None:
        with pytest.raises(ValueError, match="reserved"):
            esr.held_out_alignment({"esr_daughter_1": "A", "b": "A"}, "b")

    def test_tree(self) -> None:
        tree = load_tree(esr.held_out_tree(TREE, "Candida.albicans.IME2"))
        labels = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
        assert labels == set(TIPS) - {"Candida.albicans.IME2"} | set(esr.DUMMIES)
        daughter = tree.find_node_with_taxon_label(esr.DUMMIES[0])
        node = daughter.parent_node
        assert daughter.edge.length == 1000.0
        assert [c.taxon.label for c in node.child_nodes()] == list(esr.DUMMIES)
        assert node.edge.length == pytest.approx(0.27)
        sibling = {c.taxon.label for c in node.parent_node.child_nodes() if c.taxon}
        assert sibling == {"Candida.parapsilosis.IME2"}

    def test_node_found_through_the_daughters(self) -> None:
        written = "(a:1,(esr_daughter_1:1000,esr_daughter_2:1000)Node4:0.2,b:1)Node1;"
        assert esr.held_out_node(written) == "Node4"

    def test_internal_labels_dropped(self) -> None:
        assert "Node" not in esr.plain(TREE)
        assert esr.pendant_length(TREE, "Kluyveromyces.waltii.IME2") == 0.16


class TestPruneReinsert:
    @pytest.mark.parametrize("tip", TIPS)
    def test_round_trip(self, tip: str) -> None:
        pruned, placement = esr.prune(TREE, tip)
        assert tip not in pruned
        back = esr.reinsert(pruned, tip, placement)
        assert splits(back) == pytest.approx(splits(TREE))

    @pytest.mark.parametrize("tip", TIPS)
    @pytest.mark.parametrize("root", ["Kluyveromyces.waltii.IME2", "Candida.lusitaniae.IME2"])
    def test_round_trip_after_rerooting(self, tip: str, root: str) -> None:
        # The bug the spike found: a tip that hung off the input's top node came back in
        # the wrong place once IQ-TREE wrote the refitted tree from another root.
        if root == tip:
            return
        pruned, placement = esr.prune(TREE, tip)
        back = esr.reinsert(reroot_at(pruned, root), tip, placement)
        assert splits(back) == pytest.approx(splits(TREE))

    def test_refitted_lengths_scale_the_edge(self) -> None:
        pruned, placement = esr.prune(TREE, "Candida.albicans.IME2")
        assert placement.clade == {"Candida.parapsilosis.IME2"}
        assert placement.fraction == pytest.approx(0.49 / (0.25 + 0.49))
        doubled = load_tree(pruned)
        for node in doubled.preorder_node_iter():
            if node is not doubled.seed_node:
                node.edge.length = 2 * float(node.edge.length)
        back = splits(esr.reinsert(doubled.as_string(schema="newick"), "Candida.albicans.IME2",
                                   placement))  # fmt: skip
        assert back[frozenset({"Candida.parapsilosis.IME2"})] == pytest.approx(0.98)
        assert back[frozenset({"Candida.albicans.IME2"})] == pytest.approx(0.27)  # kept

    @pytest.mark.parametrize("tip", ["b", "a"])
    def test_zero_length_edges(self, tip: str) -> None:
        # SR238 (2026-11-02): b's parent and sibling edges both 0 divided by zero. Any
        # point of a zero-length edge is the same place: the middle. a hangs off the top
        # node, whose two other edges are 0 as well.
        tree = (
            "(a:0.3,(b:0.2,c:0):0,(d:0.1,e:0.1):0);"
            if tip == "b"
            else "(a:0.3,(b:0.2,c:0.1):0,(d:0.1,e:0.1):0);"
        )
        pruned, placement = esr.prune(tree, tip)
        assert placement.fraction == 0.5
        back = splits(esr.reinsert(pruned, tip, placement))
        assert back[frozenset({tip})] == pytest.approx(0.2 if tip == "b" else 0.3)

    def test_polytomy_refused(self) -> None:
        with pytest.raises(ValueError, match="not binary"):
            esr.prune("(a:1,(b:1,c:1,d:1):1,e:1);", "b")


class TestState:
    STATE = (
        "# comment\n"
        "Node\tSite\tState\t" + "\t".join(f"p_{a}" for a in esr.AMINO_ACIDS) + "\n"
        + "".join(
            f"{node}\t{site}\tA\t" + "\t".join(["0.9"] + ["0.005"] * 19) + "\n"
            for node in ("Node1", "Node4") for site in (1, 2, 3)
        )
    )  # fmt: skip

    def test_reads_the_node_and_drops_the_extra_column(self) -> None:
        rows = esr.read_state(self.STATE, "Node4", 2)
        assert len(rows) == 2
        assert rows[0]["A"] == 0.9 and rows[1]["V"] == 0.005

    def test_missing_columns(self) -> None:
        with pytest.raises(ValueError, match="expected 4 columns"):
            esr.read_state(self.STATE, "Node4", 4)

    def test_score_skips_gaps(self) -> None:
        probs = {a: 0.0 for a in esr.AMINO_ACIDS} | {"A": 0.7, "C": 0.3}
        records = esr.score("t", "e0", "A-C", [probs, probs, probs])
        assert [(r["site"], r["true"], r["ml"], r["pp_true"]) for r in records] == [
            (1, "A", "A", 0.7),
            (3, "C", "A", 0.3),
        ]
        result = esr.TipResult("t", "e0", 0.1, -1.0, 0.0, records)
        assert result.summary()["accuracy"] == "0.5000"
        assert result.summary()["mean_pp"] == "0.7000"


def test_codon_combinations_refused(tmp_path: Path) -> None:
    (tmp_path / "alignment.fasta").write_text(">a\nATGAAA\n>b\nATGAAG\n")
    (tmp_path / "model.txt").write_text("GY+F3X4\n")
    (tmp_path / "tree.nwk").write_text("(a:1,b:1);\n")
    with pytest.raises(ValueError, match="amino-acid"):
        esr.Combination(tmp_path, tmp_path, lambda work, args: "")


def test_unknown_level(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unknown ESR levels"):
        esr.run_combination(tmp_path, tmp_path, lambda work, args: "", levels=["e2"])


def test_reroot_helper_keeps_the_tree() -> None:
    assert splits(reroot_at(TREE, "Kluyveromyces.waltii.IME2")) == pytest.approx(splits(TREE))
    assert isinstance(load_tree(TREE), dendropy.Tree)


class FakeIQTree:
    """Answers each run with a fixture report (and a tree for a fit without -blfix),
    recording the arguments."""

    def __init__(self, report: str):
        self.report = report
        self.calls: list[list[str]] = []

    def __call__(self, work: Path, args: list[str]) -> str:
        self.calls.append(args)
        prefix = args[args.index("--prefix") + 1]
        (work / f"{prefix}.iqtree").write_text(self.report)
        if "-blfix" not in args:
            (work / f"{prefix}.treefile").write_text(reroot_at(TREE, TIPS[3]) + "\n")
        return ""


def published(folder: Path, model: str) -> Path:
    (folder / "alignment.fasta").write_text("".join(f">{t}\nMKLV\n" for t in TIPS))
    (folder / "tree.nwk").write_text(TREE + "\n")
    (folder / "model.txt").write_text(model + "\n")
    return folder


@pytest.mark.parametrize(
    ("model", "fixture", "refit"),
    [("LG+G4", "ime2_lg_g4.iqtree", False), ("LG+F+I+R4", "ime2_lg_f_i_r4.iqtree", True)],
)
def test_freerate_refits_the_lengths(tmp_path: Path, model: str, fixture: str, refit: bool) -> None:
    """IQ-TREE refuses -blfix with free +R parameters ("Fixing branch lengths not supported
    under specified site rate model"), so a +R model's first fit refits the lengths on the
    topology; the fixed model's check keeps -blfix either way."""
    report = (FIXTURES / fixture).read_text()
    run = FakeIQTree(
        report.replace("Model of substitution: LG+F+I+R4", f"Model of substitution: {model}")
    )
    work = tmp_path / "work"
    work.mkdir()
    combination = esr.Combination(published(tmp_path, model), work, run)
    first, check = run.calls
    assert ("-blfix" not in first) is refit and "-blfix" in check
    assert combination.lengths_refitted is refit
    combination.save(tmp_path / "saved")
    again = esr.Combination.load(tmp_path / "saved", work, run)
    assert again.lengths_refitted is refit and again.fixed == combination.fixed
    if refit:  # the refitted tree is the one held-out tips use
        assert splits(combination.tree) == pytest.approx(splits(TREE))


def test_lnl_tolerance() -> None:
    """0.1 on small fits; 1e-5 of |lnL| on large ones (PINK1's rounding cost 0.20 at
    lnL -77,897, 2026-11-01)."""
    assert esr.lnl_tolerance(-4412.9) == 0.1
    assert esr.lnl_tolerance(-77897.29) == pytest.approx(0.779, abs=1e-3)


class ShiftedCheck(FakeIQTree):
    """The fixed model's check run reports the fit's lnL plus `shift`."""

    def __init__(self, report: str, fitted: float, shift: float):
        super().__init__(report)
        self.fitted, self.shift = fitted, shift

    def __call__(self, work: Path, args: list[str]) -> str:
        out = super().__call__(work, args)
        prefix = args[args.index("--prefix") + 1]
        if prefix == "check":
            report = (work / "check.iqtree").read_text()
            (work / "check.iqtree").write_text(report.replace(
                f"Log-likelihood of the tree: {self.fitted:.4f}",
                f"Log-likelihood of the tree: {self.fitted + self.shift:.4f}",
            ))  # fmt: skip
        return out


@pytest.mark.parametrize(
    ("fitted", "shift", "accepted"),
    [
        (-77897.2933, 0.2018, True),  # PINK1 clustalo x RAxML-NG, Q.PFAM+F+I+R7: rounding
        (-4412.9327, 0.2, False),  # the same gap on a small fit is more than rounding
        (-77897.2933, -25.0, False),  # a misread model on a large fit is still caught
    ],
)
def test_fixed_model_check_scales(tmp_path: Path, fitted: float, shift: float,
                                  accepted: bool) -> None:  # fmt: skip
    report = (
        (FIXTURES / "ime2_lg_g4.iqtree")
        .read_text()
        .replace(
            "Log-likelihood of the tree: -4412.9327", f"Log-likelihood of the tree: {fitted:.4f}"
        )
    )
    run = ShiftedCheck(report, fitted, shift)
    work = tmp_path / "work"
    work.mkdir()
    folder = published(tmp_path, "LG+G4")
    if accepted:
        assert esr.Combination(folder, work, run).log_likelihood == pytest.approx(fitted)
    else:
        with pytest.raises(ValueError, match="gives lnL"):
            esr.Combination(folder, work, run)


def test_free_rates() -> None:
    assert esr.free_rates("LG+R4") and esr.free_rates("JTT+F+I+R6")
    assert not esr.free_rates("LG+G4") and not esr.free_rates("LG+R4{0.5,0.1,0.5,1.9}")


CHECKPOINT = """--- # IQ-TREE Checkpoint ver >= 1.6
ModelMixture21:
 Component1!ModelProtein!state_freq: 0.5, 0.25, 0.25
 prop: 0.75, 0.25, 3.7e-07
RateGammaInvar:
 gamma_shape: 1.004447982
 p_invar: 0.07693949497
"""


def test_a_mixture_needs_its_weights_fitted() -> None:
    assert esr.is_mixture("LG+C20+F+I+G4") and esr.is_mixture("LG+C20")
    assert not esr.is_mixture("LG+F+I+G4") and not esr.is_mixture("LG+CF4")


def test_jobs_still_refuse_mixtures() -> None:
    assert not esr.supported("LG+C20+F+I+G4")


def test_fixed_mixture_writes_every_value_in() -> None:
    fixed = esr.fixed_mixture("JTT+C2+F+I+G4", CHECKPOINT)
    assert fixed == (
        "MIX{JTT+F{0.5,0.25,0.25}:1:0.75,JTT+FC2pi1:1:0.25,JTT+FC2pi2:1:3.7e-07}"
        "+I{0.07693949497}+G4{1.004447982}"
    )


def test_fixed_mixture_counts_its_classes() -> None:
    with pytest.raises(ValueError, match="3 weights"):
        esr.fixed_mixture("JTT+C20+G4", CHECKPOINT)
