from pathlib import Path

import pytest

from phylobot.evaluation.sequences import calibration, identity
from phylobot.seqio import read_fasta
from phylobot.trees import clades, load_tree, robinson_foulds

LEGACY = Path(__file__).parent / "legacy" / "ime2"
EASY = Path(__file__).parent / "golden" / "ime2" / "easy" / "replicates" / "101"
YARROWIA = "Yarrowia.lipolytica.IME2"
CONTROL = "Control1.IME2"


def legacy_tree(name: str) -> str:
    return (LEGACY / "trees" / name).read_text()


def test_clades_recover_2014_ancestor_definitions() -> None:
    # The 2014 cladogram is drawn rooted at Kluyveromyces. Re-rooted on the outgroup,
    # Node10 and Node11 must be the ancestors PhyloBot named AncKLSC and AncSKCTG.
    by_clade = clades(load_tree(legacy_tree("mafft_LG.cladogram.nwk")), YARROWIA)
    anc_klsc = frozenset(
        {"Saccharomyces.cerevisiae.IME2", "Kluyveromyces.lactis.IME2", "Kluyveromyces.waltii.IME2"}
    )
    anc_skctg = frozenset(read_fasta(LEGACY / "sequences.fasta")) - {YARROWIA} | {CONTROL}
    assert by_clade[anc_klsc] == "10"
    assert by_clade[anc_skctg] == "11"


def test_clades_exclude_outgroup_root() -> None:
    by_clade = clades(load_tree(legacy_tree("mafft_LG.ml.nwk"), prune=[CONTROL]), YARROWIA)
    assert all(YARROWIA not in clade for clade in by_clade)
    assert len(by_clade) == 5  # 6 ingroup leaves, rooted binary: 5 internal nodes


def test_clades_rejects_unknown_outgroup() -> None:
    with pytest.raises(ValueError, match="outgroup"):
        clades(load_tree(legacy_tree("mafft_LG.ml.nwk")), "Not.a.taxon")


def test_2014_aligners_disagree_only_on_the_control_sequence() -> None:
    mafft, muscle = legacy_tree("mafft_LG.ml.nwk"), legacy_tree("muscle_LG.ml.nwk")
    assert robinson_foulds(mafft, muscle) > 0
    assert robinson_foulds(mafft, muscle, prune=[CONTROL]) == 0


def test_true_tree_matches_itself_and_has_all_ancestors() -> None:
    true_tree = (EASY / "true_tree.nwk").read_text()
    assert robinson_foulds(true_tree, true_tree) == 0
    ancestors = read_fasta(EASY / "true_ancestors.fasta")
    labels = set(clades(load_tree(true_tree), YARROWIA).values())
    assert labels <= set(ancestors)


def test_identity() -> None:
    assert identity("ACDE", "ACDE") == 1.0
    assert identity("ACDE", "ACDF") == 0.75
    assert identity("A-CD", "ACD") == 1.0
    assert identity("acde", "ACDE") == 1.0
    with pytest.raises(ValueError, match="length mismatch"):
        identity("ACD", "ACDE")


def test_calibration_bins() -> None:
    calls = [(0.99, True), (0.97, True), (1.0, False), (0.85, True), (0.3, False)]
    bins = {(b.low, b.high): b for b in calibration(calls)}
    top = bins[(0.95, 1.0)]
    assert top.count == 3
    assert top.accuracy == pytest.approx(2 / 3)
    assert bins[(0.8, 0.9)].count == 1
    assert bins[(0.0, 0.5)].accuracy == 0.0
    assert bins[(0.5, 0.8)].count == 0


def test_check_thresholds() -> None:
    from phylobot.evaluation.score import check

    metrics = {"rf_distance": 2, "clades_true": 5, "clades_recovered": 4, "mean_identity": 0.9}
    metrics["min_identity"] = 0.8
    assert check(metrics, {"max_rf_distance": 2, "min_mean_identity": 0.85, "_note": "x"}) == []
    failures = check(metrics, {"max_rf_distance": 0, "min_clades_recovered_fraction": 1.0})
    assert len(failures) == 2
    assert check(metrics, {"bogus": 1}) == ["unknown threshold 'bogus'"]


def test_attribution_table() -> None:
    from phylobot.evaluation.score import attribution

    def metrics(identity: float, rf: int = 0) -> dict[str, float | int]:
        return {"mean_identity": identity, "rf_distance": rf, "clades_recovered": 5,
                "residue_accuracy": identity}  # fmt: skip

    by_role = {
        "full pipeline": {**metrics(0.80, rf=2), "clades_true": 5},
        "true alignment": {**metrics(0.85, rf=2), "clades_true": 5},
        "true tree": {**metrics(0.82), "clades_true": 5},
        "true alignment + tree": {**metrics(0.90), "clades_true": 5},
    }
    table = attribution(by_role)
    assert "ancestors still miss 0.100" in table
    assert "alignment + tree inference add 0.100" in table
    assert "alignment alone +0.050" in table


def test_codon_scoring_matches_nodes_without_rooting(tmp_path: Path) -> None:
    from phylobot.evaluation.codons import check_codons, codon_identity, score_codons

    truth, results = tmp_path / "truth", tmp_path / "results"
    truth.mkdir()
    results.mkdir()
    # Truth rooted at its midpoint (root 7 has no unrooted counterpart); inferred unrooted.
    (truth / "true_tree.nwk").write_text("(((a,b)5,c)6,d)7;\n")
    (truth / "true_ancestors.fasta").write_text(">5\nATGAAA\n>6\nATGAAG\n>7\nATGCGT\n")
    (results / "tree.nwk").write_text("((a,b)Node2,c,d)Node1;\n")
    (results / "ancestors.fasta").write_text(">Node1\nMK\n>Node2\nMR\n")
    (results / "ancestors.cds.fasta").write_text(">Node1\nATGAAG\n>Node2\nATGCGT\n")
    metrics = score_codons(truth, results)
    assert (metrics["nodes_true"], metrics["nodes_compared"]) == (3, 2)
    by_node = {n["true_node"]: n for n in metrics["nodes"]}
    assert by_node["6"]["identity"] == 1.0 and by_node["6"]["codon_identity"] == 1.0
    assert by_node["5"]["identity"] == 0.5 and by_node["5"]["codon_identity"] == 0.5
    assert codon_identity("ATGAAA", "ATGAAG") == 0.5
    assert check_codons(metrics, {"min_mean_identity": 0.8, "_note": "x"}) == [
        "min_mean_identity: got 0.75, threshold 0.8"
    ]
