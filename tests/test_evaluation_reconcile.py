"""Scoring reconciliation against a simulated family (phylobot.evaluation.reconcile)."""

from phylobot.evaluation import reconcile as R

# Duplication at the root; A kept copy a, B copy b, C and D both.
TRUE = "((A_1:.1,(C_1:.1,D_1:.1)Node3:.1)Node2:.1,(B_1:.1,(C_2:.1,D_2:.1)Node5:.1)Node4:.1)Node1;"
EVENTS = (
    "event\tgene_node\tspecies_branch\tdepth\tforced\tcopies\n"
    "duplication\tNode1\tS1\t0.05\tyes\t-\n"
    "speciation\tNode2\tS1\t0.1\tno\t1a\n"
    "loss\t-\tS2\t0.2\tyes\t1a\n"
    "speciation\tNode4\tS1\t0.1\tno\t1b\n"
)
ANCESTORS = {f"Node{i}": "MKV" for i in range(1, 6)}
HEADER = "node\tevent\tspecies\tlosses\tlost_in\n"


def rec(rows: list[tuple[str, str, int]]) -> str:
    return HEADER + "".join(f"{n}\t{e}\tS\t{loss}\t\n" for n, e, loss in rows)


def test_truth_lists_visible_duplications() -> None:
    truth = R.truth(TRUE, EVENTS)
    [dup] = truth.duplications
    assert dup.node == "Node1" and len(dup.below) == 6
    assert {frozenset({"A_1", "C_1", "D_1"}), frozenset({"B_1", "C_2", "D_2"})} == set(dup.copies)


def test_the_true_tree_scores_perfectly() -> None:
    truth = R.truth(TRUE, EVENTS, ANCESTORS)
    tsv = rec([("(root)", "duplication", 0), ("Node2", "speciation", 1),
               ("Node3", "speciation", 0), ("Node4", "speciation", 1),
               ("Node5", "speciation", 0)])  # fmt: skip
    rooted = TRUE.replace(")Node1;", ");")  # GeneRax's root is unlabelled
    s = R.score(truth, TRUE, (tsv, rooted), {**ANCESTORS, "Node2": "MKA"})
    assert (s.rf, s.recall, s.precision, s.chimeric) == (0, 1.0, 1.0, [])
    assert s.clades_matched == 4 and s.losses == 2
    assert s.identity is not None and abs(s.identity - (3 + 2 / 3) / 4) < 1e-9


def test_grouping_orthologs_by_species_makes_a_chimera_and_the_flag_finds_it() -> None:
    truth = R.truth(TRUE, EVENTS)
    # A_1 beside B_1, as the species tree suggests: Node7 never existed, nor did Node6
    # (both copies of C and D, without A and B).
    wrong = "((A_1:.1,B_1:.1)Node7:.1,((C_1:.1,D_1:.1)Node3:.1,(C_2:.1,D_2:.1)Node5:.1)Node6:.1);"
    tsv = rec([("(root)", "speciation", 0), ("Node7", "speciation", 0),
               ("Node6", "duplication", 0), ("Node3", "speciation", 0),
               ("Node5", "speciation", 0)])  # fmt: skip
    s = R.score(truth, wrong, (tsv, wrong))
    assert s.chimeric == ["Node6", "Node7"] and s.recall == 0.0 and s.precision == 0.0
    assert s.rf > 0
    reconciled_tsv = rec([("(root)", "duplication", 0), ("Node2", "speciation", 0),
                          ("Node3", "speciation", 0), ("Node4", "speciation", 0),
                          ("Node5", "speciation", 0)])  # fmt: skip
    found = R.flags(truth, (tsv, wrong), (reconciled_tsv, TRUE.replace(")Node1;", ");")))
    assert found.flagged == ["Node6", "Node7"] and (found.recall, found.precision) == (1.0, 1.0)


def test_pre_duplication_and_within_copy_clades_are_not_chimeric() -> None:
    truth = R.truth(TRUE, EVENTS)
    [dup] = truth.duplications
    assert not R.crosses(dup.below, dup)  # the whole duplication: its ancestor
    assert not R.crosses(frozenset({"C_1", "D_1"}), dup)
    assert R.crosses(frozenset({"C_1", "C_2"}), dup)
