"""Path-correctness scorer: how an inferred mutational path compares with the true one.

The hand-made cases keep the inferred and true topology identical, so clade matching is
trivial and the classification logic is what is under test. The keystone is
`test_shifted_column_is_a_placement_phantom_not_a_residue_error`: a column the aligner
lumps must read as a *placement* phantom (both residues correct, columns map apart), never
a residue error.
"""

from phylobot.evaluation import paths as P
from phylobot.evaluation.paths import (
    GAIN,
    LOSS,
    MISSED,
    PHANTOM_PLACEMENT,
    PHANTOM_RESIDUE,
    SUBSTITUTION,
    TP,
    Frame,
    change_kind,
    internal_branches,
    root_to_tip,
    score_paths,
)

# A four-taxon ingroup under an outgroup; one internal branch Pnode->C to score.
SIMPLE_TREE = "(((a,b)C,c)Pnode,out);"
# Same sequences, aligned identically on both sides (no shift): placement is impossible.
SIMPLE_ALN = {"a": "MASTG", "b": "MASTG", "c": "MASTG", "out": "MASTG"}


def classes(events: list[P.PathEvent]) -> dict[str, str]:
    """{the change's label: its classification} for quick assertions."""
    out = {}
    for e in events:
        column = e.inferred_column if e.role == "inferred" else e.true_column
        out[f"{e.before}{column}{e.after}"] = e.classification
    return out


def test_change_kind() -> None:
    assert change_kind("A", "V") == SUBSTITUTION
    assert change_kind("-", "V") == GAIN
    assert change_kind("V", "-") == LOSS


def test_true_positive_residue_phantom_and_missed() -> None:
    # Pnode->C truly: A->V at col2, S->Y at col3, T lost at col4. Inferred recovers A->V
    # (TP), misses S->Y, invents G->W at col5 (its residue is wrong -> residue phantom).
    true_anc = {"Pnode": "MASTG", "C": "MVY-G", "X": "MASTG"}
    inf_anc = {"Pnode": "MASTG", "C": "MVSTW", "X": "MASTG"}
    score = score_paths(
        Frame(SIMPLE_ALN, inf_anc, SIMPLE_TREE),
        Frame(SIMPLE_ALN, true_anc, SIMPLE_TREE),
        {"t": ["Pnode", "C"]},
        "out",
    )
    assert classes(score.events) == {
        "A2V": TP,  # both endpoints correct, same true column
        "G5W": PHANTOM_RESIDUE,  # inferred child residue wrong
        "S3Y": MISSED,  # a true substitution the reconstruction lacks
        "T4-": MISSED,  # a true loss it lacks
    }
    (status,) = score.branches
    assert status.matched
    assert status.inferred == {SUBSTITUTION: 2, GAIN: 0, LOSS: 0}
    assert status.true == {SUBSTITUTION: 2, GAIN: 0, LOSS: 1}


def test_posteriors_attach_to_events() -> None:
    # The child "C" changes at columns 2 and 5 (A2V, G5W) and misses S3Y, T4-. inferred_pp
    # gives C's ML posteriors per 1-based column; each event carries the child's pp at its
    # position (a missed change reads the pp of the unchanged call at the homologous column).
    true_anc = {"Pnode": "MASTG", "C": "MVY-G", "X": "MASTG"}
    inf_anc = {"Pnode": "MASTG", "C": "MVSTW", "X": "MASTG"}
    pp = {("C", 2): 0.99, ("C", 3): 0.40, ("C", 5): 0.80}
    score = score_paths(
        Frame(SIMPLE_ALN, inf_anc, SIMPLE_TREE),
        Frame(SIMPLE_ALN, true_anc, SIMPLE_TREE),
        {"t": ["Pnode", "C"]},
        "out",
        inferred_pp=pp,
    )
    by_label = {f"{e.before}{e.inferred_column or e.true_column}{e.after}": e for e in score.events}
    assert by_label["A2V"].pp == 0.99  # TP, inferred column 2
    assert by_label["G5W"].pp == 0.80  # residue phantom, inferred column 5
    assert by_label["S3Y"].pp == 0.40  # missed: C's unchanged call at the homologous column 3


def test_true_positive_gain() -> None:
    # A residue gained on the branch, recovered exactly: a TP gain, not a phantom.
    true_anc = {"Pnode": "MA-TG", "C": "MAWTG", "X": "MA-TG"}
    inf_anc = {"Pnode": "MA-TG", "C": "MAWTG", "X": "MA-TG"}
    aln = {"a": "MAWTG", "b": "MAWTG", "c": "MAWTG", "out": "MA-TG"}
    score = score_paths(
        Frame(aln, inf_anc, SIMPLE_TREE),
        Frame(aln, true_anc, SIMPLE_TREE),
        {"t": ["Pnode", "C"]},
        "out",
    )
    (event,) = [e for e in score.events if e.role == "inferred"]
    assert event.classification == TP and event.kind == GAIN
    assert event.before == "-" and event.after == "W"


# --- the keystone: a shifted column is a placement phantom ----------------------------

# Five ingroup taxa so the parent clade can outvote the child clade. Identical topology
# both sides; the TRUE alignment separates a,b's residue (K, col1) from c,d,e's (R, col2),
# the INFERRED alignment lumps both into one column. Over {a,b} that column maps to true
# col1; over {a,b,c,d,e} the three R's carry it to true col2.
SHIFT_TREE = "(((a,b)C,(c,d,e)X)Pnode,out);"
TRUE_ALN = {"a": "MK-Z", "b": "MK-Z", "c": "M-RZ", "d": "M-RZ", "e": "M-RZ", "out": "M--Z"}
INF_ALN = {"a": "MKZ", "b": "MKZ", "c": "MRZ", "d": "MRZ", "e": "MRZ", "out": "M-Z"}


def test_shifted_column_is_a_placement_phantom_not_a_residue_error() -> None:
    # Inferred Pnode has R, C has K, at the lumped column: an R->K change. Both residues are
    # correct for their nodes (Pnode's R is true at col2, C's K is true at col1), but the
    # columns map apart -> the misalignment invented this substitution.
    true_anc = {"Pnode": "M-RZ", "C": "MK-Z", "X": "M-RZ"}
    inf_anc = {"Pnode": "MRZ", "C": "MKZ", "X": "MRZ"}
    score = score_paths(
        Frame(INF_ALN, inf_anc, SHIFT_TREE),
        Frame(TRUE_ALN, true_anc, SHIFT_TREE),
        {"t": ["Pnode", "C"]},
        "out",
    )
    (inferred,) = [e for e in score.events if e.role == "inferred"]
    assert inferred.classification == PHANTOM_PLACEMENT
    assert inferred.classification != PHANTOM_RESIDUE
    assert inferred.parent_correct and inferred.child_correct  # both residues are right
    assert (inferred.before, inferred.after) == ("R", "K")
    # In truth the branch gains K (col1) and loses R (col2) at two columns; both are missed.
    missed = {(e.kind, e.true_column) for e in score.events if e.role == MISSED}
    assert missed == {(GAIN, 2), (LOSS, 3)}


def test_placement_phantom_vanishes_under_the_true_alignment() -> None:
    # The same ancestors scored against the TRUE alignment (no shift): the change becomes a
    # real event (or disappears), never a placement phantom. This is why arm E has zero.
    true_anc = {"Pnode": "M-RZ", "C": "MK-Z", "X": "M-RZ"}
    # Give the inferred side the true alignment; ancestors on the true columns.
    inf_anc = {"Pnode": "M-RZ", "C": "MK-Z", "X": "M-RZ"}
    score = score_paths(
        Frame(TRUE_ALN, inf_anc, SHIFT_TREE),
        Frame(TRUE_ALN, true_anc, SHIFT_TREE),
        {"t": ["Pnode", "C"]},
        "out",
    )
    assert not any(e.classification == PHANTOM_PLACEMENT for e in score.events)


def test_unmatched_branch_is_reported_not_dropped() -> None:
    # A node whose clade the truth does not contain cannot be scored, but the branch is
    # still reported (its inferred changes counted) so the driver can report the fraction.
    true_tree = "(((a,c)C2,b)Pnode,out);"  # truth groups a,c, not a,b
    true_anc = {"Pnode": "MASTG", "C2": "MVSTG", "X": "MASTG"}
    inf_anc = {"Pnode": "MASTG", "C": "MVSTG", "X": "MASTG"}
    score = score_paths(
        Frame(SIMPLE_ALN, inf_anc, SIMPLE_TREE),
        Frame(SIMPLE_ALN, true_anc, true_tree),
        {"t": ["Pnode", "C"]},
        "out",
    )
    assert score.events == []  # nothing scorable
    (status,) = score.branches
    assert not status.matched
    assert status.inferred == {SUBSTITUTION: 1, GAIN: 0, LOSS: 0}
    assert status.true == {}


def test_depth_accumulates_along_a_trajectory() -> None:
    # A two-branch path Pnode->C->a: each branch keeps its own depth index.
    true_anc = {"Pnode": "MASTG", "C": "MVSTG", "X": "MASTG"}
    inf_anc = {"Pnode": "MASTG", "C": "MVSTG", "X": "MASTG"}
    tree = "((((a)C,b)Pnode,c),out);"  # a sits below C
    score = score_paths(
        Frame({"a": "MVSTG", "b": "MASTG", "c": "MASTG", "out": "MASTG"}, inf_anc, tree),
        Frame({"a": "MVSTG", "b": "MASTG", "c": "MASTG", "out": "MASTG"}, true_anc, tree),
        {"deep": ["Pnode", "C", "a"]},
        "out",
    )
    depths = sorted(b.depth for b in score.branches)
    assert depths == [0, 1]


# --- trajectory enumeration -----------------------------------------------------------


def test_internal_branches_excludes_tips_by_default() -> None:
    tree = "(((a,b)C,c)Pnode,out);"
    internal = internal_branches(tree)
    assert internal == {"Pnode->C": ["Pnode", "C"]}
    with_tips = internal_branches(tree, tips=True)
    assert with_tips["Pnode->c"] == ["Pnode", "c"]
    assert "C->a" in with_tips and "C->b" in with_tips


def test_root_to_tip_paths() -> None:
    tree = "(((a,b)C,c)Pnode,out);"
    paths = root_to_tip(tree, "Pnode")
    assert paths["a"] == ["Pnode", "C", "a"]
    assert paths["c"] == ["Pnode", "c"]
    assert "out" not in paths  # the outgroup is not below the ingroup root
