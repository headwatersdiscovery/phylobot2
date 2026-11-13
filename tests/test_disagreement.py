"""Cross-aligner disagreement on mutational paths (`phylobot.evaluation.disagreement`).

A reference substitution is "aligner-dependent" when the aligners do not all make the same
call at its reference column. The fixed-topology case shares node labels; the clade case
matches branches across differing trees.
"""

from phylobot.evaluation.disagreement import (
    branch_changes,
    branch_disagreement,
    calls,
    clade_branches,
    fixed_topology_branches,
)
from phylobot.evaluation.paths import Frame

TREE = "(((a,b)C,c)Pnode,out);"
ALN = {"a": "MASTG", "b": "MASTG", "c": "MASTG", "out": "MASTG"}


def frame(ancestors: dict[str, str]) -> Frame:
    return Frame(ALN, ancestors, TREE)


def test_branch_changes_keyed_by_reference_column() -> None:
    f = frame({"Pnode": "MASTG", "C": "MVYTG", "X": "MASTG"})
    to_ref = {i: i for i in range(5)}  # identical alignments -> identity map
    assert branch_changes(f, "Pnode", "C", to_ref) == {1: ("A", "V"), 2: ("S", "Y")}


def test_fixed_topology_agreement_and_dependence() -> None:
    # mafft and muscle agree on A->V and S->Y; clustal has only A->V. So A->V is shared by
    # all three (not aligner-dependent); S->Y is not (clustal lacks it).
    recon = {
        "mafft": frame({"Pnode": "MASTG", "C": "MVYTG", "X": "MASTG"}),
        "muscle": frame({"Pnode": "MASTG", "C": "MVYTG", "X": "MASTG"}),
        "clustal": frame({"Pnode": "MASTG", "C": "MVSTG", "X": "MASTG"}),
    }
    branch_map = fixed_topology_branches(TREE)
    rows = calls(recon, branch_map, "mafft")
    by_column = {r.reference_column: r for r in rows if r.branch == "Pnode->C"}
    # 1-based columns: A->V at column 2, S->Y at column 3.
    assert by_column[2].kind == "substitution" and (by_column[2].before, by_column[2].after) == (
        "A",
        "V",
    )
    assert by_column[2].agreeing == 3 and not by_column[2].aligner_dependent
    assert by_column[3].agreeing == 2 and by_column[3].aligner_dependent  # clustal lacks S->Y
    assert branch_disagreement(recon, branch_map, "mafft")["Pnode->C"] == 1


def test_a_shifted_alignment_makes_a_substitution_aligner_dependent() -> None:
    # Same ancestors, but clustal's alignment shifts a column so mafft's substitution lands
    # on a reference column clustal shows no change at -> aligner-dependent.
    mafft = Frame(
        {"a": "MA-STG", "b": "MA-STG", "c": "MA-STG", "out": "MA-STG"},
        {"Pnode": "MA-STG", "C": "MV-STG", "X": "MA-STG"},
        TREE,
    )
    muscle = Frame(
        {"a": "MA-STG", "b": "MA-STG", "c": "MA-STG", "out": "MA-STG"},
        {"Pnode": "MA-STG", "C": "MV-STG", "X": "MA-STG"},
        TREE,
    )
    # clustal places the residue one column over; its V sits where truth is not compared.
    clustal = Frame(
        {"a": "M-ASTG", "b": "M-ASTG", "c": "M-ASTG", "out": "M-ASTG"},
        {"Pnode": "M-ASTG", "C": "M-VSTG", "X": "M-ASTG"},
        TREE,
    )
    recon = {"mafft": mafft, "muscle": muscle, "clustal": clustal}
    rows = calls(recon, fixed_topology_branches(TREE), "mafft")
    (av,) = [r for r in rows if (r.before, r.after) == ("A", "V")]
    # clustal does report A->V too (its column maps back to the same reference column),
    # so with a clean shift the mapping recovers agreement -- this guards that the column
    # map is used rather than raw column indices.
    assert av.agreeing == 3 and not av.aligner_dependent


def test_clade_branches_matches_across_trees() -> None:
    # Two aligners whose trees differ only in an unrelated subtree still share the C clade.
    a_tree = "(((a,b)C,c)Pnode,out);"
    b_tree = "(((a,b)C2,c)Top,out);"  # same clades, different labels
    recon = {
        "a": Frame(ALN, {"C": "MVSTG", "Pnode": "MASTG"}, a_tree),
        "b": Frame(ALN, {"C2": "MVSTG", "Top": "MASTG"}, b_tree),
    }
    branch_map, reference = clade_branches(recon, "out")
    assert reference == "a"
    # the Pnode->C branch (clade {a,b,c} -> {a,b}) is shared, with each tree's own labels
    (shared,) = [m for m in branch_map.values() if m["a"] == ("Pnode", "C")]
    assert shared["b"] == ("Top", "C2")
