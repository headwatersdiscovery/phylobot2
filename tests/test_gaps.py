from pathlib import Path

import pytest

from phylobot import gaps
from phylobot.pipeline.export import place_gaps
from phylobot.seqio import read_fasta

# ((a,b)X,(c,d)Y)R  as an unrooted tree: X joins a, b and {c,d}; Y joins c, d and {a,b}.
TREE = "((a:1,b:1)X:1,(c:1,d:1)Y:1)R;"


def test_shared_presence_rule() -> None:
    alignment = {
        #      col: 0 1 2 3 4
        "a": "MKQ-A",  # col 2: insertion in a only; col 3: deletion in a only
        "b": "MK--A",  # col 4: residues in a and b only
        "c": "MK-TW",
        "d": "MK-T-",
    }
    presence = gaps.parsimony_presence(TREE, alignment)
    # X sees subtrees {a}, {b}, {c,d}; Y sees {c}, {d}, {a,b}.
    assert presence["X"] == [True, True, False, False, True]
    assert presence["Y"] == [True, True, False, True, True]
    # col 2 (a only): absent everywhere - an insertion is never pushed into ancestors.
    # col 3 (c, d only): present at Y (both its children), absent at X (only one side).
    # col 4 (a, b, c): present at both.


def test_shared_presence_single_lineage_insertion_at_the_root() -> None:
    alignment = {"a": "MKQ", "b": "MK-", "c": "MK-", "d": "MK-"}
    presence = gaps.parsimony_presence(TREE, alignment)
    assert all(not presence[node][2] for node in presence)


def test_read_binary_state_and_match_nodes(tmp_path: Path) -> None:
    state = tmp_path / "gaps.state"
    state.write_text(
        "# comment\nNode\tSite\tState\tp_0\tp_1\n"
        "N1\t1\t1\t0.1\t0.9\nN1\t2\t0\t0.8\t0.2\nN2\t1\t1\t0.3\t0.7\nN2\t2\t1\t0.4\t0.6\n"
    )
    assert gaps.read_binary_state(state) == {"N1": [0.1, 0.8], "N2": [0.3, 0.4]}
    # Same topology, different names and rooting: matched by the subtrees meeting at each
    # node. X and N1 both join {a}, {b}, {c,d}; the rooted tree's R joins only two
    # subtrees (an added root) and matches nothing. (Matching by a single split would
    # have mapped X and Y to the same node.)
    asr = "((a,b)X,(c,d)Y)R;"
    gap = "(a,b,(c,d)N2)N1;"
    assert gaps.match_nodes(asr, gap) == {"X": "N1", "Y": "N2"}


def test_ml_presence_uses_half_as_the_cut(tmp_path: Path) -> None:
    state = tmp_path / "gaps.state"
    state.write_text("Node\tSite\tState\tp_0\tp_1\nG\t1\t1\t0.5\t0.5\nG\t2\t0\t0.51\t0.49\n")
    presence, p_gap = gaps.ml_presence("((a,b)X,c,d);", "((a,b)G,c,d);", state)
    assert presence == {"X": [True, False]} and p_gap == {"X": [0.5, 0.51]}


def test_mask() -> None:
    masked = gaps.mask({"X": "MKTA", "Y": "MKTA"}, {"X": [True, False, True, False]})
    assert masked == {"X": "M-T-", "Y": "MKTA"}  # nodes without information unmasked
    with pytest.raises(ValueError, match="columns"):
        gaps.mask({"X": "MK"}, {"X": [True]})


def test_export_writes_one_masked_file_per_method(tmp_path: Path) -> None:
    alignment = {"a": "MKQ-A", "b": "MK--A", "c": "MK-TW", "d": "MK-T-"}
    ancestors = {"X": "MKQTA", "Y": "MKQTA"}
    place_gaps(ancestors, alignment, TREE, tmp_path, ("parsimony",), None, None)
    assert read_fasta(tmp_path / "ancestors.gaps-parsimony.fasta") == {"X": "MK--A", "Y": "MK-TA"}
    assert not (tmp_path / "ancestors.gaps-ml.fasta").exists()
    with pytest.raises(ValueError, match="needs --gap-state"):
        place_gaps(ancestors, alignment, TREE, tmp_path, ("ml",), None, None)
    with pytest.raises(ValueError, match="unknown gap method"):
        place_gaps(ancestors, alignment, TREE, tmp_path, ("lazarus",), None, None)


def test_external_ancestors_renamed_and_root_dropped(tmp_path: Path) -> None:
    """PRANK-style output: a rooted tree with an extra root node and #n# names."""
    ours = "(a,b,(c,d)Y)X;"
    theirs = "((a,b)#1#,(c,d)#2#)#3#;"  # #3# is an added root; #1# == our X, #2# == Y
    fasta = tmp_path / "anc.fas"
    fasta.write_text(">a\nMK-\n>#1#\nmkt\n>#2#\nM-T\n>#3#\nMKT\n")
    alignment = {"a": "MK-", "b": "MKT", "c": "M-T", "d": "MKT"}
    renamed = gaps.external_ancestors(ours, theirs, fasta, alignment)
    assert renamed == {"X": "MKT", "Y": "M-T"}  # uppercased; added root dropped
    with pytest.raises(ValueError, match="columns"):
        gaps.external_ancestors(ours, theirs, fasta, {k: v + "AA" for k, v in alignment.items()})


def test_external_ancestors_restore_all_gap_columns(tmp_path: Path) -> None:
    """PRANK drops columns where every sequence is a gap; they come back as gaps."""
    ours = "(a,b,(c,d)Y)X;"
    alignment = {"a": "M-K-T", "b": "M-KAT", "c": "M-K-T", "d": "M-KAT"}  # column 1 all gap
    fasta = tmp_path / "anc.fas"
    fasta.write_text(">X\nMKAT\n>Y\nMK-T\n")  # 4 columns: column 1 dropped
    renamed = gaps.external_ancestors(ours, "(a,b,(c,d)Y)X;", fasta, alignment)
    assert renamed == {"X": "M-KAT", "Y": "M-K-T"}
    assert gaps.all_gap_columns(alignment) == [1]
    assert gaps.restore_columns("ABC", [0, 2, 5]) == "-A-BC-"


def test_fitch_classic_cases() -> None:
    # rooted ((a,b)X,(c,d)Y)R
    one_each = {"a": "Q", "b": "-", "c": "Q", "d": "-"}  # residues in a and c only
    fitch = gaps.fitch_presence(TREE, one_each)
    shared = gaps.parsimony_presence(TREE, one_each)
    assert fitch == {"X": [False], "Y": [False], "R": [False]}  # two insertions; root tie -> gap
    assert shared == {"X": [True], "Y": [True], "R": [True]}  # residues on both sides
    both_ab = {"a": "Q", "b": "Q", "c": "-", "d": "-"}
    assert gaps.fitch_presence(TREE, both_ab) == {"X": [True], "Y": [False], "R": [False]}
    everywhere = {"a": "Q", "b": "Q", "c": "Q", "d": "-"}
    assert gaps.fitch_presence(TREE, everywhere) == {"X": [True], "Y": [True], "R": [True]}


def test_fitch_multifurcation_uses_hartigan() -> None:
    # R has three children: two with a residue, one without -> residue wins 2:1
    tree = "(a:1,b:1,(c:1,d:1)Y:1)R;"
    assert gaps.fitch_presence(tree, {"a": "Q", "b": "Q", "c": "-", "d": "-"}) == {
        "Y": [False], "R": [True]
    }  # fmt: skip


def test_fitch_skips_an_unlabelled_root() -> None:
    tree = "((a:1,b:1)X:1,(c:1,d:1)Y:1);"
    assert set(gaps.fitch_presence(tree, {"a": "Q", "b": "Q", "c": "Q", "d": "Q"})) == {"X", "Y"}
