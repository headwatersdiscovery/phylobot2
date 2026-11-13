"""phylobot.alignment_map: columns mapped between alignments of the same sequences."""

import pytest

from phylobot import alignment_map

# Two aligners on the same three sequences. They agree on columns 0 and 3 of A; A puts
# b's K and c's K in one column (1) where B splits them (1 and 2).
A = {"a": "MK-L", "b": "MKAL", "c": "MK-L"}
B = {"a": "M-K-L", "b": "MK-AL", "c": "M-K-L"}


def test_map_columns_follows_the_residues_and_scores_agreement() -> None:
    mapped = alignment_map.map_columns(A, B)
    assert mapped.target == {0: 0, 1: 2, 2: 3, 3: 4}
    # Column 1 of A holds three K's; B's column 2 holds two of them (a's, c's), column 1 b's.
    assert mapped.agreement[1] == pytest.approx(2 / 3)
    assert mapped.exact(0) and mapped.exact(2) and mapped.exact(3) and not mapped.exact(1)
    assert mapped.reference_agreement == {0: 1.0, 2: pytest.approx(2 / 3), 3: 1.0, 4: 1.0}


def test_map_columns_is_exact_between_identical_alignments() -> None:
    mapped = alignment_map.map_columns(A, A)
    assert mapped.target == {0: 0, 1: 1, 2: 2, 3: 3}
    assert all(mapped.exact(c) for c in mapped.target)


def test_mapping_needs_the_same_sequences() -> None:
    assert alignment_map.same_sequences(A, B)
    assert not alignment_map.same_sequences(A, {"a": "MKL", "b": "MKAL"})
    assert not alignment_map.same_sequences(A, {"a": "MKL", "b": "MKAL", "c": "MQL"})
    with pytest.raises(ValueError, match="different residue counts"):
        alignment_map.map_columns({"a": "MKT"}, {"a": "MK"})


def test_sequence_frame_places_columns_on_one_sequence() -> None:
    frame = alignment_map.sequence_frame("M-K-L")
    assert frame.length == 3
    assert frame.residue == [1, None, 2, None, 3]
    # Column 1 is an insertion after residue 1, column 3 after residue 2.
    assert frame.slot == [0, 1, 1, 2, 2]


def test_residue_agreement_finds_where_the_aligners_differ() -> None:
    # On sequence b: M and L agree; K and A are where the aligners part.
    scores = alignment_map.residue_agreement([A, B], "b")
    assert scores[0] == 1.0 and scores[3] == 1.0
    assert scores[1] == pytest.approx(1 / 3)  # A: {aK, bK, cK}; B: {bK}
    assert scores[2] == 1.0  # A's column 2 and B's column 3 both hold b's A alone
    assert alignment_map.residue_agreement([A], "a") == [1.0, 1.0, 1.0]
