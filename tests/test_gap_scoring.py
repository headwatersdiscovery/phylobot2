import pytest

from phylobot.evaluation.gaps import GapCalls, aligned_identity, column_map, gap_calls

# The same three sequences aligned two ways: the inferred alignment spreads an insertion
# over its own column, and shifts one residue.
TRUE = {"a": "MK-TV", "b": "MKATV", "c": "M--TV"}
INFERRED = {"a": "MKT-V", "b": "MKATV", "c": "M-T-V"}


def test_column_map_follows_the_residues() -> None:
    mapping = column_map(INFERRED, TRUE)
    # inferred col 2: residues T(a) T(c) truly at col 3, A(b) at col 2 -> majority col 3
    assert mapping == {0: 0, 1: 1, 2: 3, 3: 3, 4: 4}


def test_column_map_needs_the_same_residues() -> None:
    with pytest.raises(ValueError, match="different residue counts"):
        column_map({"a": "MKT"}, {"a": "MK"})


def test_column_map_leaves_residue_free_columns_unmapped() -> None:
    assert column_map({"a": "M-K", "b": "M-K"}, {"a": "MK", "b": "MK"}) == {0: 0, 2: 1}


def test_gap_calls_confusion() -> None:
    mapping = {0: 0, 1: 1, 2: 2, 3: 3}
    calls = gap_calls("M-KT", "M--TA", mapping)  # true col 4 holds a residue no column maps to
    assert (calls.true_residue, calls.true_gap, calls.missed_gap, calls.false_gap) == (2, 1, 1, 0)
    assert calls.lost_residues == 1
    assert calls.gap_recall == 0.5 and calls.gap_precision == 1.0 and calls.accuracy == 0.75


def test_gap_free_reconstruction_scores_zero_recall() -> None:
    calls = gap_calls("MKAT", "M--T", {0: 0, 1: 1, 2: 2, 3: 3})
    assert calls.gap_recall == 0.0 and calls.gap_precision is None  # no gaps called


def test_gap_calls_add_up() -> None:
    total = GapCalls()
    total += GapCalls(true_gap=1, true_residue=2)
    total += GapCalls(false_gap=1, lost_residues=3)
    assert (total.columns, total.lost_residues, total.accuracy) == (4, 3, 0.75)
    assert GapCalls().accuracy is None and GapCalls().gap_recall is None


@pytest.mark.parametrize(
    ("inferred", "truth", "expected"),
    [
        ("MKTAYIAK", "MKTAYIAK", 1.0),
        ("MKTAYIAK", "mktayiak", 1.0),  # case-insensitive
        ("MK-TAYIAK", "MKTAYIAK", 1.0),  # gaps ignored
        ("MKTAYIAKQQQQ", "MKTAYIAK", 8 / 12),  # extra residues lower identity
        ("MKTAYIA", "MKTAYIAK", 7 / 8),  # missing residues too
        ("MKTXYIAK", "MKTAYIAK", 7 / 8),
        ("", "MKT", 0.0),
    ],
)
def test_aligned_identity(inferred: str, truth: str, expected: float) -> None:
    assert aligned_identity(inferred, truth) == pytest.approx(expected)


def test_aligned_identity_finds_shifted_blocks() -> None:
    # A 3-residue insertion in the middle: positional comparison would collapse after it.
    assert aligned_identity("MKTAQQQYIAKLV", "MKTAYIAKLV") == pytest.approx(10 / 13)


def test_equal_length_with_different_gaps_is_not_positional() -> None:
    """Same length, gaps in different places: positional comparison would shift every
    residue after the first gap; the alignment-based identity does not."""
    inferred, truth = "MKTA-YIAKL", "MK-TAYIAKL"  # same 9 residues, different gap spot
    residues = inferred.replace("-", "")
    assert residues == truth.replace("-", "") and aligned_identity(inferred, truth) == 1.0
