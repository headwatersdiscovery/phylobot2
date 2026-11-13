"""phylobot.alignments: upload checks, the editor's edit state and its edits."""

import pytest

from phylobot import alignments, codons
from phylobot.alignments import AlignmentError, Edits, Row
from phylobot.models import Limits

LIMITS = Limits(max_sequences=10, max_length=50, max_upload_bytes=10_000, max_concurrent_jobs=1)
CODE = codons.STANDARD


def rows(*pairs: tuple[str, str]) -> list[Row]:
    return [Row(name, "", row) for name, row in pairs]


def test_ranges_round_trip_and_ignore_junk() -> None:
    assert alignments.parse_ranges("2-4, 7,9-8", 20) == {1, 2, 3, 6, 7, 8}
    assert alignments.parse_ranges("0-3,18-99,x,5-,--,-3", 20) == {0, 1, 2, 17, 18, 19}
    assert alignments.parse_ranges("", 5) == set()
    assert alignments.format_ranges({1, 2, 3, 6, 8, 9}) == "2-4,7,9-10"
    assert alignments.format_ranges([]) == ""
    assert alignments.parse_ranges("9999999999-1", 3) == set()  # too many digits: ignored


def test_parse_checks_an_upload() -> None:
    text = ">b desc here\nMK-L\n>a\nM.KL\n"
    parsed = alignments.parse(text, "protein", CODE, LIMITS)
    assert parsed == [Row("b", "desc here", "MK-L"), Row("a", "", "M-KL")]
    assert alignments.to_fasta(parsed) == ">b desc here\nMK-L\n>a\nM-KL\n"
    bad = {
        ">a\nMKL\n>b\nMK\n": "different lengths",
        ">a\nMKL\n": "at least 2 sequences",
        ">a\nMK1\n>b\nMKL\n": "invalid characters: 1",
        ">a\n---\n>b\nMKL\n": "all gaps",
        ">a\nMKL\n>a\nMKL\n": "duplicate sequence name",
        "MKL\n>a\nMKL\n": "before the first '>'",
        "".join(f">s{i}\nM\n" for i in range(11)): "at most 10 sequences",
        f">a\n{'M' * 51}\n>b\n{'M' * 51}\n": "51 residues; your limit is 50",
    }
    for text, message in bad.items():
        with pytest.raises(AlignmentError, match=message):
            alignments.parse(text, "protein", CODE, LIMITS)


def test_parse_codon_alignments() -> None:
    parsed = alignments.parse(">a\nATG---AAA\n>b\nATGCCC---\n", "codon", CODE, LIMITS)
    assert alignments.site_count(parsed, "codon") == 3
    assert alignments.amino_rows(parsed, "codon", CODE) == ["M-K", "MP-"]
    with pytest.raises(AlignmentError, match="gap inside codon"):
        alignments.parse(">a\nATG-AAAAA\n>b\nATGCCCAAA\n", "codon", CODE, LIMITS)
    with pytest.raises(AlignmentError, match="not a multiple of 3"):
        alignments.parse(">a\nATGA\n>b\nATGC\n", "codon", CODE, LIMITS)
    with pytest.raises(AlignmentError, match="choose Protein"):
        alignments.parse(">a\nATGMKL\n>b\nATGCCC\n", "codon", CODE, LIMITS)


ALN = rows(("seq10", "MKAL-"), ("seq2", "MRAL-"), ("Alpha", "M-GLV"))


def test_edits_from_the_url_never_empty_the_alignment() -> None:
    assert Edits.from_query({"cut": "1-5"}, 3, 5) == Edits()
    assert Edits.from_query({"drop": "1-3"}, 3, 5) == Edits()
    edits = Edits.from_query({"cut": "2,9", "drop": "3", "sort": "evil"}, 3, 5)
    assert edits == Edits(frozenset({1}), frozenset({2}), "")
    assert edits.query() == {"cut": "2", "drop": "3"}
    assert Edits().query() == {} and not Edits().changed


def test_apply_cuts_drops_and_sorts_naturally() -> None:
    shown = alignments.apply(ALN, Edits(frozenset({1, 4}), sort="name"), "protein")
    assert [(r.name, r.row) for r in shown] == [("Alpha", "MGL"), ("seq2", "MAL"),
                                                 ("seq10", "MAL")]  # fmt: skip
    shown = alignments.apply(ALN, Edits(drop=frozenset({0}), sort="-name"), "protein")
    assert [r.name for r in shown] == ["seq2", "Alpha"]
    codon = rows(("a", "ATGCCCAAA"), ("b", "ATG---AAA"))
    assert [r.row for r in alignments.apply(codon, Edits(frozenset({1})), "codon")] == [
        "ATGAAA", "ATGAAA"]  # fmt: skip


def test_one_off_edits_name_positions_as_shown() -> None:
    sorted_view = Edits(frozenset({0}), sort="name")  # shown: Alpha, seq2, seq10; cols 2-5
    edits, note = alignments.edit(ALN, sorted_view, {"cut_columns": "1-2"}, "protein")
    assert edits.cut == {0, 1, 2} and not note  # shown columns 1-2 were stored 2-3
    edits, _ = alignments.edit(ALN, sorted_view, {"drop_rows": "1"}, "protein")
    assert edits.drop == {2}  # Alpha, stored third
    edits, _ = alignments.edit(ALN, sorted_view, {"keep_rows": "2-3"}, "protein")
    assert edits.drop == {2}
    edits, _ = alignments.edit(ALN, sorted_view, {"keep_columns": "3-4"}, "protein")
    assert edits.cut == {0, 1, 2}
    # After removing Alpha, column 5 (and nothing else) is all gaps.
    edits, _ = alignments.edit(ALN, Edits(drop=frozenset({2})), {"cut_gaps": "1"}, "protein")
    assert edits.cut == {4}
    # Cut down to column 5, seq10 and seq2 have nothing left (Alpha is shown first).
    only_five = Edits(frozenset({0, 1, 2, 3}), sort="name")
    edits, _ = alignments.edit(ALN, only_five, {"drop_empty": "1"}, "protein")
    assert edits.drop == {0, 1}
    edits, _ = alignments.edit(ALN, Edits(), {"sort_toggle": "1"}, "protein")
    assert edits.sort == "name"
    edits, _ = alignments.edit(ALN, edits, {"sort_toggle": "1"}, "protein")
    assert edits.sort == "-name"
    edits, _ = alignments.edit(ALN, edits, {"sort_reset": "1"}, "protein")
    assert edits.sort == ""
    edits, note = alignments.edit(ALN, sorted_view, {"cut_columns": "1-4"}, "protein")
    assert edits == sorted_view and "every column" in note
    edits, note = alignments.edit(ALN, Edits(), {"drop_rows": "1-3"}, "protein")
    assert edits == Edits() and "every sequence" in note


def test_provenance_maps_back_to_the_parent() -> None:
    edits = Edits(frozenset({0, 3}), frozenset({1}), "name")
    assert alignments.provenance(ALN, edits, "protein", 3) == {
        "parent_sites": 5, "kept_sites": "2-3,5", "removed_sequences": ["seq2"], "sort": "name",
        "parent_version": 3,
    }  # fmt: skip


def test_compose_keeps_mapping_to_the_first_alignment() -> None:
    # An upload (no record yet), overwritten twice: columns map to the uploaded file's.
    first = alignments.compose({}, ALN, Edits(frozenset({0}), frozenset({2})), "protein")
    assert first == {"parent_sites": 5, "kept_sites": "2-5", "removed_sequences": ["Alpha"],
                     "sort": ""}  # fmt: skip
    now = alignments.apply(ALN, Edits(frozenset({0}), frozenset({2})), "protein")  # 2 x 4
    second = alignments.compose(first, now, Edits(frozenset({1}), sort="-name"), "protein")
    assert second == {"parent_sites": 5, "kept_sites": "2,4-5",
                      "removed_sequences": ["Alpha"], "sort": "-name"}  # fmt: skip
    # A copy's record keeps its parent and the parent's version.
    copy = alignments.provenance(ALN, Edits(frozenset({4})), "protein", 2)
    rows4 = alignments.apply(ALN, Edits(frozenset({4})), "protein")
    again = alignments.compose(copy, rows4, Edits(frozenset({0})), "protein")
    assert again["kept_sites"] == "2-4" and again["parent_version"] == 2
    # A record that does not describe the rows (4 kept, 3 here) starts afresh.
    three = alignments.apply(ALN, Edits(frozenset({0, 1})), "protein")
    assert alignments.compose(copy, three, Edits(), "protein")["parent_sites"] == 3
