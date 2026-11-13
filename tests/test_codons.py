"""Genetic codes, coding-sequence checks and codon alignments (phylobot.codons)."""

import pytest

from phylobot import codons
from phylobot.codons import CODES, STANDARD, CodonError


def test_tables_are_whole_and_known_reassignments_hold() -> None:
    assert len(codons.CODONS) == 64
    for code in codons.GENETIC_CODES:
        assert len(code.amino) == 64 and set(code.amino) <= set("ACDEFGHIKLMNPQRSTVWY*")
        assert code.iqtree or code.paml_icode is not None  # every code some tool runs
    std = STANDARD.table
    assert [std[c] for c in ("TAA", "TAG", "TGA", "ATG", "AGA")] == ["*", "*", "*", "M", "R"]
    vertebrate_mt = CODES[2].table
    assert [vertebrate_mt[c] for c in ("TGA", "AGA", "AGG", "ATA")] == ["W", "*", "*", "M"]
    assert CODES[6].table["TAA"] == "Q"  # ciliate nuclear
    assert CODES[11].amino == STANDARD.amino  # bacterial: only start codons differ


def test_tool_support_follows_iqtree_and_paml() -> None:
    assert STANDARD.tools == ["IQ-TREE", "PAML"] and STANDARD.paml_icode == 0
    assert CODES[11].paml_icode == 0  # same amino acids as the standard code
    assert CODES[15].tools == ["PAML"]  # Blepharisma: not in IQ-TREE
    assert CODES[25].tools == ["IQ-TREE"]
    assert {c.ncbi for c in codons.GENETIC_CODES if c.paml_icode is not None} == {
        1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 15,
    }  # fmt: skip
    with pytest.raises(CodonError, match="unknown genetic code"):
        codons.genetic_code(7)  # retired by NCBI
    with pytest.raises(CodonError):
        codons.genetic_code("1; rm -rf /")


def test_translation() -> None:
    assert codons.translate("ATGAAATGG", STANDARD) == "MKW"
    assert codons.translate("ATGANNTGA", STANDARD) == "MX*"
    assert codons.translate("ATGTGA", CODES[2]) == "MW"
    with pytest.raises(CodonError, match="multiple of 3"):
        codons.translate("ATGA", STANDARD)


def test_check_cds_removes_final_stop_and_refuses_internal_ones() -> None:
    assert codons.check_cds("ATGAAATAA", STANDARD) == ("ATGAAA", "MK", ["final stop codon removed"])
    assert codons.check_cds("ATGNAA", STANDARD) == ("ATGNAA", "MX", [])  # ambiguity: no note
    with pytest.raises(CodonError, match="internal stop codon TGA at codon 2"):
        codons.check_cds("ATGTGAAAA", STANDARD)
    assert codons.check_cds("ATGTGAAAA", CODES[2])[1] == "MWK"  # TGA is W in vertebrate mt
    with pytest.raises(CodonError, match="frameshift"):
        codons.check_cds("ATGAAAT", STANDARD)
    with pytest.raises(CodonError, match="no codons"):
        codons.check_cds("TAA", STANDARD)


def test_ambiguity_codes_translate_when_every_reading_agrees() -> None:
    assert codons.readings("CTR") == ["CTA", "CTG"]
    assert [codons.translate_codon(c, STANDARD) for c in ("CTR", "CTN", "ACY", "AAY")] == [
        "L", "L", "T", "N",
    ]  # fmt: skip
    assert codons.translate_codon("AMC", STANDARD) == "X"  # AAC (N) or ACC (T)
    assert codons.possible_aminos("AMC", STANDARD) == "NT"
    assert codons.translate_codon("TAR", STANDARD) == "*"  # TAA or TAG: a stop either way
    assert codons.translate_codon("YAA", STANDARD) == "X"  # TAA (stop) or CAA (Q)
    assert codons.translate_codon("NNN", STANDARD) == "X"
    with pytest.raises(CodonError, match="internal stop codon TAR"):
        codons.check_cds("ATGTARAAA", STANDARD)
    assert codons.check_cds("ATGAAATRA", STANDARD)[0] == "ATGAAA"  # final TAA/TGA removed
    found = codons.ambiguous_codons("ATGAMCCTRNNNAAA", STANDARD)
    assert [(a.number, a.codon, a.aminos, a.translation, a.unknown) for a in found] == [
        (2, "AMC", "NT", "X", False), (3, "CTR", "L", "L", False),
        (4, "NNN", "*ACDEFGHIKLMNPQRSTVWY", "X", True),
    ]  # fmt: skip
    assert codons.compatible("AAC", "AMC") and codons.compatible("AMC", "ARC")
    assert not codons.compatible("AAT", "AMC")


def test_restore_codons_undoes_prank_s_nnn() -> None:
    # PRANK -codon writes NNN for a codon with an ambiguity code; the user's comes back.
    assert codons.restore_codons("ATG---NNNAAA", "ATGAMCAAA") == "ATG---AMCAAA"
    assert codons.restore_codons("ATGNNN", "ATGNNN") == "ATGNNN"
    with pytest.raises(CodonError, match="differs"):
        codons.restore_codons("ATGNNNAAA", "ATGAACAAA")  # an unambiguous codon replaced
    with pytest.raises(CodonError, match="differs"):
        codons.restore_codons("ATG---AAA", "ATGAMCAAA")  # a codon lost


def test_aligned_rows_translate_in_whole_codons() -> None:
    assert codons.translate_row("ATG---AAA", STANDARD) == "M-K"
    with pytest.raises(CodonError, match="gap inside codon 2"):
        codons.translate_row("ATGA--AAA", STANDARD)


def test_backtranslation_puts_each_sequence_s_own_codons_back() -> None:
    assert codons.backtranslate("M--KX", "ATGAAGNNN", STANDARD) == "ATG------AAGNNN"
    assert codons.backtranslate("m-k", "ATGAAG", STANDARD) == "ATG---AAG"
    with pytest.raises(CodonError, match="residue R but codon 2 is AAG"):
        codons.backtranslate("M-R", "ATGAAG", STANDARD)  # an aligner changed a residue
    with pytest.raises(CodonError, match="left over"):
        codons.backtranslate("M--", "ATGAAG", STANDARD)
    with pytest.raises(CodonError, match="more residues"):
        codons.backtranslate("MKK", "ATGAAG", STANDARD)


def test_amino_posteriors_choose_the_amino_acid_first() -> None:
    # CTT (L) is the single most probable codon, but I (ATT + ATC) is the most probable
    # amino acid: the ancestor is I, coded by I's best codon.
    probs = {"CTT": 0.40, "ATT": 0.35, "ATC": 0.25}
    amino, codon = codons.amino_posteriors(probs, STANDARD)
    assert list(amino) == list(codons.IQTREE_AMINO_ORDER)
    assert amino["I"] == pytest.approx(0.60) and amino["L"] == pytest.approx(0.40)
    assert codon == "ATT"


def test_codon_lines_keep_codons_from_one_percent() -> None:
    lines = codons.codon_lines("Node1", 2, {"CTT": 0.4, "ATT": 0.35, "ATC": 0.245, "CTC": 0.005})
    assert lines == ["Node1\t2\tCTT\t0.40000\n", "Node1\t2\tATT\t0.35000\n",
                     "Node1\t2\tATC\t0.24500\n"]  # fmt: skip
    assert codons.codon_lines("N", 1, {"ATG": 0.004}) == ["N\t1\tATG\t0.00400\n"]  # best kept
