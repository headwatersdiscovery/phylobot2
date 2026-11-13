"""Several sequence files combined, and alignments restricted to the sequences they
share (phylobot.inputs)."""

from phylobot import inputs

A = b">a\nMKTAY\n>b\nMKTAF\n>c\nMKSAY\n"
B = b">c desc\nMK-SAY\n>d\nMRTAY\n>b\nMKTAW\n"  # c again (identical, gaps aside); b differs


def test_identical_sequences_are_used_once_and_clashes_wait_for_a_choice() -> None:
    merged = inputs.merge_sequences([("A.fasta", A), ("B.fasta", B)])
    assert [(s.name, s.sequences, s.min_length, s.max_length) for s in merged.files] == [
        ("A.fasta", 3, 5, 5), ("B.fasta", 3, 5, 5),
    ]  # fmt: skip
    assert merged.identical == [("c", "B.fasta")]
    assert merged.clashes == [inputs.Clash("b", "A.fasta", "B.fasta", "b.2")]
    assert merged.unresolved and not merged.error
    assert merged.fasta.split()[::2] == [b">a", b">b", b">c", b">d"]  # first copy kept so far


def test_clashes_renamed_or_left_to_the_first_file() -> None:
    renamed = inputs.merge_sequences([("A.fasta", A), ("B.fasta", B)], inputs.RENAME)
    assert not renamed.unresolved
    assert b">b.2\nMKTAW\n" in renamed.fasta and renamed.fasta.count(b">") == 5
    first = inputs.merge_sequences([("A.fasta", A), ("B.fasta", B)], inputs.FIRST)
    assert b"MKTAW" not in first.fasta and first.fasta.count(b">") == 4
    # A name already taken moves the suffix along; a third file gets the next one.
    taken = b">b.2\nMMMMM\n"
    third = b">b\nMKTAV\n"
    merged = inputs.merge_sequences([("A", A), ("B", B), ("C", taken + third)], inputs.RENAME)
    assert [c.renamed for c in merged.clashes] == ["b.3", "b.4"]


def test_one_type_of_sequence_per_job_and_unreadable_files() -> None:
    dna = b">x\nATGAAATGG\n>y\nATGAAGTGG\n"
    merged = inputs.merge_sequences([("prot.fasta", A), ("cds.fasta", dna)])
    assert merged.error == (
        "Incompatible files: cds.fasta is coding DNA but prot.fasta is protein. "
        "A job takes one type of sequence: untick or remove one"
    )
    assert merged.incompatible
    # ...even when one of them has another problem.
    bad = b"MKTAY\n>x\nMKTAF\n"
    both = inputs.merge_sequences([("cds.fasta", dna), ("prot.fasta", bad)])
    assert both.incompatible and both.error.startswith("Incompatible files")
    assert both.files[1].kind == "protein" and "before the first" in both.files[1].error
    assert inputs.merge_sequences([("bad.fasta", b"MKT\n")]).error.startswith("bad.fasta: line 1")
    assert inputs.merge_sequences([("bad.fasta", b"\xff\xfe")]).error == "bad.fasta: not UTF-8 text"


def test_a_name_repeated_in_one_file_is_resolved_like_one_across_files() -> None:
    text = b">a\nMKTAY\n>b\nMKTAF\n>a\nMK-TAY\n>b\nMKTAW\n"
    merged = inputs.merge_sequences([(inputs.PASTED, text)])
    assert merged.identical == [("a", inputs.PASTED)]
    assert merged.clashes == [inputs.Clash("b", inputs.PASTED, inputs.PASTED, "b.2")]
    renamed = inputs.merge_sequences([(inputs.PASTED, text)], inputs.RENAME)
    assert renamed.fasta == b">a\nMKTAY\n>b\nMKTAF\n>b.2\nMKTAW\n"
    aligned = b">a\nMK-TAY\n>b\nMKTTAF\n>b\nMKTTAW\n"
    deduped = inputs.dedupe_alignment("aln", aligned, inputs.FIRST)
    assert deduped.fasta == b">a\nMK-TAY\n>b\nMKTTAF\n" and not deduped.unresolved
    assert inputs.dedupe_alignment("aln", aligned).unresolved


def test_file_summaries() -> None:
    stats, _ = inputs.read_stats("cds.fasta", b">x\nATGAAATGG\n>y\nATGAAG\n")
    assert stats.summary == "2 coding DNA sequences, 2–3 codons"
    aligned, _ = inputs.read_stats("aln.fasta", b">x\nMK-T\n>y\nMKAT\n", aligned=True)
    assert aligned.summary == "2 protein sequences × 4 columns" and not aligned.error
    ragged, _ = inputs.read_stats("aln.fasta", b">x\nMK-T\n>y\nMKA\n", aligned=True)
    assert ragged.error == "rows have different lengths: not an alignment"


def test_restrict_to_shared_names() -> None:
    first, second = {"a", "b", "c"}, {"b", "c", "d"}
    assert inputs.shared_names([first, second]) == {"b", "c"}
    aligned = b">a\nMK-T\n>b\nM--T\n>c\nMA-T\n"
    assert inputs.restrict(aligned, {"b", "c"}, aligned=True) == b">b\nM-T\n>c\nMAT\n"
    assert inputs.restrict(b">a x\nMK\n>b\nMT\n", {"a"}, aligned=False) == b">a x\nMK\n"
