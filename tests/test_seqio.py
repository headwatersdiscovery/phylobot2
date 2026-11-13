from pathlib import Path

import pytest

from phylobot.seqio import read_fasta, write_fasta

FIXTURES = Path(__file__).parent / "fixtures" / "fasta"


def test_name_edge_cases_are_read_verbatim() -> None:
    records = read_fasta(FIXTURES / "name_edge_cases.fasta")
    assert list(records) == ["species_1", "species.2", "species__3", "spe_cies.4", "1_spe_cies.4"]
    # Two names, one sequence: duplicate detection will need to catch this.
    assert records["spe_cies.4"] == records["1_spe_cies.4"]


def test_duplicate_names_rejected(tmp_path: Path) -> None:
    path = tmp_path / "dup.fasta"
    path.write_text(">a\nMK\n>a\nMT\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_fasta(path)


def test_round_trip(tmp_path: Path) -> None:
    records = {"x": "MKT", "y": "M-T"}
    write_fasta(tmp_path / "out.fasta", records)
    assert read_fasta(tmp_path / "out.fasta") == records
