"""Per-site calibration records: the viewer's covariates, ESR tips' own."""

import csv
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Literal

import pytest

from phylobot import calibration
from phylobot.evaluation import records

TREE = "(a:0.1,b:0.2,(c:0.3,d:0.4)Node2:0.05)Node1;\n"
ALIGNMENT = {"a": "MK-", "b": "MKA", "c": "MRA", "d": "LRA"}


def posterior(**probs: float) -> dict[str, float]:
    rest = (1 - sum(probs.values())) / (20 - len(probs))
    return {a: probs.get(a, rest) for a in records.AMINO_ACIDS}


def write_tsv(path: Path, fields: list[str], rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


@pytest.fixture
def folder(tmp_path: Path) -> Path:
    (tmp_path / "alignment.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in ALIGNMENT.items()))
    (tmp_path / "tree.nwk").write_text(TREE)
    (tmp_path / "site_rates.tsv").write_text("site\trate\n1\t0.5\n2\t1.0\n3\t2.0\n")
    rows = []
    for node in ("Node1", "Node2"):
        for site, probs in enumerate(
            [posterior(M=0.9), posterior(K=0.6, R=0.3), posterior(A=0.5)], 1
        ):
            rows.append({"node": node, "site": site, **probs})
    write_tsv(tmp_path / "ancestors_posteriors.tsv", ["node", "site", *records.AMINO_ACIDS], rows)
    (tmp_path / "ancestors.gaps-parsimony.fasta").write_text(">Node1\nMK-\n>Node2\nMRA\n")
    return tmp_path


def test_ancestor_records_use_the_viewers_covariates(folder: Path) -> None:
    rows = list(records.ancestor_records(folder, "parsimony"))
    assert [(r["node"], r["site"]) for r in rows] == [
        ("Node1", 1), ("Node1", 2), ("Node1", 3), ("Node2", 1), ("Node2", 2), ("Node2", 3),
    ]  # fmt: skip
    assert rows[1]["ml"] == "K" and rows[1]["pp_ml"] == pytest.approx(0.6)
    assert rows[0]["true"] == "" and rows[0]["correct"] == ""
    assert [r["present"] for r in rows] == ["1", "1", "0", "1", "1", "1"]
    expected = calibration.covariates(
        list(ALIGNMENT.values()), TREE, "Node2", (folder / "site_rates.tsv").read_text()
    )
    assert expected is not None
    for row, x in zip(rows[3:], expected, strict=True):
        assert {k: row[k] for k in calibration.KEYS} == x


def test_without_a_gap_version(folder: Path) -> None:
    assert {r["present"] for r in records.ancestor_records(folder)} == {""}


def test_a_shifted_clade_is_scored_against_its_own_columns(tmp_path: Path) -> None:
    """The aligner shifted {a, b} one column against {c, d, e}: Node2, a's and b's
    ancestor, is right about R, but the whole alignment's map takes its column for true
    column 3 (L there). Its own descendants say true column 2."""
    folder = tmp_path / "combination"
    folder.mkdir()
    shifted = {"a": "-GRLIT", "b": "-GRLIT", "c": "GKLVT-", "d": "GKLVT-", "e": "GKLVT-"}
    shifted = {"o": "GKLVT-", **shifted}
    (folder / "alignment.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in shifted.items()))
    (folder / "tree.nwk").write_text("(o,((a:1,b:1)Node2:1,(c:1,d:1,e:1)Node3:1)Node1:1);\n")
    (folder / "site_rates.tsv").write_text(
        "site\trate\n" + "".join(f"{i}\t1.0\n" for i in range(1, 7))
    )
    rows = [{"node": node, "site": 3, **posterior(**{state: 0.99})}
            for node, state in (("Node1", "L"), ("Node2", "R"), ("Node3", "L"))]  # fmt: skip
    write_tsv(folder / "ancestors_posteriors.tsv", ["node", "site", *records.AMINO_ACIDS], rows)
    truth = tmp_path / "truth"
    truth.mkdir()
    (truth / "true_tree.nwk").write_text("(o,((a:1,b:1)T2:1,(c:1,d:1,e:1)T3:1)T1:1);\n")
    aligned = {"o": "GKLVT", "a": "GRLIT", "b": "GRLIT", "c": "GKLVT", "d": "GKLVT",
               "e": "GKLVT"}  # fmt: skip
    (truth / "true_alignment.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in aligned.items()))
    (truth / "true_ancestors.fasta").write_text(">T1\nGKLVT\n>T2\nGRLIT\n>T3\nGKLVT\n")

    def calls(columns: Literal["global", "clade"]) -> dict[str, tuple[str, int]]:
        found = records.truth_records(folder, truth, "o", columns=columns)
        return {str(r["node"]): (str(r["true"]), int(str(r["correct"]))) for r in found}

    assert calls("global")["Node2"] == ("L", 0)
    assert calls("clade")["Node2"] == ("R", 1)
    # Nodes whose descendants follow the majority score the same either way.
    assert calls("global")["Node3"] == calls("clade")["Node3"] == ("L", 1)


def test_esr_records(folder: Path, tmp_path: Path) -> None:
    sites = tmp_path / "esr_sites.tsv"
    fields = ["tip", "level", "site", "true", "ml", "pp_ml", "pp_true", *records.AMINO_ACIDS]
    first = {"tip": "c", "level": "e0", "site": 1, "true": "M", "ml": "M", "pp_ml": 0.9}
    second = {"tip": "c", "level": "e1", "site": 2, "true": "R", "ml": "K", "pp_ml": 0.6}
    write_tsv(sites, fields, [
        {**first, "pp_true": 0.9, **posterior(M=0.9)},
        {**second, "pp_true": 0.3, **posterior(K=0.6, R=0.3)},
    ])  # fmt: skip
    rows = list(records.esr_records(folder, sites))
    assert [(r["source"], r["correct"]) for r in rows] == [("esr-e0", 1), ("esr-e1", 0)]
    # c's own branch (0.3), and the nearest other sequence: a, through Node2 and Node1.
    assert rows[0]["branch_length"] == pytest.approx(math.log10(0.3 / 0.1))
    assert rows[0]["tip_distance"] == pytest.approx(math.log10(0.45 / 0.1))
    # Column statistics leave c out: column 2 is K, K, R over a, b, d.
    entropy, _ = calibration.column_statistics(["MK-", "MKA", "LRA"])
    assert rows[1]["entropy"] == pytest.approx(entropy[1])


def test_tip_place() -> None:
    place = calibration.tip_place(TREE, "a")
    assert place is not None
    assert place.branch_length == 0.1
    assert place.tip_distance == pytest.approx(0.3)  # to b, through Node1
    assert calibration.tip_place(TREE, "zz") is None


def test_write(folder: Path, tmp_path: Path) -> None:
    out = tmp_path / "records.tsv"
    assert records.write(out, records.ancestor_records(folder)) == 6
    header = out.read_text().splitlines()[0].split("\t")
    assert header == records.FIELDS


def test_truth_records(folder: Path, tmp_path: Path) -> None:
    truth = tmp_path / "truth"
    truth.mkdir()
    # Same topology, other names; the true alignment has an extra column (filled only in
    # the ancestors), so inferred column 3 maps to true column 4.
    (truth / "true_tree.nwk").write_text("(a:0.1,b:0.2,(c:0.3,d:0.4)T2:0.05)T1;\n")
    (truth / "true_alignment.fasta").write_text(">a\nMK--\n>b\nMK-A\n>c\nMR-A\n>d\nLR-A\n")
    (truth / "true_ancestors.fasta").write_text(">T1\nMKWA\n>T2\nLRWA\n")
    rows = list(records.truth_records(folder, truth, outgroup="a"))
    # Rooted on a, Node1 is the clade {b, c, d} (= T1) and Node2 is {c, d} (= T2).
    assert [(r["node"], r["site"], r["true"], r["correct"]) for r in rows] == [
        ("Node1", 1, "M", 1), ("Node1", 2, "K", 1), ("Node1", 3, "A", 1),
        ("Node2", 1, "L", 0), ("Node2", 2, "R", 0), ("Node2", 3, "A", 1),
    ]  # fmt: skip
