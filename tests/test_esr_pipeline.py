"""The ESR job's pipeline steps (phylobot.pipeline.esr): the input checks the web form
shares, renaming to safe IDs, the settings, and collecting tips into the job's files.
Fitting and reconstructing need IQ-TREE: tests/golden/test_esr.py."""

import csv
import json
from pathlib import Path

import pytest

from phylobot import esr
from phylobot.pipeline import esr as step
from phylobot.trees import load_tree

ALIGNMENT = {"a": "MKLV-", "b": "MKLAA", "c": "MRLWA", "d": "MRLW-"}
TREE = "(a:0.1,b:0.2,(c:0.3,d:0.4)Node3:0.3)Node1;"


def leaves(newick: str) -> set[str]:
    return {leaf.taxon.label for leaf in load_tree(newick).leaf_node_iter()}


class TestCheck:
    def test_good(self) -> None:
        step.check(ALIGNMENT, TREE, "LG+F+I+R4")

    @pytest.mark.parametrize(
        ("alignment", "tree", "model", "message"),
        [
            ({"a": "MK", "b": "MK", "c": "MK"}, "(a,b,c);", "LG", "at least 4"),
            (ALIGNMENT | {"d": "MRL"}, TREE, "LG", "not aligned"),
            ({k: "ACGTA" for k in ALIGNMENT}, TREE, "LG", "protein alignments only"),
            (ALIGNMENT | {"d": "MRL#-"}, TREE, "LG", "not amino acids: #"),
            (ALIGNMENT, TREE, "LG+C20+G4", "cannot fix"),
            (ALIGNMENT, TREE, "LG+G4; rm -rf /", "cannot fix"),
            (ALIGNMENT, TREE, "GY+F3X4", "cannot fix"),
            (ALIGNMENT, "(a,b,(c,e));", "LG", "not in the tree: d; not in the alignment: e"),
            (ALIGNMENT, "(a,b,(c,d,a));", "LG", "names a sequence twice"),
        ],
    )
    def test_refused(self, alignment: dict[str, str], tree: str, model: str, message: str) -> None:
        with pytest.raises(ValueError, match=message):
            step.check(alignment, tree, model)

    def test_reserved_names(self) -> None:
        alignment = ALIGNMENT | {esr.DUMMIES[0]: "MKLV-"}
        with pytest.raises(ValueError, match="reserved"):
            step.check(alignment, f"(a,b,(c,d),{esr.DUMMIES[0]});", "LG")


def test_binary() -> None:
    assert step.binary(TREE)  # IQ-TREE's trifurcating top node
    assert step.binary("((a,b),(c,d));")
    assert not step.binary("(a,b,(c,d,e)x)y;")
    assert not step.binary("(a,b,c,d);")


def test_rename_keeps_no_name() -> None:
    alignment, tree, names = step.rename({"Homo sapiens": "MK", "it's": "MR", "c": "MW"},
                                         "('Homo sapiens':1,'it''s':2,c:3);")  # fmt: skip
    assert alignment == {"t1": "MK", "t2": "MR", "t3": "MW"}
    assert leaves(tree) == {"t1", "t2", "t3"} and "Node" not in tree
    assert names == {"t1": "Homo sapiens", "t2": "it's", "t3": "c"}


def test_supported() -> None:
    for model in ("LG", "LG+G4", "WAG+I+G4", "Q.pfam+F+R4", "LG+F+I+R6", "JTT+FO+G"):
        assert esr.supported(model), model
    for model in ("", "LG+C20", "GY+F3X4+G4", "LG+H4"):
        assert not esr.supported(model), model


def staged_settings(tmp_path: Path, **values: object) -> Path:
    (tmp_path / "settings.json").write_text(json.dumps(values))
    return tmp_path


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ({"levels": ["e2"], "tips": ["a"]}, "levels"),
        ({"levels": [], "tips": ["a"]}, "levels"),
        ({"levels": ["e0"], "tips": []}, "tips"),
        ({"levels": ["e0"], "tips": [1]}, "tips"),
        ({"levels": ["e0"], "tips": ["a"], "seed": "1; rm"}, "seed"),
        ({"levels": ["e0"], "tips": ["a"], "seed": True}, "seed"),
    ],
)
def test_settings_are_checked(tmp_path: Path, values: dict[str, object], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        step._settings(staged_settings(tmp_path, **values))


def test_settings_keep_the_levels_in_order(tmp_path: Path) -> None:
    got = step._settings(staged_settings(tmp_path, levels=["e1", "e0"], tips=["a"]))
    assert got == {"levels": ["e0", "e1"], "tips": ["a"], "seed": 1, "fit_lengths": False}


def test_a_tip_id_is_checked(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="not a tip ID"):
        step.tip(tmp_path, "../x", tmp_path / "out")


def write_tip(
    folder: Path, tip: str, level: str, sites: list[tuple[str, str, float, float]]
) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    records = [{"tip": tip, "level": level, "site": i, "true": t, "ml": m, "pp_ml": p,
                "pp_true": pt, **dict.fromkeys(esr.AMINO_ACIDS, 0.05)}
               for i, (t, m, p, pt) in enumerate(sites, start=1)]  # fmt: skip
    result = esr.TipResult(tip, level, 0.1, -10.0, 0.0, records)
    for name, fields, rows in (("esr_tips.tsv", esr.TIP_FIELDS, [result.summary()]),
                               ("esr_sites.tsv", esr.SITE_FIELDS, records)):  # fmt: skip
        path = folder / name
        old = list(csv.DictReader(path.open(), delimiter="\t")) if path.exists() else []
        esr.write_tsv(path, fields, [*old, *rows])


def test_collect(tmp_path: Path) -> None:
    fit = tmp_path / "fit"
    fit.mkdir()
    (fit / "tips.txt").write_text("t2\nt1\n")
    (fit / "fit.json").write_text(json.dumps({"model": "LG", "fixed": "LG", "log_likelihood": -5}))
    (fit / "settings.json").write_text(json.dumps({"levels": ["e0", "e1"], "seed": 1}))
    (fit / "tree.nwk").write_text("(t1:1,t2:1,t3:1);\n")
    (fit / "names.tsv").write_text("id\tname\nt1\tHomo sapiens\nt2\tb\nt3\tc\n")
    write_tip(tmp_path / "tip_t1", "t1", "e0", [("M", "M", 0.95, 0.95), ("K", "R", 0.55, 0.4)])
    write_tip(tmp_path / "tip_t1", "t1", "e1", [("M", "M", 0.9, 0.9), ("K", "K", 0.5, 0.5)])
    write_tip(tmp_path / "tip_t2", "t2", "e0", [("L", "L", 0.85, 0.85)])
    out = tmp_path / "esr"
    step.collect(fit, out, [tmp_path / "tip_t1", tmp_path / "tip_t2"])
    tips = list(csv.DictReader((out / "esr_tips.tsv").open(), delimiter="\t"))
    assert [(r["tip"], r["level"]) for r in tips] == [("b", "e0"), ("Homo sapiens", "e0"),
                                                      ("Homo sapiens", "e1")]  # fmt: skip
    sites = list(csv.DictReader((out / "esr_sites.tsv").open(), delimiter="\t"))
    assert [r["tip"] for r in sites] == ["b", "Homo sapiens", "Homo sapiens", "Homo sapiens",
                                         "Homo sapiens"]  # fmt: skip
    summary = json.loads((out / "esr_summary.json").read_text())
    e0 = summary["levels"]["e0"]
    assert (e0["tips"], e0["sites"]) == (2, 3)
    assert e0["accuracy"] == pytest.approx(2 / 3) and e0["mean_pp"] == pytest.approx(2.35 / 3)
    bins = {b["low"]: b for b in e0["bins"] if b["sites"]}
    assert set(bins) == {0.5, 0.8, 0.9} and bins[0.5]["accuracy"] == 0
    assert summary["levels"]["e1"]["accuracy"] == 1 and summary["fit"]["log_likelihood"] == -5
    assert leaves((out / "tree.nwk").read_text()) == {"Homo sapiens", "b", "c"}


def test_a_posterior_of_one_lands_in_the_top_bin() -> None:
    sites = [{"level": "e0", "true": "M", "ml": "M", "pp_ml": "1.0", "pp_true": "1.0"}]
    bins = step.summarise(sites, [{"level": "e0"}])["levels"]["e0"]["bins"]
    assert bins[-1]["sites"] == 1
