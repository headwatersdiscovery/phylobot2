"""A job's cells, frame and ancestors (phylobot.uncertainty_run),
on a hand-made job whose numbers can be worked out on paper, and the pipeline's
phylobot-uncertainty around it."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from phylobot import uncertainty as u
from phylobot import uncertainty_run as r
from phylobot.pipeline import uncertainty as cli

MAFFT_LG, MAFFT_WAG, MUSCLE_LG = "mafft/iqtree/LG+G4", "mafft/iqtree/WAG+G4", "muscle/iqtree/LG+G4"
# Five sequences; MUSCLE splits the last column (e's V on its own), so reference column 3
# maps to MUSCLE's column 4 with agreement 4/5.
ALIGNMENTS = {
    "mafft": {"a": "MKL", "b": "MKL", "c": "MRL", "d": "MRL", "e": "MKV"},
    "muscle": {"a": "MK-L", "b": "MK-L", "c": "MR-L", "d": "MR-L", "e": "MKV-"},
}
# MAFFT's two trees agree (named differently); MUSCLE's has (a,b) but not (c,d) and
# not the same subtrees at the top.
TREES = {
    MAFFT_LG: "(e:0.1,(a:0.1,b:0.1)Node2:0.1,(c:0.1,d:0.1)Node3:0.1)Node1;",
    MAFFT_WAG: "(e:0.1,(c:0.1,d:0.1)Node7:0.1,(a:0.1,b:0.1)Node8:0.1)Node6;",
    MUSCLE_LG: "(e:0.1,((a:0.1,b:0.1)Node5:0.1,c:0.1)Node4:0.1,d:0.1)Node9;",
}
# Every node certain at every column; Node5 (MUSCLE's a,b ancestor) says R at column 2
# where MAFFT's say K.
STATES = {
    MAFFT_LG: {"Node1": "MKL", "Node2": "MKL", "Node3": "MRL"},
    MAFFT_WAG: {"Node6": "MKL", "Node7": "MRL", "Node8": "MKL"},
    MUSCLE_LG: {"Node9": "MKLL", "Node4": "MKLL", "Node5": "MRLL"},
}
BIC = {MAFFT_LG: 100.0, MAFFT_WAG: 102.0, MUSCLE_LG: 90.0}


def posteriors(states: dict[str, str]) -> str:
    rows = ["node\tsite\t" + "\t".join(u.AMINO)]
    for node, sequence in states.items():
        for site, residue in enumerate(sequence, start=1):
            values = ["1.00000" if a == residue else "0.00000" for a in u.AMINO]
            rows.append(f"{node}\t{site}\t" + "\t".join(values))
    return "\n".join(rows) + "\n"


def fasta(records: dict[str, str]) -> str:
    return "".join(f">{name}\n{seq}\n" for name, seq in records.items())


def combo(key: str, rank: int) -> dict[str, Any]:
    aligner, tree, model = key.split("/")
    return {
        "key": key, "aligner": aligner, "tree": tree, "model": model,
        "label": key, "aligner_label": aligner.upper(), "tree_label": "IQ-TREE",
        "model_label": model, "fit": {"model": model, "bic": str(BIC[key])}, "rank_bic": rank,
        "data": "amino",
    }  # fmt: skip


def job_files(masks: dict[str, dict[str, str]] | None = None) -> dict[str, str]:
    summary = {
        "primary": MAFFT_LG,
        "combinations": [combo(MAFFT_LG, 1), combo(MAFFT_WAG, 2), combo(MUSCLE_LG, 1)],
    }
    files = {"summary.json": json.dumps(summary)}
    for key, states in STATES.items():
        aligner = key.split("/")[0]
        files[f"{key}/alignment.fasta"] = fasta(ALIGNMENTS[aligner])
        files[f"{key}/tree.nwk"] = TREES[key]
        files[f"{key}/ancestors.fasta"] = fasta(states)
        files[f"{key}/ancestors_posteriors.tsv"] = posteriors(states)
        files[f"{key}/ancestors.gaps-parsimony.fasta"] = fasta((masks or {}).get(key, states))
    return files


def reader(files: dict[str, str]) -> r.Read:
    return files.get


def nodes_tsv(out: dict[str, str]) -> dict[tuple[str, str, str], dict[str, str]]:
    header, *lines = out["nodes.tsv"].splitlines()
    names = header.split("\t")
    rows = [dict(zip(names, line.split("\t"), strict=True)) for line in lines]
    return {(row["mode"], row["node"], row["part"]): row for row in rows}


def sites(out: dict[str, str], node: str) -> list[dict[str, str]]:
    """One row per mode, column and part, as the long form reads (the file is wide)."""
    index = dict(line.split("\t") for line in out["sites/index.tsv"].splitlines()[1:])
    header, *lines = out[f"sites/{index[node]}"].splitlines()
    names = header.split("\t")
    rows = []
    for line in lines:
        wide = dict(zip(names, line.split("\t"), strict=True))
        for part in r.PARTS:
            row = {"mode": wide["mode"], "column": wide["column"], "part": part}
            row["agreement"] = wide["agreement"]
            for k in r.SITE_VALUES:
                if part == "all":
                    value = float(wide[f"presence_{k}"]) + float(wide[f"residue_{k}"])
                else:
                    value = float(wide[f"{part}_{k}"])
                row[k] = str(value)
            rows.append(row)
    return rows


class TestTheRun:
    def test_the_ancestor_muscle_disagrees_about(self) -> None:
        out = r.run(reader(job_files()), r.Settings())
        rows = nodes_tsv(out)
        node2 = rows[("posterior", "Node2", "all")]
        assert (node2["present"], node2["of"]) == ("3", "3")
        # Pool (2/3 K, 1/3 R): total 4/9. By hand: v(A) = 4/9, v(M) = 1/9, v(A,M) = 4/9,
        # so Shapley gives the alignment 7/18 and the model 1/18.
        assert float(node2["total"]) == pytest.approx(4 / 9, abs=1e-4)
        assert float(node2["alignment"]) == pytest.approx(7 / 18, abs=1e-4)
        assert float(node2["model"]) == pytest.approx(1 / 18, abs=1e-4)
        assert float(node2["within"]) == pytest.approx(0, abs=1e-4)
        assert float(node2["tree"]) == 0 and float(node2["gaps"]) == 0  # not varied
        # Point estimates agree with the posteriors when every posterior is certain.
        assert float(rows[("map", "Node2", "all")]["total"]) == pytest.approx(4 / 9, abs=1e-4)

    def test_a_clade_missing_from_a_tree_is_counted_not_hidden(self) -> None:
        rows = nodes_tsv(r.run(reader(job_files()), r.Settings()))
        node3 = rows[("posterior", "Node3", "all")]
        assert (node3["present"], node3["of"]) == ("2", "3")
        assert float(node3["total"]) == pytest.approx(0, abs=1e-4)

    def test_nodes_are_the_reference_tree_in_preorder(self) -> None:
        out = r.run(reader(job_files()), r.Settings())
        assert [line.split("\t")[0] for line in out["sites/index.tsv"].splitlines()[1:]] == [
            "Node1", "Node2", "Node3",
        ]  # fmt: skip

    def test_columns_are_mapped_onto_the_reference(self) -> None:
        out = r.run(reader(job_files()), r.Settings())
        rows = [s for s in sites(out, "Node2") if s["mode"] == "posterior" and s["part"] == "all"]
        assert [s["column"] for s in rows] == ["1", "2", "3"]  # the reference's three columns
        assert float(rows[1]["total"]) == pytest.approx(4 / 9, abs=1e-4)
        assert float(rows[2]["agreement"]) == pytest.approx(0.8)  # MUSCLE splits it
        assert float(rows[2]["total"]) == 0  # all L there, once mapped

    def test_node_totals_are_the_sum_of_their_sites(self) -> None:
        out = r.run(reader(job_files()), r.Settings())
        rows = nodes_tsv(out)
        for mode in r.MODES:
            for part in r.PARTS:
                per_site = [
                    s for s in sites(out, "Node2") if (s["mode"], s["part"]) == (mode, part)
                ]
                total = sum(float(s["total"]) for s in per_site)
                assert float(rows[(mode, "Node2", part)]["total"]) == pytest.approx(total, abs=1e-3)

    def test_what_each_choice_says(self) -> None:
        out = r.run(reader(job_files()), r.Settings())
        index = dict(line.split("\t") for line in out["drivers/index.tsv"].splitlines()[1:])
        rows = [line.split("\t") for line in out[f"drivers/{index['Node2']}"].splitlines()[1:]]
        said = {(c, f, lv): st for c, f, lv, st, _p in rows}
        assert said[("2", "alignment", "mafft")] == "K"
        assert said[("2", "alignment", "muscle")] == "R"

    def test_pairs_and_the_summary(self) -> None:
        out = r.run(reader(job_files()), r.Settings())
        summary = json.loads(out["summary.json"])
        assert summary["reference"] == MAFFT_LG and summary["nodes"] == 3
        assert summary["columns"] == 3 and not summary["by_clade"]
        varied = {f["key"]: f["varied"] for f in summary["factors"]}
        assert varied == {"alignment": True, "tree": False, "model": True, "gaps": False}
        keys = [c["key"] for c in summary["cells"]]
        pairs = summary["pairs"]["posterior"]
        a, b = keys.index(f"{MAFFT_LG}@parsimony"), keys.index(f"{MUSCLE_LG}@parsimony")
        # Averaged over the nodes both have: only Node2, which differs at one column.
        assert pairs[a][b] == pytest.approx(1.0, abs=1e-3)


class TestGaps:
    def test_a_gap_against_a_residue_is_presence(self) -> None:
        masks = {MUSCLE_LG: {"Node9": "MKLL", "Node4": "MKLL", "Node5": "-RLL"}}
        out = r.run(reader(job_files(masks)), r.Settings())
        column1 = next(
            s for s in sites(out, "Node2")
            if (s["mode"], s["part"], s["column"]) == ("posterior", "presence", "1")
        )  # fmt: skip
        assert float(column1["total"]) == pytest.approx(2 * (1 / 3) * (2 / 3), abs=1e-4)
        assert float(column1["alignment"]) > 0

    def test_without_masks_the_original_reconstruction_is_used(self) -> None:
        files = {k: v for k, v in job_files().items() if "gaps-parsimony" not in k}
        summary = json.loads(r.run(reader(files), r.Settings())["summary.json"])
        assert {c["version"] for c in summary["cells"]} == {"original"}


class TestWeights:
    def test_flat(self) -> None:
        combos = r.eligible(json.loads(job_files()["summary.json"]), reader(job_files()))
        assert r.base_weights(combos, r.Settings()) == {k: 1.0 for k in STATES}

    def test_bic_compares_fits_within_an_alignment_only(self) -> None:
        combos = r.eligible(json.loads(job_files()["summary.json"]), reader(job_files()))
        w = r.base_weights(combos, r.Settings("bic"))
        # MAFFT's two by Schwarz weights (ΔBIC 2 → e^-1), each alignment then the same:
        # MUSCLE's lower BIC means nothing against MAFFT's.
        assert w[MAFFT_LG] / w[MAFFT_WAG] == pytest.approx(2.718281828, rel=1e-6)
        assert w[MAFFT_LG] + w[MAFFT_WAG] == pytest.approx(w[MUSCLE_LG])

    def test_a_zero_multiplier_leaves_a_level_out(self) -> None:
        settings = r.Settings(multipliers={"alignment": {"muscle": 0.0}})
        rows = nodes_tsv(r.run(reader(job_files()), settings))
        node2 = rows[("posterior", "Node2", "all")]
        assert node2["of"] == "2" and float(node2["total"]) == pytest.approx(0, abs=1e-4)

    def test_modelfinder_duplicates_count_once(self) -> None:
        files = job_files()
        summary = json.loads(files["summary.json"])
        finder = combo(MAFFT_LG, 1) | {"key": "mafft/iqtree/modelfinder", "model": "modelfinder"}
        finder["fit"] = {"model": "LG+G4", "bic": "100"}
        summary["combinations"].append(finder)
        for name in ("ancestors.fasta", "ancestors_posteriors.tsv"):
            files[f"mafft/iqtree/modelfinder/{name}"] = files[f"{MAFFT_LG}/{name}"]
        keys = [c["key"] for c in r.eligible(summary, reader(files))]
        assert MAFFT_LG not in keys and "mafft/iqtree/modelfinder" in keys


class TestSettings:
    def test_round_trip(self) -> None:
        s = r.Settings("bic", {"alignment": {"mafft": 2.0}}, ("fitch",), MAFFT_WAG)
        assert r.Settings.parse(json.loads(json.dumps(s.to_json()))) == s

    @pytest.mark.parametrize(
        "bad",
        [
            [],
            {"weighting": "aic"},
            {"surprise": 1},
            {"multipliers": {"colour": {"x": 1}}},
            {"multipliers": {"alignment": {"mafft": -1}}},
            {"multipliers": {"alignment": {"mafft": True}}},
            {"versions": ["../../etc"]},
            {"reference": 3},
        ],
    )
    def test_refused(self, bad: object) -> None:
        with pytest.raises(r.SettingsError):
            r.Settings.parse(bad)

    def test_an_unknown_reference_is_refused(self) -> None:
        with pytest.raises(r.SettingsError):
            r.run(reader(job_files()), r.Settings(reference="x/y/z"))


class TestTheCommand:
    def write(self, root: Path, files: dict[str, str]) -> None:
        for name, text in files.items():
            (root / name).parent.mkdir(parents=True, exist_ok=True)
            (root / name).write_text(text)

    def test_it_writes_the_files(self, tmp_path: Path) -> None:
        self.write(tmp_path / "published", job_files())
        (tmp_path / "settings.json").write_text("{}")
        args = [str(tmp_path / "published"), str(tmp_path / "settings.json"), str(tmp_path / "out")]
        assert cli.main(args) == 0
        out = tmp_path / "out" / "uncertainty"
        assert (out / "summary.json").is_file() and (out / "nodes.tsv").is_file()
        assert sorted(p.name for p in (out / "sites").iterdir()) == [
            "1.tsv",
            "2.tsv",
            "3.tsv",
            "index.tsv",
        ]

    def test_a_hostile_combination_key_is_refused(self, tmp_path: Path) -> None:
        files = job_files()
        summary = json.loads(files["summary.json"])
        summary["combinations"][0] |= {
            "key": "../../x/y",
            "aligner": "..",
            "tree": "..",
            "model": "x/y",
        }
        files["summary.json"] = json.dumps(summary)
        self.write(tmp_path / "published", files)
        (tmp_path / "settings.json").write_text("{}")
        args = [str(tmp_path / "published"), str(tmp_path / "settings.json"), str(tmp_path / "out")]
        assert cli.main(args) == 1
        assert not (tmp_path / "out").exists()
