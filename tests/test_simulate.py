"""The simulation engine (phylobot.simulate) without running AliSim.

AliSim itself runs in tests/golden/test_simulate.py (Docker); here the specs, the
planned commands and the collection of AliSim's outputs are checked on their own. The
strongest check is that `prepare` plans exactly the commands recorded in every
committed golden replicate's provenance.json.
"""

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from phylobot.seqio import read_fasta, write_fasta
from phylobot.simulate import SpecError, load, parse
from phylobot.simulate.engine import Plan, collect, neutral_names, prepare
from phylobot.simulate.spec import input_files

GOLDEN = Path(__file__).parent / "golden"
SPECS = sorted(GOLDEN.glob("*/*/spec.toml"))
IMAGE = "quay.io/biocontainers/iqtree:3.1.3--h8471819_0"


def replicates() -> list[Any]:
    return [
        pytest.param(spec.parent, seed, id=f"{spec.parent.relative_to(GOLDEN)}/{seed}")
        for spec in SPECS
        for seed in load(spec).seeds
        if (spec.parent / "replicates" / str(seed)).is_dir()
    ]


@pytest.mark.parametrize(("scenario", "seed"), replicates())
def test_prepare_plans_the_golden_commands(scenario: Path, seed: int, tmp_path: Path) -> None:
    spec = load(scenario / "spec.toml")
    plan = prepare(spec, seed, tmp_path, *input_files(spec, scenario))
    recorded = json.loads((scenario / "replicates" / str(seed) / "provenance.json").read_text())
    planned = plan.commands if spec.model.site_classes else plan.commands[0]
    assert planned == recorded["command"]


def aa_spec(**changes: Any) -> dict[str, Any]:
    spec: dict[str, Any] = {
        "tree": {"random": "yh", "taxa": 10, "lengths": [0.02, 0.1, 0.3]},
        "root": {"length": 200},
        "model": {"substitution": "LG+G4{0.75}"},
        "simulation": {"seeds": [5], "image": IMAGE},
    }
    for section, value in changes.items():
        spec[section] = value
    return spec


def test_random_tree_and_root_command(tmp_path: Path) -> None:
    spec = parse(aa_spec(indels={"rates": "0.03,0.03", "sizes": "GEO{3},GEO{3}"}))
    command = prepare(spec, 5, tmp_path).commands[0]
    assert command == [
        "iqtree3", "--alisim", "sim",
        "-t", "RANDOM{yh/10}", "-rlen", "0.02", "0.1", "0.3",
        "-m", "LG+G4{0.75}",
        "--length", "200",
        "--write-all", "--out-format", "fasta", "--seed", "5",
        "-redo",
        "--indel", "0.03,0.03", "--indel-size", "GEO{3},GEO{3}",
    ]  # fmt: skip
    assert not (tmp_path / "tree.nwk").exists()


def test_birth_death_codon_command(tmp_path: Path) -> None:
    spec = parse(
        aa_spec(
            tree={
                "random": "bd", "taxa": 8, "birth": 0.1, "death": 0.05,
                "lengths": [0.01, 0.05, 0.2],
            },
            root={"length": 50},
            model={"seqtype": "CODON1", "substitution": "GY{2.0,0.3}+F3X4"},
        )
    )  # fmt: skip
    command = prepare(spec, 5, tmp_path).commands[0]
    assert command[3:7] == ["-st", "CODON1", "-t", "RANDOM{bd{0.1/0.05}/8}"]
    assert command[command.index("--length") + 1] == "150"  # nucleotides: 50 codons


def test_plan_files(tmp_path: Path) -> None:
    plan = prepare(parse(aa_spec()), 5, tmp_path)
    assert Plan.read(tmp_path) == plan
    args = (tmp_path / "command0.args").read_text()
    assert args.split("\0")[:-1] == plan.commands[0][1:]  # NUL-terminated, iqtree3 dropped
    with pytest.raises(ValueError, match="iqtree3 only"):
        Plan(5, [["sh", "-c", "true"]]).write(tmp_path)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"model": {"substitution": "LG; rm -rf /"}}, "model.substitution"),
        ({"model": {"substitution": "LG", "seqtype": "DNA"}}, "model.seqtype"),
        ({"tree": {"random": "star", "taxa": 5, "lengths": [0.1, 0.2, 0.3]}}, "tree.random"),
        ({"tree": {"random": "yh", "taxa": 2, "lengths": [0.1, 0.2, 0.3]}}, "tree.taxa"),
        ({"tree": {"random": "yh", "taxa": 9999, "lengths": [0.1, 0.2, 0.3]}}, "at most"),
        ({"tree": {"random": "yh", "taxa": 5, "lengths": [0.1, 0.1, 0.1]}}, "min below max"),
        ({"tree": {"random": "bd", "taxa": 5, "birth": 0.1, "death": 0.2,
                   "lengths": [0.1, 0.2, 0.3]}}, "death"),
        ({"tree": {"file": "t.nwk", "outgroup": "A", "midpoint": True}}, "not both"),
        ({"root": {"length": 0}}, "root.length"),
        ({"root": {"file": "r.fasta", "name": "a,b"}}, "root.name"),
        ({"indels": {"rates": "0.03", "sizes": "GEO{3}; ls"}}, "indels.sizes"),
        ({"indels": {"rates": "-1,0.03", "sizes": "GEO{3}"}}, "indels.rates"),
        ({"simulation": {"seeds": []}}, "simulation.seeds"),
        ({"simulation": {"seeds": [1, 1]}}, "distinct"),
        ({"simulation": {"seeds": [1], "image": "iqtree:latest"}}, "pinned"),
        ({"simulation": {"seeds": [1], "names": "real"}}, "simulation.names"),
    ],
)  # fmt: skip
def test_bad_specs_are_refused(change: dict[str, Any], message: str) -> None:
    with pytest.raises(SpecError, match=message):
        parse(aa_spec(**change))


def test_site_classes_need_a_root_file_and_no_indels() -> None:
    m3 = load(GOLDEN / "lysin" / "codon-selection" / "spec.toml")
    assert m3.model.site_classes is not None
    data = {
        "tree": {"file": "t.nwk", "midpoint": True},
        "root": {"length": 100},
        "model": {
            "seqtype": "CODON1", "kappa": 1.5, "f3x4": [0.25] * 12,
            "site_classes": {"proportions": [0.5, 0.5], "omegas": [0.1, 2.0]},
        },
        "simulation": {"seeds": [1]},
    }  # fmt: skip
    with pytest.raises(SpecError, match="root sequence file"):
        parse(data)
    data["root"] = {"file": "r.fasta", "name": "r"}
    data["indels"] = {"rates": "0.03,0.03", "sizes": "GEO{2},GEO{2}"}
    with pytest.raises(SpecError, match="indels are not supported"):
        parse(data)


def fake_alisim(work: Path) -> None:
    """What AliSim writes for a three-leaf tree with --write-all and indels."""
    write_fasta(work / "sim.fa", {"A": "MK-L", "B": "M-QL", "C": "MKQL", "4": "MKQL"})
    write_fasta(work / "sim.unaligned.fa", {"A": "MKL", "B": "MQL", "C": "MKQL"})
    (work / "sim.full.treefile").write_text("(A:0.1,B:0.2,C:0.3)4:0;\n")


@pytest.mark.parametrize("names", ["keep", "neutral"])
def test_collect(names: str, tmp_path: Path) -> None:
    spec = parse(
        aa_spec(
            indels={"rates": "0.03,0.03", "sizes": "GEO{3},GEO{3}"},
            simulation={"seeds": [5], "names": names},
        )
    )
    work, out = tmp_path / "work", tmp_path / "out"
    work.mkdir()
    plan = prepare(spec, 5, work)
    fake_alisim(work)
    collect(spec, plan, work, out, "IQ-TREE version 3.1.3", IMAGE)

    rename = {"A": "sim001", "B": "sim002", "C": "sim003"} if names == "neutral" else {}

    def leaf(name: str) -> str:
        return rename.get(name, name)

    expected = {leaf("A"): "MKL", leaf("B"): "MQL", leaf("C"): "MKQL"}
    assert read_fasta(out / "input.fasta") == expected
    assert list(read_fasta(out / "true_alignment.fasta")) == [leaf(n) for n in "ABC"]
    assert read_fasta(out / "true_ancestors.fasta") == {"4": "MKQL"}
    tree = f"({leaf('A')}:0.1,{leaf('B')}:0.2,{leaf('C')}:0.3)4:0;\n"
    assert (out / "true_tree.nwk").read_text() == tree
    provenance = json.loads((out / "provenance.json").read_text())
    assert provenance["command"] == plan.commands[0]
    assert ("names.tsv" in provenance["sha256"]) == (names == "neutral")
    if names == "neutral":
        rows = "".join(f"sim00{i}\t{n}\n" for i, n in enumerate("ABC", 1))
        assert (out / "names.tsv").read_text() == "simulated\toriginal\n" + rows


def test_collect_keeps_hyphens_in_names(tmp_path: Path) -> None:
    """AliSim's sim.unaligned.fa loses them (MAGUK's ArabGUK1-), so it is not read."""
    spec = parse(aa_spec(indels={"rates": "0.03,0.03", "sizes": "GEO{3},GEO{3}"},
                         simulation={"seeds": [5], "names": "keep"}))  # fmt: skip
    work, out = tmp_path / "work", tmp_path / "out"
    work.mkdir()
    plan = prepare(spec, 5, work)
    write_fasta(work / "sim.fa", {"A-": "MK-L", "B-2": "M-QL", "C": "MKQL", "4": "MKQL"})
    write_fasta(work / "sim.unaligned.fa", {"A": "MKL", "B2": "MQL", "C": "MKQL"})
    (work / "sim.full.treefile").write_text("(A-:0.1,B-2:0.2,C:0.3)4:0;\n")
    collect(spec, plan, work, out, "IQ-TREE version 3.1.3", IMAGE)
    assert read_fasta(out / "input.fasta") == {"A-": "MKL", "B-2": "MQL", "C": "MKQL"}


def test_neutral_names_widen() -> None:
    names = neutral_names([f"T{i}" for i in range(1, 1001)])
    assert names["T1"] == "sim0001" and names["T1000"] == "sim1000"


def test_cli_refuses_a_seed_not_in_the_spec(tmp_path: Path) -> None:
    (tmp_path / "spec.json").write_text(json.dumps(aa_spec()))
    command = [sys.executable, "-m", "phylobot.simulate.cli", "prepare"]
    result = subprocess.run(
        [*command, str(tmp_path / "spec.json"), "6", str(tmp_path / "w")],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2 and "not one of the spec's seeds" in result.stderr
    ok = subprocess.run([*command, str(tmp_path / "spec.json"), "5", str(tmp_path / "w")])
    assert ok.returncode == 0 and (tmp_path / "w" / "plan.json").exists()
