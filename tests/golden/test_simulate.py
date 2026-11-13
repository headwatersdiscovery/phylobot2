"""The simulation engine still produces the committed golden datasets, byte for byte.

The golden replicates are frozen truth, generated before the engine moved from
tests/golden/generate.py into phylobot.simulate. Both ways of running it must reproduce
them exactly: on the host (generate.py, AliSim in Docker) and through the pipeline
(--simulate: SIM_PREPARE -> SIM_ALISIM -> SIM_COLLECT). If a test here fails, the
engine changed what it simulates: stop and report, never regenerate the goldens.
Gene families (Sim2) have no golden dataset yet; their test checks a grown family's
bundle for consistency and that the engine and the pipeline produce the same one.

Slow and needs Docker, so excluded from `make test`; run with `make golden`.
"""

import filecmp
import hashlib
import json
import tomllib
from pathlib import Path
from typing import Any

import dendropy
import pytest

from phylobot.seqio import read_fasta
from phylobot.simulate import families, load, parse
from phylobot.simulate.engine import simulate_locally

from .generate import generate
from .pipeline import run_pipeline

GOLDEN = Path(__file__).resolve().parent
ROOT = GOLDEN.parent.parent
IMAGE = "quay.io/biocontainers/iqtree:3.1.3--h8471819_0"


def replicates() -> list[Any]:
    """Every committed replicate; all but each scenario's first seed in the full tier."""
    params = []
    for spec in sorted(GOLDEN.glob("*/*/spec.toml")):
        for i, seed in enumerate(load(spec).seeds):
            if not (spec.parent / "replicates" / str(seed)).is_dir():
                continue
            marks = [pytest.mark.golden] + ([pytest.mark.golden_full] if i > 0 else [])
            name = f"{spec.parent.relative_to(GOLDEN)}/{seed}"
            params.append(pytest.param(spec.parent, seed, marks=marks, id=name))
    return params


def assert_same(produced: Path, committed: Path) -> None:
    files = sorted(p.name for p in committed.iterdir())
    assert sorted(p.name for p in produced.iterdir()) == files
    differ = [f for f in files if not filecmp.cmp(produced / f, committed / f, shallow=False)]
    assert not differ, f"{committed.relative_to(GOLDEN)}: {differ} changed"


@pytest.mark.parametrize(("scenario", "seed"), replicates())
def test_engine_reproduces_golden(scenario: Path, seed: int, tmp_path: Path) -> None:
    generate(scenario, seed, tmp_path / "out")
    assert_same(tmp_path / "out", scenario / "replicates" / str(seed))


@pytest.mark.golden
@pytest.mark.parametrize("scenario", ["ime2/indels", "lysin/codon-selection"])
def test_pipeline_reproduces_golden(scenario: str) -> None:
    """--simulate on a plain and a site-class (codeml M3) scenario's first seed."""
    directory = GOLDEN / scenario
    data = tomllib.loads((directory / "spec.toml").read_text())
    seed = data["simulation"]["seeds"][0]
    data["simulation"]["seeds"] = [seed]
    tree = directory / data["tree"].pop("file")
    root = directory / data["root"].pop("file")
    out = ROOT / "results" / "golden" / "simulate" / scenario
    out.mkdir(parents=True, exist_ok=True)
    (out / "spec.json").write_text(json.dumps(data))
    params = {
        "simulate": str(out / "spec.json"),
        "sim_tree": str(tree.resolve()),
        "sim_root": str(root.resolve()),
    }
    run = run_pipeline(out, None, params)
    assert_same(run / "replicates" / str(seed), directory / "replicates" / str(seed))


FAMILY_SPECIES = (
    "(((A:0.1,B:0.12):0.08,(C:0.1,D:0.15):0.05):0.1,((E:0.2,F:0.1):0.1,(G:0.1,H:0.1):0.2):0.1);"
)
FAMILY_SPEC = {
    "description": "Sim2: duplication and loss with hidden paralogy, then AliSim with indels",
    "tree": {},
    "root": {"length": 120},
    "model": {"substitution": "LG+G4{0.75}"},
    "indels": {"rates": "0.03,0.03", "sizes": "GEO{3},GEO{3}"},
    "family": {
        "duplication": 0.5,
        "loss": 0.5,
        "origin": 0.05,
        "duplications": [{"clade": list("ABCDEFGH")}],
        "losses": [{"clade": list("ABCD"), "copy": "1b"}, {"clade": list("EF"), "copy": "1a"}],
    },
    "simulation": {"seeds": [3], "names": "neutral", "image": IMAGE},
}


def check_family(bundle: Path) -> None:
    """The truth bundle of a grown family agrees with itself."""
    tree = dendropy.Tree.get(
        path=bundle / "true_tree.nwk", schema="newick", preserve_underscores=True
    )
    internal = [n.label for n in tree.preorder_node_iter() if n.taxon is None]
    assert list(read_fasta(bundle / "true_ancestors.fasta")) == internal  # AliSim kept labels
    rows = [r.split("\t") for r in (bundle / families.EVENTS).read_text().splitlines()[1:]]
    visible = [r[1] for r in rows if r[1] != "-"]
    assert sorted(visible) == sorted(internal)  # each ancestor is a duplication or speciation
    assert any(r[0] == "duplication" and r[4] == "yes" and r[1] != "-" for r in rows)
    genes = [r.split("\t")[0] for r in (bundle / families.SPECIES).read_text().splitlines()[1:]]
    assert genes == list(read_fasta(bundle / "input.fasta"))
    provenance = json.loads((bundle / "provenance.json").read_text())
    for file, digest in provenance["sha256"].items():
        assert hashlib.sha256((bundle / file).read_bytes()).hexdigest() == digest


@pytest.mark.golden
def test_family_through_engine_and_pipeline(tmp_path: Path) -> None:
    """Sim2: a family grown on a species tree, then AliSim, the same both ways."""
    out = ROOT / "results" / "golden" / "simulate" / "family"
    out.mkdir(parents=True, exist_ok=True)
    (out / "species.nwk").write_text(FAMILY_SPECIES + "\n")
    (out / "spec.json").write_text(json.dumps(FAMILY_SPEC))

    simulate_locally(parse(FAMILY_SPEC), 3, tmp_path / "engine", out / "species.nwk")
    check_family(tmp_path / "engine")
    params = {"simulate": str(out / "spec.json"), "sim_tree": str(out / "species.nwk")}
    run = run_pipeline(out, None, params)
    assert_same(run / "replicates" / "3", tmp_path / "engine")
