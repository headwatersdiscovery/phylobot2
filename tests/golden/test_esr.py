"""ESR against the pinned IQ-TREE, on the Ime2 golden replicate.

Three checks from the first ESR experiments:
1. The daughters carry no information: the held-out tree's lnL is the pruned tree's
   plus the daughters' own constant, 2 ln pi(L) under LG+G4, to IQ-TREE's precision.
2. Simulated under the inference model, each tip's accuracy tracks its mean PP.
3. E0 reproduces the spike's numbers, and E1 stays close to E0 (the model is right and
   the tip is one of seven).

Slow and needs Docker (about a minute), so excluded from `make test`; run with
`make golden`. Ime2 is a smoke test: nothing here is a calibration result.
"""

import csv
import json
import math
from pathlib import Path

import pytest

from phylobot import esr
from phylobot.seqio import read_fasta, write_fasta

GOLDEN = Path(__file__).resolve().parent
REPLICATE = GOLDEN / "ime2/easy/replicates/101"
IMAGE = "quay.io/biocontainers/iqtree:3.1.3--h8471819_0"
PI_LEU_LG = 0.099081

pytestmark = pytest.mark.golden


@pytest.fixture(scope="module")
def combination(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A published-looking combination: the replicate searched under LG+G4."""
    folder = tmp_path_factory.mktemp("combination")
    write_fasta(folder / "alignment.fasta", read_fasta(REPLICATE / "input.fasta"))
    esr.docker_runner(IMAGE)(folder, ["-s", "alignment.fasta", "-m", "LG+G4", "--prefix", "search"])
    (folder / "tree.nwk").write_text((folder / "search.treefile").read_text())
    (folder / "model.txt").write_text("LG+G4\n")
    return folder


@pytest.fixture(scope="module")
def rows(
    combination: Path, tmp_path_factory: pytest.TempPathFactory
) -> dict[tuple[str, str], dict[str, str]]:
    out = tmp_path_factory.mktemp("esr")
    esr.run_combination(combination, out, esr.docker_runner(IMAGE))
    with (out / "esr_tips.tsv").open() as f:
        return {(r["tip"], r["level"]): r for r in csv.DictReader(f, delimiter="\t")}


def test_daughters_carry_no_information(combination: Path, tmp_path: Path) -> None:
    run = esr.docker_runner(IMAGE)
    alignment = read_fasta(combination / "alignment.fasta")
    tree = esr.plain((combination / "tree.nwk").read_text())
    tip = "Candida.albicans.IME2"
    write_fasta(tmp_path / "full.fasta", alignment)
    (tmp_path / "full.nwk").write_text(tree)
    run(
        tmp_path,
        ["-s", "full.fasta", "-te", "full.nwk", "-blfix", "-m", "LG+G4", "--prefix", "fit"],
    )
    model = esr.fixed_model((tmp_path / "fit.iqtree").read_text())

    write_fasta(tmp_path / "held.fasta", esr.held_out_alignment(alignment, tip))
    (tmp_path / "held.nwk").write_text(esr.held_out_tree(tree, tip))
    run(
        tmp_path, ["-s", "held.fasta", "-te", "held.nwk", "-blfix", "-m", model, "--prefix", "held"]
    )
    pruned, _ = esr.prune(tree, tip)
    write_fasta(tmp_path / "pruned.fasta", {k: v for k, v in alignment.items() if k != tip})
    (tmp_path / "pruned.nwk").write_text(pruned)
    run(
        tmp_path,
        ["-s", "pruned.fasta", "-te", "pruned.nwk", "-blfix", "-m", model, "--prefix", "pruned"],
    )

    held = esr.log_likelihood((tmp_path / "held.iqtree").read_text())
    base = esr.log_likelihood((tmp_path / "pruned.iqtree").read_text())
    assert held - base == pytest.approx(2 * math.log(PI_LEU_LG), abs=2e-4)


def test_every_tip_and_level(rows: dict[tuple[str, str], dict[str, str]]) -> None:
    tips = list(read_fasta(REPLICATE / "input.fasta"))
    assert set(rows) == {(t, level) for t in tips for level in esr.LEVELS}
    for (tip, level), row in rows.items():
        assert int(row["sites"]) == 349, (tip, level)
        assert abs(float(row["accuracy"]) - float(row["mean_pp"])) < 0.05, (tip, level, row)
        e0 = rows[(tip, "e0")]
        assert abs(float(row["mean_pp"]) - float(e0["mean_pp"])) < 0.02, (tip, row)


def test_e0_matches_the_spike(rows: dict[tuple[str, str], dict[str, str]]) -> None:
    row = rows[("Candida.albicans.IME2", "e0")]
    assert (row["accuracy"], row["mean_pp"], row["log_score"]) == ("0.6218", "0.6444", "-1.2821")


def test_the_esr_job(
    combination: Path, tmp_path: Path, rows: dict[tuple[str, str], dict[str, str]]
) -> None:
    """The ESR job's workflow (main.nf --esr; phylobot.pipeline.esr) gives the engine's
    own rows, with every held-out tip its own task, and the real names put back."""
    from tests.golden.pipeline import ROOT, remove, run_pipeline

    out = ROOT / "results" / "golden" / "esr_job"
    remove(out)
    staged = out / "staged"
    staged.mkdir(parents=True)
    for name in ("alignment.fasta", "tree.nwk", "model.txt"):
        (staged / name).write_text((combination / name).read_text())
    tips = list(read_fasta(combination / "alignment.fasta"))
    settings = {"levels": list(esr.LEVELS), "tips": tips, "seed": 1}
    (staged / "settings.json").write_text(json.dumps(settings))
    run = run_pipeline(out, None, {"esr": True, "esr_dir": str(staged)})
    with (run / "esr" / "esr_tips.tsv").open() as f:
        got = {(r["tip"], r["level"]): r for r in csv.DictReader(f, delimiter="\t")}
    assert set(got) == set(rows)
    for key, row in rows.items():
        fields = ("pendant", "sites", "accuracy", "mean_pp", "log_score", "log_likelihood")
        assert {f: got[key][f] for f in fields} == {f: row[f] for f in fields}, key
    summary = json.loads((run / "esr" / "esr_summary.json").read_text())
    assert summary["levels"]["e0"]["sites"] == 349 * len(tips)
    remove(out)
