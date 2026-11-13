"""Species-tree reconciliation on the gene-family golden scenarios (reconciliation R5).

Slow and needs Docker, so excluded from `make test`; run with `make golden`.

Each replicate of tests/golden/reconcile/ runs the job the R5 validation report measured
(`reconciliation.run`: MAFFT, LG+G4, IQ-TREE's and GeneRax's trees, both reconciled
with the true species tree) and is scored against the true gene tree and its
duplication. The owner's thresholds apply to every replicate (expected_metrics.json):

- **noisy-duplication**, the pass/fail test of reconciliation: GeneRax recovers the
  gene tree and the duplication, and every chimeric ancestor of IQ-TREE's tree is
  flagged by the disagreement note.
- **hidden-paralogy**, a known limitation: GeneRax may follow the species tree into
  the wrong history, but every chimeric ancestor it gives must be flagged, and IQ-TREE's
  tree must have none.

The first seed of each scenario is the quick tier; the others run in `make golden-full`.
"""

import json
import tomllib
from pathlib import Path
from typing import Any

import pytest

from . import reconciliation

GOLDEN = Path(__file__).resolve().parent
ROOT = GOLDEN.parent.parent


def replicates() -> list[Any]:
    params = []
    for spec_path in sorted((GOLDEN / "reconcile").glob("*/spec.toml")):
        seeds = tomllib.loads(spec_path.read_text())["simulation"]["seeds"]
        for i, seed in enumerate(seeds):
            replicate = spec_path.parent / "replicates" / str(seed)
            if replicate.is_dir():
                marks = [pytest.mark.golden_full] if i > 0 else []
                name = str(replicate.relative_to(GOLDEN))
                params.append(pytest.param(replicate, marks=marks, id=name))
    return params


@pytest.mark.golden
@pytest.mark.parametrize("replicate", replicates())
def test_reconciliation(replicate: Path) -> None:
    scenario = replicate.parent.parent
    seed = int(replicate.name)
    out = ROOT / "results" / "golden" / replicate.relative_to(GOLDEN)
    newick = (replicate / "species_tree.nwk").read_text()
    genes = reconciliation.genes(replicate)
    run = reconciliation.run(out, replicate / "input.fasta", genes, newick, "protein", seed)
    scores = reconciliation.score(replicate, run)
    (out / "metrics.json").write_text(json.dumps(scores, indent=2) + "\n")
    print(f"\n{replicate.relative_to(GOLDEN)}\n{reconciliation.summary(scores)}")
    expected = json.loads((scenario / "expected_metrics.json").read_text())
    failures = reconciliation.check(scores, expected)
    assert not failures, "\n".join(failures)
