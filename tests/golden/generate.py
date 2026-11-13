"""Generate golden datasets by simulating sequence evolution with IQ-TREE's AliSim.

Each scenario is a directory holding a spec.toml (`phylobot.simulate.spec`). For every
seed in the spec, this writes replicates/<seed>/ with the truth bundle described in
`phylobot.simulate`: input.fasta (the only thing the pipeline sees), the true
alignment, ancestors and tree, true site classes for codeml M3 scenarios, and
provenance.json. The simulation itself is `phylobot.simulate.engine`, shared with the
pipeline's --simulate and with users' simulations; test_simulate.py checks it still
reproduces every committed replicate byte for byte.

Committed outputs are frozen truth; regenerating them needs the project owner's
sign-off, so existing replicates are never overwritten without --force.

    uv run python tests/golden/generate.py                 # all scenarios
    uv run python tests/golden/generate.py tests/golden/ime2/easy
"""

import argparse
from pathlib import Path

from phylobot.simulate import load
from phylobot.simulate.engine import simulate_locally
from phylobot.simulate.spec import input_files

GOLDEN = Path(__file__).resolve().parent


def generate(scenario: Path, seed: int, out: Path) -> None:
    """One replicate of a scenario into `out`."""
    spec = load(scenario / "spec.toml")
    tree, root = input_files(spec, scenario)
    simulate_locally(spec, seed, out, tree, root)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("scenarios", nargs="*", type=Path, help="scenario dirs (default: all)")
    parser.add_argument("--force", action="store_true", help="overwrite committed replicates")
    args = parser.parse_args()

    scenarios = args.scenarios or sorted(p.parent for p in GOLDEN.rglob("spec.toml"))
    for scenario in scenarios:
        scenario = scenario.resolve()
        for seed in load(scenario / "spec.toml").seeds:
            out = scenario / "replicates" / str(seed)
            label = f"{scenario.relative_to(GOLDEN)} seed {seed}"
            if out.exists() and not args.force:
                print(f"skip {label}: already generated (use --force to overwrite)")
                continue
            generate(scenario, seed, out)
            print(f"wrote {label}")


if __name__ == "__main__":
    main()
