"""The pipeline's simulation steps (SIM_PREPARE and SIM_COLLECT in main.nf).

    phylobot-simulate prepare SPEC SEED WORKDIR [--tree TREE] [--root ROOT]
    phylobot-simulate collect SPEC WORKDIR OUTDIR --version VERSION_FILE --image IMAGE

`prepare` checks the spec, writes AliSim's inputs and plan.json into WORKDIR, plus one
`command<i>.args` file per AliSim run (NUL-separated arguments, run with `xargs -0`).
`collect` reads AliSim's outputs in WORKDIR into the truth bundle in OUTDIR. The tree
and root sequence come as files, never as paths inside the spec.
"""

import argparse
from pathlib import Path

from phylobot.simulate.engine import Plan, collect, prepare
from phylobot.simulate.families import FamilyError
from phylobot.simulate.spec import SpecError, load


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    steps = parser.add_subparsers(dest="step", required=True)
    p = steps.add_parser("prepare")
    p.add_argument("spec", type=Path)
    p.add_argument("seed", type=int)
    p.add_argument("work", type=Path)
    p.add_argument("--tree", type=Path)
    p.add_argument("--root", type=Path)
    c = steps.add_parser("collect")
    c.add_argument("spec", type=Path)
    c.add_argument("work", type=Path)
    c.add_argument("out", type=Path)
    c.add_argument("--version", type=Path, required=True, help="iqtree3 --version output")
    c.add_argument("--image", required=True)
    args = parser.parse_args()

    try:
        spec = load(args.spec)
    except SpecError as error:
        parser.exit(2, f"phylobot-simulate: {error}\n")
    if args.step == "prepare":
        if args.seed not in spec.seeds:
            parser.exit(2, f"phylobot-simulate: seed {args.seed} is not one of the spec's seeds\n")
        args.work.mkdir(parents=True, exist_ok=True)
        try:
            prepare(spec, args.seed, args.work, args.tree, args.root)
        except FamilyError as error:
            parser.exit(2, f"phylobot-simulate: {error}\n")
    else:
        version = args.version.read_text().splitlines()[0]
        collect(spec, Plan.read(args.work), args.work, args.out, version, args.image)


if __name__ == "__main__":
    main()
