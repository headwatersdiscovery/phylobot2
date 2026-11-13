"""Does this machine reproduce the golden simulations, and does the CPU matter?

Regenerates golden replicates with the pinned IQ-TREE image natively and, given a static
qemu, as older or newer x86 CPUs (qemu's floating point is bit-exact, x87 included), and
says which match the committed files. The committed replicates were generated on native
x86-64; some do not regenerate byte for byte under emulation (Rosetta on Apple Silicon).

    uv run python -m tests.golden.platform_check lysin/codon-easy/301 [...]
    uv run python -m tests.golden.platform_check --qemu DIR --cpu Nehalem lysin/codon-easy/301

DIR holds qemu-x86_64-static, fetched without installing anything on the host:

    docker run --rm -v DIR:/out python:3.14.7-slim sh -c 'apt-get update -qq &&
      apt-get install -y -qq qemu-user-static && cp /usr/bin/qemu-x86_64-static /out/'
"""

import argparse
import filecmp
import hashlib
import subprocess
import tempfile
from collections.abc import Callable
from pathlib import Path

from phylobot.simulate import engine

from .generate import GOLDEN, generate


def emulated(qemu: Path, cpu: str) -> Callable[[Path, str, list[str]], None]:
    """engine.run_docker, with IQ-TREE run under qemu as `cpu`."""

    def run_docker(work: Path, image: str, command: list[str]) -> None:
        docker = ["docker", "run", "--rm", "--platform", "linux/amd64", "-v", f"{qemu}:/q"]
        docker += ["-v", f"{work}:/w", "-w", "/w", image]
        binary = ["/q/qemu-x86_64-static", "-cpu", cpu, "/usr/local/bin/iqtree3"]
        subprocess.run(docker + binary + command[1:], check=True, capture_output=True, text=True)

    return run_docker


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("replicates", nargs="+", help="scenario/seed, e.g. lysin/codon-easy/301")
    parser.add_argument("--qemu", type=Path, help="directory holding qemu-x86_64-static")
    parser.add_argument("--cpu", default="Nehalem", help="qemu CPU model (with --qemu)")
    args = parser.parse_args()
    label = f"qemu {args.cpu}" if args.qemu else "native"
    if args.qemu:
        engine.run_docker = emulated(args.qemu.resolve(), args.cpu)  # type: ignore[assignment]  # replaced on purpose
    for replicate in args.replicates:
        name, seed = replicate.rsplit("/", 1)
        committed = GOLDEN / name / "replicates" / seed
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            generate(GOLDEN / name, int(seed), out)
            differ = [
                p.name
                for p in sorted(committed.iterdir())
                if p.name != "provenance.json"  # it holds the others' checksums
                and not filecmp.cmp(out / p.name, p, shallow=False)
            ]
            digest = hashlib.sha256((out / "true_alignment.fasta").read_bytes()).hexdigest()
        verdict = "matches" if not differ else "differs: " + ", ".join(differ)
        print(f"{label:16s} {replicate:36s} alignment {digest[:12]}  {verdict}")


if __name__ == "__main__":
    main()
