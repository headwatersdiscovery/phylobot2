"""Running the real pipeline from golden tests, without filling `work/`.

Each run gets its own Nextflow work directory next to its results. It is deleted when
the run succeeds (its outputs are published to the run's results folder) and kept when
it fails, so the failing tasks' files and logs can be inspected. The shared `work/`
folder is only used by hand-run `make pipeline` (clean it with `make clean`).
"""

import json
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent.parent

#: Removes what pipeline containers wrote: on Linux they write as root, so the host user
#: cannot (Docker Desktop and OrbStack map ownership to the user, so on a Mac rmtree does).
CLEANER = "alpine:3.20.3"


def remove(path: Path) -> None:
    """Delete a folder a pipeline run wrote, root-owned files included."""
    try:
        shutil.rmtree(path)
    except FileNotFoundError:
        return
    except PermissionError:
        parent = path.resolve().parent
        command = ["docker", "run", "--rm", "-v", f"{parent}:/p", CLEANER,
                   "rm", "-rf", f"/p/{path.name}"]  # fmt: skip
        subprocess.run(command, check=True, capture_output=True)


#: Draws per ancestor for the PROPERTIES step in golden runs. The pipeline's own default
#: is 1000 (what users get); at that figure properties added 23% to the quick tier's
#: runtime, and no golden expectation reads a property file — the point of running the
#: step here is that it runs at all, on every scenario's real outputs. 50 draws keeps
#: that coverage for a fraction of the time. Scoring that needs real draws runs its own
#: pipeline, not this one.
GOLDEN_PROPERTY_DRAWS = 50


def run_pipeline(
    out: Path,
    fasta: Path | None,
    params: dict[str, Any] | None = None,
    make_args: Sequence[str] = (),
) -> Path:
    """Run `make pipeline` into out/run and return that directory. `make_args` go to
    make as they are (e.g. `-o tools-image`, or a RUNNER with environment for Nextflow)."""
    out.mkdir(parents=True, exist_ok=True)
    work = out / "work"
    command = ["make", "pipeline"]
    if fasta is not None:
        command.append(f"IN={fasta.resolve().relative_to(ROOT)}")
    command += [f"OUT={(out / 'run').relative_to(ROOT)}", f"WORK={work.relative_to(ROOT)}"]
    if params is not None:
        (out / "params.json").write_text(json.dumps(params))
        command.append(f"PARAMS={(out / 'params.json').relative_to(ROOT)}")
    command.append(f"ARGS=--property_draws {GOLDEN_PROPERTY_DRAWS}")
    command += list(make_args)
    subprocess.run(command, cwd=ROOT, check=True)  # raises (keeping work/) on failure
    remove(work)
    return out / "run"
