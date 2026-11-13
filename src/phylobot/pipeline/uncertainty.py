"""Where a finished job's ancestral uncertainty comes from.

    phylobot-uncertainty PUBLISHED SETTINGS OUTDIR

PUBLISHED is the parent job's staged results (summary.json and its combination folders,
as the web tier copied them), SETTINGS the run's settings JSON
(`uncertainty_run.Settings`). OUTDIR gets uncertainty/: summary.json, nodes.tsv, and
per ancestor sites/, drivers/ and pairs/ (numbered files, index.tsv naming the nodes).

Every path read is built from a combination key the staged summary.json names, so each
key is checked against the pipeline's own pattern before it is used.
"""

import json
import re
import sys
from pathlib import Path

from phylobot import uncertainty_run

KEY = re.compile(r"^[A-Za-z0-9_.+-]+/[A-Za-z0-9_.+-]+/[A-Za-z0-9_.+-]+$")
NAME = re.compile(r"^[A-Za-z0-9_.+-]+(/[A-Za-z0-9_.+-]+)*$")


def reader(root: Path) -> uncertainty_run.Read:
    def read(relative: str) -> str | None:
        if not NAME.match(relative) or ".." in relative.split("/"):
            raise ValueError(f"refusing to read {relative!r}")
        path = root / relative
        return path.read_text() if path.is_file() else None

    return read


def check_keys(summary_text: str | None) -> None:
    summary = json.loads(summary_text or "{}")
    for combo in summary.get("combinations", []):
        key = combo.get("key", "")
        parts = [combo.get("aligner"), combo.get("tree"), combo.get("model")]
        if not KEY.match(key) or key != "/".join(str(p) for p in parts):
            raise ValueError(f"unexpected combination key {key!r}")


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    published, settings_path, outdir = (Path(a) for a in args)
    read = reader(published)
    try:
        check_keys(read("summary.json"))
        settings = uncertainty_run.Settings.parse(json.loads(settings_path.read_text()))
        files = uncertainty_run.run(read, settings)
    except (uncertainty_run.SettingsError, ValueError) as error:
        print(f"phylobot-uncertainty: {error}", file=sys.stderr)
        return 1
    out = outdir / uncertainty_run.OUT
    for relative, text in files.items():
        path = out / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
