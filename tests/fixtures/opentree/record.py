"""Record the Open Tree of Life's answers that tests/test_species.py replays.

    uv run python tests/fixtures/opentree/record.py

Sends only public species names and taxids (the owner agreed to these, 2026-10-16):
the Ime2 and lysin examples' species and a few deliberately awkward names. Each
scenario is one JSON file of exchanges (endpoint, payload, status, body). Rerecording
after a new synthesis release changes the expected trees, so check the test diffs.
"""

import json
import sys
from importlib.resources import files
from pathlib import Path
from typing import Any

from phylobot.pipeline.validate import parse_fasta
from phylobot.species import resolve
from phylobot.species.opentree import Answer, Client, Http
from phylobot.species.table import Row

HERE = Path(__file__).parent

AWKWARD = [
    Row("misspelt", query="Saccharomyces cerevisae"),
    Row("hybrid", query="Saccharomyces pastorianus"),
    Row("genus", query="Candida"),
    Row("inside_genus", query="Candida albicans"),
    Row("nonsense", query="Notarealname foo"),
    Row("by_taxid", query="Saccharomyces cerevisiae S288C", taxid=559292),
    Row("strain_species", query="Saccharomyces cerevisiae"),
    Row("pinned", query="ott213939"),
    Row("pinned_pruned", query="ott6376729"),
    Row("gene_word", query="Candida parapsilosis ime2"),
    Row("broken", query="Crustacea"),
    Row("bird", query="Gallus gallus"),
    Row("taxid_only", taxid=4952),
]


def example_records(name: str) -> list[tuple[str, str]]:
    text = files("phylobot").joinpath("examples", f"{name}.fasta").read_text()
    return [(n, d) for n, d, _ in parse_fasta(text, name)]


class Recorder:
    def __init__(self) -> None:
        self.http = Http()
        self.exchanges: list[dict[str, Any]] = []

    def __call__(self, endpoint: str, payload: dict[str, Any]) -> Answer:
        answer = self.http(endpoint, payload)
        self.exchanges.append(
            {"endpoint": endpoint, "payload": payload, "status": answer.status, "body": answer.body}
        )
        return answer


def record(name: str, rows: list[Row]) -> None:
    recorder = Recorder()
    resolve.with_open_tree(rows, Client(recorder))
    path = HERE / f"{name}.json"
    path.write_text(json.dumps(recorder.exchanges, indent=1, sort_keys=True) + "\n")
    print(f"{path}: {len(recorder.exchanges)} exchanges", file=sys.stderr)


def main() -> None:
    record("ime2", resolve.detect(example_records("ime2")))
    record("lysin", resolve.detect(example_records("lysin")))
    record("awkward", AWKWARD)


if __name__ == "__main__":
    main()
