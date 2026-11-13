"""Packed property measurements (phylobot.property_store): what PROPERTIES writes and
the viewers read, one combination's draws in one Parquet file."""

import io
import subprocess
import sys
from pathlib import Path

from phylobot import property_store
from phylobot.property_store import Measured


def summary(node: str, version: str) -> dict[str, object]:
    return {"node": node, "gaps": version, "draws": 1000, "properties": {"charge": {"mean": 1.0}}}


def measured(nodes: int, draws: int = 1000) -> list[Measured]:
    return [
        Measured(version, f"Node{n}", {
            "charge": [float(n + d % 7) for d in range(draws)],
            "moment": [n / 100 + d / 1e6 for d in range(draws)],
        }, summary(f"Node{n}", version))
        for version in ("original", "fitch") for n in range(1, nodes + 1)
    ]  # fmt: skip


def test_one_node_round_trips(tmp_path: Path) -> None:
    property_store.write(tmp_path, measured(5, 50), ["charge", "moment"], 50)
    found = property_store.read_draws(tmp_path / "draws.parquet", "Node3", "fitch", ["charge"])
    assert found == {"charge": [float(3 + d % 7) for d in range(50)]}  # in draw order
    assert property_store.read_draws(tmp_path / "draws.parquet", "Node9", "fitch", ["charge"]) == {}
    text = (tmp_path / "summaries.json").read_text()
    assert property_store.summary(text, "Node3", "fitch") == summary("Node3", "fitch")
    assert property_store.summary(text, "Node3", "parsimony") is None
    assert property_store.measured_nodes(text, "original") == {f"Node{n}" for n in range(1, 6)}


def test_one_node_is_a_few_ranged_reads_of_a_large_file(tmp_path: Path) -> None:
    """The footer and the row group holding the node, not the file: the panel of one
    ancestor must not cost the whole combination's draws."""
    property_store.write(tmp_path, measured(120), ["charge", "moment"], 1000)
    data = (tmp_path / "draws.parquet").read_bytes()
    reads: list[int] = []

    def read_range(start: int, n: int) -> bytes:
        reads.append(n)
        return data[start : start + n]

    ranged = property_store.RangedFile(len(data), read_range)
    found = property_store.read_draws(
        io.BufferedReader(ranged, 64 * 1024), "Node77", "original", ["charge", "moment"]
    )
    assert len(found["charge"]) == len(found["moment"]) == 1000
    assert found["moment"][0] == 0.77
    assert sum(reads) < len(data) / 5, (sum(reads), len(data))


def test_a_property_with_fewer_draws_keeps_its_own(tmp_path: Path) -> None:
    rows = [Measured("original", "Node1", {"charge": [1.0, 2.0, 3.0], "moment": [0.5]}, {})]
    property_store.write(tmp_path, rows, ["charge", "moment"], 3)
    found = property_store.read_draws(
        tmp_path / "draws.parquet", "Node1", "original", ["charge", "moment", "absent"]
    )
    assert found == {"charge": [1.0, 2.0, 3.0], "moment": [0.5]}


def test_importing_it_does_not_import_pyarrow() -> None:
    code = "import sys, phylobot.property_store; print('pyarrow' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
