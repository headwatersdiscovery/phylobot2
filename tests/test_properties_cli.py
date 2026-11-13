"""`phylobot-properties`, the PROPERTIES process's command.

Builds a published combination on disk as EXPORT writes one, runs the command over it,
and reads back what the viewers read. The run on real data is the golden tier's; what is
checked here is which ancestors a run measures, which is what B5c's scopes turn into.
"""

from __future__ import annotations

import json
from pathlib import Path

from phylobot import property_store
from phylobot.pipeline import properties as cli

POSTERIORS = "node\tsite\tM\tK\tR\n" + "".join(
    f"{node}\t{site}\t0.7\t0.2\t0.1\n" for node in ("Node1", "Node2", "Node3") for site in (1, 2, 3)
)


def combination(tmp_path: Path) -> Path:
    """One combination's published files: three ancestors, posteriors, a masked copy."""
    folder = tmp_path / "mafft" / "iqtree" / "LG+G4"
    folder.mkdir(parents=True)
    (folder / "ancestors.fasta").write_text(">Node1\nMKR\n>Node2\nMKK\n>Node3\nMRK\n")
    (folder / "ancestors.gaps-fitch.fasta").write_text(">Node1\nMK-\n>Node2\nMKK\n>Node3\nMRK\n")
    (folder / "ancestors_posteriors.tsv").write_text(POSTERIORS)
    return folder


def run(folder: Path, out: Path, *extra: str) -> int:
    return cli.main([
        str(folder), str(out), "--properties", "charge", "--draws", "20", *extra
    ])  # fmt: skip


def measured(out: Path) -> set[tuple[str, str]]:
    """(node, gap version) pairs that were written."""
    if not (out / property_store.SUMMARIES).exists():
        return set()
    nodes = json.loads((out / property_store.SUMMARIES).read_text())["nodes"]
    return {(node, version) for version, found in nodes.items() for node in found}


class TestWhichAncestorsAreMeasured:
    def test_every_node_by_default(self, tmp_path: Path) -> None:
        out = tmp_path / "out"
        assert run(combination(tmp_path), out) == 0
        assert measured(out) == {(f"Node{n}", "original") for n in (1, 2, 3)}

    def test_only_the_named_nodes(self, tmp_path: Path) -> None:
        """B5c's scopes are node lists: measuring one ancestor must cost one ancestor,
        or "this ancestor" would quietly cost the same as "every ancestor"."""
        out = tmp_path / "out"
        assert run(combination(tmp_path), out, "--nodes", "Node2") == 0
        assert measured(out) == {("Node2", "original")}
        assert json.loads((out / "index.json").read_text())["asked"] == ["Node2"]

    def test_each_gap_version_asked_for(self, tmp_path: Path) -> None:
        out = tmp_path / "out"
        assert run(combination(tmp_path), out, "--nodes", "Node1", "--gap-versions", "fitch") == 0
        assert measured(out) == {("Node1", "original"), ("Node1", "fitch")}

    def test_a_node_that_is_not_here_is_reported_not_fatal(
        self, tmp_path: Path, capsys: object
    ) -> None:
        """The plan comes from the web tier's trees; if it names an ancestor these files
        do not hold, the rest are still measured and the mismatch is said out loud."""
        out = tmp_path / "out"
        assert run(combination(tmp_path), out, "--nodes", "Node2,Node9") == 0
        assert measured(out) == {("Node2", "original")}
        assert "Node9" in capsys.readouterr().err  # type: ignore[attr-defined]

    def test_the_summary_is_what_the_viewer_reads(self, tmp_path: Path) -> None:
        out = tmp_path / "out"
        run(combination(tmp_path), out, "--nodes", "Node1")
        body = property_store.summary(
            (out / property_store.SUMMARIES).read_text(), "Node1", "original"
        )
        assert body is not None and body["node"] == "Node1" and body["draws"] == 20
        assert "charge" in body["properties"] and "mean" in body["properties"]["charge"]
        draws = property_store.read_draws(
            out / property_store.DRAWS, "Node1", "original", ["charge"]
        )
        assert len(draws["charge"]) == 20  # one value per draw
        # Three files for the combination, whatever the number of nodes and versions.
        names = sorted(p.name for p in out.iterdir())
        assert names == ["draws.parquet", "index.json", "summaries.json"]
