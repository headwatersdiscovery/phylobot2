"""Write ancestor_map.tsv: every combination's ancestors by content-derived IDs.

    phylobot-ancestor-map TREES.tsv OUT.tsv

TREES.tsv has one line per combination: its key (aligner/tree/model), "rooted" or
"unrooted", and its published tree's Newick (tree.rooted.nwk with an outgroup, else
tree.nwk; sequence names restored). OUT.tsv has one row per labelled internal node:
combination, node, clade (the sequences below it, hashed; blank on an unrooted tree),
parts (the subtrees that meet at it, hashed) and leaves (how many sequences are below
it). Two combinations' nodes are the same ancestor when their clades match (with an
outgroup) or their parts do (without one): phylobot.analogs.
"""

import sys
from pathlib import Path

from phylobot.analogs import MapRow, map_rows, write_map


def build(trees: str) -> str:
    rows: list[MapRow] = []
    for line in trees.splitlines():
        if not line.strip():
            continue
        combination, kind, newick = line.split("\t", 2)
        if kind not in ("rooted", "unrooted"):
            raise ValueError(f"{combination}: unknown tree kind {kind!r}")
        rows += map_rows(combination, newick.strip(), kind == "rooted")
    return write_map(rows)


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    Path(sys.argv[2]).write_text(build(Path(sys.argv[1]).read_text()))


if __name__ == "__main__":
    main()
