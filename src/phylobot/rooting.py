"""Rooting a tree on a user-chosen outgroup.

The root is placed at the midpoint of the branch above the outgroup's most recent
common ancestor (found with the tree held by an ingroup sequence, so that ancestor is
never the root itself). Internal node labels and branch lengths are kept; the new root
is unlabelled.

The outgroup does not have to form its own clade: users often pick a few
representatives of a broad outgroup. When it does not, the report lists the ingroup
sequences that fall inside the outgroup's clade, so the job page can warn about it.

    phylobot-root TREE OUTGROUP.txt ROOTED.nwk ROOTING.json
"""

import json
import sys
from pathlib import Path
from typing import Any

from phylobot.trees import load_tree


def root_on_outgroup(newick: str, outgroup: set[str]) -> tuple[str, dict[str, Any]]:
    """(rooted Newick, report). The report has the outgroup, whether it is monophyletic,
    and the ingroup sequences inside the outgroup's clade ("intruders")."""
    tree = load_tree(newick)
    leaves = {leaf.taxon.label: leaf for leaf in tree.leaf_node_iter()}
    unknown = outgroup - set(leaves)
    if unknown:
        raise ValueError(f"outgroup sequences not in the tree: {sorted(unknown)}")
    ingroup = sorted(set(leaves) - outgroup)
    if not outgroup or not ingroup:
        raise ValueError("need at least one outgroup and one ingroup sequence")

    # Hold the tree at an ingroup sequence's parent (an existing node, so no branch is
    # split) so that the outgroup's ancestor is not the root. A scattered outgroup can
    # span the held node's other subtrees; then try the next ingroup sequence.
    ancestor = None
    for name in ingroup:
        tree.reroot_at_node(leaves[name].parent_node, update_bipartitions=True)
        candidate = tree.mrca(taxon_labels=sorted(outgroup))
        if candidate is not tree.seed_node:
            ancestor = candidate
            break
    if ancestor is None:
        # The outgroup is spread across the whole tree: no branch separates it.
        tree.reroot_at_midpoint(update_bipartitions=True)
        tree.suppress_unifurcations()
        rooted_midpoint: str = tree.as_string(
            schema="newick", suppress_rooting=True, unquoted_underscores=True
        ).strip()
        return rooted_midpoint, {
            "outgroup": sorted(outgroup),
            "monophyletic": False,
            "intruders": ingroup,
            "method": "midpoint (the outgroup is spread across the tree)",
        }
    inside = {leaf.taxon.label for leaf in ancestor.leaf_iter()}
    intruders = sorted(inside - outgroup)
    length = ancestor.edge.length or 0.0
    tree.reroot_at_edge(
        ancestor.edge, length1=length / 2, length2=length / 2, update_bipartitions=True
    )
    tree.suppress_unifurcations()
    rooted: str = tree.as_string(
        schema="newick", suppress_rooting=True, unquoted_underscores=True
    ).strip()
    report = {
        "outgroup": sorted(outgroup),
        "monophyletic": not intruders,
        "intruders": intruders,
        "method": "outgroup",
    }
    return rooted, report


def main() -> None:
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    tree, outgroup_file, rooted_path, report_path = (Path(a) for a in sys.argv[1:5])
    outgroup = {line.strip() for line in outgroup_file.read_text().splitlines() if line.strip()}
    rooted, report = root_on_outgroup(tree.read_text(), outgroup)
    rooted_path.write_text(rooted + "\n")
    report_path.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
