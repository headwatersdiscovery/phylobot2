import re

import pytest

from phylobot.rooting import root_on_outgroup
from phylobot.trees import load_tree

TREE = (
    "(s1:1.17,(s2:1.07,(s5:0.48,(s6:0.47,s7:0.12)N4:0.22)N3:0.36)N2:0.12,"
    "(s3:0.39,s4:0.25)N5:0.27)N1;"
)


def total_length(newick: str) -> float:
    return round(sum(float(x) for x in re.findall(r":([0-9.]+)", newick)), 6)


def root_children(newick: str) -> list[set[str]]:
    tree = load_tree(newick)
    return [{leaf.taxon.label for leaf in c.leaf_iter()} for c in tree.seed_node.child_nodes()]


@pytest.mark.parametrize("outgroup", [{"s1"}, {"s3", "s4"}, {"s6"}, {"s5", "s6", "s7"}])
def test_monophyletic_outgroup_is_one_side_of_the_root(outgroup: set[str]) -> None:
    rooted, report = root_on_outgroup(TREE, outgroup)
    assert report == {"outgroup": sorted(outgroup), "monophyletic": True, "intruders": [],
                      "method": "outgroup"}  # fmt: skip
    assert outgroup in root_children(rooted)
    assert total_length(rooted) == total_length(TREE)  # the root halves one branch
    assert sorted(re.findall(r"\)(N\d)", rooted)) == ["N1", "N2", "N3", "N4", "N5"]


def test_scattered_outgroup_is_reported_not_rejected() -> None:
    rooted, report = root_on_outgroup(TREE, {"s1", "s3"})
    assert report["monophyletic"] is False and report["intruders"] == ["s4"]
    assert {"s1", "s3", "s4"} in root_children(rooted)


def test_outgroup_spanning_the_held_node_tries_other_ingroup_leaves() -> None:
    _, report = root_on_outgroup(TREE, {"s2", "s3"})
    assert report["method"] == "outgroup" and report["intruders"] == ["s1", "s4"]


def test_bad_outgroups() -> None:
    with pytest.raises(ValueError, match="not in the tree"):
        root_on_outgroup(TREE, {"nope"})
    with pytest.raises(ValueError, match="at least one"):
        root_on_outgroup(TREE, set())
