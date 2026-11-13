from dataclasses import dataclass as _dataclass

from phylobot import treeview

# One near-zero branch (d) and one very long one (e); the rest about 0.1.
TREE = "(out:0.1,((a:0.1,b:0.12)Node3:0.1,(c:0.1,d:0.0001)Node4:0.08)Node2:0.1,e:5.0)Node1;"


def test_adjusted_mode_clamps_extreme_branches() -> None:
    drawing = treeview.layout(TREE, "adjusted", ancestors={"Node1", "Node2", "Node3", "Node4"})
    by_child = {b.length: b for b in drawing.branches}
    assert by_child[5.0].clamped == "long" and by_child[0.0001].clamped == "short"
    assert all(b.clamped == "" for length, b in by_child.items() if length not in (5.0, 0.0001))
    # the long branch no longer dominates: e ends near the others
    leaves = {leaf.name: leaf.x for leaf in drawing.leaves}
    assert leaves["e"] < 2 * leaves["a"]
    true = {leaf.name: leaf.x for leaf in treeview.layout(TREE, "true").leaves}
    assert true["e"] > 5 * true["a"]
    # the short branch is visible: d sits clearly right of its parent
    node4 = next(n for n in drawing.nodes if n.name == "Node4")
    assert leaves["d"] - node4.x > 10
    assert "1 short branch lengthened; 1 long branch cut (marked //)" in drawing.note
    assert drawing.scale_bar is None


def test_cladogram_aligns_extant_sequences() -> None:
    drawing = treeview.layout(TREE, "cladogram")
    assert len({leaf.x for leaf in drawing.leaves}) == 1
    assert drawing.scale_bar is None and drawing.note == ""


def test_true_lengths_have_a_scale_bar() -> None:
    drawing = treeview.layout(TREE, "true")
    assert drawing.scale_bar is not None
    x1, x2, _, label = drawing.scale_bar
    assert x2 > x1 and label == "1 substitutions per site"
    assert treeview.nice_length(0.37) == 0.2 and treeview.nice_length(7) == 5
    assert treeview.nice_length(0.1) == 0.1


def test_nodes_leaves_and_links() -> None:
    drawing = treeview.layout(
        TREE, "nonsense", ancestors={"Node2", "Node4"}, support={"Node2": "97", "Node4": "-"},
        outgroup={"out"},
    )  # fmt: skip
    assert drawing.mode == "adjusted"  # unknown modes fall back
    nodes = {n.name: n for n in drawing.nodes}
    assert nodes["Node2"].link and not nodes["Node3"].link
    assert nodes["Node2"].support == "97" and nodes["Node4"].support == ""
    assert nodes["Node2"].leaves == 4 and nodes["Node1"].leaves == 6
    assert [leaf.name for leaf in drawing.leaves if leaf.outgroup] == ["out"]
    # rows evenly spaced; each parent midway between its first and last child
    ys = [leaf.y for leaf in drawing.leaves]
    assert len({round(b - a, 3) for a, b in zip(ys, ys[1:], strict=False)}) == 1
    assert drawing.height == treeview.MARGIN["top"] + 6 * treeview.ROW + treeview.MARGIN["bottom"]
    assert drawing.stem_from is None  # labelled top node: no added root
    rooted = treeview.layout("((a:1,b:1)Node2:0.5,c:1);")
    assert rooted.stem_from is not None and rooted.nodes[0].name == ""


def test_room_for_a_labelled_root() -> None:
    plain = treeview.layout(TREE, "true")
    supported = treeview.layout(TREE, "true", support={"Node1": "100"})
    root = next(n for n in plain.nodes if n.name == "Node1")
    root_with_support = next(n for n in supported.nodes if n.name == "Node1")
    # the label ("Node1", or "Node1 100") ends at x - 6: it must start inside the image
    assert root.x - 6 - len("Node1") * treeview.NODE_CHAR_WIDTH >= treeview.MARGIN["left"] - 1
    width = len("Node1 100") * treeview.NODE_CHAR_WIDTH
    assert root_with_support.x - 6 - width >= treeview.MARGIN["left"] - 1


def test_branch_display_by_support_and_dnds() -> None:
    from phylobot.paths import BranchRate

    supports = {"Node2": {"sh_alrt": 85.0, "ufboot": 97.0}, "Node3": {"sh_alrt": 50.0},
                "Node4": {"ufboot": 99.0}}  # fmt: skip
    rates = {
        frozenset(("Node1", "Node2")): BranchRate(0.1, 0.01, 0.1, 3, 10),
        frozenset(("Node2", "Node3")): BranchRate(999, 0.02, 0.0, 6, 0),  # undefined
        frozenset(("Node4", "c")): BranchRate(10, 0.05, 0.005, 15, 1.5),
    }
    plain = treeview.layout(TREE, supports=supports, rates=rates)
    assert all(b.plain for b in plain.branches) and not plain.legends
    labels = {n.name: n.support for n in plain.nodes}
    assert labels["Node2"] == "85/97" and labels["Node3"] == "50/-" and labels["Node4"] == "-/99"
    assert plain.height == TALL  # no legends, no extra height

    display = treeview.Display(colour="dnds", width="ufboot", labels=("ufboot",), stars=True,
                               dnds_labels=True)  # fmt: skip
    drawing = treeview.layout(TREE, supports=supports, rates=rates, display=display)
    branch = {b.tooltip.split("\n")[0]: b for b in drawing.branches}
    to2, to3, to_c = branch["Node1 → Node2"], branch["Node2 → Node3"], branch["Node4 → c"]
    assert to2.stroke == treeview.dnds_colour(0.1) == "#4575b4"  # blue: dN/dS 0.1
    assert to_c.stroke == "#d73027" and to_c.label == "10.00"  # red: 10
    assert to3.stroke == treeview.MISSING and to3.dashed and to3.label == "–"  # undefined
    assert to2.stroke_width == treeview.support_width(97) and to3.stroke_width == 1.6
    assert "UFBoot 97%" in to2.tooltip and "dN/dS 0.1;" in to2.tooltip
    assert "undefined" in to3.tooltip
    labels = {n.name: n.support for n in drawing.nodes}
    assert labels["Node2"] == "***" and labels["Node4"] == "****" and labels["Node3"] == ""
    assert treeview.stars("sh_alrt", 60) == "·" and treeview.stars("sh_alrt", 80) == "**"
    titles = [g.title for g in drawing.legends]
    assert titles == ["Branch colour: dN/dS (codeml free ratios)",
                      "Branch thickness: UFBoot support (%)"]  # fmt: skip
    # two legend rows, above the tree: the tree is drawn `top` lower
    assert drawing.top == treeview.LEGEND_PAD + 2 * treeview.LEGEND_ROW
    assert drawing.height == TALL + drawing.top
    assert all(g.y < drawing.top for g in drawing.legends)
    # The ancestor viewer's palettes: support as a probability (warm weak, cool strong);
    # dN/dS cool below 1, the pale middle at 1, warm above.
    assert treeview.support_colour(40) == "#fc8d59" and treeview.support_colour(100) == "#4575b4"
    assert treeview.dnds_colour(1.0) == "#fee090"
    assert treeview.dnds_colour(0.1, "viridis") == "#440154"
    assert treeview.support_colour(80, "nonsense") == treeview.support_colour(80)
    greys = treeview.layout(TREE, supports=supports, rates=rates,
                            display=treeview.Display(colour="ufboot", palette="greys"))  # fmt: skip
    assert next(b for b in greys.branches if "UFBoot 97%" in b.tooltip).stroke == "#4b4b4b"
    # a traced path: plain branches turn orange, styled ones keep their colour
    traced = treeview.layout(TREE, path=["Node1", "Node2"], supports=supports, rates=rates,
                             display=display)  # fmt: skip
    on = next(b for b in traced.branches if b.on_path)
    assert on.drawn_stroke == on.stroke != "#e08214" and on.halo_width > on.stroke_width
    on_plain = next(b for b in treeview.layout(TREE, path=["Node1", "Node2"]).branches if b.on_path)
    assert (on_plain.drawn_stroke, on_plain.drawn_width) == ("#e08214", 4.0)


TALL = treeview.MARGIN["top"] + 6 * treeview.ROW + treeview.MARGIN["bottom"]


def test_reconciliation_events_and_losses() -> None:
    newick = "((a:1,b:1)N2:1,(c:1,d:1)N3:1,e:1)N1;"
    split = frozenset({"c", "d"})  # canonical: the side without "a"
    drawing = treeview.layout(
        newick,
        events={"N2": treeview.Event("duplication", "Saccharomycetaceae")},
        losses={split: ["Candida albicans", "Yarrowia lipolytica"]},
    )
    marked = {n.name: n.event for n in drawing.nodes}
    assert marked["N2"] == treeview.Event("duplication", "Saccharomycetaceae")
    assert marked["N3"] is None
    tips = {b.tooltip.splitlines()[0]: b.tooltip for b in drawing.branches}
    assert "2 copies lost here (in Candida albicans; Yarrowia lipolytica)" in tips["N1 → N3"]
    assert drawing.events_y is not None
    assert treeview.layout(newick).events_y is None


# Folding a clade into a triangle (collapse_clade): the member tips hide, the subtree
# becomes one short band, and the fold follows whatever maps to it.
@_dataclass
class _Style:
    label: str = ""
    fill: str = "#cfd8e6"
    border_colour: str = "#2457a6"
    border_width: float = 1.6
    alike: bool = True


def test_collapsing_a_clade_folds_it_into_a_triangle() -> None:
    full = treeview.layout(TREE, "true", ancestors={"Node1", "Node2", "Node3", "Node4"})
    folded = treeview.layout(
        TREE, "true", ancestors={"Node1", "Node2", "Node3", "Node4"},
        collapsed={"Node3": _Style(label="AB clade")},
    )  # fmt: skip
    names = {leaf.name for leaf in folded.leaves}
    assert "a" not in names and "b" not in names  # the clade's tips are hidden
    assert {"out", "c", "d", "e"} <= names
    assert [c.name for c in folded.collapsed] == ["Node3"]
    clade = folded.collapsed[0]
    assert clade.leaves == 2 and clade.label == "AB clade" and clade.alike
    assert clade.base_x > clade.apex_x  # apex at the node, base at the deepest tip
    assert clade.top < clade.y < clade.bottom
    # Node3 is drawn as a triangle, not a dot or a connector.
    assert all(n.name != "Node3" for n in folded.nodes)
    # Its incoming branch still reaches the apex; no branch leads to a hidden tip.
    assert any(b.x2 == clade.apex_x and b.y == clade.y for b in folded.branches)
    assert not any(b.tooltip.startswith("Node3 → a") for b in folded.branches)
    # Folding saves vertical space (two tips become one short band).
    assert folded.height < full.height


def test_nested_folds_keep_only_the_outer_triangle() -> None:
    folded = treeview.layout(TREE, "true", collapsed={"Node2": _Style(), "Node3": _Style()})
    assert [c.name for c in folded.collapsed] == ["Node2"]  # the inner fold is hidden inside it
    assert {leaf.name for leaf in folded.leaves} == {"out", "e"}  # a, b, c, d all hidden
    assert folded.collapsed[0].leaves == 4


def test_an_unknown_or_absent_fold_is_ignored() -> None:
    assert treeview.layout(TREE).collapsed == []
    assert treeview.layout(TREE, collapsed={"Nonsense": _Style()}).collapsed == []
    # Unfolded, the drawing is unchanged (same height as before).
    assert treeview.layout(TREE).height == TALL


def test_a_folded_triangle_carries_its_style() -> None:
    folded = treeview.layout(
        TREE, collapsed={"Node3": _Style("AB", "#dfe7d6", "#3c7a3c", 2.5, alike=False)}
    )
    clade = folded.collapsed[0]
    assert clade.fill == "#dfe7d6" and clade.border_colour == "#3c7a3c"
    assert clade.border_width == 2.5 and clade.alike is False
    # points: apex then the two base corners (one polygon string)
    assert clade.points.startswith(f"{clade.apex_x},{clade.y} {clade.base_x},")
