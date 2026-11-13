"""A mutational path followed through other combinations (phylobot.path_compare), the
tree walks it needs (paths.between, paths.branch_lengths) and its line (path_line)."""

from pathlib import Path

from phylobot import analogs, path_compare, path_line, paths
from phylobot.path_compare import ABSENT, DIFFERENT, SAME, Side

TREE = "(out:1,((a:0.1,b:0.2)Node3:0.3,(c:0.1,d:0.1)Node4:0.4)Node2:0.5)Node1;"


def test_between_walks_up_and_down_whatever_the_root() -> None:
    assert paths.between(TREE, "Node1", "c") == paths.path(TREE, "Node1", "c")
    assert paths.between(TREE, "Node3", "c") == ["Node3", "Node2", "Node4", "c"]  # over the top
    assert paths.between(TREE, "c", "Node2") == ["c", "Node4", "Node2"]  # upwards
    assert paths.between(TREE, "a", "a") == ["a"]
    assert paths.between(TREE, "Node1", "nope") == []
    # an unlabelled node on the way (the root outgroup rooting adds) is skipped
    rooted = "((a:1,b:1)Node2:1,(c:1,d:1)Node3:1);"
    assert paths.between(rooted, "Node2", "d") == ["Node2", "Node3", "d"]


def test_branch_lengths_follow_the_path() -> None:
    names = paths.path(TREE, "Node1", "a")
    assert names == ["Node1", "Node2", "Node3", "a"]
    assert paths.branch_lengths(TREE, names) == [0.5, 0.3, 0.1]
    rooted = "((a:1,b:1)Node2:1,(c:1,d:1)Node3:2);"
    assert paths.branch_lengths(rooted, ["Node2", "Node3"]) == [3.0]  # through the root


def side(key: str, newick: str, alignment: dict[str, str], ancestors: dict[str, str]) -> Side:
    return Side(key, key, alignment, ancestors, newick, analogs.places(newick))


def main_side() -> Side:
    #                 123456
    alignment = {"a": "MKL-AQ", "b": "MRL-AQ", "c": "MKLWAE", "d": "MKLWAE"}
    return side(
        "main", "((a,b)Node2,(c,d)Node3)Node1;", alignment,
        {"Node1": "MKLWAE", "Node2": "MKL-AQ", "Node3": "MKLWAE"},
    )  # fmt: skip


def test_same_alignment_another_reconstruction() -> None:
    main = main_side()
    # the same columns and tree, nodes numbered otherwise; it keeps W at Node2's twin
    # and reads Q where the traced path reads E at the top
    other = side(
        "other", "((a,b)NodeB,(c,d)NodeC)NodeA;", main.alignment,
        {"NodeA": "MKLWAQ", "NodeB": "MKLWAQ", "NodeC": "MKLWAE"},
    )  # fmt: skip
    route = ["Node1", "Node2", "a"]
    found = path_compare.compare(main, route, [other], by_clade=True)
    assert found.changes == [{3: ("W", "-"), 5: ("E", "Q")}, {}]
    track = found.followed[0]
    assert track.twins == {"Node1": "NodeA", "Node2": "NodeB", "a": "a"}
    assert track.route == ["NodeA", "NodeB", "a"] and track.missing == []
    # NodeA -> NodeB changes nothing: both of the traced branch's changes are absent
    assert {c: s.verdicts[0].state for c, s in found.support[0].items()} == {3: ABSENT, 5: ABSENT}
    first, second = found.branches
    assert (first[0].comparable, first[0].direct, first[0].absent, first[0].extra) == (
        True, True, 2, 0,
    )  # fmt: skip
    assert second[0].extra == 1  # NodeB -> a loses the W there instead
    assert found.levels(0) == {"all": 0, "most": 0, "few": 0, "none": 2, "na": 0}


def test_another_alignment_is_laid_on_the_traced_columns() -> None:
    main = main_side()
    # the gap one column earlier: residues land in other columns, the same changes
    alignment = {"a": "MK-LAQ", "b": "MR-LAQ", "c": "MKWLAE", "d": "MKWLAE"}
    other = side(
        "shifted", "((a,b)X,(c,d)Y)Z;", alignment,
        {"Z": "MKWLAE", "X": "MK-LAE", "Y": "MKWLAE"},
    )  # fmt: skip
    route = ["Node1", "Node2", "a"]
    found = path_compare.compare(main, route, [other], by_clade=True)
    columns = found.followed[0].columns
    assert columns is not None and columns[5] == 5 and columns[4] == 4
    # W lost on the first branch in both, but at another residue's column there: the
    # shifted alignment loses its column 2, which holds the traced column 3's residues' L...
    states = {c: s.verdicts[0].state for c, s in found.support[0].items()}
    assert states[5] == ABSENT  # E -> Q happens a branch later there
    assert found.support[1] == {}  # the traced path changes nothing into a
    assert found.branches[1][0].extra == 1  # ...and there it is: E -> Q into a
    # scrolling together: its columns against the traced ones, never backwards
    sync = path_compare.sync_columns(columns, 6, 6)
    assert sync == sorted(sync) and len(sync) == 6 and sync[-1] == 5
    assert path_compare.sync_on(alignment["a"], main.alignment["a"]) == [0, 1, 2, 2, 4, 5]
    assert path_compare.sync_on("MK", "MR") is None  # not the same sequence


def test_support_levels_and_verdicts() -> None:
    main = main_side()

    def twin(key: str, top: str, below: str) -> Side:
        return side(
            key, "((a,b)Node2,(c,d)Node3)Node1;", main.alignment,
            {"Node1": top, "Node2": below, "Node3": top},
        )  # fmt: skip

    others = [
        twin("same", "MKLWAE", "MKL-AQ"),  # both changes
        twin("half", "MKLWAE", "MKL-AD"),  # the loss, and E -> D instead of E -> Q
        twin("none", "MKLWAE", "MKLWAE"),  # neither
    ]
    found = path_compare.compare(main, ["Node1", "Node2"], others, by_clade=True)
    loss, swap = found.support[0][3], found.support[0][5]
    assert (loss.agreeing, loss.compared, loss.level) == (2, 3, "most")
    assert (swap.agreeing, swap.compared, swap.level) == (1, 3, "few")
    assert [v.state for v in swap.verdicts] == [SAME, DIFFERENT, ABSENT]
    assert swap.verdicts[1].call == ("E", "D") and swap.sides(ABSENT) == ["none"]
    assert found.followed[0].shared == {"Node2": {3: (0, 3), 5: (0, 5)}}
    assert found.followed[1].shared == {"Node2": {3: (0, 3)}}
    assert path_compare.Support().level == "na"
    assert path_compare.Support((path_compare.Verdict("x", SAME),)).level == "all"


def test_a_tree_without_the_branch_gives_no_verdict() -> None:
    main = main_side()
    # another topology: no (a,b) clade, so Node2 has no twin
    other = side(
        "other", "((a,c)P,(b,d)Q)R;", main.alignment,
        {"R": "MKLWAE", "P": "MKLWAE", "Q": "MKLWAE"},
    )  # fmt: skip
    found = path_compare.compare(main, ["Node1", "Node2", "a"], [other], by_clade=True)
    track = found.followed[0]
    assert track.missing == ["Node2"] and track.twins == {"Node1": "R", "a": "a"}
    assert track.route == ["R", "P", "a"]  # its own way from the top to a
    assert all(not m.comparable for row in found.branches for m in row)
    assert found.support[0][3].level == "na" and found.levels(0)["na"] == 2
    # another sequence set: nothing can be placed
    alien = side("alien", other.newick, {"a": "MMMM", "b": "MMMM", "c": "MM-M", "d": "MMMM"},
                 {"R": "MMMM", "P": "MMMM", "Q": "MMMM"})  # fmt: skip
    found = path_compare.compare(main, ["Node1", "Node2"], [alien], by_clade=True)
    assert found.followed[0].columns is None and not found.branches[0][0].comparable


def test_merged_rows_keep_each_ancestor_beside_itself() -> None:
    def followed(twins: dict[str, str], route: list[str]) -> path_compare.Followed:
        return path_compare.Followed(main_side(), twins, route, None)

    route = ["N1", "N2", "N3", "tip"]
    one = followed({"N1": "A", "N3": "C", "tip": "tip"}, ["A", "X", "C", "tip"])
    two = followed({"N1": "P", "N2": "Q", "N3": "R", "tip": "tip"}, ["P", "Q", "Y", "R", "tip"])
    assert path_compare.merged_rows(route, [one, two]) == [
        ["N1", "A", "P"],
        ["", "X", ""],  # only the first combination's tree has it
        ["N2", "", "Q"],
        ["", "", "Y"],
        ["N3", "C", "R"],
        ["tip", "tip", "tip"],
    ]
    assert path_compare.merged_rows(route, []) == [[n] for n in route]


def test_line_is_to_scale_until_stops_would_touch() -> None:
    names = ["N1", "N2", "N3", "tip"]
    levels = [{"all": 2, "none": 1}, {}, {"na": 4}]
    drawn = path_line.line(
        names, names, {"tip"}, ["u1", "u2", "u3", ""], [0.3, 0.0, 0.1], levels, ["a", "b", "c"]
    )
    xs = [s.x for s in drawn.stations]
    assert xs == sorted(xs) and xs[0] == path_line.PAD and xs[-1] == path_line.WIDTH - path_line.PAD
    first, second, third = drawn.stretches
    assert second.width >= path_line.MIN_STRETCH - 0.2 and drawn.stretched == 1
    assert abs(first.width / third.width - 3) < 0.05  # the real ones stay in proportion
    assert [(b.level, b.count) for b in first.bands] == [("all", 2), ("none", 1)]
    assert abs(first.bands[0].width - 2 * first.bands[1].width) < 0.1
    assert second.bands == [] and second.total == 0 and third.total == 4
    assert drawn.scale_width > 0 and "substitutions per site" in drawn.scale_label
    assert [level for level, _ in drawn.levels] == ["all", "none", "na"]
    assert drawn.stations[-1].extant and not drawn.stations[0].extant
    equal = path_line.line(
        names, names, {"tip"}, [""] * 4, [0.3, 0.0, 0.1], levels, [""] * 3, equal=True
    )
    widths = {round(s.width) for s in equal.stretches}
    assert len(widths) == 1 and equal.scale_width == 0 and equal.equal and not equal.stretched
    # names that would overlap go up a row
    close = path_line.line(
        ["A" * 20, "B" * 20], ["A" * 20, "B" * 20], set(), ["", ""], [0.0], [{}], [""], equal=True
    )
    assert {s.lane for s in close.stations} == {0}
    crowded = path_line.line(
        names,
        ["Ancestor-of-everything"] * 4,
        set(),
        [""] * 4,
        [1.0, 0.001, 0.001],
        levels,
        [""] * 3,
    )
    assert max(s.lane for s in crowded.stations) >= 1


def test_level_colours_are_the_same_in_the_stylesheet() -> None:
    """The chips' rules (site.css) repeat path_line.LEVEL_COLOURS, which the line, the
    bars and the viewers take from Python: two lists of the same colours drift."""
    css = (Path(path_line.__file__).parent / "static/phylobot/site.css").read_text()
    for level, colour in path_line.LEVEL_COLOURS.items():
        assert f".change.lvl-{level} {{ border-left-color: {colour}; }}" in css
    assert tuple(path_line.LEVEL_COLOURS) == path_compare.LEVELS
