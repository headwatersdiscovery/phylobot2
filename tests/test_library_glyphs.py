"""phylobot.library_glyphs: the library page's site map, sparkline and overview."""

import dataclasses

import pytest

from phylobot import ancestors, library_glyphs


def dists(*columns: dict[str, float]) -> list[list[tuple[str, float]]]:
    return [sorted(c.items(), key=lambda sp: (-sp[1], sp[0])) for c in columns]


COLUMNS = dists(
    {"M": 0.95, "L": 0.05},
    {"K": 0.4, "R": 0.35, "Q": 0.25},  # two plausible alternatives, not just one
    {"A": 1.0},  # masked in the template below
    {"Q": 0.7, "E": 0.3},
)


def test_site_map_skips_gaps_and_stacks_every_plausible_alternative() -> None:
    m = library_glyphs.site_map("MK-Q", COLUMNS, threshold=0.2)
    tiles = m.tiles
    assert m.sequence == "MKQ" and m.length == 3 and m.ambiguous == 2
    assert [(t.residue, t.column, t.state) for t in tiles] == [
        (1, 1, "M"),
        (2, 2, "K"),
        (3, 4, "Q"),
    ]
    # Shown at the threshold: the third state too, which ML and AltAll cannot cover.
    assert [a.state for a in tiles[1].shown] == ["R", "Q"]
    # Kept down to the lowest threshold a user may pick, so the page can move without
    # the server: M's alternative L (0.05) is carried, not shown.
    assert [a.state for a in tiles[0].alternatives] == ["L"] and tiles[0].shown == []
    assert (tiles[0].top, tiles[1].top) == (5, 35) and tiles[1].code == "R35 Q25"
    assert tiles[1].title == "residue 2 · column 2 · K 0.40" and tiles[1].pct == 40
    assert m.mean_pp == pytest.approx((0.95 + 0.4 + 0.7) / 3)
    assert library_glyphs.site_map("MK-Q", COLUMNS, 0.3).ambiguous == 2
    assert library_glyphs.site_map("MK-Q", COLUMNS, 0.36).ambiguous == 0
    assert library_glyphs.site_map("MK-Q", COLUMNS, 0.05).ambiguous == 3


def test_percentages_agree_with_the_threshold_on_its_grid() -> None:
    # 0.29 * 100 is 28.999... in floating point: rounded down naively it would be 28.
    assert library_glyphs.percent(0.29) == 29 and library_glyphs.percent(0.2) == 20
    assert library_glyphs.percent(0.199) == 19 and library_glyphs.percent(1.0) == 100
    for step in range(1, 51):
        t = step / 100
        for pp in (t - 0.004, t, t + 0.004, 0.29, 0.57, 0.199):
            assert (library_glyphs.percent(pp) >= library_glyphs.threshold_percent(t)) == (
                pp >= t - 1e-12
            ), (pp, t)


def test_colour_classes_cover_every_percentage() -> None:
    classes = library_glyphs.colour_classes("viridis")
    assert [i for i, _ in classes] == list(range(101))
    assert classes[40][1] == ancestors.colour(0.4, "viridis")


def test_site_map_without_posteriors_and_line_wrapping() -> None:
    m = library_glyphs.site_map("MK-Q", None, 0.2)
    assert m.mean_pp is None and m.ambiguous == 0
    assert all(t.pp is None and not t.uncertain(0.01) for t in m.tiles)
    long = library_glyphs.site_map("A" * 130, None, 0.2)
    assert [len(line) for line in long.lines] == [60, 60, 10]
    assert long.lines[1][0].residue == 61
    assert long.fasta("Node1").splitlines() == [">Node1", "A" * 60, "A" * 60, "A" * 10]
    with pytest.raises(ValueError, match="columns differ"):
        library_glyphs.site_map("MK", COLUMNS, 0.2)


def test_sparkline_ticks_follow_the_threshold() -> None:
    m = library_glyphs.site_map("MK-Q", COLUMNS, 0.2)
    spark = library_glyphs.sparkline(m.tiles, width=100, height=20)
    assert [top for _, top in spark.ticks] == [5, 35, 30] and not spark.empty
    assert spark.count(20) == 2 and spark.count(5) == 3 and spark.count(36) == 0
    assert spark.tick_path(20).count("M") == 2 and spark.tick_path(36) == ""
    assert spark.code.count(":") == 3
    for point in spark.points.split():
        x, y = map(float, point.split(","))
        assert 0 <= x <= 100 and 0 <= y <= 20
    big = library_glyphs.site_map("K" * 900, dists(*[{"K": 0.9, "R": 0.1}] * 900), 0.2)
    assert len(library_glyphs.sparkline(big.tiles, width=120).points.split()) == 120
    assert library_glyphs.sparkline([]).empty


# Two alignments of the same sequences (a, b, c), as in test_alignment_map: they agree on
# M and L; A puts the three K's in one column where B splits b's from a's and c's.
ALN_A = {"a": "MK-L", "b": "MKAL", "c": "MK-L"}
ALN_B = {"a": "M-K-L", "b": "MK-AL", "c": "M-K-L"}
CERTAIN = [[("M", 1.0)], [("K", 0.6), ("R", 0.4)], [("A", 1.0)], [("L", 0.9), ("I", 0.1)]]
CERTAIN_B = [[("M", 1.0)], [("K", 1.0)], [("K", 0.6), ("R", 0.4)], [("A", 1.0)], [("L", 1.0)]]


def entry(name: str, aligner: str, template: str, alignment: dict[str, str], dists=None):  # type: ignore[no-untyped-def]
    return library_glyphs.Entry(
        name, f"target-{name}", name, template, dists, aligner, aligner.upper(), alignment
    )


def test_layouts_offered_and_the_default() -> None:
    one = [entry("N1", "mafft", "MKAL", ALN_A, CERTAIN)]
    assert [k for k, _ in library_glyphs.alignment_options(one)] == ["mafft", "sequence"]
    assert library_glyphs.extant_options(one) == ["a", "b", "c"]
    assert library_glyphs.resolve_layout(one, "", []).alignment == "mafft"
    two = one + [entry("N2", "muscle", "MKKAL", ALN_B, CERTAIN_B)]
    assert [k for k, _ in library_glyphs.alignment_options(two)][0] == library_glyphs.OWN
    # With several alignments the default keeps them apart, never on one frame.
    assert library_glyphs.resolve_layout(two, "", []).alignment == library_glyphs.OWN


def test_resolve_layout_replaces_anything_unknown() -> None:
    two = [
        entry("N1", "mafft", "MKAL", ALN_A, CERTAIN),
        entry("N2", "muscle", "MKKAL", ALN_B, CERTAIN_B),
    ]
    for bad in ["../x", "<script>", "", "nobody"]:
        assert library_glyphs.resolve_layout(two, bad, [bad]).alignment == library_glyphs.OWN
    # Several extant sequences are kept on an alignment, in the order given, without repeats.
    layout = library_glyphs.resolve_layout(two, "mafft", ["c", "a", "c", "nobody"])
    assert layout.extant == ("c", "a")
    # On a sequence exactly one is the frame: the first chosen, else the first sequence.
    assert library_glyphs.resolve_layout(two, "sequence", ["b", "c"]).extant == ("b",)
    assert library_glyphs.resolve_layout(two, "sequence", []).extant == ("a",)


def test_own_columns_keep_each_alignment_apart() -> None:
    two = [
        entry("N1", "mafft", "MKAL", ALN_A, CERTAIN),
        entry("N2", "muscle", "MKKAL", ALN_B, CERTAIN_B),
    ]
    frames = library_glyphs.frames(two, library_glyphs.Layout(library_glyphs.OWN))
    assert [(f.title, f.units, len(f.rows)) for f in frames] == [
        ("MAFFT alignment · 4 columns", 4, 1), ("MUSCLE alignment · 5 columns", 5, 1),
    ]  # fmt: skip
    assert all(not f.tracks for f in frames)


def test_extant_sequences_are_rows_above_the_ancestors() -> None:
    one = [entry("N1", "mafft", "MKAL", ALN_A, CERTAIN)]
    (frame,) = library_glyphs.frames(one, library_glyphs.Layout("mafft", ("c", "a")))
    assert [(r.name, r.ancestor) for r in frame.rows] == [("c", False), ("a", False), ("N1", True)]
    assert [r.states for r in frame.rows] == ["MK-L", "MK-L", "MKAL"]
    # An extant row has no posteriors, so nothing to colour or call uncertain.
    assert set(frame.rows[0].pp) == {-1} and set(frame.rows[0].top) == {0}
    data = frame.data("coolwarm", 20)
    assert data["extant"] == 2 and data["scheme"] == "posterior" and data["cell"] == "fit"
    assert data["threshold"] == 20 and data["anchors"][2] == "target-N1"


def test_on_an_alignment_others_are_mapped_and_fade_where_aligners_disagree() -> None:
    two = [
        entry("N1", "mafft", "MKAL", ALN_A, CERTAIN),
        entry("N2", "muscle", "MKKAL", ALN_B, CERTAIN_B),
    ]
    (frame,) = library_glyphs.frames(two, library_glyphs.Layout("mafft"))
    first, mapped = frame.rows
    assert (first.name, mapped.name) == ("MAFFT N1", "MUSCLE N2")
    assert first.states == "MKAL" and not first.faded
    # B's columns 1 and 2 both map to A's column 1; the one agreeing better is shown there,
    # the other is not drawn, and the column is faded.
    assert mapped.states == "MKAL" and mapped.faded == [1]
    assert "share a column with another" in frame.note
    (track,) = frame.tracks
    assert track["label"] == "Aligners agree"
    assert track["series"][0]["values"] == [1.0, 0.667, 1.0, 1.0]  # rounded for the page


def test_on_a_sequence_columns_become_its_residues() -> None:
    two = [
        entry("N1", "mafft", "MKAL", ALN_A, CERTAIN),
        entry("N2", "muscle", "MKKAL", ALN_B, CERTAIN_B),
    ]
    (frame,) = library_glyphs.frames(two, library_glyphs.Layout("sequence", ("a",)))
    assert (frame.title, frame.units, frame.unit) == ("On a · 3 residues", 3, "residue")
    assert [r.name for r in frame.rows] == ["a", "MAFFT N1", "MUSCLE N2"]
    assert [r.states for r in frame.rows] == ["MKL", "MKL", "MKL"]
    # Residues the ancestors have where a has none are counted, not drawn: b's A after a's
    # second residue, and MUSCLE's extra K after the first (B splits the K column).
    inserted = next(t for t in frame.tracks if t["key"] == "inserted")
    assert inserted["series"][0]["values"] == [1, 1, 0]
    agreement = next(t for t in frame.tracks if t["key"] == "agreement")
    assert agreement["series"][0]["values"] == [1.0, 0.667, 1.0]


def test_targets_that_cannot_share_a_frame_keep_their_own_columns() -> None:
    other = {"x": "MK-L", "y": "MKAL"}  # another sequence set
    entries = [
        entry("N1", "mafft", "MKAL", ALN_A, CERTAIN),
        entry("Z", "elsewhere", "MKAL", other, CERTAIN),
    ]
    first, rest = library_glyphs.frames(entries, library_glyphs.Layout("mafft"))
    assert [r.name for r in first.rows] == ["N1"]
    assert rest.note == "Not the same sequences: shown on their own columns."
    _, rest = library_glyphs.frames(entries, library_glyphs.Layout("sequence", ("a",)))
    assert rest.note == "Their alignments do not contain a: shown on their own columns."


def test_rows_carry_every_alternative_and_the_strongest() -> None:
    one = [entry("N1", "mafft", "MKAL", ALN_A, CERTAIN)]
    (frame,) = library_glyphs.frames(one, library_glyphs.Layout("mafft"))
    (row,) = frame.rows
    assert row.pp == [100, 60, 100, 90] and row.top == [0, 40, 0, 10]
    assert row.alts == ["", "R40", "", "I10"]


def test_parse_zoom() -> None:
    assert library_glyphs.parse_zoom("fit") == "fit"
    assert library_glyphs.parse_zoom("3.5") == 3.5
    assert library_glyphs.parse_zoom("999") == library_glyphs.MAX_CELL
    assert library_glyphs.parse_zoom("0") == library_glyphs.MIN_CELL
    for junk in ["", "abc", None, "<script>"]:
        assert library_glyphs.parse_zoom(junk) == "fit"


def member(of: library_glyphs.Entry, name: str, sequence: str, cut: bool = False):  # type: ignore[no-untyped-def]
    """A member drawn on its ancestor's columns, as `library_views.member_viewers` makes it."""
    pairs = zip(of.template, sequence, strict=True)
    marked = frozenset(i for i, (a, b) in enumerate(pairs) if a != b)
    return dataclasses.replace(
        of, label=name, anchor=f"member-{name}", template=sequence, member=True,
        marked=marked, cut=cut,
    )  # fmt: skip


def test_members_are_rows_outlined_where_they_leave_the_ml_sequence() -> None:
    """L5c: a member's row shows its own states, coloured by its ancestor's posterior for
    them, and outlines its changes rather than the ancestor's uncertain sites."""
    n1 = entry("N1", "mafft", "MKAL", ALN_A, CERTAIN)
    rows = [member(n1, "ML", "MKAL"), member(n1, "ALT", "MRAI"), member(n1, "S1", "MRAL", cut=True)]
    (frame,) = library_glyphs.frames(rows, library_glyphs.Layout("mafft"))
    assert [r.name for r in frame.rows] == ["ML", "ALT", "S1"] and frame.members == 3
    ml, alt, cut = frame.rows
    assert ml.marked == [] and alt.marked == [1, 3] and cut.marked == [1]
    assert alt.pp == [100, 40, 100, 10]  # the posterior of the state it carries
    assert set(alt.top) == {0}  # no uncertainty outlines on a member
    assert (ml.dim, cut.dim) == (False, True)
    data = dataclasses.replace(frame, group="m").data("coolwarm", 20)
    assert data["member"] == [True] * 3 and data["marked"][1] == [1, 3]
    assert data["dim"] == [False, False, True] and data["anchors"][1] == "member-ALT"


def test_members_follow_their_ancestor_onto_another_layout() -> None:
    """Changes are carried through a mapping like everything else on the row, and a member's
    name is prefixed with its aligner once several combinations share one frame, so members
    from different aligners are told apart."""
    n1 = entry("N1", "mafft", "MKAL", ALN_A, CERTAIN)
    n2 = entry("N2", "muscle", "MKKAL", ALN_B, CERTAIN_B)
    rows = [member(n1, "A1", "MRAL"), member(n2, "B1", "MKRAL")]
    (frame,) = library_glyphs.frames(rows, library_glyphs.Layout("mafft"))
    assert [r.name for r in frame.rows] == ["MAFFT A1", "MUSCLE B1"]
    # B's column 3 (the R) maps to A's column 2.
    assert frame.rows[0].marked == [1] and frame.rows[1].marked == [1]
    (frame,) = library_glyphs.frames(rows, library_glyphs.Layout("sequence", ("a",)))
    assert [r.marked for r in frame.rows[1:]] == [[1], [1]]  # a's second residue


# --- Length against the extant sequences -------------------------------------------------

IME2_LENGTHS = [369, 512, 540, 601, 655, 712, 800]  # the shape of the Ime2 example


def test_extant_range_flags_lengths_outside_it() -> None:
    extant = library_glyphs.extant_range(IME2_LENGTHS)
    assert extant is not None and extant.check(369) is None and extant.check(800) is None
    long = extant.check(1052)
    assert long is not None and long.beyond == 1
    assert long.text.startswith("252 residues longer than the longest extant sequence.")
    assert "369–800 residues (7)" in long.text
    short = extant.check(300)
    assert short is not None and short.beyond == -1
    assert short.text.startswith("69 residues shorter than the shortest")
    assert library_glyphs.extant_range([]) is None


def test_extant_lengths_count_each_protein_once_across_alignments() -> None:
    mafft = {"a": "MK--Q", "b": "M-KLQ"}
    muscle = {"a": "MKQ--", "b": "MKLQ-", "c": "MA"}
    assert sorted(library_glyphs.extant_lengths([mafft, muscle])) == [2, 3, 4]


def test_length_chart_draws_every_ancestor_on_one_histogram() -> None:
    chart = library_glyphs.length_chart(
        IME2_LENGTHS,
        [
            ("Node3", "Node3 (original)", 1052),
            ("Node3 · fitch", "Node3 (fitch)", 684),
            ("Node4", "Node4 (original)", 684),
        ],
    )
    assert chart.stats == "Extant (n = 7): mean 598 · median 601 · sd 142"
    assert sum(int(b.title.split(": ")[1].split()[0]) for b in chart.bars) == 7
    # The axis stretches to the outlier, so nothing is pinned; ancestors at one length
    # share a line and a tooltip.
    assert [m.names for m in chart.marks] == [
        ["Node3 (fitch)", "Node4 (original)"],
        ["Node3 (original)"],
    ]
    inside, outside = chart.marks
    assert not inside.beyond
    assert inside.title == "Node3 (fitch): 684 residues\nNode4 (original): 684 residues"
    assert inside.label == "Node3 · fitch, Node4" and outside.label == "Node3"
    assert outside.beyond and "252 longer than the longest extant" in outside.title
    assert chart.left < inside.x < outside.x < chart.right
    assert all(chart.left <= x <= chart.right for x, _ in chart.x_ticks)
    # Labels sit above the plot, each line rising to its own; kept on the figure.
    assert all(m.label_y < m.top <= chart.top for m in chart.marks)
    assert outside.label_x + len(outside.label) * library_glyphs.LABEL_CHAR <= chart.width
    assert library_glyphs.length_chart([], [("Node1", "Node1", 10)]).empty


def test_length_chart_stacks_labels_that_would_clash() -> None:
    near = [(f"Node{i}", f"Node{i}", 500 + 6 * i) for i in range(1, 4)]  # ~8 px apart
    chart = library_glyphs.length_chart(IME2_LENGTHS, near)
    assert len({m.label_y for m in chart.marks}) == 3
    apart = [("Node1", "Node1", 400), ("Node2", "Node2", 700)]
    one_row = library_glyphs.length_chart(IME2_LENGTHS, apart)
    assert len({m.label_y for m in one_row.marks}) == 1
    assert one_row.top < chart.top  # the band grows only with the rows it needs


def test_label_rows_keep_lines_out_of_labels() -> None:
    # Side by side, the labels would overlap: one goes up a row, and its line rises
    # beside the other's label, not through it.
    rows = library_glyphs.label_rows([(100, 98, 138), (115, 113, 153)])
    assert rows[0] != rows[1]
    spans = [(100.0, 98.0, 138.0), (115.0, 113.0, 153.0)]
    for i, span in enumerate(spans):
        others = [s for j, s in enumerate(spans) if j != i]
        other_rows = [r for j, r in enumerate(rows) if j != i]
        assert library_glyphs.label_clashes(span, rows[i], others, other_rows) == 0
