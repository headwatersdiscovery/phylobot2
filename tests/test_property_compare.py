"""One ancestor's properties across combinations (`phylobot.property_compare`).

Pure: every test builds candidates by hand and checks the geometry and numbers the page
would draw. The page itself is covered by `tests/test_views.py`.
"""

from typing import Any

from phylobot import property_compare as pc


def candidate(key: str, values: list[float], *, shown: bool = True, **fields: Any) -> pc.Candidate:
    summary = {
        "properties": {
            "charge": {
                "label": "Net charge", "units": "e", "level": 0.95,
                "mean": sum(values) / len(values), "sd": 1.0,
                "low": min(values), "high": max(values), "ml": values[0],
                "standardised_shift": 0.3,
            } | fields,
        }
    }  # fmt: skip
    return pc.Candidate(
        key=key, node="Node1", label=key, aligner=key.split("/")[0],
        shown=shown, values={"charge": values}, summary=summary,
    )  # fmt: skip


def panel_of(*candidates: pc.Candidate) -> pc.Panel:
    panel = pc.panel("charge", "Net charge", "e", list(candidates))
    assert panel is not None
    return panel


class TestSharedAxis:
    def test_every_series_is_drawn_on_the_same_axis(self) -> None:
        """The point of the page: two distributions can only be compared on one scale."""
        panel = panel_of(
            candidate("a/iqtree/LG", [1.0, 2.0, 3.0]),
            candidate("b/iqtree/LG", [10.0, 11.0, 12.0]),
        )
        assert panel.axis.low <= 1.0 and panel.axis.high >= 12.0
        # Both outlines span the same x range, so their positions mean the same thing.
        xs = [[float(p.split(",")[0]) for p in s.outline.split()] for s in panel.series]
        assert min(xs[0]) == min(xs[1]) and max(xs[0]) == max(xs[1])

    def test_hiding_a_series_does_not_move_the_axis(self) -> None:
        """Otherwise two readings of the same page would not be comparable — the one
        thing this page exists to make possible."""
        both = panel_of(
            candidate("a/iqtree/LG", [1.0, 2.0, 3.0]),
            candidate("b/iqtree/LG", [10.0, 11.0, 12.0]),
        )
        one_off = panel_of(
            candidate("a/iqtree/LG", [1.0, 2.0, 3.0]),
            candidate("b/iqtree/LG", [10.0, 11.0, 12.0], shown=False),
        )
        assert (both.axis.low, both.axis.high) == (one_off.axis.low, one_off.axis.high)
        assert both.series[0].outline == one_off.series[0].outline

    def test_a_hidden_series_stays_in_the_table_but_not_the_plot(self) -> None:
        panel = panel_of(
            candidate("a/iqtree/LG", [1.0, 2.0]),
            candidate("b/iqtree/LG", [3.0, 4.0], shown=False),
        )
        assert [s.key for s in panel.series] == ["a/iqtree/LG", "b/iqtree/LG"]
        assert [s.key for s in panel.drawn] == ["a/iqtree/LG"]


class TestDistributions:
    def test_series_are_proportions_not_counts(self) -> None:
        """Combinations may have different numbers of draws; raw counts would make the
        bigger one look like a taller distribution rather than a different one."""
        panel = panel_of(
            candidate("a/iqtree/LG", [5.0] * 100),
            candidate("b/iqtree/LG", [5.0] * 10),
        )
        tops = [min(float(p.split(",")[1]) for p in s.outline.split()) for s in panel.series]
        assert tops[0] == tops[1]  # same height: both put everything in one bin

    def test_the_area_closes_to_the_baseline(self) -> None:
        panel = panel_of(candidate("a/iqtree/LG", [1.0, 2.0, 3.0]))
        s = panel.series[0]
        ys = [float(p.split(",")[1]) for p in s.area.split()]
        assert ys[0] == panel.axis.baseline and ys[-1] == panel.axis.baseline

    def test_a_count_property_gets_integer_bins(self) -> None:
        panel = panel_of(candidate("a/iqtree/LG", [3.0, 4.0, 4.0, 5.0]))
        assert panel.axis.count == 3 and panel.axis.low == 2.5

    def test_identical_draws_do_not_divide_by_zero(self) -> None:
        panel = panel_of(candidate("a/iqtree/LG", [2.0] * 5))
        assert panel.series[0].outline


class TestNumbers:
    def test_the_lead_states_the_spread_of_means(self) -> None:
        panel = panel_of(
            candidate("a/iqtree/LG", [1.0, 1.0], mean=1.0),
            candidate("b/iqtree/LG", [4.0, 4.0], mean=4.0),
        )
        assert panel.spread == "1 to 4"

    def test_one_series_has_no_spread_to_report(self) -> None:
        assert panel_of(candidate("a/iqtree/LG", [1.0, 2.0])).spread == ""

    def test_a_hidden_series_is_left_out_of_the_spread(self) -> None:
        panel = panel_of(
            candidate("a/iqtree/LG", [1.0, 1.0], mean=1.0),
            candidate("b/iqtree/LG", [9.0, 9.0], mean=9.0, shown=False),
        )
        assert panel.spread == ""

    def test_the_offset_keeps_its_sign(self) -> None:
        panel = panel_of(candidate("a/iqtree/LG", [1.0, 2.0], standardised_shift=-0.4))
        assert panel.series[0].shift == "-0.40"

    def test_colours_cycle_rather_than_running_out(self) -> None:
        panel = panel_of(*(candidate(f"a{i}/iqtree/LG", [float(i), 1.0]) for i in range(7)))
        assert [s.colour for s in panel.series] == [1, 2, 3, 4, 5, 1, 2]
        assert all(1 <= s.colour <= pc.COLOURS for s in panel.series)


class TestNothingToDraw:
    def test_a_property_nobody_measured_is_no_panel(self) -> None:
        assert pc.panel("moment", "Moment", "", [candidate("a/iqtree/LG", [1.0])]) is None

    def test_no_candidates_is_no_panel(self) -> None:
        assert pc.panel("charge", "Net charge", "e", []) is None


class TestMeanMarks:
    def test_each_shown_series_gets_one(self) -> None:
        """Labelled with the combination, not the value: the question these lines answer
        is which combination sits where, and the values are in the table below."""
        panel = panel_of(
            candidate("a/iqtree/LG", [1.0, 2.0], mean=1.5),
            candidate("b/iqtree/LG", [3.0, 4.0], mean=3.5),
        )
        assert [m.text for m in panel.means] == ["a/iqtree/LG", "b/iqtree/LG"]
        assert all(panel.axis.left <= m.x <= panel.axis.right for m in panel.means)

    def test_a_hidden_series_has_no_mark(self) -> None:
        panel = panel_of(
            candidate("a/iqtree/LG", [1.0, 2.0], mean=1.5),
            candidate("b/iqtree/LG", [3.0, 4.0], mean=3.5, shown=False),
        )
        assert [m.text for m in panel.means] == ["a/iqtree/LG"]

    def test_they_can_be_switched_off(self) -> None:
        panel = pc.panel(
            "charge", "Net charge", "e",
            [candidate("a/iqtree/LG", [1.0, 2.0], mean=1.5)], show_means=False,
        )  # fmt: skip
        assert panel is not None and panel.means == []
        assert panel.axis.headroom == 0  # and the plot takes the room back

    def test_crowded_labels_stack_instead_of_shifting(self) -> None:
        """Means cluster — that is the usual finding — so labels at nearly the same x
        are normal. Each must stay over its own line, so they stack."""
        panel = panel_of(
            *(candidate(f"a{i}/iqtree/LG", [0.0, 100.0], mean=50.0 + i * 0.1) for i in range(4))
        )
        assert sorted(m.lane for m in panel.means) == [0, 1, 2, 3]
        # Every label sits over its own line (x unchanged), not pushed aside.
        assert all(abs(m.label_x - m.x) < 1 for m in panel.means)
        # Lanes are stacked upwards, and the plot leaves room for them.
        ys = sorted(m.y for m in panel.means)
        assert all(b - a >= pc.LANE_HEIGHT - 0.01 for a, b in zip(ys, ys[1:], strict=False))
        assert panel.axis.headroom >= 4 * pc.LANE_HEIGHT
        assert min(ys) > 0  # not off the top of the viewBox

    def test_far_apart_labels_share_one_lane(self) -> None:
        panel = panel_of(
            candidate("a/iqtree/LG", [0.0, 100.0], mean=5.0),
            candidate("b/iqtree/LG", [0.0, 100.0], mean=95.0),
        )
        assert [m.lane for m in panel.means] == [0, 0]

    def test_a_label_at_the_edge_stays_inside_the_plot(self) -> None:
        panel = panel_of(candidate("a/iqtree/LG", [0.0, 100.0], mean=100.0))
        mark = panel.means[0]
        assert mark.label_x < mark.x  # pulled in, so the text is not clipped
        assert mark.label_x <= panel.axis.right


class TestPalettes:
    def test_three_are_offered(self) -> None:
        assert list(pc.PALETTES) == ["bright", "okabe", "ink"]
        assert all(label and note for label, note in pc.PALETTES.values())

    def test_unknown_palettes_fall_back(self) -> None:
        """The value comes from the query string, so it is raw request text."""
        assert pc.palette_of("okabe") == "okabe"
        assert pc.palette_of("<script>") == pc.DEFAULT_PALETTE
        assert pc.palette_of(None) == pc.DEFAULT_PALETTE
        assert pc.palette_of(["ink"]) == pc.DEFAULT_PALETTE


class TestLabelLanes:
    """Mean labels stack rather than slide: a label that moved sideways would stop
    naming the line beneath it (`property_compare.lanes`)."""

    def test_a_label_stays_over_its_own_line(self) -> None:
        placed = pc.lanes([(0, 400.0, "MAFFT · LG+G4"), (1, 405.0, "MUSCLE · WAG+G4")])
        for _, x, _, label_x, _ in placed:
            assert label_x == x

    def test_overlapping_labels_go_to_different_lanes(self) -> None:
        placed = pc.lanes([(0, 400.0, "MAFFT · LG+G4"), (1, 405.0, "MUSCLE · WAG+G4")])
        assert sorted(lane for *_, lane in placed) == [0, 1]

    def test_far_apart_labels_share_a_lane(self) -> None:
        placed = pc.lanes([(0, 80.0, "a"), (1, 600.0, "b")])
        assert [lane for *_, lane in placed] == [0, 0]

    def test_a_label_at_the_edge_is_pulled_inside(self) -> None:
        ((_, _, _, label_x, _),) = pc.lanes([(0, 0.0, "a long label here")])
        assert label_x > 0

    def test_the_margins_are_not_clobbered_between_marks(self) -> None:
        """Each label is clamped to the *plot's* edges. The label's own extent once
        shared those names, so every mark after the first was clamped against the
        previous label and slid away from its line."""
        marks = [(0, 241.5, "posterior mean"), (1, 246.9, "f(ML)"), (2, 24.4, "f(AltAll)")]
        placed = pc.lanes(marks, left=8.0, right=8.0, width=460)
        by_colour = {colour: label_x for colour, _, _, label_x, _ in placed}
        assert by_colour[0] == 241.5 and by_colour[1] == 246.9
