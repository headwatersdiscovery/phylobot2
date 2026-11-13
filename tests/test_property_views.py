"""The ancestor viewer's property panel (`phylobot.property_views`).

Pure: every test here builds a summary dict as the pipeline writes it and checks what the
panel would draw. The panel's presence in the real viewer is covered by
`tests/test_ancestor_viewer.py`.
"""

from typing import Any

from phylobot import about, properties, property_store, property_views


def summary(**properties: dict[str, Any]) -> dict[str, Any]:
    """A stored summary.json, with the fields the pipeline always writes."""
    return {
        "node": "Node3",
        "gaps": "fitch",
        "draws": 1000,
        "seed": "1:Node3:fitch",
        "sampling": "per-column marginal posteriors, drawn independently",
        "weighted": False,
        "properties": properties,
    }


def entry(**fields: Any) -> dict[str, Any]:
    """One property's stored entry; `charge`-shaped unless overridden."""
    return {
        "label": "Net charge", "units": "e", "n": 1000, "level": 0.95, "weighted": False,
        "mean": -3.2, "sd": 4.0, "low": -11.0, "high": 4.6, "ml": -1.5,
        "width": 15.6, "shift": 1.7, "standardised_shift": 0.425, "cv": -1.25,
    } | fields  # fmt: skip


class TestNumbers:
    def test_integers_stay_integers(self) -> None:
        assert property_views.number(3.0) == "3"
        assert property_views.number(-12.0) == "-12"

    def test_significant_figures_not_decimals(self) -> None:
        # Properties span counts, charge and the hydrophobic moment, so a fixed number of
        # decimals is wrong for at least one of them.
        assert property_views.number(0.035612) == "0.0356"
        assert property_views.number(-3.2449) == "-3.24"

    def test_negative_zero_is_zero(self) -> None:
        assert property_views.number(-0.0) == "0"


class TestPanel:
    def test_it_reads_the_stored_fields(self) -> None:
        panel = property_views.of(summary(charge=entry()))
        assert panel is not None
        assert panel.draws == 1000 and panel.level == "95"
        assert panel.sampling.startswith("per-column marginal")
        assert panel.weighted is False
        row = panel.rows[0]
        assert row.key == "charge" and row.label == "Net charge" and row.units == "e"
        assert (row.mean, row.sd, row.low, row.high, row.ml) == ("-3.2", "4", "-11", "4.6", "-1.5")

    def test_the_offset_carries_its_sign(self) -> None:
        # The offset's direction is not a fixed bias, so the panel shows a signed distance
        # and never a correction.
        assert property_views.of(summary(c=entry())) is not None
        plus = property_views.of(summary(c=entry(standardised_shift=0.425)))
        minus = property_views.of(summary(c=entry(standardised_shift=-0.425)))
        assert plus is not None and minus is not None
        assert plus.rows[0].shift == "+0.42" and minus.rows[0].shift == "-0.42"
        assert plus.rows[0].shift_notable and minus.rows[0].shift_notable

    def test_a_small_offset_is_not_notable(self) -> None:
        panel = property_views.of(summary(c=entry(standardised_shift=0.04)))
        assert panel is not None and not panel.rows[0].shift_notable

    def test_several_properties_keep_their_order(self) -> None:
        panel = property_views.of(
            summary(
                charge=entry(),
                moment=entry(label="Hydrophobic moment", units=""),
                hydropathy_regions=entry(label="Hydrophobic regions", units=""),
            )
        )
        assert panel is not None
        assert [r.key for r in panel.rows] == ["charge", "moment", "hydropathy_regions"]

    def test_nothing_to_draw_is_none(self) -> None:
        assert property_views.of(None) is None
        assert property_views.of(summary()) is None

    def test_a_malformed_entry_is_skipped_not_drawn_as_zero(self) -> None:
        panel = property_views.of(summary(charge=entry(low=None), moment=entry()))
        assert panel is not None and [r.key for r in panel.rows] == ["moment"]

    def test_a_level_other_than_95(self) -> None:
        panel = property_views.of(summary(c=entry(level=0.8)))
        assert panel is not None and panel.level == "80"


def counts(values: list[float]) -> int:
    """Total count over a plot's bars: binning must never lose or duplicate a draw."""
    panel = property_views.of(summary(c=entry()), {"c": values})
    assert panel is not None and panel.rows[0].plot is not None
    return sum(bar.count for bar in panel.rows[0].plot.bars)


class TestHistogram:
    def test_every_draw_lands_in_exactly_one_bar(self) -> None:
        assert counts([float(i % 37) / 3 for i in range(1000)]) == 1000

    def test_the_tallest_bar_fills_the_plot(self) -> None:
        panel = property_views.of(summary(c=entry()), {"c": [1.0] * 50 + [2.0] * 10})
        assert panel is not None and (p := panel.rows[0].plot) is not None
        assert max(bar.height for bar in p.bars) == panel.baseline - panel.top
        assert all(panel.top <= bar.y <= panel.baseline for bar in p.bars)

    def test_a_count_property_gets_one_bar_per_integer(self) -> None:
        # hydropathy_regions is a count: half-integer bins would misrepresent it.
        panel = property_views.of(summary(c=entry()), {"c": [9.0, 10, 10, 11, 12, 12, 12]})
        assert panel is not None and (p := panel.rows[0].plot) is not None
        assert p.discrete and len(p.bars) == 4
        assert [bar.count for bar in p.bars] == [1, 2, 1, 3]

    def test_a_continuous_property_is_binned(self) -> None:
        panel = property_views.of(summary(c=entry()), {"c": [i / 100 for i in range(500)]})
        assert panel is not None and (p := panel.rows[0].plot) is not None
        assert not p.discrete and len(p.bars) == property_views.BINS

    def test_identical_draws_do_not_divide_by_zero(self) -> None:
        panel = property_views.of(summary(c=entry()), {"c": [3.0] * 20})
        assert panel is not None and (p := panel.rows[0].plot) is not None
        assert sum(bar.count for bar in p.bars) == 20

    def test_no_draws_means_no_plot_but_still_a_row(self) -> None:
        panel = property_views.of(summary(c=entry()), {})
        assert panel is not None and panel.rows[0].plot is None
        assert panel.rows[0].mean == "-3.2"  # the numbers still come from the summary


class TestMarkers:
    def test_mean_and_ml_are_always_drawn(self) -> None:
        panel = property_views.of(summary(c=entry()), {"c": [-3.0, -3.5, -2.0, -4.0]})
        assert panel is not None and (p := panel.rows[0].plot) is not None
        assert [m.kind for m in p.markers] == ["mean", "ml"]
        assert all(m.inside for m in p.markers)

    def test_altall_is_drawn_when_the_viewer_offers_one(self) -> None:
        panel = property_views.of(
            summary(c=entry()), {"c": [-3.0, -3.5, -2.0]}, {"c": -2.5}, altall_threshold="0.2"
        )
        assert panel is not None and (p := panel.rows[0].plot) is not None
        assert [m.kind for m in p.markers] == ["mean", "ml", "altall"]
        assert panel.rows[0].altall == "-2.5" and panel.altall_threshold == "0.2"

    def test_no_altall_leaves_the_caption_empty(self) -> None:
        panel = property_views.of(summary(c=entry()), {"c": [-3.0]}, altall_threshold="0.2")
        assert panel is not None and panel.altall_threshold == ""
        assert panel.rows[0].altall == ""

    def test_a_predictor_without_an_altall_still_draws_beside_one_that_has_it(self) -> None:
        """The web tier carries no torch, so f(AltAll) is computed only for the
        pure-Python properties (`views.property_panel`). A panel mixing the two must
        drop the marker on the predictor's row, not break or imply agreement."""
        panel = property_views.of(
            summary(charge=entry(), disorder=entry(label="Disordered fraction")),
            {"charge": [-3.0, -3.5], "disorder": [0.4, 0.5]},
            {"charge": -2.5},  # no "disorder": the viewer could not compute it
            altall_threshold="0.2",
        )
        assert panel is not None
        charge, disorder = panel.rows
        assert charge.altall == "-2.5" and disorder.altall == ""
        assert (p := disorder.plot) is not None
        assert [m.kind for m in p.markers] == ["mean", "ml"]
        assert panel.altall_threshold == "0.2"  # still captioned, for the row that has one

    def test_a_marker_outside_the_draws_is_still_on_the_axis(self) -> None:
        # f(ML) far from the posterior is exactly what the panel exists to show, so it
        # must never be clipped away: an absent marker would read as agreement.
        panel = property_views.of(summary(c=entry(ml=40.0)), {"c": [-3.0, -3.5, -2.0]})
        assert panel is not None and (p := panel.rows[0].plot) is not None
        ml = next(m for m in p.markers if m.kind == "ml")
        assert panel.left <= ml.x <= panel.right
        assert ml.x > max(bar.x + bar.width for bar in p.bars)


class TestMethods:
    def test_each_registered_property_has_a_method_note(self) -> None:
        """The panel's (?) affordance shows it, and the About page cites it."""
        for key in properties.REGISTRY:
            assert key in about.PROPERTY_METHODS, key
            method = about.PROPERTY_METHODS[key]
            assert method.name and len(method.note) > 80

    def test_the_cited_methods_carry_resolvable_references(self) -> None:
        for key in ("moment", "hydropathy_regions"):
            refs = about.PROPERTY_METHODS[key].refs
            assert refs, key
            assert all(ref.doi and ref.authors and ref.title for ref in refs)

    def test_charge_is_described_but_not_cited(self) -> None:
        """It is our own counting model, so citing anything would misattribute it."""
        assert not about.PROPERTY_METHODS["charge"].refs
        assert "pKa" in about.PROPERTY_METHODS["charge"].note

    def test_the_panel_carries_the_method_through(self) -> None:
        panel = property_views.of(summary(moment=entry()), methods=about.PROPERTY_METHODS)
        assert panel is not None and panel.rows[0].method is about.PROPERTY_METHODS["moment"]


class TestKeys:
    def test_the_paths_match_what_the_pipeline_publishes(self) -> None:
        # PROPERTIES writes into the combination's properties/ (main.nf publishes it there).
        assert (
            property_store.summaries_path("mafft/iqtree/LG+G4")
            == "mafft/iqtree/LG+G4/properties/summaries.json"
        )
        assert property_views.has(
            {"mafft/iqtree/LG+G4/properties/summaries.json"}, "mafft/iqtree/LG+G4"
        )
