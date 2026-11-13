"""Measuring properties after a job has finished (`phylobot.properties_run`). Pure: no
storage, no Docker, no web."""

from typing import Any

import pytest

from phylobot import properties
from phylobot import properties_run as pr


def summary(*keys: str) -> dict[str, Any]:
    """A stored summary.json that measured `keys`."""
    return {"draws": 1000, "properties": {k: {"mean": 1.0} for k in keys}}


class TestWhatIsMissing:
    def test_a_job_that_measured_nothing_is_offered_everything(self) -> None:
        assert pr.missing(None) == [o.key for o in pr.offers(None)]
        assert set(pr.missing(None)) == set(properties.REGISTRY)

    def test_what_is_there_is_not_offered_again(self) -> None:
        assert "charge" not in pr.missing(summary("charge"))
        assert "moment" in pr.missing(summary("charge"))

    def test_nothing_missing_is_an_empty_offer(self) -> None:
        assert pr.missing(summary(*properties.REGISTRY)) == []

    def test_a_property_the_registry_no_longer_knows_is_ignored(self) -> None:
        """A summary from an older version can carry keys we have since dropped; nothing
        can be offered about them, and they must not look 'done'."""
        assert pr.measured(summary("charge", "gone_away")) == {"charge"}

    def test_the_done_ones_are_still_listed(self) -> None:
        """The list doubles as the answer to "what do I already have?", so a measured
        property stays on it, marked rather than hidden."""
        listed = pr.offers(summary("charge"))
        assert [o.key for o in listed] == [o.key for o in pr.offers(None)]
        done = next(o for o in listed if o.key == "charge")
        assert done.done and not done.offered and done.cost == 0


class TestOrder:
    def test_cheapest_tier_first_then_cheapest(self) -> None:
        listed = pr.offers(None)
        assert [o.tier for o in listed] == sorted(o.tier for o in listed)
        assert listed[0].key == "charge"  # the cheapest of the three
        keys = [o.key for o in listed]
        # Sliding windows make the moment the dearest pure-Python property; the
        # predictors cost more still, ESM-2 (a model call per draw) the most.
        assert keys.index("moment") > keys.index("hydropathy_regions")
        assert keys[-2:] == ["disorder", "esm_likelihood"]

    def test_the_order_does_not_change_between_requests(self) -> None:
        assert [o.key for o in pr.offers(None)] == [o.key for o in pr.offers(None)]


class TestCost:
    def test_it_scales_with_targets_draws_and_length(self) -> None:
        one = pr.cost_of("charge", targets=1, draws=1000, length=400)
        assert pr.cost_of("charge", targets=4, draws=1000, length=400) == pytest.approx(4 * one)
        assert pr.cost_of("charge", targets=1, draws=2000, length=400) == pytest.approx(2 * one)
        assert pr.cost_of("charge", targets=1, draws=1000, length=800) == pytest.approx(2 * one)

    def test_nothing_to_measure_costs_nothing(self) -> None:
        assert pr.cost_of("charge", targets=0, draws=1000, length=400) == 0
        assert pr.cost_of("not_a_property", targets=9, draws=1000, length=400) == 0
        assert pr.total_cost([], 10, 1000, 400) == 0

    def test_the_total_is_the_sum_of_the_ticked(self) -> None:
        keys = ["charge", "moment"]
        each = [pr.cost_of(k, 3, 1000, 400) for k in keys]
        assert pr.total_cost(keys, 3, 1000, 400) == pytest.approx(sum(each))


class TestScopes:
    ANALOGS = {"a/iqtree/LG": "Node1", "b/iqtree/LG": "Node7"}
    NODES = {"a/iqtree/LG": ["Node1", "Node2"], "b/iqtree/LG": ["Node7", "Node8"]}

    def expand(self, scope: str) -> list[pr.Target]:
        return pr.expand(scope, "a/iqtree/LG", "Node1", self.ANALOGS, self.NODES)

    def test_this_ancestor_only(self) -> None:
        assert self.expand("node") == [pr.Target("a/iqtree/LG", "Node1")]

    def test_this_ancestor_everywhere_it_was_recovered(self) -> None:
        """Each combination contributes its *own* name for the node: names belong to a
        tree, so Node1 here is not Node1 there."""
        assert self.expand("analogs") == [
            pr.Target("a/iqtree/LG", "Node1"),
            pr.Target("b/iqtree/LG", "Node7"),
        ]

    def test_every_ancestor_in_this_combination(self) -> None:
        assert self.expand("combination") == [
            pr.Target("a/iqtree/LG", "Node1"),
            pr.Target("a/iqtree/LG", "Node2"),
        ]

    def test_everything(self) -> None:
        targets = self.expand("job")
        assert len(targets) == 4
        assert pr.Target("b/iqtree/LG", "Node8") in targets
        # Never another combination's ancestors under this one's node names.
        assert pr.Target("b/iqtree/LG", "Node1") not in targets

    def test_scopes_grow(self) -> None:
        sizes = [len(self.expand(s)) for s in pr.SCOPES]
        assert sizes == sorted(sizes), dict(zip(pr.SCOPES, sizes, strict=True))

    def test_an_unknown_scope_falls_back(self) -> None:
        assert pr.scope_of("job") == "job"
        assert pr.scope_of("../../etc") == pr.DEFAULT_SCOPE
        assert pr.scope_of(None) == pr.DEFAULT_SCOPE


class TestTheChecklistsFields:
    """The checklist is part of the viewer's GET form, so what it submits decides what is
    ticked. An unticked box submits nothing, so "none" and "first load" are
    told apart by the field's presence, not by its contents."""

    OFFERED = ["charge", "moment", "hydropathy_regions"]

    def test_a_first_load_ticks_everything_missing(self) -> None:
        assert pr.chosen(None, self.OFFERED) == self.OFFERED

    def test_a_first_load_never_ticks_a_predictor(self) -> None:
        """Offered, but a minutes-per-ancestor model is never queued by default (owner,
        2026-10-21)."""
        offered = ["charge", "disorder", "esm_likelihood"]
        assert pr.chosen(None, offered) == ["charge"]
        assert pr.chosen(None, ["disorder", "esm_likelihood"]) == []
        assert pr.chosen(["", "esm_likelihood"], offered) == ["esm_likelihood"]

    def test_the_hidden_sentinel_alone_means_none(self) -> None:
        assert pr.chosen([""], self.OFFERED) == []

    def test_ticked_boxes_keep_checklist_order(self) -> None:
        assert pr.chosen(["", "hydropathy_regions", "charge"], self.OFFERED) == [
            "charge",
            "hydropathy_regions",
        ]

    def test_anything_not_offered_is_ignored(self) -> None:
        """Raw request text: a key that is not on offer (already measured, or invented)
        can never be smuggled into the run."""
        assert pr.chosen(["charge", "disorder", "../../etc"], self.OFFERED) == ["charge"]


class TestCostForAReader:
    def test_it_grows_with_the_work(self) -> None:
        costs = [pr.cost_of("moment", t, 1000, 400) for t in (1, 10, 100)]
        assert costs == sorted(costs) and costs[2] == pytest.approx(costs[0] * 100)

    def test_units_follow_the_size(self) -> None:
        """Hours, minutes and seconds, never raw CPU-seconds: "7986 CPU-s" is what the
        viewer showed before the predictors made the numbers big (owner, 2026-10-21)."""
        assert pr.cost_text(0.4) == "under a second"
        assert pr.cost_text(39) == "39 s"
        assert pr.cost_text(90) == "1 min 30 s"
        assert pr.cost_text(600) == "10 min"
        assert pr.cost_text(3600) == "1 h"
        assert pr.cost_text(7986) == "2 h 13 min"

    def test_a_row_carries_its_own_wording_and_a_unit_cost(self) -> None:
        """The list shows the duration; the dialog rescales `unit_cost` in the browser
        as the reader changes the draws or the scope."""
        row = next(
            o for o in pr.offers(None, targets=4, draws=1000, length=400) if o.key == "moment"
        )
        assert row.cost_text == pr.cost_text(row.cost)
        assert row.unit_cost == pytest.approx(pr.cost_of("moment", 1, 1000, 400))
        assert row.cost == pytest.approx(row.unit_cost * 4)

    def test_a_measured_row_costs_nothing_but_can_still_be_priced(self) -> None:
        """`cost` is what it would take to have it, which is nothing when you already
        do; `unit_cost` prices recomputing it, which costs the same as measuring did."""
        done = next(o for o in pr.offers(summary("moment"), length=400) if o.key == "moment")
        assert done.done and done.cost == 0
        assert done.unit_cost == pytest.approx(pr.cost_of("moment", 1, 1000, 400))

    def test_a_measured_row_carries_the_draws_it_was_measured_on(self) -> None:
        rows = pr.offers({"properties": {"moment": {"draws": 250}}})
        assert next(o for o in rows if o.key == "moment").measured_draws == 250

    def test_nothing_chosen_costs_nothing(self) -> None:
        assert pr.total_cost([], 100, 1000, 400) == 0
        assert pr.cost_text(0) == "—"
