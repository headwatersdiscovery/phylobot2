"""phylobot.properties: the property-posterior primitives."""

import math
import random
from collections.abc import Sequence

import pytest

from phylobot import properties
from phylobot.library import Source, Target


def target(template: str, dists: list[list[tuple[str, float]]]) -> Target:
    source = Source(combination=0, node="Node1", gap="parsimony")
    return Target(source=source, template=template, dists=dists)


def certain(template: str) -> Target:
    """An ancestor whose every column is certain, so every draw is the template."""
    return target(template, [[(c, 1.0)] for c in template])


class TestPropertyFunctions:
    def test_gaps_and_case_do_not_count(self) -> None:
        assert properties.residues("ka-R e") == "KARE"

    def test_net_charge_counts_sidechains(self) -> None:
        assert properties.net_charge("KR") == pytest.approx(2.0)
        assert properties.net_charge("DE") == pytest.approx(-2.0)
        assert properties.net_charge("KRDE") == pytest.approx(0.0)
        assert properties.net_charge("H") == pytest.approx(0.1)

    def test_net_charge_ignores_gaps_so_it_tracks_the_protein(self) -> None:
        assert properties.net_charge("K--R") == properties.net_charge("KR")

    def test_net_charge_is_extensive(self) -> None:
        assert properties.net_charge("KK") > properties.net_charge("K")

    def test_hydrophobic_moment_is_zero_without_a_full_window(self) -> None:
        assert properties.hydrophobic_moment("AAA") == 0.0

    def test_a_uniform_sequence_has_almost_no_moment(self) -> None:
        """The per-residue vectors spread around the helix and cancel, however
        hydrophobic the residues are. The moment measures asymmetry, not hydrophobicity."""
        assert properties.hydrophobic_moment("L" * 21) < 0.05

    def test_an_amphipathic_heptad_has_a_large_moment(self) -> None:
        """Hydrophobics at i, i+3, i+4 put them on one face, which is what the moment is
        for; the same residues spread uniformly do not."""
        heptad = properties.hydrophobic_moment("LKKLLKK" * 3)
        assert heptad > properties.hydrophobic_moment("L" * 21)
        assert heptad > 0.2

    def test_hydrophobic_moment_is_intensive(self) -> None:
        once = properties.hydrophobic_moment("LLKKKLLKKKLLKKKLLKKK")
        twice = properties.hydrophobic_moment("LLKKKLLKKKLLKKKLLKKK" * 2)
        assert twice == pytest.approx(once, rel=0.2)

    def test_hydropathy_regions_counts_runs_not_windows(self) -> None:
        assert properties.hydropathy_regions("K" * 20) == 0.0
        assert properties.hydropathy_regions("I" * 20) == 1.0
        assert properties.hydropathy_regions("I" * 10 + "K" * 10 + "I" * 10) == 2.0

    def test_one_substitution_can_remove_a_region(self) -> None:
        """Why this property is here: it is discontinuous. Alanine windows sit just above
        the threshold at 1.8, so a single glycine drops the mean to 1.49 and the region
        disappears. A draw differing at one site can land on the other side."""
        assert properties.hydropathy_regions("A" * 7) == 1.0
        assert properties.hydropathy_regions("AAAGAAA") == 0.0

    def test_one_substitution_can_split_a_region_in_two(self) -> None:
        assert properties.hydropathy_regions("A" * 21) == 1.0
        assert properties.hydropathy_regions("A" * 10 + "G" + "A" * 10) == 2.0

    def test_registry_keys_are_stable_and_callable(self) -> None:
        assert set(properties.REGISTRY) == {
            "charge",
            "moment",
            "hydropathy_regions",
            "disorder",
            "esm_likelihood",
        }
        for prop in properties.REGISTRY.values():
            assert prop.key and prop.label
            if prop.image == properties.TOOLS:  # the others need phylobot-ml (make ml-test)
                assert isinstance(prop.fn("ACDEFGHIKLMNPQRSTVWY"), float)

    def test_the_defaults_are_the_pure_python_ones(self) -> None:
        """The predictors are opt-in (owner, 2026-10-21); the first three stay on."""
        assert properties.DEFAULTS == ("charge", "moment", "hydropathy_regions")
        for key, prop in properties.REGISTRY.items():
            assert prop.default == (prop.image == properties.TOOLS), key

    def test_importing_the_registry_never_imports_torch(self) -> None:
        """The web tier and phylobot-tools carry no torch; the predictors import it only
        when they compute something."""
        import subprocess
        import sys

        code = "import sys, phylobot.properties; print('torch' in sys.modules)"
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        assert out.stdout.strip() == "False", out.stderr

    def test_a_run_needs_phylobot_ml_only_for_a_predictor(self) -> None:
        assert properties.image_for(["charge", "moment"]) == properties.TOOLS
        assert properties.image_for(["charge", "disorder"]) == properties.ML
        assert properties.image_for([]) == properties.TOOLS


class TestValues:
    def test_each_distinct_sequence_is_computed_once(self) -> None:
        seen: list[list[str]] = []

        def batch(sequences: Sequence[str]) -> list[float]:
            seen.append(list(sequences))
            return [float(len(s)) for s in sequences]

        prop = properties.Property("n", "n", lambda s: 0.0, "", "", batch=batch)
        assert prop.values(["MK", "M", "MK", "MK"]) == [2.0, 1.0, 2.0, 2.0]
        assert seen == [["MK", "M"]]

    def test_without_a_batch_it_uses_fn(self) -> None:
        calls: list[str] = []

        def fn(sequence: str) -> float:
            calls.append(sequence)
            return 1.0

        prop = properties.Property("n", "n", fn, "", "")
        assert prop.values(["A", "A", "C"]) == [1.0, 1.0, 1.0]
        assert calls == ["A", "C"]


class TestDrawing:
    def test_certain_columns_reproduce_the_template(self) -> None:
        drawn = properties.draw_sequences(certain("MKKLE"), draws=5, seed="s")
        assert drawn == ["MKKLE"] * 5

    def test_masked_columns_stay_gaps(self) -> None:
        t = target("M-K", [[("M", 1.0)], [("-", 1.0)], [("K", 0.5), ("R", 0.5)]])
        for sequence in properties.draw_sequences(t, draws=10, seed="s"):
            assert sequence[1] == "-"

    def test_the_seed_reproduces_the_draws(self) -> None:
        t = target("MK", [[("M", 1.0)], [("K", 0.5), ("R", 0.5)]])
        assert properties.draw_sequences(t, 20, seed="a") == properties.draw_sequences(t, 20, "a")

    def test_different_seeds_differ(self) -> None:
        t = target("MK", [[("M", 1.0)], [("K", 0.5), ("R", 0.5)]])
        assert properties.draw_sequences(t, 40, "a") != properties.draw_sequences(t, 40, "b")

    def test_draw_n_is_stable_however_many_were_asked_for(self) -> None:
        """As for library members: asking for more draws must not change the earlier ones,
        or a stored seed would not reproduce a published figure."""
        t = target("MK", [[("M", 1.0)], [("K", 0.5), ("R", 0.5)]])
        assert properties.draw_sequences(t, 50, "s")[:10] == properties.draw_sequences(t, 10, "s")

    def test_draws_follow_the_posterior(self) -> None:
        t = target("MK", [[("M", 1.0)], [("K", 0.8), ("R", 0.2)]])
        drawn = properties.draw_sequences(t, 2000, seed="s")
        share = sum(s[1] == "K" for s in drawn) / len(drawn)
        assert 0.75 < share < 0.85

    def test_draws_must_be_positive(self) -> None:
        with pytest.raises(ValueError, match="at least 1"):
            properties.draw_sequences(certain("MK"), draws=0)


class TestSummary:
    def test_summarises_a_known_sample(self) -> None:
        summary = properties.summarise([1.0, 2.0, 3.0, 4.0, 5.0], ml=3.0)
        assert summary.n == 5
        assert summary.mean == pytest.approx(3.0)
        assert summary.sd == pytest.approx(math.sqrt(2.5))
        assert summary.shift == pytest.approx(0.0)

    def test_interval_brackets_the_bulk_and_has_a_width(self) -> None:
        summary = properties.summarise([float(i) for i in range(101)], ml=50.0)
        assert summary.low == pytest.approx(2.5)
        assert summary.high == pytest.approx(97.5)
        assert summary.width == pytest.approx(95.0)

    def test_covers_asks_whether_the_truth_is_inside(self) -> None:
        summary = properties.summarise([float(i) for i in range(101)], ml=50.0)
        assert summary.covers(50.0)
        assert not summary.covers(-1.0)

    def test_shift_is_signed_and_standardised(self) -> None:
        summary = properties.summarise([0.0, 1.0, 2.0], ml=4.0)
        assert summary.shift == pytest.approx(3.0)
        assert summary.standardised_shift == pytest.approx(3.0 / summary.sd)

    def test_a_constant_property_has_no_spread_and_no_standardised_shift(self) -> None:
        """Such nodes are excluded from the gate; they must not raise or divide by zero."""
        summary = properties.summarise([2.0] * 10, ml=2.0)
        assert summary.sd == 0.0
        assert summary.cv == 0.0
        assert summary.standardised_shift == 0.0

    def test_single_draw_is_allowed(self) -> None:
        summary = properties.summarise([7.0], ml=7.0)
        assert summary.sd == 0.0
        assert summary.low == summary.high == 7.0

    def test_weights_move_the_mean_and_are_recorded(self) -> None:
        values = [0.0, 10.0]
        heavy_low = properties.summarise(values, ml=0.0, weights=[9.0, 1.0])
        assert heavy_low.mean == pytest.approx(1.0)
        assert heavy_low.weighted
        assert not properties.summarise(values, ml=0.0).weighted

    def test_rejects_bad_input(self) -> None:
        with pytest.raises(ValueError, match="no values"):
            properties.summarise([], ml=0.0)
        with pytest.raises(ValueError, match="between 0 and 1"):
            properties.summarise([1.0], ml=1.0, level=1.0)
        with pytest.raises(ValueError, match="differ in length"):
            properties.summarise([1.0, 2.0], ml=1.0, weights=[1.0])
        with pytest.raises(ValueError, match="sum to zero"):
            properties.summarise([1.0, 2.0], ml=1.0, weights=[0.0, 0.0])

    def test_summarise_target_keeps_ml_apart_from_the_draws(self) -> None:
        t = target("KK", [[("K", 1.0)], [("K", 0.5), ("D", 0.5)]])
        summary = properties.summarise_target(t, properties.REGISTRY["charge"], draws=400, seed="s")
        assert summary.ml == pytest.approx(2.0)  # the template is KK
        assert 0.0 < summary.mean < 2.0  # draws mix in the negative D
        assert summary.shift > 0.0


class TestDelta:
    def test_difference_adds_the_uncertainties(self) -> None:
        older = properties.summarise([0.0, 1.0, 2.0], ml=1.0)
        newer = properties.summarise([10.0, 11.0, 12.0], ml=11.0)
        change = properties.delta(older, newer)
        assert change.mean == pytest.approx(10.0)
        assert change.sd == pytest.approx(math.hypot(newer.sd, older.sd))
        assert change.ml == pytest.approx(10.0)
        assert change.excludes_zero
        assert change.sign == 1

    def test_an_overlapping_change_does_not_exclude_zero(self) -> None:
        a = properties.summarise([0.0, 5.0, 10.0], ml=5.0)
        b = properties.summarise([1.0, 6.0, 11.0], ml=6.0)
        assert not properties.delta(a, b).excludes_zero

    def test_sign_is_negative_when_the_property_falls(self) -> None:
        high = properties.summarise([10.0, 11.0], ml=10.0)
        low = properties.summarise([1.0, 2.0], ml=1.0)
        assert properties.delta(high, low).sign == -1

    def test_normal_quantile_matches_known_values(self) -> None:
        assert properties._normal_quantile(0.975) == pytest.approx(1.959964, abs=1e-5)
        assert properties._normal_quantile(0.5) == pytest.approx(0.0, abs=1e-9)
        assert properties._normal_quantile(0.001) == pytest.approx(-3.090232, abs=1e-5)


class TestComponents:
    def test_identical_groups_put_all_variance_within(self) -> None:
        split = properties.components({"a": [1.0, 2.0, 3.0], "b": [1.0, 2.0, 3.0]})
        assert split.between == pytest.approx(0.0, abs=1e-9)
        assert split.within > 0.0
        assert split.between_fraction == pytest.approx(0.0)

    def test_separated_groups_put_variance_between(self) -> None:
        split = properties.components({"a": [0.0, 0.1, 0.2], "b": [10.0, 10.1, 10.2]})
        assert split.between > split.within
        assert split.between_fraction > 0.9

    def test_a_negative_estimate_is_clamped_and_flagged(self) -> None:
        """Method of moments can go below zero when the true between variance is ~0."""
        split = properties.components({"a": [0.0, 10.0], "b": [0.0, 10.0], "c": [0.0, 10.0]})
        assert split.between == 0.0
        assert split.negative

    def test_counts_groups_and_draws(self) -> None:
        split = properties.components({"a": [1.0, 2.0], "b": [3.0, 4.0, 5.0]})
        assert (split.groups, split.n) == (2, 5)
        assert split.grand_mean == pytest.approx(3.0)

    def test_needs_two_usable_groups(self) -> None:
        with pytest.raises(ValueError, match="at least two groups"):
            properties.components({"a": [1.0, 2.0], "b": [3.0]})


class TestSignConsistency:
    def test_all_agreeing_is_total(self) -> None:
        falls = [properties.Delta(-1.0, 0.1, -1.2, -0.8, -1.0) for _ in range(5)]
        assert properties.sign_consistency(falls) == (5, 5, 1.0)

    def test_disagreement_reports_the_majority(self) -> None:
        deltas = [properties.Delta(1.0, 0.1, 0.8, 1.2, 1.0) for _ in range(3)]
        deltas += [properties.Delta(-1.0, 0.1, -1.2, -0.8, -1.0)]
        agreeing, total, fraction = properties.sign_consistency(deltas)
        assert (agreeing, total) == (3, 4)
        assert fraction == pytest.approx(0.75)

    def test_a_property_that_does_not_move_cannot_look_robust(self) -> None:
        flat = [properties.Delta(0.0, 0.1, -0.2, 0.2, 0.0) for _ in range(4)]
        assert properties.sign_consistency(flat) == (0, 4, 0.0)

    def test_no_deltas_is_empty_not_an_error(self) -> None:
        assert properties.sign_consistency([]) == (0, 0, 0.0)


class TestVaries:
    """`varies`, which decides which nodes a calibration check excludes. Guarded here
    because the scripts that use it are not run by `make test`, so a dangling reference
    from them went unnoticed once already."""

    def test_a_spread_property_varies(self) -> None:
        assert properties.summarise([1.0, 2.0, 3.0], ml=2.0).varies

    def test_an_identical_property_does_not(self) -> None:
        assert not properties.summarise([2.0] * 10, ml=2.0).varies

    def test_it_holds_for_a_property_that_crosses_zero(self) -> None:
        """cv is useless here (the mean is near zero), which is why `varies` exists."""
        summary = properties.summarise([-5.0, 0.1, 5.0], ml=0.0)
        assert summary.varies
        assert summary.cv > 10  # the number that made the cv floor meaningless


class TestPIT:
    def test_truth_below_everything_is_zero(self) -> None:
        assert properties.pit([1.0, 2.0, 3.0], truth=0.0) == 0.0

    def test_truth_above_everything_is_one(self) -> None:
        assert properties.pit([1.0, 2.0, 3.0], truth=9.0) == 1.0

    def test_truth_at_the_centre_is_a_half(self) -> None:
        assert properties.pit([0.0, 1.0, 2.0, 3.0], truth=1.5) == pytest.approx(0.5)

    def test_ties_are_split(self) -> None:
        assert properties.pit([1.0, 2.0, 2.0, 3.0], truth=2.0) == pytest.approx(0.5)

    def test_a_calibrated_sample_is_about_uniform(self) -> None:
        """The property PIT exists for: draws from the same law as the truth give PIT
        values spread evenly, so a histogram over many ancestors is flat."""
        rng = random.Random(7)
        draws = [rng.gauss(0, 1) for _ in range(400)]
        pits = [properties.pit(draws, rng.gauss(0, 1)) for _ in range(400)]
        assert 0.42 < sum(pits) / len(pits) < 0.58
        assert 0.2 < sum(1 for p in pits if p < 0.5) / len(pits) < 0.8

    def test_too_wide_a_posterior_humps_in_the_middle(self) -> None:
        rng = random.Random(7)
        wide = [rng.gauss(0, 5) for _ in range(400)]
        middle = sum(1 for _ in range(400) if 0.25 < properties.pit(wide, rng.gauss(0, 1)) < 0.75)
        assert middle / 400 > 0.75  # far more central than the uniform 0.5

    def test_empty_is_an_error(self) -> None:
        with pytest.raises(ValueError, match="no values"):
            properties.pit([], truth=0.0)


class TestCRPS:
    def test_a_point_forecast_is_the_absolute_error(self) -> None:
        assert properties.crps([3.0], truth=5.0) == pytest.approx(2.0)

    def test_it_is_zero_for_a_perfect_sharp_forecast(self) -> None:
        assert properties.crps([4.0] * 10, truth=4.0) == pytest.approx(0.0)

    def test_being_wrong_costs_more_than_being_right(self) -> None:
        near = properties.crps([0.0, 1.0, 2.0], truth=1.0)
        far = properties.crps([0.0, 1.0, 2.0], truth=9.0)
        assert far > near

    def test_it_punishes_vagueness_unlike_coverage(self) -> None:
        """A wide interval covers everything, so coverage cannot tell these apart; CRPS
        can, which is why the arms are compared on it."""
        rng = random.Random(3)
        sharp = [rng.gauss(0, 1) for _ in range(500)]
        vague = [rng.gauss(0, 20) for _ in range(500)]
        assert properties.crps(sharp, 0.0) < properties.crps(vague, 0.0)

    def test_empty_is_an_error(self) -> None:
        with pytest.raises(ValueError, match="no values"):
            properties.crps([], truth=0.0)
