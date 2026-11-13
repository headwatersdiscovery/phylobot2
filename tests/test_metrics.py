"""Calibration metrics: exact small cases, and the behaviour that matters on
synthetic calls with a known calibration."""

import math

import numpy as np
import pytest

from phylobot.evaluation import metrics


def synthetic(n: int, distort: float, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """PPs uniform on [0.05, 1]; right with probability pp**distort (1 = calibrated,
    above 1 = overconfident)."""
    rng = np.random.default_rng(seed)
    pp = rng.uniform(0.05, 1.0, n)
    return pp, (rng.uniform(size=n) < pp**distort).astype(float)


class TestExact:
    def test_ece_one_bin(self) -> None:
        assert metrics.ece([0.9, 0.9], [1, 0], n_bins=1) == pytest.approx(0.4)

    def test_ece_weights_bins_by_count(self) -> None:
        # bin [0.5, 0.6): one call, right (|1 - 0.55| = 0.45); bin [0.9, 1]: three, all right.
        value = metrics.ece([0.55, 0.95, 0.95, 0.95], [1, 1, 1, 1], n_bins=10)
        assert value == pytest.approx((0.45 + 3 * 0.05) / 4)

    def test_pp_of_one_lands_in_the_last_bin(self) -> None:
        bins = metrics.reliability([1.0, 0.0], [1, 0], n_bins=2)
        assert [b.count for b in bins] == [1, 1]

    def test_brier_and_log_score(self) -> None:
        uniform = np.full((1, 20), 0.05)
        assert metrics.brier(uniform, [3]) == pytest.approx(0.95**2 + 19 * 0.05**2)
        assert metrics.log_score(uniform, [3]) == pytest.approx(math.log(20))
        sure = np.eye(20)[[3]]
        assert metrics.brier(sure, [3]) == 0.0
        assert metrics.classwise_ece(np.eye(20)[[1, 2, 3]], [1, 2, 3]) == 0.0
        floored = -math.log(metrics.LOG_FLOOR)  # IQ-TREE's printed precision
        assert metrics.log_score(np.eye(20)[[1]], [2]) == pytest.approx(floored)

    def test_murphy_sums_to_the_brier_score(self) -> None:
        pp, y = synthetic(5000, 1.5)
        m = metrics.murphy(pp, y)
        assert m.reliability - m.resolution + m.uncertainty + m.within == pytest.approx(m.brier)
        # Constant within bins: the within-bin term vanishes.
        m = metrics.murphy([0.15, 0.15, 0.85, 0.85], [0, 1, 1, 1])
        assert m.within == pytest.approx(0.0)
        assert m.reliability - m.resolution + m.uncertainty == pytest.approx(m.brier)

    def test_bad_input(self) -> None:
        with pytest.raises(ValueError):
            metrics.ece([1.2], [1])
        with pytest.raises(ValueError):
            metrics.ece([0.5], [2])
        with pytest.raises(ValueError):
            metrics.ece([], [])
        with pytest.raises(ValueError):
            metrics.brier(np.eye(3), [0, 1, 5])
        with pytest.raises(ValueError, match="strategy"):
            metrics.bin_edges([0.5], 10, "log")


class TestBehaviour:
    def test_calibrated_calls_score_near_zero(self) -> None:
        pp, y = synthetic(20000, 1.0)
        assert metrics.ece(pp, y) < 0.015
        assert metrics.ece(pp, y, strategy=metrics.COUNT) < 0.015

    def test_overconfident_calls_do_not(self) -> None:
        pp, y = synthetic(20000, 2.0)
        assert metrics.ece(pp, y) > 0.1
        bins = metrics.reliability(pp, y)
        assert all(b.accuracy < b.mean_pp for b in bins if b.count and b.mean_pp < 0.95)

    def test_equal_count_bins(self) -> None:
        pp, y = synthetic(1000, 1.0)
        counts = [b.count for b in metrics.reliability(pp, y, 5, metrics.COUNT)]
        assert counts == [200] * 5

    def test_intervals_cover_the_estimate_and_widen_with_clustering(self) -> None:
        pp, y = synthetic(2000, 1.0)
        bins = metrics.reliability(pp, y, n_bins=5, resamples=200)
        for b in bins:
            assert b.interval is not None
            assert b.interval[0] <= b.accuracy <= b.interval[1]
        # Ten clusters whose calls share an outcome: far less information than 2000 sites.
        clustered = np.repeat(np.arange(10), 200)
        y_shared = np.repeat(np.random.default_rng(1).uniform(size=10) < 0.5, 200).astype(float)
        p_flat = np.full(2000, 0.5)
        lo_i, hi_i = metrics.ece_interval(p_flat, y_shared, resamples=200)
        lo_c, hi_c = metrics.ece_interval(p_flat, y_shared, groups=clustered, resamples=200)
        assert hi_c - lo_c > hi_i - lo_i

    def test_bootstrap_is_seeded(self) -> None:
        pp, y = synthetic(500, 1.0)
        assert metrics.ece_interval(pp, y, resamples=50) == metrics.ece_interval(
            pp, y, resamples=50
        )

    def test_strata(self) -> None:
        s = metrics.strata(np.arange(100), 4)
        assert np.bincount(s).tolist() == [25, 25, 25, 25]
        assert s[0] == 0 and s[-1] == 3
