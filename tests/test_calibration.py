"""Calibrated Posterior Probabilities (phylobot.calibration): temperature scaling keeps the
ML state, the covariates are what their definitions say, and T stays in bounds."""

import math

import pytest

from phylobot import calibration
from phylobot.calibration import Coefficients


def test_temperature_scaling_keeps_the_order() -> None:
    # With T = 1.5, (0.80, 0.15, 0.05) becomes (0.67, 0.22, 0.11).
    scaled = calibration.scale([("M", 0.80), ("K", 0.15), ("R", 0.05)], 1.5)
    assert [s for s, _ in scaled] == ["M", "K", "R"]
    assert [round(p, 2) for _, p in scaled] == [0.67, 0.22, 0.11]
    assert math.isclose(sum(p for _, p in scaled), 1.0)
    sharper = calibration.scale([("M", 0.80), ("K", 0.15), ("R", 0.05)], 0.5)
    assert sharper[0][1] > 0.8 and max(sharper, key=lambda sp: sp[1])[0] == "M"
    assert calibration.scale([("M", 0.6), ("K", 0.4)], 1.0) == [("M", 0.6), ("K", 0.4)]
    assert calibration.scale([("M", 1.0), ("K", 0.0)], 2.0) == [("M", 1.0), ("K", 0.0)]


def test_the_ml_state_never_changes_at_any_temperature() -> None:
    dist = [("A", 0.41), ("C", 0.39), ("D", 0.2)]
    for t in (0.1, 0.5, 1.0, 2.0, 10.0):
        scaled = calibration.scale(dist, t)
        assert max(scaled, key=lambda sp: sp[1])[0] == "A"


def test_temperatures_from_the_formula_are_bounded() -> None:
    zero = dict.fromkeys(calibration.KEYS, 0.0)
    # The shipped built-in is the fitted global T = 1.24 (no covariate terms), not identity.
    assert not calibration.BUILT_IN.identity and calibration.BUILT_IN.fitted
    assert math.isclose(calibration.BUILT_IN.temperature(zero), 1.24, abs_tol=1e-3)
    assert all(calibration.BUILT_IN.betas[k] == 0.0 for k in calibration.KEYS)
    assert Coefficients().identity and Coefficients().temperature(zero) == 1.0
    flatter = Coefficients(intercept=math.log(1.5))
    assert math.isclose(flatter.temperature(zero), 1.5) and not flatter.identity
    by_branch = Coefficients(betas={**zero, "branch_length": 0.5})
    assert math.isclose(by_branch.temperature({**zero, "branch_length": 1.0}), math.exp(0.5))
    assert Coefficients(intercept=5).temperature(zero) == calibration.T_MAX
    assert Coefficients(intercept=-5).temperature(zero) == calibration.T_MIN


def test_column_entropy_and_gaps() -> None:
    entropy, gaps = calibration.column_statistics(["MKA-", "MRA-", "MKC-", "MRDA"])
    assert entropy[0] == 0.0  # every sequence agrees
    assert math.isclose(entropy[1], 1 / math.log2(20))  # K/R half and half: one bit
    assert entropy[3] == 0.0 and gaps == [0, 0, 0, 0.75]
    with pytest.raises(ValueError, match="differ"):
        calibration.column_statistics(["MK", "M"])


def test_node_place_is_rooting_free() -> None:
    newick = "((a:0.1,b:0.3)Node2:0.2,(c:1.0,d:2.0)Node3:0.5,e:0.4)Node1;"
    top = calibration.node_place(newick, "Node1")
    assert top is not None
    assert math.isclose(top.branch_length, (0.2 + 0.5 + 0.4) / 3)
    assert math.isclose(top.tip_distance, 0.3)  # via Node2 to a
    deep = calibration.node_place(newick, "Node3")
    assert deep is not None
    assert math.isclose(deep.branch_length, (1.0 + 2.0 + 0.5) / 3)
    assert math.isclose(deep.tip_distance, 0.8)  # up through Node1 to a, not down to c
    assert calibration.node_place(newick, "Node9") is None
    assert calibration.node_place(newick, "a") is None  # extant, not an ancestor


def test_covariates_per_column() -> None:
    newick = "((a:0.1,b:0.1)Node2:0.1,c:0.1)Node1;"
    rates = "site\trate\tcategory\tcategory_rate\n1\t10.0\t4\t2.4\n2\t1.0\t2\t1.0\n"
    rows = calibration.covariates(["MK", "MR", "M-"], newick, "Node2", rates)
    assert rows is not None and len(rows) == 2
    assert math.isclose(rows[0]["branch_length"], 0.0, abs_tol=1e-12)  # 0.1, the reference
    assert math.isclose(rows[0]["tip_distance"], 0.0, abs_tol=1e-12)
    assert math.isclose(rows[0]["site_rate"], 1.0) and rows[1]["site_rate"] == 0.0
    assert rows[0]["entropy"] == 0.0 and rows[1]["entropy"] > 0
    assert rows[0]["gap_fraction"] == 0.0 and math.isclose(rows[1]["gap_fraction"], 1 / 3)
    # Without site_rates.tsv (no rate heterogeneity) every rate is 1.
    plain = calibration.covariates(["MK", "MR", "M-"], newick, "Node2", None)
    assert plain is not None and {r["site_rate"] for r in plain} == {0.0}
    assert calibration.covariates(["MK"], newick, "Node7", None) is None


def test_every_covariate_has_a_coefficient() -> None:
    assert set(calibration.BUILT_IN.betas) == set(calibration.KEYS)
    with pytest.raises(ValueError, match="columns differ"):
        calibration.apply([[("M", 1.0)]], [1.0, 2.0])
    s = calibration.summary([2.0, 1.0, 3.0, 4.0])
    assert s is not None and (s.low, s.median, s.high, s.sites) == (1.0, 2.5, 4.0, 4)
