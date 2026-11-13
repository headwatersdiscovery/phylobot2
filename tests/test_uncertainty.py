"""The decomposition (phylobot.uncertainty): hand-made cells whose
shares can be worked out on paper, and random ones for the identities that must hold
whatever the input."""

from collections.abc import Sequence
from itertools import permutations

import numpy as np
import pytest

from phylobot import uncertainty as u

TWO = ("alignment", "tree")


def one_hot(*states: str) -> np.ndarray:
    """[cells, 1 site, 20] point masses on the given residues."""
    out = np.zeros((len(states), 1, 20))
    for c, state in enumerate(states):
        out[c, 0, u.AMINO.index(state)] = 1.0
    return out


def residues_only(
    q: np.ndarray,
    weights: list[float],
    levels: Sequence[Sequence[str]],
    factors: tuple[str, ...] = TWO,
) -> u.Split:
    return u.categorical(np.zeros(q.shape[:2]), q, np.asarray(weights, float), levels, factors)


def random_cells(
    rng: np.random.Generator, cells: int, sites: int, gaps: bool = True
) -> tuple[np.ndarray, np.ndarray]:
    q = rng.dirichlet(np.full(20, 0.3), size=(cells, sites))
    g = (
        rng.uniform(size=(cells, sites)) * (rng.uniform(size=(cells, sites)) < 0.4)
        if gaps
        else np.zeros((cells, sites))
    )
    return g, q


def gini21(gap: np.ndarray, q: np.ndarray, w: np.ndarray) -> np.ndarray:
    """1 − Σp̄² over the 21 states of the weighted pool, computed directly."""
    w = w / w.sum()
    full = np.concatenate([(1 - gap)[..., None] * q, gap[..., None]], axis=-1)
    pooled = np.einsum("c,csk->sk", w, full)
    return np.asarray(1 - (pooled**2).sum(axis=-1))


class TestByHand:
    def test_one_cell_is_all_posterior(self) -> None:
        q = np.zeros((1, 1, 20))
        q[0, 0, u.AMINO.index("K")] = 0.5
        q[0, 0, u.AMINO.index("R")] = 0.5
        split = residues_only(q, [1], [("mafft", "iqtree")])
        assert split.all.total[0] == pytest.approx(0.5)
        assert split.all.within[0] == pytest.approx(0.5)
        assert all(s[0] == pytest.approx(0) for s in split.all.shares.values())

    def test_two_aligners_disagreeing_is_all_alignment(self) -> None:
        split = residues_only(
            one_hot("K", "R"), [1, 1], [("mafft", "iqtree"), ("muscle", "iqtree")]
        )
        shares = split.all.summed()
        assert shares["total"] == pytest.approx(0.5)  # two random draws differ half the time
        assert shares["within"] == pytest.approx(0)
        assert shares["alignment"] == pytest.approx(0.5)
        assert shares["tree"] == pytest.approx(0)  # one level: not varied

    def test_crossed_design_where_only_the_aligner_matters(self) -> None:
        levels = [
            ("mafft", "iqtree"),
            ("mafft", "raxml"),
            ("muscle", "iqtree"),
            ("muscle", "raxml"),
        ]
        split = residues_only(one_hot("K", "K", "R", "R"), [1] * 4, levels)
        shares = split.all.summed()
        assert shares["alignment"] == pytest.approx(0.5)
        assert shares["tree"] == pytest.approx(0)

    def test_a_pure_interaction_is_shared_evenly(self) -> None:
        # K where aligner and tree "agree", R where they don't: neither alone explains
        # anything (v(A) = v(T) = 0), together everything (v(A,T) = 0.5).
        levels = [
            ("mafft", "iqtree"),
            ("muscle", "raxml"),
            ("mafft", "raxml"),
            ("muscle", "iqtree"),
        ]
        split = residues_only(one_hot("K", "K", "R", "R"), [1] * 4, levels)
        shares = split.all.summed()
        assert shares["alignment"] == pytest.approx(0.25)
        assert shares["tree"] == pytest.approx(0.25)

    def test_weights_move_the_pool(self) -> None:
        # Weights 3:1 on K and R: 1 − (0.75² + 0.25²) = 0.375.
        split = residues_only(
            one_hot("K", "R"), [3, 1], [("mafft", "iqtree"), ("muscle", "iqtree")]
        )
        assert split.all.total[0] == pytest.approx(0.375)
        assert split.all.shares["alignment"][0] == pytest.approx(0.375)

    def test_a_zero_weight_cell_says_nothing(self) -> None:
        levels = [("mafft", "iqtree"), ("muscle", "iqtree"), ("clustalo", "iqtree")]
        split = residues_only(one_hot("K", "K", "R"), [1, 1, 0], levels)
        assert split.all.total[0] == pytest.approx(0)

    def test_a_site_nobody_speaks_for_is_zero(self) -> None:
        split = residues_only(one_hot("K", "R"), [0, 0], [("a", "t"), ("b", "t")])
        assert split.all.total[0] == 0 and split.all.within[0] == 0

    def test_cells_sharing_every_level_count_as_other(self) -> None:
        split = residues_only(one_hot("K", "R"), [1, 1], [("mafft", "iqtree"), ("mafft", "iqtree")])
        assert split.all.other[0] == pytest.approx(0.5)
        assert split.all.shares["alignment"][0] == pytest.approx(0)


class TestPresenceAndResidue:
    def test_the_example_that_breaks_subtraction_stays_non_negative(self) -> None:
        # K-or-R pooled with a gap. 21-state total 0.625.
        q = np.zeros((2, 1, 20))
        q[0, 0, u.AMINO.index("K")] = q[0, 0, u.AMINO.index("R")] = 0.5
        q[1, 0, u.AMINO.index("A")] = 1.0  # ignored: this cell puts a gap there
        gap = np.array([[0.0], [1.0]])
        split = u.categorical(gap, q, np.ones(2), [("mafft", "t"), ("muscle", "t")], TWO)
        assert split.all.total[0] == pytest.approx(0.625)
        assert split.presence.total[0] == pytest.approx(0.5)  # 2 · ½ · ½
        assert split.residue.total[0] == pytest.approx(0.125)  # ¼ both residues · ½ differ
        assert split.residue.within[0] == pytest.approx(0.125)
        assert split.residue.shares["alignment"][0] == pytest.approx(0)
        assert split.presence.shares["alignment"][0] == pytest.approx(0.5)

    def test_a_gap_in_every_cell_has_no_residue_part(self) -> None:
        gap = np.ones((2, 1))
        split = u.categorical(gap, one_hot("K", "R"), np.ones(2), [("a", "t"), ("b", "t")], TWO)
        assert split.all.total[0] == pytest.approx(0)

    def test_ml_gap_probability_is_posterior_presence_uncertainty(self) -> None:
        gap = np.array([[0.5]])
        split = u.categorical(gap, one_hot("K"), np.ones(1), [("a", "t")], TWO)
        assert split.presence.within[0] == pytest.approx(0.5)
        assert split.all.total[0] == pytest.approx(1 - 0.25 - 0.25)


class TestIdentities:
    @pytest.mark.parametrize("seed", range(5))
    def test_parts_add_up_to_the_21_state_variance(self, seed: int) -> None:
        rng = np.random.default_rng(seed)
        gap, q = random_cells(rng, 6, 30)
        w = rng.uniform(0.2, 2, size=6)
        levels = [(a, t, m) for a in "xy" for t in "pqr" for m in "m"]
        split = u.categorical(gap, q, w, levels, ("alignment", "tree", "model"))
        np.testing.assert_allclose(split.all.total, gini21(gap, q, w), atol=1e-10)
        for part in split.parts().values():
            between = part.total - part.within
            assert np.all(part.within >= -1e-12)
            assert all(np.all(s >= -1e-12) for s in part.shares.values())
            np.testing.assert_allclose(sum(part.shares.values()) + part.other, between, atol=1e-10)

    def test_shapley_matches_averaging_over_every_order(self) -> None:
        rng = np.random.default_rng(7)
        # Unbalanced: a crossed 2x3x2 design with three cells missing.
        levels = [(a, t, m) for a in "xy" for t in "pqr" for m in "LW"][:-3]
        factors = ("alignment", "tree", "model")
        means = rng.dirichlet(np.ones(4), size=(len(levels), 5))
        within = rng.uniform(size=(len(levels), 5))
        w = rng.uniform(0.5, 1.5, size=len(levels))
        found = u.decompose(means, within, w, levels, factors)

        def v(chosen: list[int]) -> np.ndarray:
            if not chosen:
                return np.zeros(5)
            ww = w / w.sum()
            pooled = np.einsum("c,csk->sk", ww, means)
            out = np.zeros(5)
            for key in {tuple(c[i] for i in chosen) for c in levels}:
                mine = [k for k, c in enumerate(levels) if tuple(c[i] for i in chosen) == key]
                weight = ww[mine].sum()
                mean = np.einsum("c,csk->sk", ww[mine], means[mine]) / weight
                out += weight * ((mean - pooled) ** 2).sum(-1)
            return out

        orders = list(permutations(range(3)))
        for i, factor in enumerate(factors):
            expected = np.zeros(5)
            for order in orders:
                before = sorted(order[: order.index(i)])
                expected += v(before + [i]) - v(before)
            np.testing.assert_allclose(found.shares[factor], expected / len(orders), atol=1e-10)

    def test_property_input_is_the_law_of_total_variance(self) -> None:
        # The engine is generic (§3.4): scalar draws in, means and variances.
        rng = np.random.default_rng(3)
        draws = rng.normal(size=(4, 200)) + np.array([[0], [0], [2], [2]])
        means = draws.mean(axis=1)[:, None, None]
        within = draws.var(axis=1)[:, None]
        levels = [
            ("mafft", "iqtree"),
            ("mafft", "raxml"),
            ("muscle", "iqtree"),
            ("muscle", "raxml"),
        ]
        found = u.decompose(means, within, np.ones(4), levels, TWO)
        assert found.total[0] == pytest.approx(draws.var())
        assert found.shares["alignment"][0] > 10 * found.shares["tree"][0]


class TestModesAndPairs:
    def test_most_probable_has_no_posterior_part(self) -> None:
        rng = np.random.default_rng(1)
        gap, q = random_cells(rng, 4, 20)
        point_gap, point = u.most_probable(gap, q)
        split = u.categorical(
            point_gap, point, np.ones(4), [("a", "t"), ("b", "t"), ("a", "s"), ("b", "s")], TWO
        )
        np.testing.assert_allclose(split.all.within, 0, atol=1e-12)
        assert set(np.unique(point_gap)) <= {0.0, 1.0}

    def test_pair_differences(self) -> None:
        q = one_hot("K", "R")
        pairs = u.pair_differences(np.zeros((2, 1)), q)
        np.testing.assert_allclose(pairs, [[0, 1], [1, 0]])
        # The diagonal is a cell's own uncertainty, summed over sites.
        rng = np.random.default_rng(2)
        gap, q = random_cells(rng, 3, 10)
        own = u.categorical(gap[:1], q[:1], np.ones(1), [("a", "t")], TWO).all.total.sum()
        assert u.pair_differences(gap, q)[0, 0] == pytest.approx(own)

    def test_level_states_say_what_each_choice_reconstructs(self) -> None:
        levels = [("mafft", "iqtree"), ("mafft", "raxml"), ("muscle", "iqtree")]
        gap = np.array([[0.0], [0.0], [1.0]])
        states = u.level_states(gap, one_hot("K", "K", "R"), np.ones(3), levels, 0)
        assert u.STATES[states["mafft"][0][0]] == "K" and states["mafft"][1][0] == pytest.approx(1)
        assert u.STATES[states["muscle"][0][0]] == u.GAP
