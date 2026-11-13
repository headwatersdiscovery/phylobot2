"""Where an ancestor's uncertainty comes from: the decomposition.

Each *cell* (a combination and gap version) gives, per site, a mean vector and a
within-cell variance. Pooled over the cells with weights, the variance splits exactly
into the mean within-cell variance (the posterior's own uncertainty) and the variance
between cells, and that between term splits over the factors that make a cell (aligner,
tree, model, gap version) by Shapley values: each factor's marginal contribution
averaged over every order in which the factors could be added (§1.2). The shares are
non-negative and add up to the between term, whatever the design's balance.

`decompose` knows nothing about sequences: a categorical state passes its one-hot mean
(the distribution) and 1 − Σp²; a property would pass the mean and variance of its draws
(§3.4). `categorical` is the sequence case, split into presence (residue or gap?) and
residue (which one?) so that both parts stay variances (§1.3). Pure numpy; no file
reading here.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import combinations as subsets
from math import factorial

import numpy as np

# IQ-TREE's column order in ancestors_posteriors.tsv; the gap is the 21st state.
AMINO = "ARNDCQEGHILKMFPSTWYV"
GAP = "-"
STATES = AMINO + GAP

FACTORS = ("alignment", "tree", "model", "gaps")
TINY = 1e-12  # weights below this are no weight (a site no cell speaks for)


@dataclass(frozen=True)
class Components:
    """One decomposition, per site: total = within + Σ shares + other."""

    total: np.ndarray  # [sites]
    within: np.ndarray  # [sites]: the mean within-cell variance
    shares: dict[str, np.ndarray]  # factor -> [sites]: its Shapley share of the between term
    other: np.ndarray  # [sites]: between cells that share every level (zero when cells are unique)

    def scaled(self, by: np.ndarray) -> Components:
        return Components(
            self.total * by,
            self.within * by,
            {f: s * by for f, s in self.shares.items()},
            self.other * by,
        )

    def __add__(self, other: Components) -> Components:
        return Components(
            self.total + other.total,
            self.within + other.within,
            {f: s + other.shares[f] for f, s in self.shares.items()},
            self.other + other.other,
        )

    def summed(self) -> dict[str, float]:
        """Every component summed over sites (a node's totals, in residues)."""
        out = {"total": float(self.total.sum()), "within": float(self.within.sum())}
        out |= {f: float(s.sum()) for f, s in self.shares.items()}
        out["other"] = float(self.other.sum())
        return out


def _normalised(weights: np.ndarray) -> np.ndarray:
    """Weights [cells, sites] made to sum to 1 at every site that has any."""
    weights = np.clip(weights, 0.0, None)
    total = weights.sum(axis=0)
    return np.asarray(np.divide(weights, total, out=np.zeros_like(weights), where=total > TINY))


def _grouping(levels: Sequence[Sequence[str]], chosen: Sequence[int]) -> np.ndarray:
    """[groups, cells] membership, cells grouped by their levels of the chosen factors."""
    ids: dict[tuple[str, ...], int] = {}
    member = [ids.setdefault(tuple(cell[i] for i in chosen), len(ids)) for cell in levels]
    grouping = np.zeros((len(ids), len(levels)))
    grouping[member, np.arange(len(levels))] = 1.0
    return grouping


def _between(
    means: np.ndarray, w: np.ndarray, grouping: np.ndarray, pooled_sq: np.ndarray
) -> np.ndarray:
    """v(S) per site: Σ_g W_g ‖m_g − p̄‖², written as Σ_g ‖Σ_{c∈g} w_c μ_c‖² / W_g − ‖p̄‖²
    (the weights sum to 1, so the cross terms cancel)."""
    group_weight = grouping @ w  # [groups, sites]
    group_sum = np.einsum("gc,csk->gsk", grouping, w[..., None] * means)
    squared = np.einsum("gsk,gsk->gs", group_sum, group_sum)
    spread = np.divide(squared, group_weight, out=np.zeros_like(squared), where=group_weight > TINY)
    return np.asarray(np.clip(spread.sum(axis=0) - pooled_sq, 0.0, None))


def decompose(
    means: np.ndarray,
    within: np.ndarray,
    weights: np.ndarray,
    levels: Sequence[Sequence[str]],
    factors: Sequence[str] = FACTORS,
) -> Components:
    """Total, within and each factor's Shapley share of the between-cell variance.

    `means` is [cells, sites, dims], `within` [cells, sites], `weights` [cells, sites]
    (or [cells], the same at every site; normalised per site, and a site with no weight
    is all zeros), `levels` one tuple per cell naming its level of each factor.
    """
    cells = means.shape[0]
    if within.shape != means.shape[:2] or len(levels) != cells:
        raise ValueError("means, within and levels disagree on the cells or sites")
    if any(len(cell) != len(factors) for cell in levels):
        raise ValueError("every cell needs one level per factor")
    weights = np.asarray(weights, dtype=float)
    if weights.ndim == 1:
        weights = np.repeat(weights[:, None], means.shape[1], axis=1)
    w = _normalised(weights)
    pooled = np.einsum("cs,csk->sk", w, means)
    pooled_sq = np.einsum("sk,sk->s", pooled, pooled)
    within_mean = np.einsum("cs,cs->s", w, within)
    between = _between(means, w, np.eye(cells), pooled_sq)  # every cell its own group

    n = len(factors)
    value: dict[frozenset[int], np.ndarray] = {frozenset(): np.zeros(means.shape[1])}
    for size in range(1, n + 1):
        for chosen in subsets(range(n), size):
            value[frozenset(chosen)] = _between(means, w, _grouping(levels, chosen), pooled_sq)
    shares = {}
    for i, factor in enumerate(factors):
        share = np.zeros(means.shape[1])
        rest = [j for j in range(n) if j != i]
        for size in range(n):
            coefficient = factorial(size) * factorial(n - size - 1) / factorial(n)
            for chosen in subsets(rest, size):
                s = frozenset(chosen)
                share += coefficient * (value[s | {i}] - value[s])
        shares[factor] = share
    explained = value[frozenset(range(n))]
    other = np.clip(between - explained, 0.0, None)
    return Components(within_mean + between, within_mean, shares, other)


@dataclass(frozen=True)
class Split:
    """The sequence decomposition in its two parts (§1.3); `all` is their sum, which is
    the 21-state decomposition's total exactly."""

    presence: Components
    residue: Components
    gap: np.ndarray = field(repr=False)  # [sites]: the pooled P(gap)

    @property
    def all(self) -> Components:
        return self.presence + self.residue

    def parts(self) -> dict[str, Components]:
        return {"all": self.all, "presence": self.presence, "residue": self.residue}


def categorical(
    gap: np.ndarray,
    residues: np.ndarray,
    weights: np.ndarray,
    levels: Sequence[Sequence[str]],
    factors: Sequence[str] = FACTORS,
) -> Split:
    """Decompose ancestral states: `gap` [cells, sites] is each cell's P(gap), `residues`
    [cells, sites, 20] its residue distribution *given* a residue (rows summing to 1),
    `weights` [cells] (zero for a cell that lacks the node).

    presence: the 2-state vector (g, 1 − g), variance 2g(1 − g) within a cell.
    residue: the conditional distribution q, weighted by w(1 − g) (a cell with a gap there
    says nothing about which residue), scaled by (1 − ḡ)², the chance two draws from the
    pool are both residues. presence + residue = 1 − Σp̄² over the 21 states.
    """
    gap = np.clip(gap, 0.0, 1.0)
    weights = np.asarray(weights, dtype=float)
    per_site = np.repeat(weights[:, None], gap.shape[1], axis=1) if weights.ndim == 1 else weights
    w = _normalised(per_site)
    pooled_gap = np.einsum("cs,cs->s", w, gap)
    two = np.stack([gap, 1.0 - gap], axis=-1)
    presence = decompose(two, 2.0 * gap * (1.0 - gap), per_site, levels, factors)
    gini = 1.0 - np.einsum("csk,csk->cs", residues, residues)
    residue = decompose(residues, gini, per_site * (1.0 - gap), levels, factors)
    return Split(presence, residue.scaled((1.0 - pooled_gap) ** 2), pooled_gap)


def most_probable(gap: np.ndarray, residues: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Each cell's point estimate as distributions (§1.4): a gap where P(gap) > 0.5, and
    the one-hot of its most probable residue (the first, on a tie, in AMINO order)."""
    point = np.zeros_like(residues)
    best = residues.argmax(axis=-1)
    np.put_along_axis(point, best[..., None], 1.0, axis=-1)
    return (gap > 0.5).astype(float), point


def pair_differences(gap: np.ndarray, residues: np.ndarray) -> np.ndarray:
    """[cells, cells]: the expected positions at which a draw from one cell and a draw
    from the other differ, Σ_sites 1 − p_c·p_d over the 21 states. The diagonal is a
    cell's own posterior uncertainty."""
    present = 1.0 - gap
    same_residue = np.einsum("csk,dsk->cds", residues, residues)
    same = (
        gap[:, None, :] * gap[None, :, :] + present[:, None, :] * present[None, :, :] * same_residue
    )
    return np.asarray((1.0 - same).sum(axis=-1))


def level_states(
    gap: np.ndarray,
    residues: np.ndarray,
    weights: np.ndarray,
    levels: Sequence[Sequence[str]],
    factor: int,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """What each level of one factor says at every site: {level: (state index into
    STATES, its probability)} from the level's cells pooled with their weights ([cells]
    or [cells, sites]). Levels whose cells all have zero weight are left out; at a site
    where they have none, the state is the gap with probability 0."""
    weights = np.asarray(weights, dtype=float)
    if weights.ndim == 1:
        weights = np.repeat(weights[:, None], gap.shape[1], axis=1)
    full = np.concatenate([(1.0 - gap)[..., None] * residues, gap[..., None]], axis=-1)
    out = {}
    for level in dict.fromkeys(cell[factor] for cell in levels):
        mine = np.array([cell[factor] == level for cell in levels])
        if not (weights[mine] > TINY).any():
            continue
        pooled = np.einsum("cs,csk->sk", _normalised(weights[mine]), full[mine])
        best = pooled.argmax(axis=-1)
        best[pooled.sum(axis=-1) <= TINY] = len(STATES) - 1
        out[level] = (best, np.take_along_axis(pooled, best[:, None], axis=-1)[:, 0])
    return out
