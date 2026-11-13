"""Calibration metrics, over the records of
`phylobot.evaluation.records` or any (probability, outcome) pairs.

Two views of the same calls:
- **top-label**: the ML state's PP against whether it was right. Reliability bins, ECE,
  and the Murphy decomposition of its Brier score. This is what a user reads off a site.
- **full posterior**: every state's probability against the true state. Classwise ECE
  (Kull et al. 2019), the multi-class Brier score and the log score. These see the
  alternatives that AltAll and property draws use.

Sites on one tree are not independent, so intervals come from a **cluster bootstrap**:
resample whole groups (a node, a held-out tip, a replicate), never single sites, when
`groups` is given.
"""

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

WIDTH, COUNT = "width", "count"
# IQ-TREE prints posteriors to 5 decimals, so a printed 0 means "below 0.000005": the true
# state's PP is floored there, never at an arbitrary tiny number that would let one site
# dominate a mean (one such site in Thioredoxin's ESR scored -690 nats under a 1e-300
# floor, -6.4 a site over its sequence; 2026-11-02).
LOG_FLOOR = 5e-6

Floats = NDArray[np.float64]


@dataclass(frozen=True)
class Bin:
    low: float
    high: float
    count: int
    mean_pp: float
    accuracy: float
    interval: tuple[float, float] | None = None  # bootstrap interval for the accuracy


@dataclass(frozen=True)
class Murphy:
    """Brier = reliability - resolution + uncertainty (+ within-bin term, which is 0 when
    every PP in a bin is the same and small otherwise; kept so the sum is exact)."""

    brier: float
    reliability: float
    resolution: float
    uncertainty: float
    within: float


def _arrays(pp: ArrayLike, correct: ArrayLike) -> tuple[Floats, Floats]:
    p = np.asarray(pp, dtype=float)
    y = np.asarray(correct, dtype=float)
    if p.shape != y.shape or p.ndim != 1:
        raise ValueError("pp and correct must be 1-d and the same length")
    if p.size == 0:
        raise ValueError("no calls")
    if np.any((p < 0) | (p > 1)) or not np.all(np.isin(y, (0.0, 1.0))):
        raise ValueError("pp must lie in [0, 1] and correct be 0 or 1")
    return p, y


def bin_edges(pp: ArrayLike, n_bins: int = 10, strategy: str = WIDTH) -> Floats:
    """Equal-width bins on [0, 1], or equal-count bins from the PPs' quantiles (ties can
    leave a bin with more than its share; duplicate edges are merged)."""
    if strategy == WIDTH:
        return np.linspace(0.0, 1.0, n_bins + 1)
    if strategy == COUNT:
        q = np.quantile(np.asarray(pp, dtype=float), np.linspace(0.0, 1.0, n_bins + 1))
        q[0], q[-1] = 0.0, 1.0
        return np.unique(q)
    raise ValueError(f"unknown binning strategy {strategy!r}")


def _assign(p: Floats, edges: Floats) -> NDArray[np.intp]:
    """Bin index per call; the last bin is closed on the right."""
    return np.clip(np.searchsorted(edges, p, side="right") - 1, 0, len(edges) - 2)


def _bootstrap(
    statistic: Callable[[NDArray[np.intp]], Floats],
    n: int,
    groups: ArrayLike | None,
    resamples: int,
    seed: int,
) -> Floats:
    """`statistic` (of row indices) over cluster-bootstrap resamples: one row per draw."""
    rng = np.random.default_rng(seed)
    if groups is None:
        members = [np.array([i]) for i in range(n)]
    else:
        g = np.asarray(groups)
        if g.shape != (n,):
            raise ValueError("groups must give one group per call")
        members = [np.flatnonzero(g == k) for k in np.unique(g)]
    draws = []
    for _ in range(resamples):
        chosen = rng.integers(0, len(members), len(members))
        draws.append(statistic(np.concatenate([members[c] for c in chosen])))
    return np.array(draws)


def reliability(
    pp: ArrayLike,
    correct: ArrayLike,
    n_bins: int = 10,
    strategy: str = WIDTH,
    groups: ArrayLike | None = None,
    resamples: int = 0,
    level: float = 0.95,
    seed: int = 1,
) -> list[Bin]:
    """Reliability-diagram bins: the mean PP and the fraction right in each. With
    `resamples`, a percentile interval for each bin's accuracy (empty bins in a resample
    are left out of that bin's interval)."""
    p, y = _arrays(pp, correct)
    edges = bin_edges(p, n_bins, strategy)
    k = len(edges) - 1
    which = _assign(p, edges)

    def accuracies(rows: NDArray[np.intp]) -> Floats:
        counts = np.bincount(which[rows], minlength=k)
        hits = np.bincount(which[rows], weights=y[rows], minlength=k)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(counts > 0, hits / np.maximum(counts, 1), np.nan)

    intervals: list[tuple[float, float] | None] = [None] * k
    if resamples:
        draws = _bootstrap(accuracies, len(p), groups, resamples, seed)
        tail = (1 - level) / 2 * 100
        for b in range(k):
            column = draws[:, b][~np.isnan(draws[:, b])]
            if column.size:
                lo, hi = np.percentile(column, [tail, 100 - tail])
                intervals[b] = (float(lo), float(hi))
    bins = []
    for b in range(k):
        mask = which == b
        n = int(mask.sum())
        bins.append(
            Bin(
                low=float(edges[b]),
                high=float(edges[b + 1]),
                count=n,
                mean_pp=float(p[mask].mean()) if n else float("nan"),
                accuracy=float(y[mask].mean()) if n else float("nan"),
                interval=intervals[b] if n else None,
            )
        )
    return bins


def ece(pp: ArrayLike, correct: ArrayLike, n_bins: int = 10, strategy: str = WIDTH) -> float:
    """Expected calibration error: the count-weighted mean |accuracy - mean PP| over bins."""
    p, y = _arrays(pp, correct)
    edges = bin_edges(p, n_bins, strategy)
    which = _assign(p, edges)
    total = 0.0
    for b in range(len(edges) - 1):
        mask = which == b
        if mask.any():
            total += mask.sum() * abs(y[mask].mean() - p[mask].mean())
    return float(total / len(p))


def ece_interval(
    pp: ArrayLike,
    correct: ArrayLike,
    groups: ArrayLike | None = None,
    n_bins: int = 10,
    strategy: str = WIDTH,
    resamples: int = 1000,
    level: float = 0.95,
    seed: int = 1,
) -> tuple[float, float]:
    """A cluster-bootstrap percentile interval for the ECE."""
    p, y = _arrays(pp, correct)

    def statistic(rows: NDArray[np.intp]) -> Floats:
        return np.array([ece(p[rows], y[rows], n_bins, strategy)])

    draws = _bootstrap(statistic, len(p), groups, resamples, seed)[:, 0]
    tail = (1 - level) / 2 * 100
    lo, hi = np.percentile(draws, [tail, 100 - tail])
    return float(lo), float(hi)


def _posteriors(
    posteriors: ArrayLike, true: Sequence[int] | ArrayLike
) -> tuple[Floats, NDArray[np.intp]]:
    q = np.asarray(posteriors, dtype=float)
    t = np.asarray(true, dtype=np.intp)
    if q.ndim != 2 or t.shape != (q.shape[0],):
        raise ValueError("posteriors must be calls x states, true one state index per call")
    if q.shape[0] == 0:
        raise ValueError("no calls")
    if np.any((t < 0) | (t >= q.shape[1])):
        raise ValueError("a true state index is out of range")
    return q, t


def classwise_ece(
    posteriors: ArrayLike, true: Sequence[int] | ArrayLike, n_bins: int = 10
) -> float:
    """Mean over states of the ECE of P(state) against "the state is the true one"
    (equal-width bins). States never predicted above 0 and never true add 0."""
    q, t = _posteriors(posteriors, true)
    onehot = np.zeros_like(q)
    onehot[np.arange(len(t)), t] = 1.0
    return float(np.mean([ece(q[:, k], onehot[:, k], n_bins) for k in range(q.shape[1])]))


def brier(posteriors: ArrayLike, true: Sequence[int] | ArrayLike) -> float:
    """Multi-class Brier score: mean over calls of sum over states of (p - outcome)^2."""
    q, t = _posteriors(posteriors, true)
    onehot = np.zeros_like(q)
    onehot[np.arange(len(t)), t] = 1.0
    return float(np.mean(np.sum((q - onehot) ** 2, axis=1)))


def log_score(posteriors: ArrayLike, true: Sequence[int] | ArrayLike) -> float:
    """Mean negative log probability of the true state (lower is better)."""
    q, t = _posteriors(posteriors, true)
    return float(-np.mean(np.log(np.maximum(q[np.arange(len(t)), t], LOG_FLOOR))))


def murphy(pp: ArrayLike, correct: ArrayLike, n_bins: int = 10, strategy: str = WIDTH) -> Murphy:
    """The top-label Brier score, (pp - correct)^2 averaged, decomposed over bins."""
    p, y = _arrays(pp, correct)
    edges = bin_edges(p, n_bins, strategy)
    which = _assign(p, edges)
    base = y.mean()
    rel = res = within = 0.0
    for b in range(len(edges) - 1):
        mask = which == b
        if not mask.any():
            continue
        w = mask.sum() / len(p)
        f, o = p[mask].mean(), y[mask].mean()
        rel += w * (f - o) ** 2
        res += w * (o - base) ** 2
        # (p - y)^2 = (p - f)^2 + 2(p - f)(f - y) + (f - y)^2 within the bin
        within += w * (np.mean((p[mask] - f) ** 2) + 2 * np.mean((p[mask] - f) * (f - y[mask])))
    return Murphy(
        brier=float(np.mean((p - y) ** 2)),
        reliability=float(rel),
        resolution=float(res),
        uncertainty=float(base * (1 - base)),
        within=float(within),
    )


def strata(values: ArrayLike, n: int = 4) -> NDArray[np.intp]:
    """A stratum index per call, by quantiles of a covariate (0 = lowest)."""
    v = np.asarray(values, dtype=float)
    edges = np.unique(np.quantile(v, np.linspace(0.0, 1.0, n + 1)))
    return np.clip(np.searchsorted(edges, v, side="right") - 1, 0, len(edges) - 2)
