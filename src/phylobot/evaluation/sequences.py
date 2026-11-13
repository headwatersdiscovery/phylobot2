"""Ancestral sequence accuracy and calibration."""

from collections.abc import Iterable
from dataclasses import dataclass


def identity(inferred: str, truth: str) -> float:
    """Fraction of residues identical between two ungapped sequences of equal length.

    Only valid when no indels were simulated; indel scenarios will need a pairwise
    alignment first.
    """
    a, b = inferred.replace("-", ""), truth.replace("-", "")
    if len(a) != len(b):
        raise ValueError(f"length mismatch ({len(a)} vs {len(b)}); align before comparing")
    if not a:
        raise ValueError("empty sequences")
    return sum(x.upper() == y.upper() for x, y in zip(a, b, strict=True)) / len(a)


@dataclass(frozen=True)
class CalibrationBin:
    low: float
    high: float
    count: int
    mean_pp: float
    accuracy: float


def calibration(
    calls: Iterable[tuple[float, bool]], edges: tuple[float, ...] = (0.0, 0.5, 0.8, 0.9, 0.95, 1.0)
) -> list[CalibrationBin]:
    """Bin site calls by posterior probability and report how often each bin was right.

    Well-calibrated output has accuracy close to mean_pp in every bin.
    """
    bins: list[list[tuple[float, bool]]] = [[] for _ in edges[:-1]]
    for pp, correct in calls:
        for i in range(len(edges) - 1):
            last = i == len(edges) - 2
            if edges[i] <= pp < edges[i + 1] or (last and pp == edges[-1]):
                bins[i].append((pp, correct))
                break
    return [
        CalibrationBin(
            low=edges[i],
            high=edges[i + 1],
            count=len(b),
            mean_pp=sum(p for p, _ in b) / len(b) if b else 0.0,
            accuracy=sum(c for _, c in b) / len(b) if b else 0.0,
        )
        for i, b in enumerate(bins)
    ]
