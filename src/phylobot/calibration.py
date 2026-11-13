"""Calibrated posterior probabilities: covariate-dependent temperature scaling.

IQ-TREE's marginal
posteriors are taken as they are and reshaped site by site: each state's probability is
raised to the power 1/T and the distribution renormalised. T = 1 changes nothing, T > 1
flattens (confidence moves to the alternatives), T < 1 sharpens. The order of the states
never changes, so **the ML ancestor is the same** whatever T is; only the stated
confidence moves.

T depends on the site's situation through a linear formula in covariates that any job
has, without knowing the truth:

    ln T = intercept + sum over covariates of (coefficient x covariate)

Each covariate is 0 at a reference site (`Covariate.reference`), so the intercept is ln T
there. The coefficients are fitted once, on a benchmark with known ground truth, and
shipped with their provenance (`Coefficients.source`). The Settings page overrides them
(`site_settings`), and each user may set their own (My calibration).

Pure Python over published files (the alignment, tree.nwk, site_rates.tsv), applied per
request like the gap versions are per file: nothing is written into a job, and
`ancestors_posteriors.tsv` stays IQ-TREE's.
"""

import heapq
import math
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import dendropy

from phylobot.trees import load_tree

RAW, CALIBRATED = "raw", "calibrated"
# T is kept inside these bounds whatever the coefficients say, so a mistyped coefficient
# cannot turn every site into a certainty (T -> 0) or into a uniform distribution.
T_MIN, T_MAX = 0.1, 10.0
FLOOR = 0.001  # lengths and rates below this count as this (logarithms of 0)
AMINO_STATES = 20
NOT_RESIDUES = frozenset("-.?XxBbZzJj*")


@dataclass(frozen=True)
class Covariate:
    key: str
    label: str
    meaning: str  # what it measures and in what unit, for the Settings page
    reference: str  # the value at which it is 0


COVARIATES = (
    Covariate(
        "branch_length",
        "Branch length at the node",
        "log10 of the mean length of the branches meeting at the node (substitutions "
        "per site), relative to 0.1: +1 at 1.0, -1 at 0.01. Long branches are where "
        "posteriors are known to be overconfident.",
        "0.1 substitutions per site",
    ),
    Covariate(
        "tip_distance",
        "Distance to the nearest sequence",
        "log10 of the path length from the node to its nearest extant sequence, "
        "relative to 0.1. Deep ancestors are far from any data.",
        "0.1 substitutions per site",
    ),
    Covariate(
        "site_rate",
        "Site rate",
        "log10 of IQ-TREE's posterior mean rate for the column (site_rates.tsv; 1 for "
        "every column without rate heterogeneity): +1 at ten times the mean rate.",
        "the mean rate (1)",
    ),
    Covariate(
        "entropy",
        "Column entropy",
        "Shannon entropy of the extant residues in the column (gaps left out), as a "
        "fraction of the most possible (log2 20 bits): 0 when every sequence agrees, "
        "1 when all twenty amino acids are equally common.",
        "a conserved column",
    ),
    Covariate(
        "gap_fraction",
        "Gap fraction",
        "The fraction of extant sequences with a gap in the column, 0 to 1.",
        "no gaps",
    ),
)
KEYS = tuple(c.key for c in COVARIATES)


@dataclass(frozen=True)
class Coefficients:
    intercept: float = 0.0
    betas: Mapping[str, float] = field(default_factory=lambda: dict.fromkeys(KEYS, 0.0))
    source: str = "placeholder"  # where they come from, shown with every calibrated view
    fitted: bool = False  # fitted on the benchmark (not placeholders, not hand-set)

    @property
    def identity(self) -> bool:
        """T = 1 at every site: calibrated equals raw."""
        return self.intercept == 0 and all(self.betas.get(k, 0.0) == 0 for k in KEYS)

    def temperature(self, x: Mapping[str, float]) -> float:
        log_t = self.intercept + sum(self.betas.get(k, 0.0) * x[k] for k in KEYS)
        return min(max(math.exp(log_t), T_MIN), T_MAX)


# The built-in coefficients: a global MAP-state temperature fitted on a benchmark, with no
# covariate terms. ln(1.24) = 0.2151. Users may override the intercept and every covariate
# on My calibration; the site may override them in Settings.
BUILT_IN = Coefficients(
    intercept=0.2151,  # ln(1.24): global MAP-state temperature
    source="Fitted on a benchmark: global MAP-state T = 1.24",
    fitted=True,
)

Distribution = list[tuple[str, float]]


def scale(dist: Sequence[tuple[str, float]], t: float) -> Distribution:
    """`dist` with each probability raised to 1/t and renormalised; same order."""
    if t == 1.0:
        return list(dist)
    powered = [(s, p ** (1 / t) if p > 0 else 0.0) for s, p in dist]
    total = sum(p for _, p in powered)
    if total <= 0:
        return list(dist)
    return [(s, p / total) for s, p in powered]


def apply(
    dists: Sequence[Sequence[tuple[str, float]]], temps: Sequence[float]
) -> list[Distribution]:
    """Every column's distribution scaled by its own temperature."""
    if len(dists) != len(temps):
        raise ValueError(f"columns differ: posteriors {len(dists)}, temperatures {len(temps)}")
    return [scale(d, t) for d, t in zip(dists, temps, strict=True)]


def _log10_rel(value: float, reference: float) -> float:
    return math.log10(max(value, FLOOR) / reference)


def column_statistics(rows: Sequence[str]) -> tuple[list[float], list[float]]:
    """(entropy, gap fraction) per alignment column over the extant `rows`."""
    if not rows:
        return [], []
    width = len(rows[0])
    if any(len(r) != width for r in rows):
        raise ValueError("alignment rows differ in length")
    most = math.log2(AMINO_STATES)
    entropy, gaps = [], []
    for column in zip(*rows, strict=True):
        residues = Counter(c.upper() for c in column if c not in NOT_RESIDUES)
        n = sum(residues.values())
        h = -sum(k / n * math.log2(k / n) for k in residues.values()) if n else 0.0
        entropy.append(min(h / most, 1.0))
        gaps.append(sum(1 for c in column if c in "-.") / len(column))
    return entropy, gaps


@dataclass(frozen=True)
class NodePlace:
    branch_length: float  # mean length of the branches meeting at the node
    tip_distance: float  # path length to the nearest extant sequence


def node_place(newick: str, node: str) -> NodePlace | None:
    """Where `node` sits in the tree, rooting-free (the unrooted tree.nwk). None if the
    tree has no node of that name."""
    tree = load_tree(newick)
    internal = (n for n in tree.preorder_node_iter() if n.label == node and not n.is_leaf())
    target = next(internal, None)
    return None if target is None else _place(tree, target)


def tip_place(newick: str, tip: str) -> NodePlace | None:
    """The same for an extant sequence held out and reconstructed (ESR, `phylobot.esr`):
    its one branch, and the path to the nearest *other* sequence, as if its own were
    unknown, which is what it is while it is reconstructed."""
    tree = load_tree(newick)
    target = tree.find_node_with_taxon_label(tip)
    return None if target is None else _place(tree, target)


def _place(tree: dendropy.Tree, target: dendropy.Node) -> NodePlace:
    neighbours: dict[object, list[tuple[object, float]]] = {}
    for n in tree.preorder_node_iter():
        if n.parent_node is not None:
            length = float(n.edge.length or 0.0)
            neighbours.setdefault(n, []).append((n.parent_node, length))
            neighbours.setdefault(n.parent_node, []).append((n, length))
    lengths = [length for _, length in neighbours.get(target, [])]
    # Dijkstra to the nearest leaf: through the parent as well as the children.
    best: dict[int, float] = {id(target): 0.0}
    queue: list[tuple[float, int, object]] = [(0.0, 0, target)]
    tie = 1
    nearest = math.inf
    while queue:
        d, _, n = heapq.heappop(queue)
        if d > best.get(id(n), math.inf):
            continue
        if n.is_leaf() and n is not target:  # type: ignore[attr-defined]
            nearest = d
            break
        for other, length in neighbours.get(n, []):
            if d + length < best.get(id(other), math.inf):
                best[id(other)] = d + length
                heapq.heappush(queue, (d + length, tie, other))
                tie += 1
    mean = sum(lengths) / len(lengths) if lengths else 0.0
    return NodePlace(mean, nearest if nearest < math.inf else 0.0)


def site_rates(tsv: str | None, width: int) -> list[float]:
    """Per column, IQ-TREE's rate from site_rates.tsv; 1 everywhere without one (IQ-TREE
    writes none without rate heterogeneity, where every site's rate is 1)."""
    if not tsv:
        return [1.0] * width
    rates = [1.0] * width
    for line in tsv.splitlines()[1:]:
        fields = line.split("\t")
        if len(fields) >= 2 and fields[0].isdigit() and 1 <= int(fields[0]) <= width:
            rates[int(fields[0]) - 1] = float(fields[1])
    return rates


def covariates(
    extants: Sequence[str], newick: str, node: str, rates_tsv: str | None
) -> list[dict[str, float]] | None:
    """Per alignment column, the node's covariates (each 0 at its reference). None if
    the node is not in the tree."""
    place = node_place(newick, node)
    return None if place is None else covariate_rows(place, extants, rates_tsv)


def covariate_rows(
    place: NodePlace, extants: Sequence[str], rates_tsv: str | None
) -> list[dict[str, float]]:
    """Per alignment column, the covariates of a node at `place`."""
    entropy, gaps = column_statistics(extants)
    rates = site_rates(rates_tsv, len(entropy))
    branch = _log10_rel(place.branch_length, 0.1)
    tips = _log10_rel(place.tip_distance, 0.1)
    return [
        {
            "branch_length": branch,
            "tip_distance": tips,
            "site_rate": _log10_rel(rate, 1.0),
            "entropy": h,
            "gap_fraction": g,
        }
        for h, g, rate in zip(entropy, gaps, rates, strict=True)
    ]


def temperatures(coefficients: Coefficients, rows: Sequence[Mapping[str, float]]) -> list[float]:
    return [coefficients.temperature(x) for x in rows]


@dataclass(frozen=True)
class Summary:
    """T over an ancestor's sites, for the viewer's note."""

    low: float
    median: float
    high: float
    sites: int


def summary(temps: Sequence[float]) -> Summary | None:
    if not temps:
        return None
    ordered = sorted(temps)
    n = len(ordered)
    middle = ordered[n // 2] if n % 2 else (ordered[n // 2 - 1] + ordered[n // 2]) / 2
    return Summary(ordered[0], middle, ordered[-1], n)
