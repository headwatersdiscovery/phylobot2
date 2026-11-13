"""Property posteriors: a sequence property's distribution over an ancestor's posterior.

The primitive is

    node posterior -> N draws -> map(property function) -> distribution + interval

so the deliverable is never one number: it is a credible interval, the shift between the
property of the ML ancestor and the posterior mean, and (across a job's combinations) a
decomposition separating reconstruction noise from the choice of aligner, tree and model.

This module is pure: no storage, no Django, no containers. It reuses
`library.posterior_sample` and `library.Target` as the sampler and the per-ancestor data
shape, since both are general (masked columns stay gaps, codons work at `width=3`, and a
stored seed reproduces a draw months later). **That reuse is not the Library concept.** A
Library is a budget-constrained, diversity-maximising, control-reserving set of sequences
to order into wells, so it is deliberately a *biased* sample of the posterior and a
credible interval computed from one would be wrong. Properties draw their own unbiased
sample.

The first three property functions are deliberately cheap, dependency-free and chosen to
span sensitivity classes rather than to be the best available predictors. Real
predictors are registered here too, but computed by `phylobot.properties_ml` in the
phylobot-ml image; importing this module never imports torch.
"""

import math
import random
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

from phylobot.library import Source, Target, posterior_sample

GAP = "-"
DEFAULT_DRAWS = 1000  # owner decision 2026-10-20; user-configurable within DRAWS
#: What a user may ask for. The floor is where an interval starts meaning something (a
#: 95% interval from 100 draws is read off its 3rd and 98th values); the ceiling keeps a
#: tick box from queueing an hour of work. Cost is linear in the count either way.
MIN_DRAWS = 100
MAX_DRAWS = 10000


def draws_in_range(asked: object, default: int = DEFAULT_DRAWS) -> int:
    """A draw count from raw request text: the number asked for, clamped, or the default
    for anything that is not a number. The form validates and reports; this is for the
    paths where a silent fallback beats a 404 (the viewer's own controls)."""
    try:
        wanted = int(str(asked))
    except TypeError, ValueError:
        return default
    return max(MIN_DRAWS, min(MAX_DRAWS, wanted))


DEFAULT_LEVEL = 0.95

# Eisenberg et al. 1984 consensus hydrophobicity, used for the hydrophobic moment.
EISENBERG = {
    "A": 0.62,
    "R": -2.53,
    "N": -0.78,
    "D": -0.90,
    "C": 0.29,
    "Q": -0.85,
    "E": -0.74,
    "G": 0.48,
    "H": -0.40,
    "I": 1.38,
    "L": 1.06,
    "K": -1.50,
    "M": 0.64,
    "F": 1.19,
    "P": 0.12,
    "S": -0.18,
    "T": -0.05,
    "W": 0.81,
    "Y": 0.26,
    "V": 1.08,
}
# Kyte & Doolittle 1982 hydropathy, used for the thresholded window count.
KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}
# Net charge at pH 7: histidine's pKa is near 6, so it carries little charge there.
CHARGES = {"K": 1.0, "R": 1.0, "H": 0.1, "D": -1.0, "E": -1.0}

HELIX_DEGREES = 100.0  # residues per turn in an alpha helix, for the moment's phase
MOMENT_WINDOW = 11
HYDROPATHY_WINDOW = 7
HYDROPATHY_THRESHOLD = 1.5


def residues(sequence: str) -> str:
    """The sequence without gaps, upper-cased. Properties are of the protein, not of its
    alignment, so an ancestor's gap columns never contribute."""
    return "".join(c for c in sequence.upper() if c != GAP and c.isalpha())


def net_charge(sequence: str) -> float:
    """Net charge at pH 7. Smooth, additive and *extensive*: it scales with length, so it
    is sensitive to the indel errors gap placement already measures, and blind to swaps
    among like-charged residues."""
    return sum(CHARGES.get(c, 0.0) for c in residues(sequence))


_phase_cache: dict[int, tuple[tuple[float, ...], tuple[float, ...]]] = {}


def _phases(window: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Sines and cosines of each offset's helical phase, cached: the moment is computed
    over hundreds of windows per sequence and millions of draws per study, and recomputing
    these trigonometric values inside the loop dominated the first version."""
    if window not in _phase_cache:
        delta = math.radians(HELIX_DEGREES)
        _phase_cache[window] = (
            tuple(math.sin(delta * k) for k in range(window)),
            tuple(math.cos(delta * k) for k in range(window)),
        )
    return _phase_cache[window]


def hydrophobic_moment(sequence: str, window: int = MOMENT_WINDOW) -> float:
    """Mean Eisenberg hydrophobic moment over sliding windows.

    Smooth, positional and *intensive*: it measures how well hydrophobicity follows the
    helical period, so a substitution matters according to where it falls rather than only
    what it is. This is the class that would notice an A30P-style helix disruption. A
    uniform sequence scores near zero however hydrophobic it is, because the per-residue
    vectors spread around the helix and cancel: this measures asymmetry.
    """
    seq = residues(sequence)
    if len(seq) < window:
        return 0.0
    values = [EISENBERG.get(c, 0.0) for c in seq]
    sines, cosines = _phases(window)
    spans = len(values) - window + 1
    total = 0.0
    for start in range(spans):
        sin_sum = cos_sum = 0.0
        for offset in range(window):
            h = values[start + offset]
            sin_sum += h * sines[offset]
            cos_sum += h * cosines[offset]
        total += math.hypot(sin_sum, cos_sum)
    return total / spans / window


def hydropathy_regions(
    sequence: str,
    window: int = HYDROPATHY_WINDOW,
    threshold: float = HYDROPATHY_THRESHOLD,
) -> float:
    """The number of distinct runs of windows whose mean Kyte-Doolittle hydropathy exceeds
    `threshold`.

    Deliberately **non-smooth**: one substitution can tip a window across the threshold and
    split or merge a run, so the distribution over draws may be multimodal. It stands in
    for a real aggregation predictor, chosen for that *shape* rather than for its biology,
    and it is what tells us whether the machinery survives discontinuous properties.
    """
    seq = residues(sequence)
    if len(seq) < window:
        return 0.0
    values = [KYTE_DOOLITTLE.get(c, 0.0) for c in seq]
    regions = 0
    inside = False
    for start in range(len(values) - window + 1):
        mean = sum(values[start : start + window]) / window
        if mean > threshold:
            if not inside:
                regions += 1
                inside = True
        else:
            inside = False
    return float(regions)


#: Tier 0 is pure Python over the draws (milliseconds per ancestor); the higher tiers
#: need a model, a structure or a simulation and are gated on compute credits. The tier
#: orders the checklist a user picks from when measuring properties after a job has
#: finished, so it is a property of the registry rather than of the page.
TIER_0 = 0


#: Which image computes a property. "tools" is pure Python (phylobot-tools, and the web
#: tier can compute it too, as the ancestor viewer does for f(AltAll)); "ml" needs torch
#: and model weights, which only phylobot-ml carries (`docker/ml.Dockerfile`), so the
#: pipeline runs PROPERTIES there whenever one is asked for.
TOOLS = "tools"
ML = "ml"


@dataclass(frozen=True)
class Property:
    key: str  # short and stable: it names output paths and columns
    label: str
    fn: Callable[[str], float]
    units: str
    note: str
    tier: int = TIER_0
    #: Reference CPU-seconds for 1000 draws of a 400-residue ancestor, from the figures
    #: behind `workload`. Only used to order and price the checklist, so an
    #: order-of-magnitude value is enough; measure before any of it reaches a bill.
    cost: float = 1.0
    image: str = TOOLS
    #: Measured by default: on New Job's "Measure properties" and in the pipeline's own
    #: `properties` default. The neural predictors are opt-in (owner decision),
    #: since they cost minutes to hours where the first three cost seconds; the ancestor
    #: viewer still offers every property, ticking only these.
    default: bool = True
    #: Many sequences in one call, for the predictors whose cost is a model call. None
    #: means `fn` per sequence.
    batch: Callable[[Sequence[str]], list[float]] | None = None
    #: What computed it, written into each summary so a figure names its predictor.
    method: str = ""

    def values(self, sequences: Sequence[str]) -> list[float]:
        """The property of every sequence, computing each distinct one once.

        Draws from a confident ancestor repeat, and a model call is the expensive part,
        so duplicates are folded before the call rather than after."""
        distinct = list(dict.fromkeys(sequences))
        if self.batch is not None:
            computed = dict(zip(distinct, self.batch(distinct), strict=True))
        else:
            computed = {s: self.fn(s) for s in distinct}
        return [computed[s] for s in sequences]


def _disorder_batch(sequences: Sequence[str]) -> list[float]:
    from phylobot import properties_ml  # torch: only in phylobot-ml, only when asked

    return properties_ml.disorder_fraction(sequences)


def _esm_batch(sequences: Sequence[str]) -> list[float]:
    from phylobot import properties_ml

    return properties_ml.esm_log_likelihood(sequences)


def disorder_fraction(sequence: str) -> float:
    """The share of residues in IDRs called by metapredict V3 (`phylobot.properties_ml`)."""
    return _disorder_batch([sequence])[0]


def esm_log_likelihood(sequence: str) -> float:
    """ESM-2 35M's mean per-residue log-probability, one unmasked pass
    (`phylobot.properties_ml`)."""
    return _esm_batch([sequence])[0]


REGISTRY: dict[str, Property] = {
    p.key: p
    for p in (
        Property(
            "charge",
            "Net charge (pH 7)",
            net_charge,
            "e",
            "Additive and extensive; sensitive to indel errors.",
            cost=0.6,
        ),
        Property(
            "moment",
            "Hydrophobic moment",
            hydrophobic_moment,
            "",
            "Positional and intensive; sensitive to helical periodicity.",
            cost=2.4,  # sliding windows over every draw: the dearest of the three
        ),
        Property(
            "hydropathy_regions",
            "Hydrophobic regions",
            hydropathy_regions,
            "count",
            "Thresholded and non-smooth; may be multimodal over draws.",
            cost=1.2,
        ),
        Property(
            "disorder",
            "Disordered fraction",
            disorder_fraction,
            "fraction",
            "Share of residues in predicted disordered regions (metapredict).",
            # Measured 2026-10-21, amd64 emulated on the M2: 36 s for 1000 draws of a
            # 1050-residue ancestor, so ~14 s at the 400-residue reference.
            cost=14.0,
            image=ML,
            default=False,
            batch=_disorder_batch,
            method="metapredict 3.0.2 (V3), IDR calls",
        ),
        Property(
            "esm_likelihood",
            "ESM-2 log-likelihood",
            esm_log_likelihood,
            "nats/residue",
            "How natural the sequence looks to a protein language model; higher is more plausible.",
            # ~1.1 s per 400-residue draw under the M2's amd64 emulation, ~0.25 s native
            # (2026-10-21): the reference figures are the laptop's, so the emulated one.
            cost=1100.0,
            image=ML,
            default=False,
            batch=_esm_batch,
            method="ESM-2 35M (esm2_t12_35M_UR50D), one unmasked pass, mean log p",
        ),
    )
}

#: What "Measure properties" measures unless asked otherwise.
DEFAULTS: tuple[str, ...] = tuple(k for k, p in REGISTRY.items() if p.default)


def image_for(keys: Sequence[str]) -> str:
    """The image a run measuring `keys` needs: phylobot-ml if any of them needs it."""
    return ML if any(REGISTRY[k].image == ML for k in keys if k in REGISTRY) else TOOLS


def stream(seed: str | int, source: Source) -> random.Random:
    """One random stream per seed and ancestor, as `library.generate` does, so a target's
    draws neither depend on nor disturb another's and draw *n* is stable however many were
    asked for."""
    return random.Random(f"{seed}:{source.combination}:{source.node}:{source.gap}")


def draw_sequences(target: Target, draws: int, seed: str | int = "") -> list[str]:
    """`draws` independent draws from the ancestor's per-column marginal posteriors.

    IQ-TREE reports no joint distribution, so columns are drawn independently. That is a
    real limitation, not a shortcut: independent draws can combine individually plausible
    residues into a sequence that never occurs, which is why draws are reweighted by
    plausibility once that exists. Until then an unweighted
    distribution must be reported as such.
    """
    if draws < 1:
        raise ValueError(f"draws must be at least 1, got {draws}")
    rng = stream(seed, target.source)
    return [posterior_sample(target, rng) for _ in range(draws)]


def quantile(sorted_values: Sequence[float], q: float) -> float:
    """The q-quantile by linear interpolation between order statistics."""
    if not sorted_values:
        raise ValueError("no values")
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = q * (len(sorted_values) - 1)
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return sorted_values[low]
    return sorted_values[low] + (sorted_values[high] - sorted_values[low]) * (position - low)


@dataclass(frozen=True)
class Summary:
    """A property's posterior distribution for one ancestor.

    `width` is reported with `covers` everywhere, because coverage alone can be bought by
    being vague: an interval of (-inf, +inf) covers everything and says nothing.
    """

    n: int
    mean: float
    sd: float
    low: float
    high: float
    ml: float
    level: float
    weighted: bool

    @property
    def width(self) -> float:
        return self.high - self.low

    @property
    def shift(self) -> float:
        """f(ML) - posterior mean. Non-zero means the published single ancestor's property
        is not the posterior's centre."""
        return self.ml - self.mean

    @property
    def standardised_shift(self) -> float:
        """The shift in units of the posterior sd, so it compares across properties."""
        return self.shift / self.sd if self.sd > 0 else 0.0

    @property
    def cv(self) -> float:
        """Coefficient of variation.

        Only meaningful for a property that stays well away from zero. Net charge crosses
        zero, so its cv explodes: 13.1 on a real Ime2 ancestor whose posterior mean was
        -0.28. Nothing is therefore *excluded* on cv; see `varies`.
        """
        return self.sd / abs(self.mean) if self.mean else 0.0

    @property
    def varies(self) -> bool:
        """Whether there is any uncertainty to propagate.

        A node whose property is identical in every draw says nothing about calibration,
        and keeping it would flatter coverage for a trivial reason, so the smoke test
        excludes it.
        """
        return self.sd > 0.0

    def covers(self, truth: float) -> bool:
        return self.low <= truth <= self.high


def summarise(
    values: Sequence[float],
    ml: float,
    level: float = DEFAULT_LEVEL,
    weights: Sequence[float] | None = None,
) -> Summary:
    """Summarise draws of a property. `ml` is the property of the ML ancestor, kept apart
    from the draws so the point-estimate shift is always available."""
    if not values:
        raise ValueError("no values to summarise")
    if not 0.0 < level < 1.0:
        raise ValueError(f"level must be between 0 and 1, got {level}")
    if weights is not None and len(weights) != len(values):
        raise ValueError("weights and values differ in length")
    tail = (1.0 - level) / 2.0
    if weights is None:
        mean = sum(values) / len(values)
        spread = sum((v - mean) ** 2 for v in values)
        variance = spread / (len(values) - 1) if len(values) > 1 else 0.0
        ordered = sorted(values)
        low, high = quantile(ordered, tail), quantile(ordered, 1.0 - tail)
    else:
        total = sum(weights)
        if total <= 0:
            raise ValueError("weights sum to zero")
        mean = sum(v * w for v, w in zip(values, weights, strict=True)) / total
        variance = sum(w * (v - mean) ** 2 for v, w in zip(values, weights, strict=True)) / total
        low, high = _weighted_interval(values, weights, total, tail)
    return Summary(
        n=len(values),
        mean=mean,
        sd=math.sqrt(variance),
        low=low,
        high=high,
        ml=ml,
        level=level,
        weighted=weights is not None,
    )


def _weighted_interval(
    values: Sequence[float], weights: Sequence[float], total: float, tail: float
) -> tuple[float, float]:
    """The weighted quantiles bounding the central interval."""
    pairs = sorted(zip(values, weights, strict=True))
    low = high = pairs[0][0]
    running = 0.0
    seen_low = False
    for value, weight in pairs:
        running += weight / total
        if not seen_low and running >= tail:
            low, seen_low = value, True
        if running >= 1.0 - tail:
            high = value
            break
    return low, high


def summarise_target(
    target: Target,
    prop: Property,
    draws: int = DEFAULT_DRAWS,
    seed: str | int = "",
    level: float = DEFAULT_LEVEL,
) -> Summary:
    """Draw from one ancestor and summarise `prop` over the draws."""
    values = [prop.fn(seq) for seq in draw_sequences(target, draws, seed)]
    return summarise(values, prop.fn(target.template), level)


def draw_values(
    target: Target,
    props: Sequence[Property],
    draws: int = DEFAULT_DRAWS,
    seed: str | int = "",
) -> dict[str, list[float]]:
    """Every property's values over one set of draws, keyed by property.

    Drawing is the expensive half, so several properties share one set of draws rather
    than redrawing per property. They therefore see the *same* sequences, which is also
    what makes properties of one ancestor comparable to each other.
    """
    sequences = draw_sequences(target, draws, seed)
    return {prop.key: prop.values(sequences) for prop in props}


def pit(values: Sequence[float], truth: float) -> float:
    """The probability integral transform of `truth` under the draws' distribution.

    The share of draws at or below the truth, with ties split. If the posterior is
    calibrated, PIT values across many ancestors are **uniform on [0, 1]**; a hump in the
    middle means intervals that are too wide, and mass at both ends means too narrow.

    Unlike "did the interval cover it", this uses *every* ancestor rather than only those
    where two arms disagree, which is what gives a calibration study its power.
    """
    if not values:
        raise ValueError("no values")
    below = sum(1 for v in values if v < truth)
    equal = sum(1 for v in values if v == truth)
    return (below + 0.5 * equal) / len(values)


def crps(values: Sequence[float], truth: float) -> float:
    """Continuous ranked probability score of the draws against the truth, lower better.

    A proper scoring rule that rewards being both *right* and *sharp*, so it cannot be
    gamed by widening the interval the way coverage can. Computed from the ensemble as
    mean|X - y| - 0.5 mean|X - X'|, the second term in O(n log n) from the sorted draws.
    """
    if not values:
        raise ValueError("no values")
    n = len(values)
    if n == 1:
        return abs(values[0] - truth)
    ordered = sorted(values)
    accuracy = sum(abs(v - truth) for v in ordered) / n
    # sum_{i<j} (x_j - x_i) = sum_i x_i * (2i - n + 1), i zero-based on the sorted values.
    weighted = sum(v * (2 * i - n + 1) for i, v in enumerate(ordered))
    spread = weighted / (n * n)
    return accuracy - spread


@dataclass(frozen=True)
class Delta:
    """A property's change between two ancestors, with the uncertainty of both.

    The draws of two nodes are independent, so the variance of the difference is the sum
    of theirs; there is no pairing to exploit and pretending otherwise would understate
    the interval.
    """

    mean: float
    sd: float
    low: float
    high: float
    ml: float

    @property
    def excludes_zero(self) -> bool:
        return self.low > 0.0 or self.high < 0.0

    @property
    def sign(self) -> int:
        return (self.mean > 0) - (self.mean < 0)


def delta(older: Summary, newer: Summary, level: float = DEFAULT_LEVEL) -> Delta:
    """`newer` minus `older`: the property's change along the path between two nodes."""
    mean = newer.mean - older.mean
    sd = math.hypot(newer.sd, older.sd)
    # The normal quantile for the central interval; adequate for a difference of means and
    # it needs no dependency. 1.959964 is the 97.5th percentile.
    z = _normal_quantile((1.0 + level) / 2.0)
    return Delta(mean=mean, sd=sd, low=mean - z * sd, high=mean + z * sd, ml=newer.ml - older.ml)


def _normal_quantile(p: float) -> float:
    """The standard normal quantile, by Acklam's rational approximation (absolute error
    below 1e-9 over the range we use). Pure Python keeps `phylobot-tools` dependency-free."""
    if not 0.0 < p < 1.0:
        raise ValueError(f"p must be between 0 and 1, got {p}")
    a = (
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    )
    b = (
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    )
    c = (
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464968e00,
        2.938163982698783e00,
    )
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00, 3.754408661907416e00)
    low, high = 0.02425, 1.0 - 0.02425
    if p < low:
        q = math.sqrt(-2.0 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )
    if p > high:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )
    q = p - 0.5
    r = q * q
    return (
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5])
        * q
        / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    )


@dataclass(frozen=True)
class Components:
    """A property's variance split between and within groups of draws.

    Groups are a job's combinations, or one factor of them (aligner, tree source, model),
    so `between` is how much the property depends on that methodological choice and
    `within` is reconstruction noise given it.

    This is a method-of-moments decomposition, which is honest but not REML: on an
    unbalanced design (UFBoot only on IQ-TREE's own trees, JTTDCMut unavailable on
    RAxML-NG) the components are approximate and `between` can come out slightly negative
    when the true between-group variance is near zero, in which case it is reported as
    zero and `negative` is set. REML would need another dependency.
    """

    groups: int
    n: int
    grand_mean: float
    between: float
    within: float
    negative: bool

    @property
    def total(self) -> float:
        return self.between + self.within

    @property
    def between_fraction(self) -> float:
        return self.between / self.total if self.total > 0 else 0.0


def components(groups: Mapping[str, Sequence[float]]) -> Components:
    """Split a property's variance between and within groups of draws."""
    usable = {key: values for key, values in groups.items() if len(values) > 1}
    if len(usable) < 2:
        raise ValueError("need at least two groups with more than one draw each")
    counts = {key: len(values) for key, values in usable.items()}
    means = {key: sum(values) / counts[key] for key, values in usable.items()}
    n = sum(counts.values())
    grand = sum(sum(values) for values in usable.values()) / n
    k = len(usable)
    within_ss = sum(sum((v - means[key]) ** 2 for v in values) for key, values in usable.items())
    within = within_ss / (n - k)
    between_ss = sum(counts[key] * (means[key] - grand) ** 2 for key in usable)
    between_ms = between_ss / (k - 1)
    # The expected between-group mean square carries a within-group term; n0 is the
    # effective group size for an unbalanced design.
    n0 = (n - sum(c * c for c in counts.values()) / n) / (k - 1)
    raw = (between_ms - within) / n0 if n0 > 0 else 0.0
    return Components(
        groups=k,
        n=n,
        grand_mean=grand,
        between=max(raw, 0.0),
        within=within,
        negative=raw < 0.0,
    )


def sign_consistency(deltas: Sequence[Delta]) -> tuple[int, int, float]:
    """How often a property's change points the same way: (agreeing, total, fraction).

    The headline statistic for whether an evolutionary claim survives method choice
    Zero-mean deltas count in the total but agree with neither direction, so a property
    that does not move cannot look robust.
    """
    signs = [d.sign for d in deltas]
    total = len(signs)
    if not total:
        return 0, 0, 0.0
    agreeing = max(signs.count(1), signs.count(-1))
    return agreeing, total, agreeing / total
