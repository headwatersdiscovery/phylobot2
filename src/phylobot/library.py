"""Ancestral library design: candidate sequences from a finished job's reconstructions.

Pure Python over the files the pipeline already wrote — no Django, no object storage, no
compute — so the views, a future API and the tests all call the same functions. It
grows in this order:

1. `alt_all`, the standard robustness control of the ASR literature (Eick et al. 2017,
   MBE 34: 247): every site whose second-most-probable state is plausible is switched to
   it, giving one sequence that is "wrong everywhere it could be wrong".
2. `generate` / `reduce` / `select` / `coverage`, the library proper: every sequence the
   ticked axes imply, reduced to distinct proteins, cut to what the user can order.

Distributions come from `ancestors.distributions` (per alignment column, the states and
their probabilities, most probable first). Thresholds are never constants in the caller:
they are controls on the page, saved in the library's spec, and their defaults live here.
"""

from __future__ import annotations

import random
import re
import string
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from itertools import combinations
from math import comb, inf

# Owner decision 2026-10-04: AltAll's default threshold. The literature uses 0.2 or 0.3;
# the user can change it, so this is only where the page starts.
ALT_ALL_THRESHOLD = 0.2
# What a user may ask for. Above 0.5 an "alternative" would have to beat the ML state,
# which never happens, and the sequence would be the ML one with extra steps.
MIN_THRESHOLD = 0.01
MAX_THRESHOLD = 0.5


def parse_threshold(text: object, default: float = ALT_ALL_THRESHOLD, places: int = 2) -> float:
    """A probability threshold from request text: the default for anything unparsable,
    clamped to [MIN_THRESHOLD, MAX_THRESHOLD] and rounded, so it round-trips through a
    URL and names a file reproducibly."""
    try:
        value = float(str(text))
    except TypeError, ValueError:
        return default
    if value != value or value in (float("inf"), float("-inf")):  # NaN, infinities
        return default
    return round(min(max(value, MIN_THRESHOLD), MAX_THRESHOLD), places)


@dataclass(frozen=True)
class Swap:
    """One site where AltAll departs from the ML reconstruction."""

    site: int  # 1-based alignment column
    ml: str
    alt: str
    pp: float  # the alternative's posterior probability
    ml_pp: float


@dataclass(frozen=True)
class AltAll:
    """The AltAll sequence for one ancestor, and the sites it switched."""

    sequence: str  # aligned: the template's gaps are kept
    threshold: float
    swaps: list[Swap] = field(default_factory=list)

    @property
    def ambiguous(self) -> int:
        return len(self.swaps)

    @property
    def unaligned(self) -> str:
        return self.sequence.replace("-", "")

    def sites_tsv(self) -> str:
        """The switched sites, one per row: the ambiguous sites of this ancestor."""
        rows = ["site\tml\tml_pp\talt\talt_pp"]
        rows += [f"{s.site}\t{s.ml}\t{s.ml_pp:.5f}\t{s.alt}\t{s.pp:.5f}" for s in self.swaps]
        return "\n".join(rows) + "\n"


def split_columns(template: str, width: int = 1) -> list[str]:
    """An aligned row as its columns: one character each for residues, three for codons."""
    if width < 1 or len(template) % width:
        raise ValueError(f"row of {len(template)} is not whole columns of {width}")
    return [template[i : i + width] for i in range(0, len(template), width)]


def alternatives(
    template: str,
    dists: list[list[tuple[str, float]]],
    width: int = 1,
) -> list[Swap]:
    """Every unmasked column's best alternative to the reconstructed state, whatever its
    probability. Independent of any threshold, so one pass answers "how many sites would
    AltAll switch at t?" for every t — which is what the viewer's slider reports.
    """
    columns = split_columns(template, width)
    if len(columns) != len(dists):
        raise ValueError(f"columns differ: ancestor {len(columns)}, posteriors {len(dists)}")
    gap = "-" * width
    found: list[Swap] = []
    for site, (ml, dist) in enumerate(zip(columns, dists, strict=True), start=1):
        if ml == gap:  # masked: no state to have an alternative to
            continue
        others = [(state, pp) for state, pp in dist if state != ml]
        if not others:
            continue
        # Highest probability, ties by state, so the choice never depends on input order.
        alt, pp = min(others, key=lambda sp: (-sp[1], sp[0]))
        ml_pp = next((p for state, p in dist if state == ml), 0.0)
        found.append(Swap(site=site, ml=ml, alt=alt, pp=pp, ml_pp=ml_pp))
    return found


def count_at(found: list[Swap], threshold: float) -> int:
    """How many of `alternatives` AltAll would take at `threshold`."""
    return sum(1 for swap in found if swap.pp >= threshold)


def alt_all(
    template: str,
    dists: list[list[tuple[str, float]]],
    threshold: float = ALT_ALL_THRESHOLD,
    width: int = 1,
    found: list[Swap] | None = None,
) -> AltAll:
    """AltAll for the ancestor `template` (an aligned row: `ancestors.fasta`, one of its
    gap-masked copies, or the codon row at `width=3`) given its per-column distributions.

    At every column the most probable state other than the template's is taken if its
    probability is at least `threshold`. Masked columns (a gap in the template) stay
    gaps: AltAll is an alternative *reconstruction*, not an alternative gap placement,
    which is a separate axis of the library. `found` reuses an `alternatives` call.
    """
    columns = split_columns(template, width)
    if found is None:
        found = alternatives(template, dists, width)
    swaps = [swap for swap in found if swap.pp >= threshold]
    taken = {swap.site: swap.alt for swap in swaps}
    sequence = "".join(taken.get(site, ml) for site, ml in enumerate(columns, start=1))
    return AltAll(sequence=sequence, threshold=threshold, swaps=swaps)


# --- The library proper: generator -> reducer -> selector ----------------------------
#
# The axes explode (four aligners x two tree sources x three models x four gap versions
# x twenty samples is 1,920 sequences for one node) and most of what they give is the
# same sequence, so identity is the signal: generate everything the ticked axes imply,
# reduce to distinct proteins carrying how many choices produced each, then select a set
# that fits what the user can order.
#
# Axes that vary *within* one ancestor (AltAll, posterior samples, single-site
# alternates) live in `Axes`. Axes that *multiply* ancestors (alignment, tree, model, gap
# version) are the caller's, which passes one `Target` for each, so this module never
# learns the file layout and the cross-combination axes need no change of shape.

# One letter per gap version in member names; the original reconstruction adds none.
# The ancestor viewer's version menu also holds PAML codeml's reconstructions, which are
# not gap versions but reconstructions in their own right: see `reconstruction_segment`.
GAP_LETTERS = {
    "": "",
    "none": "",
    "original": "",  # the ancestor viewer's key for the unmasked reconstruction
    "parsimony": "P",
    "fitch": "F",
    "ml": "M",
    "prank": "K",
    "grasp": "G",
}
# codeml's reconstructions are keyed "codeml-<model>" (`ancestors.CODEML_PREFIX`; not
# imported, so this module stays free of the file layout) and named "D<model>" (dN/dS),
# e.g. DM3, DFR: one segment per model, so two models' ancestors never share a name.
CODEML_VERSION = "codeml-"
CODEML_LETTER = "D"
# Member names go on a vendor's order form. Owner decision 2026-10-07; the cap is still
# to be checked against the vendors' forms, so a name over it is flagged, never cut.
NAME_LIMIT = 24
CODE_LIMIT = 6
NAME_CHARACTERS = frozenset(string.ascii_uppercase + string.digits + "-")
# The coverage readout's threshold: a control on the page, this is where it starts.
COVERAGE_THRESHOLD = ALT_ALL_THRESHOLD
# `max_substitutions` above 1 is a binomial in the ambiguous sites (40 sites: 780 pairs,
# 9,880 triples). `count_candidates` is cheap, so the page refuses before generating.
CANDIDATE_LIMIT = 20_000
STRATEGIES = ("supported", "spanning", "per_node")


@dataclass(frozen=True)
class Strategy:
    """How a fill strategy is presented: its name, the one-line hint beside it, and
    whether the page keeps it behind "Advanced".

    `advanced` is an editorial judgement, not a capability: every strategy stays
    selectable, because which of them serves a bench best is an open question (only accuracy
    was measured, on one simulated family). What that settled is that
    `supported` is the safe default, so it is the one the page loads with and the others
    are a disclosure away (owner decision 2026-10-13).
    """

    key: str
    name: str
    hint: str
    advanced: bool = False


STRATEGY_CHOICES = (
    Strategy("supported", "Most supported", "the choices most settings agree on"),
    Strategy(
        "spanning", "Spanning", "cover the disagreement, not near-identical variants", advanced=True
    ),  # fmt: skip
    Strategy("per_node", "Per ancestor", "equal shares, for a lineage", advanced=True),
)


def strategy_choices(chosen: str) -> tuple[tuple[Strategy, ...], tuple[Strategy, ...], bool]:
    """(the strategies shown, those behind "Advanced", whether to open it).

    The disclosure opens when the library's saved strategy is one of the advanced ones,
    so a spec never hides the setting it is actually using.
    """
    plain = tuple(s for s in STRATEGY_CHOICES if not s.advanced)
    advanced = tuple(s for s in STRATEGY_CHOICES if s.advanced)
    return plain, advanced, any(s.key == chosen for s in advanced)


# How a sequence reached the library, best name first: ML > ALT > single-site > sample.
KINDS = ("ml", "altall", "single", "sample")
CONTROLS = ("ml", "altall")
# `Axes.samples` has no default of its own — the generator must be asked — but the page
# has one: 10 (owner decision 2026-10-07). Measured on the real Ime2 library, ten samples
# take the whole coverage gain (81.9% -> 97.7% at budget 24) and a hundred add 180 members
# for nothing; 0 would leave every strategy returning one set, since `spanning` alone ever
# picks a sample. The presets sit beside the field.
DEFAULT_SAMPLES = 10
SAMPLE_PRESETS = (0, 10, 100)
DEFAULT_BUDGET = 24
DEFAULT_STRATEGY = "supported"
MAX_BUDGET = 10_000  # a number of wells, not a limit anyone should reach
MAX_SUBSTITUTIONS = 2  # above 2 the binomial runs away; `count_candidates` refuses first


@dataclass(frozen=True)
class Source:
    """The ancestor a candidate came from: one combination, node and gap version.

    `combination` is the index in `summary.json`'s order (stable for a finished job),
    which `names.tsv` maps back to aligner / tree / model. `weight` is the combination's
    `w_bic` from `model_comparison.tsv`. `label` is for display only and never reaches a
    member name.
    """

    combination: int
    node: str
    gap: str = ""
    weight: float = 0.0
    label: str = ""

    @property
    def key(self) -> tuple[int, str, str]:
        """What makes two sources the same choice; weight and label only ride along."""
        return (self.combination, self.node, self.gap)


@dataclass(frozen=True)
class Target:
    """One ancestor to draw candidates from: its aligned ML row (`ancestors.fasta`, one
    of its gap versions, or the codon row at `width=3`) and its per-column
    distributions, as `ancestors.distributions` gives them."""

    source: Source
    template: str
    dists: list[list[tuple[str, float]]]
    width: int = 1

    @property
    def gap_column(self) -> str:
        return "-" * self.width


@dataclass(frozen=True)
class Axes:
    """The sequence-space depth controls, saved in the library's spec.

    `max_substitutions` m gives every combination of up to m single-site alternates
    (those at or above `threshold`); 0 is off, 1 the single-site alternates themselves.
    `samples` is off at 0 and has no module default (see above).
    """

    alt_all: bool = True
    threshold: float = ALT_ALL_THRESHOLD
    samples: int = 0
    max_substitutions: int = 0


@dataclass(frozen=True)
class Candidate:
    """One generated sequence, before deduplication."""

    sequence: str  # aligned, so the preview can draw it as an alignment row
    variant: str  # ML, ALT, S42, K300R (or K300RL22M for several alternates at once)
    kind: str  # one of KINDS
    source: Source

    @property
    def unaligned(self) -> str:
        return self.sequence.replace("-", "")

    @property
    def priority(self) -> tuple[int, int, str, str, int, str]:
        """Which of several candidates names a sequence they share: by kind, then by
        where it came from. Sample numbers compare as numbers, so S2 precedes S10."""
        number = int(self.variant[1:]) if self.kind == "sample" else 0
        source = self.source
        return (
            KINDS.index(self.kind),
            source.combination,
            source.node,
            source.gap,
            number,
            self.variant,
        )


def swap_label(swap: Swap) -> str:
    """A single-site alternate as the literature writes it: ML state, site, alternative
    (`K300R`), which is also `paths.Change.label`. The site is the 1-based alignment
    column (owner decision 2026-10-07): the gap axis makes members of one library differ
    in length, so only the column is comparable between them."""
    return f"{swap.ml}{swap.site}{swap.alt}"


def with_swaps(target: Target, swaps: Iterable[Swap]) -> str:
    """The target's reconstruction with `swaps` applied."""
    taken = {swap.site: swap.alt for swap in swaps}
    columns = split_columns(target.template, target.width)
    return "".join(taken.get(site, ml) for site, ml in enumerate(columns, start=1))


def posterior_sample(target: Target, rng: random.Random) -> str:
    """One draw from the ancestor's per-column *marginal* posteriors.

    IQ-TREE gives no joint distribution, so every column is drawn on its own — which the
    page and the methods paragraph must say. Masked columns stay gaps: a sample is an
    alternative reconstruction, not an alternative gap placement. Exactly one draw is
    taken per unmasked column, so sample n is the same sequence however many samples were
    asked for, which is what makes a stored seed reproduce a member months later.
    """
    columns = split_columns(target.template, target.width)
    if len(columns) != len(target.dists):
        raise ValueError(f"columns differ: ancestor {len(columns)}, posteriors {len(target.dists)}")
    gap = target.gap_column
    drawn: list[str] = []
    for ml, dist in zip(columns, target.dists, strict=True):
        if ml == gap:
            drawn.append(ml)
            continue
        total = sum(pp for _, pp in dist)
        point = rng.random() * total
        running = 0.0
        state = ml  # a distribution that sums to nothing leaves the ML state standing
        for candidate_state, pp in dist:
            running += pp
            if point < running:
                state = candidate_state
                break
        drawn.append(state)
    return "".join(drawn)


def ambiguous_swaps(target: Target, threshold: float) -> list[Swap]:
    """The target's alternatives at or above `threshold`: its ambiguous sites."""
    found = alternatives(target.template, target.dists, target.width)
    return [swap for swap in found if swap.pp >= threshold]


def count_candidates(targets: Sequence[Target], axes: Axes) -> int:
    """How many sequences the ticked axes imply, without generating one: the live count
    the library page shows before the user commits to a budget (the New Job page's
    honest-arithmetic pattern). Deduplication only ever reduces it."""
    total = 0
    for target in targets:
        ambiguous = len(ambiguous_swaps(target, axes.threshold))
        total += 1 + (1 if axes.alt_all else 0) + max(axes.samples, 0)
        total += sum(comb(ambiguous, k) for k in range(1, max(axes.max_substitutions, 0) + 1))
    return total


def generate(targets: Sequence[Target], axes: Axes, seed: str | int = "") -> list[Candidate]:
    """Every sequence the ticked axes imply, one `Candidate` each (duplicates included:
    `reduce` is what collapses them, and how many choices agree is the signal).

    Deterministic: posterior samples draw from `seed`, which the library's spec keeps, so
    the same spec over a finished job's unchanging files reproduces the same members.
    """
    if (asked := count_candidates(targets, axes)) > CANDIDATE_LIMIT:
        raise ValueError(
            f"{asked} sequences asked for, over the {CANDIDATE_LIMIT} limit: "
            "lower the substitutions from ML, the samples, or the number of targets"
        )
    made: list[Candidate] = []
    for target in targets:
        source = target.source
        made.append(Candidate(target.template, "ML", "ml", source))
        swaps = ambiguous_swaps(target, axes.threshold)
        if axes.alt_all:
            design = alt_all(
                target.template, target.dists, axes.threshold, target.width, found=swaps
            )
            made.append(Candidate(design.sequence, "ALT", "altall", source))
        for size in range(1, max(axes.max_substitutions, 0) + 1):
            for chosen in combinations(swaps, size):
                # Several alternates run together, since the character set has no
                # separator to spare ("K300R/L22M" is how a bench would write it).
                variant = "".join(swap_label(swap) for swap in chosen)
                made.append(Candidate(with_swaps(target, chosen), variant, "single", source))
        if axes.samples > 0:
            # One stream per target, keyed by the source, so a target's samples do not
            # shift when another target is added to the library.
            rng = random.Random(f"{seed}:{source.combination}:{source.node}:{source.gap}")
            for number in range(1, axes.samples + 1):
                drawn = posterior_sample(target, rng)
                made.append(Candidate(drawn, f"S{number}", "sample", source))
    return made


def short_code(code: str) -> str:
    """A library's code as a name segment: upper case, letters and digits, `CODE_LIMIT`
    long. Anything unusable falls back to LIB, so a name is always well formed."""
    kept = "".join(c for c in str(code).upper() if c.isascii() and c.isalnum())
    return kept[:CODE_LIMIT] or "LIB"


def node_tag(node: str) -> str:
    """The node as a name segment: `Node15` -> `N15`, anything else `N` + its
    characters, so a tree whose nodes are named otherwise still gives a usable name."""
    if match := re.fullmatch(r"node(\d+)", str(node), re.IGNORECASE):
        return f"N{match.group(1)}"
    kept = "".join(c for c in str(node).upper() if c.isascii() and c.isalnum())
    return f"N{kept}"


def member_name(code: str, source: Source, variant: str) -> str:
    """`<CODE>-C<k>-N<node>-<variant>[-<gap>]`, e.g. `IME2-C1-N15-ALT-F` (owner decision
    2026-10-07).

    Every segment is always present, never dropped because it does not vary, so names
    stay stable as an editable library grows and match across export receipts; only the
    gap letter is left off the original reconstruction. A very wide alignment can push a
    single-site name past `NAME_LIMIT`, so the page checks `name_fits` and warns rather
    than renaming silently — a name that changed between exports would be worse.
    """
    parts = [short_code(code), f"C{source.combination}", node_tag(source.node), variant.upper()]
    if segment := reconstruction_segment(source.gap):
        parts.append(segment)
    return "-".join(parts)


def reconstruction_segment(version: str) -> str:
    """The name segment for an ancestor viewer version: a gap method's letter, or
    `CODEML_LETTER` + the codeml model ("" for the original reconstruction).

    An unknown version is refused rather than abbreviated: a guessed letter once gave
    every codeml model "C", so ancestors of different models shared one name.
    """
    if version in GAP_LETTERS:
        return GAP_LETTERS[version]
    if version.startswith(CODEML_VERSION):
        model = short_code(version[len(CODEML_VERSION) :])
        if model != "LIB":  # short_code's fallback: no usable model name
            return CODEML_LETTER + model
    raise ValueError(f"no member-name segment for the reconstruction {version!r}")


SWAP_PATTERN = re.compile(r"([A-Z-])(\d+)([A-Z-])")


def parse_swaps(variant: str) -> list[tuple[str, int, str]]:
    """The swaps a single-site variant's name describes: (ML state, alignment column,
    alternative), in the order they were written. "K300R" gives one; "K300RL22M", two.

    The column is the 1-based *alignment* column, as the naming decision requires; an
    export carries the residue number beside it, which only the ancestor's own row knows.
    """
    return [(m[1], int(m[2]), m[3]) for m in SWAP_PATTERN.finditer(variant)]


def residue_number(template: str, column: int, width: int = 1) -> int:
    """Which residue of this ancestor sits at that alignment column, counting from 1.
    0 when the column is a gap there, which an export shows rather than inventing one."""
    columns = split_columns(template, width)
    if not 1 <= column <= len(columns):
        return 0
    gap = "-" * width
    if columns[column - 1] == gap:
        return 0
    return sum(1 for c in columns[:column] if c != gap)


def variant_label(variant: str, kind: str) -> str:
    """The variant as a table reads it. The member *name* keeps the short form (`S3`),
    which is what goes on an order form and must not drift; this is for the eye only."""
    if kind == "ml":
        return "Most likely"
    if kind == "altall":
        return "AltAll"
    if kind == "sample":
        # "Draw", not "sample": a reader takes a sample to be a specimen of the library,
        # while this is one draw from the per-site posteriors (owner, 2026-10-08).
        return f"Draw {variant[1:]}" if variant[1:].isdigit() else variant
    return variant  # a single-site alternate is already its own label, e.g. K300R


def name_fits(name: str) -> bool:
    """Whether a member name is short enough and in the character set a vendor's order
    form accepts."""
    return bool(name) and len(name) <= NAME_LIMIT and set(name) <= NAME_CHARACTERS


@dataclass(frozen=True)
class Member:
    """One distinct protein, with every choice that produced it.

    Members are deduplicated on the *unaligned* sequence — the protein an
    experimentalist would order — while `sequence` keeps the representative candidate's
    aligned row, so the preview can draw the set as alignment-viewer rows.
    """

    sequence: str
    kind: str
    variant: str
    source: Source  # the representative candidate's
    candidates: list[Candidate] = field(default_factory=list)
    width: int = 1

    @property
    def unaligned(self) -> str:
        return self.sequence.replace("-", "")

    @property
    def length(self) -> int:
        """Residues (codons at `width=3`), which is what goes in the order."""
        return len(self.unaligned) // self.width

    @property
    def sources(self) -> list[Source]:
        """The distinct choices that gave this sequence, in a stable order."""
        seen: dict[tuple[int, str, str], Source] = {}
        for candidate in self.candidates:
            seen.setdefault(candidate.source.key, candidate.source)
        return [seen[key] for key in sorted(seen)]

    @property
    def support(self) -> int:
        """How many choices agree on this sequence. Several variants of one source (its
        ML sequence and a sample that landed on it) count once: the question is how many
        reasonable choices gave this protein, not how many paths reached it."""
        return len(self.sources)

    @property
    def weight(self) -> float:
        """The summed `w_bic` of the distinct combinations behind it, so a sequence that
        the well-fitting models agree on outranks one with the same support from poor
        ones."""
        by_combination = {source.combination: source.weight for source in self.sources}
        return sum(by_combination.values())

    @property
    def gaps(self) -> list[str]:
        """The gap versions this sequence came from ("" for the original)."""
        return sorted({source.gap for source in self.sources})

    @property
    def is_control(self) -> bool:
        return self.kind in CONTROLS

    def name(self, code: str) -> str:
        return member_name(code, self.source, self.variant)

    def differences(self, reference: str) -> list[str]:
        """This member against an aligned reference row (the ML ancestor), as change
        labels in the mutational paths' convention. Gaps take part, so a gap version's
        losses read `K300-`."""
        mine = split_columns(self.sequence, self.width)
        theirs = split_columns(reference, self.width)
        if len(mine) != len(theirs):
            raise ValueError(f"rows of {len(theirs)} and {len(mine)} columns cannot be compared")
        return [
            f"{was}{site}{now}"
            for site, (was, now) in enumerate(zip(theirs, mine, strict=True), start=1)
            if was != now
        ]

    def mean_posterior(self, target: Target) -> float:
        """The mean posterior probability of the states this member carries, over its
        unmasked columns: how probable the whole sequence's residues are, which separates
        the ML ancestor (highest by construction) from a sample (lowest)."""
        columns = split_columns(self.sequence, self.width)
        if len(columns) != len(target.dists):
            raise ValueError(
                f"columns differ: member {len(columns)}, posteriors {len(target.dists)}"
            )
        gap = "-" * self.width
        found = [
            next((pp for state, pp in dist if state == mine), 0.0)
            for mine, dist in zip(columns, target.dists, strict=True)
            if mine != gap
        ]
        return sum(found) / len(found) if found else 0.0


def reduce(candidates: Iterable[Candidate], width: int = 1) -> list[Member]:
    """Distinct proteins, each carrying the choices that produced it.

    Identity is on the unaligned sequence, so the same protein reached through different
    alignments or gap versions is one member. Its name comes from its highest-priority
    candidate (ML > ALT > single-site > sample) and the manifest lists them all. The
    order is support, then summed weight, then the representative's own order, so the
    members table and every export are stable for a given spec.
    """
    groups: dict[str, list[Candidate]] = {}
    for candidate in candidates:
        groups.setdefault(candidate.unaligned, []).append(candidate)
    members = [
        Member(
            sequence=(best := min(group, key=lambda c: c.priority)).sequence,
            kind=best.kind,
            variant=best.variant,
            source=best.source,
            candidates=sorted(group, key=lambda c: c.priority),
            width=width,
        )
        for group in groups.values()
    ]
    members.sort(key=lambda m: (-m.support, -m.weight, m.candidates[0].priority))
    return members


def distance(a: str, b: str, width: int = 1) -> int:
    """Columns at which two aligned rows differ: the spanning strategy's distance.

    Rows of different widths cannot be compared, which is why v1 spans within one
    alignment; the cross-alignment axes (alternative alignments, trees and models) are
    what need a length-tolerant distance.
    """
    if len(a) != len(b):
        raise ValueError(f"rows of {len(a)} and {len(b)} columns cannot be compared")
    mine, theirs = split_columns(a, width), split_columns(b, width)
    return sum(1 for x, y in zip(mine, theirs, strict=True) if x != y)


@dataclass(frozen=True)
class Selection:
    """The set to order, and what it was chosen from."""

    members: list[Member]
    budget: int
    strategy: str
    considered: int

    @property
    def dropped(self) -> int:
        return max(self.considered - len(self.members), 0)


# The ancestor viewer's keys for the unmasked reconstruction, which carries every
# insertion seen in any sequence: on an indel dataset the worst of a target's versions.
UNMASKED = ("", "none", "original")


def gap_rank(member: Member, preferred: Sequence[str] = ()) -> int:
    """How early a control's gap version should be ordered: the job's default gap method
    (`preferred`) first, then the other masks and reconstructions, the unmasked original
    last. A member several versions agree on takes its best rank."""

    def rank(gap: str) -> int:
        if gap in preferred:
            return 0
        return 2 if gap in UNMASKED else 1

    return min((rank(gap) for gap in member.gaps), default=2)


def control_key(member: Member, preferred_gaps: Sequence[str] = ()) -> tuple[object, ...]:
    """The order `select` reserves controls in: ML before AltAll, then the best-fitting
    combination (summed `w_bic`), then the default gap version.

    The first wells are what a small library orders, so they must hold the ancestor the
    user would have ordered without one. The old order, `(kind, combination, node, gap)`,
    put `gap=""` — the unmasked reconstruction — first: over the ten `ime2/indels`
    replicates well 1 was -0.188 identity against that ancestor; this order gives -0.004.
    """
    return (
        KINDS.index(member.kind),
        -member.weight,
        gap_rank(member, preferred_gaps),
        member.source.combination,
        member.source.node,
        member.source.gap,
    )


def select(
    members: Sequence[Member],
    budget: int,
    strategy: str = "supported",
    controls: bool = True,
    preferred_gaps: Sequence[str] = (),
) -> Selection:
    """The members to order, at most `budget` of them.

    The controls go in first — the ML ancestor and AltAll, the ASR literature's reference
    and its robustness control (Eick et al. 2017), in `control_key`'s order, with the
    job's default gap method (`preferred_gaps`) ahead of the other versions — then the
    strategy fills the rest:

    - `supported`: the choices that most reasonable settings agree on (`reduce`'s order).
    - `spanning`: greedy max-min, so the wells cover the disagreement instead of holding
      24 near-identical variants.
    - `per_node`: equal shares per node, round robin, for a library over a lineage. The
      share is the *representative's* node, so one protein reconstructed identically at
      two nodes counts for the first of them — it is still one well to order.

    Members are taken in the order given, which is `reduce`'s ranking; `select` never
    re-ranks them. `controls=False` only stops places being reserved: a control that is
    also the best supported member is still picked on its merits.

    The budget is never exceeded, because it is a number of wells: if the controls alone
    overrun it they are cut too, and the page says so.
    """
    if strategy not in STRATEGIES:
        raise ValueError(f"unknown strategy {strategy!r}: expected one of {', '.join(STRATEGIES)}")
    chosen: list[Member] = []
    taken: set[str] = set()

    def take(member: Member) -> None:
        if member.unaligned not in taken:
            taken.add(member.unaligned)
            chosen.append(member)

    if controls:
        for member in sorted(
            (m for m in members if m.is_control),
            key=lambda m: control_key(m, preferred_gaps),
        ):
            take(member)
    rest = [member for member in members if member.unaligned not in taken]
    room = max(budget - len(chosen), 0)

    if strategy == "spanning":
        # Only rows of one alignment can be compared (`distance` refuses the rest), so
        # span within each comparable group and take from the groups in turn, the
        # largest first. A library over two alignments used to raise from here.
        # Spanning *across* alignments needs `phylobot.alignment_map` to put their
        # columns on one frame.
        frames: dict[tuple[int, int], list[Member]] = {}
        for member in rest:
            frames.setdefault((len(member.sequence), member.width), []).append(member)
        queues = sorted(frames.values(), key=lambda q: (-len(q), len(q[0].sequence)))
        # Each queue carries the distance from every member still in it to the nearest
        # thing already chosen from its frame. Recomputing that at every pick is what made
        # a 629-member library take fifteen seconds; one pass against the member just
        # taken is the same greedy choice for a fraction of the work.
        # Nothing chosen yet means "unreached", not zero: a zero would stick, because the
        # update below only ever lowers the distance. The first pick is then decided by
        # the tie-break, which is the incoming ranking, exactly as before.
        nearest: list[list[float]] = []
        for queue in queues:
            width, columns = queue[0].width, len(queue[0].sequence)
            picked = [m.sequence for m in chosen if m.width == width and len(m.sequence) == columns]
            nearest.append(
                [
                    min((distance(m.sequence, seq, width) for seq in picked), default=inf)
                    for m in queue
                ]
            )
        while len(chosen) < budget and any(queues):
            for queue, far in zip(queues, nearest, strict=True):
                if not queue or len(chosen) >= budget:
                    continue
                # The member furthest from everything comparable picked so far. Ties —
                # and the first pick, when no control was reserved — go to the better
                # ranked member, so the incoming order stays the one ranking here.
                furthest = max(range(len(queue)), key=lambda i: (far[i], -i))
                latest = queue.pop(furthest)  # not `taken`: that name is the chosen set
                far.pop(furthest)
                take(latest)
                for i, member in enumerate(queue):
                    step = distance(member.sequence, latest.sequence, member.width)
                    if step < far[i]:
                        far[i] = step
    elif strategy == "per_node":
        by_node: dict[str, list[Member]] = {}
        for member in rest:
            by_node.setdefault(member.source.node, []).append(member)
        queues = [by_node[node] for node in sorted(by_node)]
        while len(chosen) < budget and any(queues):  # round robin, a node at a time
            for queue in queues:
                if queue and len(chosen) < budget:
                    take(queue.pop(0))
    else:
        for member in rest[:room]:
            take(member)

    return Selection(
        members=chosen[:budget], budget=budget, strategy=strategy, considered=len(members)
    )


@dataclass(frozen=True)
class Coverage:
    """What a selected set achieves, stated so it can be checked rather than believed."""

    threshold: float
    sites: int  # ambiguous columns considered
    total: int  # plausible states at them
    covered: int
    missing: list[tuple[str, int, str]] = field(default_factory=list)  # node, column, state

    @property
    def fraction(self) -> float:
        return self.covered / self.total if self.total else 1.0

    @property
    def complete(self) -> bool:
        return self.covered == self.total

    @property
    def sentence(self) -> str:
        """The readout the library page shows above the budget."""
        if not self.total:
            return (
                f"No column has two states at pp >= {self.threshold:.2f}: "
                "nothing is ambiguous at this threshold."
            )
        if self.complete:
            return (
                f"Every state with pp >= {self.threshold:.2f} at all {self.sites} ambiguous "
                f"sites appears in at least one member."
            )
        return (
            f"{self.covered} of {self.total} states with pp >= {self.threshold:.2f} at "
            f"{self.sites} ambiguous sites appear in a member; {self.total - self.covered} "
            "do not."
        )


def plausible_states(
    targets: Sequence[Target], threshold: float
) -> list[tuple[int, int, str, Target]]:
    """Every state a library could be asked to carry: at each ambiguous column of each
    target, one entry per state at or above `threshold`. `coverage` counts how many of
    these a set holds; the curve below counts them as members are added."""
    out: list[tuple[int, int, str, Target]] = []
    for index, target in enumerate(targets):
        columns = split_columns(target.template, target.width)
        gap = target.gap_column
        for site, (ml, dist) in enumerate(zip(columns, target.dists, strict=True)):
            if ml == gap:
                continue
            states = [state for state, pp in dist if pp >= threshold]
            if len(states) < 2:
                continue
            out.extend((index, site, state, target) for state in states)
    return out


def coverage_curve(
    ordered: Sequence[Member], targets: Sequence[Target], threshold: float = COVERAGE_THRESHOLD
) -> list[int]:
    """How many plausible states the first k members carry, for k = 0, 1, ... in the
    order `select` chose them.

    Every strategy's selection at budget k is the first k of its picks — `supported`
    takes a prefix, `per_node` goes round robin, `spanning` adds one at a time — so one
    walk gives the whole curve, and the page can show what another budget would buy
    without running the selector again.
    """
    wanted = plausible_states(targets, threshold)
    if not wanted:
        return [0] * (len(ordered) + 1)
    rows: dict[int, list[str]] = {}
    covered = [False] * len(wanted)
    curve = [0]
    found = 0
    for member in ordered:
        for i, (index, site, state, target) in enumerate(wanted):
            if covered[i]:
                continue
            if index not in rows:
                fits = len(member.sequence) == len(target.template) and member.width == target.width
                rows[index] = split_columns(member.sequence, target.width) if fits else []
            row = rows[index]
            if row and row[site] == state:
                covered[i] = True
                found += 1
        rows.clear()
        curve.append(found)
    return curve


def coverage(
    selected: Sequence[Member],
    targets: Sequence[Target],
    threshold: float = COVERAGE_THRESHOLD,
) -> Coverage:
    """Whether the chosen set carries every plausible state at every ambiguous site.

    A site is ambiguous when two or more states reach `threshold`; a state is covered
    when some member has it in that column. Masked columns are skipped: there is no
    state to cover where the ancestor has no residue. Members whose row does not match
    the target's alignment are left out of its check, which is how a library over
    several alignments (later work) will not claim coverage it cannot have.
    """
    sites = total = covered = 0
    missing: list[tuple[str, int, str]] = []
    for target in targets:
        rows = [
            split_columns(member.sequence, target.width)
            for member in selected
            if len(member.sequence) == len(target.template) and member.width == target.width
        ]
        columns = split_columns(target.template, target.width)
        gap = target.gap_column
        for index, (ml, dist) in enumerate(zip(columns, target.dists, strict=True), start=1):
            if ml == gap:
                continue
            plausible = [state for state, pp in dist if pp >= threshold]
            if len(plausible) < 2:
                continue
            sites += 1
            for state in plausible:
                total += 1
                if any(row[index - 1] == state for row in rows):
                    covered += 1
                else:
                    missing.append((target.source.node, index, state))
    return Coverage(threshold=threshold, sites=sites, total=total, covered=covered, missing=missing)


# --- The saved spec (the `Library` record's `spec` field) ----------------------------
#
# Members are never stored: a finished job's files do not change and the generator is
# deterministic, so the spec reproduces them, and sequences stay out of Postgres. The
# spec is plain JSON. These functions are its only writers, so its shape and defaults
# live here and the views stay thin. Axes, thresholds, budget and strategy join it with
# the library page; a spec without them means their defaults.

SPEC_VERSION = 1


@dataclass(frozen=True)
class Design:
    """What the library page's Design section holds, as saved in the spec: the axes, how
    many wells there are, and how to fill them."""

    axes: Axes
    budget: int = DEFAULT_BUDGET
    strategy: str = DEFAULT_STRATEGY

    @property
    def entry(self) -> dict[str, object]:
        """The design as stored. The threshold rides with the axes although it is a View
        control, so an export reproduces the set it was taken from."""
        return {
            "axes": {
                "alt_all": bool(self.axes.alt_all),
                "threshold": float(self.axes.threshold),
                "samples": int(self.axes.samples),
                "max_substitutions": int(self.axes.max_substitutions),
            },
            "budget": int(self.budget),
            "strategy": self.strategy,
        }


def default_design(threshold: float = ALT_ALL_THRESHOLD) -> Design:
    return Design(axes=Axes(threshold=threshold, samples=DEFAULT_SAMPLES))


def design_from(spec: dict[str, object], threshold: float | None = None) -> Design:
    """The saved design, or the default where a spec has none (every library before L5a).
    `threshold` overrides the stored one: it is the page's View control, so the live
    slider governs what counts as ambiguous."""
    stored = spec.get("axes")
    axes_entry = stored if isinstance(stored, dict) else {}
    saved_threshold = axes_entry.get("threshold")
    use = threshold if threshold is not None else saved_threshold
    axes = Axes(
        alt_all=bool(axes_entry.get("alt_all", True)),
        threshold=float(use if isinstance(use, int | float) else ALT_ALL_THRESHOLD),
        samples=clamp_int(axes_entry.get("samples", DEFAULT_SAMPLES), 0, 10_000),
        max_substitutions=clamp_int(axes_entry.get("max_substitutions", 0), 0, MAX_SUBSTITUTIONS),
    )
    strategy = spec.get("strategy")
    return Design(
        axes=axes,
        budget=clamp_int(spec.get("budget", DEFAULT_BUDGET), 1, MAX_BUDGET),
        strategy=strategy if strategy in STRATEGIES else DEFAULT_STRATEGY,
    )


def clamp_int(value: object, low: int, high: int) -> int:
    """A stored or posted number brought into range; anything unusable gives `low`."""
    try:
        number = int(float(str(value)))
    except TypeError, ValueError:
        return low
    return max(low, min(high, number))


@dataclass(frozen=True)
class View:
    """How the library is shown: the choices that change nothing about which sequences
    exist, only how they are drawn. Saved beside the design so reopening a library brings
    back the view it was left in (owner, 2026-10-08).

    The keys are checked by the page against its own catalogues of palettes, schemes and
    alignments, not here: this module knows nothing of either.
    """

    palette: str = ""
    scheme: str = ""
    zoom: str = ""
    alignment: str = ""
    extant: tuple[str, ...] = ()
    members: str = ""  # which members the members' viewer draws

    @property
    def entry(self) -> dict[str, object]:
        return {
            "palette": self.palette,
            "scheme": self.scheme,
            "zoom": self.zoom,
            "alignment": self.alignment,
            "extant": list(self.extant[:NAMES_SAVED]),
            "members": self.members,
        }


NAMES_SAVED = 200  # extant sequences remembered in a view; a guard, not a real limit


def view_from(spec: dict[str, object]) -> View:
    """The saved view, or an empty one. Every field is read as text and left to the page
    to recognise, so an unknown palette or aligner falls back there rather than here."""
    stored = spec.get("view")
    entry = stored if isinstance(stored, dict) else {}
    names = entry.get("extant")
    return View(
        palette=str(entry.get("palette", "") or ""),
        scheme=str(entry.get("scheme", "") or ""),
        zoom=str(entry.get("zoom", "") or ""),
        alignment=str(entry.get("alignment", "") or ""),
        extant=tuple(str(n) for n in names[:NAMES_SAVED]) if isinstance(names, list) else (),
        members=str(entry.get("members", "") or ""),
    )


def saved_at(spec: dict[str, object]) -> str:
    """When the design was last saved, as stored (ISO 8601), or "" if it never was."""
    when = spec.get("saved_at")
    return str(when) if isinstance(when, str) else ""


def with_design(
    spec: dict[str, object], design: Design, view: View | None = None, when: str = ""
) -> dict[str, object]:
    """The spec with this design — and the view it was saved in — written into it. Targets
    and seed are untouched."""
    out = dict(spec) | design.entry
    if view is not None:
        out["view"] = view.entry
    if when:
        out["saved_at"] = when
    return out


def new_spec(seed: int) -> dict[str, object]:
    """An empty library's spec. The seed is drawn once, when the library is created, so
    its posterior samples are the same at every visit and in every export."""
    return {"targets": [], "seed": int(seed)}


def target_entry(
    combination: str,
    node: str,
    version: str,
    partition: Iterable[Iterable[str]],
    below: int,
) -> dict[str, object]:
    """One target as stored. `version` is the ancestor viewer's key ("original", a gap
    method, or a codeml reconstruction).

    `partition` is the subtrees (as leaf sets) that meet at the node, from
    `gaps.node_partitions`: a node's identity independent of rooting. v1 matches targets
    by node name only, but storing it lets the cross-tree axes find the same ancestor
    under another topology without migrating saved specs. A single split would not do:
    it names a *branch*, and so fits the nodes at both of its ends equally. `below` is
    for display only — the sequences below the node on the tree as the viewers draw it.
    """
    return {
        "combination": combination,
        "node": node,
        "version": version,
        "partition": sorted(sorted(part) for part in partition),
        "below": int(below),
    }


def target_key(entry: dict[str, object]) -> tuple[str, str, str]:
    """What makes two targets the same: an ancestor in one reconstruction."""
    return (str(entry["combination"]), str(entry["node"]), str(entry["version"]))


def spec_targets(spec: dict[str, object]) -> list[dict[str, object]]:
    targets = spec.get("targets")
    return [t for t in targets if isinstance(t, dict)] if isinstance(targets, list) else []


def add_targets(
    spec: dict[str, object], entries: Iterable[dict[str, object]]
) -> tuple[dict[str, object], int]:
    """The spec with `entries` appended, skipping any already in it (so a second click,
    or a lineage overlapping a node already added, changes nothing), and how many were
    new. Returns a new dict; the caller saves it."""
    targets = spec_targets(spec)
    present = {target_key(t) for t in targets}
    added = 0
    for entry in entries:
        if (key := target_key(entry)) not in present:
            present.add(key)
            targets.append(entry)
            added += 1
    return {**spec, "targets": targets}, added


def remove_target(spec: dict[str, object], key: tuple[str, str, str]) -> dict[str, object]:
    """The spec without the target `key` (by key, never by position, so a page that is
    out of date cannot remove the wrong one)."""
    return {**spec, "targets": [t for t in spec_targets(spec) if target_key(t) != key]}


def default_code(job_name: str, number: int) -> str:
    """A new library's short code: from the job's name, else LIB<n>."""
    code = "".join(c for c in job_name.upper() if c.isascii() and c.isalnum())[:CODE_LIMIT]
    return code or f"LIB{number}"[:CODE_LIMIT]
