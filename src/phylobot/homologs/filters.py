"""The homolog filter stack, cheapest first, with the count left after each filter.

Filters only remove hits; nothing is reordered. `apply` returns every stage's count so
the review page can show the effect of each filter as the user adjusts it.
"""

from collections.abc import Callable
from dataclasses import dataclass, field

from phylobot.homologs.hits import Hit

# Placeholder bands pending the owner's sign-off. "Recent"
# keeps close homologs, so the ancestors are young and well supported; "ancient" reaches
# down to the twilight zone's edge for deep ancestors.
PRESETS = {
    "recent": {"min_identity": 50.0, "max_identity": 90.0},
    "ancient": {"min_identity": 30.0, "max_identity": 90.0},
}


@dataclass(frozen=True)
class Settings:
    max_evalue: float = 1e-5
    min_coverage: float = 70.0  # % of the query covered by the hit's alignments
    min_identity: float = 30.0
    max_identity: float = 90.0

    def __post_init__(self) -> None:
        if not self.max_evalue >= 0:
            raise ValueError("the E-value cut-off must be 0 or more")
        for name in ("min_coverage", "min_identity", "max_identity"):
            if not 0 <= getattr(self, name) <= 100:
                raise ValueError(f"{name.replace('_', ' ')} must be between 0 and 100")
        if self.min_identity > self.max_identity:
            raise ValueError("minimum identity is above maximum identity")

    @classmethod
    def preset(cls, name: str, **overrides: float) -> Settings:
        return cls(**{**PRESETS[name], **overrides})


@dataclass(frozen=True)
class Stage:
    name: str
    description: str
    remaining: int


def stages(settings: Settings) -> list[tuple[str, str, Callable[[Hit], bool]]]:
    s = settings
    return [
        ("evalue", f"E-value ≤ {s.max_evalue:g}", lambda h: h.evalue <= s.max_evalue),
        (
            "coverage",
            f"covers ≥ {s.min_coverage:g}% of the query",
            lambda h: h.coverage >= s.min_coverage,
        ),
        (
            "identity",
            f"identity {s.min_identity:g}-{s.max_identity:g}%",
            lambda h: s.min_identity <= h.identity <= s.max_identity,
        ),
    ]


def apply(hits: list[Hit], settings: Settings) -> tuple[list[Hit], list[Stage]]:
    """The hits that pass every filter, and the count remaining after each one."""
    kept = list(hits)
    report = [Stage("search", "search results", len(kept))]
    for name, description, keep in stages(settings):
        kept = [hit for hit in kept if keep(hit)]
        report.append(Stage(name, description, len(kept)))
    return kept, report


# The query appears among the hits as a full-length, fully identical match (blastp finds
# the query's own sequence in the database). The identity filter's upper bound is meant to
# drop near-identical *redundant* homologs, not the query itself, so a hit this close is
# treated as the query and kept by default — the reconstruction needs the sequence it
# started from. The user can still untick it.
QUERY_IDENTITY = 99.95
QUERY_COVERAGE = 99.0


def is_query_hit(hit: Hit) -> bool:
    return hit.identity >= QUERY_IDENTITY and hit.coverage >= QUERY_COVERAGE


def query_ids(hits: list[Hit]) -> frozenset[str]:
    return frozenset(h.id for h in hits if is_query_hit(h))


MAX_CHOICES = 5000  # hand-picked includes + excludes honoured in one selection


BASES = ("filters", "all", "none")  # what the selection starts from before hand ticks


@dataclass(frozen=True)
class Selection:
    """What goes into the job or download: the hits the filters keep, plus the ones the
    user included by hand, minus the ones they excluded by hand. "Select all" and "Select
    none" change the base the hand ticks apply to (`base`), so a choice of every hit is
    one URL parameter, never a thousand IDs."""

    by_filters: list[Hit]  # kept by the filters
    selected: list[Hit]  # after the user's own choices, in search order
    included: list[Hit]  # selected although the filters remove them
    excluded: list[Hit]  # not selected although the filters keep them
    stages: list[Stage]
    base: str = "filters"
    default: frozenset[str] = frozenset()  # the base's hits, before hand ticks
    include_ids: list[str] = field(default_factory=list)  # ticks against the base
    exclude_ids: list[str] = field(default_factory=list)
    always: frozenset[str] = frozenset()  # the query hits kept by default (is_query_hit)


def select(
    hits: list[Hit],
    settings: Settings,
    include: list[str],
    exclude: list[str],
    base: str = "filters",
    always: frozenset[str] = frozenset(),
) -> Selection:
    """Apply the filters, then the user's choices: from the filters' hits (or every hit,
    or none, as `base` says), add `include` and remove `exclude`. IDs not among the
    hits are ignored, as are ticks that agree with the base (they change nothing).
    `always` names hits kept by default on top of the filters (the query, `query_ids`);
    they go into the base like filter-kept hits, so a tick removing one is an exclude and
    the user can still drop it. "Select none" starts from nothing, query included."""
    kept, stages = apply(hits, settings)
    kept_ids = {h.id for h in kept}
    if base not in BASES:
        base = "filters"
    always = frozenset(always) & {h.id for h in hits}
    base_default = {"filters": kept_ids, "all": {h.id for h in hits}, "none": set()}[base]
    default = base_default | (always if base != "none" else frozenset())
    include_ids = set(include[:MAX_CHOICES]) - default
    exclude_ids = (set(exclude[:MAX_CHOICES]) & default) - include_ids
    chosen = (default | include_ids) - exclude_ids
    selected = [h for h in hits if h.id in chosen]
    # `included` is hand-picks only (the query, auto-kept, is reported on its own, not as a
    # choice the user made); `excluded` is filter-kept hits the user dropped.
    included = [h for h in selected if h.id not in kept_ids and h.id not in always]
    excluded = [h for h in kept if h.id not in chosen]
    query_kept = any(h.id in always and h.id not in kept_ids for h in selected)
    if included or excluded:
        stages = [*stages, Stage("choices", "after your own choices", len(selected))]
    elif query_kept:
        stages = [*stages, Stage("choices", "with the query kept", len(selected))]
    return Selection(
        kept, selected, included, excluded, stages, base, frozenset(default),
        [h.id for h in hits if h.id in include_ids],
        [h.id for h in hits if h.id in exclude_ids],
        always,
    )  # fmt: skip


def reason(hit: Hit, settings: Settings) -> str:
    """Which filter removes this hit ("" if none), for the not-selected list."""
    for _, description, keep in stages(settings):
        if not keep(hit):
            return f"fails {description}"
    return ""
