"""What the libraries landing page says about each library, computed from stored data.

Pure functions over a library's spec, its job's params and its export receipts: nothing
here reads object storage or generates members. That is the point — reading a library's
posteriors to count its members costs about a third of a second each, and a user may
have thirty (owner decision 2026-10-07). The card shows what the database already knows,
and the member count appears once the library has been designed and exported.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

# A library's state, read from the spec and the export receipts rather than stored: the
# three steps of the flow (collect ancestors -> design the set -> export).
DRAFT = "draft"
DESIGNED = "designed"
EXPORTED = "exported"
STATES = (DRAFT, DESIGNED, EXPORTED)
STATE_LABELS = {
    DRAFT: "Draft",
    DESIGNED: "Designed",
    EXPORTED: "Exported",
}
# The spec keys L5a writes when the design is saved. A spec holding any of them has been
# designed, however the page's controls grow later.
DESIGN_KEYS = ("axes", "budget", "strategy", "filters")


def designed(spec: dict[str, Any]) -> bool:
    """Whether the design half has been saved, as opposed to targets alone."""
    return any(spec.get(key) not in (None, {}, [], "") for key in DESIGN_KEYS)


def state(spec: dict[str, Any], exports: int) -> str:
    """Which of the three steps a library has reached."""
    if exports:
        return EXPORTED
    return DESIGNED if designed(spec) else DRAFT


def state_label(spec: dict[str, Any], exports: int) -> str:
    """The badge's text: "Exported x3" counts them, the others do not."""
    key = state(spec, exports)
    if key == EXPORTED:
        return f"{STATE_LABELS[key]} ×{exports}" if exports > 1 else STATE_LABELS[key]
    return STATE_LABELS[key]


@dataclass(frozen=True)
class Axis:
    """One axis a library can draw on, and how much of it the job actually ran.

    `count` is the number of choices the job has; an axis with fewer than two gives the
    design nothing to compare, which is what the robustness preset exists to fix.
    """

    key: str
    label: str  # the plural, for "3 aligners"
    singular: str
    count: int

    @property
    def thin(self) -> bool:
        return self.count < 2

    @property
    def sentence(self) -> str:
        return f"{self.count} {self.label if self.count != 1 else self.singular}"


def model_count(params: dict[str, Any]) -> int:
    """How many substitution models the job's combinations use: the fixed ones it was
    given, plus one for ModelFinder's pick (one model per alignment, whatever the size of
    the candidate set it chose from)."""
    fixed = list(params.get("models") or []) if params.get("fixed_models") else []
    codon = list(params.get("codon_models") or [])
    picked = 1 if params.get("modelfinder") else 0
    picked += 1 if params.get("codon_modelfinder") else 0
    return len(fixed) + len(codon) + picked


def axes(params: dict[str, Any]) -> list[Axis]:
    """The job's axes, from its stored parameters — never from `summary.json`, which
    would be a storage read per card."""
    aligners = [a for a in params.get("aligners") or [] if a != "none"]
    trees = [t for t in params.get("tree_builders") or [] if t != "none"]
    gaps = list(params.get("gap_methods") or [])
    return [
        Axis("aligners", "aligners", "aligner", len(aligners)),
        Axis("trees", "tree searches", "tree search", len(trees)),
        Axis("models", "models", "model", model_count(params)),
        Axis("gaps", "gap methods", "gap method", len(gaps)),
    ]


def inventory(params: dict[str, Any]) -> str:
    """The one-line axes inventory: "3 aligners · 2 tree searches · 4 models"."""
    return " · ".join(axis.sentence for axis in axes(params) if axis.count)


@dataclass(frozen=True)
class JobEntry:
    """One job a library draws on, as the card's sub-item lists it.

    A library belongs to one job today, so this is a list of one; collections of libraries
    across jobs (planned) are what make it a list, and the axes inventory is per job
    because a library over two jobs has two of them.
    """

    job: Any  # phylobot.models.Job
    inventory: str


def job_entry(job: Any) -> JobEntry:
    """A job as the card lists it. There is deliberately no "this job is thin" notice: the
    inventory beside it already says "1 aligner", and in v1 a library compares one
    combination at a time, so a single aligner costs nothing yet. When L5a's axes controls
    arrive, a disabled axis can say why in the one place it matters."""
    return JobEntry(job=job, inventory=inventory(job.params or {}))


@dataclass(frozen=True)
class Card:
    """One library on the landing page. Everything here is stored, so a page of thirty
    costs one query and no storage reads."""

    library: Any  # phylobot.models.Library
    exports: int
    last_export: Any | None  # the newest LibraryExport, for "Download last export"
    targets: int
    nodes: list[str]  # the ancestors' node names, for a recognisable line of chips
    distinct_nodes: int
    state: str
    state_label: str
    jobs: list[JobEntry]

    @property
    def more_nodes(self) -> int:
        """Ancestors whose nodes the card's chips do not show."""
        return max(self.distinct_nodes - len(self.nodes), 0)

    @property
    def members(self) -> int | None:
        """The last export's member count, or None where nothing has been exported: the
        only honest member count a page that never generates can show."""
        return self.last_export.member_count if self.last_export else None


NODES_SHOWN = 6


def card(library: Any, exports: Sequence[Any]) -> Card:
    """One library's card. `exports` is its receipts, newest first (the model's own
    ordering), passed in so the view can prefetch them in one query."""
    spec = library.spec or {}
    targets = library.targets
    nodes: list[str] = []
    for target in targets:
        node = str(target.get("node", ""))
        if node and node not in nodes:
            nodes.append(node)
    return Card(
        library=library,
        exports=len(exports),
        last_export=exports[0] if exports else None,
        targets=len(targets),
        nodes=nodes[:NODES_SHOWN],
        distinct_nodes=len(nodes),
        state=state(spec, len(exports)),
        state_label=state_label(spec, len(exports)),
        jobs=[job_entry(library.job)],
    )
