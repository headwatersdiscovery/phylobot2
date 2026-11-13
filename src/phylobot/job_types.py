"""What kind of job each Job is, and the My Jobs filters.

A job's types come from its params, never stored separately. A New Job analysis is
cumulative (alignments, then trees, then ancestors; `stop_after`), so it carries every
stage it reached, plus "dN/dS" when codeml ran. A branch-site test (a follow-up from
the tree viewer) is "dN/dS" too: both are codeml selection analyses on a finished tree,
and the job's name says which test it is. A simulation is one type, and so is an ESR job
(extant-sequence reconstruction, the ESR page). Each type has a
filter on the database (`Q`), so the filters apply before the page's row limit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date
from typing import Any
from urllib.parse import urlencode

from django.db.models import Case, F, IntegerField, OuterRef, Q, QuerySet, Subquery, Value, When
from django.db.models.functions import Lower, NullIf
from django.http import QueryDict

from phylobot.models import Job, Tag


@dataclass(frozen=True)
class JobType:
    key: str
    label: str
    description: str
    group: str  # colour: the work estimate's groups (align, trees, ancestors, codeml) or sim

    @property
    def stage(self) -> bool:
        """One of the cumulative analysis stages (drawn as a chain)."""
        return self.key in STAGES


TYPES = (
    JobType("alignments", "Alignment", "Aligns the sequences (or uses uploaded ones)", "align"),
    JobType("trees", "Trees", "Builds phylogenies and fits substitution models", "trees"),
    JobType("ancestors", "Ancestors", "Reconstructs ancestral sequences", "ancestors"),
    JobType("dnds", "dN/dS", "Tests for selection with codeml (dN/dS) on coding DNA: "
            "models run in a job, or a branch-site test from the tree viewer", "codeml"),
    JobType("simulation", "Simulation", "Simulates sequences along a tree (Simulations page)",
            "sim"),
    JobType("properties", "Ancestral properties", "Measures sequence properties (charge, "
            "hydrophobicity and the rest) over an ancestor's posterior, on ancestors a "
            "finished job already reconstructed", "ancestors"),
    JobType("uncertainty", "Uncertainty", "Splits the variation of a finished job's "
            "ancestors into the posterior's own and that of each choice: aligner, tree, "
            "model, gap placement", "ancestors"),
    JobType("esr", "ESR", "Extant-sequence reconstruction: holds out sequences, reconstructs "
            "each as if it were an ancestor and scores the posterior against the real one "
            "(ESR page)", "ancestors"),
    JobType("reconciled", "Reconciled", "Reconciled with a species tree: GeneRax marks gene "
            "duplications and losses on its trees", "recon"),
)  # fmt: skip
BY_KEY = {t.key: t for t in TYPES}
STAGES = ("alignments", "trees", "ancestors")


def types_of(params: dict[str, Any]) -> list[JobType]:
    """The job's types, in pipeline order."""
    if params.get("simulate"):
        return [BY_KEY["simulation"]]
    if params.get("esr"):
        return [BY_KEY["esr"]]
    if params.get("branch_site"):
        return [BY_KEY["dnds"]]
    if params.get("properties_only"):  # a follow-up from the ancestor viewer (B5c)
        return [BY_KEY["properties"]]
    if params.get("uncertainty_only"):  # a follow-up from the Uncertainty page
        return [BY_KEY["uncertainty"]]
    stop = params.get("stop_after") or "ancestors"
    reached = STAGES[: STAGES.index(stop) + 1] if stop in STAGES else STAGES
    keys = [*reached, *(["dnds"] if params.get("codeml_models") else [])]
    keys += ["reconciled"] if params.get("reconcile") else []
    return [BY_KEY[k] for k in keys]


def _isnt(key: str, value: Any) -> Q:
    """params[key] is missing or not `value` (a plain negation would drop missing keys,
    since comparing a missing JSON key gives NULL)."""
    return ~Q(params__has_key=key) | ~Q(**{f"params__{key}": value})


def type_filter(key: str) -> Q:
    """Jobs that include this type (the database's view of `types_of`)."""
    stop_unset = ~Q(params__has_key="stop_after")
    analysis = (
        _isnt("simulate", True)
        & _isnt("branch_site", True)
        & _isnt("properties_only", True)
        & _isnt("uncertainty_only", True)
        & _isnt("esr", True)
    )
    match key:
        case "simulation":
            return Q(params__simulate=True)
        case "alignments":
            return analysis
        case "trees":
            return analysis & (stop_unset | Q(params__stop_after__in=["trees", "ancestors"]))
        case "ancestors":
            return analysis & (stop_unset | Q(params__stop_after="ancestors"))
        case "dnds":
            codeml = Q(params__has_key="codeml_models") & ~Q(params__codeml_models=[])
            return analysis & codeml | Q(params__branch_site=True) & _isnt("simulate", True)
        case "properties":
            return Q(params__properties_only=True)
        case "uncertainty":
            return Q(params__uncertainty_only=True)
        case "esr":
            return Q(params__esr=True)
        case "reconciled":
            return analysis & Q(params__reconcile=True)
    raise ValueError(key)


@dataclass(frozen=True)
class Composition:
    """What a job was made of, for My Jobs' second line under its name: how many
    sequences, alignments (aligners and the user's own) and models, and tree sources."""

    sequences: int | None
    alignments: list[str]
    models: list[str]
    trees: list[str]

    @property
    def text(self) -> str:
        return " · ".join(self.parts)

    @property
    def parts(self) -> list[str]:
        """The line's phrases (each kept on one line when it wraps)."""
        parts = [f"{self.sequences} sequences"] if self.sequences else []
        if self.alignments:
            n = len(self.alignments)
            aligned = f"{n} alignment{'s' if n != 1 else ''}"
            if self.models:
                m = len(self.models)
                aligned += f" × {m} model{'s' if m != 1 else ''}"
            parts.append(aligned)
        if len(self.trees) > 1:
            parts.append(f"{len(self.trees)} tree sources")
        return parts

    @property
    def details(self) -> str:
        lines = []
        if self.alignments:
            lines.append("Alignments: " + ", ".join(self.alignments))
        if self.models:
            lines.append("Models: " + ", ".join(self.models))
        if self.trees:
            lines.append("Trees: " + ", ".join(self.trees))
        return "\n".join(lines)


ALIGNER_LABELS = {
    "mafft": "MAFFT", "muscle": "MUSCLE", "clustalo": "Clustal Omega", "prank": "PRANK",
    "prank_codon": "PRANK (codon model)",
}  # fmt: skip
TREE_LABELS = {"iqtree": "IQ-TREE", "raxml": "RAxML-NG", "generax": "GeneRax (reconciled)"}


def composition(params: dict[str, Any]) -> Composition | None:
    """A job's composition from its params; None for simulations, ESR and branch-site tests,
    which are described by their type alone."""
    if params.get("simulate") or params.get("branch_site") or params.get("esr"):
        return None
    alignments = [ALIGNER_LABELS.get(a, a) for a in params.get("aligners") or [] if a != "none"]
    alignments += [f"your {label}" for label in params.get("alignment_labels") or []]
    models: list[str] = []
    trees: list[str] = []
    if (params.get("stop_after") or "ancestors") != "alignments":
        if params.get("modelfinder"):
            models.append("ModelFinder")
        if params.get("fixed_models"):
            models += list(params.get("models") or [])
        if params.get("codon_modelfinder"):
            models.append("codon ModelFinder")
        models += list(params.get("codon_models") or [])
        trees = [TREE_LABELS.get(t, t) for t in params.get("tree_builders") or ["iqtree"]
                 if t != "none"]  # fmt: skip
        trees += [f"your {label}" for label in params.get("tree_labels") or []]
    count = params.get("sequence_count")
    return Composition(count if isinstance(count, int) else None, alignments, models, trees)


DATA = {"protein": "Protein", "codon": "Coding DNA"}
ORIGIN = {"real": "Real sequences", "simulated": "Simulated data"}


def _date(value: str | None) -> date | None:
    try:
        return date.fromisoformat(value) if value else None
    except ValueError:
        return None


@dataclass
class Filters:
    """The My Jobs filters, read from the query string (unknown values are dropped)."""

    q: str = ""
    types: list[str] = field(default_factory=list)
    statuses: list[str] = field(default_factory=list)
    data: str = ""
    origin: str = ""
    tag: str = ""
    submitted_from: date | None = None
    submitted_to: date | None = None
    finished_from: date | None = None
    finished_to: date | None = None

    @classmethod
    def parse(cls, query: QueryDict) -> Filters:
        return cls(
            q=query.get("q", "").strip()[:100],
            types=[t for t in query.getlist("type") if t in BY_KEY],
            statuses=[s for s in query.getlist("status") if s in Job.Status.values],
            data=d if (d := query.get("data", "")) in DATA else "",
            origin=o if (o := query.get("origin", "")) in ORIGIN else "",
            tag=query.get("tag", "").lstrip("#").lower(),
            submitted_from=_date(query.get("submitted_from")),
            submitted_to=_date(query.get("submitted_to")),
            finished_from=_date(query.get("finished_from")),
            finished_to=_date(query.get("finished_to")),
        )

    @property
    def active(self) -> bool:
        return any(getattr(self, name) for name in self.__dataclass_fields__)

    @property
    def count(self) -> int:
        """How many filters are set (a date window counts once), for the Filters button."""
        return sum(map(bool, (
            self.q, self.types, self.statuses, self.data, self.origin, self.tag,
            self.submitted_from or self.submitted_to, self.finished_from or self.finished_to,
        )))  # fmt: skip

    def apply(self, rows: QuerySet[Job]) -> QuerySet[Job]:
        if self.q:
            rows = rows.filter(search(self.q))
        if self.types:
            wanted = Q()
            for key in self.types:
                wanted |= type_filter(key)
            rows = rows.filter(wanted)
        if self.statuses:
            rows = rows.filter(status__in=self.statuses)
        if self.data == "codon":
            rows = rows.filter(params__seqtype="codon")
        elif self.data == "protein":
            rows = rows.filter(_isnt("seqtype", "codon"))
        if self.origin == "simulated":
            rows = rows.filter(Q(params__simulated=True) | Q(params__simulate=True))
        elif self.origin == "real":
            rows = rows.filter(_isnt("simulated", True) & _isnt("simulate", True))
        if self.tag:
            rows = rows.filter(tags__name=self.tag)
        for when, start, end in (
            ("created_at", self.submitted_from, self.submitted_to),
            ("finished_at", self.finished_from, self.finished_to),
        ):
            if start:
                rows = rows.filter(**{f"{when}__date__gte": start})
            if end:
                rows = rows.filter(**{f"{when}__date__lte": end})
        return rows


def search(text: str) -> Q:
    """Jobs whose name contains `text`, or whose ID starts with it (as pages show IDs)."""
    match = Q(name__icontains=text)
    if _hex(text):
        match |= Q(id__startswith=text.lower())
    return match


def _hex(text: str) -> bool:
    """Whether the search text could be the start of a job ID (as My Jobs shows them)."""
    return all(c in "0123456789abcdef-" for c in text.lower())


# My Jobs sorting: column key -> (label, first click descending?). Dates start newest
# first, CPU time and storage largest first, the rest A to Z; jobs without a value (no
# name, tag, finish or measurement) always go last.
COLUMNS = {
    "name": ("Job", False),
    "type": ("Type", False),
    "tags": ("Tags", False),
    "submitted": ("Submitted", True),
    "status": ("Status", False),
    "finished": ("Finished", True),
    "cpu": ("CPU", True),  # these two share one header cell ("Usage")
    "disk": ("Storage", True),
}
DEFAULT_SORT = ("submitted", True)
STATUS_ORDER = ("queued", "running", "succeeded", "failed", "cancelled")  # lifecycle


def _sort_value(key: str) -> Any:
    match key:
        case "name":
            return NullIf(Lower("name"), Value(""))
        case "type":  # as the pills read: alignment-only first, simulations last
            return Case(
                When(type_filter("simulation"), then=6),
                When(type_filter("esr"), then=5),
                When(type_filter("dnds"), then=4),
                When(type_filter("ancestors"), then=3),
                When(type_filter("trees"), then=2),
                default=1,
                output_field=IntegerField(),
            )
        case "tags":  # by the job's first tag, alphabetically
            first = Tag.objects.filter(jobs=OuterRef("pk")).order_by("name").values("name")[:1]
            return Subquery(first)
        case "submitted":
            return F("created_at")
        case "status":
            return Case(
                *(When(status=s, then=i) for i, s in enumerate(STATUS_ORDER)),
                output_field=IntegerField(),
            )
        case "finished":
            return F("finished_at")
        case "cpu":
            return F("cpu_seconds")
        case "disk":
            return F("storage_bytes")
    raise ValueError(key)


@dataclass(frozen=True)
class Sort:
    """The My Jobs table's order (`sort` and `dir` in the query string)."""

    key: str = DEFAULT_SORT[0]
    descending: bool = DEFAULT_SORT[1]

    @classmethod
    def parse(cls, query: QueryDict) -> Sort:
        key = query.get("sort", "")
        if key not in COLUMNS:
            return cls()
        direction = query.get("dir", "")
        return cls(key, direction == "desc" if direction in ("asc", "desc") else COLUMNS[key][1])

    @property
    def default(self) -> bool:
        return (self.key, self.descending) == DEFAULT_SORT

    @property
    def params(self) -> dict[str, str]:
        """The query string's sort fields ({} for the default order)."""
        if self.default:
            return {}
        return {"sort": self.key, "dir": "desc" if self.descending else "asc"}

    @property
    def query(self) -> str:
        """ "?sort=...&dir=..." to keep this order on a link, or "" for the default."""
        return "?" + urlencode(self.params) if self.params else ""

    def apply(self, rows: QuerySet[Job]) -> QuerySet[Job]:
        value = _sort_value(self.key)
        order = value.desc(nulls_last=True) if self.descending else value.asc(nulls_last=True)
        return rows.order_by(order, "-created_at", "id")

    def headers(self, query: QueryDict) -> list[dict[str, Any]]:
        """One entry per column: its label, the link that sorts by it (or reverses the
        current order), and aria-sort for the current column; other parameters kept."""
        out = []
        for key, (label, first_desc) in COLUMNS.items():
            current = key == self.key
            descending = not self.descending if current else first_desc
            link = query.copy()
            for name in ("sort", "dir"):
                link.pop(name, None)
            link.update(Sort(key, descending).params)
            state = ("descending" if self.descending else "ascending") if current else ""
            out.append({"key": key, "label": label, "query": link.urlencode(), "sort": state})
        return out
