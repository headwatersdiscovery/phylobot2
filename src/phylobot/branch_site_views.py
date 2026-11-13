"""The branch-site test as a follow-up job (phylobot.branch_site).

In the tree viewer of a coding-DNA job the user picks a foreground (a node, and either
the branch leading to it or its whole clade); the tree highlights it, and "Test for
positive selection here" queues a new Job: its own accounting record, on that
combination's tree and the job's codon alignment for that aligner. Its inputs are
written here from the parent's results (alignment1.fasta, tree1.nwk with codeml's #1
marks), and its params say what it tests (`branch_site_of`); the pipeline's
`--branch_site` mode does the rest. The job page shows the result (`results`).
"""

import json
import logging
from dataclasses import dataclass, field
from typing import Any
from urllib.parse import urlencode

from django.db import transaction
from django.http import Http404, HttpRequest, HttpResponse
from django.shortcuts import get_object_or_404, redirect
from django.urls import reverse
from django.views.decorators.http import require_POST

from phylobot import ancestors, branch_site, jobs, manifest, storage, trees, workload
from phylobot.models import Job

logger = logging.getLogger(__name__)

SCOPES = list(branch_site.SCOPES)


def codon_alignment(job: Job, aligner: str, keys: set[str]) -> str | None:
    """The job's codon alignment for `aligner` (published for coding DNA), if any."""
    path = f"{aligner}/alignment.cds.fasta"
    return path if job.params.get("seqtype") == "codon" and path in keys else None


@dataclass
class NodeOption:
    name: str
    label: str
    extant: bool


@dataclass
class Panel:
    """The tree viewer's branch-site panel: which node and scope are chosen, the
    branches they mean on the drawn tree, and the tests already run here."""

    available: bool
    nodes: list[NodeOption] = field(default_factory=list)
    node: str = ""
    scope: str = "branch"
    foreground: set[str] = field(default_factory=set)  # drawn nodes whose branch is marked
    error: str = ""
    tests: list[Job] = field(default_factory=list)  # on the same tree, any combination
    existing: Job | None = None  # this very test, already run or running
    query: dict[str, str] = field(default_factory=dict)  # for the viewer's URL


def node_options(newick: str) -> list[NodeOption]:
    """Ancestral nodes root first, then the extant sequences."""
    tree = trees.load_tree(newick)
    internal = [
        NodeOption(n.label, f"{n.label} · {len(n.leaf_nodes())} sequences", False)
        for n in tree.preorder_internal_node_iter()
        if n.label and n is not tree.seed_node
    ]
    extant = sorted(leaf.taxon.label for leaf in tree.leaf_node_iter())
    return internal + [NodeOption(name, name, True) for name in extant]


def drawn_foreground(newick: str, splits: set[branch_site.Split]) -> set[str]:
    """The names of the drawn tree's nodes whose branch is in `splits` (both of the
    root's branches when the outgroup root splits a foreground branch)."""
    tree = trees.load_tree(newick)
    everyone = branch_site.leaves_below(tree.seed_node)
    return {
        branch_site.name(node)
        for node in tree.preorder_node_iter()
        if node is not tree.seed_node
        and trees.canonical_split(branch_site.leaves_below(node), everyone) in splits
    }


def tests_on(owner: Any, tree: str) -> list[Job]:
    """The user's branch-site tests on this tree (`branch_site.tree_key`: from any of
    their jobs and combinations with the same codon alignment and topology), newest first."""
    return list(
        Job.objects.for_user(owner)
        .filter(params__branch_site_of__tree_key=tree)
        .order_by("-created_at")
    )


def existing_test(owner: Any, key: str) -> Job | None:
    """The user's test of exactly this analysis (`branch_site.test_key`), if one has
    succeeded or is queued or running: the oldest, so every link leads to the same job.
    Failed and cancelled tests do not count, so they can be run again."""
    return (
        Job.objects.for_user(owner)
        .filter(params__branch_site_of__key=key)
        .exclude(status__in=[Job.Status.FAILED, Job.Status.CANCELLED])
        .order_by("created_at")
        .first()
    )


def panel(job: Job, combination: dict[str, Any], keys: set[str], newick: str, query: Any) -> Panel:
    """The panel for the tree viewer: `newick` is the tree as drawn."""
    alignment = codon_alignment(job, combination["aligner"], keys)
    if alignment is None:
        return Panel(available=False)
    codons = storage.read_bytes(f"{job.storage_prefix}results/{alignment}")
    tree = branch_site.tree_key(codons, newick, job.params.get("genetic_code") or 1)
    options = node_options(newick)
    out = Panel(available=True, nodes=options, tests=tests_on(job.owner, tree))
    node = query.get("fg", "")
    if node not in {o.name for o in options}:
        return out
    scope = query.get("fg_scope", "branch")
    out.scope = scope if scope in branch_site.SCOPES else "branch"
    out.node = node
    out.query = {"fg": node, "fg_scope": out.scope}
    try:
        splits = branch_site.foreground_splits(newick, node, out.scope)
    except branch_site.MarkError as exc:
        out.error = str(exc).capitalize() + "."
        return out
    out.foreground = drawn_foreground(newick, splits)
    out.existing = existing_test(job.owner, branch_site.test_key(tree, splits))
    return out


@require_POST
def submit(request: HttpRequest, job_id: str) -> HttpResponse:
    """Queue a branch-site test on one combination of the user's job. Every value is
    checked against the job's own results; nothing from the request reaches a file
    except as a choice among them."""
    parent = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if parent.status != Job.Status.SUCCEEDED:
        raise Http404
    prefix = f"{parent.storage_prefix}results/"
    keys = manifest.paths(parent)
    if "summary.json" not in keys:
        raise Http404
    summary = json.loads(storage.read_bytes(prefix + "summary.json"))
    combos = {c["key"]: c for c in summary.get("combinations", [])}
    combination = combos.get(request.POST.get("combination", ""))
    if combination is None or f"{combination['key']}/tree.nwk" not in keys:
        raise Http404
    alignment = codon_alignment(parent, combination["aligner"], keys)
    if alignment is None:
        raise Http404
    key = combination["key"]
    rooted = f"{key}/tree.rooted.nwk" in keys
    drawn = storage.read_bytes(prefix + f"{key}/{'tree.rooted.nwk' if rooted else 'tree.nwk'}")
    unrooted = storage.read_bytes(prefix + f"{key}/tree.nwk").decode()
    node, scope = request.POST.get("fg", ""), request.POST.get("fg_scope", "")
    if node not in {o.name for o in node_options(drawn.decode())} or scope not in SCOPES:
        raise Http404
    try:
        splits = branch_site.foreground_splits(drawn.decode(), node, scope)
        marked = branch_site.mark_tree(unrooted, splits)
    except branch_site.MarkError:
        raise Http404 from None
    codons = storage.read_bytes(prefix + alignment)
    code = parent.params.get("genetic_code") or 1
    tree = branch_site.tree_key(codons, drawn.decode(), code)
    test = branch_site.test_key(tree, splits)
    of = {
        "job": str(parent.id),
        "combination": key,
        "aligner": combination["aligner"],
        "label": combination.get("label", key),
        "node": node,
        "scope": scope,
        "branches": len(splits),
        "tree_key": tree,
        "key": test,
    }
    params = {
        "branch_site": True,
        "branch_site_of": of,
        "seqtype": "codon",
        "genetic_code": code,
        "aligners": ["none"],
        "tree_builders": ["none"],
        "alignment_labels": [f"{combination['aligner']} codon alignment"],
        "tree_labels": [f"{key} tree"],
    }
    size = workload.size_of(codons.replace(b"-", b""), codons=True) or workload.IME2_SIZE
    name = f"Branch-site test: {node} ({scope}) in {parent.name or str(parent.id)[:8]}"
    with transaction.atomic():
        # One submission at a time per parent job, so two clicks queue one test.
        Job.objects.select_for_update().filter(id=parent.id).first()
        if (existing := existing_test(parent.owner, test)) is not None:
            return redirect("job_detail", job_id=existing.id)  # the same analysis: its result
        job = Job.objects.create(
            owner=parent.owner,
            name=name[:100],
            params=params,
            work_units=workload.estimate(params, size).units,
            input_size=workload.size_record(size),
        )
        job.tags.set(parent.tags.all())
        storage.put_bytes(f"{job.storage_prefix}inputs/alignment1.fasta", codons)
        storage.put_bytes(f"{job.storage_prefix}inputs/tree1.nwk", (marked + "\n").encode())
        transaction.on_commit(lambda: jobs.enqueue(job))
    return redirect("job_detail", job_id=job.id)


@dataclass(frozen=True)
class Site:
    """One codon column: Model A's BEB posterior of each site class (0, 1, 2a, 2b) and
    of positive selection on the foreground (2a + 2b)."""

    column: int
    residue: str  # in the alignment's first sequence ("-" for a gap)
    classes: tuple[float, float, float, float]
    positive: float


SIGNIFICANT = 0.95  # the posterior codeml (and most papers) star
# The ancestor viewer's trace, with room above for the labels of sites near 1.
TRACE_FRAME = ancestors.Frame.of((1000, 280), (66, 14, 32, 58))
LABEL_GAP = 34  # SVG units between two labelled sites' columns; closer ones go unlabelled
HIST_FRAME = ancestors.Frame.of(ancestors.HIST_SIZE, ancestors.HIST_MARGIN)


@dataclass(frozen=True)
class SitePlots:
    """The job page's figures: P(positive selection on the foreground) along the
    alignment, and its distribution over sites."""

    sites: int
    points: str  # SVG polyline
    marks: list[tuple[float, float, int, float, bool]]  # x, y, column, P, labelled; P > 0.95
    threshold_y: float
    x_ticks: list[tuple[float, int]]
    trace_frame: ancestors.Frame
    histogram: list[ancestors.Bar]
    hist_frame: ancestors.Frame
    hist_ticks: list[tuple[float, str]]
    tooltip: list[list[object]]  # per site: column, column, then label/P pairs (ancestors.js)


def site_plots(sites: list[Site]) -> SitePlots | None:
    if not sites:
        return None
    f, n = TRACE_FRAME, len(sites)

    def x(i: int) -> float:
        return f.left + (i / (n - 1) if n > 1 else 0.5) * f.plot_width

    points = " ".join(f"{x(i):.1f},{f.y(s.positive):.1f}" for i, s in enumerate(sites))
    marks = []
    last: float = -LABEL_GAP  # the previous label's x: labels never overlap
    for i, site in enumerate(sites):
        if site.positive > SIGNIFICANT:
            labelled = x(i) - last >= LABEL_GAP
            last = x(i) if labelled else last
            marks.append((round(x(i), 1), round(f.y(site.positive), 1), site.column,
                          site.positive, labelled))  # fmt: skip
    step = ancestors.nice_step(n)
    ticks = [(round(x(k - 1), 1), sites[k - 1].column) for k in [1, *range(step, n + 1, step)]]
    h = HIST_FRAME
    counts = [0] * 10
    for site in sites:
        counts[min(int(site.positive * 10), 9)] += 1
    slot = h.plot_width / 10
    bars = []
    for i, count in enumerate(counts):
        share = count / n
        height = share * h.plot_height
        bars.append(ancestors.Bar(
            x=round(h.left + i * slot + 1, 1), y=round(h.bottom - height, 1),
            width=round(slot - 2, 1), height=round(height, 1), proportion=share,
            label=f"{i / 10:.1f}–{(i + 1) / 10:.1f}: {count} site{'s' * (count != 1)}",
            fill="#2457a6",
        ))  # fmt: skip
    hist_ticks = [(round(h.left + k * 2 * slot, 1), f"{k * 0.2:.1f}") for k in range(6)]
    tooltip: list[list[object]] = [
        [s.column, s.column, f"{s.residue} · positive", round(s.positive, 3),
         *(v for label, p in zip(("class 2a", "class 2b"), s.classes[2:], strict=True)
           for v in (label, round(p, 3)))]
        for s in sites
    ]  # fmt: skip
    return SitePlots(
        n, points, marks, round(f.y(SIGNIFICANT), 1), ticks, f, bars, h, hist_ticks, tooltip
    )


@dataclass
class Result:
    """A finished branch-site job's results, for the job page."""

    statistic: float
    df: int
    p_value: float
    p_mixture: float
    models: list[dict[str, str]]
    sites: list[tuple[int, str, float]]  # BEB: site, residue, P (all codeml lists, > 0.5)
    omega2: float
    omega2_at_bound: bool  # codeml's upper limit (999): the foreground ω2 is not estimable
    every_site: list[Site]  # sites.tsv (jobs before it existed: none)
    plots: SitePlots | None

    @property
    def significant(self) -> list[tuple[int, str, float]]:
        """Sites above 0.95, the posterior conventionally reported."""
        return [s for s in self.sites if s[2] > SIGNIFICANT]

    @property
    def p_text(self) -> str:
        return p_text(self.p_value)

    @property
    def p_mixture_text(self) -> str:
        return p_text(self.p_mixture)

    @property
    def highest(self) -> Site | None:
        return max(self.every_site, key=lambda s: s.positive, default=None)


def p_text(p: float) -> str:
    """A p-value for the page: four decimals, or powers of ten below 0.0001."""
    if p >= 1e-4:
        return f"{p:.4f}".rstrip("0").rstrip(".") if p < 1 else "1"
    if p <= 0:
        return "< 10⁻³⁰⁰"
    mantissa, exponent = f"{p:.1e}".split("e")
    superscript = str(int(exponent)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{mantissa} × 10{superscript}"


def read_tsv(text: str) -> list[dict[str, str]]:
    lines = text.strip().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def results(read: Any, keys: set[str]) -> Result | None:
    """The job's branch_site/ files, if it finished."""
    needed = [f"branch_site/{name}" for name in ("test.tsv", "models.tsv", "beb.tsv")]
    if not all(path in keys for path in needed):
        return None
    [test] = read_tsv(read("branch_site/test.tsv"))
    models = read_tsv(read("branch_site/models.tsv"))
    sites = [(int(r["site"]), r["residue"], float(r["p_positive"]))
             for r in read_tsv(read("branch_site/beb.tsv"))]  # fmt: skip
    omega2 = float(models[0]["omega2"])
    every_site = []
    if "branch_site/sites.tsv" in keys:
        every_site = [
            Site(int(r["site"]), r["residue"],
                 (float(r["p0"]), float(r["p1"]), float(r["p2a"]), float(r["p2b"])),
                 float(r["p_positive"]))
            for r in read_tsv(read("branch_site/sites.tsv"))
        ]  # fmt: skip
    return Result(
        statistic=float(test["statistic"]),
        df=int(test["df"]),
        p_value=float(test["p_value"]),
        p_mixture=float(test["p_mixture"]),
        models=models,
        sites=sites,
        omega2=omega2,
        omega2_at_bound=omega2 >= 998,
        every_site=every_site,
        plots=site_plots(every_site),
    )


def context(job: Job, keys: list[str], read: Any) -> dict[str, Any]:
    """The job page's branch-site part: what was tested, links back, and the result."""
    of = job.params.get("branch_site_of") or {}
    parent_id = str(of.get("job") or "")  # written by `submit`, so a valid UUID if set
    parent = Job.objects.for_user(job.owner).filter(id=parent_id).first() if parent_id else None
    tree_url = ""
    if parent is not None:
        query = {"combination": of.get("combination", ""), "fg": of.get("node", ""),
                 "fg_scope": of.get("scope", "branch")}  # fmt: skip
        tree_url = f"{reverse('tree_viewer', args=[parent.id])}?{urlencode(query)}"
    return {
        "branch_site_of": of,
        "branch_site_foreground": branch_site.describe(of),
        "branch_site_parent": parent,
        "branch_site_tree_url": tree_url,
        "branch_site": results(read, set(keys)) if job.status == Job.Status.SUCCEEDED else None,
    }
