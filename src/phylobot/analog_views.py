"""The same ancestor across the job's combinations, in the viewers (`phylobot.analogs`).

The tree viewer shows each node's analogs in a popover; the ancestor viewer lists them
with their identity to the ancestor shown, draws them on one frame with the library's
alignment viewer (`library_glyphs.frames`), and adds them all to a library at once.
Matched on the job's `ancestor_map.tsv` (written by the pipeline; an extension's joined
with its parent's), one read instead of every combination's tree: parsing 285 trees
took 2.7 s per page.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Any
from urllib.parse import urlencode

from django.http import Http404, HttpRequest, HttpResponse
from django.urls import reverse
from django.views.decorators.http import require_POST

from phylobot import (
    analogs,
    ancestors,
    catalog,
    library,
    library_glyphs,
    library_views,
    results,
)

GAP = analogs.GAP


def by_clade(selection: Any) -> bool:
    """Match by clade on the outgroup-rooted trees, or (no outgroup) by partition."""
    return bool(selection.job.params.get("outgroup"))


def combinations(selection: Any) -> list[dict[str, Any]]:
    """The job's combinations with ancestors, in the viewers' menu order."""
    return sorted(
        (c for c in selection.summary.get("combinations", [])
         if f"{c['key']}/ancestors.fasta" in selection.keys),
        key=results.combination_order(selection.summary.get("primary", "")),
    )  # fmt: skip


def trees(selection: Any) -> dict[str, dict[str, analogs.Located]]:
    """{combination: its nodes' IDs}, for every combination with ancestors, from the
    job's ancestor map. Matching by clade needs the rooted tree; a combination without
    one is left out rather than matched on a root that means nothing. Cached on the
    selection's Results (one request)."""
    files = selection.results  # an extension's map is joined with its parent's
    rows = files.memo.get(analogs.MAP)
    if rows is None:
        rows = analogs.read_map(files.read(analogs.MAP)) if analogs.MAP in files.keys else []
        files.memo[analogs.MAP] = rows
    placed = analogs.placed(rows, by_clade(selection))
    return {c["key"]: placed[c["key"]] for c in combinations(selection) if c["key"] in placed}


def short_label(combo: dict[str, Any], several_trees: bool) -> str:
    """Tree (when the job has several) and model: what tells a combination's ancestor
    from the others on the same alignment."""
    fitted = combo.get("fit", {}).get("model", "")
    model = (
        f"ModelFinder → {fitted}"
        if catalog.is_modelfinder(combo["model"]) and fitted
        else combo.get("model_label", combo["model"])
    )
    return f"{combo.get('tree_label', combo['tree'])} · {model}" if several_trees else model


def row_name(combo: dict[str, Any], several_trees: bool) -> str:
    """`short_label` for the viewer's row names, which are cut at about 40 characters:
    ModelFinder's pick alone, without saying it was ModelFinder's."""
    model = combo.get("fit", {}).get("model", "") if catalog.is_modelfinder(combo["model"]) else ""
    model = model or combo.get("model_label", combo["model"])
    return f"{combo.get('tree_label', combo['tree'])} · {model}" if several_trees else model


def compact(combo: dict[str, Any], several_trees: bool) -> str:
    """`label` with ModelFinder's pick alone (`row_name`), for the tree's popovers."""
    return f"{combo.get('aligner_label', combo['aligner'])} · {row_name(combo, several_trees)}"


def label(combo: dict[str, Any], several_trees: bool) -> str:
    """The combination in full: alignment, then `short_label`."""
    return f"{combo.get('aligner_label', combo['aligner'])} · {short_label(combo, several_trees)}"


def match_note(clade: bool) -> str:
    if clade:
        return (
            "Matched by clade: the node with the same sequences below it on each "
            "combination's outgroup-rooted tree."
        )
    return (
        "No outgroup, so the trees are unrooted: matched where the same subtrees meet, "
        "which needs the same branching on every side of the node."
    )


# Branch supports a popover lists, as supports.tsv names them. Agreement (k/n) is left
# out: it is counted against the tree on screen, so means nothing for another tree's node.
POPOVER_SUPPORTS = (("sh_alrt", "SH-aLRT"), ("ufboot", "UFBoot"))


@dataclass(frozen=True)
class Link:
    label: str
    node: str
    url: str
    alike: bool
    supports: tuple[str, ...] = ()  # one per Popover.measures column ("–" where none)


@dataclass(frozen=True)
class Popover:
    """What the tree viewer shows on clicking a node: the way to its ancestor, the same
    ancestor in other combinations, and (where the job's ancestor map places it) its
    name, which can be changed there."""

    node: str
    url: str  # this ancestor
    found: list[Link]
    missing: list[str]  # combinations whose trees lack the node
    nameable: bool = False  # in the ancestor map, so it can be named
    called: str = ""  # the user's name for this ancestor (set by the view)
    measures: tuple[str, ...] = ()  # support columns: those any of its rows has
    own: tuple[tuple[str, str], ...] = ()  # (measure, value): this node's own supports

    @property
    def others(self) -> int:
        return len(self.found) + len(self.missing)


def popovers(selection: Any, nodes: set[str]) -> list[Popover]:
    """A popover for each of `nodes` (the drawing's linked nodes) in the tree viewer's
    combination, in the reconstruction's order; the same ancestor elsewhere only where
    the job's ancestor map has this combination."""
    shown = selection.combination["key"]
    placed = trees(selection)
    combos = {c["key"]: c for c in combinations(selection)}
    several = len({c["tree"] for c in combos.values()}) > 1
    clade = by_clade(selection)
    viewer = reverse("ancestor_viewer", args=[selection.job.id])
    matched = analogs.index(placed, shown, clade) if shown in placed else {}
    # Each listed combination's supports.tsv, read together: a few hundred small files
    # on a large job.
    needed = {shown} | {a.combination for found in matched.values() for a in found}
    names = [f"{k}/supports.tsv" for k in sorted(needed) if f"{k}/supports.tsv" in selection.keys]
    supports = {
        name.rsplit("/", 1)[0]: {row["node"]: row for row in results.read_tsv(text)}
        for name, text in selection.results.read_many(names).items()
    }

    def support(key: str, node: str) -> dict[str, str]:
        return supports.get(key, {}).get(node, {})

    def url(key: str, node: str) -> str:
        query = {"aligner": combos[key]["aligner"], "combination": key, "node": node}
        return f"{viewer}?{urlencode(query)}"

    out = []
    for node in selection.node_names:
        if node not in nodes or shown not in combos:
            continue
        found = matched.get(node, [])
        own = support(shown, node)
        rows = [support(a.combination, a.node) for a in found]
        columns = [(c, n) for c, n in POPOVER_SUPPORTS if any(r.get(c) for r in rows)]
        links = [
            Link(
                compact(combos[a.combination], several), a.node, url(a.combination, a.node),
                a.alike, tuple(row.get(c) or "–" for c, _ in columns),
            )
            for a, row in zip(found, rows, strict=True)
        ]  # fmt: skip
        lacking = (
            [label(combos[k], several) for k in analogs.missing(placed, shown, found)]
            if node in matched
            else []
        )
        out.append(
            Popover(
                node, url(shown, node), links, lacking, nameable=node in placed.get(shown, {}),
                measures=tuple(n for _, n in columns),
                own=tuple((n, own[c]) for c, n in POPOVER_SUPPORTS if own.get(c)),
            )
        )  # fmt: skip
    return out


@dataclass(frozen=True)
class Row:
    """One ancestor in the comparison table."""

    label: str
    node: str
    url: str
    this: bool
    alike: bool
    rank: int | None
    length: int | None  # None where the combination lacks this reconstruction
    identity: float | None  # percent; None where it could not be placed beside this one
    anchor: str


COMPARE_PARAM = "compare_on"


def comparison(
    selection: Any, node: str, query: Any, palette: str, threshold: float, extant: str
) -> dict[str, Any] | None:
    """The ancestor viewer's comparison: this ancestor and its analogs, the same
    reconstruction version (gap method) throughout. None when the job has a single
    combination with ancestors, where there is nothing to compare."""
    combos = {c["key"]: c for c in combinations(selection)}
    if len(combos) < 2:
        return None
    shown = selection.combination["key"]
    placed = trees(selection)
    clade = by_clade(selection)
    found = analogs.find(placed, shown, node, clade) if shown in placed else []
    lacking = [combos[k] for k in analogs.missing(placed, shown, found)]
    lacking += [c for k, c in combos.items() if k not in placed and k != shown]
    several = len({c["tree"] for c in combos.values()}) > 1
    version = selection.version
    members = [(shown, node, True)] + [(a.combination, a.node, a.alike) for a in found]
    targets = [{"combination": k, "node": n, "version": version} for k, n, _ in members]
    loaded = library_views.load_targets(selection.job, targets, selection)

    entries: list[library_glyphs.Entry] = []
    for index, ((key, member, _), data) in enumerate(zip(members, loaded, strict=True)):
        if data is None:
            continue
        combo = combos[key]
        this = " ◂" if index == 0 else ""  # the ancestor shown (the table's bold row)
        entries.append(
            library_glyphs.Entry(
                f"{member} · {row_name(combo, several)}{this}", f"analog-{index}",
                f"{member} · {label(combo, several)}", data.template, data.dists,
                combo["aligner"], combo.get("aligner_label", combo["aligner"]), data.alignment,
            )
        )  # fmt: skip

    # Identity is always measured on this ancestor's alignment, whatever the layout shown,
    # so the numbers do not change with the menu.
    identities: dict[str, float] = {}
    if entries and entries[0].anchor == "analog-0":
        reference = library_glyphs.Layout(entries[0].aligner)
        frame = library_glyphs.frames(entries, reference)[0]
        rows = {r.anchor: r for r in frame.rows if r.ancestor}
        mine = rows.get("analog-0")
        lengths = {e.anchor: len(e.template.replace(GAP, "")) for e in entries}
        for anchor, row in rows.items():
            if mine is not None:
                identities[anchor] = analogs.identity(
                    mine.states, row.states, lengths["analog-0"], lengths[anchor]
                ).percent

    viewer = reverse("ancestor_viewer", args=[selection.job.id])
    table = []
    for index, ((key, member, alike), data) in enumerate(zip(members, loaded, strict=True)):
        combo = combos[key]
        here = {"aligner": combo["aligner"], "combination": key, "node": member, "gaps": version}
        anchor = f"analog-{index}"
        length = None if data is None else len(data.template.replace(GAP, ""))
        table.append(
            Row(
                label(combo, several), member, f"{viewer}?{urlencode(here)}", index == 0, alike,
                combo.get("rank_bic"), length, identities.get(anchor), anchor,
            )
        )  # fmt: skip

    layout = None
    frames: list[tuple[library_glyphs.Frame, dict[str, Any]]] = []
    options = library_glyphs.alignment_options(entries) if entries else []
    if entries:
        chosen = str(query.get(COMPARE_PARAM, ""))
        if chosen not in dict(options):
            chosen = entries[0].aligner
        layout = library_glyphs.resolve_layout(entries, chosen, [extant])
        least = max(library_glyphs.threshold_percent(threshold), library_glyphs.FLOOR)
        cell = library_glyphs.parse_zoom(query.get("zoom", "fit"))
        drawn = library_glyphs.frames(entries, layout)
        frames = [
            (f, f.data(palette, least, cell))
            for f in (dataclasses.replace(f, index=i, group="c") for i, f in enumerate(drawn))
        ]
    addable = sum(1 for data in loaded if data is not None)
    return {
        "rows": table,
        "found": len(found),
        "others": len(combos) - 1,
        "missing": [label(c, several) for c in lacking],
        "by_clade": clade,
        "add_label": f"Add all {addable} to" if addable > 1 else "Add to",
        "note": match_note(clade),
        "frames": frames,
        "layout": layout.alignment if layout else "",
        "layouts": [
            (
                key,
                "Each alignment separately"
                if key == library_glyphs.OWN
                else f"On {extant}'s residues"
                if key == library_glyphs.ON_SEQUENCE
                else f"On the {text} alignment",
            )
            for key, text in options
        ],  # fmt: skip
        "param": COMPARE_PARAM,
        "zoom": query.get("zoom", "fit"),
        "version_label": ancestors.version_label(version),
    }


@require_POST
def add_to_library(request: HttpRequest, job_id: str) -> HttpResponse:
    """Add an ancestor and every analog of it (POST combination, gaps, node) to a library,
    each in the same reconstruction version where its combination has it. The analogs are
    found again here, never taken from the page."""
    from phylobot.views import select_ancestor  # views imports this module

    job = library_views.owned_job(request, job_id)
    selection = select_ancestor(job, request.POST)
    if selection is None:
        raise Http404
    if selection.combination["key"] != request.POST.get(
        "combination"
    ) or selection.version != request.POST.get("gaps"):
        raise Http404
    node = request.POST.get("node", "")
    if node not in set(selection.node_names) & set(selection.ancestors):
        raise Http404
    shown = selection.combination["key"]
    found = analogs.find(trees(selection), shown, node, by_clade(selection))
    entries = []
    for key, member in [(shown, node)] + [(a.combination, a.node) for a in found]:
        other = selection.other(key, selection.version)
        if other is None or member not in other.ancestors:
            continue  # this combination lacks the version (codeml's, say)
        identity = library_views.node_identity(other, member)
        entries.append(library.target_entry(key, member, selection.version, *identity))
    named = f"{node} and its {len(entries) - 1} analog{'s' if len(entries) != 2 else ''}"
    return library_views.save_targets(request, job, entries, named)
