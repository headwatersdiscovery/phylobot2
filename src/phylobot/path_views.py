"""The mutational path page (`/jobs/<id>/path/`).

A path traced on the tree viewer opens here, a page of its own: the path drawn as a
straight line (`phylobot.path_line`), its changes branch by branch, the ancestors along
it as alignment viewers, and all of it compared with the job's other alignment x model
combinations (`phylobot.path_compare`), because a path read from one alignment can be
wrong where alignments differ, and the changes several alignments agree about are the ones to
trust. A path can be pinned (`phylobot.saved_paths`) to come back to.

Like the other viewers it is one GET form, so every view has an address; the POST
buttons (library, pin) go through page-level forms. The path's two ends are node names
of the chosen combination's tree: when the combination changes, `of` says which tree
they were names in and they are carried over by clade (`translate`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import pairwise
from typing import Any
from urllib.parse import parse_qs, urlencode

from django.contrib import messages
from django.core.cache import cache
from django.http import Http404, HttpRequest, HttpResponse, QueryDict
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.views.decorators.http import require_GET, require_POST

from phylobot import (
    analog_views,
    ancestor_names,
    ancestors,
    gaps,
    library_glyphs,
    path_compare,
    path_line,
    paths,
    results,
    saved_paths,
    views,
)
from phylobot.models import Job

# How the alignment viewers show the path's changes.
MUTATIONS = {
    "only": "Colour only the changes",
    "support": "Colour the changes by how widely they are shared",
    "red": "Colour every site, outline the changes in red",
    "white": "Colour every site, outline the changes in white",
}
DEFAULT_MUTATIONS = "only"
LAYOUTS = {"stack": "Stacked", "side": "Side by side"}  # the keys name their icons too
VIEWER_LABELS = {"all": "in every alignment × model combo", "most": "in most combos"}
SHOWN = {
    "all": "Every change",
    "shared": "Made by every compared combination",
    "disputed": "Not made by every one",
}
# Which combinations a change is checked against: each alignment's best by BIC.
TOPS = {"1": "Each alignment's best model", "3": "Each alignment's best 3 models"}
EVERY = "all"
DEFAULT_TOP = "1"
MAX_COMPARED = 40  # "every combination" is offered up to this many (285 on the GFP job)
MAX_ALSO = 6  # alignments shown beside the traced one
MAX_EXTANT = 20
DEFAULT_ZOOM = 14.0
NAMES_WIDTH = {"stack": 190, "side": 130}
ROW = 18
CACHE_SECONDS = 3600
BACK_LIMIT = 2000
NOT_ON_PATH = len(path_compare.LEVELS)  # a followed combination's change the traced path lacks
LEVEL_COLOURS = [*(path_line.LEVEL_COLOURS[level] for level in path_compare.LEVELS), "#6b7280"]


def _choice(query: Any, name: str, options: Any, default: str) -> str:
    value = str(query.get(name, ""))
    return value if value in options else default


def translate(selection: Any, query: Any) -> tuple[str, str, list[str]]:
    """The path's two ends as the chosen combination's tree names them. `of` is the
    combination the names in the address belong to; where that is another one, each
    ancestor is carried over by its clade (the job's ancestor map). Returns the ends and
    the names that could not be carried over (the chosen tree has no such ancestor)."""
    start, end = query.get("path_from", ""), query.get("path_to", "")
    key = selection.combination["key"]
    of = query.get("of", "")
    if not of or of == key:
        return start, end, []
    placed = analog_views.trees(selection)
    theirs, mine = placed.get(of, {}), placed.get(key, {})
    by_clade = analog_views.by_clade(selection)
    lookup = {(p.below if by_clade else p.parts): node for node, p in mine.items()}
    lost = []

    def carried(name: str) -> str:
        if not name or name not in theirs:
            return name  # an extant sequence (or nothing): the same name everywhere
        place = theirs[name]
        twin = lookup.get(place.below if by_clade else place.parts, "")
        if not twin:
            lost.append(name)
        return twin

    return carried(start), carried(end), lost


@dataclass
class Choice:
    """The page's view, every value checked against the job."""

    base: dict[str, Any]  # views.path_context
    route: list[str]
    lost: list[str]
    top: str
    compared: list[str]  # combinations each change is checked against
    also: list[str]  # combinations shown beside the traced one
    layout: str
    extant: list[str]
    mutations: str
    palette: str
    zoom: float | str
    lengths: str
    shown: str
    back: str  # the tree viewer's own settings, as a query string
    combos: dict[str, dict[str, Any]] = field(default_factory=dict)

    def view(self, selection: Any, back: bool = True) -> dict[str, Any]:
        """The canonical address of this view."""
        combo = selection.combination
        base = self.base
        view: dict[str, Any] = {"aligner": combo["aligner"], "combination": combo["key"]}
        view["gaps"] = selection.version
        if base["path_from"]:
            view["path_from"] = base["path_from"]
        if base["path_to"]:
            view["path_to"] = base["path_to"]
        view |= {"top": self.top, "lengths": self.lengths, "show": self.shown}
        if self.also:
            view["also"] = self.also
        view |= {"layout": self.layout, "mutations": self.mutations, "palette": self.palette}
        view["zoom"] = "fit" if self.zoom == "fit" else f"{self.zoom:g}"
        if self.extant:
            view["extant"] = self.extant
        view["cols"] = "changes" if base["only_changes"] else "all"
        if base["has_codons"]:
            view |= {"silent": "show" if base["show_silent"] else "hide", "cells": base["cells"]}
        if back and self.back:
            view["back"] = self.back
        return view


def best_per_alignment(combos: list[dict[str, Any]], top: int | None) -> list[str]:
    """Each alignment's `top` combinations by BIC (all of them for None), in menu order."""
    if top is None:
        return [c["key"] for c in combos]
    kept = []
    for aligner in dict.fromkeys(c["aligner"] for c in combos):
        mine = [c for c in combos if c["aligner"] == aligner]
        ranked = sorted(mine, key=lambda c: (c.get("rank_bic") or 10**6, mine.index(c)))
        kept += [c["key"] for c in ranked[:top]]
    wanted = set(kept)
    return [c["key"] for c in combos if c["key"] in wanted]


def clean_back(text: str) -> str:
    """The tree viewer's display settings out of a "back" value (a query string): only
    its own known keys, re-encoded, so the link built from it can lead nowhere else."""
    parsed = parse_qs(str(text)[:BACK_LIMIT].lstrip("?"), keep_blank_values=False)
    kept = {k: v[-1] for k, v in parsed.items() if k in saved_paths.TREE_VIEW_KEYS}
    return urlencode(kept)


def choose(selection: Any, query: Any, table: bool = True) -> Choice:
    """Resolve the address against the job. The one-off buttons are applied here: `add`
    (a combination), `add_best` (every other alignment's best) and `drop`."""
    key = selection.combination["key"]
    start, end, lost = translate(selection, query)
    asked = query.copy()
    asked["path_from"], asked["path_to"], asked["trace"] = start, end, "1"
    names = list(ancestors.parse_fasta(selection.read("alignment.fasta")))
    extant = list(dict.fromkeys(n for n in query.getlist("extant") if n in names))[:MAX_EXTANT]
    asked.setlist("compare", extant)
    base = views.path_context(selection, asked, selection.tree(), with_table=table)

    combos = {c["key"]: c for c in analog_views.combinations(selection)}
    others = [c for k, c in combos.items() if k != key]
    top = _choice(query, "top", {**TOPS, EVERY: ""}, DEFAULT_TOP)
    if top == EVERY and len(others) > MAX_COMPARED:
        top = DEFAULT_TOP
    pool = best_per_alignment(list(combos.values()), None if top == EVERY else int(top))
    compared = [k for k in pool if k != key]

    also = list(dict.fromkeys(k for k in query.getlist("also") if k in combos and k != key))
    add = query.get("add", "")
    if query.get("add_best"):  # each other alignment's best, in one go
        mine = selection.combination["aligner"]
        elsewhere = [c for c in others if c["aligner"] != mine]
        also += [k for k in best_per_alignment(elsewhere, 1) if k not in also]
    elif add in combos and add != key and add not in also:
        also.append(add)
    dropped = query.get("drop", "")
    also = [k for k in also if k != dropped][:MAX_ALSO]
    return Choice(
        base=base,
        route=base["path"],
        lost=lost,
        top=top,
        compared=compared,
        also=also,
        layout=_choice(query, "layout", LAYOUTS, "stack"),
        extant=extant,
        mutations=_choice(query, "mutations", MUTATIONS, DEFAULT_MUTATIONS),
        palette=base["palette"],
        zoom=library_glyphs.parse_zoom(query.get("zoom", DEFAULT_ZOOM)),
        lengths=_choice(query, "lengths", ("real", "equal"), "real"),
        shown=_choice(query, "show", SHOWN, "all"),
        back=clean_back(query.get("back", "")),
        combos=combos,
    )


def load_sides(selection: Any, keys: list[str]) -> dict[str, path_compare.Side]:
    """Each combination's reconstruction in the selection's gap version, for those that
    have it. Read together, and kept for an hour: a finished job's files never change,
    and every control on the page asks for the same ones again."""
    version = selection.version
    files = selection.results
    placed = analog_views.trees(selection)
    combos = {c["key"]: c for c in analog_views.combinations(selection)}
    found: dict[str, tuple[dict[str, str], dict[str, str], str]] = {}
    wanted: dict[str, tuple[str, str, str]] = {}
    for key in keys:
        fasta = views.version_files(key).get(version, "")
        alignment = f"{key}/alignment.fasta"
        rooted = f"{key}/tree.rooted.nwk"
        tree = rooted if rooted in selection.keys else f"{key}/tree.nwk"
        if key not in combos or not {fasta, alignment, tree} <= selection.keys:
            continue
        cached = cache.get(f"path-side:v1:{selection.job.id}:{key}:{version}")
        if cached is not None:
            found[key] = cached
        else:
            wanted[key] = (fasta, alignment, tree)
    if wanted:
        texts = files.read_many(list(dict.fromkeys(n for names in wanted.values() for n in names)))
        for key, (fasta, alignment, tree) in wanted.items():
            found[key] = (
                ancestors.parse_fasta(texts[alignment]),
                ancestors.parse_fasta(texts[fasta]),
                texts[tree],
            )
            cache.set(f"path-side:v1:{selection.job.id}:{key}:{version}", found[key], CACHE_SECONDS)
    return {
        key: path_compare.Side(key, combos[key]["aligner"], *found[key], places=placed.get(key, {}))
        for key in keys
        if key in found
    }


def probabilities(selection: Any, nodes: dict[str, list[str]]) -> dict[str, dict[str, list[float]]]:
    """{combination: {node: its ML residues' posterior probabilities}} for the nodes asked
    for, where the gap version keeps IQ-TREE's residues (PRANK and GRASP have none)."""
    version = selection.version
    if version not in (ancestors.ORIGINAL, *gaps.MASKING):
        return {}
    table = ancestors.version_folder(version) + "ancestors_pp.tsv"
    names = {key: f"{key}/{table}" for key in nodes if f"{key}/{table}" in selection.keys}
    texts = selection.results.read_many(list(names.values()))
    return {
        key: paths.ml_probabilities(texts[name], set(nodes[key])) for key, name in names.items()
    }


def change_label(call: path_compare.Call, column: int) -> str:
    return f"{call[0]}{column + 1}{call[1]}"


def frame(
    side: path_compare.Side,
    rows: list[str],
    labels: dict[str, str],
    pp: dict[str, list[float]],
    extant: list[str],
    choice: Choice,
    levels: dict[str, dict[int, int]],
    notes: dict[str, dict[int, str]],
    sync: list[int] | None,
    match: list[int] | None = None,
) -> dict[str, Any]:
    """One alignment viewer's data (alignment.js): `rows` are the path's nodes in this
    combination, oldest first ("" for a row kept empty beside another combination's
    ancestor); the chosen extant sequences follow. A cell is marked where it differs from
    the path row above it."""
    width = len(next(iter(side.alignment.values()), ""))
    data: dict[str, list[Any]] = {
        k: [] for k in ("names", "ancestor", "rows", "pp", "marked", "onPath", "blank",
                        "levels", "notes")
    }  # fmt: skip
    previous = ""
    previous_name = ""
    for name in rows:
        states = side.row(name) if name else None
        if states is None:
            states = path_compare.GAP * width
        found = pp.get(name)
        usable = found is not None and len(found) == width and name in side.ancestors
        marked = []
        row_notes = dict(notes.get(name, {}))
        if name and previous:
            for column, call in path_compare.differences(previous, states).items():
                marked.append(column)
                row_notes.setdefault(
                    column,
                    f"{change_label(call, column)} from {labels.get(previous_name, previous_name)}",
                )
        data["names"].append(labels.get(name, name) if name else "–")
        data["ancestor"].append(bool(name) and name in side.ancestors)
        data["rows"].append(states)
        data["pp"].append(
            [library_glyphs.percent(p) if s != path_compare.GAP else -1
             for p, s in zip(found or [], states, strict=True)] if usable else []
        )  # fmt: skip
        data["marked"].append(marked)
        data["onPath"].append(bool(name))
        data["blank"].append(not name)
        data["levels"].append({str(c): v for c, v in levels.get(name, {}).items()})
        data["notes"].append({str(c): v for c, v in row_notes.items()})
        if name:
            previous, previous_name = states, name
    for name in extant:
        if name not in side.alignment or name in rows:
            continue
        data["names"].append(name)
        data["ancestor"].append(False)
        data["rows"].append(side.alignment[name].upper())
        for k in ("pp", "marked"):
            data[k].append([])
        data["onPath"].append(False)
        data["blank"].append(False)
        for k in ("levels", "notes"):
            data[k].append({})
    outlined = choice.mutations in ("red", "white")
    return data | {
        "amino": [],
        "codon": False,
        "extant": sum(1 for a in data["ancestor"] if not a),
        "reference": -1,
        "scheme": "posterior" if outlined else "changes",
        "cell": choice.zoom,
        "row": ROW,
        "height": len(data["rows"]),
        "tracks": [],  # none: side by side the rows must line up, and stacked they add little
        "ancestorViewer": None,
        "palette": library_glyphs.palette_stops(choice.palette),
        "threshold": 101,
        "markColour": {"red": "#d7191c", "white": "#ffffff"}.get(choice.mutations, ""),
        "markWidth": 2.5 if outlined else 0,
        "levelColours": LEVEL_COLOURS if choice.mutations == "support" else [],
        "namesWidth": NAMES_WIDTH[choice.layout],
        "sync": sync,
        # Per column, the traced alignment's column holding the same residues (-1: none);
        # path.js boxes the same site in every viewer as the pointer moves.
        "match": match,
        "hoverEvents": True,
    }


def support_title(found: path_compare.Support, labels: dict[str, str]) -> str:
    """A change's support in words, for its tooltip."""
    if not found.compared:
        return "not compared: no other combination has this branch"
    parts = [f"the same in {found.agreeing} of {found.compared}"]
    for state, word in ((path_compare.SAME, "same"), (path_compare.DIFFERENT, "different"),
                        (path_compare.ABSENT, "no change")):  # fmt: skip
        here = [v for v in found.verdicts if v.state == state]
        if not here or (state == path_compare.SAME and len(here) == found.compared):
            continue
        shown = []
        for v in here[:6]:
            text = labels.get(v.side, v.side)
            shown.append(f"{text} ({v.call[0]}→{v.call[1]})" if v.call else text)
        more = f" and {len(here) - 6} more" if len(here) > 6 else ""
        parts.append(f"{word}: {', '.join(shown)}{more}")
    return " · ".join(parts)


def level_words(counts: dict[str, int]) -> str:
    """ "9 in every one, 2 in most, 1 only here"."""
    return ", ".join(
        f"{counts[level]} {path_compare.LEVEL_LABELS[level]}"
        for level in path_compare.LEVELS
        if counts.get(level)
    )


def context_for(request: HttpRequest, selection: Any, choice: Choice) -> dict[str, Any]:
    """Everything the page shows for a traced path."""
    job = selection.job
    key = selection.combination["key"]
    route = choice.route
    combos = choice.combos
    several = len({c["tree"] for c in combos.values()}) > 1
    label = {k: analog_views.label(c, several) for k, c in combos.items()}
    short = {k: analog_views.compact(c, several) for k, c in combos.items()}
    by_clade = analog_views.by_clade(selection)

    needed = list(dict.fromkeys([key, *choice.compared, *choice.also]))
    sides = load_sides(selection, needed)
    main = sides[key]
    others = [sides[k] for k in needed[1:] if k in sides]
    comparison = path_compare.compare(main, route, others, by_clade)
    followed = {f.side.key: f for f in comparison.followed}
    slot = {f.side.key: i for i, f in enumerate(comparison.followed)}
    judged = [k for k in choice.compared if k in sides]  # the ones the counts are over

    # Support counted over the compared combinations only: one shown beside the path but
    # outside that set is followed too, and must not move the numbers.
    def narrowed(found: path_compare.Support) -> path_compare.Support:
        return path_compare.Support(tuple(v for v in found.verdicts if v.side in set(judged)))

    support = [{c: narrowed(s) for c, s in branch.items()} for branch in comparison.support]
    level_index = {level: i for i, level in enumerate(path_compare.LEVELS)}

    called = {k: ancestor_names.for_selection(selection, k) for k in [key, *choice.also]}

    def named(node: str, side: str = key) -> str:
        return ancestor_names.label(node, called.get(side, {}))

    # --- The changes, branch by branch (the path table's, each with its support).
    table = choice.base["path_table"]
    totals = dict.fromkeys(path_compare.LEVELS, 0)
    for_excel: dict[tuple[int, int], tuple[int, int, str]] = {}
    branches = []
    for i, branch in enumerate(table.branches):
        counts = dict.fromkeys(path_compare.LEVELS, 0)
        chips = []
        for change in branch.changes:
            found = support[i].get(change.column - 1)
            # A silent or uncertain codon change is no amino-acid change: not compared.
            if change.kind not in ("substitution", "gain", "loss"):
                found = None
            level = found.level if found is not None else ""
            title = support_title(found, short) if found is not None else ""
            if found is not None:
                counts[level] += 1
                for_excel[i, change.column] = (found.agreeing, found.compared, title)
            hidden = (choice.shown == "shared" and level != path_compare.ALL) or (
                choice.shown == "disputed" and level in ("", path_compare.ALL, path_compare.NA)
            )
            chips.append({
                "change": change,
                "level": level,
                "agreeing": found.agreeing if found is not None else None,
                "compared": found.compared if found is not None else None,
                "title": title,
                "hidden": hidden,
            })  # fmt: skip
        for level, n in counts.items():
            totals[level] += n
        judged_here = sum(counts.values())
        branches.append({
            "branch": branch,
            "chips": [c for c in chips if not c["hidden"]],
            "hidden": sum(1 for c in chips if c["hidden"]),
            "counts": counts,
            "words": level_words(counts),
            "bar": [
                (level, path_line.LEVEL_COLOURS[level], round(100 * n / judged_here, 2), n)
                for level, n in counts.items() if n
            ],
        })  # fmt: skip

    # --- The line.
    lengths = paths.branch_lengths(main.newick, route)
    viewer = reverse("ancestor_viewer", args=[job.id])
    here = {"aligner": selection.combination["aligner"], "combination": key,
            "gaps": selection.version}  # fmt: skip
    urls = [
        f"{viewer}?{urlencode(here | {'node': n})}" if n in main.ancestors else "" for n in route
    ]
    titles = []
    for i, b in enumerate(branches):
        n = len(table.branches[i].changes)
        words = f": {b['words']}" if b["words"] and judged else ""
        titles.append(
            f"{named(route[i])} → {named(route[i + 1])} · {lengths[i]:.4g} substitutions per "
            f"site · {n} change{'s' if n != 1 else ''}{words}"
        )
    line = path_line.line(
        route, [named(n) for n in route], set(main.alignment), urls, lengths,
        [b["counts"] if judged else {"na": sum(b["counts"].values())} for b in branches],
        titles, equal=choice.lengths == "equal",
    )  # fmt: skip

    # --- The alignment viewers.
    shown = [k for k in choice.also if k in sides]
    lanes = path_compare.merged_rows(route, [followed[k] for k in shown])
    side_by_side = choice.layout == "side" and bool(shown)
    wanted = {key: route} | {k: followed[k].route for k in shown}
    pp = probabilities(selection, wanted)
    main_levels: dict[str, dict[int, int]] = {}
    main_notes: dict[str, dict[int, str]] = {}
    for i, (parent, child) in enumerate(pairwise(route)):
        for column, call in comparison.changes[i].items():
            found = support[i][column]
            main_levels.setdefault(child, {})[column] = level_index[found.level]
            text = f"{change_label(call, column)} on {named(parent)} → {named(child)}"
            if found.compared:
                text += f"\nthe same in {found.agreeing} of {found.compared} compared"
            main_notes.setdefault(child, {})[column] = text
    canonical = choice.view(selection)
    page = reverse("path_viewer", args=[job.id])

    def url(**changes: Any) -> str:
        view = {k: v for k, v in (canonical | changes).items() if v not in ("", None, [])}
        return f"{page}?{urlencode(view, doseq=True)}"

    frames = [{
        "key": key,
        "title": label[key],
        "main": True,
        "columns": len(main.row(route[0]) or ""),
        "json_id": "aln-data-p0",
        "data": frame(
            main, [r[0] for r in lanes] if side_by_side else route,
            {n: named(n) for n in route}, pp.get(key, {}), choice.extant, choice,
            main_levels, main_notes, None,
        ),
    }]  # fmt: skip
    anchor = next((n for n in choice.extant if n in main.alignment), "")
    for index, other_key in enumerate(shown, start=1):
        track = followed[other_key]
        side = track.side
        levels: dict[str, dict[int, int]] = {}
        notes: dict[str, dict[int, str]] = {}
        for upper, lower in pairwise(track.route):
            for column, call in path_compare.differences(
                side.row(upper) or "", side.row(lower) or ""
            ).items():
                hit = track.shared.get(lower, {}).get(column)
                if hit is None:
                    levels.setdefault(lower, {})[column] = NOT_ON_PATH
                    continue
                branch, main_column = hit
                levels.setdefault(lower, {})[column] = level_index[
                    support[branch][main_column].level
                ]
                notes.setdefault(lower, {})[column] = (
                    f"{change_label(call, column)}: on the traced path too "
                    f"({change_label(call, main_column)})"
                )
        width = len(next(iter(side.alignment.values()), ""))
        sync = None
        if anchor and anchor in side.alignment:
            sync = path_compare.sync_on(side.alignment[anchor], main.alignment[anchor])
        if sync is None and track.columns is not None:
            sync = path_compare.sync_columns(track.columns, width, frames[0]["columns"])
        matches = [b[slot[other_key]] for b in comparison.branches]
        total = sum(m.total for m in matches if m.comparable)
        frames.append({
            "key": other_key,
            "title": label[other_key],
            "main": False,
            "columns": width,
            "json_id": f"aln-data-p{index}",
            "trace_url": url(combination=other_key, aligner=side.aligner, of=key,
                             also=[k for k in [key, *choice.also] if k != other_key][:MAX_ALSO]),
            "missing": [named(n) for n in track.missing],
            "placed": track.columns is not None,
            "same": sum(m.same for m in matches if m.comparable),
            "total": total,
            "data": frame(
                side, [r[index] for r in lanes] if side_by_side else track.route,
                {n: ancestor_names.label(n, called.get(other_key, {})) for n in track.route},
                pp.get(other_key, {}), choice.extant, choice, levels, notes, sync,
                None if track.columns is None
                else [track.columns.get(c, -1) for c in range(width)],
            ),
        })  # fmt: skip

    # --- The support table: a branch per row, a compared combination per column.
    matrix = []
    for i, b in enumerate(branches):
        cells = []
        for other_key in judged:
            m = comparison.branches[i][slot[other_key]]
            ratio = m.same / m.total if m.comparable and m.total else None
            cells.append({
                "match": m,
                "style": results.tint(ratio) if ratio is not None else "",
                "title": (
                    f"{short[other_key]}: {m.same} the same, {m.different} different, "
                    f"{m.absent} not there; {m.extra} more only there"
                    + ("" if m.direct else " (other ancestors lie between the two in its tree)")
                ) if m.comparable else (
                    f"{short[other_key]}: its tree has no "
                    + " or ".join(named(n) for n, twin in
                                  ((route[i], m.parent), (route[i + 1], m.child)) if not twin)
                    if not (m.parent and m.child)
                    else f"{short[other_key]}: its alignment holds other sequences"
                ),
            })  # fmt: skip
        matrix.append({"branch": b["branch"], "cells": cells, "index": i})
    column_totals = []
    for other_key in judged:
        matches = [b[slot[other_key]] for b in comparison.branches]
        comparable = [m for m in matches if m.comparable]
        same, total = sum(m.same for m in comparable), sum(m.total for m in comparable)
        column_totals.append({
            "same": same, "total": total, "branches": len(comparable),
            "style": results.tint(same / total) if total else "",
        })  # fmt: skip

    saved = saved_paths.find(job, key, route[0], route[-1])
    compared_note = (
        f"{len(judged)} other combination{'s' if len(judged) != 1 else ''}"
        if judged else "no other combination"
    )  # fmt: skip
    aligners, models = views.combination_menus(selection)
    tree = reverse("tree_viewer", args=[job.id])
    back = dict(parse_qs(choice.back)) | here | {
        "path_from": route[0], "path_to": route[-1], "trace": "1"}  # fmt: skip
    judged_total = sum(totals.values()) - totals[path_compare.NA]
    return choice.base | {
        "choice": choice,
        "route_labels": [named(n) for n in route],
        "branches": branches,
        "line": line,
        "frames": frames,
        "side_by_side": side_by_side,
        "matrix": matrix,
        "matrix_columns": [
            {"key": k, "label": short[k], "total": t}
            for k, t in zip(judged, column_totals, strict=True)
        ],
        "totals": totals,
        "total_words": level_words(totals),
        "judged_total": judged_total,
        "compared_note": compared_note,
        "compared": len(judged),
        "levels": [
            (level, path_line.LEVEL_COLOURS[level], path_compare.LEVEL_LABELS[level])
            for level in path_compare.LEVELS
        ],
        # The viewers' legend spells the first two out; the rest read on from them.
        "viewer_levels": [
            (level, path_line.LEVEL_COLOURS[level], VIEWER_LABELS.get(level, text))
            for level, text in path_compare.LEVEL_LABELS.items()
        ],
        "by_clade": by_clade,
        "not_on_path_colour": LEVEL_COLOURS[NOT_ON_PATH],
        "add_options": [
            (
                aligner_label,
                [
                    (k, analog_views.short_label(c, several))
                    for k, c in combos.items()
                    if c["aligner"] == aligner and k != key and k not in choice.also
                ],
            )
            for aligner, aligner_label in aligners.items()
        ],
        "can_add": len(choice.also) < MAX_ALSO and len(combos) > 1 + len(choice.also),
        "unshown_best": [
            k
            for k in best_per_alignment([c for k, c in combos.items() if k != key], 1)
            if k not in choice.also and combos[k]["aligner"] != selection.combination["aligner"]
        ],
        "sequence_names": list(main.alignment),
        "extant": choice.extant,
        "anchor": anchor,
        "aligners": aligners,
        "aligner": selection.combination["aligner"],
        "models": models,
        "combination": key,
        "tops": [*TOPS.items()]
        + (
            [(EVERY, f"Every combination ({len(combos) - 1})")]
            if len(combos) - 1 <= MAX_COMPARED
            else []
        ),
        "mutation_choices": list(MUTATIONS.items()),
        "layouts": list(LAYOUTS.items()),
        "shown_choices": list(SHOWN.items()),
        "library_path_nodes": [n for n in route if n in selection.node_names],
        "saved": saved,
        "path_support": for_excel if judged else None,
        "compared_with": "; ".join(label[k] for k in judged),
        "view_query": urlencode(choice.view(selection, back=False), doseq=True),
        "canonical_url": url(),
        "path_excel_url": url(export="xlsx"),
        "tree_url": f"{tree}?{urlencode(back, doseq=True)}",
        "match_note": analog_views.match_note(by_clade),
    }


@require_GET
def page(request: HttpRequest, job_id: str) -> HttpResponse:
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    selection = views.select_ancestor(job, request.GET)
    htmx = bool(request.headers.get("HX-Request"))
    # With HTMX the way back to the tree comes again out of band (_path_back.html).
    context: dict[str, Any] = {"job": job, "selection": selection, "oob": htmx}
    if selection is not None:
        choice = choose(selection, request.GET)
        context["saved_paths"] = saved_paths.listing(job)
        tree = reverse("tree_viewer", args=[job.id])
        here = {"aligner": selection.combination["aligner"],
                "combination": selection.combination["key"]}  # fmt: skip
        context |= {"choice": choice, "tree_url": f"{tree}?{urlencode(here)}",
                    "lost": choice.lost}  # fmt: skip
        if choice.route:
            context |= context_for(request, selection, choice)
            if request.GET.get("export") == "xlsx":
                return views.path_workbook(context)
            moved = request.GET.get("of", context["combination"]) != context["combination"]
            one_off = moved or any(request.GET.get(k) for k in ("add", "add_best", "drop"))
            if one_off and not htmx:
                return redirect(context["canonical_url"])
        else:
            aligners, models = views.combination_menus(selection)
            context |= {"aligners": aligners, "models": models,
                        "aligner": here["aligner"], "combination": here["combination"]}  # fmt: skip
    template = "_path_viewer.html" if htmx else "path.html"
    response = render(request, f"phylobot/{template}", context)
    if htmx and context.get("canonical_url"):
        response["HX-Push-Url"] = context["canonical_url"]
    return response


def _posted_view(request: HttpRequest, job: Job) -> tuple[Any, Choice]:
    """The path a POST names (its `view`: the page's own query), checked like any address."""
    query = QueryDict(request.POST.get("view", "")[:8000])
    selection = views.select_ancestor(job, query)
    if selection is None:
        raise Http404
    choice = choose(selection, query, table=False)
    if not choice.route:
        raise Http404
    return selection, choice


@require_POST
def save(request: HttpRequest, job_id: str) -> HttpResponse:
    """Pin the path shown (or rename it): POST `view` (the page's query) and `name`."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    selection, choice = _posted_view(request, job)
    key = selection.combination["key"]
    problem = saved_paths.save(
        job, request.POST.get("name", ""), key, choice.route[0], choice.route[-1],
        choice.view(selection, back=False),
    )  # fmt: skip
    if problem:
        messages.error(request, problem)
    else:
        messages.success(request, "Path saved: it is listed here and on the tree viewer.")
    page_url = reverse("path_viewer", args=[job.id])
    return redirect(f"{page_url}?{urlencode(choice.view(selection), doseq=True)}")


@require_POST
def unsave(request: HttpRequest, job_id: str, saved_id: int) -> HttpResponse:
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    saved_paths.remove(job, saved_id)
    target = request.POST.get("next", "")
    if not (target.startswith(f"/jobs/{job.id}/") and "//" not in target and "\\" not in target):
        target = reverse("tree_viewer", args=[job.id])
    return redirect(target)
