"""The structure page (`/jobs/<id>/structure/`): ancestors on protein structures.

A job has structures attached (`Structure`: a PDB entry, an AlphaFold model, or an
uploaded file); the page draws them in the browser (3Dmol.js, `structure.js`) painted
with the job's ancestors.

The page has two sequences at most, as a branch or a mutational path has two ends: an
ancestor, and optionally a second one to compare it with (another ancestor, or an extant
sequence). It shows them on one to four **panels**, each a structure painted with one of
the two, so sequences and structures can be mixed: both ends on one structure, each end on
the structure closest to it, one ancestor on every structure. Panels of different
structures are superposed through the job's alignment (`structures.superpose`) and share a
camera.

Like the other viewers the page is one GET form, so every view has an address. Each
panel's 3D stage is kept across the form's HTMX swaps (`hx-preserve`), so changing the
ancestor repaints the structures without reloading them or moving the camera. Attaching,
uploading and removing structures are POST forms outside it.

The web tier does no modelling: it reads a structure's residues, aligns each chain with
the job's extant sequence most like it (`structures.mapping`) and looks up what an
ancestor has at each residue's alignment column. Coordinates are never changed, only
moved rigidly to lie on the first panel's.
"""

from __future__ import annotations

import io
import re
import uuid
import zipfile
from dataclasses import dataclass, field
from itertools import pairwise
from pathlib import PurePosixPath
from typing import Any
from urllib.parse import parse_qs, urlencode, urlsplit

from django.contrib import messages
from django.core.cache import cache
from django.http import Http404, HttpRequest, HttpResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.views.decorators.http import require_GET, require_POST

from phylobot import (
    analog_views,
    analogs,
    ancestor_names,
    ancestors,
    gaps,
    paths,
    storage,
    views,
)
from phylobot.forms import SequenceLookupForm
from phylobot.homologs import lookup
from phylobot.homologs.remote import SearchError
from phylobot.models import Job, Structure
from phylobot.structures import colours, mapping, pdb, remote, superpose

MAX_PER_JOB = 20
MAX_AT_ONCE = 8  # structures attached by one submission
MAX_PANELS = 4  # each is a WebGL context in the browser
MAX_ZIP_BYTES = 60 * 1024 * 1024
MAX_ZIP_MEMBERS = 500
CACHE_SECONDS = 3600
A, B = "a", "b"  # the page's two sequences, as a panel names the one it paints
STYLES = {
    "cartoon": "Cartoon",
    "trace": "Backbone trace",
    "stick": "Sticks",
    "sphere": "Spheres",
}
SIDECHAINS = {
    "marked": "At the coloured differences",
    "all": "Everywhere",
    "none": "Nowhere",
}
BACKGROUNDS = {"auto": "Follow the skin", "white": "White", "black": "Black"}
PRESETS = {
    "one": "One structure, the changes painted on it",
    "shared": "Both on one structure",
    "closest": "Each on the structure closest to it",
    "every": "One sequence on every structure",
}
# The page's own query keys: what a POST carries back to the page it came from.
VIEW_KEYS = (
    "aligner", "combination", "gaps", "node", "to", "shown", "ps", "pw", "pc", "scheme",
    "palette", "style", "sidechains", "surface", "ligands", "labels", "spin", "background",
    "link", "scroll",
)  # fmt: skip
PANEL_KEYS = ("ps", "pw", "pc")
MAX_CHANGES = 400  # rows of the changes table
MAX_LABELS = 60  # labels drawn in a 3D view
STRUCTURE_FILES = (".cif", ".mmcif", ".pdb", ".ent")


def _choice(query: Any, name: str, options: Any, default: str) -> str:
    value = str(query.get(name, ""))
    return value if value in options else default


def _on(query: Any, name: str, default: bool) -> bool:
    """A checkbox with a hidden "off" before it (an unticked box submits nothing)."""
    value = query.get(name, "")
    return default if value not in ("on", "off") else value == "on"


def attached(job: Job) -> list[Structure]:
    return list(Structure.objects.filter(job=ancestor_names.holder(job)))


def parsed(structure: Structure) -> pdb.Parsed:
    """The structure's chains. Its file never changes, so the reading is kept."""
    key = f"structure-parsed:v2:{structure.id}"
    found = cache.get(key)
    if found is None:
        found = pdb.parse(storage.read_bytes(structure.storage_key))
        cache.set(key, found, CACHE_SECONDS)
    assert isinstance(found, pdb.Parsed)
    return found


def chain_maps(
    structure: Structure, chains: pdb.Parsed, aligner: str, alignment: dict[str, str]
) -> dict[str, mapping.ChainMap]:
    """Each chain's reference sequence and alignment with it, kept per alignment (the
    job's alignments hold the same sequences, but an extension's may not)."""
    key = f"structure-map:v2:{structure.id}:{structure.job_id}:{aligner}"
    found = cache.get(key)
    if found is None:
        found = mapping.chain_maps(chains, mapping.ungapped(alignment))
        cache.set(key, found, CACHE_SECONDS)
    assert isinstance(found, dict)
    return found


def carried(selection: Any, name: str, came_from: str) -> str:
    """An ancestor's name in the chosen combination's tree, where the address names it
    in another's (`shown`): the same ancestor by clade, as the ancestor viewer does."""
    key = selection.combination["key"]
    if not name or came_from == key or f"{came_from}/ancestors.fasta" not in selection.keys:
        return name
    placed = analog_views.trees(selection)
    if name not in placed.get(came_from, {}):
        return name  # an extant sequence: the same name everywhere
    same = [
        a.node
        for a in analogs.find(placed, came_from, name, analog_views.by_clade(selection))
        if a.combination == key
    ]
    return same[0] if same else ""


def node_pp(selection: Any, node: str, row: str) -> list[float] | None:
    """The posterior probability of `node`'s state at every column, where the gap version
    keeps IQ-TREE's residues (one ranged read of the posteriors table)."""
    if selection.version not in (ancestors.ORIGINAL, *gaps.MASKING):
        return None
    name = "ancestors_posteriors.tsv"
    folder = ancestors.version_folder(selection.version)
    if f"{selection.combination['key']}/{folder}{name}" not in selection.keys:
        return None
    dists = ancestors.distributions(selection.read_node(name, node), node)
    if len(dists) != len(row):
        return None
    return [dict(dist).get(state.upper(), 0.0) for dist, state in zip(dists, row, strict=True)]


@dataclass(frozen=True)
class Shown:
    """One of the page's two sequences."""

    which: str  # A or B
    name: str
    label: str
    row: str
    pp: list[float] | None
    extant: bool

    @property
    def residues(self) -> str:
        """The sequence itself, gaps removed (what a structure predictor is given)."""
        return self.row.replace(mapping.GAP, "").upper()


@dataclass
class Loaded:
    """An attached structure read for this request: its chains and their maps."""

    structure: Structure
    chains: pdb.Parsed | None = None
    maps: dict[str, mapping.ChainMap] = field(default_factory=dict)
    error: str = ""

    @property
    def id(self) -> str:
        return str(self.structure.id)

    def columns(self, chain: pdb.Chain, alignment: dict[str, str]) -> list[int | None]:
        found = self.maps.get(chain.id)
        if found is None or found.reference not in alignment:
            return [None] * len(chain.residues)
        return mapping.columns(found, alignment[found.reference])


def load(structure: Structure, aligner: str, alignment: dict[str, str]) -> Loaded:
    try:
        chains = parsed(structure)
    except Exception as exc:  # a file gone from storage, or one no longer readable
        return Loaded(structure, error=f"This structure's file could not be read ({exc}).")
    return Loaded(structure, chains, chain_maps(structure, chains, aligner, alignment))


@dataclass(frozen=True)
class Fit:
    """How well one of the page's sequences sits on a structure's best chain for it."""

    chain: str
    identity: float  # of the chain's residues that the sequence has a residue for
    coverage: float  # share of the sequence's residues that have a place on the chain
    compared: int

    @property
    def score(self) -> float:
        return self.identity * self.coverage

    @property
    def is_model(self) -> bool:
        """The structure is this very sequence: a model of the ancestor (or, for an extant
        sequence, its own structure)."""
        return self.identity >= 0.995 and self.coverage >= 0.95


def fit_of(loaded: Loaded, alignment: dict[str, str], shown: Shown) -> Fit | None:
    if loaded.chains is None:
        return None
    length = len(shown.residues)
    best: Fit | None = None
    for chain in loaded.chains.chains:
        if chain.id not in loaded.maps:
            continue
        same = compared = 0
        for residue, column in zip(chain.residues, loaded.columns(chain, alignment), strict=True):
            if column is None or column >= len(shown.row) or shown.row[column] == mapping.GAP:
                continue
            compared += 1
            same += shown.row[column].upper() == residue.aa
        if not compared or not length:
            continue
        found = Fit(chain.id, same / compared, compared / length, compared)
        if best is None or found.score > best.score:
            best = found
    return best


@dataclass
class Panel:
    """A structure painted with one of the page's two sequences."""

    loaded: Loaded
    which: str = A
    chain: str = ""  # the chain asked for ("": the one that fits the sequence best)


def closest(loaded: list[Loaded], fits: dict[str, dict[str, Fit | None]], which: str) -> Loaded:
    """The structure the sequence sits on best (the first, if it sits on none)."""

    def score(one: Loaded) -> float:
        found = fits[one.id].get(which)
        return found.score if found else -1.0

    return max(loaded, key=score)


def choose_panels(
    query: Any, loaded: list[Loaded], fits: dict[str, dict[str, Fit | None]], has_b: bool
) -> list[Panel]:
    """The page's panels from its address: `ps`, `pw` and `pc` name each panel's
    structure, sequence and chain, in order. The one-off buttons are applied here:
    `preset` (a ready-made arrangement), `add` and `drop`."""
    by_id = {one.id: one for one in loaded}
    structures, which, chains = (query.getlist(key) for key in PANEL_KEYS)
    panels = []
    for i, structure_id in enumerate(structures[:MAX_PANELS]):
        if structure_id in by_id:
            sequence = which[i] if i < len(which) and which[i] in (A, B) and has_b else A
            panels.append(
                Panel(by_id[structure_id], sequence, chains[i] if i < len(chains) else "")
            )
    if not panels:
        asked = by_id.get(query.get("structure", ""))  # a link that names one structure
        panels = [Panel(asked or closest(loaded, fits, A))]

    preset = query.get("preset", "")
    first = panels[0]
    if preset == "one":
        panels = [Panel(first.loaded, A)]
    elif preset == "shared" and has_b:
        panels = [Panel(first.loaded, A), Panel(first.loaded, B)]
    elif preset == "closest" and has_b:
        panels = [Panel(closest(loaded, fits, A), A), Panel(closest(loaded, fits, B), B)]
    elif preset == "every":
        panels = [Panel(one, first.which) for one in loaded[:MAX_PANELS]]
    elif query.get("add") and len(panels) < MAX_PANELS:
        shown = {p.loaded.id for p in panels}
        unshown = next((one for one in loaded if one.id not in shown), None)
        if has_b and all(p.which == A for p in panels):
            panels.append(Panel(panels[-1].loaded, B))  # the other end, on the same structure
        else:
            panels.append(Panel(unshown or panels[-1].loaded, panels[-1].which))
    drop = query.get("drop", "")
    if drop.isdigit() and int(drop) < len(panels) and len(panels) > 1:
        del panels[int(drop)]
    return panels


def tip(site: mapping.Site, painted: Shown, other: Shown | None, note: str = "") -> str:
    """The text shown when a residue is pointed at."""
    text = f"{site.residue.name} {site.residue.label} · chain {site.chain}"
    if site.column is None:
        return f"{text} · not in the alignment"
    text += f" · column {site.column + 1}"

    def one(shown: Shown, state: str, pp: float | None) -> str:
        if state == mapping.GAP:
            return f"{shown.label}: no residue"
        return f"{shown.label}: {state}" + (f" (probability {pp:.2f})" if pp is not None else "")

    text += " · " + one(painted, site.state, site.pp)
    if other is not None:
        text += " · " + one(other, site.other, site.other_pp)
    return text + (f" · {note}" if note else "")


def change_label(a: str, column: int, b: str) -> str:
    """As the mutational path page writes a change: K300R, K300-, -300R."""
    return f"{a}{column + 1}{b}"


def path_steps(
    selection: Any, alignment: dict[str, str], first: Shown, second: Shown | None, called: Any
) -> tuple[colours.Steps | None, dict[int, list[str]]]:
    """The route between the page's two sequences, where the second is below the first
    on the tree and more than one branch away: for every column that changes along it,
    the branch it last changed on (`colours.Steps`) and each change in words."""
    if second is None or first.extant:
        return None, {}
    route = paths.path(selection.tree(), first.name, second.name)
    if len(route) < 3:
        return None, {}
    rows = []
    for name in route:
        row = selection.ancestors.get(name) or alignment.get(name)
        if row is None:
            return None, {}
        rows.append(row.upper())
    labels = [ancestor_names.label(name, called) for name in route]
    of_column: dict[int, int] = {}
    notes: dict[int, list[str]] = {}
    branches = []
    for step, ((parent, child), (above, below)) in enumerate(
        zip(pairwise(rows), pairwise(labels), strict=True)
    ):
        branches.append(f"{above} → {below}")
        for column, (a, b) in enumerate(zip(parent, child, strict=True)):
            if a != b:
                of_column[column] = step
                notes.setdefault(column, []).append(
                    f"{change_label(a, column, b)} on {above} → {below}"
                )
    return colours.Steps(branches, of_column), notes


@dataclass(frozen=True)
class ChangeRow:
    label: str
    column: int  # from 1
    kind: str  # sub | lost | gained
    where: list[str]  # per panel: "A LYS 295", or "" where its structure has no residue
    pp: float | None
    other_pp: float | None
    when: str = ""  # the branch(es) of the path it changed on

    @property
    def on_structure(self) -> bool:
        return any(self.where)


def view_context(job: Job, selection: Any, query: Any, structures: list[Structure]) -> Any:
    """Everything the viewer form shows: the choices, and the panels painted."""
    combo, key = selection.combination, selection.combination["key"]
    aligners, models = views.combination_menus(selection)
    support, _ = views.node_support(selection)
    called = ancestor_names.for_selection(selection)
    nodes = ancestors.node_choices(selection.tree(), selection.node_names, support, called)
    alignment = ancestors.parse_fasta(selection.read("alignment.fasta"))
    came_from = query.get("shown", "")
    node = carried(selection, query.get("node", ""), came_from)
    if node not in selection.ancestors:
        node = next((n.name for n in nodes if n.name in selection.ancestors), "")
    if not node:
        return {"no_ancestors": True}
    to = carried(selection, query.get("to", ""), came_from)
    if to == node or (to not in selection.ancestors and to not in alignment):
        to = ""

    def shown(which: str, name: str) -> Shown:
        if name in selection.ancestors:
            row = selection.ancestors[name]
            return Shown(which, name, ancestor_names.label(name, called), row,
                         node_pp(selection, name, row), False)  # fmt: skip
        return Shown(which, name, name, alignment[name], None, True)

    first = shown(A, node)
    second = shown(B, to) if to else None
    palette = _choice(query, "palette", ancestors.PALETTES, ancestors.DEFAULT_PALETTE)
    view: dict[str, Any] = {
        "aligner": combo["aligner"], "combination": key, "gaps": selection.version,
        "node": node,
    }  # fmt: skip
    if to:
        view["to"] = to
    context: dict[str, Any] = {
        "aligners": aligners, "aligner": combo["aligner"], "models": models,
        "combination": key, "nodes": nodes, "node": node, "to": to,
        "extant_names": list(alignment), "first": first, "second": second,
        "versions": [(v, ancestors.version_label(v)) for v in selection.versions],
        "version": selection.version,
        "version_label": ancestors.version_label(selection.version),
        "palettes": [(k, label) for k, (label, _) in ancestors.PALETTES.items()],
        "palette": palette, "structures": structures,
    }  # fmt: skip
    if structures:
        loaded = [load(s, combo["aligner"], alignment) for s in structures]
        steps = path_steps(selection, alignment, first, second, called)
        context |= board(job, query, loaded, alignment, first, second, palette, steps)
        view |= context["view"]
    context["canonical_url"] = (
        reverse("structure_viewer", args=[job.id]) + "?" + urlencode(view, doseq=True)
    )
    return context


def board(
    job: Job,
    query: Any,
    loaded: list[Loaded],
    alignment: dict[str, str],
    first: Shown,
    second: Shown | None,
    palette: str,
    path: tuple[colours.Steps | None, dict[int, list[str]]],
) -> dict[str, Any]:
    """The page's panels painted: for each, its structure's chains tied to the alignment,
    every residue's colour in the chosen scheme, its sequence strip; and for the page the
    superposition that lets the panels share a camera, and the table of changes."""
    sequences = {A: first} | ({B: second} if second else {})
    fits = {one.id: {w: fit_of(one, alignment, s) for w, s in sequences.items()} for one in loaded}
    panels = choose_panels(query, loaded, fits, second is not None)
    steps, path_notes = path
    has_pp = any(sequences[p.which].pp is not None for p in panels)
    schemes = colours.available(second is not None, has_pp, steps is not None)
    default = "changes" if second is not None else "pp" if has_pp else "differs"
    scheme = _choice(query, "scheme", schemes, default)
    link = _on(query, "link", True)
    display = {
        "scheme": scheme,
        "palette": palette,
        "style": _choice(query, "style", STYLES, "cartoon"),
        "sidechains": _choice(query, "sidechains", SIDECHAINS, "marked"),
        "surface": "on" if _on(query, "surface", False) else "off",
        "ligands": "on" if _on(query, "ligands", True) else "off",
        "labels": "on" if _on(query, "labels", False) else "off",
        "spin": "on" if _on(query, "spin", False) else "off",
        "background": _choice(query, "background", BACKGROUNDS, "auto"),
        "link": "on" if link else "off",
        "scroll": "on" if _on(query, "scroll", True) else "off",
    }

    shown_panels: list[dict[str, Any]] = []
    data_panels: list[dict[str, Any]] = []
    focused: list[list[mapping.Site]] = []  # each panel's focus chain
    for index, panel in enumerate(panels):
        one = panel.loaded
        painted = sequences[panel.which]
        other = next((s for w, s in sequences.items() if w != panel.which), None)
        found_fit = fits[one.id].get(panel.which)
        entry: dict[str, Any] = {
            "index": index, "number": index + 1, "structure": one.structure,
            "which": panel.which, "painted": painted, "fit": found_fit, "error": one.error,
            "chains": list(one.maps), "focus": "", "strip": [],
        }  # fmt: skip
        shown_panels.append(entry)
        if one.chains is None:
            focused.append([])
            data_panels.append({"file": "", "keys": [], "error": True})
            continue
        focus = panel.chain if panel.chain in one.maps else ""
        focus = focus or (found_fit.chain if found_fit else next(iter(one.maps), ""))
        reverse_ = panel.which == B
        sites: list[mapping.Site] = []
        mine: list[mapping.Site] = []
        for chain in one.chains.chains:
            found = mapping.sites(
                chain, one.columns(chain, alignment), painted.row, painted.pp,
                other.row if other else "", other.pp if other else None,
            )  # fmt: skip
            sites += found
            if chain.id == focus:
                mine = found
        focused.append(mine)
        paint = colours.paint(scheme, sites, palette, one.structure.predicted, reverse_, steps)
        if scheme == "changes":
            marked = [i for i, s in enumerate(sites) if s.change]
        elif scheme == "differs":
            marked = [i for i, s in enumerate(sites) if s.differs]
        elif scheme == "steps" and steps is not None:
            marked = [i for i, s in enumerate(sites) if s.column in steps.of_column]
        else:
            marked = []
        labels = {}
        for i in marked[:MAX_LABELS]:
            s = sites[i]
            if s.column is None:
                continue
            if scheme == "differs":
                labels[i] = change_label(s.residue.aa, s.column, s.state)
            else:  # always named from the first sequence to the second
                a, b = (s.other, s.state) if reverse_ else (s.state, s.other)
                labels[i] = change_label(a, s.column, b)
        between = mapping.between(mine, painted.row, other.row if other else "")
        notes = {
            s.key: "; ".join(path_notes.get(s.column, [])) for s in sites if s.column is not None
        }
        colour_of = dict(zip((s.key for s in sites), paint.colours, strict=True))
        after = {b.after: b for b in between}
        for site in mine:
            colour = colour_of[site.key]
            style = f"background: {colour}; color: {colours.text_colour(colour)}"
            states = {panel.which: site.state, (other.which if other else ""): site.other}
            entry["strip"].append(
                {
                    "site": site,
                    "tip": tip(site, painted, other, notes.get(site.key, "")),
                    "tick": site.residue.number % 10 == 0 and not site.residue.icode,
                    "between": after.get(site.key),
                    "a": states.get(A, ""),
                    "b": states.get(B, ""),
                    # the painted sequence's cell has the colour; the other's where they differ
                    "a_style": style if panel.which == A or site.change else "",
                    "b_style": style if panel.which == B or site.change else "",
                }
            )
        entry |= {
            "focus": focus,
            "focus_map": one.maps.get(focus),
            "residues": len(mine),
            "legend": paint.legend,
            "legend_note": paint.note,
            "counts": {
                "sub": sum(s.differs == "sub" for s in mine),
                "absent": sum(s.differs == "absent" for s in mine),
                "extra": sum(b.count for b in between),
            },
        }
        data_panels.append(
            {
                "file": reverse("structure_file", args=[job.id, one.structure.id]),
                "format": one.chains.format,
                "title": f"{painted.label} on {one.structure.label}",
                "keys": [s.key for s in sites],
                "cols": [s.column for s in sites],
                "focus": focus,
                "colours": paint.colours,
                "tips": [tip(s, painted, other, notes.get(s.key, "")) for s in sites],
                "marked": marked,
                "labels_text": {str(i): text for i, text in labels.items()},
                "between": [
                    {"after": b.after, "before": b.before, "tip": between_tip(b, painted, other)}
                    for b in between
                ],
                "matrix": None,
                "linked": link,
            }
        )

    # One frame for every panel: each other structure is laid on the first panel's, by the
    # C-alphas of the residues the two have in the same alignment columns.
    frame = {s.column: s.residue for s in focused[0] if s.column is not None}
    for index, panel in enumerate(panels):
        entry, data = shown_panels[index], data_panels[index]
        entry["superposed"] = ""
        if not link or len(panels) == 1 or data.get("error"):
            continue
        if index == 0:
            entry["superposed"] = "frame"
        elif panel.loaded.id == panels[0].loaded.id:
            entry["superposed"] = "same"
        else:
            pairs = [(s.residue, frame[s.column]) for s in focused[index] if s.column in frame]
            fitted = superpose.fit(
                [(r.x, r.y, r.z) for r, _ in pairs], [(r.x, r.y, r.z) for _, r in pairs]
            )
            if fitted is None:
                entry["superposed"] = "free"
                data["linked"] = False
            else:
                entry["superposed"] = "fitted"
                entry["superposition"] = fitted
                data["matrix"] = fitted.matrix

    changes = change_rows(first, second, focused, path_notes) if second else []
    view: dict[str, Any] = {
        "ps": [p.loaded.id for p in panels],
        "pw": [p.which for p in panels],
        "pc": [entry["focus"] for entry in shown_panels],
    } | display
    fit_rows = [
        {"structure": one.structure, "a": fits[one.id].get(A), "b": fits[one.id].get(B)}
        for one in loaded
    ]
    offered = ["one"]
    if second is not None:
        offered.append("shared")
        if len(loaded) > 1:
            offered.append("closest")
    if len(loaded) > 1:
        offered.append("every")
    return {
        "view": view,
        "display": display,
        "panels": shown_panels,
        "several": len(panels) > 1,
        "can_add": len(panels) < MAX_PANELS,
        "presets": [(k, PRESETS[k]) for k in offered],
        "structure_data": {"panels": data_panels, "betweenColour": colours.INSERTION, **display},
        "legend": shown_panels[0].get("legend", []),
        "legend_note": shown_panels[0].get("legend_note", ""),
        "legend_per_panel": scheme == "confidence" and len(panels) > 1,
        "schemes": [(k, colours.SCHEMES[k]) for k in schemes],
        "scheme_help": [(colours.SCHEMES[k], colours.ABOUT[k]) for k in schemes],
        "styles": STYLES,
        "sidechain_choices": SIDECHAINS,
        "backgrounds": BACKGROUNDS,
        "fit_rows": fit_rows,
        "steps": steps,
        "changes": changes[:MAX_CHANGES],
        "changes_total": len(changes),
        "changes_on_structure": sum(c.on_structure for c in changes),
        "insertion_colour": colours.INSERTION,
    }


def change_rows(
    first: Shown,
    second: Shown,
    focused: list[list[mapping.Site]],
    path_notes: dict[int, list[str]],
) -> list[ChangeRow]:
    """Every column where the two sequences differ, and the residue of each panel's
    structure (its focus chain) that the column falls on."""
    at = [{s.column: s for s in sites if s.column is not None} for sites in focused]
    rows = []
    for column, (a, b) in enumerate(zip(first.row.upper(), second.row.upper(), strict=True)):
        if a == b:
            continue
        kind = "lost" if b == mapping.GAP else "gained" if a == mapping.GAP else "sub"
        where = []
        for sites in at:
            s = sites.get(column)
            where.append(f"{s.chain} {s.residue.name} {s.residue.label}" if s else "")
        when = "; ".join(n.split(" on ", 1)[-1] for n in path_notes.get(column, []))
        rows.append(
            ChangeRow(
                change_label(a, column, b),
                column + 1,
                kind,
                where,
                first.pp[column] if first.pp and a != mapping.GAP else None,
                second.pp[column] if second.pp and b != mapping.GAP else None,
                when,
            )
        )
    return rows


def between_tip(b: mapping.Between, painted: Shown, other: Shown | None) -> str:
    text = f"Columns {b.first + 1} to {b.last + 1}: {painted.label} has {b.count} residue"
    text += "" if b.count == 1 else "s"
    if other is not None:
        text += f", {other.label} has {b.other}"
    return text + " the structure has no residue for"


def manage_context(job: Job, selection: Any, structures: list[Structure]) -> dict[str, Any]:
    """The Attached Structures section: what is attached, and the ways to attach more."""
    names: list[str] = []
    if selection is not None:
        names = list(ancestors.parse_fasta(selection.read("alignment.fasta")))
    have = {(s.source, s.accession) for s in structures}
    return {
        "structures": structures,
        "suggested": [
            (name, accession)
            for name, accession in remote.accessions_in(names).items()
            if (remote.ALPHAFOLD, accession) not in have
        ][:200],
        "lookup_form": SequenceLookupForm(initial={"reviewed_only": False}),
        "model_species": lookup.MODEL_SPECIES,
        "sources": remote.SOURCES,
        "max_structures": MAX_PER_JOB,
        "max_at_once": MAX_AT_ONCE,
        "can_attach": len(structures) < MAX_PER_JOB,
        "max_mb": pdb.MAX_BYTES // (1024 * 1024),
    }


@require_GET
def page(request: HttpRequest, job_id: str) -> HttpResponse:
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    selection = views.select_ancestor(job, request.GET)
    structures = attached(job)
    context: dict[str, Any] = {"job": job, "selection": selection}
    if selection is not None:
        context |= view_context(job, selection, request.GET, structures)
    htmx = bool(request.headers.get("HX-Request"))
    # The preset, add and drop buttons are one-off actions: the address gets the
    # resulting view, never the button, so reloading does not apply it again.
    one_off = any(key in request.GET for key in ("preset", "add", "drop"))
    if one_off and not htmx and "canonical_url" in context:
        return redirect(context["canonical_url"])
    if htmx:
        response = render(request, "phylobot/_structure_viewer.html", context)
        if "canonical_url" in context:
            response["HX-Push-Url"] = context["canonical_url"]
        return response
    context |= manage_context(job, selection, structures)
    return render(request, "phylobot/structure.html", context)


@require_GET
def find(request: HttpRequest, job_id: str) -> HttpResponse:
    """Structures of a protein found by name: UniProt's PDB entries and AlphaFold model
    for each match (`structures.remote.search`), each with a box to tick. Only the name
    and species leave the machine, cleaned, after the consent tick."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    form = SequenceLookupForm(request.GET or None)
    found: list[remote.Protein] | None = None
    error = ""
    if form.is_valid():
        try:
            found = remote.search(form.cleaned_data["protein"], form.cleaned_data["species"])
            if form.cleaned_data["reviewed_only"]:
                found = [p for p in found if p.reviewed]
        except SearchError as exc:
            error = str(exc)
    elif request.GET:
        error = str(next((e for errors in form.errors.values() for e in errors), "Fill in a name."))
    have = {(s.source, s.accession) for s in attached(job)}
    context = {"job": job, "found": found, "find_error": error, "max_at_once": MAX_AT_ONCE,
               "attached_pdb": {a for s, a in have if s == remote.RCSB},
               "attached_models": {a for s, a in have if s == remote.ALPHAFOLD}}  # fmt: skip
    return render(request, "phylobot/_structure_lookup_results.html", context)


def back_to(request: HttpRequest, job: Job, **changed: Any) -> HttpResponse:
    """Back to the page the POST came from, as its address described it: only the
    page's own keys are carried, re-encoded, so the redirect can lead nowhere else."""
    asked = parse_qs(urlsplit(request.headers.get("Referer", "")).query)
    view: dict[str, Any] = {}
    for key, values in asked.items():
        if key in PANEL_KEYS:
            view[key] = [v[:300] for v in values[:MAX_PANELS]]
        elif key in VIEW_KEYS:
            view[key] = values[-1][:300]
    view |= changed
    view = {k: v for k, v in view.items() if v}
    url = reverse("structure_viewer", args=[job.id])
    return redirect(f"{url}?{urlencode(view, doseq=True)}" if view else url)


def showing(structure: Structure) -> dict[str, Any]:
    """The page's panels replaced by one, of this structure."""
    return {"ps": [str(structure.id)], "pw": [A], "pc": []}


def store(
    job: Job,
    source: str,
    accession: str,
    label: str,
    title: str,
    data: bytes,
    predicted: bool = False,
) -> Structure:
    """Check a structure file and attach it; StructureError if it is not one."""
    chains = pdb.parse(data)
    structure = Structure(
        id=uuid.uuid4(),
        job=ancestor_names.holder(job),
        label=label[:120],
        source=source,
        accession=accession,
        title=(title or chains.title)[:300],
        chains=[{"id": c.id, "residues": len(c.residues)} for c in chains.chains],
        size=len(data),
        format=chains.format,
        predicted=predicted or source == remote.ALPHAFOLD or pdb.looks_predicted(data),
    )
    storage.put_bytes(structure.storage_key, data)
    structure.save()
    return structure


def report(
    request: HttpRequest, done: list[Structure], already: list[str], failed: list[str]
) -> None:
    if done:
        messages.success(request, "Attached " + ", ".join(s.label for s in done) + ".")
    if already:
        messages.info(request, "Already attached: " + ", ".join(already) + ".")
    for text in failed:
        messages.error(request, f"Not attached: {text}.")


@require_POST
def attach(request: HttpRequest, job_id: str) -> HttpResponse:
    """Attach PDB entries and AlphaFold models, fetched by their identifiers: those
    ticked in the lookup's results (`pick`, "rcsb:2SRC"), or those typed (`identifier`)."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if request.POST.get("consent") != "on":
        messages.error(request, "Tick the box to fetch a structure: its identifier is sent to "
                       "the database it comes from.")  # fmt: skip
        return back_to(request, job)
    shapes = {remote.RCSB: remote.PDB_ID, remote.ALPHAFOLD: remote.ACCESSION}
    try:
        picks = request.POST.getlist("pick")
        if picks:
            wanted: list[tuple[str, str]] = []
            for pick in picks[: MAX_AT_ONCE + 1]:
                source, _, accession = pick.partition(":")
                shape = shapes.get(source)
                if shape is None or not shape.fullmatch(accession):
                    raise SearchError("not an identifier of that database")
                if (source, accession.upper()) not in wanted:
                    wanted.append((source, accession.upper()))
            if len(wanted) > MAX_AT_ONCE:
                raise SearchError(f"at most {MAX_AT_ONCE} structures at a time")
        else:
            wanted = remote.clean_ids(request.POST.get("identifier", ""), MAX_AT_ONCE)
    except SearchError as exc:
        messages.error(request, f"Nothing was attached: {exc}.")
        return back_to(request, job)
    have = {(s.source, s.accession): s for s in attached(job)}
    done: list[Structure] = []
    already: list[Structure] = []
    failed: list[str] = []
    for source, accession in wanted:
        label = accession if source == remote.RCSB else f"AF-{accession}"
        if (source, accession) in have:
            already.append(have[(source, accession)])
        elif len(have) + len(done) >= MAX_PER_JOB:
            failed.append(f"{label} (a job holds at most {MAX_PER_JOB} structures)")
        else:
            try:
                fetched = remote.fetch(source, accession)
                done.append(store(job, source, accession, label, fetched.title, fetched.data))
            except (SearchError, pdb.StructureError) as exc:
                failed.append(f"{label} ({exc})")
    report(request, done, [s.label for s in already], failed)
    one = (done or already)[0] if len(done) + len(already) == 1 else None
    return back_to(request, job, **(showing(one) if one else {}))


def from_zip(data: bytes) -> list[tuple[str, bytes]]:
    """The structure files in a ZIP, as (name, content): what a structure predictor's
    download is. A predictor writes several models; where the archive names its best
    (AlphaFold Server's `..._model_0.cif`, AlphaFold 3's `<job>_model.cif` beside its
    per-sample folders) only that one is taken. Hostile input: nothing is extracted to
    disk, names are only read, and no member is read beyond the size a structure may be."""
    try:
        archive = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile as exc:
        raise pdb.StructureError("the ZIP could not be read") from exc
    members = [
        info
        for info in archive.infolist()[:MAX_ZIP_MEMBERS]
        if not info.is_dir()
        and info.filename.lower().endswith(STRUCTURE_FILES)
        and info.file_size <= pdb.MAX_BYTES
    ]
    best = [m for m in members if re.search(r"_model_0\.cif$", m.filename, re.IGNORECASE)] or [
        m
        for m in members
        if m.filename.lower().endswith("_model.cif") and "seed-" not in m.filename.lower()
    ]
    found = []
    for info in (best or members)[:MAX_AT_ONCE]:
        try:
            with archive.open(info) as member:
                content = member.read(pdb.MAX_BYTES + 1)
        except zipfile.BadZipFile, RuntimeError, NotImplementedError, OSError:
            continue  # encrypted, or damaged
        if len(content) <= pdb.MAX_BYTES:
            found.append((PurePosixPath(info.filename).name, content))
    if not found:
        raise pdb.StructureError("the ZIP holds no PDB or mmCIF file")
    return found


def clean_name(name: str) -> str:
    kept = "".join(c for c in name if c.isalnum() or c in " ._-()")
    return " ".join(kept.split())[:80]


@require_POST
def upload(request: HttpRequest, job_id: str) -> HttpResponse:
    """Attach the user's own files, PDB format or mmCIF, or a ZIP of them: structures not
    in a database, or models of this job's ancestors predicted elsewhere. Hostile input:
    size-capped before it is read, and stored only if `structures.pdb` finds a protein
    chain in it."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    sent = request.FILES.getlist("file")
    if not sent:
        messages.error(request, "Choose a PDB or mmCIF file (or a ZIP of them) to upload.")
        return back_to(request, job)
    files: list[tuple[str, bytes]] = []
    failed: list[str] = []
    for one in sent[:MAX_AT_ONCE]:
        name = clean_name(one.name or "") or "Uploaded"
        zipped = name.lower().endswith(".zip")
        limit = MAX_ZIP_BYTES if zipped else pdb.MAX_BYTES
        try:
            if one.size is None or one.size > limit:
                raise pdb.StructureError(f"the file is over {limit // (1024 * 1024)} MB")
            data = one.read(limit + 1)
            files += from_zip(data) if zipped else [(name, data)]
        except pdb.StructureError as exc:
            failed.append(f"{name} ({exc})")
    if len(sent) > MAX_AT_ONCE or len(files) > MAX_AT_ONCE:
        failed.append(f"more than {MAX_AT_ONCE} files at a time (the rest were left)")
    typed = " ".join(request.POST.get("label", "").split())[:80]
    predicted = request.POST.get("predicted") == "on"
    count = len(attached(job))
    done: list[Structure] = []
    for name, data in files[:MAX_AT_ONCE]:
        label = typed if typed and len(files) == 1 else clean_name(PurePosixPath(name).stem) or name
        if count + len(done) >= MAX_PER_JOB:
            failed.append(f"{label} (a job holds at most {MAX_PER_JOB} structures)")
            continue
        try:
            done.append(store(job, remote.UPLOAD, "", label, "", data, predicted))
        except pdb.StructureError as exc:
            failed.append(f"{label} ({exc})")
    report(request, done, [], failed)
    return back_to(request, job, **(showing(done[0]) if len(done) == 1 else {}))


@require_POST
def drop(request: HttpRequest, job_id: str, structure_id: str) -> HttpResponse:
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    structure = get_object_or_404(Structure, id=structure_id, job=ancestor_names.holder(job))
    storage.delete_key(structure.storage_key)
    label = structure.label
    structure.delete()
    messages.success(request, f"Removed {label}.")
    # The panels are chosen afresh: one of them may have shown the structure.
    return back_to(request, job, ps=[], pw=[], pc=[])


@require_GET
def file(request: HttpRequest, job_id: str, structure_id: str) -> HttpResponse:
    """The structure's file, for structure.js (and as a download): plain text, never
    anything a browser would render."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    structure = get_object_or_404(Structure, id=structure_id, job=ancestor_names.holder(job))
    try:
        data = storage.read_bytes(structure.storage_key)
    except Exception as exc:
        raise Http404("the structure's file is gone") from exc
    response = HttpResponse(data, content_type="text/plain; charset=us-ascii")
    response["X-Content-Type-Options"] = "nosniff"
    response["Cache-Control"] = "private, max-age=3600"
    if request.GET.get("download"):
        name = structure.accession or str(structure.id)
        response["Content-Disposition"] = f'attachment; filename="{name}.{structure.format}"'
    return response
