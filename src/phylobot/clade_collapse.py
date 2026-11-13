"""Clades the user folds into triangles in the tree viewer, to save vertical space in a
publication figure (`phylobot.treeview` draws them; `tree.js` and `_tree_viewer.html` are
the UI).

A fold belongs to the ancestor's **clade**, the same content-derived ID an ancestor name
uses (`phylobot.analogs`, `ancestor_map.tsv`): the clade ID with an outgroup, else the
parts ID. So a clade folded on one combination is folded on every other combination that
has it, its custom label and triangle style come along, and the fold persists the next
time the tree viewer is opened. Like `phylobot.ancestor_names`, the state is held on the
topmost job an extension inherits from (`holder`), and the label and colours are user text
applied only where the web tier draws the tree, never in a pipeline file or storage path.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from django.db import IntegrityError, transaction

from phylobot import analog_views
from phylobot.ancestor_names import holder, match_key
from phylobot.models import CladeCollapse, Job

LABEL_MAX = 60
# A label is free text (clade names carry spaces, "/", parentheses), only kept printable
# and bounded; Django escapes it into the SVG. Control characters and angle brackets out.
_CONTROL = re.compile(r"[\x00-\x1f\x7f<>]")

DEFAULT_FILL = "#cfd8e6"
DEFAULT_BORDER = "#2457a6"
DEFAULT_BORDER_WIDTH = 1.6
# Offered in the popover's "border thickness" menu (value → label).
BORDER_WIDTHS = ((1.0, "Thin"), (1.6, "Medium"), (2.5, "Thick"), (4.0, "Heavy"))
# Quick-pick swatches in the popover (the colour input gives any value): muted, print-safe
# fills and darker borders to match. "none" is an open (unfilled) triangle.
FILL_SWATCHES = ("#cfd8e6", "#dfe7d6", "#f0e2cc", "#ead6e0", "#d6e6e6", "#e6e6e6", "none")
BORDER_SWATCHES = ("#2457a6", "#3c7a3c", "#a6682a", "#8e44ad", "#2a8a8a", "#1d2330")
_HEX = re.compile(r"#[0-9a-fA-F]{6}")


@dataclass(frozen=True)
class Style:
    """A folded clade's drawing, resolved for one combination's node."""

    label: str
    fill: str
    border_colour: str
    border_width: float
    alike: bool  # the same subtrees meet here as where it was folded (else approximate)


def clean_label(text: str) -> str:
    """A safe, bounded clade label ("" if blank)."""
    return _CONTROL.sub("", str(text)).strip()[:LABEL_MAX]


def clean_colour(text: str, fallback: str, *, allow_none: bool = False) -> str:
    """A `#rrggbb` colour, or `fallback`; "none" is kept when `allow_none`."""
    value = str(text).strip().lower()
    if allow_none and value == "none":
        return "none"
    return value if _HEX.fullmatch(value) else fallback


def clean_width(text: Any) -> float:
    try:
        value = float(text)
    except TypeError, ValueError:
        return DEFAULT_BORDER_WIDTH
    allowed = [w for w, _ in BORDER_WIDTHS]
    return min(allowed, key=lambda w: abs(w - value))


def save(
    job: Job, clade: str, parts: str, label: str, fill: str, border_colour: str, border_width: Any
) -> None:
    """Fold a clade (or restyle one already folded). Values are cleaned here."""  # fmt: skip
    owner = holder(job)
    defaults = {
        "parts": parts,
        "label": clean_label(label),
        "fill": clean_colour(fill, DEFAULT_FILL, allow_none=True),
        "border_colour": clean_colour(border_colour, DEFAULT_BORDER),
        "border_width": clean_width(border_width),
    }
    try:
        with transaction.atomic():
            CladeCollapse.objects.update_or_create(job=owner, clade=clade, defaults=defaults)
    except IntegrityError:  # folded meanwhile
        pass


def remove(job: Job, clade: str) -> None:
    """Expand the clade again."""
    CladeCollapse.objects.filter(job=holder(job), clade=clade).delete()


def cleaned_entry(entry: dict[str, Any]) -> dict[str, Any]:
    """Cleaned style fields for one imported collapse (phylobot.job_import, hostile input);
    its clade/parts structure is checked there."""
    return {
        "label": clean_label(entry.get("label", "")),
        "fill": clean_colour(entry.get("fill", ""), DEFAULT_FILL, allow_none=True),
        "border_colour": clean_colour(entry.get("border_colour", ""), DEFAULT_BORDER),
        "border_width": clean_width(entry.get("border_width", "")),
    }


def for_nodes(job: Job, nodes: dict[str, Any], by_clade: bool) -> dict[str, Style]:
    """{node: Style} for the folded clades among one combination's nodes
    ({node: analogs.Located}); {} when none are folded."""
    rows = {r.clade: r for r in CladeCollapse.objects.filter(job=holder(job))}
    if not rows:
        return {}
    out: dict[str, Style] = {}
    for node, place in nodes.items():
        row = rows.get(match_key(place, by_clade))
        if row is not None:
            out[node] = Style(
                row.label, row.fill, row.border_colour, row.border_width,
                row.parts == str(place.parts),
            )  # fmt: skip
    return out


def for_selection(selection: Any, combination: str = "") -> dict[str, Style]:
    """{node: Style} for a viewer selection's combination, from the job's ancestor map."""
    placed = analog_views.trees(selection)
    nodes = placed.get(combination or selection.combination["key"], {})
    return for_nodes(selection.job, nodes, analog_views.by_clade(selection)) if nodes else {}
