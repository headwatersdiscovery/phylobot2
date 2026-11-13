"""Measuring properties after a job has finished.

Answers the three questions the offer is made of:

- which registry properties are **missing** for an ancestor, and which are already there;
- what a **scope** expands to — this ancestor, this ancestor everywhere it was recovered,
  every ancestor here, everything;
- what the chosen properties over that scope would **cost**.

The point of the feature is not "run what you skipped" — with the pipeline default on,
a new job already has every tier-0 property. It is **run what did not exist yet**: the
registry grows, and no job run before a property existed could otherwise ever have it,
while re-running a whole pipeline for one descriptor is absurd when the ancestors have
not changed. So the offer is built from the *difference* between the registry and what
was stored.

Pure: callers pass what they have read. Matching the same ancestor across combinations is
`phylobot.analogs`, and staging and submission are the web tier's.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from phylobot import properties

#: The scopes offered, in increasing cost, with what each means to a reader.
SCOPES: dict[str, str] = {
    "node": "This ancestor, in this combination",
    "analogs": "This ancestor, in every combination that recovered it",
    "combination": "Every ancestor in this combination",
    "job": "Every ancestor in every combination",
}
DEFAULT_SCOPE = "analogs"


def scope_of(asked: object) -> str:
    """The chosen scope, or the default for anything unknown (raw request text)."""
    return asked if isinstance(asked, str) and asked in SCOPES else DEFAULT_SCOPE


def measured(summary: dict[str, Any] | None) -> set[str]:
    """The registry properties a stored `summary.json` already holds. Keys it carries
    that the registry no longer knows are ignored: they are another version's, and
    nothing can be offered about them."""
    if not summary:
        return set()
    stored = summary.get("properties")
    if not isinstance(stored, dict):
        return set()
    return {k for k in stored if k in properties.REGISTRY}


@dataclass(frozen=True)
class Offer:
    """One property on the checklist."""

    key: str
    label: str
    units: str
    note: str
    tier: int
    #: Already measured for this ancestor: shown ticked and disabled, since it is the
    #: answer to "what do I already have?" as much as an offer.
    done: bool
    #: Reference CPU-seconds for the chosen scope and draws, or 0 when already done.
    cost: float
    #: Reference CPU-seconds for *one* ancestor at 1000 draws. The measure dialog scales
    #: this in the browser as the reader changes the draws or the scope, rather than
    #: re-requesting the viewer and closing itself mid-decision (`property_measure.js`).
    unit_cost: float = 0.0
    #: For a row already measured: the draws it was measured on, so the one list can
    #: offer "recompute" beside "measure" and say what each would replace. `cost` is 0
    #: for such a row (nothing to spend to have what you have) while `unit_cost` is not,
    #: because recomputing it costs the same per ancestor as measuring it did.
    measured_draws: int = 0

    @property
    def offered(self) -> bool:
        return not self.done

    @property
    def cost_text(self) -> str:
        """The estimate as the list shows it: hours and minutes, not CPU-seconds."""
        return cost_text(self.cost)


def cost_of(key: str, targets: int, draws: int, length: int) -> float:
    """Reference CPU-seconds to measure one property over `targets` ancestors.

    Linear in the draws and in the sequence length, from the per-property figures in the
    registry (1000 draws of a 400-residue ancestor). An estimate: it orders the
    checklist and warns before a long run; it is not yet a measured cost model.
    """
    prop = properties.REGISTRY.get(key)
    if prop is None or targets <= 0:
        return 0.0
    return prop.cost * targets * (draws / 1000) * (max(length, 1) / 400)


def offers(
    summary: dict[str, Any] | None,
    targets: int = 1,
    draws: int = 1000,
    length: int = 400,
) -> list[Offer]:
    """The checklist: every registry property, cheapest tier first, with what it would
    cost over `targets` ancestors and whether it is already there."""
    have = measured(summary)  # read once: `offers` is called per request
    stored = (summary or {}).get("properties")
    stored = stored if isinstance(stored, dict) else {}

    def already(key: str) -> int:
        entry = stored.get(key)
        try:
            return int((entry or {}).get("draws") or 0)
        except TypeError, ValueError:
            return 0

    out = [
        Offer(
            key=prop.key,
            label=prop.label,
            units=prop.units,
            note=prop.note,
            tier=prop.tier,
            done=prop.key in have,
            cost=0.0 if prop.key in have else cost_of(prop.key, targets, draws, length),
            unit_cost=cost_of(prop.key, 1, 1000, length),
            measured_draws=already(prop.key) if prop.key in have else 0,
        )
        for prop in properties.REGISTRY.values()
    ]
    # Cheapest tier first, then cheapest within it; a stable key breaks ties so the list
    # does not reorder itself between requests.
    return sorted(out, key=lambda o: (o.tier, o.cost, o.key))


#: The checklist's field name in the viewer's GET form. A hidden field of this name is
#: submitted before the boxes, so an all-unticked checklist is distinguishable from a
#: first load (an unticked checkbox submits nothing).
FIELD = "measure"
SCOPE_FIELD = "measure_scope"


def chosen(asked: list[str] | None, offered: list[str]) -> list[str]:
    """The ticked properties, in checklist order.

    `asked` is the submitted values, or None when the field is absent — a first load,
    where every missing *default* property is ticked. With the field present, a key that
    is not in it is unticked, including all of them: that is what the hidden sentinel
    buys.
    """
    if asked is None:
        # Only the cheap defaults: a predictor that costs minutes per ancestor is offered
        # but never pre-ticked (owner, 2026-10-21), so the button never queues an hour by
        # accident.
        return [key for key in offered if key in properties.DEFAULTS]
    wanted = set(asked)
    return [key for key in offered if key in wanted]


def cost_text(seconds: float) -> str:
    """A cost as a reader thinks of one: hours, minutes and seconds (owner, 2026-10-21).

    Raw CPU-seconds were unreadable once the predictors arrived — "7986 CPU-s" is a real
    figure the viewer showed. Two units at most, because these are estimates from
    per-property reference figures and a third would claim more than we know. The caller
    says the word "estimated"; this returns only the quantity.
    """
    if seconds <= 0:
        return "—"
    if seconds < 1:
        return "under a second"
    hours, rest = divmod(int(round(seconds)), 3600)
    minutes, secs = divmod(rest, 60)
    if hours:
        return f"{hours} h {minutes} min" if minutes else f"{hours} h"
    if minutes:
        return f"{minutes} min {secs} s" if secs else f"{minutes} min"
    return f"{secs} s"


def missing(summary: dict[str, Any] | None) -> list[str]:
    """The registry properties this ancestor does not have, in checklist order: what the
    panel offers (`chosen` ticks the defaults among them), and empty when there is
    nothing to offer."""
    return [o.key for o in offers(summary) if o.offered]


@dataclass(frozen=True)
class Target:
    """One ancestor to measure: a combination and a node of its tree."""

    combination: str
    node: str


def expand(
    scope: str,
    combination: str,
    node: str,
    analogs: dict[str, str],
    nodes: dict[str, list[str]],
) -> list[Target]:
    """The ancestors a scope covers.

    `analogs` is {combination: its name for this node}, already matched by
    `phylobot.analogs` on the web tier, where the trees are — the pipeline is handed an
    explicit list and never repeats clade matching.

    `nodes` is {combination: its own ancestral node names}. It has to be per combination:
    node names belong to a tree, so `Node3` of one combination is not `Node3` of another,
    and reusing one list across all of them would name ancestors that do not exist.
    """
    if scope == "node":
        return [Target(combination, node)]
    if scope == "analogs":
        here = [Target(combination, node)]
        return here + [Target(c, n) for c, n in sorted(analogs.items()) if c != combination]
    if scope == "combination":
        return [Target(combination, n) for n in nodes.get(combination, [])]
    return [Target(c, n) for c in sorted(nodes) for n in nodes[c]]


def total_cost(keys: list[str], targets: int, draws: int, length: int) -> float:
    """What the whole run would cost, for the figure shown before the button."""
    return sum(cost_of(key, targets, draws, length) for key in keys)
