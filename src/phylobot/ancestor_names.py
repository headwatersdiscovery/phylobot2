"""Names users give ancestors ("Anc.GFP"), shown on every combination that has them.

A node's label (Node4) is one tree's numbering; the same ancestor is Node9 in another
combination. A name therefore belongs to the ancestor's **clade**: the content-derived ID
`ancestor_map.tsv` gives it (`phylobot.analogs`; the clade ID with an outgroup, else the
ID of the subtrees meeting at it). Every combination whose node has that ID shows the
name; where the clade matches but the branching inside differs from the node it was
named on, the name is marked approximate.

Names are user text, so they are held to a character set that is safe as a FASTA
identifier and in a file name, never reach a pipeline file, a storage path or a Nextflow
`meta`, and are applied only where the web tier displays and exports (`for_selection`).
An extension shares its parent's names: they are stored on the topmost job it inherits
from (`holder`).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from django.db import IntegrityError, transaction

from phylobot import analog_views, joined
from phylobot.models import AncestorName, AncestorNote, Job

MAX_LENGTH = 40
PATTERN = re.compile(rf"[A-Za-z][A-Za-z0-9._-]{{0,{MAX_LENGTH - 1}}}")
NODE_LABEL = re.compile(r"node\d+", re.IGNORECASE)  # the trees' own labels: confusing
RULES = f"up to {MAX_LENGTH} characters: letters, digits, '.', '_' and '-', starting with a letter"


@dataclass(frozen=True)
class Named:
    name: str
    alike: bool  # the same subtrees meet here as where it was named, not only the clade


def holder(job: Job) -> Job:
    """The job whose names this one shows: the topmost job it inherits results from."""
    while (parent := joined.inheritance(job)[0]) is not None:
        job = parent
    return job


def problem(name: str, sequences: set[str]) -> str:
    """Why `name` cannot be used ("" if it can). `sequences`: the job's sequence names,
    which a name must not repeat (a FASTA of both would be ambiguous)."""
    if not PATTERN.fullmatch(name):
        return f"A name is {RULES}."
    if NODE_LABEL.fullmatch(name):
        return "That is how the trees label nodes: choose a name of your own."
    if name.casefold() in {s.casefold() for s in sequences}:
        return "That is one of the job's sequence names."
    return ""


def save(job: Job, clade: str, parts: str, name: str) -> str:
    """Name a clade (renaming it if it has a name); "" or why not. Validate first."""
    owner = holder(job)
    clash = AncestorName.objects.filter(job=owner, name__iexact=name).exclude(clade=clade)
    if clash.exists():
        return "Another ancestor of this job already has that name."
    try:
        with transaction.atomic():
            AncestorName.objects.update_or_create(
                job=owner, clade=clade, defaults={"name": name, "parts": parts}
            )
    except IntegrityError:  # named meanwhile
        return "Another ancestor of this job already has that name."
    return ""


def exist(job: Job) -> bool:
    """Whether any ancestor of the job (or what it inherits from) has a name."""
    return AncestorName.objects.filter(job=holder(job)).exists()


def remove(job: Job, clade: str) -> None:
    AncestorName.objects.filter(job=holder(job), clade=clade).delete()


def match_key(place: Any, by_clade: bool) -> str:
    """The ID a name is stored under, for a node of the ancestor map."""
    return str(place.below if by_clade else place.parts)


def for_nodes(job: Job, nodes: dict[str, Any], by_clade: bool) -> dict[str, Named]:
    """{node: its name} for one combination's nodes ({node: analogs.Located})."""
    names = {n.clade: n for n in AncestorName.objects.filter(job=holder(job))}
    if not names:
        return {}
    out = {}
    for node, place in nodes.items():
        found = names.get(match_key(place, by_clade))
        if found is not None:
            out[node] = Named(found.name, found.parts == place.parts)
    return out


def for_selection(selection: Any, combination: str = "") -> dict[str, Named]:
    """{node: its name} for a viewer selection's combination (or another of the job's),
    from the job's ancestor map; {} where it has no map."""
    placed = analog_views.trees(selection)
    nodes = placed.get(combination or selection.combination["key"], {})
    return for_nodes(selection.job, nodes, analog_views.by_clade(selection)) if nodes else {}


def label(node: str, names: dict[str, Named]) -> str:
    """ "Anc.GFP (Node4)", or the node alone."""
    found = names.get(node)
    return f"{found.name} ({node})" if found else node


def fasta_id(node: str, names: dict[str, Named]) -> str:
    """The identifier a downloaded FASTA record carries: the name if it has one."""
    found = names.get(node)
    return found.name if found else node


def save_note(job: Job, clade: str, parts: str, text: str) -> None:
    """Describe a clade; empty text removes the description."""
    owner = holder(job)
    text = text.replace("\r\n", "\n").strip()[: AncestorNote.LIMIT]
    if not text:
        AncestorNote.objects.filter(job=owner, clade=clade).delete()
        return
    AncestorNote.objects.update_or_create(
        job=owner, clade=clade, defaults={"text": text, "parts": parts}
    )


def note_for(job: Job, place: Any, by_clade: bool) -> str:
    """The description of the ancestor at `place` (an ancestor-map node), or ""."""
    found = AncestorNote.objects.filter(job=holder(job), clade=match_key(place, by_clade)).first()
    return found.text if found else ""
