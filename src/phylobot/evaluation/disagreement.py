"""Cross-aligner disagreement on mutational paths.

Where `evaluation.paths` scores a path against the *truth*, this scores how far several
aligners *agree* about it, which is all one can measure on real data; on simulated data
it can be set against true error. Several reconstructions of the same
sequences -- one per aligner -- are laid on one reference column frame (each aligner's
columns mapped to the reference's by `alignment_map`, through the extant residues they
share). On a branch, each aligner reports a set of changes; a change is **aligner-
dependent** when the aligners do not all make the same call at that reference column
(a different before/after, or one aligner shows no change there at all).

Two ways the "same branch" is identified across aligners:

- **Fixed topology**: every aligner reconstructs on the *same* tree, so a branch is
  the same (parent, child) labels in each -- node matching is by construction, and the
  disagreement is purely alignment-driven.
- **Clade matching**: each aligner searched its own tree, so a branch is matched
  by its (parent clade, child clade) on the outgroup-rooted trees (`trees.clades`), the
  way `phylobot.analogs` matches nodes. A branch only counts where every aligner has it;
  the disagreement then mixes alignment with tree, which the report states.

Disagreement is **not error**: four aligners can agree and all be wrong (the simulation
only a simulation can speak to accuracy). Pure over in-memory frames; columns 0-based
internally, 1-based in the public rows.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phylobot.alignment_map import column_map
from phylobot.evaluation.paths import GAP, Frame, change_kind, internal_branches
from phylobot.trees import clades, load_tree

Change = tuple[str, str]  # (before, after), each a residue or "-"
Branch = tuple[str, str]  # a node-label pair, or used as a canonical id


@dataclass(frozen=True)
class Call:
    """One reference substitution (or indel) on one branch, and whether the aligners agree
    about it."""

    branch: str  # a readable id for the branch (reference labels, or clade sizes)
    reference_column: int  # 1-based, in the reference aligner's alignment
    kind: str  # substitution | gain | loss
    before: str
    after: str
    aligner_dependent: bool  # not every aligner makes this same call at this column
    agreeing: int  # how many aligners make exactly this call
    aligners: int  # how many aligners were compared
    branch_length: float  # the reference tree's length for this branch


def _to_reference(frame: Frame, reference: Frame) -> dict[int, int]:
    """Each column of `frame` -> the reference column where most of its residues sit."""
    return column_map(dict(frame.extant), dict(reference.extant))


def branch_changes(
    frame: Frame, parent: str, child: str, to_reference: dict[int, int]
) -> dict[int, Change]:
    """{reference column: (before, after)} for every column the branch's two endpoints
    differ at, in the chosen gap version. A column that maps to no reference column (an
    insertion the reference alignment lacks) is dropped: there is nothing to compare it to."""
    a, b = frame.row(parent).upper(), frame.row(child).upper()
    out: dict[int, Change] = {}
    for column, (before, after) in enumerate(zip(a, b, strict=True)):
        if before == after:
            continue
        reference_column = to_reference.get(column)
        if reference_column is not None:
            out[reference_column] = (before, after)
    return out


def _branch_lengths(newick: str) -> dict[Branch, float]:
    """{(parent label, child label): branch length} as the file roots the tree."""
    tree = load_tree(newick)
    out: dict[Branch, float] = {}
    for node in tree.preorder_internal_node_iter():
        if not node.label:
            continue
        for kid in node.child_nodes():
            name = kid.taxon.label if kid.is_leaf() else (kid.label or "")
            if name:
                out[(node.label, name)] = float(kid.edge.length or 0.0)
    return out


def fixed_topology_branches(reference_tree: str) -> dict[str, dict[str, Branch]]:
    """Every internal branch of the shared tree, as {branch id: {aligner: (parent, child)}}
    where the labels are the same for every aligner (the fixed topology). The aligner keys
    are filled in by `calls` from the reconstructions it is given."""
    branches = internal_branches(reference_tree)
    # One placeholder "" key; `calls` copies the shared labels to every aligner.
    return {name: {"": (pair[0], pair[1])} for name, pair in branches.items()}


def clade_branches(
    reconstructions: Mapping[str, Frame], anchor: str
) -> tuple[dict[str, dict[str, Branch]], str]:
    """Match internal branches across aligners by (parent clade, child clade) on the
    outgroup-rooted trees. Returns {branch id: {aligner: (parent label, child label)}} for
    the branches every aligner recovers, and the reference aligner (the first), whose column
    frame the calls are laid on."""
    aligners = list(reconstructions)
    CladeBranch = tuple[frozenset[str], frozenset[str]]
    per_aligner_branches: dict[str, dict[CladeBranch, Branch]] = {}
    for aligner in aligners:
        tree = reconstructions[aligner].tree
        label_clade = {label: clade for clade, label in clades(load_tree(tree), anchor).items()}
        pairs: dict[CladeBranch, Branch] = {}
        for p_label, c_label in internal_branches(tree).values():
            p_clade, c_clade = label_clade.get(p_label), label_clade.get(c_label)
            if p_clade is not None and c_clade is not None:
                pairs[(p_clade, c_clade)] = (p_label, c_label)
        per_aligner_branches[aligner] = pairs
    shared = (
        set.intersection(*(set(b) for b in per_aligner_branches.values())) if aligners else set()
    )
    out: dict[str, dict[str, Branch]] = {}
    for key in shared:
        p_clade, c_clade = key
        name = f"{len(p_clade)}>{len(c_clade)}:{min(c_clade)}"  # a stable, readable id
        out[name] = {a: per_aligner_branches[a][key] for a in aligners}
    return out, aligners[0]


def calls(
    reconstructions: Mapping[str, Frame],
    branch_map: Mapping[str, Mapping[str, Branch]],
    reference: str,
    *,
    kinds: tuple[str, ...] = ("substitution", "gain", "loss"),
) -> list[Call]:
    """For every branch in `branch_map`, every change the *reference* aligner reports (of
    the given kinds), with whether the aligners agree about it. `branch_map` maps a branch
    id to each aligner's (parent, child) labels (from `fixed_topology_branches`, its aligner
    keys replaced, or `clade_branches`). The reference aligner's column frame and branch
    lengths are used."""
    aligners = list(reconstructions)
    to_reference = {
        a: _to_reference(reconstructions[a], reconstructions[reference]) for a in aligners
    }
    lengths = _branch_lengths(reconstructions[reference].tree)
    rows: list[Call] = []
    for name, by_aligner in branch_map.items():
        if reference not in by_aligner:
            # fixed-topology maps carry one placeholder key; the labels are shared.
            labels = next(iter(by_aligner.values()))
            by_aligner = dict.fromkeys(aligners, labels)
        per_aligner = {
            a: branch_changes(
                reconstructions[a], by_aligner[a][0], by_aligner[a][1], to_reference[a]
            )
            for a in aligners
        }
        ref_labels = by_aligner[reference]
        length = lengths.get((ref_labels[0], ref_labels[1]), 0.0)
        for column, (before, after) in per_aligner[reference].items():
            kind = change_kind(before, after)
            if kind not in kinds:
                continue
            here = {a: per_aligner[a].get(column) for a in aligners}
            agreeing = sum(1 for call in here.values() if call == (before, after))
            rows.append(Call(
                name, column + 1, kind, before, after,
                aligner_dependent=agreeing < len(aligners), agreeing=agreeing,
                aligners=len(aligners), branch_length=length,
            ))  # fmt: skip
    return rows


def branch_disagreement(
    reconstructions: Mapping[str, Frame],
    branch_map: Mapping[str, Mapping[str, Branch]],
    reference: str,
) -> dict[str, int]:
    """{branch id: how many of its reference changes are aligner-dependent}. The per-branch
    disagreement score that can be set against true error."""
    score: dict[str, int] = {}
    for call in calls(reconstructions, branch_map, reference):
        score[call.branch] = score.get(call.branch, 0) + int(call.aligner_dependent)
    return score


__all__ = [
    "Call",
    "GAP",
    "branch_changes",
    "branch_disagreement",
    "calls",
    "clade_branches",
    "fixed_topology_branches",
]
