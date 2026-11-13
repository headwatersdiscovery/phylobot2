"""Scoring mutational trajectories against the truth.

The calibration work scores *residue identity at a matched column* (is this node's
residue right?). A user reading a **mutational trajectory** — the substitutions along an
ancestor->descendant path, as `phylobot.paths` shows them (K300R, -300R, K300-) — feels a
different error: a column the aligner shifted manufactures a *phantom substitution* on one
branch and a compensating one elsewhere, invisible to residue-at-matched-column accuracy
but a real error in the path. This module scores the path.

**The atom is a branch.** A reconstruction has nodes only at tree vertices, so on a branch
parent->child we see only the two endpoints: the net change `parent[col]` -> `child[col]`.
Unobserved history within a branch (a site that went A->B->A on one edge) is invisible on
both sides, by construction, to inferred and true alike. This is the *net readable path*
scope — the same thing the trajectory viewer shows and a resurrection
experiment would test — applied identically to both. A trajectory is just an ordered list
of nodes; its branches are the consecutive pairs, and errors accumulate with depth.

**Node matching** is the golden scorer's (`evaluation.records.truth_records`): internal
nodes by clade on the outgroup-rooted trees (`trees.clades`); extant tips by name. A branch
whose endpoints both match the truth is scored; one that does not is reported, not dropped
(hard nodes the aligner shifted too far to recover exist — report their fraction).

**Columns** map inferred->true the way the *clade* scoring does: each node's own
extant descendants vote (`alignment_map.majority_columns` over the clade), with the whole
alignment's map as a fallback for a column the clade does not cover. Parent and child each
carry their own map; a placement phantom is exactly where the two disagree.

**Per inferred event** (a column where the branch's two endpoints differ), using each
endpoint's residue in the chosen gap version (a residue, or `-` for a gap) and whether it
matches the truth at that node and mapped column:

- `tp`: both endpoints correct and map to the *same* true column, which truly changed.
- `phantom_placement`: both endpoints correct, but their columns map apart (or nowhere) --
  a change the misalignment invented.
- `phantom_residue`: an endpoint's residue is itself wrong -- reconstruction error.

A **true** change between the matched ancestors that no inferred event covers is `missed`
(a substitution the user never sees). Precision = tp / (tp + phantoms) over inferred events;
recall = covered true changes / all true changes. Gains (`-`->X) and losses (X->`-`) are
classified the same way; the driver filters by kind.

Pure over in-memory data (the frames below); reading files is the caller's job. Columns
are 1-based in the public rows (as `phylobot.paths.Change`), 0-based
internally.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

from phylobot.alignment_map import majority_columns
from phylobot.trees import clades, load_tree

GAP = "-"
SUBSTITUTION, GAIN, LOSS = "substitution", "gain", "loss"
TP, PHANTOM_PLACEMENT, PHANTOM_RESIDUE, MISSED = (
    "tp", "phantom_placement", "phantom_residue", "missed",
)  # fmt: skip
PHANTOMS = frozenset({PHANTOM_PLACEMENT, PHANTOM_RESIDUE})


def change_kind(before: str, after: str) -> str:
    """substitution / gain / loss from two differing characters (a residue or `-`)."""
    if before == GAP:
        return GAIN
    if after == GAP:
        return LOSS
    return SUBSTITUTION


@dataclass(frozen=True)
class Frame:
    """One side of the comparison. `extant` are the aligned rows of the observed sequences
    (identical residues on both sides -- only the gaps differ); `ancestors` are the internal
    nodes' rows on the same columns (the inferred side in the chosen gap version); `tree` is
    the Newick whose labels key `ancestors`. Both frames must hold the same `extant`
    sequences (`alignment_map.same_sequences`)."""

    extant: Mapping[str, str]
    ancestors: Mapping[str, str]
    tree: str

    def row(self, node: str) -> str:
        """A node's aligned row, whether it is an ancestor or an extant tip."""
        return self.ancestors[node] if node in self.ancestors else self.extant[node]


@dataclass(frozen=True)
class PathEvent:
    """One scored change on one branch of one trajectory."""

    trajectory: str
    depth: int  # branches from the trajectory's start to this one (0 = first)
    parent: str  # inferred node names
    child: str
    role: str  # "inferred" (a change the reconstruction shows) or "missed" (a true one it lacks)
    classification: str  # tp | phantom_placement | phantom_residue | missed
    kind: str  # substitution | gain | loss
    before: str
    after: str
    inferred_column: int | None  # 1-based; None for a missed change
    true_column: int | None  # 1-based true column (the match's, or a missed change's own)
    parent_correct: bool | None  # residue right for its node (None where the column maps nowhere)
    child_correct: bool | None
    pp: float | None = None  # the child's ML-state posterior at this position (when `inferred_pp`
    # was given): for an inferred change, at its inferred column; for a missed one, at the
    # inferred column homologous to the true change — the confidence of the no-change the
    # reconstruction kept there. None where the position has no residue PP (a gap) or no map.


@dataclass(frozen=True)
class BranchStatus:
    """Every branch the scorer looked at, matched or not, with its change counts (so the
    driver can report the net substitution-count error and the unmatched fraction)."""

    trajectory: str
    depth: int
    parent: str
    child: str
    matched: bool  # both endpoints recovered in the truth
    inferred: dict[str, int] = field(default_factory=dict)  # kind -> count (always present)
    true: dict[str, int] = field(default_factory=dict)  # kind -> count (empty when unmatched)


@dataclass(frozen=True)
class PathScore:
    events: list[PathEvent]
    branches: list[BranchStatus]


class _Columns:
    """Inferred->true column maps, built per node from its own clade and cached, over a
    shared whole-alignment fallback. Mirrors `records.truth_records(columns="clade")`."""

    def __init__(self, inferred: Frame, truth: Frame, label_clade: Mapping[str, frozenset[str]]):
        self._inf = inferred
        self._true = truth
        self._clade = label_clade
        self._global = majority_columns(dict(inferred.extant), dict(truth.extant))
        self._cache: dict[str, dict[int, int]] = {}

    def _node_map(self, node: str) -> dict[int, int]:
        if node not in self._cache:
            clade = self._clade.get(node, frozenset({node}))  # a tip votes with itself
            source = {n: self._inf.extant[n] for n in clade}
            reference = {n: self._true.extant[n] for n in clade}
            self._cache[node] = majority_columns(source, reference)
        return self._cache[node]

    def true_column(self, node: str, column: int) -> int | None:
        """The true column the node places its inferred `column` at (0-based both), by its
        clade, falling back to the whole alignment; None if neither maps it."""
        own = self._node_map(node).get(column)
        return own if own is not None else self._global.get(column)


def _matching(
    inferred_tree: str, true_tree: str, anchor: str
) -> tuple[dict[str, str], dict[str, frozenset[str]]]:
    """(`inferred label -> true label`, `inferred label -> its clade`) for every ingroup
    internal node both trees recover, by clade on the anchor-rooted trees."""
    inferred = clades(load_tree(inferred_tree), anchor)
    true = clades(load_tree(true_tree), anchor)
    match = {label: true[clade] for clade, label in inferred.items() if clade in true}
    # Extant tips are the observed data: they match by name, and are their own "ancestor".
    for tip in load_tree(inferred_tree).leaf_node_iter():
        match[tip.taxon.label] = tip.taxon.label
    label_clade = {label: clade for clade, label in inferred.items()}
    return match, label_clade


def _true_changes(parent_row: str, child_row: str) -> dict[int, tuple[str, str]]:
    """{0-based true column: (before, after)} for every column the two true ancestors
    differ at (the net change on the branch between them)."""
    changes = {}
    for column, (before, after) in enumerate(zip(parent_row, child_row, strict=True)):
        before, after = before.upper(), after.upper()
        if before != after:
            changes[column] = (before, after)
    return changes


def _score_branch(
    trajectory: str,
    depth: int,
    parent: str,
    child: str,
    inferred: Frame,
    truth: Frame,
    match: Mapping[str, str],
    columns: _Columns,
    inferred_pp: Mapping[tuple[str, int], float] | None = None,
) -> tuple[list[PathEvent], BranchStatus]:
    """One branch: classify its inferred changes, then find the true changes it missed."""
    a, b = inferred.row(parent).upper(), inferred.row(child).upper()

    def pp_at(node: str, column0: int) -> float | None:
        return inferred_pp.get((node, column0 + 1)) if inferred_pp is not None else None

    inferred_counts: dict[str, int] = {SUBSTITUTION: 0, GAIN: 0, LOSS: 0}
    for before, after in zip(a, b, strict=True):
        if before != after:
            inferred_counts[change_kind(before, after)] += 1

    true_parent, true_child = match.get(parent), match.get(child)
    if true_parent is None or true_child is None:
        status = BranchStatus(trajectory, depth, parent, child, False, inferred_counts, {})
        return [], status

    true_changes = _true_changes(truth.row(true_parent), truth.row(true_child))
    true_counts: dict[str, int] = {SUBSTITUTION: 0, GAIN: 0, LOSS: 0}
    for before, after in true_changes.values():
        true_counts[change_kind(before, after)] += 1

    events: list[PathEvent] = []
    covered: set[int] = set()
    for column, (before, after) in enumerate(zip(a, b, strict=True)):
        if before == after:
            continue
        tcol_parent = columns.true_column(parent, column)
        tcol_child = columns.true_column(child, column)
        parent_ok = (
            None if tcol_parent is None else before == truth.row(true_parent)[tcol_parent].upper()
        )
        child_ok = (
            None if tcol_child is None else after == truth.row(true_child)[tcol_child].upper()
        )
        true_column: int | None = None
        if parent_ok and child_ok and tcol_parent is not None and tcol_parent == tcol_child:
            classification = TP
            true_column = tcol_parent
            covered.add(tcol_parent)
        elif parent_ok and child_ok:
            classification = PHANTOM_PLACEMENT  # correct residues, columns map apart
        elif parent_ok is None and child_ok is None:
            classification = PHANTOM_PLACEMENT  # an inserted column mapping to no true column
        else:
            classification = PHANTOM_RESIDUE
        events.append(PathEvent(
            trajectory, depth, parent, child, "inferred", classification,
            change_kind(before, after), before, after, column + 1,
            None if true_column is None else true_column + 1, parent_ok, child_ok,
            pp_at(child, column),
        ))  # fmt: skip

    # For the missed changes, map each true column back to the child's inferred column, to read
    # the confidence of the (unchanged) call the reconstruction kept there.
    to_inferred: dict[int, int] = {}
    if inferred_pp is not None and true_changes:
        for column in range(len(b)):
            mapped = columns.true_column(child, column)
            if mapped is not None:
                to_inferred.setdefault(mapped, column)
    for column, (before, after) in sorted(true_changes.items()):
        if column in covered:
            continue
        inferred_col = to_inferred.get(column)
        missed_pp = pp_at(child, inferred_col) if inferred_col is not None else None
        events.append(PathEvent(
            trajectory, depth, parent, child, "missed", MISSED,
            change_kind(before, after), before, after, None, column + 1, None, None, missed_pp,
        ))  # fmt: skip

    status = BranchStatus(trajectory, depth, parent, child, True, inferred_counts, true_counts)
    return events, status


def score_paths(
    inferred: Frame,
    truth: Frame,
    trajectories: Mapping[str, Sequence[str]],
    anchor: str,
    inferred_pp: Mapping[tuple[str, int], float] | None = None,
) -> PathScore:
    """Score each trajectory's branches against the truth. `trajectories` maps a name to an
    ordered node sequence (oldest first), as `phylobot.paths.path` returns; `anchor` is one
    outgroup sequence, which roots both trees for clade matching. `inferred_pp` ({(node,
    1-based column): the node's ML-state posterior}, from ancestors_pp.tsv) fills each event's
    `pp`. Returns every scored change (`PathEvent`) and every branch's status (`BranchStatus`)."""
    match, label_clade = _matching(inferred.tree, truth.tree, anchor)
    columns = _Columns(inferred, truth, label_clade)
    events: list[PathEvent] = []
    branches: list[BranchStatus] = []
    for name, nodes in trajectories.items():
        for depth, (parent, child) in enumerate(zip(nodes, nodes[1:], strict=False)):
            branch_events, status = _score_branch(
                name, depth, parent, child, inferred, truth, match, columns, inferred_pp
            )
            events.extend(branch_events)
            branches.append(status)
    return PathScore(events, branches)


# --- enumerating trajectories from a tree ---------------------------------------------


def internal_branches(newick: str, tips: bool = False) -> dict[str, list[str]]:
    """Every (parent, child) edge whose parent is a labelled internal node, as the file
    roots the tree, each as its own one-branch trajectory `"parent->child"`. The independent
    datum of P1. By default `child` must be a labelled internal node (an internal branch);
    `tips=True` also includes edges to extant tips."""
    tree = load_tree(newick)
    out: dict[str, list[str]] = {}
    for node in tree.preorder_internal_node_iter():
        if not node.label:
            continue
        for kid in node.child_nodes():
            if kid.is_leaf():
                if tips:
                    out[f"{node.label}->{kid.taxon.label}"] = [node.label, kid.taxon.label]
            elif kid.label:
                out[f"{node.label}->{kid.label}"] = [node.label, kid.label]
    return out


def root_to_tip(newick: str, root: str) -> dict[str, list[str]]:
    """The path from internal node `root` down to each extant tip below it, keyed by tip.
    The depth curve's trajectories; overlapping ones share backbone branches (not
    independent -- the branch-level numbers carry the statistics)."""
    from phylobot.paths import descendants, path

    out: dict[str, list[str]] = {}
    for leaf in descendants(newick, root):
        if leaf.extant:
            names = path(newick, root, leaf.name)
            if names:
                out[leaf.name] = names
    return out


def branch_trajectories(names: Iterable[str]) -> dict[str, list[str]]:
    """Helper: wrap an explicit ordered node list as one named trajectory per adjacent pair,
    for callers that already have a path and want each branch scored independently."""
    seq = list(names)
    return {f"{a}->{b}": [a, b] for a, b in zip(seq, seq[1:], strict=False)}
