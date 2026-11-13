"""Gap placement in ancestral sequences: decide which ancestral columns were absent.

IQ-TREE reconstructs a residue at every alignment column, so ancestors come out as long
as the whole alignment. Each method here decides, per ancestral node and column,
whether the ancestor had a residue; masking replaces the others with '-'. Masking never
changes the alignment, the tree, or any residue or its probability.

Masking methods (residues unchanged):
- parsimony, "shared presence": at a node, a column holds a residue if residues occur
  in at least two of the subtrees that meet there. On a rooted tree that means "on both
  sides of the node", so an insertion in one lineage is never pushed back into its
  ancestors. Rooting-free.
- fitch, "Fitch parsimony": the classic two-pass algorithm on the gap/residue pattern
  of each column. It needs a root (the outgroup-rooted tree if there is one, otherwise
  the tree's top node); nodes with more than two children use Hartigan's generalisation.
  Ties are resolved toward a gap at the root and toward the parent's state below it.
- ml: the alignment recoded as gap (0) / residue (1), reconstructed by IQ-TREE on the
  same tree (binary GTR2, branch lengths rescaled by one fitted factor); a column is a
  gap where P(gap) > 0.5. The probabilities are kept (gaps-ml.tsv).

Alternative reconstructions (their own residues too; same alignment columns):
- prank: PRANK's ancestors on the combination's alignment (kept as is) and its tree
  rooted on the outgroup (their indel models depend on the root, so both need one).
- grasp: GRASP's joint reconstruction with indel inference (BEP), on the same inputs.
Their node names are matched to ours by the subtrees that meet at each node.
"""

from pathlib import Path
from typing import Any

from phylobot.seqio import read_fasta
from phylobot.trees import load_tree

GAP = "-"
METHODS = ("parsimony", "fitch", "ml", "prank", "grasp")
MASKING = ("parsimony", "fitch", "ml")  # the rest reconstruct ancestors themselves
NEED_OUTGROUP = ("prank", "grasp")  # their indel models depend on the root


def default_methods(has_outgroup: bool) -> list[str]:
    """The default gap placement (owner decision, from the ime2/indels replicates): Fitch
    when an outgroup roots the tree, shared presence (rooting-free) otherwise. Keep in
    step with main.nf and the New Job page's script."""
    return ["fitch"] if has_outgroup else ["parsimony"]


def parsimony_presence(newick: str, alignment: dict[str, str]) -> dict[str, list[bool]]:
    """{internal node label: [residue present? per column]} by the shared-presence rule.

    `alignment` holds the leaves' aligned rows, keyed by the tree's leaf names.
    """
    tree = load_tree(newick)
    columns = len(next(iter(alignment.values())))
    everyone = frozenset(alignment)
    # For every column, the set of leaves with a residue there.
    present_in: list[frozenset[str]] = [
        frozenset(name for name, row in alignment.items() if row[col] != GAP)
        for col in range(columns)
    ]
    result: dict[str, list[bool]] = {}
    for node in tree.postorder_internal_node_iter():
        # The subtrees meeting at this node: each child's leaves, plus (unless this is
        # the root) everything outside the node.
        below = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        subtrees = [
            frozenset(leaf.taxon.label for leaf in child.leaf_iter())
            for child in node.child_nodes()
        ]
        if below != everyone:
            subtrees.append(everyone - below)
        result[node.label or ""] = [
            sum(1 for subtree in subtrees if subtree & present) >= 2 for present in present_in
        ]
    return result


def fitch_presence(rooted_newick: str, alignment: dict[str, str]) -> dict[str, list[bool]]:
    """{labelled internal node: [residue present? per column]} by Fitch parsimony on the
    gap (False) / residue (True) pattern, using the tree's root as given. An unlabelled
    root (e.g. one added by outgroup rooting) is used but not reported."""
    tree = load_tree(rooted_newick)
    states = fitch_nodes(tree, alignment)
    return {
        node.label: states[id(node)]
        for node in reversed(list(tree.postorder_node_iter()))
        if not node.is_leaf() and node.label
    }


def fitch_nodes(tree: Any, alignment: dict[str, str]) -> dict[int, list[bool]]:
    """{id(internal node): [residue present? per column]}, every internal node of `tree`
    (a dendropy tree), labelled or not."""
    columns = len(next(iter(alignment.values())))
    result: dict[int, list[bool]] = {}
    nodes = list(tree.postorder_node_iter())
    for col in range(columns):
        # Up pass: candidate state sets (Hartigan: states in the most children's sets).
        sets: dict[int, frozenset[bool]] = {}
        for node in nodes:
            if node.is_leaf():
                sets[id(node)] = frozenset({alignment[node.taxon.label][col] != GAP})
                continue
            counts = {state: 0 for state in (False, True)}
            for child in node.child_nodes():
                for state in sets[id(child)]:
                    counts[state] += 1
            best = max(counts.values())
            sets[id(node)] = frozenset(state for state, n in counts.items() if n == best)
        # Down pass: the root prefers a gap on a tie; below, keep the parent's state
        # when it is a candidate, otherwise prefer a gap on a tie.
        chosen: dict[int, bool] = {}
        for node in reversed(nodes):  # preorder: parents before children
            if node.is_leaf():
                continue
            candidates = sets[id(node)]
            parent = node.parent_node
            if parent is not None and chosen[id(parent)] in candidates:
                chosen[id(node)] = chosen[id(parent)]
            else:
                chosen[id(node)] = candidates == frozenset({True})
            result.setdefault(id(node), []).append(chosen[id(node)])
    return result


def read_binary_state(path: Path) -> dict[str, list[float]]:
    """{node: [P(gap) per column]} from IQ-TREE's .state for a gap/residue (0/1) alignment."""
    probabilities: dict[str, list[float]] = {}
    for line in path.read_text().splitlines():
        if line.startswith("#") or line.startswith("Node\t") or not line.strip():
            continue
        node, _site, _state, p_gap, _p_residue = line.split("\t")
        probabilities.setdefault(node, []).append(float(p_gap))
    return probabilities


def node_partitions(newick: str) -> dict[frozenset[frozenset[str]], str]:
    """{the subtrees meeting at an internal node (as leaf sets): its label}.

    This identifies a node in an unrooted sense, independent of where a tree is rooted:
    a root added by splitting a branch (as PRANK does) meets only two subtrees and
    matches no node of an unrooted tree.
    """
    tree = load_tree(newick)
    everyone = frozenset(leaf.taxon.label for leaf in tree.leaf_node_iter())
    result: dict[frozenset[frozenset[str]], str] = {}
    for node in tree.postorder_internal_node_iter():
        below = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        parts = {frozenset(leaf.taxon.label for leaf in c.leaf_iter()) for c in node.child_nodes()}
        if below != everyone:
            parts.add(everyone - below)
        result[frozenset(parts)] = node.label or ""
    return result


def match_nodes(target_newick: str, source_newick: str) -> dict[str, str]:
    """{target node label: source node label} for internal nodes where the same subtrees
    meet. Both trees must have the same leaves; names and rooting may differ."""
    source = node_partitions(source_newick)
    return {
        label: source[parts]
        for parts, label in node_partitions(target_newick).items()
        if parts in source
    }


def all_gap_columns(alignment: dict[str, str]) -> list[int]:
    """Columns where every sequence has a gap (e.g. a residue lost in every lineage)."""
    columns = len(next(iter(alignment.values())))
    return [c for c in range(columns) if all(row[c] == GAP for row in alignment.values())]


def restore_columns(row: str, dropped: list[int]) -> str:
    """Re-insert gaps at `dropped` column positions (ascending) into `row`."""
    chars = list(row)
    for column in dropped:
        chars.insert(column, GAP)
    return "".join(chars)


def external_ancestors(
    asr_newick: str, their_newick: str, their_fasta: Path, alignment: dict[str, str]
) -> dict[str, str]:
    """Another program's ancestors (PRANK, GRASP), renamed to our node names.

    Their rows must be aligned to the same alignment. Some programs drop all-gap columns
    (PRANK does, even with -keep); rows short by exactly those get them back as gaps, the
    only defensible call where no sequence has a residue. Any other length mismatch is
    an error. Nodes of theirs that match none of ours (e.g. an added root) are dropped.
    """
    columns = len(next(iter(alignment.values())))
    dropped = all_gap_columns(alignment)
    rows = {name.split()[0]: row.upper() for name, row in read_fasta(their_fasta).items()}
    renamed = {}
    for ours, theirs in match_nodes(asr_newick, their_newick).items():
        if theirs not in rows:
            continue
        row = rows[theirs]
        if dropped and len(row) == columns - len(dropped):
            row = restore_columns(row, dropped)
        if len(row) != columns:
            raise ValueError(f"{theirs}: {len(row)} columns, alignment has {columns}")
        renamed[ours] = row
    return renamed


def ml_presence(
    asr_newick: str, gap_newick: str, gap_state: Path
) -> tuple[dict[str, list[bool]], dict[str, list[float]]]:
    """Presence per residue-ASR node from the gap reconstruction, and P(gap) per column.
    The gap reconstruction's node names are matched to the residue ASR's by split."""
    p_gap = read_binary_state(gap_state)
    matched = match_nodes(asr_newick, gap_newick)
    probabilities = {node: p_gap[src] for node, src in matched.items() if src in p_gap}
    presence = {node: [p <= 0.5 for p in probs] for node, probs in probabilities.items()}
    return presence, probabilities


def mask(ancestors: dict[str, str], presence: dict[str, list[bool]]) -> dict[str, str]:
    """Replace residues with '-' where the ancestor had none. Nodes without presence
    information are left unmasked."""
    masked = {}
    for node, row in ancestors.items():
        present = presence.get(node)
        if present is None:
            masked[node] = row
            continue
        if len(present) != len(row):
            raise ValueError(f"{node}: {len(row)} columns but presence for {len(present)}")
        masked[node] = "".join(c if keep else GAP for c, keep in zip(row, present, strict=True))
    return masked
