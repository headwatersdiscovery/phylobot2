"""Mutational paths: the sequences along the tree from an ancestor to a descendant.

The path runs down the tree (away from the root; the outgroup-rooted tree when there
is one), through every ancestor in between, and ends at an ancestor or an extant
sequence. Each ancestor is shown as its ML sequence in the chosen gap version, with the
posterior probability of each residue (ancestors_pp.tsv) for colouring; an extant
sequence is shown as aligned. A cell that differs from the row above has changed on
that branch: a substitution, or a residue gained or lost (the gap version decides where
an ancestor has a gap).

For codon reconstructions (coding DNA under IQ-TREE's codon models or codeml) every row
also has its codons (the version's coding ancestors, the extant codon alignment), so a
substitution is non-synonymous (the amino acid changed) and a codon change that keeps
the amino acid is a synonymous ("silent") change, listed as e.g. L300L. Silent changes
can be hidden; the table can show codons instead of amino acids. An extant codon with
an ambiguity code (e.g. AMC: AAC or ACC) is no change where it may be the codon above
(`codons.compatible`), and a change it may or may not be is "uncertain" (N300X). When codeml's
free-ratio model ran on the combination, each branch also gets its dN/dS estimate
(codeml/FR/branches.tsv; `branch_rates`).
"""

from dataclasses import dataclass, field
from typing import Any

from phylobot import codons as codon_table
from phylobot.ancestors import GAP, colour
from phylobot.trees import load_tree


@dataclass(frozen=True)
class PathNode:
    name: str
    extant: bool
    leaves: int  # extant sequences below it (1 for an extant sequence)


def descendants(newick: str, ancestor: str) -> list[PathNode]:
    """Everything below `ancestor` (labelled internal nodes and extant sequences), in
    preorder. Empty if `ancestor` is not an internal node of the tree."""
    tree = load_tree(newick)
    start = next((n for n in tree.preorder_internal_node_iter() if n.label == ancestor), None)
    if start is None:
        return []
    below = []
    for node in start.preorder_iter():
        if node is start:
            continue
        if node.is_leaf():
            below.append(PathNode(node.taxon.label, True, 1))
        elif node.label:
            below.append(PathNode(node.label, False, len(node.leaf_nodes())))
    return below


def path(newick: str, ancestor: str, descendant: str) -> list[str]:
    """Names from `ancestor` down to `descendant`, both included; [] if `descendant` is
    not below `ancestor`. Unlabelled nodes in between (none in IQ-TREE trees) are skipped."""
    tree = load_tree(newick)
    target = next(
        (n for n in tree.preorder_node_iter()
         if (n.taxon.label if n.is_leaf() else n.label) == descendant),
        None,
    )  # fmt: skip
    names: list[str] = []
    node = target
    while node is not None:
        name = node.taxon.label if node.is_leaf() else (node.label or "")
        if name:
            names.append(name)
        if name == ancestor and not node.is_leaf():
            return names[::-1]
        node = node.parent_node
    return []


def _named(newick: str) -> dict[str, Any]:
    """{name: node} for a tree's extant sequences and labelled internal nodes."""
    found = {}
    for node in load_tree(newick).preorder_node_iter():
        name = node.taxon.label if node.is_leaf() else (node.label or "")
        if name:
            found[name] = node
    return found


def _ancestry(node: Any) -> list[Any]:
    """A node and everything above it, the root last."""
    out = []
    while node is not None:
        out.append(node)
        node = node.parent_node
    return out


def between(newick: str, first: str, last: str) -> list[str]:
    """Names along the tree from `first` to `last`, both included, whichever way the file
    roots the tree: up from `first` to where the two meet, then down to `last`. For an
    ancestor and one of its descendants this is `path`. Used to follow a path through
    another combination's tree, whose root (without an outgroup) may lie anywhere.
    Unlabelled nodes are skipped; [] if either name is not in the tree."""
    nodes = _named(newick)
    if first not in nodes or last not in nodes:
        return []
    up, down = _ancestry(nodes[first]), _ancestry(nodes[last])
    above_last = {id(n) for n in down}
    meet = next(i for i, n in enumerate(up) if id(n) in above_last)
    descent = down[: next(i for i, n in enumerate(down) if n is up[meet])]
    names = {id(node): name for name, node in nodes.items()}
    return [names[id(n)] for n in [*up[: meet + 1], *reversed(descent)] if id(n) in names]


def branch_lengths(newick: str, names: list[str]) -> list[float]:
    """The tree's length between each two consecutive `names` of a path (`path` or
    `between`), summed over any unlabelled node in between; 0 where the file has none."""
    nodes = _named(newick)

    def climb(node: Any, top: Any) -> float:
        total = 0.0
        while node is not top:
            total += float(node.edge.length or 0.0)
            node = node.parent_node
        return total

    lengths = []
    for a, b in zip(names, names[1:], strict=False):
        above_b = {id(n) for n in _ancestry(nodes[b])}
        meet = next(n for n in _ancestry(nodes[a]) if id(n) in above_b)
        lengths.append(climb(nodes[a], meet) + climb(nodes[b], meet))
    return lengths


def ml_probabilities(pp_tsv: str, nodes: set[str]) -> dict[str, list[float]]:
    """{node: [posterior probability of the ML state, per site]} from ancestors_pp.tsv."""
    by_node: dict[str, list[tuple[int, float]]] = {}
    for line in pp_tsv.splitlines()[1:]:
        node, site, _state, pp = line.split("\t")
        if node in nodes:
            by_node.setdefault(node, []).append((int(site), float(pp)))
    return {node: [pp for _, pp in sorted(sites)] for node, sites in by_node.items()}


@dataclass(frozen=True)
class PathCell:
    char: str  # the residue, or the codon when the table shows codons
    pp: float | None = None
    style: str = ""
    changed: bool = False  # amino acid differs from the row above (path rows only)
    silent: bool = False  # same amino acid as the row above, another codon (shown silent)
    codon: str = ""  # the cell's codon, when codons are known


@dataclass
class PathRow:
    name: str
    extant: bool
    on_path: bool  # False for extant sequences added for comparison
    cells: list[PathCell] = field(default_factory=list)


@dataclass(frozen=True)
class Change:
    column: int  # 1-based alignment column
    before: str  # "-" for a gain
    after: str  # "-" for a loss
    pp_before: float | None
    pp_after: float | None
    codon_before: str = ""  # with codons known; "" otherwise
    codon_after: str = ""

    @property
    def label(self) -> str:
        """E.g. K300R, -300R (gained), K300- (lost), L300L (silent); numbered by
        alignment column."""
        return f"{self.before}{self.column}{self.after}"

    @property
    def kind(self) -> str:
        """gain, loss, substitution (the amino acid changed: non-synonymous when codons
        are known), synonymous (same amino acid, another codon) or uncertain (a codon
        with an ambiguity code, translated as X: it may or may not change the amino acid)."""
        if self.before == GAP:
            return "gain"
        if self.after == GAP:
            return "loss"
        ambiguous = any(codon_table.is_ambiguous(c) for c in (self.codon_before, self.codon_after))
        if ambiguous and "X" in (self.before, self.after):
            return "uncertain"
        return "synonymous" if self.before == self.after else "substitution"

    @property
    def codons(self) -> str:
        """E.g. "AAA → AGA"; "" without codons."""
        if not (self.codon_before or self.codon_after):
            return ""
        return f"{self.codon_before or '---'} → {self.codon_after or '---'}"


@dataclass(frozen=True)
class BranchRate:
    """codeml's free-ratio estimate for one branch: dN/dS, dN, dS, and the expected
    numbers of non-synonymous (N*dN) and synonymous (S*dS) changes on it."""

    omega: float
    dn: float
    ds: float
    n_dn: float
    s_ds: float

    @property
    def defined(self) -> bool:
        """Whether dN/dS means anything: at least one synonymous change expected (with
        dS near 0 codeml reports 999, or another huge ratio)."""
        return self.ds > 0 and self.s_ds >= 1


def branch_rates(tsv: str) -> dict[frozenset[str], BranchRate]:
    """{the branch's two end nodes: its estimate} from codeml/FR/branches.tsv. Keyed by
    the unordered pair, since codeml's tree and ours may be rooted differently."""
    rows = [line.split("\t") for line in tsv.splitlines()]
    header = rows[0]
    rates = {}
    for row in rows[1:]:
        value = dict(zip(header, row, strict=True))
        rates[frozenset((value["parent"], value["child"]))] = BranchRate(
            *(float(value[k]) for k in ("dN_dS", "dN", "dS", "N_dN", "S_dS"))
        )
    return rates


@dataclass(frozen=True)
class BranchChanges:
    parent: str
    child: str
    changes: list[Change]  # shown ones: silent changes only when shown
    hidden_silent: int = 0  # silent changes on this branch left out of `changes`
    rate: BranchRate | None = None  # codeml's free-ratio estimate, when it ran

    def count(self, kind: str) -> int:
        return sum(1 for c in self.changes if c.kind == kind)

    @property
    def substitutions(self) -> int:
        return self.count("substitution")

    @property
    def gains(self) -> int:
        return self.count("gain")

    @property
    def losses(self) -> int:
        return self.count("loss")

    @property
    def synonymous(self) -> int:
        return self.count("synonymous") + self.hidden_silent

    @property
    def uncertain(self) -> int:
        return self.count("uncertain")


@dataclass
class PathTable:
    columns: list[int]  # 1-based alignment columns shown
    first_numbers: list[int | None]  # residue number in the path's first sequence, per shown column
    rows: list[PathRow]
    branches: list[BranchChanges]
    changed_columns: int
    codons: bool = False  # codons known: changes are classed (non-)synonymous
    show_silent: bool = True
    cells: str = "amino"  # what the cells show: "amino" acids or "codons"

    @property
    def uncertain(self) -> bool:
        """Some change is uncertain (an extant codon with an ambiguity code)."""
        return any(b.uncertain for b in self.branches)


def trace(
    names: list[str],
    sequences: dict[str, str],
    probabilities: dict[str, list[float]],
    extant: set[str],
    compare: list[str],
    palette: str,
    only_changes: bool = True,
    codons: dict[str, str] | None = None,
    codon_probabilities: dict[str, list[float]] | None = None,
    show_silent: bool = True,
    cells: str = "amino",
    rates: dict[frozenset[str], BranchRate] | None = None,
) -> PathTable:
    """The path table. `names` is the path (oldest first); `sequences` holds every row's
    aligned sequence (ancestors in the chosen gap version); `probabilities` the ML
    residues' posterior probabilities of ancestors (none for PRANK/GRASP versions);
    `compare` extant sequences shown beneath for comparison. `codons`, for codon
    reconstructions, holds every row's coding sequence (three nucleotides per column):
    changes are then classed (non-)synonymous, silent ones listed if `show_silent`, and
    `cells="codons"` shows codons coloured by `codon_probabilities` (the ML codon's).
    `rates` (branch_rates) gives each branch codeml's free-ratio dN/dS where it ran
    (a branch split by the outgroup root has none)."""
    columns = len(sequences[names[0]])
    for name in [*names, *compare]:
        if len(sequences[name]) != columns:
            raise ValueError(f"{name}: {len(sequences[name])} columns, expected {columns}")
        if codons is not None and len(codons[name]) != 3 * columns:
            raise ValueError(f"{name}: {len(codons[name])} nucleotides, expected {3 * columns}")
    if codons is None:
        cells = "amino"
    shown_pp = (codon_probabilities or {}) if cells == "codons" else probabilities
    # Probabilities that do not cover every column are not used (cells stay uncoloured).
    shown_pp = {k: v for k, v in shown_pp.items() if len(v) == columns}
    probabilities = {k: v for k, v in probabilities.items() if len(v) == columns}

    def upper(name: str) -> str:
        return sequences[name].upper()

    def codon(name: str, c: int) -> str:
        """The row's codon at column c ("" without codons or at a gap)."""
        if codons is None:
            return ""
        triplet = codons[name][3 * c : 3 * c + 3].upper()
        return "" if GAP in triplet else triplet

    def same_codon(ca: str, cb: str) -> bool:
        """Codons that may be the same (one with an ambiguity code, e.g. AAC and AMC)."""
        return bool(ca and cb) and ca != cb and codon_table.compatible(ca, cb)

    branches = []
    changed: set[int] = set()
    for parent, child in zip(names, names[1:], strict=False):
        a, b = upper(parent), upper(child)
        pa, pb = probabilities.get(parent), probabilities.get(child)
        diffs = []
        silent = 0
        for c in range(columns):
            ca, cb = codon(parent, c), codon(child, c)
            if (a[c] == b[c] and ca == cb) or same_codon(ca, cb):
                continue
            if a[c] == b[c] and not show_silent:
                silent += 1
                continue
            diffs.append(Change(c + 1, a[c], b[c], pa[c] if pa else None,
                                pb[c] if pb else None, ca, cb))  # fmt: skip
        changed |= {d.column - 1 for d in diffs}
        rate = (rates or {}).get(frozenset((parent, child)))
        branches.append(BranchChanges(parent, child, diffs, silent, rate))
    shown = sorted(changed) if only_changes else list(range(columns))

    def row_cells(name: str, previous: str | None) -> list[PathCell]:
        row, pp = upper(name), shown_pp.get(name)
        out = []
        for c in shown:
            char = row[c]
            triplet = codon(name, c)
            differs = (
                previous is not None and previous[c] != char
                and not same_codon(codon(previous_name, c), triplet)
            )  # fmt: skip
            quiet = (
                show_silent and previous is not None and not differs
                and triplet != codon(previous_name, c)
            )  # fmt: skip
            text = triplet if cells == "codons" and triplet else char
            if char == GAP:
                out.append(PathCell("–", changed=differs))
            elif pp is not None and name not in extant:
                out.append(PathCell(text, pp[c], colour(pp[c], palette), differs, quiet, triplet))
            else:
                out.append(PathCell(text, changed=differs, silent=quiet, codon=triplet))
        return out

    rows = []
    previous = None
    previous_name = ""
    for name in names:
        rows.append(PathRow(name, name in extant, True, row_cells(name, previous)))
        previous, previous_name = upper(name), name
    for name in compare:
        rows.append(PathRow(name, True, False, row_cells(name, None)))

    numbers: list[int | None] = []
    count = 0
    first = upper(names[0])
    by_column: list[int | None] = []
    for char in first:
        if char == GAP:
            by_column.append(None)
        else:
            count += 1
            by_column.append(count)
    numbers = [by_column[c] for c in shown]
    return PathTable(
        [c + 1 for c in shown], numbers, rows, branches, len(changed),
        codons is not None, show_silent, cells,
    )  # fmt: skip
