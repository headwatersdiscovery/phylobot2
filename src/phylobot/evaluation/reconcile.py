"""Score a reconciled combination against a simulated gene family.

The truth is a simulated family (`phylobot.simulate.families`): the rooted true gene tree
(internal nodes Node1, ... as in true_ancestors.fasta) and true_events.tsv, whose
duplications name the gene node they became ("-" when pruned away). A combination is a
pipeline result folder with reconciliation.nwk (GeneRax's rooting of the combination's
tree, our node labels) and reconciliation.tsv (each node's event, plus the "(root)"
row: GeneRax's root lies on one of our branches).

Inferred clades are taken in GeneRax's rooting, as the job's pages take them
(`reconcile.paralogy`), so a wrong root shows up as wrong clades, which is what a user
would see. Measures:

- `rf`: unrooted Robinson-Foulds distance of the combination's tree to the true tree;
- duplication recall and precision: true and inferred duplications matched by the
  leaves below them;
- `chimeric`: inferred ancestors whose clade crosses a true duplication, i.e. holds
  genes from two of its copies without holding the whole duplication. No such gene ever
  existed; this is the count reconciliation exists to reduce;
- `identity`: mean identity of the inferred ancestors whose clade is a true node's,
  against that node's true sequence;
- the paralog-mixing flag (`reconcile.paralogy`, comparing a tree with the reconciled
  tree of the same alignment and model): which ancestors it flags, against `chimeric`.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field

from phylobot import reconcile
from phylobot.evaluation.gaps import aligned_identity
from phylobot.evaluation.sequences import identity
from phylobot.trees import load_tree, robinson_foulds

Leaves = frozenset[str]


@dataclass(frozen=True)
class Duplication:
    node: str  # the true gene node it became
    below: Leaves
    copies: tuple[Leaves, ...]  # the leaves of each copy (its children)


@dataclass
class Truth:
    newick: str
    clades: dict[Leaves, str]  # every internal node's leaves -> its name
    duplications: list[Duplication]
    ancestors: Mapping[str, str] = field(default_factory=dict)  # name -> true sequence


def truth(newick: str, events_tsv: str, ancestors: Mapping[str, str] | None = None) -> Truth:
    """The true gene tree and its visible duplications."""
    tree = load_tree(newick)
    clades: dict[Leaves, str] = {}
    children: dict[str, tuple[Leaves, ...]] = {}
    for node in tree.preorder_internal_node_iter():
        below = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        clades[below] = node.label or ""
        children[node.label or ""] = tuple(
            frozenset(x.taxon.label for x in c.leaf_iter()) for c in node.child_nodes()
        )
    by_name = {name: below for below, name in clades.items()}
    header, *rows = events_tsv.splitlines()
    columns = header.split("\t")
    duplications = []
    for line in rows:
        row = dict(zip(columns, line.split("\t"), strict=True))
        node = row["gene_node"]
        if row["event"] == "duplication" and node != "-":
            duplications.append(Duplication(node, by_name[node], children[node]))
    return Truth(newick, clades, duplications, dict(ancestors or {}))


def crosses(clade: Leaves, duplication: Duplication) -> bool:
    """A clade holding genes of two copies of a duplication, but not all of it."""
    touched = sum(1 for copy in duplication.copies if clade & copy)
    return touched >= 2 and not duplication.below <= clade


def chimeric_nodes(clades: Mapping[str, Leaves], true: Truth) -> set[str]:
    """The named inferred ancestors that mix paralogs of a true duplication."""
    return {
        name
        for name, below in clades.items()
        if name != reconcile.ROOT and any(crosses(below, d) for d in true.duplications)
    }


@dataclass
class Score:
    rf: int
    ancestors: int  # inferred ancestors scored (named internal nodes)
    duplications_true: int
    duplications_inferred: int
    duplications_found: int  # inferred duplications on a true duplication's clade
    chimeric: list[str]
    clades_matched: int  # inferred ancestors whose clade is a true node's
    identity: float | None  # their mean identity to the truth (None without sequences)
    losses: int

    @property
    def recall(self) -> float | None:
        return self.duplications_found / self.duplications_true if self.duplications_true else None

    @property
    def precision(self) -> float | None:
        if not self.duplications_inferred:
            return None
        return self.duplications_found / self.duplications_inferred


def sequence_identity(inferred: str, true: str) -> float:
    if "-" not in inferred + true and len(inferred) == len(true):
        return identity(inferred, true)
    return aligned_identity(inferred, true)


def score(
    true: Truth,
    tree_newick: str,
    reconciliation: tuple[str, str],
    ancestors: Mapping[str, str] | None = None,
) -> Score:
    """One combination: its tree.nwk, its (reconciliation.tsv, reconciliation.nwk) and,
    optionally, its ancestors by node name."""
    tsv, rooted = reconciliation
    labels, _ = reconcile.rooted_clades(rooted)
    events = reconcile.event_names(tsv)
    inferred_dups = {labels[name] for name, event in events.items()
                     if event == "duplication" and name in labels}  # fmt: skip
    true_dups = {d.below for d in true.duplications}
    named = {name: below for name, below in labels.items() if name != reconcile.ROOT}
    matched = {name: true.clades[below] for name, below in named.items() if below in true.clades}
    identities = [
        sequence_identity(ancestors[name], true.ancestors[node])
        for name, node in matched.items()
        if ancestors and name in ancestors and node in true.ancestors
    ]
    losses = 0
    for line in tsv.splitlines()[1:]:
        cells = line.split("\t")
        if len(cells) == 5 and cells[3].isdigit():
            losses += int(cells[3])
    return Score(
        rf=robinson_foulds(true.newick, tree_newick),
        ancestors=len(named),
        duplications_true=len(true_dups),
        duplications_inferred=len(inferred_dups),
        duplications_found=len(true_dups & inferred_dups),
        chimeric=sorted(chimeric_nodes(named, true)),
        clades_matched=len(matched),
        identity=sum(identities) / len(identities) if identities else None,
        losses=losses,
    )


@dataclass(frozen=True)
class FlagScore:
    flagged: list[str]
    chimeric: list[str]

    @property
    def true_positives(self) -> int:
        return len(set(self.flagged) & set(self.chimeric))

    @property
    def recall(self) -> float | None:
        return self.true_positives / len(self.chimeric) if self.chimeric else None

    @property
    def precision(self) -> float | None:
        return self.true_positives / len(self.flagged) if self.flagged else None


def flags(true: Truth, own: tuple[str, str], reconciled: tuple[str, str]) -> FlagScore:
    """The paralog-mixing flag on every ancestor of one combination (its reconciliation
    files), against the reconciled tree of the same alignment and model."""
    labels, _ = reconcile.rooted_clades(own[1])
    named = {name: below for name, below in labels.items() if name != reconcile.ROOT}
    flagged = []
    for name in named:
        found = reconcile.paralogy(name, own, reconciled)
        if found is not None and found.mixes:
            flagged.append(name)
    return FlagScore(sorted(flagged), sorted(chimeric_nodes(named, true)))
