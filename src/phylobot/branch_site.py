"""The branch-site test of positive selection (codeml Model A against its null; Zhang,
Nielsen & Yang 2005), on foreground branches the user marks.

Model A (model = 2, NSsites = 2) has four site classes: 0 (dN/dS < 1 everywhere), 1
(dN/dS = 1 everywhere), 2a and 2b (dN/dS < 1 or = 1 on the background, ω2 >= 1 on the
foreground). The null fixes ω2 = 1. Twice the difference in log-likelihood is tested
against χ² with 1 degree of freedom (conservative; the null distribution is a 50:50
mixture of 0 and χ²₁, whose p-value is half as large). Sites with a high posterior of
being in class 2a or 2b on the foreground come from codeml's Bayes empirical Bayes.

Foreground branches are written with codeml's own marks: `#1` after a node marks the
branch leading to it, `$1` every branch of its clade (the branch above it included).
Branches are identified here by their split (the leaves on one side), so a mark made on
an outgroup-rooted tree carries over to the unrooted tree codeml analyses.

Used by VALIDATE (`encode_marks`, `strip_mark`, `splits_of`, `mark_tree`,
`write_controls`), BRANCH_SITE_EXPORT (`parse_mlc`,
`test`) and the tree viewer (`foreground_splits`, `mark_tree`).
"""

import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import dendropy

from phylobot import trees
from phylobot.codeml import CONTROL, FLOAT, LNL, floats
from phylobot.pipeline.export import quote

Split = frozenset[str]  # a branch: the leaves on one side (trees.canonical_split)

ALTERNATIVE, NULL = "A", "A_null"
# Model A from two starting values of ω2 (it has local optima, as M2a and M8 do); the
# null has ω2 fixed at 1.
STARTS = {ALTERNATIVE: (1.5, 4.0), NULL: (1.0,)}
SCOPES = {
    "branch": "the branch leading to it",
    "clade": "every branch in its clade (the branch leading to it included)",
}
MARK = {"#": "branch", "$": "clade"}
SUFFIX = {"branch": "__phylobot_fg_branch__", "clade": "__phylobot_fg_clade__"}


class MarkError(ValueError):
    pass


def name(node: dendropy.Node) -> str:
    """A node's name: its taxon's for a parsed leaf, else its label."""
    return str(node.taxon.label if node.taxon is not None else node.label or "")


def leaves_below(node: dendropy.Node) -> frozenset[str]:
    return frozenset(name(leaf) for leaf in node.leaf_iter())


def branch_splits(tree: dendropy.Tree) -> dict[Split, dendropy.Node]:
    """Every branch of the tree as an unrooted split, with a node it leads to (the root's
    two branches of a rooted tree are one unrooted branch)."""
    everyone = leaves_below(tree.seed_node)
    out: dict[Split, dendropy.Node] = {}
    for node in tree.preorder_node_iter():
        if node is not tree.seed_node:
            out.setdefault(trees.canonical_split(leaves_below(node), everyone), node)
    return out


def node_named(tree: dendropy.Tree, label: str) -> dendropy.Node | None:
    for node in tree.preorder_node_iter():
        if name(node) == label:
            return node
    return None


def splits_of(tree: dendropy.Tree, node: dendropy.Node, scope: str) -> set[Split]:
    """The foreground branches a mark on `node` means: its own branch, or its clade's."""
    if node is tree.seed_node and scope == "branch":
        raise MarkError("the tree's top node has no branch leading to it")
    everyone = leaves_below(tree.seed_node)
    nodes = [node] if scope == "branch" else list(node.preorder_iter())
    return {
        trees.canonical_split(leaves_below(n), everyone) for n in nodes if n is not tree.seed_node
    }


def check_background(foreground: set[Split], every: set[Split]) -> None:
    if not foreground:
        raise MarkError("no foreground branch is marked")
    if foreground >= every:
        raise MarkError("every branch is foreground: the test needs background branches too")


def foreground_splits(newick: str, node: str, scope: str) -> set[Split]:
    """The foreground chosen in the tree viewer: `node` (an ancestor's label or an extant
    sequence's name) on the tree as drawn (outgroup-rooted or not), and the scope."""
    if scope not in SCOPES:
        raise MarkError(f"unknown scope {scope!r}")
    tree = trees.load_tree(newick)
    found = node_named(tree, node)
    if found is None:
        raise MarkError(f"no node {node!r} in the tree")
    foreground = splits_of(tree, found, scope)
    check_background(foreground, set(branch_splits(tree)))
    return foreground


def mark_tree(newick: str, foreground: set[Split]) -> str:
    """The tree for codeml: unrooted, leaf names and branch lengths kept, no internal
    labels, ` #1` after each foreground branch's node."""
    tree = trees.load_tree(newick)
    if len(tree.seed_node.child_nodes()) == 2:
        tree.collapse_basal_bifurcation()  # codeml analyses an unrooted tree
    by_split = branch_splits(tree)
    missing = foreground - set(by_split)
    if missing:
        raise MarkError(f"{len(missing)} foreground branches are not in the tree")
    check_background(foreground, set(by_split))
    marked = {id(by_split[s]) for s in foreground}

    def write(node: dendropy.Node) -> str:
        if node.is_leaf():
            text = quote(name(node))
        else:
            text = "(" + ",".join(write(child) for child in node.child_nodes()) + ")"
        if id(node) in marked:
            text += " #1"
        if node is not tree.seed_node and node.edge.length is not None:
            text += f":{node.edge.length:.10g}"
        return text

    return write(tree.seed_node) + ";"


def encode_marks(text: str) -> str:
    """A user's Newick with codeml marks (`#1`, `$1`) -> the same tree with each mark
    folded into the label before it (a leaf's name, an internal node's label, often
    empty), so a Newick parser keeps it; `strip_mark` takes it out again."""
    for suffix in SUFFIX.values():
        if suffix in text:
            raise MarkError("the tree contains a reserved name")
    if any(n != "1" for n in re.findall(r"[#$](\d+)", text)):
        raise MarkError("only #1 and $1 marks are supported (one foreground class)")
    if re.search(r"[#$](?!1)", text):
        raise MarkError("'#' and '$' are only allowed in codeml's #1 and $1 marks")
    # After a quoted label the suffix goes inside the quotes.
    return re.sub(
        r"(')?\s*([#$])1",
        lambda m: SUFFIX[MARK[m.group(2)]] + (m.group(1) or ""),
        text,
    )


def strip_mark(label: str | None) -> tuple[str | None, str | None]:
    """A label as parsed -> (the label without its mark, the mark's scope or None)."""
    if label is None:
        return None, None
    for scope, suffix in SUFFIX.items():
        if label.endswith(suffix):
            return label[: -len(suffix)] or None, scope
    return label, None


@dataclass
class SiteClasses:
    """Model A's site classes (0, 1, 2a, 2b): proportions and dN/dS on the background
    and the foreground."""

    proportions: list[float]
    background: list[float]
    foreground: list[float]


@dataclass
class Fit:
    lnl: float
    parameters: int
    kappa: float
    classes: SiteClasses
    beb: list[tuple[int, str, float]]  # site, residue, P(foreground ω > 1)


def parse_mlc(text: str) -> Fit:
    """One Model A (or null) run's estimates from mlc."""
    match = LNL.search(text)
    if match is None:
        raise ValueError("codeml output has no lnL line")
    kappa = re.search(r"^kappa \(ts/tv\) =\s*" + FLOAT, text, re.MULTILINE)
    table = re.search(
        r"^proportion\s+(.*)\n^background w\s+(.*)\n^foreground w\s+(.*)$", text, re.MULTILINE
    )
    if kappa is None or table is None:
        raise ValueError("codeml output is not from the branch-site model")
    beb = []
    if "Bayes Empirical Bayes (BEB)" in text:
        block = text.split("Bayes Empirical Bayes (BEB)", 1)[1].split("The grid", 1)[0]
        beb = [
            (int(m.group(1)), m.group(2), float(m.group(3)))
            for m in re.finditer(r"^\s*(\d+) (\S)\s+([\d.]+)\**\s*$", block, re.MULTILINE)
        ]
    return Fit(
        lnl=float(match.group(3)),
        parameters=int(match.group(2)),
        kappa=float(kappa.group(1)),
        classes=SiteClasses(*(floats(table.group(i)) for i in (1, 2, 3))),
        beb=beb,
    )


SITE_ROW = re.compile(r"^\s*(\d+) (\S)\s+" + r"\s+".join([FLOAT] * 4) + r"\s+\(\s*\d+\)")


def parse_site_classes(rst: str) -> list[tuple[int, str, tuple[float, float, float, float]]]:
    """Model A's BEB posteriors of the four site classes (0, 1, 2a, 2b) at every site,
    from rst (written even without ancestors). The residue is the first sequence's."""
    if "Bayes Empirical Bayes (BEB) probabilities for 4 classes" not in rst:
        raise ValueError("codeml rst has no BEB site-class probabilities")
    block = rst.split("Bayes Empirical Bayes (BEB) probabilities for 4 classes", 1)[1]
    sites = []
    for line in block.splitlines()[1:]:
        if match := SITE_ROW.match(line):
            site, residue, *probabilities = match.groups()
            p0, p1, p2a, p2b = map(float, probabilities)
            sites.append((int(site), residue, (p0, p1, p2a, p2b)))
        elif sites and not line.strip():
            break  # the block ends at its first blank line after the rows
    if not sites:
        raise ValueError("codeml's BEB site-class table is empty")
    return sites


@dataclass(frozen=True)
class Test:
    statistic: float  # 2 (lnL_A - lnL_null), 0 if negative (optimisation noise)
    p_value: float  # χ² with 1 degree of freedom
    p_mixture: float  # the 50:50 mixture of 0 and χ²₁


def test(alternative: float, null: float) -> Test:
    statistic = max(2 * (alternative - null), 0.0)
    p = math.erfc(math.sqrt(statistic / 2))  # P(χ²₁ >= statistic)
    return Test(statistic, p, p / 2 if statistic > 0 else 1.0)


def write_controls(outdir: Path, icode: int) -> None:
    """codeml/A/start<n>.ctl and codeml/A_null/start1.ctl (VALIDATE's control files, as
    for the site models: CODEML runs each)."""
    for key, starts in STARTS.items():
        folder = outdir / key
        folder.mkdir(parents=True, exist_ok=True)
        for n, omega in enumerate(starts, 1):
            text = CONTROL.format(model=2, nssites=2, icode=icode, omega=omega, ncatg=1)
            fixed = "1" if key == NULL else "0"
            text = text.replace("fix_omega = 0", f"fix_omega = {fixed}")
            # Ancestors are not wanted here; CODEML still expects an rst file.
            text = text.replace("RateAncestor = 1", "RateAncestor = 0")
            (folder / f"start{n}.ctl").write_text(text)


def tree_key(alignment: bytes, newick: str, genetic_code: int) -> str:
    """What makes two trees the same for the branch-site test: the codon alignment, the
    genetic code and the unrooted topology. Branch lengths and the tree's model are left
    out: codeml re-estimates every length, and its starting values do not change the
    optimum (lysin, Node10: the same lnL to 6 decimals from the codon tree's and the
    amino-acid tree's lengths, 2026-10-02). Node names and rooting are left out too."""
    topology = sorted(sorted(split) for split in branch_splits(trees.load_tree(newick)))
    identity = {
        "alignment": hashlib.sha256(alignment).hexdigest(),
        "genetic_code": int(genetic_code),
        "topology": topology,
    }
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()


def test_key(tree: str, foreground: set[Split]) -> str:
    """The same analysis: the same tree (`tree_key`) and the same foreground branches."""
    identity = {"tree": tree, "foreground": sorted(sorted(split) for split in foreground)}
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()


def describe(of: dict[str, Any]) -> str:  # the job's params["branch_site_of"]
    """'the branch leading to Node10', 'every branch in Node10's clade'."""
    node = of.get("node", "?")
    if of.get("scope") == "clade":
        return f"every branch in {node}'s clade (the branch leading to it included)"
    return f"the branch leading to {node}"
