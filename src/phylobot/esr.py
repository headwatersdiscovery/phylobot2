"""Extant-sequence reconstruction (ESR; Sennett & Theobald 2024): hold out a tip,
reconstruct it as if it were an ancestor, and score the posterior against the real
sequence: a ground truth for ancestral posteriors on real data.

IQ-TREE reports posteriors for internal nodes only, so the held-out tip becomes
internal: two daughters on 1000-length branches, the tip's own sequence out of the
alignment, `-te -blfix` with every model parameter fixed, `--ancestral`. The daughters'
only residue sits in one extra column that every real sequence has as `?`: IQ-TREE
refuses all-`?` sequences, and S&T's poly-Leu is not neutral under +I (a Leu breaks
every constant column's invariant class). The extra column is exact under every model.

Levels:
- **e0**: S&T's sequence-wise ESR. The combination's tree and branch lengths, the model
  fitted on them. The tip informs all of it, its own pendant length included.
- **e1**: topology fixed; branch lengths and model refitted without the tip; the tip
  put back on the same split, at the same point, with its e0 pendant length (which the
  data without it cannot estimate).

Works from a combination's published files (alignment.fasta, tree.nwk, model.txt), or
from an uploaded alignment and tree (`fit_lengths`: the tree's branch lengths are fitted
first, as the pipeline's FIT_IQTREE fits a provided tree). A FreeRate (+R) model refits
the branch lengths too: IQ-TREE will not fit free +R parameters with -blfix. A site-profile
mixture (+C20, +C20+F) works in this engine, fixed from the fit's checkpoint, but `supported`
still refuses it, so no job offers one: only code calling the engine directly can. Amino-acid
combinations only for now. A fit is saved (`save`) so that each held-out tip can be its
own pipeline task (`phylobot.pipeline.esr`) without fitting the model again. IQ-TREE
re-roots and renames what it writes, so nodes are found by their relation to the
daughters and tips placed by split, never by name or position.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import subprocess
import tempfile
import time
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import dendropy

from .evaluation.metrics import LOG_FLOOR as LOG_FLOOR  # one floor for every log score
from .seqio import read_fasta, write_fasta
from .trees import load_tree

AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"
LEVELS = ("e0", "e1")
DUMMIES = ("esr_daughter_1", "esr_daughter_2")
DUMMY_LENGTH = 1000.0
# The fixed model must reproduce the fit it came from. The .iqtree report rounds every
# parameter to 4 decimals, and the cost grows with the alignment and the parameter count:
# about 0.05 under +R4 on Ime2 (spike), 0.20-0.22 under Q.PFAM+F+I+R7 on PINK1 (lnL
# about -77,900; 2026-11-01). So the allowance is 0.1 or 1e-5 of |lnL|, whichever is
# larger; a misread model string is off by tens to thousands.
LNL_TOLERANCE = 0.1
LNL_RELATIVE_TOLERANCE = 1e-5


def lnl_tolerance(log_likelihood: float) -> float:
    """How far a fixed model's lnL may be from the fit's: rounding, not a misread."""
    return max(LNL_TOLERANCE, LNL_RELATIVE_TOLERANCE * abs(log_likelihood))


TIP_FIELDS = ["tip", "level", "pendant", "sites", "accuracy", "mean_pp", "log_score",
              "log_likelihood", "seconds"]  # fmt: skip
SITE_FIELDS = ["tip", "level", "site", "true", "ml", "pp_ml", "pp_true", *AMINO_ACIDS]

# Runs IQ-TREE with these arguments in `work` and returns its stdout.
Runner = Callable[[Path, list[str]], str]


def docker_runner(image: str, seed: int = 1) -> Runner:
    def run(work: Path, args: list[str]) -> str:
        command = ["docker", "run", "--rm", "--platform", "linux/amd64"]
        command += ["-v", f"{work}:/w", "-w", "/w", image, "iqtree3", *args]
        command += ["-redo", "--quiet", "--seed", str(seed), "-T", "1"]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode:
            errors = [x for x in (result.stdout + result.stderr).splitlines() if "ERROR" in x]
            raise RuntimeError(f"IQ-TREE failed: {'; '.join(errors) or result.returncode}")
        return result.stdout

    return run


def local_runner(seed: int = 1) -> Runner:
    """IQ-TREE on this machine's PATH (inside the pipeline's container)."""

    def run(work: Path, args: list[str]) -> str:
        command = ["iqtree3", *args, "-redo", "--quiet", "--seed", str(seed), "-T", "1"]
        result = subprocess.run(command, cwd=work, capture_output=True, text=True, check=False)
        if result.returncode:
            errors = [x for x in (result.stdout + result.stderr).splitlines() if "ERROR" in x]
            raise RuntimeError(f"IQ-TREE failed: {'; '.join(errors) or result.returncode}")
        return result.stdout

    return run


# Model components `fixed_model` can fix; anything else (mixtures, +C20, +H) is refused
# before a job runs (`supported`). `fixed_mixture` fixes mixtures, for the benchmark only.
FIXABLE = re.compile(r"^(F|FO|FQ|I|G\d*|R\d+)$")
MIXTURE = re.compile(r"\+C\d+(\+|$)")


def is_mixture(model: str) -> bool:
    """Whether the model is a site-profile mixture (needs -mwopt, or the weights are not
    fitted: IQ-TREE 3.1.3 keeps the published ones)."""
    return MIXTURE.search(model) is not None


def free_rates(model: str) -> bool:
    """Whether the model has a FreeRate component (+R4, not yet fixed as +R4{...}), which
    IQ-TREE will not fit with branch lengths fixed (-blfix)."""
    return any(re.fullmatch(r"R\d+", part) for part in model.split("+")[1:])


def supported(model: str) -> bool:
    """Whether ESR can fix this model's parameters (an amino-acid matrix plus +F, +I,
    +G, +R)."""
    parts = model.split("+")
    return bool(parts[0]) and all(FIXABLE.match(p) for p in parts[1:])


# --- the model ------------------------------------------------------------------------


def _grab(pattern: str, text: str) -> str:
    m = re.search(pattern, text)
    if m is None:
        raise ValueError(f"not in the IQ-TREE report: {pattern}")
    return m.group(1)


def model_name(report: str) -> str:
    return _grab(r"Model of substitution: (\S+)", report)


def fixed_model(report: str) -> str:
    """The report's model with every parameter fixed: +F{20 frequencies}, +I{p},
    +G{alpha}, +R{w1,r1,...}. The exchangeabilities are the named matrix's."""
    name = model_name(report)
    parts = name.split("+")
    out = [parts[0]]
    for p in parts[1:]:
        if p.startswith("F"):
            freqs = [_grab(rf"pi\({a}\) = ([\d.]+)", report) for a in AMINO_ACIDS]
            out.append("F{" + ",".join(freqs) + "}")
        elif p.startswith("I"):
            out.append(p + "{" + _grab(r"Proportion of invariable sites: ([\d.]+)", report) + "}")
        elif p.startswith("G"):
            out.append(p + "{" + _grab(r"Gamma shape alpha: ([\d.]+)", report) + "}")
        elif p.startswith("R"):
            listed = _grab(r"Site proportion and rates:(.*)", report)
            pairs = re.findall(r"\(([\d.]+),([\d.]+)\)", listed)
            out.append(p + "{" + ",".join(f"{w},{r}" for w, r in pairs) + "}")
        else:
            raise ValueError(f"cannot fix model component +{p} of {name}")
    return "+".join(out)


def _floats(pattern: str, text: str) -> list[float]:
    return [float(x) for x in _grab(pattern, text).split(",")]


def fixed_mixture(model: str, checkpoint: str) -> str:
    """A site-profile mixture with every parameter fixed, for example
    MIX{LG+F{...}:1:w,LG+FC20pi1:1:w,...}+I{p}+G4{alpha}. From the checkpoint: the report
    prints too few digits to give back the fit's likelihood. A written-in +F class survives
    the held-out alignment, whose extra column would change an empirical one; IQ-TREE's
    report of such a model calls that class plain `LG` but honours its frequencies."""
    parts = model.split("+")
    matrix, rest = parts[0], parts[1:]
    profiles = next(p for p in rest if re.fullmatch(r"C\d+", p))
    weights = _floats(r"\n\s*prop: ([^\n]+)", checkpoint)
    classes = [f"{matrix}+F{profiles}pi{k}" for k in range(1, int(profiles[1:]) + 1)]
    if "F" in rest:
        pi = _floats(r"Component1!ModelProtein!state_freq: ([^\n]+)", checkpoint)
        classes.insert(0, f"{matrix}+F{{" + ",".join(f"{x:.10g}" for x in pi) + "}")
    if len(classes) != len(weights):
        raise ValueError(f"{len(weights)} weights in the checkpoint for {len(classes)} classes")
    mix = ",".join(f"{c}:1:{w:.10g}" for c, w in zip(classes, weights, strict=True))
    out = [f"MIX{{{mix}}}"]
    for p in rest:
        if p.startswith("I"):
            out.append(f"{p}{{{_grab(r'p_invar: ([0-9.e-]+)', checkpoint)}}}")
        elif p.startswith("G"):
            out.append(f"{p}{{{_grab(r'gamma_shape: ([0-9.e-]+)', checkpoint)}}}")
        elif p != "F" and p != profiles:
            raise ValueError(f"cannot fix model component +{p} of {model}")
    return "+".join(out)


def log_likelihood(report: str) -> float:
    return float(_grab(r"Log-likelihood of the tree: (-?[\d.]+)", report))


# --- trees ----------------------------------------------------------------------------


@dataclass(frozen=True)
class Placement:
    """Where a pruned tip sat: on the edge above `clade`, `fraction` of the edge's
    length from the clade's end, on a pendant branch of length `pendant`."""

    clade: frozenset[str]
    fraction: float
    pendant: float


def _newick(tree: dendropy.Tree) -> str:
    tree.seed_node.edge.length = None
    text: str = tree.as_string(
        schema="newick",
        suppress_rooting=True,
        suppress_internal_node_labels=True,
        unquoted_underscores=True,
    )
    return text.strip() + "\n"


def _unroot(tree: dendropy.Tree) -> None:
    """A two-child top node splits one edge in two; make it one edge again, so every
    edge is a node's edge whichever way the tree was rooted when written."""
    top = tree.seed_node
    children = top.child_nodes()
    if len(children) != 2:
        return
    keep, dissolve = sorted(children, key=lambda c: c.is_internal())
    if not dissolve.is_internal():
        return  # two tips: nothing to dissolve
    keep.edge.length = float(keep.edge.length) + float(dissolve.edge.length)
    top.remove_child(dissolve)
    for grandchild in dissolve.child_nodes():
        top.add_child(grandchild)


def _leaves(node: dendropy.Node) -> frozenset[str]:
    return frozenset(leaf.taxon.label for leaf in node.leaf_iter())


def _tip(tree: dendropy.Tree, tip: str) -> dendropy.Node:
    node = tree.find_node_with_taxon_label(tip)
    if node is None:
        raise ValueError(f"{tip!r} is not a tip of the tree")
    return node


def plain(newick: str) -> str:
    """The tree with internal labels dropped (IQ-TREE would read them as supports)."""
    return _newick(load_tree(newick))


def pendant_length(newick: str, tip: str) -> float:
    return float(_tip(load_tree(newick), tip).edge.length)


def held_out_tree(newick: str, tip: str) -> str:
    """The tree with `tip` made internal: two daughters on 1000-length branches."""
    tree = load_tree(newick)
    node = _tip(tree, tip)
    for name in DUMMIES:
        node.new_child(taxon=tree.taxon_namespace.require_taxon(label=name),
                       edge_length=DUMMY_LENGTH)  # fmt: skip
    node.taxon = None
    return _newick(tree)


def _fraction(part: float, total: float) -> float:
    """`part`'s share of an edge of length `total`; anywhere on a zero-length edge is the
    same place, so its middle (SR238: a tip whose parent and sibling edges were both 0,
    2026-11-02)."""
    return part / total if total > 0 else 0.5


def prune(newick: str, tip: str) -> tuple[str, Placement]:
    """The tree without `tip`, and where it sat. Binary trees with an optional
    trifurcating top node (as IQ-TREE writes them) only."""
    tree = load_tree(newick)
    leaf = _tip(tree, tip)
    parent = leaf.parent_node
    pendant = float(leaf.edge.length)
    parent.remove_child(leaf)
    children = parent.child_nodes()
    if parent is not tree.seed_node and len(children) == 1:
        # The parent is now of degree 2: merge its edge into its only child's.
        only = children[0]
        above, below = float(parent.edge.length), float(only.edge.length)
        placement = Placement(_leaves(only), _fraction(below, above + below), pendant)
        grandparent = parent.parent_node
        grandparent.remove_child(parent)
        grandparent.add_child(only)
        only.edge.length = above + below
    elif parent is tree.seed_node and len(children) == 2:
        # It hung off the top node; the two remaining edges become one.
        x, y = (float(c.edge.length) for c in children)
        placement = Placement(_leaves(children[0]), _fraction(x, x + y), pendant)
        _unroot(tree)
    else:
        raise ValueError(f"cannot hold out {tip!r}: the tree is not binary around it")
    return _newick(tree), placement


def reinsert(newick: str, tip: str, placement: Placement) -> str:
    """Put `tip` back on the edge its placement names, in a tree that may be rooted
    anywhere (IQ-TREE re-roots the trees it writes)."""
    tree = load_tree(newick)
    _unroot(tree)
    everyone = _leaves(tree.seed_node)
    fraction = placement.fraction
    below = None
    for want, f in ((placement.clade, fraction), (everyone - placement.clade, 1 - fraction)):
        below = next((n for n in tree.preorder_node_iter()
                      if n is not tree.seed_node and _leaves(n) == want), None)  # fmt: skip
        if below is not None:
            fraction = f
            break
    if below is None:
        raise ValueError(f"the tree has no edge for {tip!r}'s split")
    parent = below.parent_node
    length = float(below.edge.length)
    parent.remove_child(below)
    middle = parent.new_child(edge_length=length * (1 - fraction))
    middle.add_child(below)
    below.edge.length = length * fraction
    middle.new_child(taxon=tree.taxon_namespace.require_taxon(label=tip),
                     edge_length=placement.pendant)  # fmt: skip
    return _newick(tree)


def held_out_node(treefile: str) -> str:
    """IQ-TREE's name for the held-out node in a tree it wrote: the daughters' parent."""
    tree = dendropy.Tree.get(data=treefile, schema="newick", preserve_underscores=True)
    daughter = tree.find_node_with_taxon_label(DUMMIES[0])
    if daughter is None or daughter.parent_node is None:
        raise ValueError("no held-out daughter in IQ-TREE's tree")
    label: str = daughter.parent_node.label
    return label


# --- alignments and posteriors --------------------------------------------------------


def held_out_alignment(alignment: dict[str, str], tip: str) -> dict[str, str]:
    """Without `tip`, plus one column holding the daughters' only residue."""
    clash = set(DUMMIES) & set(alignment)
    if clash:
        raise ValueError(f"sequence names reserved for ESR: {sorted(clash)}")
    length = len(alignment[tip])
    out = {name: seq + "?" for name, seq in alignment.items() if name != tip}
    for name in DUMMIES:
        out[name] = "?" * length + "L"
    return out


def read_state(text: str, node: str, length: int) -> list[dict[str, float]]:
    """The node's posterior at columns 1..length of an IQ-TREE .state file (the extra
    column's row dropped)."""
    rows: dict[int, dict[str, float]] = {}
    header: list[str] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "Node":
            header = [h.removeprefix("p_") for h in fields[3:]]
            continue
        if fields[0] == node and int(fields[1]) <= length:
            rows[int(fields[1])] = dict(zip(header, map(float, fields[3:]), strict=True))
    if sorted(rows) != list(range(1, length + 1)):
        raise ValueError(f"{node}: expected {length} columns, found {len(rows)}")
    return [rows[i] for i in range(1, length + 1)]


@dataclass
class TipResult:
    tip: str
    level: str
    pendant: float
    log_likelihood: float
    seconds: float
    sites: list[dict[str, object]]

    def summary(self) -> dict[str, object]:
        n = len(self.sites)
        correct = sum(s["true"] == s["ml"] for s in self.sites)
        pp_ml = sum(float(str(s["pp_ml"])) for s in self.sites)
        logs = sum(math.log(max(float(str(s["pp_true"])), LOG_FLOOR)) for s in self.sites)
        return {
            "tip": self.tip, "level": self.level, "pendant": f"{self.pendant:.6g}", "sites": n,
            "accuracy": f"{correct / n:.4f}" if n else "",
            "mean_pp": f"{pp_ml / n:.4f}" if n else "",
            "log_score": f"{logs / n:.4f}" if n else "",
            "log_likelihood": f"{self.log_likelihood:.4f}", "seconds": f"{self.seconds:.1f}",
        }  # fmt: skip


def score(
    tip: str, level: str, truth: str, posterior: Sequence[dict[str, float]]
) -> list[dict[str, object]]:
    """One record per column where the held-out sequence has an amino acid."""
    records: list[dict[str, object]] = []
    for site, (true, probs) in enumerate(zip(truth, posterior, strict=True), start=1):
        if true not in AMINO_ACIDS:
            continue
        ml = max(AMINO_ACIDS, key=lambda a: probs[a])
        records.append({"tip": tip, "level": level, "site": site, "true": true, "ml": ml,
                        "pp_ml": probs[ml], "pp_true": probs[true],
                        **{a: probs[a] for a in AMINO_ACIDS}})  # fmt: skip
    return records


# --- running --------------------------------------------------------------------------


class Combination:
    """One combination's published files, fitted once with its branch lengths fixed so
    every held-out tip reuses the same parameters."""

    def __init__(self, folder: Path, work: Path, run: Runner, fit_lengths: bool = False):
        self.work, self.run = work, run
        self.alignment = read_fasta(folder / "alignment.fasta")
        residues = set("".join(self.alignment.values()).upper()) - set("-?.")
        if residues <= set("ACGTUN"):
            raise ValueError("ESR handles amino-acid combinations only (codon: not yet)")
        self.model = (folder / "model.txt").read_text().strip()
        self.tree = plain((folder / "tree.nwk").read_text())
        write_fasta(work / "full.fasta", self.alignment)
        (work / "full.nwk").write_text(self.tree)
        # IQ-TREE refuses -blfix with a free FreeRate model (its rates and the branch
        # lengths share a scale), so a +R model refits the lengths on the fixed topology,
        # as an uploaded tree does. Its fixed parameters take -blfix (the check below).
        self.lengths_refitted = fit_lengths or free_rates(self.model)
        lengths = [] if self.lengths_refitted else ["-blfix"]
        self.run(work, ["-s", "full.fasta", "-te", "full.nwk", *lengths, "-m", self.model,
                        *self.weights, "--prefix", "full"])  # fmt: skip
        if self.lengths_refitted:
            self.tree = plain((work / "full.treefile").read_text())
            (work / "full.nwk").write_text(self.tree)
        report = (work / "full.iqtree").read_text()
        self.fixed = self.fixed_from("full")
        self.run(work, ["-s", "full.fasta", "-te", "full.nwk", "-blfix", "-m", self.fixed,
                        "--prefix", "check"])  # fmt: skip
        fitted, fixed = log_likelihood(report), log_likelihood((work / "check.iqtree").read_text())
        if abs(fitted - fixed) > lnl_tolerance(fitted):
            raise ValueError(f"fixed model {self.fixed} gives lnL {fixed}, the fit {fitted}")
        self.log_likelihood = fitted

    @property
    def weights(self) -> list[str]:
        return ["-mwopt"] if is_mixture(self.model) else []

    def fixed_from(self, prefix: str) -> str:
        """The model this fit gave, every parameter fixed."""
        if is_mixture(self.model):
            with gzip.open(self.work / f"{prefix}.ckp.gz", "rt") as f:
                return fixed_mixture(self.model, f.read())
        return fixed_model((self.work / f"{prefix}.iqtree").read_text())

    def save(self, folder: Path) -> None:
        """What `load` needs: the alignment, the tree (with its fitted lengths) and the
        model, named and fixed."""
        folder.mkdir(parents=True, exist_ok=True)
        write_fasta(folder / "alignment.fasta", self.alignment)
        (folder / "tree.nwk").write_text(self.tree)
        fit = {"model": self.model, "fixed": self.fixed, "log_likelihood": self.log_likelihood,
               "lengths_refitted": self.lengths_refitted}  # fmt: skip
        (folder / "fit.json").write_text(json.dumps(fit, indent=2) + "\n")

    @classmethod
    def load(cls, folder: Path, work: Path, run: Runner) -> Combination:
        """A saved fit, ready for held-out tips (no IQ-TREE run)."""
        self = cls.__new__(cls)
        self.work, self.run = work, run
        self.alignment = read_fasta(folder / "alignment.fasta")
        self.tree = (folder / "tree.nwk").read_text()
        fit = json.loads((folder / "fit.json").read_text())
        self.model, self.fixed = fit["model"], fit["fixed"]
        self.log_likelihood = float(fit["log_likelihood"])
        self.lengths_refitted = bool(fit.get("lengths_refitted"))
        return self

    def reconstruct(self, prefix: str, tip: str, tree: str, model: str, level: str,
                    pendant: float) -> TipResult:  # fmt: skip
        start = time.monotonic()
        write_fasta(self.work / f"{prefix}.fasta", held_out_alignment(self.alignment, tip))
        (self.work / f"{prefix}.nwk").write_text(held_out_tree(tree, tip))
        self.run(self.work, ["-s", f"{prefix}.fasta", "-te", f"{prefix}.nwk", "-blfix",
                             "-m", model, "--ancestral", "--prefix", prefix])  # fmt: skip
        node = held_out_node((self.work / f"{prefix}.treefile").read_text())
        truth = self.alignment[tip]
        posterior = read_state((self.work / f"{prefix}.state").read_text(), node, len(truth))
        lnl = log_likelihood((self.work / f"{prefix}.iqtree").read_text())
        sites = score(tip, level, truth, posterior)
        return TipResult(tip, level, pendant, lnl, time.monotonic() - start, sites)

    def e0(self, i: int, tip: str) -> TipResult:
        pendant = pendant_length(self.tree, tip)
        return self.reconstruct(f"e0_{i}", tip, self.tree, self.fixed, "e0", pendant)

    def e1(self, i: int, tip: str) -> TipResult:
        start = time.monotonic()
        pruned, placement = prune(self.tree, tip)
        prefix = f"e1fit_{i}"
        write_fasta(self.work / f"{prefix}.fasta",
                    {k: v for k, v in self.alignment.items() if k != tip})  # fmt: skip
        (self.work / f"{prefix}.nwk").write_text(pruned)
        self.run(self.work, ["-s", f"{prefix}.fasta", "-te", f"{prefix}.nwk", "-m", self.model,
                             *self.weights, "--prefix", prefix])  # fmt: skip
        refit = self.fixed_from(prefix)
        tree = reinsert((self.work / f"{prefix}.treefile").read_text(), tip, placement)
        result = self.reconstruct(f"e1_{i}", tip, tree, refit, "e1", placement.pendant)
        result.seconds = time.monotonic() - start
        return result


def write_tsv(path: Path, fields: list[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_combination(
    folder: Path,
    out: Path,
    run: Runner,
    levels: Sequence[str] = LEVELS,
    tips: Sequence[str] | None = None,
) -> list[TipResult]:
    """ESR on one combination: esr_tips.tsv (one row per tip and level) and
    esr_sites.tsv (one row per scored column) in `out`."""
    if unknown := set(levels) - set(LEVELS):
        raise ValueError(f"unknown ESR levels: {sorted(unknown)}")
    out = out.resolve()  # mounted into the container
    out.mkdir(parents=True, exist_ok=True)
    results: list[TipResult] = []
    with tempfile.TemporaryDirectory(prefix="esr-", dir=out) as tmp:
        combination = Combination(folder, Path(tmp), run)
        chosen = list(combination.alignment) if tips is None else list(tips)
        for i, tip in enumerate(chosen, start=1):
            if tip not in combination.alignment:
                raise ValueError(f"{tip!r} is not in the alignment")
            for level in levels:
                results.append(getattr(combination, level)(i, tip))
    write_tsv(out / "esr_tips.tsv", TIP_FIELDS, (r.summary() for r in results))
    write_tsv(out / "esr_sites.tsv", SITE_FIELDS, (s for r in results for s in r.sites))
    return results


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("combination", type=Path, help="a published combination folder")
    parser.add_argument("out", type=Path)
    parser.add_argument("--image", required=True, help="the pinned IQ-TREE image")
    parser.add_argument("--levels", default=",".join(LEVELS))
    parser.add_argument("--tips", help="comma-separated; default every tip")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args(argv)
    tips = args.tips.split(",") if args.tips else None
    run_combination(args.combination, args.out, docker_runner(args.image, args.seed),
                    args.levels.split(","), tips)  # fmt: skip


if __name__ == "__main__":
    main()
