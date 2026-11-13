"""Convert IQ-TREE output into PhyloBot's result files, restoring original names.

    phylobot-export NAMES.tsv ALIGNMENT.fasta ASR_TREE SUPPORT_TREE STATEFILE MODEL.txt \
        REPORT OUTDIR

ASR_TREE has IQ-TREE node names (Node1...) but no support values; SUPPORT_TREE
has branch supports but no names. Supports are matched to names by split.
--supports lists the measures in SUPPORT_TREE's node labels, in IQ-TREE's order
(sh_alrt,ufboot when both ran: labels 'SH-aLRT/UFBoot'); supports.tsv has one column
per measure (ufboot, sh_alrt), blank for a measure that did not run.
REPORT is IQ-TREE's .iqtree file from the tree step (search or fit); for a ModelFinder
combination, --modelfinder-report is MODELFINDER's (the candidate ranking).
Optional gap placement (see phylobot.gaps): --gap-methods parsimony,ml,prank,grasp
writes ancestors.gaps-<method>.fasta for each method (ancestors.fasta stays unmasked).
ml needs --gap-state/--gap-tree (IQ-TREE's reconstruction of the gap/residue alignment)
and writes gaps-ml.tsv; prank and grasp need --<method>-ancestors/--<method>-tree
(that program's ancestors and node-labelled tree). With an outgroup, --rooted-tree and
--rooting (from ROOT_TREE) are published as tree.rooted.nwk and rooting.json, and Fitch
uses the rooted tree (otherwise the tree's top node).
For a codon model (coding DNA), --codon-alignment and --genetic-code: STATEFILE holds
codon posteriors, summed into amino-acid ones (codons.amino_posteriors) for every file
above, so the viewers read codon and amino-acid combinations alike; also written:
alignment.cds.fasta, ancestors.cds.fasta (each site's best codon for its best amino
acid), ancestors_codon_pp.tsv (that codon's posterior), ancestors_codon_posteriors.tsv
(every codon with posterior 0.01 or more: node, site, codon, p) and, per masking gap method,
ancestors.gaps-<method>.cds.fasta. ALIGNMENT is then the protein alignment (one column
per codon).
Writes to OUTDIR: alignment.fasta, tree.nwk, supports.tsv, ancestors.fasta,
ancestors_pp.tsv, ancestors_posteriors.tsv, model.txt, fit.tsv (final-tree fit
statistics), and modelfinder.tsv (ModelFinder's candidate ranking, if it ran).
--site-rates and --site-lh (IQ-TREE's -wsr and -wsl, from the reconstruction on the
final tree) give site_lh.tsv (each column's log-likelihood) and, for a model with rate
heterogeneity (+G, +R; otherwise IQ-TREE lists none), site_rates.tsv (posterior mean
rate, most probable category: 0 invariable, 1 slowest..., and that category's rate).

Jobs that stop early (stop_after) use the smaller exports:

    phylobot-export-alignment NAMES.tsv ALIGNMENT.fasta OUTDIR [CODON_ALIGNMENT.fasta]
    phylobot-export-tree NAMES.tsv ALIGNMENT.fasta TREE MODEL.txt REPORT OUTGROUP.txt OUTDIR \
        [--supports ...] [--modelfinder-report ...]

export-tree names the tree's internal nodes Node1, Node2, ... (preorder), as the
reconstruction would, so trees look the same with or without ancestors; supports go
to supports.tsv. With an outgroup (OUTGROUP.txt non-empty) it also roots the tree.
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import dendropy

from phylobot import ancestors, codons, gaps, rooting
from phylobot.seqio import read_fasta, write_fasta
from phylobot.trees import canonical_split, internal_splits, load_tree

NEWICK_SPECIAL = set("()[]':;, \t")
SAFE_ID = re.compile(r"\bs\d{4,}\b")


def quote(name: str) -> str:
    """Quote a Newick label if it contains characters with syntactic meaning."""
    if any(c in NEWICK_SPECIAL for c in name):
        return "'" + name.replace("'", "''") + "'"
    return name


# fit.tsv columns: statistics of the final ML tree under the combination's model.
FIT_FIELDS = (
    "model", "log_likelihood", "free_parameters", "aic", "aicc", "bic", "sites", "tree_length",
)  # fmt: skip
FIT_PATTERNS = {
    "model": r"^Model of substitution: (\S+)",
    "log_likelihood": r"^Log-likelihood of the tree: (-?[\d.]+)",
    "free_parameters": r"^Number of free parameters.*: (\d+)",
    "aic": r"^Akaike information criterion \(AIC\) score: (-?[\d.]+)",
    "aicc": r"^Corrected Akaike information criterion \(AICc\) score: (-?[\d.]+)",
    "bic": r"^Bayesian information criterion \(BIC\) score: (-?[\d.]+)",
    "sites": r"^Input data: \d+ sequences with (\d+) ",
    "tree_length": r"^Total tree length \(sum of branch lengths\): ([\d.]+)",
}
MODELFINDER_FIELDS = (
    "model", "log_likelihood", "aic", "aic_in_95", "w_aic", "aicc", "aicc_in_95", "w_aicc",
    "bic", "bic_in_95", "w_bic",
)  # fmt: skip


def parse_fit(report: str, mixture: str = "") -> dict[str, str]:
    """Final-tree statistics from an IQ-TREE .iqtree report (the last match wins: in a
    ModelFinder report the final tree's section follows the model selection). A mixture's
    report says "Mixture model of substitution: MIX{...}" with every class listed: its row
    carries `mixture`, the catalog name, instead."""
    fit: dict[str, str] = {}
    for field, pattern in FIT_PATTERNS.items():
        matches = re.findall(pattern, report, re.MULTILINE)
        if field == "model" and not matches and mixture:
            matches = re.findall(r"^Mixture model of substitution: \S+", report, re.MULTILINE)
            matches = [mixture] if matches else []
        if not matches:
            raise ValueError(f"IQ-TREE report has no {field}")
        fit[field] = matches[-1]
    return fit


def parse_modelfinder(report: str) -> list[dict[str, str]]:
    """ModelFinder's table ('List of models sorted by BIC scores'); empty if it did not run.
    Scores come from ModelFinder's initial tree, not the final one; '+' marks models in
    the 95% confidence set for that criterion."""
    lines = report.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.startswith("List of models sorted"))
    except StopIteration:
        return []
    rows = []
    for line in lines[start + 3 :]:
        parts = line.split()
        if len(parts) != len(MODELFINDER_FIELDS):
            break
        row = dict(zip(MODELFINDER_FIELDS, parts, strict=True))
        for flag in ("aic_in_95", "aicc_in_95", "bic_in_95"):
            row[flag] = "yes" if row[flag] == "+" else "no"
        rows.append(row)
    return rows


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    lines = ["\t".join(fields)] + ["\t".join(row[f] for f in fields) for row in rows]
    path.write_text("\n".join(lines) + "\n")


def read_names(path: Path) -> dict[str, str]:
    rows = path.read_text().splitlines()[1:]
    return {row.split("\t")[0]: row.split("\t")[1] for row in rows}


def supports_by_node(asr_newick: str, support_newick: str) -> list[tuple[str, str]]:
    """Match each named node of the ASR tree to its support value, by split."""
    taxa = dendropy.TaxonNamespace()
    named = internal_splits(load_tree(asr_newick, taxa))
    supported = internal_splits(load_tree(support_newick, taxa))
    everyone = frozenset(taxon.label for taxon in taxa)
    by_split = {
        canonical_split(leaves, everyone): label
        for label, leaves in supported
        if leaves != everyone
    }
    return [(node, by_split.get(canonical_split(leaves, everyone), "")) for node, leaves in named]


# supports.tsv columns, after "node". Node labels list the measures that ran in IQ-TREE's
# order (SH-aLRT first), so the labels are split by what the tree step was asked for.
SUPPORT_COLUMNS = ("ufboot", "sh_alrt")


def index_by_node(source: Path) -> int:
    """Write the byte offsets of each node's rows beside a posteriors table
    (<name>.index.tsv), so the viewers read one node with a ranged read instead of the
    whole table, which at the pilot limits is hundreds of megabytes per combination.
    Returns how many nodes."""
    index = ancestors.index_table(source)
    (source.parent / ancestors.index_name(source.name)).write_text(ancestors.write_index(index))
    return len(index.nodes)


def write_supports(path: Path, supports: list[tuple[str, str]], kinds: tuple[str, ...]) -> None:
    unknown = set(kinds) - set(SUPPORT_COLUMNS)
    if unknown:
        raise ValueError(f"unknown support measures {sorted(unknown)}")
    rows = ["\t".join(("node", *SUPPORT_COLUMNS))]
    for node, label in supports:
        values = label.split("/") if label else [""] * len(kinds)
        if len(values) != len(kinds):
            raise ValueError(f"support label {label!r} at {node} does not match {kinds}")
        by_kind = dict(zip(kinds, values, strict=True))
        rows.append("\t".join((node, *(by_kind.get(c, "") for c in SUPPORT_COLUMNS))))
    path.write_text("\n".join(rows) + "\n")


def write_alignment(
    names: dict[str, str], alignment: Path, outdir: Path, filename: str = "alignment.fasta"
) -> None:
    with open(outdir / filename, "w") as handle:
        for safe, seq in read_fasta(alignment).items():
            handle.write(f">{names[safe]}\n{seq}\n")


def restorer(names: dict[str, str]) -> Any:
    """Newick text with safe IDs swapped for quoted original names."""
    return lambda text: SAFE_ID.sub(lambda m: quote(names[m.group(0)]), text)


def write_model_and_fit(
    model_path: Path, report: Path, outdir: Path, ranking_report: Path | None = None
) -> None:
    """model.txt, fit.tsv (the tree step's final fit) and, for a ModelFinder combination,
    modelfinder.tsv (the ranking in MODELFINDER's report)."""
    model = model_path.read_text().strip()
    (outdir / "model.txt").write_text(model + "\n")
    write_tsv(outdir / "fit.tsv", FIT_FIELDS, [parse_fit(report.read_text(), model)])
    if ranking_report and ranking_report.suffix == ".tsv":
        # A pick an extension reused from the job it extends (main.nf, modelfinder_reused):
        # that job's modelfinder.tsv, already a ranking, kept as it is.
        rows = ranking_report.read_text().splitlines()
        if rows and rows[0].split("\t") == list(MODELFINDER_FIELDS):
            (outdir / "modelfinder.tsv").write_text("\n".join(rows) + "\n")
    elif ranking_report and (ranking := parse_modelfinder(ranking_report.read_text())):
        write_tsv(outdir / "modelfinder.tsv", MODELFINDER_FIELDS, ranking)


def write_rooted(
    restore: Any, names: dict[str, str], rooted_newick: str, report: dict[str, Any], outdir: Path
) -> None:
    (outdir / "tree.rooted.nwk").write_text(restore(rooted_newick) + "\n")
    for key in ("outgroup", "intruders"):
        report[key] = [names[safe] for safe in report[key]]
    (outdir / "rooting.json").write_text(json.dumps(report, indent=2) + "\n")


def export_alignment(
    names_path: Path, alignment: Path, outdir: Path, codon_alignment: Path | None = None
) -> None:
    """An alignment with its original names (every job; alignment-only jobs end here),
    and for coding DNA its codon form, alignment.cds.fasta."""
    outdir.mkdir(parents=True, exist_ok=True)
    names = read_names(names_path)
    write_alignment(names, alignment, outdir)
    if codon_alignment is not None:
        write_alignment(names, codon_alignment, outdir, "alignment.cds.fasta")


def name_nodes(newick: str) -> str:
    """The tree with internal nodes named Node1, Node2, ... in preorder (support values,
    which IQ-TREE writes as internal labels, are dropped: they go to supports.tsv)."""
    tree = load_tree(newick)
    for index, node in enumerate(tree.preorder_internal_node_iter(), 1):
        node.label = f"Node{index}"
    named: str = tree.as_string(
        schema="newick", suppress_rooting=True, unquoted_underscores=True
    ).strip()
    return named


def export_tree(
    names_path: Path,
    alignment: Path,
    tree: Path,
    model_path: Path,
    report: Path,
    outgroup_path: Path,
    outdir: Path,
    supports: tuple[str, ...] = ("ufboot",),
    ranking_report: Path | None = None,
) -> None:
    """A combination's tree without ancestors (jobs that stop after trees)."""
    names = read_names(names_path)
    outdir.mkdir(parents=True, exist_ok=True)
    write_alignment(names, alignment, outdir)
    restore = restorer(names)
    supported = tree.read_text().strip()
    named = name_nodes(supported)
    (outdir / "tree.nwk").write_text(restore(named) + "\n")
    write_supports(outdir / "supports.tsv", supports_by_node(named, supported), supports)
    outgroup = {line.strip() for line in outgroup_path.read_text().splitlines() if line.strip()}
    if outgroup:
        rooted, rooting_report = rooting.root_on_outgroup(named, outgroup)
        write_rooted(restore, names, rooted, rooting_report, outdir)
    write_model_and_fit(model_path, report, outdir, ranking_report)


def export(
    names_path: Path,
    alignment: Path,
    asr_tree: Path,
    support_tree: Path,
    statefile: Path,
    model_path: Path,
    report: Path,
    outdir: Path,
    gap_methods: tuple[str, ...] = (),
    gap_state: Path | None = None,
    gap_tree: Path | None = None,
    external: dict[str, tuple[Path, Path]] | None = None,
    rooted_tree: Path | None = None,
    rooting_json: Path | None = None,
    support_kinds: tuple[str, ...] = ("ufboot",),
    ranking_report: Path | None = None,
    codon_alignment: Path | None = None,
    genetic_code: int = 1,
) -> None:
    names = read_names(names_path)
    outdir.mkdir(parents=True, exist_ok=True)
    write_alignment(names, alignment, outdir)
    code = codons.genetic_code(genetic_code)
    if codon_alignment is not None:
        write_alignment(names, codon_alignment, outdir, "alignment.cds.fasta")

    # Rewrite IQ-TREE's own Newick text, swapping safe IDs for quoted original names,
    # rather than re-serializing: the tree is published exactly as IQ-TREE wrote it.
    newick = asr_tree.read_text().strip()
    supports = supports_by_node(newick, support_tree.read_text())
    restore = restorer(names)
    (outdir / "tree.nwk").write_text(restore(newick) + "\n")
    rooted_newick = rooted_tree.read_text().strip() if rooted_tree else None
    if rooted_newick is not None and rooting_json is not None:
        report_json = json.loads(rooting_json.read_text())
        write_rooted(restore, names, rooted_newick, report_json, outdir)
    write_supports(outdir / "supports.tsv", supports, support_kinds)

    states: dict[str, list[str]] = {}
    coding: dict[str, list[str]] = {}  # codon models: each node's codons
    with (
        open(statefile) as src,
        open(outdir / "ancestors_pp.tsv", "w") as pp,
        open(outdir / "ancestors_posteriors.tsv", "w") as post,
    ):
        pp.write("node\tsite\tstate\tpp\n")
        codon_pp = open(outdir / "ancestors_codon_pp.tsv", "w") if codon_alignment else None
        codon_post = (
            open(outdir / "ancestors_codon_posteriors.tsv", "w") if codon_alignment else None
        )
        if codon_pp and codon_post:
            codon_pp.write("node\tsite\tcodon\tpp\n")
            codon_post.write("node\tsite\tcodon\tp\n")
        header: list[str] = []
        for line in src:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "Node":
                header = [h[2:] for h in fields[3:]]
                columns = list(codons.IQTREE_AMINO_ORDER) if codon_pp else header
                post.write("node\tsite\t" + "\t".join(columns) + "\n")
                continue
            node, site, state, *probs = fields
            if codon_pp:
                by_codon = {c: float(p) for c, p in zip(header, probs, strict=True)}
                amino, codon = codons.amino_posteriors(by_codon, code)
                state = max(amino, key=lambda a: amino[a])
                probs = [f"{amino[a]:.5f}" for a in codons.IQTREE_AMINO_ORDER]
                coding.setdefault(node, []).append(codon)
                codon_pp.write(f"{node}\t{site}\t{codon}\t{by_codon[codon]:.5f}\n")
                assert codon_post is not None
                codon_post.writelines(codons.codon_lines(node, int(site), by_codon))
            best = max(float(p) for p in probs)
            states.setdefault(node, []).append(state)
            pp.write(f"{node}\t{site}\t{state}\t{best:.5f}\n")
            post.write(f"{node}\t{site}\t" + "\t".join(probs) + "\n")
        if codon_pp and codon_post:
            codon_pp.close()
            codon_post.close()
    # Where each node's rows are, for the viewers' ranged reads.
    index_by_node(outdir / "ancestors_posteriors.tsv")
    if codon_alignment:
        index_by_node(outdir / "ancestors_codon_posteriors.tsv")

    with open(outdir / "ancestors.fasta", "w") as handle:
        for node, residues in states.items():
            handle.write(f">{node}\n{''.join(residues)}\n")
    if coding:
        write_fasta(outdir / "ancestors.cds.fasta", {n: "".join(c) for n, c in coding.items()})
    place_gaps(
        {node: "".join(residues) for node, residues in states.items()},
        read_fasta(alignment),
        newick,
        outdir,
        gap_methods,
        gap_state,
        gap_tree,
        external or {},
        rooted_newick,
    )
    if coding:
        mask_coding_ancestors(outdir, gap_methods, {n: "".join(c) for n, c in coding.items()})

    write_model_and_fit(model_path, report, outdir, ranking_report)


def mask_coding_ancestors(outdir: Path, methods: tuple[str, ...], cds: dict[str, str]) -> None:
    """ancestors.gaps-<method>.cds.fasta for each masking method: a codon becomes '---'
    where the masked protein ancestor has a gap (PRANK and GRASP reconstruct proteins
    of their own, so they get no coding version)."""
    for method in methods:
        if method in ("prank", "grasp"):
            continue
        masked = read_fasta(outdir / f"ancestors.gaps-{method}.fasta")
        rows = {
            node: "".join(
                "---" if residue == "-" else cds[node][3 * i : 3 * i + 3]
                for i, residue in enumerate(protein)
            )
            for node, protein in masked.items()
        }
        write_fasta(outdir / f"ancestors.gaps-{method}.cds.fasta", rows)


def place_gaps(
    ancestors: dict[str, str],
    alignment: dict[str, str],
    asr_newick: str,
    outdir: Path,
    methods: tuple[str, ...],
    gap_state: Path | None,
    gap_tree: Path | None,
    external: dict[str, tuple[Path, Path]] | None = None,
    rooted_newick: str | None = None,
) -> None:
    """Write ancestors.gaps-<method>.fasta for each requested method: masked copies for
    parsimony, fitch and ml, other programs' reconstructions (renamed to our nodes) for
    prank and grasp. Fitch uses `rooted_newick` if given, else the tree as rooted."""
    external = external or {}
    for method in methods:
        if method in ("prank", "grasp"):
            if method not in external:
                raise ValueError(f"gap method {method} needs --{method}-ancestors/--{method}-tree")
            fasta, tree = external[method]
            rows = gaps.external_ancestors(asr_newick, tree.read_text(), fasta, alignment)
            write_fasta(outdir / f"ancestors.gaps-{method}.fasta", rows)
            continue
        if method == "parsimony":
            presence = gaps.parsimony_presence(asr_newick, alignment)
        elif method == "fitch":
            presence = gaps.fitch_presence(rooted_newick or asr_newick, alignment)
        elif method == "ml":
            if gap_state is None or gap_tree is None:
                raise ValueError("gap method ml needs --gap-state and --gap-tree")
            presence, p_gap = gaps.ml_presence(asr_newick, gap_tree.read_text(), gap_state)
            with open(outdir / "gaps-ml.tsv", "w") as handle:
                handle.write("node\tsite\tp_gap\n")
                for node, probabilities in p_gap.items():
                    handle.writelines(
                        f"{node}\t{site}\t{p:.5f}\n" for site, p in enumerate(probabilities, 1)
                    )
        else:
            raise ValueError(f"unknown gap method {method!r}")
        write_fasta(outdir / f"ancestors.gaps-{method}.fasta", gaps.mask(ancestors, presence))


def main() -> None:
    parser = argparse.ArgumentParser(usage=__doc__)
    for name in ("names", "alignment", "asr_tree", "support_tree", "statefile", "model",
                 "report", "outdir"):  # fmt: skip
        parser.add_argument(name, type=Path)
    parser.add_argument("--gap-methods", default="")
    parser.add_argument("--supports", default="ufboot")
    parser.add_argument("--modelfinder-report", type=Path)
    parser.add_argument("--gap-state", type=Path)
    parser.add_argument("--gap-tree", type=Path)
    parser.add_argument("--rooted-tree", type=Path)
    parser.add_argument("--rooting", type=Path)
    for program in ("prank", "grasp"):
        parser.add_argument(f"--{program}-ancestors", type=Path)
        parser.add_argument(f"--{program}-tree", type=Path)
    parser.add_argument("--codon-alignment", type=Path)
    parser.add_argument("--genetic-code", type=int, default=1)
    parser.add_argument("--site-rates", type=Path)
    parser.add_argument("--site-lh", type=Path)
    args = parser.parse_args()
    methods = tuple(m for m in args.gap_methods.split(",") if m)
    external = {
        program: (getattr(args, f"{program}_ancestors"), getattr(args, f"{program}_tree"))
        for program in ("prank", "grasp")
        if getattr(args, f"{program}_ancestors") is not None
    }
    export(
        args.names, args.alignment, args.asr_tree, args.support_tree, args.statefile,
        args.model, args.report, args.outdir, methods, args.gap_state, args.gap_tree,
        external, args.rooted_tree, args.rooting, measures(args.supports),
        args.modelfinder_report, args.codon_alignment, args.genetic_code,
    )  # fmt: skip
    write_site_statistics(args.site_rates, args.site_lh, args.outdir)


def write_site_statistics(rates: Path | None, site_lh: Path | None, outdir: Path) -> None:
    """site_rates.tsv and site_lh.tsv, one row per alignment column (see the docstring)."""
    if rates is not None:
        rows = [
            line.split()
            for line in rates.read_text().splitlines()
            if line.strip() and not line.startswith("#") and not line.startswith("Site")
        ]
        if rows:  # none without rate heterogeneity
            write_tsv(outdir / "site_rates.tsv", ("site", "rate", "category", "category_rate"),
                      [dict(zip(("site", "rate", "category", "category_rate"), r, strict=True))
                       for r in rows])  # fmt: skip
    if site_lh is not None:
        lines = site_lh.read_text().split("\n", 1)[1].split()
        if not lines or lines[0] != "Site_Lh":
            raise ValueError(f"{site_lh} is not IQ-TREE's site log-likelihood file")
        lh_rows = [{"site": str(i), "log_likelihood": v} for i, v in enumerate(lines[1:], 1)]
        write_tsv(outdir / "site_lh.tsv", ("site", "log_likelihood"), lh_rows)


def measures(text: str) -> tuple[str, ...]:
    """--supports: comma-separated, in label order; '' when none ran."""
    return tuple(m for m in text.split(",") if m)


def main_alignment() -> None:
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    export_alignment(*(Path(a) for a in sys.argv[1:]))


def main_tree() -> None:
    parser = argparse.ArgumentParser(usage=__doc__)
    for name in ("names", "alignment", "tree", "model", "report", "outgroup", "outdir"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--supports", default="ufboot")
    parser.add_argument("--modelfinder-report", type=Path)
    args = parser.parse_args()
    export_tree(
        args.names, args.alignment, args.tree, args.model, args.report, args.outgroup,
        args.outdir, measures(args.supports), args.modelfinder_report,
    )  # fmt: skip


if __name__ == "__main__":
    main()
