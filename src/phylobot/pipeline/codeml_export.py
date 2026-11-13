"""codeml's results for one combination, in PhyloBot's files.

    phylobot-codeml-export ASR_TREE GENETIC_CODE RUNS_DIR OUTDIR NAMES

RUNS_DIR holds <model>.start<n>.mlc and .rst for every codeml run of this combination
(one per site model and starting omega). Per model the run with the best likelihood is
kept. OUTDIR gets:

  codeml.tsv      one row per model: lnL, parameters, kappa, tree length, the dN/dS
                  estimates, and every start's lnL
  selection.tsv   likelihood-ratio tests of positive selection (M1a vs M2a, M7 vs M8),
                  for the pairs that ran
  beb.tsv         Bayes empirical Bayes sites of M2a and M8 (posterior P(dN/dS > 1))
  <model>/        ancestors.fasta, ancestors.cds.fasta, ancestors_pp.tsv,
                  ancestors_posteriors.tsv, ancestors_codon_pp.tsv,
                  ancestors_codon_posteriors.tsv: as for IQ-TREE's
                  codon reconstructions (codon posteriors summed into amino acids, the
                  coding ancestor the best codon for the best amino acid), on the same
                  node names (ASR_TREE: the combination's IQ-TREE tree), so the viewers
                  read them alike
  FR/branches.tsv the free-ratio model's estimates per branch (parent, child: our node
                  names, extant sequences under their original names from NAMES, the
                  job's names.tsv), from mlc's per-branch table

Every model's per-branch table is in mlc, but only the free-ratio model gives each
branch its own dN/dS, so only its table is exported.
"""

import re
import sys
from pathlib import Path

from phylobot import codeml, codons, gaps
from phylobot.pipeline.export import index_by_node, read_names
from phylobot.seqio import write_fasta

RUN = re.compile(r"^(M\w+|FR)\.start(\d+)\.mlc$")


def export(
    asr_tree: Path,
    code: codons.GeneticCode,
    runs_dir: Path,
    outdir: Path,
    names: dict[str, str] | None = None,
) -> None:
    runs: dict[str, list[tuple[int, codeml.Fit, Path]]] = {}
    for path in sorted(runs_dir.iterdir()):
        if (match := RUN.match(path.name)) and match.group(1) in codeml.MODELS:
            fit = codeml.parse_mlc(path.read_text())
            runs.setdefault(match.group(1), []).append((int(match.group(2)), fit, path))
    if not runs:
        raise ValueError(f"no codeml runs in {runs_dir}")
    outdir.mkdir(parents=True, exist_ok=True)
    ours = asr_tree.read_text().strip()
    rows, lnl, beb = [], {}, []
    for key in [m.key for m in codeml.SITE_MODELS if m.key in runs]:
        start, fit, mlc = max(runs[key], key=lambda run: run[1].lnl)
        model = codeml.MODELS[key]
        lnl[key] = fit.lnl
        starts = "; ".join(f"{model.starts[n - 1]:g}: {f.lnl:.4f}" for n, f, _ in sorted(runs[key]))
        rows.append([key, f"{fit.lnl:.6f}", str(fit.parameters), f"{fit.kappa:.5g}",
                     f"{fit.tree_length:.5g}", fit.summary(key), starts])  # fmt: skip
        beb += [
            [key, str(s), r, f"{p:.3f}", f"{m:.3f}", f"{e:.3f}"] for s, r, p, m, e in fit.beb or []
        ]
        reconstruction = codeml.parse_rst(mlc.with_suffix(".rst").read_text())
        write_ancestors(ours, reconstruction, code, outdir / key)
        if key == codeml.FREE_RATIOS:
            branches = codeml.parse_branches(mlc.read_text())
            write_branches(ours, reconstruction, branches, names or {}, outdir / key)
    write_table(outdir / "codeml.tsv",
                ["model", "log_likelihood", "free_parameters", "kappa", "tree_length",
                 "estimates", "starts"], rows)  # fmt: skip
    tests = [
        [f"{t.null} vs {t.alternative}", f"{t.statistic:.4f}", "2", f"{t.p_value:.4g}"]
        for t in codeml.likelihood_ratio_tests(lnl)
    ]
    write_table(outdir / "selection.tsv", ["test", "statistic", "df", "p_value"], tests)
    write_table(
        outdir / "beb.tsv", ["model", "site", "residue", "p_positive", "mean_omega", "se"], beb
    )


def write_ancestors(
    ours: str, reconstruction: codeml.Reconstruction, code: codons.GeneticCode, folder: Path
) -> None:
    """codeml's ancestors under our node names (matched by the subtrees meeting at each
    node, which is independent of rooting)."""
    names = gaps.match_nodes(ours, reconstruction.tree)
    if len(names) != len(reconstruction.posteriors):
        raise ValueError("codeml's tree does not match the combination's tree")
    folder.mkdir(parents=True, exist_ok=True)
    proteins, coding = {}, {}
    with (
        open(folder / "ancestors_pp.tsv", "w") as pp,
        open(folder / "ancestors_posteriors.tsv", "w") as post,
        open(folder / "ancestors_codon_pp.tsv", "w") as codon_pp,
        open(folder / "ancestors_codon_posteriors.tsv", "w") as codon_post,
    ):
        pp.write("node\tsite\tstate\tpp\n")
        post.write("node\tsite\t" + "\t".join(codons.IQTREE_AMINO_ORDER) + "\n")
        codon_pp.write("node\tsite\tcodon\tpp\n")
        codon_post.write("node\tsite\tcodon\tp\n")
        for node, label in sorted(names.items(), key=lambda item: node_order(item[0])):
            residues, triplets = [], []
            for site, probs in enumerate(reconstruction.posteriors[label], 1):
                amino, codon = codons.amino_posteriors(probs, code)
                best = max(amino, key=lambda a: amino[a])
                residues.append(best)
                triplets.append(codon)
                pp.write(f"{node}\t{site}\t{best}\t{amino[best]:.5f}\n")
                values = "\t".join(f"{amino[a]:.5f}" for a in codons.IQTREE_AMINO_ORDER)
                post.write(f"{node}\t{site}\t{values}\n")
                codon_pp.write(f"{node}\t{site}\t{codon}\t{probs[codon]:.5f}\n")
                codon_post.writelines(codons.codon_lines(node, site, probs))
            proteins[node], coding[node] = "".join(residues), "".join(triplets)
    write_fasta(folder / "ancestors.fasta", proteins)
    write_fasta(folder / "ancestors.cds.fasta", coding)
    # Where each node's rows are, as for IQ-TREE's reconstruction.
    index_by_node(folder / "ancestors_posteriors.tsv")
    index_by_node(folder / "ancestors_codon_posteriors.tsv")


BRANCH_FIELDS = ["parent", "child", "t", "N", "S", "dN_dS", "dN", "dS", "N_dN", "S_dS"]


def write_branches(
    ours: str,
    reconstruction: codeml.Reconstruction,
    branches: list[codeml.Branch],
    names: dict[str, str],
    folder: Path,
) -> None:
    """branches.tsv: codeml's branches under our names (ancestors by the subtrees meeting
    at each node, extant sequences by codeml's tip numbers, then their original names)."""
    by_number = {
        int(label[1:]): node for node, label in gaps.match_nodes(ours, reconstruction.tree).items()
    }
    by_number |= {n: names.get(safe, safe) for n, safe in reconstruction.tips.items()}
    rows = []
    for b in branches:
        if b.parent not in by_number or b.child not in by_number:
            raise ValueError(f"codeml branch {b.parent}..{b.child} has no node in our tree")
        rows.append([by_number[b.parent], by_number[b.child], *(
            f"{v:g}" for v in (b.t, b.n, b.s, b.omega, b.dn, b.ds, b.n_dn, b.s_ds)
        )])  # fmt: skip
    write_table(folder / "branches.tsv", BRANCH_FIELDS, rows)


def node_order(name: str) -> tuple[int, str]:
    """Node1, Node2, ... in number order (IQ-TREE's names), anything else after."""
    digits = re.sub(r"\D", "", name)
    return (int(digits) if digits else 1 << 30, name)


def write_table(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.write_text("\n".join("\t".join(row) for row in [header, *rows]) + "\n")


def main() -> None:
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    tree, code, runs, outdir, names = sys.argv[1:]
    try:
        export(
            Path(tree), codons.genetic_code(code), Path(runs), Path(outdir),
            read_names(Path(names)),
        )  # fmt: skip
    except (ValueError, codons.CodonError) as exc:
        print(f"codeml export error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
