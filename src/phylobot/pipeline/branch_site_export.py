"""The branch-site test's results, in PhyloBot's files.

    phylobot-branch-site-export TREE RUNS_DIR OUTDIR NAMES [ALIGNMENT]

TREE is codeml's tree (safe IDs, #1 on the foreground branches; VALIDATE's), RUNS_DIR
holds A.start<n>.mlc and A_null.start1.mlc (with their .ctl control files, and .rst
files, unused), NAMES the job's names.tsv, ALIGNMENT the codon alignment codeml read.
Per model the run with the best likelihood is kept. OUTDIR gets:

  models.tsv      Model A and its null: lnL, parameters, kappa, the site classes
                  (proportions, background and foreground dN/dS), every start's lnL
  test.tsv        the likelihood-ratio test: 2 dlnL, df, p (chi-squared, 1 df) and p
                  under the 50:50 mixture of 0 and chi-squared 1 df
  beb.tsv         Model A's Bayes empirical Bayes sites: posterior probability of being
                  under positive selection on the foreground (codeml lists P > 0.5)
  sites.tsv       every site's BEB posteriors of the four site classes (0, 1, 2a, 2b,
                  from Model A's rst) and of positive selection (2a + 2b)
  foreground.nwk  codeml's tree with the original sequence names and the #1 marks
  summary.json    {"branch_site": true, ...}, for the job page
  codeml/         codeml's own files, to check or rerun them as they ran (README.txt)
"""

import json
import re
import shutil
import sys
from pathlib import Path

from phylobot import branch_site
from phylobot.pipeline.codeml_export import write_table
from phylobot.pipeline.export import read_names, restorer

RUN = re.compile(r"^(A|A_null)\.start(\d+)\.mlc$")
CLASSES = ("0", "1", "2a", "2b")


README = """\
codeml's own files for this branch-site test (PAML 4.10.10), exactly as they ran:

  aln.fasta, tree.nwk   the codon alignment and the tree (#1 on the foreground
                        branches), with PhyloBot's sequence IDs: names.tsv maps them
                        to the original names
  <run>.ctl, .mlc, .rst each run's control file, codeml's main output and its
                        supplementary output (per-site class probabilities): A.start<n>
                        is Model A from the n-th starting omega, A_null.start1 its null
                        (omega2 fixed at 1); PhyloBot keeps each model's best start

To rerun one: codeml A.start1.ctl (it writes mlc in this folder).
"""


def write_sites(rst: Path, outdir: Path) -> None:
    """sites.tsv: every site's BEB class posteriors under Model A (its best start) and
    their sum for the foreground's positively selected classes (2a + 2b)."""
    rows = [
        [str(site), residue, *(f"{p:.5f}" for p in classes), f"{classes[2] + classes[3]:.5f}"]
        for site, residue, classes in branch_site.parse_site_classes(rst.read_text())
    ]
    write_table(outdir / "sites.tsv", ["site", "residue", "p0", "p1", "p2a", "p2b",
                                       "p_positive"], rows)  # fmt: skip


def write_codeml_files(
    runs_dir: Path, tree: Path, alignment: Path | None, names: Path | None, folder: Path
) -> None:
    """codeml/: the runs' control files and mlc, and the inputs they name."""
    folder.mkdir(parents=True, exist_ok=True)
    for path in sorted(runs_dir.iterdir()):
        if path.suffix in (".ctl", ".mlc", ".rst"):
            shutil.copyfile(path, folder / path.name)
    newick = tree.read_text().strip()
    sequences = len(re.findall(r"\bs\d{4,}\b", newick))
    (folder / "tree.nwk").write_text(f"{sequences} 1\n{newick}\n")  # codeml's header line
    if alignment is not None:
        shutil.copyfile(alignment, folder / "aln.fasta")
    if names is not None:
        shutil.copyfile(names, folder / "names.tsv")
    (folder / "README.txt").write_text(README)


def export(
    tree: Path,
    runs_dir: Path,
    outdir: Path,
    names: dict[str, str],
    alignment: Path | None = None,
    names_file: Path | None = None,
) -> None:
    runs: dict[str, list[tuple[int, branch_site.Fit]]] = {}
    for path in sorted(runs_dir.iterdir()):
        if match := RUN.match(path.name):
            fit = branch_site.parse_mlc(path.read_text())
            runs.setdefault(match.group(1), []).append((int(match.group(2)), fit))
    if set(runs) != {branch_site.ALTERNATIVE, branch_site.NULL}:
        raise ValueError(f"codeml runs of Model A and its null are missing from {runs_dir}")
    outdir.mkdir(parents=True, exist_ok=True)
    best: dict[str, branch_site.Fit] = {}
    rows = []
    for key in (branch_site.ALTERNATIVE, branch_site.NULL):
        start, fit = max(runs[key], key=lambda run: run[1].lnl)
        best[key] = fit
        if key == branch_site.ALTERNATIVE:
            best_start = start
        c = fit.classes
        starts = branch_site.STARTS[key]
        tried = "; ".join(f"{starts[n - 1]:g}: {f.lnl:.4f}" for n, f in sorted(runs[key]))
        rows.append([
            key, f"{fit.lnl:.6f}", str(fit.parameters), f"{fit.kappa:.5g}",
            *(f"{p:.5g}" for p in c.proportions),
            f"{c.background[0]:.5g}", f"{c.foreground[2]:.5g}", tried,
        ])  # fmt: skip
    write_table(
        outdir / "models.tsv",
        ["model", "log_likelihood", "free_parameters", "kappa",
         *(f"p{k}" for k in CLASSES), "omega0", "omega2", "starts"],
        rows,
    )  # fmt: skip
    result = branch_site.test(best[branch_site.ALTERNATIVE].lnl, best[branch_site.NULL].lnl)
    write_table(
        outdir / "test.tsv",
        ["statistic", "df", "p_value", "p_mixture"],
        [[f"{result.statistic:.4f}", "1", f"{result.p_value:.4g}", f"{result.p_mixture:.4g}"]],
    )
    write_sites(runs_dir / f"{branch_site.ALTERNATIVE}.start{best_start}.rst", outdir)
    beb = best[branch_site.ALTERNATIVE].beb
    write_table(
        outdir / "beb.tsv",
        ["site", "residue", "p_positive"],
        [[str(s), r, f"{p:.3f}"] for s, r, p in beb],
    )
    (outdir / "foreground.nwk").write_text(restorer(names)(tree.read_text().strip()) + "\n")
    write_codeml_files(runs_dir, tree, alignment, names_file, outdir / "codeml")
    summary = {"branch_site": True, "stop_after": "branch_site", "alignments": [],
               "combinations": []}  # fmt: skip
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main() -> None:
    if len(sys.argv) not in (5, 6):
        sys.exit(__doc__)
    tree, runs, outdir, names, *alignment = sys.argv[1:]
    try:
        export(
            Path(tree), Path(runs), Path(outdir), read_names(Path(names)),
            Path(alignment[0]) if alignment else None, Path(names),
        )  # fmt: skip
    except ValueError as exc:
        print(f"branch-site export error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
