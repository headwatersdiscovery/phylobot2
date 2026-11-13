"""PAML codeml: codon site models as a second codon reconstruction, and the free-ratio
branch model (one dN/dS per branch).

codeml runs on fixed trees (each alignment's best codon combination by default; its
own tree search is not competitive) with the site models below. The model definitions
live here only: VALIDATE writes one control file per model and starting omega
(`write_controls`), CODEML runs them, CODEML_EXPORT parses the output (`parse_mlc`,
`parse_rst`) and tests M1a against M2a and M7 against M8 (`likelihood_ratio_tests`).
The free-ratio model's per-branch estimates come from mlc's "dN & dS for each branch"
table (`parse_branches`).

Checked with PAML 4.10.10 on lysin (2026-10-01): codeml reads FASTA and a Newick tree
without internal labels; it writes marginal (not joint) ancestors under every site model.
"""

import math
import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class SiteModel:
    key: str
    nssites: int
    ncatg: int
    label: str
    starts: tuple[float, ...]  # starting omegas; the best likelihood is kept
    default: bool  # ticked on the New Job page (owner decision 2026-10-01)
    model: int = 0  # codeml's branch model: 0 one dN/dS for all branches, 1 free ratios


# Free ratios: one dN/dS per branch (codeml model = 1, NSsites = 0), the branch
# estimates for mutational paths (owner decision: before the branch-site test). Many
# parameters: estimates, not a test.
FREE_RATIOS = "FR"


# M2a and M8 have known local optima (Yang 2000; Wong et al. 2004): two starting omegas.
SITE_MODELS: tuple[SiteModel, ...] = (
    SiteModel("M0", 0, 1, "M0: one dN/dS for all sites", (0.4,), True),
    SiteModel("M1a", 1, 2, "M1a: nearly neutral (dN/dS < 1 or = 1)", (0.4,), False),
    SiteModel("M2a", 2, 3, "M2a: positive selection (adds dN/dS > 1)", (0.5, 2.0), False),
    SiteModel("M3", 3, 3, "M3: three dN/dS classes", (0.4,), True),
    SiteModel("M7", 7, 10, "M7: beta-distributed dN/dS (0-1)", (0.4,), True),
    SiteModel("M8", 8, 10, "M8: beta + dN/dS > 1", (0.5, 2.0), True),
    SiteModel(FREE_RATIOS, 0, 1, "Free ratios: one dN/dS per branch", (0.4,), False, 1),
)
MODELS = {m.key: m for m in SITE_MODELS}
DEFAULT_MODELS = [m.key for m in SITE_MODELS if m.default]
# Nested pairs tested for positive selection (2 degrees of freedom each).
TESTS = (("M1a", "M2a"), ("M7", "M8"))
WHERE = ("best", "all")  # each alignment's best codon combination, or every one

CONTROL = """\
      seqfile = aln.fasta
     treefile = tree.nwk
      outfile = mlc
        noisy = 0
      verbose = 2
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        model = {model}
      NSsites = {nssites}
        icode = {icode}
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = {omega}
        ncatG = {ncatg}
    cleandata = 0
 RateAncestor = 1
"""


def write_controls(outdir: Path, models: list[str], icode: int) -> None:
    """codeml/<model>/start<n>.ctl for each model and starting omega (F3X4 codon
    frequencies, CodonFreq = 2, as codeml's own examples)."""
    for key in models:
        model = MODELS[key]
        folder = outdir / key
        folder.mkdir(parents=True, exist_ok=True)
        for n, omega in enumerate(model.starts, 1):
            text = CONTROL.format(
                model=model.model, nssites=model.nssites, icode=icode, omega=omega,
                ncatg=model.ncatg,
            )  # fmt: skip
            (folder / f"start{n}.ctl").write_text(text)


def chi2_sf_2df(statistic: float) -> float:
    """P(chi-squared with 2 degrees of freedom >= statistic): exp(-x/2), exactly."""
    return math.exp(-max(statistic, 0.0) / 2)


@dataclass(frozen=True)
class Test:
    null: str
    alternative: str
    statistic: float  # 2 (lnL1 - lnL0)
    p_value: float


def likelihood_ratio_tests(lnl: dict[str, float]) -> list[Test]:
    """M1a vs M2a and M7 vs M8, for the pairs that both ran."""
    tests = []
    for null, alternative in TESTS:
        if null in lnl and alternative in lnl:
            statistic = 2 * (lnl[alternative] - lnl[null])
            tests.append(Test(null, alternative, statistic, chi2_sf_2df(statistic)))
    return tests


LNL = re.compile(r"^lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\):\s+(-?[\d.]+)", re.MULTILINE)


def lnl_of(mlc: str) -> float:
    match = LNL.search(mlc)
    if match is None:
        raise ValueError("codeml output has no lnL line")
    return float(match.group(3))


FLOAT = r"(-?[\d.]+(?:e-?\d+)?)"


def floats(line: str) -> list[float]:
    return [float(x) for x in re.findall(FLOAT, line)]


@dataclass
class Fit:
    """One codeml run's estimates (mlc)."""

    lnl: float
    parameters: int
    kappa: float
    tree_length: float
    omega: float | None = None  # M0
    proportions: list[float] | None = None  # site classes (M1a, M2a, M3, M7, M8)
    omegas: list[float] | None = None
    beta: tuple[float, float] | None = None  # M7, M8: beta(p, q)
    beb: list[tuple[int, str, float, float, float]] | None = None  # site, residue, P, mean, se

    def summary(self, model: str) -> str:
        """The dN/dS estimates in one line, for tables."""
        p, w = self.proportions or [], self.omegas or []
        if model == FREE_RATIOS:
            return "one dN/dS per branch (branches.tsv)"
        if model == "M0" and self.omega is not None:
            return f"ω = {self.omega:.3g}"
        if model in ("M7", "M8") and self.beta:
            text = f"beta(p = {self.beta[0]:.3g}, q = {self.beta[1]:.3g})"
            if model == "M8" and p and w:
                text = f"{1 - p[-1]:.3g} {text}; {p[-1]:.3g} at ω = {w[-1]:.3g}"
            return text
        return "; ".join(f"{pi:.3g} at ω = {wi:.3g}" for pi, wi in zip(p, w, strict=True))


def parse_mlc(text: str) -> Fit:
    match = LNL.search(text)
    if match is None:
        raise ValueError("codeml output has no lnL line")
    fit = Fit(
        lnl=float(match.group(3)),
        parameters=int(match.group(2)),
        kappa=_value(text, r"^kappa \(ts/tv\) =\s*" + FLOAT),
        tree_length=_value(text, r"^tree length =\s*" + FLOAT),
    )
    if omega := re.search(r"^omega \(dN/dS\) =\s*" + FLOAT, text, re.MULTILINE):
        fit.omega = float(omega.group(1))
    if classes := re.search(r"site classes \(K=\d+\)\s*\n\s*\np:(.*)\nw:(.*)", text):
        fit.proportions, fit.omegas = floats(classes.group(1)), floats(classes.group(2))
    if beta := re.search(r"\bp =\s*" + FLOAT + r"\s+q =\s*" + FLOAT, text):
        fit.beta = (float(beta.group(1)), float(beta.group(2)))
    if "Bayes Empirical Bayes (BEB)" in text:
        block = text.split("Bayes Empirical Bayes (BEB)", 1)[1]
        fit.beb = [
            (int(m.group(1)), m.group(2), float(m.group(3)), float(m.group(4)), float(m.group(5)))
            for m in re.finditer(
                r"^\s*(\d+) (\S)\s+([\d.]+)\**\s+([\d.]+) \+- ([\d.]+)", block, re.MULTILINE
            )
        ]
    return fit


def _value(text: str, pattern: str) -> float:
    match = re.search(pattern, text, re.MULTILINE)
    if match is None:
        raise ValueError(f"codeml output lacks {pattern!r}")
    return float(match.group(1))


@dataclass(frozen=True)
class Branch:
    """One row of mlc's "dN & dS for each branch" table: codeml's node numbers (tips
    1..n in alignment order, then the ancestors), the branch length t, the synonymous
    and non-synonymous sites S and N, dN/dS, dN, dS, and the expected numbers of
    non-synonymous (N*dN) and synonymous (S*dS) changes on the branch."""

    parent: int
    child: int
    t: float
    n: float
    s: float
    omega: float
    dn: float
    ds: float
    n_dn: float
    s_ds: float


BRANCH_ROW = re.compile(r"^\s*(\d+)\.\.(\d+)" + r"\s+" + r"\s+".join([FLOAT] * 8) + r"\s*$")


def parse_branches(mlc: str) -> list[Branch]:
    """The per-branch estimates (every codeml model prints them; under the free-ratio
    model each branch has its own dN/dS). codeml writes 999 for dN/dS when dS is 0."""
    if "dN & dS for each branch" not in mlc:
        raise ValueError("codeml output has no per-branch dN & dS table")
    block = mlc.split("dN & dS for each branch", 1)[1].split("tree length for dN", 1)[0]
    branches = []
    for line in block.splitlines():
        if match := BRANCH_ROW.match(line):
            parent, child, *values = match.groups()
            branches.append(Branch(int(parent), int(child), *map(float, values)))
    if not branches:
        raise ValueError("codeml's per-branch table is empty")
    return branches


@dataclass
class Reconstruction:
    """codeml's marginal reconstruction (rst, verbose = 2): its tree with node numbers as
    labels, per node and site the codon posteriors (codons below 0.001 left out), and
    codeml's number for each tip."""

    tree: str
    posteriors: dict[str, list[dict[str, float]]]  # node label -> sites
    tips: dict[int, str] = field(default_factory=dict)  # codeml's number -> sequence name


CODON_PROB = re.compile(r"([ACGT]{3})\(([\d.]+)\)")


def parse_rst(text: str) -> Reconstruction:
    lines = text.splitlines()
    try:
        start = lines.index("tree with node labels for Rod Page's TreeView")
    except ValueError:
        raise ValueError("codeml rst has no labelled tree") from None
    labelled = lines[start + 1].strip()
    tips = {int(n): name for n, name in re.findall(r"(?:^|[(,\s])(\d+)_([^\s,()]+)", labelled)}
    newick = re.sub(r"\b\d+_(s\d+)", r"\1", labelled)  # "1_s0001" -> "s0001"
    newick = re.sub(r"\)\s*(\d+)\s*", r")N\1", newick)  # node numbers as labels
    if (root := re.search(r"Nodes (\d+) to \d+ are ancestral", text)) and newick.endswith(");"):
        newick = newick[:-1] + f"N{root.group(1)};"  # the root is left unnumbered
    posteriors: dict[str, list[dict[str, float]]] = {}
    node = None
    for line in lines:
        if line.startswith("Prob distribution at node "):
            node = "N" + line.split()[4].rstrip(",")
            posteriors[node] = []
        elif line.startswith("Prob of best state"):
            node = None  # the distributions are over
        elif node is not None and ":" in line and line.strip()[:1].isdigit():
            probs = {c: float(p) for c, p in CODON_PROB.findall(line.split(":", 1)[1])}
            posteriors[node].append({c: p for c, p in probs.items() if p >= 0.001})
    return Reconstruction(newick.replace(" ", ""), posteriors, tips)
