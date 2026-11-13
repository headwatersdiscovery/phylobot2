"""PAML codeml: models, control files, output parsing, tests and export (phylobot.codeml)."""

import math
from pathlib import Path

import pytest

from phylobot import codeml, codons
from phylobot.pipeline.codeml_export import export

FIXTURES = Path(__file__).parent / "fixtures" / "codeml"


def test_models_and_defaults() -> None:
    # Owner decision 2026-10-01: M0, M3, M7, M8 ticked; M1a/M2a offered, unticked.
    assert codeml.DEFAULT_MODELS == ["M0", "M3", "M7", "M8"]
    assert [m.key for m in codeml.SITE_MODELS] == ["M0", "M1a", "M2a", "M3", "M7", "M8", "FR"]
    assert codeml.MODELS["M8"].starts == (0.5, 2.0)  # local optima: two starts
    assert codeml.TESTS == (("M1a", "M2a"), ("M7", "M8"))


def test_control_files(tmp_path: Path) -> None:
    codeml.write_controls(tmp_path, ["M0", "M8"], icode=1)
    assert sorted(str(p.relative_to(tmp_path)) for p in tmp_path.rglob("*.ctl")) == [
        "M0/start1.ctl", "M8/start1.ctl", "M8/start2.ctl",
    ]  # fmt: skip
    m8 = (tmp_path / "M8" / "start2.ctl").read_text()
    assert "NSsites = 8" in m8 and "ncatG = 10" in m8 and "omega = 2.0" in m8
    assert "icode = 1" in m8 and "verbose = 2" in m8 and "RateAncestor = 1" in m8


def test_likelihood_ratio_tests() -> None:
    assert codeml.chi2_sf_2df(0) == 1.0
    assert codeml.chi2_sf_2df(5.991) == pytest.approx(0.05, abs=1e-4)  # chi2(2) 95% point
    [test] = codeml.likelihood_ratio_tests({"M7": -4524.684, "M8": -4464.620, "M0": -4682.37})
    assert (test.null, test.alternative) == ("M7", "M8")
    assert test.statistic == pytest.approx(120.128) and test.p_value < 1e-20
    assert codeml.likelihood_ratio_tests({"M2a": -1.0}) == []  # M1a didn't run


def test_parse_mlc_on_lysin() -> None:
    m0 = codeml.parse_mlc((FIXTURES / "lysin_M0.mlc").read_text())
    assert (m0.lnl, m0.parameters, m0.kappa, m0.omega) == (-4682.3731, 49, 1.56127, 0.92956)
    assert m0.summary("M0") == "ω = 0.93" and m0.beb is None
    m8 = codeml.parse_mlc((FIXTURES / "lysin_M8.mlc").read_text())
    assert m8.lnl == -4464.620176 and m8.beta == (0.3103, 0.29107)
    assert m8.proportions and m8.omegas and len(m8.proportions) == 11
    assert m8.summary("M8") == "0.731 beta(p = 0.31, q = 0.291); 0.269 at ω = 2.95"
    assert m8.beb and m8.beb[0] == (4, "H", 1.0, 2.872, 0.483)  # BEB, not the NEB table


RST = """\
Supplemental results for CODEML (seqf: aln.fasta  treef: tree.nwk)

tree with node labels for Rod Page's TreeView
((1_s0001, 2_s0002) 6 , 3_s0003, 4_s0004);

Nodes 5 to 6 are ancestral

Prob distribs at nodes, those with p < 0.001 not listed

Prob distribution at node 5, by site

   site  Freq   Data

      1      1   ATG (M) ATG (M) ATG (M) ATG (M) : ATG(1.000) AAA(0.000)
      2      1   CTT (L) ATT (I) ATC (I) CTT (L) : CTT(0.400) ATT(0.350) ATC(0.250)

Prob distribution at node 6, by site

   site  Freq   Data

      1      1   ATG (M) ATG (M) ATG (M) ATG (M) : ATG(1.000)
      2      1   CTT (L) ATT (I) ATC (I) CTT (L) : CTT(0.900) ATT(0.100)

Prob of best state at each node, listed by site
"""


def test_parse_rst() -> None:
    reconstruction = codeml.parse_rst(RST)
    assert reconstruction.tree == "((s0001,s0002)N6,s0003,s0004)N5;"
    assert reconstruction.posteriors["N5"][1] == {"CTT": 0.4, "ATT": 0.35, "ATC": 0.25}
    assert reconstruction.posteriors["N6"][0] == {"ATG": 1.0}  # 0.000 entries left out


def test_export(tmp_path: Path) -> None:
    runs = tmp_path / "runs"
    runs.mkdir()
    m0 = (FIXTURES / "lysin_M0.mlc").read_text()
    m8 = (FIXTURES / "lysin_M8.mlc").read_text()
    (runs / "M0.start1.mlc").write_text(m0)
    (runs / "M8.start1.mlc").write_text(m8.replace("-4464.620176", "-4470.000000"))
    (runs / "M8.start2.mlc").write_text(m8)  # the better start
    for name in ("M0.start1", "M8.start1", "M8.start2"):
        (runs / f"{name}.rst").write_text(RST)
    tree = tmp_path / "asr.treefile"
    tree.write_text("((s0001:0.1,s0002:0.1)Node2:0.1,s0003:0.1,s0004:0.1)Node1;\n")
    export(tree, codons.STANDARD, runs, tmp_path / "out")
    out = tmp_path / "out"
    rows = [line.split("\t") for line in (out / "codeml.tsv").read_text().splitlines()]
    assert [r[0] for r in rows[1:]] == ["M0", "M8"]
    assert rows[2][1] == "-4464.620176" and rows[2][6] == "0.5: -4470.0000; 2: -4464.6202"
    assert (out / "selection.tsv").read_text().splitlines() == ["test\tstatistic\tdf\tp_value"]
    assert (out / "beb.tsv").read_text().splitlines()[1].startswith("M8\t4\tH\t1.000")
    # Ancestors under our node names; I wins at site 2 of Node1 (0.6 summed) with ATT.
    assert (out / "M8" / "ancestors.fasta").read_text().split() == [">Node1", "MI", ">Node2", "ML"]
    assert (out / "M8" / "ancestors.cds.fasta").read_text().split()[1] == "ATGATT"
    pp = (out / "M8" / "ancestors_pp.tsv").read_text().splitlines()
    assert pp[2] == "Node1\t2\tI\t0.60000"
    assert math.isclose(float(pp[4].split("\t")[3]), 0.9)
    codon_rows = (out / "M8" / "ancestors_codon_posteriors.tsv").read_text().splitlines()
    assert "Node1\t2\tATT\t0.35000" in codon_rows


def test_parse_branches_on_lysin() -> None:
    # Free ratios on lysin (PAML 4.10.10, 2026-10-02): 47 branches; 999 where dS is 0.
    mlc = (FIXTURES / "lysin_FR.mlc").read_text()
    branches = codeml.parse_branches(mlc)
    assert len(branches) == 47
    assert branches[0] == codeml.Branch(26, 1, 0.107, 305.7, 99.3, 999.0, 0.0472, 0.0, 14.4, 0.0)
    assert branches[-1].omega == 0.3705 and branches[-1].s_ds == 21.0
    assert codeml.parse_mlc(mlc).parameters == 95
    with pytest.raises(ValueError, match="per-branch"):
        codeml.parse_branches("lnL(ntime: 1  np: 2):  -1.0")


FR_MLC = """\
lnL(ntime:  5  np: 11):  -100.000000      +0.000000
tree length =   0.50000
kappa (ts/tv) =  2.00000

dN & dS for each branch

 branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS

   5..6      0.100     4.5     1.5  0.5000  0.0200  0.0400   0.1   0.1
   6..1      0.100     4.5     1.5 999.0000  0.0300  0.0000   0.1   0.0
   6..2      0.100     4.5     1.5  1.0000  0.0300  0.0300   0.1   0.0
   5..3      0.100     4.5     1.5  2.0000  0.0400  0.0200   0.2   0.0
   5..4      0.100     4.5     1.5  0.2500  0.0100  0.0400   0.0   0.1

tree length for dN:       0.1300
"""


def test_export_free_ratio_branches(tmp_path: Path) -> None:
    runs = tmp_path / "runs"
    runs.mkdir()
    (runs / "FR.start1.mlc").write_text(FR_MLC)
    (runs / "FR.start1.rst").write_text(RST)
    tree = tmp_path / "asr.treefile"
    tree.write_text("((s0001:0.1,s0002:0.1)Node2:0.1,s0003:0.1,s0004:0.1)Node1;\n")
    names = {"s0001": "alpha", "s0002": "beta", "s0003": "gamma", "s0004": "delta"}
    export(tree, codons.STANDARD, runs, tmp_path / "out", names)
    out = tmp_path / "out"
    rows = [line.split("\t") for line in (out / "FR" / "branches.tsv").read_text().splitlines()]
    assert rows[0] == ["parent", "child", "t", "N", "S", "dN_dS", "dN", "dS", "N_dN", "S_dS"]
    assert [r[:2] for r in rows[1:]] == [
        ["Node1", "Node2"], ["Node2", "alpha"], ["Node2", "beta"], ["Node1", "gamma"],
        ["Node1", "delta"],
    ]  # fmt: skip
    assert rows[2][5] == "999" and rows[1][5] == "0.5"
    codeml_rows = (out / "codeml.tsv").read_text().splitlines()
    assert codeml_rows[1].startswith("FR\t-100.000000\t11")
    assert "one dN/dS per branch" in codeml_rows[1]
