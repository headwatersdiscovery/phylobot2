"""The pipeline's Python steps (phylobot.pipeline), run as the pipeline runs them."""

import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from phylobot import ancestors, catalog
from phylobot.pipeline import export
from phylobot.pipeline.export import parse_fit, parse_modelfinder, quote, supports_by_node

FIXTURES = Path(__file__).parent / "fixtures" / "fasta"
PIPELINE = Path(__file__).parent.parent / "pipeline"


def module(step: str) -> list[str]:
    return [sys.executable, "-m", f"phylobot.pipeline.{step}"]


LIMITS = {"max_sequences": 500, "max_length": 5000, "max_provided": 5}
FOUR = "".join(f">seq{i}\nMKTAYIAK{'LV'[i % 2]}Q\n" for i in range(4))


def validate(
    tmp_path: Path,
    sequences: str | None = None,
    alignments: list[str] = [],  # noqa: B006 - never mutated
    trees: list[str] = [],  # noqa: B006
    labels: dict[str, list[str]] | None = None,
    models: list[str] | None = None,
    raxml: bool = False,
    candidates: list[str] | None = None,
    seqtype: str | None = None,
    genetic_code: int | str | None = None,
    codon_models: list[str] | None = None,
    keys: list[str] | None = None,
    keep_lengths: list[bool] | None = None,
    **limits: int,
) -> str:
    """Stage files and a manifest as the workflow does; return stderr ('' on success)."""
    labels = labels or {}
    if sequences is not None:
        (tmp_path / "input.fasta").write_text(sequences)
    for i, text in enumerate(alignments, 1):
        (tmp_path / f"alignment{i}.fasta").write_text(text)
    for i, text in enumerate(trees, 1):
        (tmp_path / f"tree{i}.nwk").write_text(text)
    manifest = {
        "sequences": sequences is not None,
        "alignments": [
            {
                "filename": f"aln{i}.fasta",
                "label": (labels.get("alignments") or [""] * 9)[i],
                **({"key": keys[i]} if keys else {}),
            }
            for i in range(len(alignments))
        ],
        "trees": [
            {
                "filename": f"t{i}.nwk",
                "label": (labels.get("trees") or [""] * 9)[i],
                **({"keep_lengths": True} if keep_lengths and keep_lengths[i] else {}),
            }
            for i in range(len(trees))
        ],
        **LIMITS,
        **limits,
        "models": models,
        "raxml": raxml,
        **({"candidates": candidates} if candidates is not None else {}),
        **({"seqtype": seqtype} if seqtype is not None else {}),
        **({"genetic_code": genetic_code} if genetic_code is not None else {}),
        **({"codon_models": codon_models} if codon_models is not None else {}),
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    proc = subprocess.run(
        [*module("validate"), "manifest.json"], cwd=tmp_path, capture_output=True, text=True
    )
    return proc.stderr if proc.returncode else ""


def test_edge_case_names_get_safe_ids(tmp_path: Path) -> None:
    assert validate(tmp_path, (FIXTURES / "name_edge_cases.fasta").read_text()) == ""
    names = (tmp_path / "names.tsv").read_text().splitlines()[1:]
    pairs = [n.split("\t")[:2] for n in names]
    assert pairs[:2] == [["s0001", "species_1"], ["s0002", "species.2"]]
    sanitized = (tmp_path / "sanitized.fasta").read_text()
    assert ">s0005\nMANEDVLL\n" in sanitized
    report = json.loads((tmp_path / "validation.json").read_text())
    assert report["identical_sequence_groups"] == [["spe_cies.4", "1_spe_cies.4"]]


def test_header_description_is_split_off(tmp_path: Path) -> None:
    fasta = "".join(f">seq{i} some description\nMKTAYIAK\n" for i in range(4))
    assert validate(tmp_path, fasta) == ""
    assert "seq0\tsome description" in (tmp_path / "names.tsv").read_text()


@pytest.mark.parametrize(
    ("fasta", "message"),
    [
        (">a\nMK\n>b\nMK\n>c\nMK\n", "at least 4"),
        (">a\nMK\n>a\nMT\n>b\nMK\n>c\nMK\n", "name 'a' (two different sequences: rename one)"),
        (">a\nMK1\n>b\nMK\n>c\nMK\n>d\nMK\n", "invalid characters: 1"),
        (">a\n\n>b\nMK\n>c\nMK\n>d\nMK\n", "is empty"),
        ("MK\n>a\nMK\n", "before the first"),
        (">a\nACGTACGT\n>b\nACGTACGA\n>c\nACGTTCGT\n>d\nACGAACGT\n", "--seqtype codon"),
        (">\nMK\n", "empty sequence name"),
    ],
)
def test_bad_input_rejected(tmp_path: Path, fasta: str, message: str) -> None:
    assert message in validate(tmp_path, fasta)


def test_guardrails(tmp_path: Path) -> None:
    fasta = "".join(f">s{i}\nMKTAYIAK\n" for i in range(6))
    assert "at most 5 sequences" in validate(tmp_path, fasta, max_sequences=5)
    assert "limit is 4" in validate(tmp_path, fasta, max_length=4)
    many = [FOUR] * 3
    assert "at most 2 provided alignments" in validate(tmp_path, FOUR, many, max_provided=2)


def aligned(rows: dict[str, str]) -> str:
    return "".join(f">{name}\n{row}\n" for name, row in rows.items())


ROWS = {"seq0": "MKTAYIAKLQ", "seq1": "MKTAYIAKVQ", "seq2": "MKTAYIAKLQ", "seq3": "MKTAYIAKVQ"}


def test_provided_alignment_is_renamed_and_kept_in_sequence_order(tmp_path: Path) -> None:
    gapped = {"seq3": "MKTAY-IAKVQ", "seq0": "MKTAY-IAKLQ", "seq1": "MKTAY-IAKVQ"}
    gapped["seq2"] = "mktay.iaklq"  # lowercase and '.' gaps: accepted, with a warning
    assert validate(tmp_path, FOUR, [aligned(gapped)], labels={"alignments": ["Mine"]}) == ""
    written = (tmp_path / "provided/alignments/provided1.fasta").read_text()
    assert written.splitlines()[:2] == [">s0001", "MKTAY-IAKLQ"]
    assert ">s0003\nMKTAY-IAKLQ" in written
    report = json.loads((tmp_path / "validation.json").read_text())
    assert "lowercase" in report["warnings"][0]
    sources = (tmp_path / "sources.tsv").read_text().splitlines()
    assert sources[1] == "provided1\talignment\tMine\taln0.fasta"


def test_inherited_alignments_keep_their_aligner_keys(tmp_path: Path) -> None:
    # A job extending another: its MAFFT alignment stays `mafft`; the user's own are
    # numbered among themselves, and only they count towards max_provided.
    rows = aligned({n: "MKTAY-IAK" + r[-2:] for n, r in ROWS.items()})
    assert validate(
        tmp_path, FOUR, [rows, rows, rows], keys=["mafft", "", "muscle"], max_provided=1,
        labels={"alignments": ["MAFFT", "Mine", "MUSCLE"]},
    ) == ""  # fmt: skip
    assert {p.name for p in (tmp_path / "provided/alignments").iterdir()} == {
        "mafft.fasta", "provided1.fasta", "muscle.fasta",
    }  # fmt: skip
    sources = (tmp_path / "sources.tsv").read_text().splitlines()
    assert sources[1:] == [
        "mafft\talignment\tMAFFT\taln0.fasta",
        "provided1\talignment\tMine\taln1.fasta",
        "muscle\talignment\tMUSCLE\taln2.fasta",
    ]


@pytest.mark.parametrize(
    ("keys", "message"),
    [
        (["mafft", "mafft"], "used twice"),
        (["mafft", "../x"], "unknown aligner key"),
        (["prank_codon", ""], "needs coding DNA"),
        (["", ""], "at most 1 provided alignments"),
    ],
)
def test_bad_alignment_keys_rejected(tmp_path: Path, keys: list[str], message: str) -> None:
    rows = aligned(ROWS)
    assert message in validate(tmp_path, FOUR, [rows, rows], keys=keys, max_provided=1)


def test_sequences_can_come_from_the_alignment(tmp_path: Path) -> None:
    gapped = {name: row[:5] + "-" + row[5:] for name, row in ROWS.items()}
    assert validate(tmp_path, None, [aligned(gapped)]) == ""
    assert ">s0001\nMKTAYIAKLQ\n" in (tmp_path / "sanitized.fasta").read_text()


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ({**ROWS, "seq3": "MKTAYIAKV"}, "different lengths"),
        ({k: v for k, v in ROWS.items() if k != "seq3"}, "missing 'seq3'"),
        ({**ROWS, "extra": "MKTAYIAKLQ"}, "not in the sequences: 'extra'"),
        ({**ROWS, "seq0": "MKTAYIAKLL"}, "differ from the sequences: 'seq0'"),
        ({**ROWS, "seq0": "MKTAYIAK*Q"}, "invalid characters"),
    ],
)
def test_bad_alignment_rejected(tmp_path: Path, rows: dict[str, str], message: str) -> None:
    assert message in validate(tmp_path, FOUR, [aligned(rows)])


def test_provided_tree_renamed_with_labels_stripped(tmp_path: Path) -> None:
    tree = "((seq0:0.1,seq1:0.2)95:0.3,seq2:0.1,seq3:0.4)root;"
    assert validate(tmp_path, FOUR, trees=[tree]) == ""
    written = (tmp_path / "provided/trees/provided1.nwk").read_text().strip()
    assert written == "((s0001:0.1,s0002:0.2):0.3,s0003:0.1,s0004:0.4);"
    assert "provided1\ttree\tt0.nwk\tt0.nwk" in (tmp_path / "sources.tsv").read_text()


def test_kept_branch_lengths_are_written_exactly(tmp_path: Path) -> None:
    """--tree_branch_lengths keep: IQ-TREE gets the lengths as given (a simulation's true
    tree), digit for digit; a rooted tree's root branch may have none."""
    tree = "((seq0:0.6618865364,seq1:0.0000001)9:0.25,(seq2:0.1,seq3:1e-3):0.75);"
    assert validate(tmp_path, FOUR, trees=[tree], keep_lengths=[True]) == ""
    written = (tmp_path / "provided/trees/provided1.nwk").read_text().strip()
    assert written == "((s0001:0.6618865364,s0002:1e-07):0.25,(s0003:0.1,s0004:0.001):0.75);"


@pytest.mark.parametrize(
    ("tree", "message"),
    [
        ("(seq0,seq1,(seq2,seq3));", "5 branches without a length"),
        ("(seq0:0.1,seq1:0.1,(seq2:0.1,seq3):0.2);", "1 branch without a length"),
        ("(seq0:0.1,seq1:-0.1,(seq2:0.1,seq3:0.1):0.2);", "1 negative branch length"),
    ],
)
def test_kept_branch_lengths_must_all_be_there(tmp_path: Path, tree: str, message: str) -> None:
    error = validate(tmp_path, FOUR, trees=["(seq0,seq1,(seq2,seq3));", tree],
                     keep_lengths=[False, True])  # fmt: skip
    assert f"tree 2: {message}" in error
    # Re-estimated lengths (the default) need none.
    assert validate(tmp_path, FOUR, trees=[tree.replace("-", "")]) == ""


def test_tree_names_are_exact(tmp_path: Path) -> None:
    fasta = "".join(f">Candida_{i}\nMKTAYIAK\n" for i in range(4))
    underscores = "(Candida_0,Candida_1,(Candida_2,Candida_3));"
    assert validate(tmp_path, fasta, trees=[underscores]) == ""
    wrong_case = "(candida_0,Candida_1,(Candida_2,Candida_3));"
    assert "missing 'Candida_0'" in validate(tmp_path, fasta, trees=[wrong_case])


def test_polytomy_accepted(tmp_path: Path) -> None:
    assert validate(tmp_path, FOUR, trees=["(seq0,seq1,seq2,seq3);"]) == ""


def test_quoted_names_and_comments_accepted(tmp_path: Path) -> None:
    fasta = FOUR.replace(">seq0", ">it's,odd")
    tree = "('it''s,odd':0.1[&note=x],seq1,(seq2,seq3)[comment, here]);"
    assert validate(tmp_path, fasta, trees=[tree]) == ""


@pytest.mark.parametrize(
    ("tree", "message"),
    [
        ("(seq0,seq1,(seq2,seq2));", "duplicate leaf names 'seq2'"),
        ("(seq0,seq1,seq2);", "missing 'seq3'"),
        ("(seq0,seq1,(seq2,seq3);", "not valid Newick"),
        ("(" * 50 + "seq0" + ")" * 50 + ";", "nesting is deeper"),
        ("(seq0,seq1,(seq2,seq3,));", "empty element"),
        ("(seq0,,seq1,(seq2,seq3));", "empty element"),
        ("(seq0,seq1,(seq2,seq3):0.1:0.2);", "more than one branch length"),
        ("(seq0,seq1,(seq2,seq3));(seq0,seq1,(seq2,seq3));", "single tree"),
        ("(seq0,seq1,(seq2,seq3))", "single tree"),
        ("(seq0,seq1,('seq2,x',seq3));", "missing 'seq2'"),
    ],
)
def test_bad_tree_rejected(tmp_path: Path, tree: str, message: str) -> None:
    assert message in validate(tmp_path, FOUR, trees=[tree])


def test_models_checked_and_candidates_written(tmp_path: Path) -> None:
    assert validate(tmp_path, FOUR) == ""
    assert (tmp_path / "modelfinder_args.txt").read_text().strip() == ""  # IQ-TREE default
    assert "unknown substitution model 'LG+X'" in validate(tmp_path, FOUR, models=["LG+X"])
    assert validate(tmp_path, FOUR, models=["LG+G4", "WAG+R4"]) == ""
    args = (tmp_path / "modelfinder_args.txt").read_text().strip()
    assert args == "-mset LG,WAG -mrate G,R -mfreq FU"
    # RAxML-NG lacks JTTDCMut: refused if ticked, and left out of the default candidates
    refused = validate(tmp_path, FOUR, models=["JTTDCMut+G4", "LG+G4"], raxml=True)
    assert "RAxML-NG does not have the JTTDCMut matrix" in refused
    assert validate(tmp_path, FOUR, raxml=True) == ""
    flags = (tmp_path / "modelfinder_args.txt").read_text().split()
    assert flags[0] == "-mset" and "JTTDCMut" not in flags[1] and len(flags) == 2


def test_mixtures_are_full_searches_only(tmp_path: Path) -> None:
    # A mixture is no ModelFinder candidate: left out when the models are the fallback...
    assert validate(tmp_path, FOUR, models=["LG+C20+G4", "LG+G4"]) == ""
    assert (tmp_path / "modelfinder_args.txt").read_text().split() == [
        "-mset", "LG", "-mrate", "G", "-mfreq", "FU"]  # fmt: skip
    # ...refused when named, and RAxML-NG has none.
    refused = validate(tmp_path, FOUR, models=["LG"], candidates=["LG+C20+G4"])
    assert "site-profile mixtures are not ModelFinder candidates" in refused
    refused = validate(tmp_path, FOUR, models=["LG+C20+F+G4"], raxml=True)
    assert "RAxML-NG has no site-profile mixtures" in refused


def test_a_mixture_report_gives_the_catalog_name() -> None:
    report = """Input data: 7 sequences with 400 amino-acid sites
Mixture model of substitution: MIX{LG+F+F,LG+FC20pi1,LG+FC20pi2}+G4
Log-likelihood of the tree: -100.5 (s.e. 1.0)
Number of free parameters (#branches + #model parameters): 31
Akaike information criterion (AIC) score: 263.0
Corrected Akaike information criterion (AICc) score: 264.0
Bayesian information criterion (BIC) score: 387.0
Total tree length (sum of branch lengths): 1.5
"""
    assert parse_fit(report, "LG+C20+F+G4")["model"] == "LG+C20+F+G4"
    with pytest.raises(ValueError, match="no model"):
        parse_fit(report)


def test_every_step_that_fits_a_mixture_fits_its_weights() -> None:
    """The catalog and mixtures.nf agree, and each IQ-TREE step fitting a fixed model passes
    -mwopt: --ancestral and the topology test refit the model, so one left out would silently
    use C20's published weights."""
    rule = (PIPELINE / "modules/mixtures.nf").read_text()
    pattern = re.search(r"model ==~ /(.+)/", rule)
    assert pattern
    for model in catalog.all_fixed_models() + catalog.all_codon_models():
        assert bool(re.fullmatch(pattern[1], model)) == catalog.is_mixture(model), model
    for module in ("tree_iqtree", "fit_iqtree", "ancestral_iqtree", "topology_test"):
        text = (PIPELINE / f"modules/{module}.nf").read_text()
        assert "mixtureFlag(meta)" in text and "${weights} " in text, module


def test_main_accepts_every_catalog_model() -> None:
    pattern = r"[A-Za-z][A-Za-z0-9.]*(\+[A-Z][0-9]{0,2})*"
    assert pattern in (PIPELINE / "main.nf").read_text()
    assert all(re.fullmatch(pattern, m) for m in catalog.all_fixed_models())


def test_candidates_apart_from_full_search_models(tmp_path: Path) -> None:
    # ModelFinder's candidates come from their own list, not from the full-search models...
    assert validate(tmp_path, FOUR, models=["JTT"], candidates=["LG+G4", "WAG+R4"]) == ""
    args = (tmp_path / "modelfinder_args.txt").read_text().strip()
    assert args == "-mset LG,WAG -mrate G,R -mfreq FU"
    # ...both are checked...
    assert "unknown substitution model" in validate(tmp_path, FOUR, candidates=["LG+X"])
    refused = validate(tmp_path, FOUR, models=["LG"], candidates=["JTTDCMut+G4"], raxml=True)
    assert "RAxML-NG does not have the JTTDCMut matrix" in refused
    # ...and the whole grid means IQ-TREE's default set (no flags).
    assert validate(tmp_path, FOUR, models=["JTT"], candidates=catalog.all_models()) == ""
    assert (tmp_path / "modelfinder_args.txt").read_text().strip() == ""


def test_labels_checked(tmp_path: Path) -> None:
    same = {"alignments": ["Mine", "mine"]}
    assert "used twice" in validate(tmp_path, FOUR, [aligned(ROWS)] * 2, labels=same)
    control = {"alignments": ["bad\tlabel"]}
    assert "control characters" in validate(tmp_path, FOUR, [aligned(ROWS)], labels=control)
    hostile = {"alignments": ["it's $(rm -rf /); ok"]}
    assert validate(tmp_path, FOUR, [aligned(ROWS)], labels=hostile) == ""


def test_support_labels_split_into_measures(tmp_path: Path) -> None:
    from phylobot.pipeline.export import write_supports

    out = tmp_path / "supports.tsv"
    write_supports(out, [("Node1", ""), ("Node2", "85.5/97")], ("sh_alrt", "ufboot"))
    assert out.read_text() == "node\tufboot\tsh_alrt\nNode1\t\t\nNode2\t97\t85.5\n"
    write_supports(out, [("Node2", "85.5")], ("sh_alrt",))
    assert out.read_text() == "node\tufboot\tsh_alrt\nNode2\t\t85.5\n"
    write_supports(out, [("Node2", "")], ())  # no support measure ran
    assert out.read_text() == "node\tufboot\tsh_alrt\nNode2\t\t\n"
    with pytest.raises(ValueError, match="does not match"):  # labels and request disagree
        write_supports(out, [("Node2", "85.5/97")], ("sh_alrt",))


def test_supports_matched_to_node_names_by_split() -> None:
    # Same topology, rooted differently: matching must not depend on rooting.
    named = "(s1:1,(s2:1,(s3:1,s4:1)Node3:1)Node2:1,(s5:1,s6:1)Node4:1)Node1;"
    supported = "((s3:1,s4:1)100:1,s2:1,(s1:1,(s5:1,s6:1)96:1)55:1);"
    assert supports_by_node(named, supported) == [
        ("Node3", "100"),
        ("Node2", "55"),
        ("Node4", "96"),
        ("Node1", ""),
    ]


def test_quote_newick_names() -> None:
    assert quote("Candida.albicans") == "Candida.albicans"
    assert quote("has space") == "'has space'"
    assert quote("it's(odd)") == "'it''s(odd)'"


FIXTURE_REPORTS = Path(__file__).parent / "fixtures" / "iqtree"


def test_parse_fit_from_fixed_model_report() -> None:
    fit = parse_fit((FIXTURE_REPORTS / "fixed_LG+G4.iqtree").read_text())
    assert fit == {
        "model": "LG+G4",
        "log_likelihood": "-4412.9327",
        "free_parameters": "12",
        "aic": "8849.8655",
        "aicc": "8850.7941",
        "bic": "8896.1263",
        "sites": "349",
        "tree_length": "4.8148",
    }
    assert parse_modelfinder((FIXTURE_REPORTS / "fixed_LG+G4.iqtree").read_text()) == []


def test_parse_modelfinder_report() -> None:
    report = (FIXTURE_REPORTS / "modelfinder.iqtree").read_text()
    ranking = parse_modelfinder(report)
    assert ranking[0]["model"] == "LG+G4" and ranking[0]["bic_in_95"] == "yes"
    assert ranking[0]["w_bic"] == "0.973" and ranking[1]["bic_in_95"] == "no"
    assert len(ranking) > 10
    # The final-tree section comes after model selection; the final fit wins.
    assert parse_fit(report)["bic"] == "8896.1607" and parse_fit(report)["sites"] == "350"


def summarize(
    tmp_path: Path,
    fits: list[str],
    primary: str,
    sources: str = "",
    alignments: str = "",
    stop_after: str = "ancestors",
) -> dict[str, Any]:
    (tmp_path / "sources.tsv").write_text("key\tkind\tlabel\tfilename\n" + sources)
    (tmp_path / "fits.tsv").write_text("\n".join(fits) + "\n")
    command = [*module("summary"), "sources.tsv", "fits.tsv", primary, alignments, stop_after]
    subprocess.run(command, cwd=tmp_path, check=True)
    result: dict[str, Any] = json.loads((tmp_path / "summary.json").read_text())
    return result


def fit_row(aligner: str, tree: str, key: str, model: str, bic: float) -> str:
    """aligner, tree, model key, then model, logL, params, AIC, AICc, BIC, sites, length."""
    return "\t".join([aligner, tree, key, model, "-4400", "12", str(bic - 40), str(bic - 39),
                      str(bic), "349", "4.8"])  # fmt: skip


def test_summary_marks_primary_and_labels(tmp_path: Path) -> None:
    fits = [fit_row("mafft", "iqtree", "modelfinder", "LG+G4", 8896.1),
            fit_row("mafft", "iqtree", "LG+G4", "LG+G4", 8896.1)]  # fmt: skip
    summary = summarize(tmp_path, fits, "mafft/iqtree/modelfinder")
    assert summary["primary"] == "mafft/iqtree/modelfinder"
    keys = [c["key"] for c in summary["combinations"]]
    assert keys == ["mafft/iqtree/LG+G4", "mafft/iqtree/modelfinder"]
    labels = [c["label"] for c in summary["combinations"]]
    assert labels == ["MAFFT × IQ-TREE × LG+G4", "MAFFT × IQ-TREE × ModelFinder"]
    assert summary["combinations"][1]["fit"]["bic"] == "8896.1"


def test_summary_falls_back_when_primary_missing(tmp_path: Path) -> None:
    sources = "provided1\talignment\tMy alignment\tmine.fasta\n"
    fits = [fit_row("provided1", "iqtree", "LG", "LG", 9000)]
    summary = summarize(tmp_path, fits, "mafft/iqtree/LG", sources=sources)
    assert summary["primary"] == "provided1/iqtree/LG"
    assert summary["combinations"][0]["label"] == "My alignment × IQ-TREE × LG"


def test_summary_of_a_job_that_stops_after_alignments(tmp_path: Path) -> None:
    sources = "provided1\talignment\tCurated\tmine.fasta\n"
    summary = summarize(tmp_path, [], "", sources, "mafft,provided1", "alignments")
    assert summary["stop_after"] == "alignments" and summary["combinations"] == []
    assert summary["primary"] == ""
    assert summary["alignments"] == [
        {"key": "mafft", "label": "MAFFT", "path": "mafft/alignment.fasta"},
        {"key": "provided1", "label": "Curated", "path": "provided1/alignment.fasta"},
    ]
    assert (tmp_path / "model_comparison.tsv").read_text().count("\n") == 1  # header only


def test_model_comparison_ranks_within_each_alignment(tmp_path: Path) -> None:
    fits = [
        fit_row("mafft", "iqtree", "modelfinder", "LG+G4", 8896.0),  # same fit as LG+G4 below
        fit_row("mafft", "iqtree", "LG+G4", "LG+G4", 8896.0),
        fit_row("mafft", "iqtree", "WAG+R4", "WAG+R4", 8900.0),
        fit_row("provided1", "iqtree", "LG+G4", "LG+G4", 1.0),  # other data: not compared
    ]
    summarize(tmp_path, fits, "mafft/iqtree/modelfinder")
    lines = (tmp_path / "model_comparison.tsv").read_text().splitlines()
    rows = [dict(zip(lines[0].split("\t"), line.split("\t"), strict=True)) for line in lines[1:]]
    mafft = [r for r in rows if r["aligner"] == "mafft"]
    assert [r["model_key"] for r in mafft] == ["LG+G4", "modelfinder", "WAG+R4"]
    assert [r["rank_bic"] for r in mafft] == ["1", "1", "2"]  # duplicate fit shares a rank
    assert [r["delta_bic"] for r in mafft] == ["0.0000", "0.0000", "4.0000"]
    # weights over distinct fits: exp(0) and exp(-2), normalised
    assert float(mafft[0]["w_bic"]) == float(mafft[1]["w_bic"]) == pytest.approx(0.8808, abs=1e-4)
    provided = [r for r in rows if r["aligner"] == "provided1"]
    assert provided[0]["rank_bic"] == "1" and provided[0]["w_bic"] == "1"


def test_export_alignment_and_tree_for_early_stops(tmp_path: Path) -> None:
    from phylobot.pipeline.export import export_alignment, export_tree

    (tmp_path / "names.tsv").write_text(
        "id\tname\tdescription\n" + "".join(f"s000{i}\tseq {i}\t\n" for i in range(1, 6))
    )
    (tmp_path / "aln.fasta").write_text("".join(f">s000{i}\nMK-L\n" for i in range(1, 6)))
    export_alignment(tmp_path / "names.tsv", tmp_path / "aln.fasta", tmp_path / "a")
    assert (tmp_path / "a" / "alignment.fasta").read_text().startswith(">seq 1\nMK-L\n")

    # IQ-TREE writes UFBoot supports as internal labels; export names the nodes instead
    (tmp_path / "tree.nwk").write_text(
        "(s0001:0.1,(s0002:0.2,s0003:0.2)95:0.1,(s0004:0.1,s0005:0.1)80:0.3);\n"
    )
    (tmp_path / "model.txt").write_text("LG+G4\n")
    (tmp_path / "outgroup.txt").write_text("s0001\n")
    export_tree(tmp_path / "names.tsv", tmp_path / "aln.fasta", tmp_path / "tree.nwk",
                tmp_path / "model.txt", FIXTURE_REPORTS / "fixed_LG+G4.iqtree",
                tmp_path / "outgroup.txt", tmp_path / "t")  # fmt: skip
    out = tmp_path / "t"
    tree = (out / "tree.nwk").read_text()
    assert "Node1" in tree and "Node3" in tree and "95" not in tree and "'seq 2'" in tree
    supports = [line.split("\t") for line in (out / "supports.tsv").read_text().splitlines()]
    assert sorted(row[1] for row in supports) == ["", "80", "95", "ufboot"]
    assert {row[2] for row in supports} == {"sh_alrt", ""}
    assert (out / "tree.rooted.nwk").is_file()
    assert json.loads((out / "rooting.json").read_text())["outgroup"] == ["seq 1"]
    assert (out / "fit.tsv").is_file() and (out / "model.txt").read_text() == "LG+G4\n"
    assert not (out / "modelfinder.tsv").exists()  # a fixed model: no ranking

    # A ModelFinder combination: the ranking comes from MODELFINDER's own report
    export_tree(tmp_path / "names.tsv", tmp_path / "aln.fasta", tmp_path / "tree.nwk",
                tmp_path / "model.txt", FIXTURE_REPORTS / "fixed_LG+G4.iqtree",
                tmp_path / "outgroup.txt", tmp_path / "mf", ("ufboot",),
                FIXTURE_REPORTS / "modelfinder.iqtree")  # fmt: skip
    ranking = (tmp_path / "mf" / "modelfinder.tsv").read_text().splitlines()
    assert ranking[0].startswith("model\t") and len(ranking) > 2
    assert not (out / "ancestors.fasta").exists()
    # no outgroup: no rooted tree
    (tmp_path / "outgroup.txt").write_text("")
    export_tree(tmp_path / "names.tsv", tmp_path / "aln.fasta", tmp_path / "tree.nwk",
                tmp_path / "model.txt", FIXTURE_REPORTS / "fixed_LG+G4.iqtree",
                tmp_path / "outgroup.txt", tmp_path / "u")  # fmt: skip
    assert not (tmp_path / "u" / "tree.rooted.nwk").exists()


# Coding DNA (--seqtype codon)

CDS = {  # translations MKWL, MKWV, MRWL, MKFL; seq0 ends with a stop codon
    "seq0": "ATGAAATGGCTGTAA",
    "seq1": "ATGAAGTGGGTG",
    "seq2": "ATGCGTTGGCTT",
    "seq3": "ATGAAATTTCTG",
}


def fasta(records: dict[str, str]) -> str:
    return "".join(f">{name}\n{seq}\n" for name, seq in records.items())


def test_coding_sequences_translated_and_kept(tmp_path: Path) -> None:
    assert validate(tmp_path, fasta(CDS), seqtype="codon") == ""
    assert (tmp_path / "sanitized.fasta").read_text().split()[1::2] == [
        "MKWL", "MKWV", "MRWL", "MKFL",
    ]  # fmt: skip
    cds = (tmp_path / "sanitized.cds.fasta").read_text().split()
    assert cds[:2] == [">s0001", "ATGAAATGGCTG"]  # final stop codon removed
    report = json.loads((tmp_path / "validation.json").read_text())
    assert report["seqtype"] == "codon"
    assert report["genetic_code"] == {"ncbi": 1, "name": "Standard"}
    assert report["max_length"] == 4  # codons
    assert report["warnings"] == ["final stop codon removed: 1 sequence ('seq0')"]


def test_protein_jobs_report_their_type_and_write_no_codon_files(tmp_path: Path) -> None:
    assert validate(tmp_path, FOUR) == ""
    assert json.loads((tmp_path / "validation.json").read_text())["seqtype"] == "protein"
    assert not (tmp_path / "sanitized.cds.fasta").exists()


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"seq1": "ATGAAGTGGGT"}, "'seq1': length 11 is not a multiple of 3 (a frameshift"),
        ({"seq2": "ATGTGATGGCTT"}, "'seq2': internal stop codon TGA at codon 2"),
        ({"seq3": "ATGAAJTTTCTG"}, "not nucleotide (IUPAC) codes: J"),
        ({"seq3": "AUGAAATTTCTG"}, "not nucleotide (IUPAC) codes: U"),
        ({"seq0": "TAA"}, "no codons besides a stop codon"),
    ],
)
def test_bad_coding_sequences_rejected(
    tmp_path: Path, change: dict[str, str], message: str
) -> None:
    assert message in validate(tmp_path, fasta(CDS | change), seqtype="codon")


def test_ambiguity_codes_accepted_and_reported(tmp_path: Path) -> None:
    # AAR is K either way; AMC (AAC or ACC) is N or T, so X; lowercase is fine.
    ambiguous = CDS | {"seq1": "ATGAARTGGGTG", "seq3": "ATGamcTTTCTG"}
    assert validate(tmp_path, fasta(ambiguous), seqtype="codon") == ""
    assert (tmp_path / "sanitized.fasta").read_text().split()[3::2] == ["MKWV", "MRWL", "MXFL"]
    assert (tmp_path / "sanitized.cds.fasta").read_text().split()[-1] == "ATGAMCTTTCTG"
    report = json.loads((tmp_path / "validation.json").read_text())
    found = report["ambiguous_codons"]
    assert (found["codons"], found["sequences"], found["unknown"], found["translated"]) == (
        2, 2, 0, 1,
    )  # fmt: skip
    assert found["listed"] == [
        {"sequence": "seq1", "number": 2, "codon": "AAR", "aminos": "K", "translation": "K"},
        {"sequence": "seq3", "number": 2, "codon": "AMC", "aminos": "NT", "translation": "X"},
    ]
    assert report["warnings"][-1].startswith("ambiguous bases: 2 codons in 2 sequences")
    # The same codons from an alignment of them: the same fingerprint (the New Job page's
    # confirmation holds whichever way the sequences came).
    rows = {n: s if n == "seq0" else s + "---" for n, s in ambiguous.items()}
    other = tmp_path / "aln"
    other.mkdir()
    assert validate(other, None, [fasta(rows)], seqtype="codon") == ""
    again = json.loads((other / "validation.json").read_text())["ambiguous_codons"]
    assert again["fingerprint"] == found["fingerprint"]
    # Another genetic code translates differently: confirm again.
    third = tmp_path / "mt"
    third.mkdir()
    assert validate(third, fasta(ambiguous), seqtype="codon", genetic_code=2) == ""
    mt = json.loads((third / "validation.json").read_text())["ambiguous_codons"]
    assert mt["fingerprint"] != found["fingerprint"]
    plain = tmp_path / "plain"
    plain.mkdir()
    assert validate(plain, fasta(CDS), seqtype="codon") == ""
    assert "ambiguous_codons" not in json.loads((plain / "validation.json").read_text())


def test_coding_input_checks(tmp_path: Path) -> None:
    assert "look like protein" in validate(tmp_path, FOUR, seqtype="codon")
    assert "unknown sequence type" in validate(tmp_path, fasta(CDS), seqtype="rna")
    assert "unknown genetic code" in validate(tmp_path, fasta(CDS), seqtype="codon", genetic_code=7)
    assert "is 4 codons; limit is 3" in validate(
        tmp_path, fasta(CDS), seqtype="codon", max_length=3
    )


def test_genetic_code_decides_what_is_a_stop(tmp_path: Path) -> None:
    mito = CDS | {"seq2": "ATGTGATGGCTT"}  # TGA: stop in the standard code, W in mt
    assert validate(tmp_path, fasta(mito), seqtype="codon", genetic_code=2) == ""
    assert (tmp_path / "sanitized.fasta").read_text().split()[5] == "MWWL"


CODON_ROWS = {
    "seq0": "ATGAAATGGCTGTAA",
    "seq1": "ATGAAGTGGGTG---",
    "seq2": "ATGCGTTGGCTT---",
    "seq3": "ATGAAATTTCTG---",
}


def test_provided_codon_alignment_kept_and_translated(tmp_path: Path) -> None:
    gapped = {
        "seq0": "ATG---AAATGGCTG",
        "seq1": "ATG---AAGTGGGTG",
        "seq2": "ATG---CGTTGGCTT",
        "seq3": "---ATGAAATTTCTG",
    }
    assert validate(tmp_path, fasta(CDS), [fasta(gapped)], seqtype="codon") == ""
    protein = (tmp_path / "provided/alignments/provided1.fasta").read_text().split()
    assert protein[1::2] == ["M-KWL", "M-KWV", "M-RWL", "-MKFL"]
    codon = (tmp_path / "provided/codon_alignments/provided1.fasta").read_text().split()
    assert codon[1] == "ATG---AAATGGCTG"


def test_codon_alignment_final_stops_trimmed_and_sequences_taken_from_it(
    tmp_path: Path,
) -> None:
    assert validate(tmp_path, None, [fasta(CODON_ROWS)], seqtype="codon") == ""
    codon = (tmp_path / "provided/codon_alignments/provided1.fasta").read_text().split()
    assert codon[1] == "ATGAAATGGCTG"  # the stop column, left all gaps, is dropped
    assert (tmp_path / "sanitized.cds.fasta").read_text().split()[1] == "ATGAAATGGCTG"


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"seq1": "ATGAA-GTGGGTG--"}, "gap inside codon 2"),
        ({"seq1": "ATGTAATGGGTG---"}, "internal stop codon TAA in codon column 2"),
        ({"seq1": "ATGAAATGGGTG---"}, "rows differ from the sequences: 'seq1'"),
    ],
)
def test_bad_codon_alignment_rejected(tmp_path: Path, change: dict[str, str], message: str) -> None:
    assert message in validate(tmp_path, fasta(CDS), [fasta(CODON_ROWS | change)], seqtype="codon")


def codon_alignment(tmp_path: Path, form: str, rows: dict[str, str]) -> str:
    (tmp_path / "aligned.fasta").write_text(fasta(rows))
    (tmp_path / "cds.fasta").write_text(fasta({n: s[:12] for n, s in CDS.items()}))
    proc = subprocess.run(
        [*module("codon_alignment"), "mafft", form, "aligned.fasta", "cds.fasta", "1", "out"],
        cwd=tmp_path, capture_output=True, text=True,
    )  # fmt: skip
    return proc.stderr if proc.returncode else ""


def test_protein_alignment_back_translated(tmp_path: Path) -> None:
    rows = {"seq0": "MKWL-", "seq1": "MKWV-", "seq2": "MRWL-", "seq3": "-MKFL"}
    assert codon_alignment(tmp_path, "protein", rows) == ""
    codon = (tmp_path / "out/alignment.cds.fasta").read_text().split()
    assert codon[1] == "ATGAAATGGCTG---" and codon[7] == "---ATGAAATTTCTG"
    assert (tmp_path / "out/alignment.fasta").read_text().split()[1::2] == list(rows.values())


def test_codon_aligner_output_checked_against_the_input(tmp_path: Path) -> None:
    rows = {n: s[:12] for n, s in CDS.items()}
    assert codon_alignment(tmp_path, "codon", rows) == ""
    assert (tmp_path / "out/alignment.fasta").read_text().split()[1] == "MKWL"
    # PRANK -codon replaces a stop codon with NNN and exits 0: caught here.
    changed = rows | {"seq2": "ATGNNNTGGCTT"}
    assert "seq2: the aligned sequence differs" in codon_alignment(tmp_path, "codon", changed)
    # ...but where the input codon has an ambiguity code, the user's codon comes back.
    (tmp_path / "cds.fasta").write_text(fasta({n: s[:12] for n, s in CDS.items()}
                                              | {"seq2": "ATGCGYTGGCTT"}))  # fmt: skip
    (tmp_path / "aligned.fasta").write_text(fasta(changed))
    proc = subprocess.run(
        [*module("codon_alignment"), "prank_codon", "codon", "aligned.fasta", "cds.fasta",
         "1", "out"], cwd=tmp_path, capture_output=True, text=True,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "out/alignment.cds.fasta").read_text().split()[5] == "ATGCGYTGGCTT"
    assert (tmp_path / "out/alignment.fasta").read_text().split()[5] == "MRWL"  # CGY: R
    altered = {"seq0": "MKWL", "seq1": "MKWV", "seq2": "MKWL", "seq3": "MKFL"}
    assert "residue K but codon 2 is CGT" in codon_alignment(tmp_path, "protein", altered)
    (tmp_path / "empty.fasta").write_text("")
    proc = subprocess.run(
        [*module("codon_alignment"), "prank_codon", "codon", "empty.fasta", "cds.fasta", "1",
         "out"], cwd=tmp_path, capture_output=True, text=True,
    )  # fmt: skip
    assert "prank_codon wrote no alignment" in proc.stderr


def test_summary_ranks_codon_and_amino_acid_models_apart(tmp_path: Path) -> None:
    fits = [fit_row("mafft", "iqtree", "modelfinder", "LG+G4", 6000),
            fit_row("mafft", "iqtree", "WAG+G4", "WAG+G4", 6010),
            fit_row("mafft", "iqtree", "GY+F3X4+G4", "GY+F3X4+G4", 9000)]  # fmt: skip
    summary = summarize(tmp_path, fits, "mafft/iqtree/modelfinder")
    ranks = {c["model"]: (c["data"], c["rank_bic"]) for c in summary["combinations"]}
    assert ranks == {"modelfinder": ("amino", 1), "WAG+G4": ("amino", 2),
                     "GY+F3X4+G4": ("codon", 1)}  # fmt: skip
    rows = (tmp_path / "model_comparison.tsv").read_text().splitlines()
    assert rows[0].endswith("\tdata") and rows[-1].endswith("\tcodon")


def test_summary_labels_codon_fits_on_amino_acid_trees(tmp_path: Path) -> None:
    fits = [fit_row("mafft", "iqtree", "LG+G4", "LG+G4", 6000),
            fit_row("mafft", "iqtree", "modelfinder", "LG+G4", 6000),
            fit_row("mafft", "iqtree", "WAG", "WAG", 6100),
            fit_row("mafft", "iqtree_aa", "GY+F3X4+G4", "GY+F3X4+G4", 9000),
            fit_row("mafft", "raxml_aa-modelfinder", "GY+F3X4+G4", "GY+F3X4+G4", 9001),
            fit_row("mafft", "iqtree_aa-Q.pfam+G4", "GY+F3X4+G4", "GY+F3X4+G4", 9002)]  # fmt: skip
    summary = summarize(tmp_path, fits, "mafft/iqtree/LG+G4")
    labels = {c["tree"]: c["tree_label"] for c in summary["combinations"]}
    assert labels == {
        "iqtree": "IQ-TREE",
        "iqtree_aa": "IQ-TREE amino-acid tree (best: LG+G4)",  # ties: first model key
        "raxml_aa-modelfinder": "RAxML-NG amino-acid tree (ModelFinder)",
        "iqtree_aa-Q.pfam+G4": "IQ-TREE amino-acid tree (Q.pfam+G4)",
    }
    assert [c["data"] for c in summary["combinations"]].count("codon") == 3


def test_codon_models_checked(tmp_path: Path) -> None:
    assert validate(tmp_path, fasta(CDS), seqtype="codon", codon_models=["GY+F3X4+G4"]) == ""
    assert "unknown codon model 'MG+F'" in validate(
        tmp_path, fasta(CDS), seqtype="codon", codon_models=["MG+F"]
    )
    assert "codon models need coding DNA" in validate(tmp_path, FOUR, codon_models=["GY+F3X4"])


def test_export_sums_codon_posteriors_into_amino_acids(tmp_path: Path) -> None:
    from phylobot.codons import CODES, IQTREE_AMINO_ORDER
    from phylobot.pipeline.export import export

    sense = [c for c in codons_in_order() if CODES[1].table[c] != "*"]
    rows = []
    for node in ("Node1", "Node2"):
        # site 1: ATG certain; site 2: CTT 0.4 but ATT 0.35 + ATC 0.25 (I wins)
        for site, probs in ((1, {"ATG": 1.0}), (2, {"CTT": 0.4, "ATT": 0.35, "ATC": 0.25})):
            best = max(probs, key=lambda c: probs[c])
            values = "\t".join(f"{probs.get(c, 0.0):.5f}" for c in sense)
            rows.append(f"{node}\t{site}\t{best}\t{values}")
    header = "Node\tSite\tState\t" + "\t".join(f"p_{c}" for c in sense)
    (tmp_path / "asr.state").write_text("# IQ-TREE\n" + header + "\n" + "\n".join(rows) + "\n")
    (tmp_path / "names.tsv").write_text(
        "id\tname\tdescription\n" + "".join(f"s000{i}\tseq{i}\t\n" for i in range(1, 5))
    )
    (tmp_path / "aln.fasta").write_text("".join(f">s000{i}\nMI\n" for i in range(1, 5)))
    (tmp_path / "cds.fasta").write_text("".join(f">s000{i}\nATGATT\n" for i in range(1, 5)))
    tree = "((s0001,s0002)Node2,s0003,s0004)Node1;"
    (tmp_path / "asr.treefile").write_text(tree + "\n")
    (tmp_path / "support.treefile").write_text("((s0001,s0002),s0003,s0004);\n")
    (tmp_path / "model.txt").write_text("GY+F3X4+G4\n")
    report = (FIXTURE_REPORTS / "fixed_LG+G4.iqtree").read_text()
    (tmp_path / "fit.iqtree").write_text(report)
    out = tmp_path / "out"
    export(
        tmp_path / "names.tsv", tmp_path / "aln.fasta", tmp_path / "asr.treefile",
        tmp_path / "support.treefile", tmp_path / "asr.state", tmp_path / "model.txt",
        tmp_path / "fit.iqtree", out, ("parsimony",), support_kinds=(),
        codon_alignment=tmp_path / "cds.fasta", genetic_code=1,
    )  # fmt: skip
    assert (out / "ancestors.fasta").read_text().split() == [">Node1", "MI", ">Node2", "MI"]
    assert (out / "ancestors.cds.fasta").read_text().split()[1] == "ATGATT"
    posteriors = (out / "ancestors_posteriors.tsv").read_text().splitlines()
    assert posteriors[0].split("\t")[2:] == list(IQTREE_AMINO_ORDER)
    pp = (out / "ancestors_pp.tsv").read_text().splitlines()
    assert pp[2] == "Node1\t2\tI\t0.60000"
    codon_pp = (out / "ancestors_codon_pp.tsv").read_text().splitlines()
    assert codon_pp[2] == "Node1\t2\tATT\t0.35000"
    assert (out / "alignment.cds.fasta").read_text().split()[:2] == [">seq1", "ATGATT"]
    assert (out / "ancestors.gaps-parsimony.cds.fasta").exists()
    codon_rows = (out / "ancestors_codon_posteriors.tsv").read_text().splitlines()
    assert codon_rows[0] == "node\tsite\tcodon\tp"
    assert codon_rows[2:5] == ["Node1\t2\tCTT\t0.40000", "Node1\t2\tATT\t0.35000",
                               "Node1\t2\tATC\t0.25000"]  # fmt: skip


def codons_in_order() -> list[str]:
    from phylobot.codons import CODONS

    return sorted(CODONS)  # IQ-TREE's codon columns: alphabetical (AAA, AAC, ...)


def test_codon_modelfinder_ranks_with_codon_models(tmp_path: Path) -> None:
    fits = [fit_row("mafft", "iqtree", "modelfinder", "LG+G4", 6000),
            fit_row("mafft", "iqtree", "modelfinder_codon", "MG+F1X4+G4", 9100),
            fit_row("mafft", "iqtree", "GY+F3X4+G4", "GY+F3X4+G4", 9200)]  # fmt: skip
    summary = summarize(tmp_path, fits, "mafft/iqtree/modelfinder")
    ranks = {c["model"]: (c["data"], c["rank_bic"]) for c in summary["combinations"]}
    assert ranks["modelfinder_codon"] == ("codon", 1) and ranks["GY+F3X4+G4"] == ("codon", 2)
    labels = {c["model"]: c["model_label"] for c in summary["combinations"]}
    assert labels["modelfinder_codon"] == "ModelFinder (codon)"


def test_codon_modelfinder_flags_written(tmp_path: Path) -> None:
    ticked = ["MG+F1X4+G4", "GY+F3X4"]
    (tmp_path / "input.fasta").write_text(fasta(CDS))
    manifest = {"sequences": True, **LIMITS, "seqtype": "codon",
                "codon_modelfinder": True, "codon_candidates": ticked}  # fmt: skip
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    subprocess.run([*module("validate"), "manifest.json"], cwd=tmp_path, check=True)
    flags = (tmp_path / "codon_modelfinder_args.txt").read_text().split()
    assert flags == ["-mset", "MG,GY", "-mrate", "E,G", "-mfreq", "F3X4,F1X4"]
    manifest["seqtype"] = "protein"
    (tmp_path / "input.fasta").write_text(FOUR)
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    proc = subprocess.run([*module("validate"), "manifest.json"], cwd=tmp_path,
                          capture_output=True, text=True)  # fmt: skip
    assert "codon models need coding DNA" in proc.stderr


def test_codeml_models_checked_and_controls_written(tmp_path: Path) -> None:
    def run(seqtype: str, code: int, models: list[str], records: str) -> str:
        (tmp_path / "input.fasta").write_text(records)
        manifest = {"sequences": True, **LIMITS, "seqtype": seqtype, "genetic_code": code,
                    "codeml_models": models}  # fmt: skip
        (tmp_path / "manifest.json").write_text(json.dumps(manifest))
        proc = subprocess.run([*module("validate"), "manifest.json"], cwd=tmp_path,
                              capture_output=True, text=True)  # fmt: skip
        return proc.stderr

    assert run("codon", 1, ["M0", "M8"], fasta(CDS)) == ""
    assert (tmp_path / "codeml" / "M8" / "start2.ctl").exists()
    assert "icode = 0" in (tmp_path / "codeml" / "M0" / "start1.ctl").read_text()
    assert "unknown codeml models 'M9'" in run("codon", 1, ["M9"], fasta(CDS))
    assert "codon models need coding DNA" in run("protein", 1, ["M0"], FOUR)
    assert "PAML does not have genetic code 25" in run("codon", 25, ["M0"], fasta(CDS))


def test_site_rates_and_log_likelihoods(tmp_path: Path) -> None:
    from phylobot.pipeline.export import write_site_statistics

    # IQ-TREE 3.1.3's -wsr and -wsl files (lysin, Q.PFAM+G4; the first columns).
    rates = tmp_path / "asr.rate"
    rates.write_text(
        "# Site-specific subtitution rates determined by empirical Bayesian method\n"
        "#   Site:   Alignment site ID\nSite\tRate\tCat\tC_Rate\n"
        "1\t0.22921\t1\t0.10717\n2\t1.71435\t4\t2.49342\n"
    )
    site_lh = tmp_path / "asr.sitelh"
    site_lh.write_text("1 2\nSite_Lh    -4.22469 -33.444\n")
    write_site_statistics(rates, site_lh, tmp_path)
    assert (tmp_path / "site_rates.tsv").read_text() == (
        "site\trate\tcategory\tcategory_rate\n1\t0.22921\t1\t0.10717\n2\t1.71435\t4\t2.49342\n"
    )
    assert (tmp_path / "site_lh.tsv").read_text() == (
        "site\tlog_likelihood\n1\t-4.22469\n2\t-33.444\n"
    )
    # Without rate heterogeneity IQ-TREE writes the header only: no site_rates.tsv.
    (tmp_path / "site_rates.tsv").unlink()
    rates.write_text("# comments\nSite\tRate\tCat\tC_Rate\n")
    write_site_statistics(rates, site_lh, tmp_path)
    assert not (tmp_path / "site_rates.tsv").exists()


def test_index_by_node_finds_each_nodes_rows(tmp_path: Path) -> None:
    """EXPORT's byte index: each node's rows are ranges of the table itself."""
    source = tmp_path / "ancestors_posteriors.tsv"
    source.write_text(
        "node\tsite\tM\tK\n"
        "Node1\t1\t0.9\t0.1\nNode1\t2\t0.8\t0.2\n"
        "Node2\t1\t0.5\t0.5\n"
        "Node1\t3\t0.7\t0.3\n"  # out of order: a second run, never lost
    )
    assert export.index_by_node(source) == 2
    index = ancestors.read_index((tmp_path / "ancestors_posteriors.index.tsv").read_text())
    assert len(index.nodes["Node1"]) == 2 and len(index.nodes["Node2"]) == 1
    node1 = index.text("Node1", ancestors.file_range(source))
    assert node1 is not None and node1.splitlines()[0] == "node\tsite\tM\tK"  # with header
    assert [row.split("\t")[1] for row in node1.splitlines()[1:]] == ["1", "2", "3"]
    node2 = index.text("Node2", ancestors.file_range(source))
    assert node2 == "node\tsite\tM\tK\nNode2\t1\t0.5\t0.5\n"
    assert index.text("Node9", ancestors.file_range(source)) is None
    assert ancestors.read_index(ancestors.write_index(index)) == index
    assert source.read_text().count("\n") == 5  # the whole table is untouched
