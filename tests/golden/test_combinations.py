"""Every alignment x tree combination must receive its own inputs.

Runs the real pipeline with MAFFT, two provided alignments, IQ-TREE and RAxML-NG searches
and a provided tree (9 combinations), where each alignment has a distinct width and only
IQ-TREE's searched trees carry UFBoot support. Nextflow finishes tasks in arbitrary order,
so a pairing mistake in main.nf shows up here as a mismatched alignment or tree.
"""

import json
from pathlib import Path

import pytest

from phylobot.seqio import read_fasta, write_fasta

from .pipeline import run_pipeline

GOLDEN = Path(__file__).resolve().parent
ROOT = GOLDEN.parent.parent
TRUTH = GOLDEN / "ime2" / "easy" / "replicates" / "101"


@pytest.mark.golden
def test_each_combination_gets_its_own_inputs() -> None:
    out = ROOT / "results" / "golden" / "combinations"
    out.mkdir(parents=True, exist_ok=True)
    true_alignment = read_fasta(TRUTH / "true_alignment.fasta")
    widened = {name: row + "-" for name, row in true_alignment.items()}  # one extra column
    write_fasta(out / "widened.fasta", widened)
    params = {
        "alignment": [
            str((TRUTH / "true_alignment.fasta").relative_to(ROOT)),
            str((out / "widened.fasta").relative_to(ROOT)),
        ],
        "alignment_labels": ["True alignment", "True alignment + gap column"],
        "tree": [str((TRUTH / "true_tree.nwk").relative_to(ROOT))],
        "tree_labels": ["True tree"],
        "tree_builders": "iqtree,raxml",
    }
    run = run_pipeline(out, TRUTH / "input.fasta", params)
    mafft = read_fasta(run / "mafft" / "iqtree" / "modelfinder" / "alignment.fasta")
    expected = {"mafft": mafft, "provided1": true_alignment, "provided2": widened}

    summary = json.loads((run / "summary.json").read_text())
    keys = sorted(c["key"] for c in summary["combinations"])
    assert keys == [
        f"{a}/{t}/modelfinder" for a in sorted(expected) for t in ("iqtree", "provided1", "raxml")
    ]
    labels = {c["key"]: c["label"] for c in summary["combinations"]}
    assert (
        labels["provided2/provided1/modelfinder"]
        == "True alignment + gap column × True tree × ModelFinder"
    )

    for combination in summary["combinations"]:
        path = run / combination["path"]
        alignment = read_fasta(path / "alignment.fasta")
        assert alignment == expected[combination["aligner"]], combination["key"]
        width = len(next(iter(alignment.values())))
        pp_rows = (path / "ancestors_pp.tsv").read_text().splitlines()[1:]
        assert max(int(row.split("\t")[1]) for row in pp_rows) == width, combination["key"]
        supports = [r.split("\t")[1] for r in (path / "supports.tsv").read_text().splitlines()[1:]]
        searched = combination["tree"] == "iqtree"  # UFBoot: IQ-TREE's own search only
        assert any(supports) == searched, combination["key"]
    # One ModelFinder run per alignment: every tree source used the same pick
    for aligner in expected:
        picks = {(run / aligner / t / "modelfinder" / "model.txt").read_text()
                 for t in ("iqtree", "provided1", "raxml")}  # fmt: skip
        assert len(picks) == 1, aligner


@pytest.mark.golden
def test_each_model_is_its_own_combination() -> None:
    """ModelFinder plus two fixed models: three combinations, each run under its model."""
    params = {"models": ["LG+G4", "WAG+R4"], "fixed_models": True, "modelfinder": True}
    run = run_pipeline(ROOT / "results" / "golden" / "models", TRUTH / "input.fasta", params)
    summary = json.loads((run / "summary.json").read_text())
    assert summary["primary"] == "mafft/iqtree/modelfinder"
    keys = sorted(c["key"] for c in summary["combinations"])
    assert keys == ["mafft/iqtree/LG+G4", "mafft/iqtree/WAG+R4", "mafft/iqtree/modelfinder"]
    for model in ("LG+G4", "WAG+R4"):
        assert (run / "mafft" / "iqtree" / model / "model.txt").read_text().strip() == model
    # ModelFinder chose among the ticked cells' cross product: LG,WAG x G,R.
    chosen = (run / "mafft" / "iqtree" / "modelfinder" / "model.txt").read_text().strip()
    assert chosen.split("+")[0] in {"LG", "WAG"}
    assert "+G" in chosen or "+R" in chosen
    validation = json.loads((run / "validation.json").read_text())
    assert sorted(validation["modelfinder_candidates"]) == [
        "LG+G4", "LG+R4", "WAG+G4", "WAG+R4"
    ]  # fmt: skip


ALIGNERS = ("mafft", "muscle", "clustalo", "prank")


@pytest.mark.golden
def test_every_aligner_gets_its_own_alignment() -> None:
    """All four aligners on real, gapped data (the 2014 Ime2 sequences): four alignments,
    each flowing through its own tree and reconstruction, compared separately."""
    run = run_pipeline(
        ROOT / "results" / "golden" / "aligners",
        ROOT / "src" / "phylobot" / "examples" / "ime2.fasta",
        {"aligners": ",".join(ALIGNERS)},
    )
    summary = json.loads((run / "summary.json").read_text())
    keys = sorted(c["key"] for c in summary["combinations"])
    assert keys == sorted(f"{a}/iqtree/modelfinder" for a in ALIGNERS)
    labels = {c["aligner"]: c["aligner_label"] for c in summary["combinations"]}
    assert labels == {"mafft": "MAFFT", "muscle": "MUSCLE", "clustalo": "Clustal Omega",
                      "prank": "PRANK"}  # fmt: skip
    alignments = {
        a: read_fasta(run / a / "iqtree" / "modelfinder" / "alignment.fasta") for a in ALIGNERS
    }
    sequences = read_fasta(ROOT / "src" / "phylobot" / "examples" / "ime2.fasta")
    for aligner, alignment in alignments.items():
        # same sequences, names restored, residues untouched
        assert {n: r.replace("-", "") for n, r in alignment.items()} == sequences, aligner
        width = len(next(iter(alignment.values())))
        rows = (run / aligner / "iqtree" / "modelfinder" / "ancestors_pp.tsv").read_text()
        assert max(int(r.split("\t")[1]) for r in rows.splitlines()[1:]) == width, aligner
    widths = {a: len(next(iter(aln.values()))) for a, aln in alignments.items()}
    assert len(set(widths.values())) == len(ALIGNERS), widths  # real data: all differ
    comparison = (run / "model_comparison.tsv").read_text().splitlines()
    assert {line.split("\t")[0] for line in comparison[1:]} == set(ALIGNERS)
    versions = dict(line.split("\t") for line in (run / "versions.tsv").read_text().splitlines())
    assert versions["muscle"] == "5.3" and versions["clustalo"] == "1.2.4"
    assert versions["prank"] == "250331"  # as the binary reports itself


@pytest.mark.golden
def test_stop_after_alignments() -> None:
    """An alignment-only job: each alignment (an aligner's and the user's own) with its
    original names, and a summary without combinations. No trees, no ancestors."""
    params = {
        "stop_after": "alignments",
        "aligners": "mafft,muscle",
        "alignment": [str((TRUTH / "true_alignment.fasta").relative_to(ROOT))],
        "alignment_labels": ["True alignment"],
    }
    run = run_pipeline(ROOT / "results" / "golden" / "stop-alignments", TRUTH / "input.fasta",
                       params)  # fmt: skip
    summary = json.loads((run / "summary.json").read_text())
    assert summary["stop_after"] == "alignments" and summary["combinations"] == []
    assert [a["key"] for a in summary["alignments"]] == ["mafft", "muscle", "provided1"]
    assert summary["alignments"][2]["label"] == "True alignment"
    sequences = read_fasta(TRUTH / "input.fasta")
    for key in ("mafft", "muscle"):
        alignment = read_fasta(run / key / "alignment.fasta")
        assert {n: r.replace("-", "") for n, r in alignment.items()} == sequences, key
    assert read_fasta(run / "provided1" / "alignment.fasta") == read_fasta(
        TRUTH / "true_alignment.fasta"
    )
    assert not list(run.glob("*/*/*/tree.nwk")) and not list(run.glob("**/ancestors.fasta"))


@pytest.mark.golden
def test_stop_after_trees() -> None:
    """A trees-only job with the user's own tree: IQ-TREE's tree and the fitted user tree,
    each with named nodes, supports (searched tree only) and fit statistics, rooted on the
    outgroup; no ancestral reconstruction."""
    outgroup = next(iter(read_fasta(TRUTH / "input.fasta")))
    params = {
        "stop_after": "trees",
        "tree": [str((TRUTH / "true_tree.nwk").relative_to(ROOT))],
        "tree_labels": ["True tree"],
        "outgroup": [outgroup],
        "gap_methods": ["ml"],  # ignored: no ancestors
    }
    run = run_pipeline(ROOT / "results" / "golden" / "stop-trees", TRUTH / "input.fasta", params)
    summary = json.loads((run / "summary.json").read_text())
    assert summary["stop_after"] == "trees"
    keys = sorted(c["key"] for c in summary["combinations"])
    assert keys == ["mafft/iqtree/modelfinder", "mafft/provided1/modelfinder"]
    for combination in summary["combinations"]:
        path = run / combination["path"]
        assert "Node1" in (path / "tree.nwk").read_text()
        assert (path / "tree.rooted.nwk").is_file() and (path / "fit.tsv").is_file()
        supports = [r.split("\t")[1] for r in (path / "supports.tsv").read_text().splitlines()[1:]]
        assert any(supports) == (combination["tree"] == "iqtree"), combination["key"]
        assert not (path / "ancestors.fasta").exists() and not (path / "gaps-ml.tsv").exists()
    assert (run / "mafft" / "alignment.fasta").is_file()


@pytest.mark.golden
def test_codon_models_fitted_on_the_amino_acid_trees() -> None:
    """C3d: codon ModelFinder's pick fitted on IQ-TREE's and RAxML-NG's amino-acid trees,
    with no codon search. One amino-acid model (ModelFinder, the New Job page's default)
    is the case groupTuple(size: 1) silently dropped. The fits keep the amino-acid
    topology but none of its UFBoot values (FIT_IQTREE removes internal labels)."""
    lysin = ROOT / "src" / "phylobot" / "examples" / "lysin.fasta"
    params = {
        "seqtype": "codon",
        "stop_after": "trees",
        "tree_builders": ["iqtree", "raxml"],
        "codon_modelfinder": True,
        "codon_candidates": ["GY+F3X4", "MG+F3X4"],
        "codon_trees": ["amino"],
        "bootstrap": 1000,
        "alrt": 0,
    }
    run = run_pipeline(ROOT / "results" / "golden" / "codon-amino-trees", lysin, params)
    summary = json.loads((run / "summary.json").read_text())
    combos = {c["key"]: c for c in summary["combinations"]}
    assert sorted(combos) == [
        "mafft/iqtree/modelfinder",
        "mafft/iqtree_aa/modelfinder_codon",
        "mafft/raxml/modelfinder",
        "mafft/raxml_aa/modelfinder_codon",
    ]
    for source, label in (("iqtree", "IQ-TREE"), ("raxml", "RAxML-NG")):
        codon = combos[f"mafft/{source}_aa/modelfinder_codon"]
        assert codon["data"] == "codon"
        assert codon["tree_label"] == f"{label} amino-acid tree (best: ModelFinder)"
        supports = (run / codon["path"] / "supports.tsv").read_text().splitlines()[1:]
        assert supports and not any(row.split("\t")[1] for row in supports)  # no UFBoot
    # The amino-acid IQ-TREE tree had UFBoot, so stripping its labels was needed.
    searched = (run / "mafft/iqtree/modelfinder/supports.tsv").read_text().splitlines()[1:]
    assert any(row.split("\t")[1] for row in searched)
