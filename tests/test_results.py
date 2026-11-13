import pytest

from phylobot import results


def test_tint_runs_red_to_blue() -> None:
    assert results.tint(1.0) == "background-color: rgb(161, 201, 225)"  # best: light blue
    assert results.tint(0.5) == "background-color: rgb(251, 251, 251)"  # neutral
    assert results.tint(0.0) == "background-color: rgb(234, 176, 166)"  # worst: light red
    assert results.tint(7) == results.tint(1.0)  # clamped


@pytest.mark.parametrize(
    ("column", "value", "values", "expected"),
    [
        ("bic", 100.0, [100.0, 105.0, 130.0], 1.0),  # best row
        ("bic", 105.0, [100.0, 105.0, 130.0], 0.5),  # 5 units behind
        ("bic", 130.0, [100.0, 105.0, 130.0], 0.0),  # >= 10 behind: capped
        ("log_likelihood", -100.0, [-100.0, -102.5], 1.0),
        ("log_likelihood", -102.5, [-100.0, -102.5], 0.5),  # 2*dlogL = 5
        ("w_bic", 0.4, [0.4, 0.1], 1.0),  # the table's best weight
        ("w_bic", 0.1, [0.4, 0.1], 0.25),
        ("w_bic", 0.0, [0.0, 0.0], 0.0),
        ("free_parameters", 12.0, [12.0, 31.0], None),
    ],
)
def test_goodness(column: str, value: float, values: list[float], expected: float | None) -> None:
    assert results.goodness(column, value, values) == expected


def test_long_modelfinder_tables_are_folded() -> None:
    header = "\t".join(
        ["model", "log_likelihood", "aic", "aic_in_95", "w_aic", "aicc", "aicc_in_95", "w_aicc",
         "bic", "bic_in_95", "w_bic"]
    )  # fmt: skip
    rows = [
        f"M{i}\t-{100 + i}\t{200 + i}\tno\t0\t{201 + i}\tno\t0\t{210 + i}\tno\t0" for i in range(40)
    ]
    table = results.modelfinder_table("\n".join([header, *rows]))
    assert len(table.rows) == results.MODELFINDER_SHOWN and table.total == 40
    assert table.title == "ModelFinder candidates"


def test_one_modelfinder_table_per_alignment() -> None:
    header = ("model\tlog_likelihood\taic\taic_in_95\tw_aic\taicc\taicc_in_95\tw_aicc\tbic"
              "\tbic_in_95\tw_bic")  # fmt: skip
    ranking = f"{header}\nLG+G4\t-100\t200\tyes\t1\t201\tyes\t1\t210\tyes\t1\n"
    combos = [{"key": f"mafft/{tree}/modelfinder", "aligner": "mafft", "tree": tree,
               "model": "modelfinder", "path": f"mafft/{tree}/modelfinder"}
              for tree in ("iqtree", "raxml")]  # fmt: skip
    summary = {"combinations": combos}
    shared = {c["key"]: ranking for c in combos}  # one MODELFINDER run, both sources
    [result] = results.alignment_results(summary, None, shared, lambda path: path)
    assert [t.title for t in result.modelfinder] == ["ModelFinder candidates"]
    # older jobs: a provided tree ran its own ModelFinder, with its own ranking
    own = shared | {"mafft/raxml/modelfinder": ranking.replace("LG+G4", "WAG+G4")}
    [result] = results.alignment_results(summary, None, own, lambda path: path)
    assert [t.title for t in result.modelfinder] == [
        "ModelFinder candidates · tree: iqtree",
        "ModelFinder candidates · tree: raxml",
    ]


def test_alignment_stats() -> None:
    fasta = ">a\nMK-T\n>b\nM--T\n>c\nMKAT\n"  # 3 gaps in 12 cells
    stats = results.alignment_stats(fasta)
    assert (stats.sequences, stats.columns) == (3, 4)
    assert stats.gap_percent == pytest.approx(25.0)
    wrapped = ">a\nMK\n-T\n>b\nM-\n-T\n"  # line-wrapped FASTA
    assert results.alignment_stats(wrapped).columns == 4
    assert results.alignment_stats("").gap_percent == 0.0


def test_modelfinder_pick_that_was_also_ticked_is_one_row() -> None:
    fields = "aligner tree model_key model log_likelihood free_parameters aic aicc bic "
    fields += "delta_bic w_bic w_aic w_aicc rank_bic"
    rows = [
        "mafft iqtree modelfinder LG+G4 -100 10 220 221 250 0 0.9 0.9 0.9 1",
        "mafft iqtree LG+G4 LG+G4 -100 10 220 221 250 0 0.9 0.9 0.9 1",
        "mafft iqtree WAG+G4 WAG+G4 -102 10 224 225 254 4 0.1 0.1 0.1 2",
    ]
    comparison = "\n".join(line.replace(" ", "\t") for line in [fields, *rows])
    combos = [
        {"key": f"mafft/iqtree/{m}", "aligner": "mafft", "tree": "iqtree", "model": m}
        for m in ("modelfinder", "LG+G4", "WAG+G4")
    ]
    summary = {"primary": "mafft/iqtree/LG+G4", "combinations": combos}
    [result] = results.alignment_results(summary, comparison, {}, lambda path: path)
    models = [row[1] for row in result.table.rows]
    assert [cell.text for cell in models] == ["ModelFinder\u00a0→ LG+G4 (also ticked)", "WAG+G4"]
    # without statistics (older jobs) the fit in summary.json decides
    picks = ("LG+G4", "LG+G4", "WAG+G4")
    fitted = [{**c, "fit": {"model": m}} for c, m in zip(combos, picks, strict=True)]
    older = {"primary": "", "combinations": fitted}
    [result] = results.alignment_results(older, None, {}, lambda path: path)
    assert len(result.table.rows) == 2


def test_codon_and_amino_acid_models_get_separate_tables() -> None:
    fields = "aligner tree model_key model log_likelihood free_parameters aic aicc bic "
    fields += "delta_bic w_bic w_aic w_aicc rank_bic data"
    rows = [
        "mafft iqtree modelfinder LG+G4 -100 10 220 221 250 0 1 1 1 1 amino",
        "mafft iqtree GY+F3X4+G4 GY+F3X4+G4 -300 60 720 721 900 0 0.8 0.8 0.8 1 codon",
        "mafft provided1 GY+F3X4+G4 GY+F3X4+G4 -310 60 740 741 920 20 0.2 0.2 0.2 2 codon",
    ]
    comparison = "\n".join(line.replace(" ", "\t") for line in [fields, *rows])
    keys = ["mafft/iqtree/modelfinder", "mafft/iqtree/GY+F3X4+G4", "mafft/provided1/GY+F3X4+G4"]
    combos = [
        {"key": k, "aligner": "mafft", "tree": k.split("/")[1], "model": k.split("/")[2]}
        for k in keys
    ]
    summary = {"primary": keys[0], "combinations": combos}
    available = {"mafft/alignment.cds.fasta", f"{keys[1]}/ancestors.cds.fasta"}
    amino, codon = results.alignment_results(
        summary, comparison, {}, lambda path: path, available=available
    )
    assert (amino.data, codon.data) == ("amino", "codon")
    assert amino.table.title == "Alignment: mafft · amino-acid models"
    assert codon.table.title == "Alignment: mafft · codon models"
    assert (amino.first, codon.first) == (True, False)  # the overview lists mafft once
    assert codon.best == keys[1] and codon.anchor == "alignment-mafft-codon"
    assert [row[0].text for row in codon.table.rows] == ["1", "2"]  # ranked among codons
    assert "coding DNA" in codon.table.rows[0][-1].title  # in the row's ZIP
    assert codon.codon_href == "mafft/alignment.cds.fasta"


def test_codeml_section() -> None:
    codeml_tsv = "\n".join([
        "model\tlog_likelihood\tfree_parameters\tkappa\ttree_length\testimates\tstarts",
        "M7\t-4546.038\t50\t1.33\t8.6\tbeta(p = 0.278, q = 0.209)\t0.4: -4546.0380",
        "M8\t-4479.934\t52\t1.55\t9.0\t0.737 beta; 0.263 at ω = 3.08\t0.5: -4479.9",
    ])  # fmt: skip
    selection = "test\tstatistic\tdf\tp_value\nM7 vs M8\t132.2081\t2\t1.956e-29\n"
    beb = ("model\tsite\tresidue\tp_positive\tmean_omega\tse\n"
           "M8\t4\tH\t1.000\t3.392\t0.310\nM8\t6\tV\t0.668\t2.553\t1.177\n"
           "M8\t7\tE\t0.998\t3.388\t0.325\n")  # fmt: skip
    combo = {"key": "mafft/iqtree/modelfinder_codon", "aligner": "mafft", "tree": "iqtree",
             "model": "modelfinder_codon", "tree_label": "IQ-TREE",
             "model_label": "ModelFinder (codon)"}  # fmt: skip
    result = results.codeml_result(
        combo, codeml_tsv, selection, beb, lambda path: path,
        lambda key, model: f"view?{key}&{model}",
    )  # fmt: skip
    assert result.table.title == "PAML codeml on IQ-TREE · ModelFinder (codon)"
    m8 = result.table.rows[1]
    assert m8[0].text == "M8" and m8[4].text == f"{2 * 52 + 2 * 4479.934:.1f}"  # AIC
    links = dict(m8[-1].links)
    assert links["view"] == "view?mafft/iqtree/modelfinder_codon&M8"
    assert links["coding DNA"] == "mafft/iqtree/modelfinder_codon/codeml/M8/ancestors.cds.fasta"
    assert result.tests == [
        "M7 vs M8: 2ΔlnL = 132.2, p = 2e-29 (χ², 2 df): positive selection supported at 5%."
    ]
    assert result.sites[0].startswith("M8: 2 sites with P(dN/dS > 1) > 0.95")
    assert "4H, 7E" in result.sites[0] and "6V" not in result.sites[0]


def test_reconciled_combinations_show_duplications_and_losses() -> None:
    combos = [
        {"key": f"mafft/{t}/LG+G4", "aligner": "mafft", "tree": t, "model": "LG+G4"}
        for t in ("iqtree", "raxml")
    ]
    summary = {"primary": "mafft/iqtree/LG+G4", "combinations": combos}
    events = {"mafft/iqtree/LG+G4": {"duplications": 3, "losses": 12}}
    [result] = results.alignment_results(
        summary, None, {}, lambda path: path, reconciliation=events
    )
    column = result.table.columns.index("Duplications / losses")
    cells = {row[1].text: row[column] for row in result.table.rows}
    assert cells["iqtree"].text == "3 / 12"
    events["mafft/iqtree/LG+G4"]["weak"] = 2
    [weak] = results.alignment_results(summary, None, {}, lambda path: path, reconciliation=events)
    row = next(r for r in weak.table.rows if r[1].text == "iqtree")
    assert row[column].text == "3 / 12 · 2 doubtful"
    assert "doubtful duplication" in weak.table.subtitle
    assert cells["iqtree"].href == "mafft/iqtree/LG+G4/reconciliation.tsv"
    assert cells["raxml"].text == "–"
    [plain] = results.alignment_results(summary, None, {}, lambda path: path)
    assert "Duplications / losses" not in plain.table.columns


def test_reconciled_tree_rows_show_the_au_test() -> None:
    combos = [
        {"key": f"mafft/{t}/LG+G4", "aligner": "mafft", "tree": t, "model": "LG+G4"}
        for t in ("iqtree", "raxml", "generax")
    ]
    summary = {"primary": "mafft/iqtree/LG+G4", "combinations": combos}
    tsv = (
        "tree\tlogL\tdeltaL\tbp_rell\tp_kh\tp_sh\tc_elw\tp_au\n"
        "generax\t-5767.06\t1.17\t0.368\t0.329\t0.329\t0.38\t0.329\n"
        "iqtree\t-5765.90\t0\t0.632\t0.671\t1\t0.62\t1\n"
        "raxml\t-5790.00\t24.1\t0\t0.01\t0.02\t0.0\t0.0123\n"
    )
    topology = results.topology_rows("mafft", "LG+G4", tsv)
    [result] = results.alignment_results(summary, None, {}, lambda path: path, topology=topology)
    column = result.table.columns.index("AU\u00a0p")
    cells = {row[1].text: row[column] for row in result.table.rows}
    assert cells["generax"].text == "0.33"
    assert cells["generax"].href == "mafft/topology_tests/LG+G4.tsv"
    assert cells["raxml"].text == "0.012 · rejected"
    assert "AU test" in result.table.subtitle
    [plain] = results.alignment_results(summary, None, {}, lambda path: path)
    assert "AU\u00a0p" not in plain.table.columns
