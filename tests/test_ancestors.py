import statistics

import pytest

from phylobot import ancestors

POSTERIORS = (
    "node\tsite\tA\tR\tN\n"
    "Node1\t1\t0.10\t0.85\t0.05\n"
    "Node1\t2\t0.995\t0.004\t0.001\n"
    "Node1\t3\t0.40\t0.40\t0.20\n"
    "Node2\t1\t1.0\t0.0\t0.0\n"
    "Node2\t2\t1.0\t0.0\t0.0\n"
    "Node2\t3\t1.0\t0.0\t0.0\n"
)


def test_distributions_sorted_most_probable_first() -> None:
    dists = ancestors.distributions(POSTERIORS, "Node1")
    assert [d[0] for d in dists] == [("R", 0.85), ("A", 0.995), ("A", 0.40)]
    assert dists[2][:2] == [("A", 0.40), ("R", 0.40)]  # ties: alphabetical
    assert ancestors.posterior_lines(POSTERIORS, "Node2")[0].startswith("node\tsite")
    assert len(ancestors.posterior_lines(POSTERIORS, "Node2")) == 4
    assert ancestors.posterior_lines(POSTERIORS, "Node") == ["node\tsite\tA\tR\tN"]


def test_tiles_skip_tiny_states_and_number_the_extant() -> None:
    dists = ancestors.distributions(POSTERIORS, "Node1")
    rows = ancestors.site_rows("original", "RAA", dists, ["M-K"])
    assert [r.column for r in rows] == [1, 2, 3]
    assert [r.extants for r in rows] == [[("M", 1)], [("–", None)], [("K", 2)]]
    two = ancestors.site_rows("original", "RAA", dists, ["M-K", "-QR"])
    assert [r.extants for r in two][1] == [("–", None), ("Q", 1)]
    with pytest.raises(ValueError, match="columns differ"):
        ancestors.site_rows("original", "RAA", dists, ["M-K", "QR"])
    assert [t.state for t in rows[0].tiles] == ["R", "A", "N"]
    assert [t.state for t in rows[1].tiles] == ["A"]  # 0.004 and 0.001 < MIN_SHOWN
    assert all(t.style.startswith("background-color") for r in rows for t in r.tiles)


def test_masked_and_external_versions() -> None:
    dists = ancestors.distributions(POSTERIORS, "Node1")
    masked = ancestors.site_rows("fitch", "R-A", dists, ["MKK"])
    assert [t.kind for t in masked[1].tiles] == ["gap"]
    assert masked[0].tiles[0].pp == 0.85
    prank = ancestors.site_rows("prank", "k-A", dists, ["MKK"])
    assert [(t.state, t.pp, t.kind) for t in prank[0].tiles] == [("K", None, "external")]
    assert prank[1].tiles[0].kind == "gap"
    with pytest.raises(ValueError, match="columns differ"):
        ancestors.site_rows("original", "RA", dists, ["MKK"])


def test_palettes_run_warm_to_cool() -> None:
    for palette in ancestors.PALETTES:
        stops = ancestors.PALETTES[palette][1]
        assert stops[0][0] == 0.0 and stops[-1][0] == 1.0
        assert ancestors.colour(1.0, palette).startswith(
            "background-color: rgb({}, {}, {})".format(*stops[-1][1])
        )
        assert ancestors.colour(0.0, palette).startswith(
            "background-color: rgb({}, {}, {})".format(*stops[0][1])
        )
    assert "rgb(69, 117, 180); color: #fff" in ancestors.colour(1.0, "coolwarm")
    assert ancestors.colour(0.5, "nope") == ancestors.colour(0.5, ancestors.DEFAULT_PALETTE)


def test_node_choices_root_first_with_descendants() -> None:
    newick = "(out,((a,b)Node3,(c,d,e,f)Node4)Node2)Node1;"
    choices = ancestors.node_choices(newick, ["Node4", "Node3", "Node2", "Node1"], {"Node3": "97"})
    assert [c.name for c in choices] == ["Node1", "Node2", "Node3", "Node4"]
    assert choices[2].label == "Node3 · 2 sequences · support 97"  # no names: trees grow
    assert choices[3].label == "Node4 · 4 sequences"
    assert choices[3].clade == {"c", "d", "e", "f"}


def test_parse_fasta() -> None:
    assert ancestors.parse_fasta(">x desc\nAB\nC\n>y\n\nD\n") == {"x": "ABC", "y": "D"}


def test_ml_gap_probabilities_become_tiles() -> None:
    dists = ancestors.distributions(POSTERIORS, "Node1")
    gaps_tsv = "node\tsite\tp_gap\nNode1\t3\t0.9\nNode1\t1\t0.001\nNode1\t2\t0.2\n"
    p_gap = ancestors.gap_probabilities(gaps_tsv, "Node1")
    assert p_gap == [0.001, 0.2, 0.9]  # in site order
    rows = ancestors.site_rows("ml", "RA-", dists, ["MKK"], p_gap=p_gap)
    assert [t.kind for t in rows[0].tiles] == ["state"] * 3  # P(gap) 0.001: no tile
    assert [(t.state, t.pp) for t in rows[1].tiles] == [("A", 0.995), ("–", 0.2)]
    assert rows[1].tiles[1].style.startswith("background-color")
    assert [(t.state, t.pp, t.kind) for t in rows[2].tiles] == [("–", 0.9, "gap")]
    with pytest.raises(ValueError, match="gap probabilities"):
        ancestors.site_rows("ml", "RA-", dists, ["MKK"], p_gap=[0.1])


def test_node_plots_summarise_the_ml_residues() -> None:
    dists = ancestors.distributions(POSTERIORS, "Node1")  # ML pp 0.85, 0.995, 0.40
    plots = ancestors.node_plots("original", "RAA", dists, "coolwarm")
    assert plots is not None and plots.sites == 3
    assert round(plots.mean, 4) == round((0.85 + 0.995 + 0.40) / 3, 4)
    assert plots.sd == pytest.approx(statistics.pstdev([0.85, 0.995, 0.40]))
    shares = {bar.label: bar.proportion for bar in plots.histogram if bar.proportion}
    assert shares == {"0.4–0.5": 1 / 3, "0.8–0.9": 1 / 3, "0.9–1.0": 1 / 3}
    assert plots.histogram[9].fill.startswith("rgb(")
    assert [t.rank for t in plots.traces] == ["ML residue", "2nd best", "3rd best"]
    assert all(len(t.points.split()) == 3 for t in plots.traces)
    assert plots.traces[0].colour == ancestors.LINE_STYLES["standard"][1][0][0]
    dashed = ancestors.node_plots("original", "RAA", dists, "coolwarm", "dashes")
    assert dashed is not None and [t.dash for t in dashed.traces] == ["", "6 4", "1.5 3"]
    fallback = ancestors.node_plots("original", "RAA", dists, "coolwarm", "nope")
    assert fallback is not None and fallback.traces == plots.traces
    # tooltips: site in the ancestor, alignment column, then the top three states
    assert plots.tooltip[1] == [2, 2, "A", 0.995, "R", 0.004, "N", 0.001]
    # a masked version plots only the ancestor's own sites, numbered along it
    masked = ancestors.node_plots("fitch", "R-A", dists, "coolwarm")
    assert masked is not None and masked.sites == 2
    assert [t[:2] for t in masked.tooltip] == [[1, 1], [2, 3]]
    assert ancestors.node_plots("prank", "RAA", dists, "coolwarm") is None
    assert ancestors.node_plots("fitch", "---", dists, "coolwarm") is None


def test_nice_step() -> None:
    assert [ancestors.nice_step(n) for n in (3, 10, 11, 45, 461, 1052, 20000)] == [
        1, 1, 2, 5, 50, 200, 2000,
    ]  # fmt: skip


def test_codon_states() -> None:
    tsv = ("node\tsite\tcodon\tp\nN1\t1\tATG\t1.00000\nN1\t2\tCTT\t0.40000\n"
           "N1\t2\tATT\t0.35000\nN1\t2\tATC\t0.25000\nN2\t1\tATG\t1.00000\n")  # fmt: skip
    dists = ancestors.codon_distributions(tsv, "N1")
    assert dists == [[("ATG", 1.0)], [("CTT", 0.4), ("ATT", 0.35), ("ATC", 0.25)]]
    notes = {"ATG": "M", "CTT": "L", "ATT": "I", "ATC": "I"}
    rows = ancestors.site_rows("original", "MI", dists, [], notes=notes)
    assert [(t.state, t.note) for t in rows[1].tiles] == [("CTT", "L"), ("ATT", "I"), ("ATC", "I")]
    assert ancestors.coding_file("original") == "ancestors.cds.fasta"
    assert ancestors.coding_file("fitch") == "ancestors.gaps-fitch.cds.fasta"
    assert ancestors.coding_file("codeml-M8") == "codeml/M8/ancestors.cds.fasta"


def test_short_version_labels_name_every_version() -> None:
    """The library table shows these; the menus' explanatory labels wrap a row to three
    lines. Every version needs one, including a codeml reconstruction."""
    for version in ancestors.VERSIONS:
        short = ancestors.short_version_label(version)
        assert short and len(short) <= len(ancestors.version_label(version))
    assert ancestors.short_version_label("original") == "None"
    assert ancestors.short_version_label("codeml-M3") == "codeml M3"
