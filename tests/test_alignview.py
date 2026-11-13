import math
from collections import Counter

import pytest

from phylobot import alignview
from phylobot.alignview import AMINO, BLOSUM62


def test_blosum62_is_symmetric_with_known_values() -> None:
    assert all(BLOSUM62[a, b] == BLOSUM62[b, a] for a in AMINO for b in AMINO)
    assert (BLOSUM62["W", "W"], BLOSUM62["C", "C"], BLOSUM62["A", "A"]) == (11, 9, 4)
    assert (BLOSUM62["I", "V"], BLOSUM62["W", "D"]) == (3, -4)
    # Capra & Singh's rounded values, used as published (they sum to 1.002)
    assert sum(alignview.BACKGROUND.values()) == pytest.approx(1.0, abs=0.005)


def test_entropy_and_sum_of_pairs() -> None:
    assert alignview.entropy(Counter("AAAA"), 20) == 0
    assert alignview.entropy(Counter(AMINO), 20) == pytest.approx(1.0)
    assert alignview.sum_of_pairs("AA") == 4
    assert alignview.sum_of_pairs("AR-") == -1  # gaps left out
    assert alignview.sum_of_pairs("AAR") == pytest.approx((4 - 1 - 1) / 3)
    assert alignview.sum_of_pairs("A--") is None


def test_jensen_shannon_conservation() -> None:
    weights = [0.25] * 4
    conserved = alignview.jensen_shannon("WWWW", weights)
    mixed = alignview.jensen_shannon("ALVK", weights)
    assert conserved is not None and mixed is not None and 0 <= mixed < conserved <= 1
    # a rare residue conserved scores higher than a common one conserved
    common = alignview.jensen_shannon("LLLL", weights)
    assert common is not None and common < conserved
    assert alignview.jensen_shannon("WW--", weights) is None  # over 30% gaps
    # one gap of five: scaled by the weighted non-gap fraction
    full = alignview.jensen_shannon("WWWWW", [0.2] * 5)
    one_gap = alignview.jensen_shannon("WWWW-", [0.2] * 5)
    assert full and one_gap and one_gap < full
    assert alignview.windowed([1.0, None, 0.0], 1) == [1.0, None, 0.0]  # no scored neighbours
    assert alignview.windowed([1.0, 0.0, 0.0], 1) == [0.5, 0.25, 0.0]


def test_henikoff_weights() -> None:
    rows = ["AAA", "AAA", "CCC"]
    weights = alignview.henikoff_weights(alignview.columns(rows), 3)
    assert sum(weights) == pytest.approx(1.0)
    assert weights[2] == pytest.approx(0.5) and weights[0] == pytest.approx(0.25)


def test_composition_tracks() -> None:
    rows = ["MK-D", "MR-E", "MK-W"]
    tracks = alignview.composition_tracks(rows)
    assert set(tracks) == set(alignview.COMPOSITION)
    assert tracks["consensus"].text == ["M", "K", "-", "D"]
    values = {k: t.series[0].values for k, t in tracks.items()}
    assert values["consensus"][:3] == [1.0, pytest.approx(2 / 3), None]
    assert values["occupancy"] == [1.0, 1.0, 0.0, 1.0]
    assert values["distinct"] == [1.0, 2.0, 0.0, 3.0]
    assert values["entropy"][0] == 0 and values["entropy"][2] is None
    assert values["entropy_gaps"][2] == 0  # all gaps: one state
    assert values["charge"][1] == 1.0 and values["charge"][3] == pytest.approx(-2 / 3)
    assert values["hydropathy"][0] == pytest.approx(1.9)
    assert values["jsd"][2] is None and values["jsd"][0] is not None
    profile = values["hydropathy_profile"]  # of the first row: M K D (gap skipped)
    assert profile[2] is None and profile[0] == pytest.approx((1.9 - 3.9 - 3.5) / 3)


def test_tree_order_ladderized_with_outgroup_last() -> None:
    newick = "((out,x)N2,((a,b)N4,c)N3)N1;"
    assert alignview.tree_order(newick) == ["out", "x", "c", "a", "b"]
    assert alignview.tree_order(newick, ["out"]) == ["c", "a", "b", "out", "x"]
    assert alignview.tree_order(newick, ["out"], internal=True) == [
        "N1", "N3", "c", "N4", "a", "b", "N2", "out", "x"
    ]  # fmt: skip


def test_substitutions_mean_pp_and_gap_shares() -> None:
    sequences = {"N1": "MKA", "N2": "MRA", "a": "MR-", "b": "LKA", "c": "MKA"}
    counts = alignview.substitutions("((a,b)N2,c)N1;", sequences, 3)
    assert counts == [1.0, 2.0, 0.0]  # b: M->L; N1->N2 K->R, N2->b R->K; a's gap skipped
    pp = "node\tsite\tstate\tpp\nN1\t1\tM\t0.5\nN2\t1\tM\t1.0\nN1\t2\tK\t0.8\n"
    assert alignview.mean_pp(pp, 3) == [0.75, 0.8, None]
    shares = alignview.gap_shares({"Fitch": {"N1": "M-A", "N2": "M--"}}, 3)
    assert shares[0].label == "Fitch" and shares[0].values == [0.0, 1.0, 0.5]
    tsv = "site\trate\n1\t0.5\n3\tx\n9\t2\n"
    assert alignview.tsv_column(tsv, "rate", 3) == [0.5, None, None]


def test_clustal_and_tracks_tsv() -> None:
    text = alignview.clustal({"a": "MKV", "b": "MRI"}, block=2)
    lines = text.splitlines()
    assert lines[0].startswith("CLUSTAL W")
    assert lines[3:6] == ["a      MK", "b      MR", "       *:"]
    assert lines[7:10] == ["a      V", "b      I", "       :"]
    assert "*" not in alignview.clustal({"a": "MK", "b": "MK"}, protein=False)
    track = alignview.Track("x", "X", "", [alignview.Series("X", [1.0, None])])
    assert alignview.tracks_tsv([track], 2) == "column\tX\n1\t1\n2\t\n"
    assert math.isclose(track.json()["series"][0]["values"][0], 1.0)
