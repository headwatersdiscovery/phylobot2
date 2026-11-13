"""Cross-combination disagreement on the library page (phylobot.library_disagreement):
sites where the same clade's ML residue differs across combinations."""

from phylobot.library_disagreement import Overview, State, analyse, compare, disputed_columns


def test_identical_reconstructions_have_no_disputes() -> None:
    extant = {"x": "ACDE", "y": "ACDF"}
    a = State("mafft/iqtree/WAG", "ACDE", extant)
    b = State("muscle/iqtree/LG", "ACDE", extant)
    dispute = compare([a, b])
    assert dispute is not None and dispute.combinations == 2 and dispute.disputed == 0


def test_a_differing_site_is_disputed() -> None:
    extant = {"x": "ACDE", "y": "ACDF"}
    a = State("mafft/iqtree/WAG", "ACDE", extant)
    b = State("muscle/iqtree/LG", "ACDQ", extant)  # differs at the 4th column
    dispute = compare([a, b])
    assert dispute is not None and dispute.disputed == 1
    # Exactly one of the two combinations is counted as off the (tied) majority.
    assert sum(dispute.per_key.values()) == 1


def test_three_combinations_count_the_minority() -> None:
    extant = {"x": "ACDE", "y": "ACDF", "z": "ACDG"}
    a = State("a", "ACDE", extant)
    b = State("b", "ACDE", extant)
    c = State("c", "ACDQ", extant)  # the odd one out at column 4
    dispute = compare([a, b, c])
    assert dispute is not None and dispute.disputed == 1
    assert dispute.per_key == {"a": 0, "b": 0, "c": 1}  # c differs from the A majority


def test_fewer_than_two_combinations_is_none() -> None:
    extant = {"x": "ACDE"}
    assert compare([State("a", "ACDE", extant)]) is None
    # The same combination twice (e.g. two gap versions) is not aligner disagreement.
    assert compare([State("a", "ACDE", extant), State("a", "ACDQ", extant)]) is None


def test_disputed_columns_mark_the_frame() -> None:
    # The in-viewer track: 0-based reference columns where the combinations differ.
    extant = {"x": "ACDE", "y": "ACDF"}
    a = State("a", "ACDE", extant)
    b = State("b", "ACDQ", extant)  # differs at column 3 (0-based)
    assert disputed_columns([a, b], extant) == {3}
    assert disputed_columns([a], extant) == set()  # one combination: nothing to compare
    assert disputed_columns([a, b], {}) == set()  # no reference frame


def test_analyse_groups_clades_and_sums_disjoint_sites() -> None:
    extant = {"x": "ACDE", "y": "ACDF"}
    groups = {
        "clade1": [State("a", "ACDE", extant), State("b", "ACDQ", extant)],  # 1 dispute
        "clade2": [State("a", "ACDE", extant), State("b", "AKDE", extant)],  # 1 dispute (col 2)
        "clade3": [State("a", "ACDE", extant)],  # single combination: skipped
    }
    disputes, overview = analyse(groups)
    assert set(disputes) == {"clade1", "clade2"}
    assert overview == Overview(clades=2, disputed_sites=2)
