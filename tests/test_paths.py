import pytest

from phylobot import paths

TREE = "(out,((a,b)Node3,(c,d)Node4)Node2)Node1;"


def test_descendants_and_path() -> None:
    below = paths.descendants(TREE, "Node2")
    assert [(n.name, n.extant, n.leaves) for n in below] == [
        ("Node3", False, 2), ("a", True, 1), ("b", True, 1),
        ("Node4", False, 2), ("c", True, 1), ("d", True, 1),
    ]  # fmt: skip
    assert paths.descendants(TREE, "a") == [] and paths.descendants(TREE, "Nope") == []
    assert paths.path(TREE, "Node1", "c") == ["Node1", "Node2", "Node4", "c"]
    assert paths.path(TREE, "Node2", "Node3") == ["Node2", "Node3"]
    assert paths.path(TREE, "Node3", "c") == []  # not downstream
    assert paths.path(TREE, "Node4", "Node2") == []  # upstream
    assert paths.path(TREE, "a", "a") == []  # starts must be ancestors


def test_ml_probabilities() -> None:
    tsv = "node\tsite\tstate\tpp\nN1\t2\tK\t0.5\nN1\t1\tM\t0.9\nN2\t1\tM\t1.0\n"
    assert paths.ml_probabilities(tsv, {"N1"}) == {"N1": [0.9, 0.5]}


def test_trace_marks_changes_branch_by_branch() -> None:
    sequences = {"Node1": "MKLA-", "Node2": "MRLA-", "c": "MRL-Q", "x": "AKLAQ"}
    probs = {"Node1": [1.0, 0.6, 1.0, 0.9, 1.0], "Node2": [1.0, 0.8, 1.0, 0.95, 1.0]}
    table = paths.trace(["Node1", "Node2", "c"], sequences, probs, {"c", "x"}, ["x"], "coolwarm")
    assert [b.child for b in table.branches] == ["Node2", "c"]
    first, second = table.branches
    assert [c.label for c in first.changes] == ["K2R"] and first.substitutions == 1
    assert first.changes[0].pp_before == 0.6 and first.changes[0].pp_after == 0.8
    assert [c.label for c in second.changes] == ["A4-", "-5Q"]
    assert (second.losses, second.gains, second.substitutions) == (1, 1, 0)
    # only the changing columns (2, 4, 5), numbered along the first sequence
    assert table.columns == [2, 4, 5] and table.first_numbers == [2, 4, None]
    assert table.changed_columns == 3
    rows = {row.name: row for row in table.rows}
    assert [c.changed for c in rows["Node1"].cells] == [False, False, False]
    assert [c.changed for c in rows["Node2"].cells] == [True, False, False]
    assert [c.char for c in rows["c"].cells] == ["R", "–", "Q"]
    assert rows["Node2"].cells[0].style.startswith("background-color")
    assert rows["c"].cells[0].style == ""  # extant: no posterior probability
    assert not rows["x"].on_path and not any(c.changed for c in rows["x"].cells)
    full = paths.trace(["Node1", "Node2"], sequences, probs, set(), [], "coolwarm", False)
    assert full.columns == [1, 2, 3, 4, 5]
    partial = paths.trace(["Node1", "Node2"], sequences, {"Node1": [1.0]}, set(), [], "x")
    assert all(cell.style == "" for row in partial.rows for cell in row.cells)
    with pytest.raises(ValueError, match="columns"):
        paths.trace(["Node1", "Node2"], sequences | {"Node2": "MR"}, probs, set(), [], "x")


def test_codon_paths_class_changes_synonymous_or_not() -> None:
    # Column 1: Node1 CTG (L) -> Node2 CTA (L): silent; column 2: AAA (K) -> AGA (R).
    sequences = {"Node1": "LKM", "Node2": "LRM", "c": "LRM"}
    codons = {"Node1": "CTGAAAATG", "Node2": "CTAAGAATG", "c": "CTAAGGATG"}
    probs = {"Node1": [1.0, 0.6, 1.0], "Node2": [1.0, 0.8, 1.0]}
    codon_pp = {"Node1": [0.7, 0.5, 1.0], "Node2": [0.6, 0.4, 1.0]}
    route = ["Node1", "Node2", "c"]
    table = paths.trace(route, sequences, probs, {"c"}, [], "coolwarm", codons=codons)
    first, second = table.branches
    assert [(c.label, c.kind, c.codons) for c in first.changes] == [
        ("L1L", "synonymous", "CTG → CTA"), ("K2R", "substitution", "AAA → AGA"),
    ]  # fmt: skip
    assert (first.substitutions, first.synonymous) == (1, 1)
    assert [c.label for c in second.changes] == ["R2R"]  # AGA -> AGG, silent
    assert table.columns == [1, 2] and table.codons
    node2 = next(r for r in table.rows if r.name == "Node2").cells
    assert (node2[0].silent, node2[0].changed) == (True, False)
    assert (node2[1].silent, node2[1].changed) == (False, True)
    assert node2[0].codon == "CTA" and node2[0].char == "L"
    # Hidden: silent changes leave the list and the columns, but are still counted.
    hidden = paths.trace(route, sequences, probs, {"c"}, [], "x", codons=codons,
                         show_silent=False)  # fmt: skip
    assert [c.label for c in hidden.branches[0].changes] == ["K2R"]
    assert hidden.branches[0].synonymous == 1 and hidden.branches[1].hidden_silent == 1
    assert hidden.columns == [2]
    assert not any(cell.silent for row in hidden.rows for cell in row.cells)
    # Codon cells, coloured by the ML codon's probability.
    shown = paths.trace(route, sequences, probs, {"c"}, [], "coolwarm", codons=codons,
                        codon_probabilities=codon_pp, cells="codons")  # fmt: skip
    node2 = next(r for r in shown.rows if r.name == "Node2").cells
    assert [(c.char, c.pp) for c in node2] == [("CTA", 0.6), ("AGA", 0.4)]
    assert [c.char for c in shown.rows[-1].cells] == ["CTA", "AGG"]  # extant
    # Without codons nothing changes: no codon cells, no silent changes.
    plain = paths.trace(route, sequences, probs, {"c"}, [], "x", cells="codons")
    assert plain.cells == "amino" and not plain.codons and plain.columns == [2]
    with pytest.raises(ValueError, match="nucleotides"):
        paths.trace(route, sequences, probs, {"c"}, [], "x", codons=codons | {"c": "CTA"})


def test_branch_rates_keyed_by_both_ends() -> None:
    tsv = ("parent\tchild\tt\tN\tS\tdN_dS\tdN\tdS\tN_dN\tS_dS\n"
           "Node1\tNode2\t0.1\t300\t100\t0.5\t0.01\t0.02\t3\t2\n"
           "Node2\tc\t0.1\t300\t100\t999\t0.02\t0\t6\t0\n")  # fmt: skip
    rates = paths.branch_rates(tsv)
    assert rates[frozenset(("Node2", "Node1"))].omega == 0.5
    assert rates[frozenset(("Node1", "Node2"))].defined
    assert not rates[frozenset(("c", "Node2"))].defined  # dS = 0: codeml's 999
    sequences = {"Node1": "MK", "Node2": "MR", "c": "MR"}
    table = paths.trace(["Node1", "Node2", "c"], sequences, {}, {"c"}, [], "x", rates=rates)
    assert table.branches[0].rate and table.branches[0].rate.omega == 0.5
    assert table.branches[1].rate and not table.branches[1].rate.defined
    assert paths.trace(["Node1", "Node2"], sequences, {}, set(), [], "x").branches[0].rate is None


def test_ambiguous_extant_codons_change_only_where_they_must() -> None:
    # Column 1: Node1 AAC (N) -> extant AMC (AAC or ACC, X): may be the same codon, so no
    # change. Column 2: AAT (N) -> AMC: not AAT whatever it is, but N or T: uncertain.
    # Column 3: CTG (L) -> CTY (CTC or CTT, L): another codon for sure, the same amino
    # acid: silent.
    sequences = {"Node1": "NNL", "c": "XXL"}
    codons = {"Node1": "AACAATCTG", "c": "AMCAMCCTY"}
    probs = {"Node1": [1.0, 1.0, 1.0]}
    table = paths.trace(["Node1", "c"], sequences, probs, {"c"}, [], "x", codons=codons)
    (branch,) = table.branches
    assert [(c.label, c.kind) for c in branch.changes] == [
        ("N2X", "uncertain"), ("L3L", "synonymous"),
    ]  # fmt: skip
    assert (branch.substitutions, branch.uncertain, branch.synonymous) == (0, 1, 1)
    assert table.uncertain and table.columns == [2, 3]
    full = paths.trace(["Node1", "c"], sequences, probs, {"c"}, [], "x", codons=codons,
                       only_changes=False)  # fmt: skip
    assert [cell.changed for cell in full.rows[-1].cells] == [False, True, False]
