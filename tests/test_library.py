"""phylobot.library: AltAll and the thresholds behind it (L1), and the library proper —
generator, reducer, selector and coverage readout (L3)."""

import dataclasses
from math import inf

import pytest

from phylobot import library


def dists(*columns: dict[str, float]) -> list[list[tuple[str, float]]]:
    """Per-column distributions as `ancestors.distributions` gives them: most probable
    first, ties by state."""
    return [sorted(c.items(), key=lambda sp: (-sp[1], sp[0])) for c in columns]


def test_parse_threshold_falls_back_clamps_and_rounds() -> None:
    assert library.parse_threshold("0.3") == 0.3
    assert library.parse_threshold(0.25) == 0.25
    assert library.parse_threshold("0.333") == 0.33
    # Anything unparsable or out of range becomes a usable number, never an error.
    junk: list[object] = ["", "abc", None, "nan", "inf", "-inf", [], "0.2; DROP TABLE"]
    for value in junk:
        assert library.parse_threshold(value) == library.ALT_ALL_THRESHOLD
    assert library.parse_threshold("0.9") == library.MAX_THRESHOLD
    assert library.parse_threshold("-1") == library.MIN_THRESHOLD
    assert library.parse_threshold("", default=0.4) == 0.4


def test_alt_all_switches_sites_above_the_threshold() -> None:
    design = library.alt_all(
        "MKK",
        dists(
            {"M": 0.9, "K": 0.1},  # alternative too improbable
            {"K": 0.6, "R": 0.4},  # switched
            {"K": 0.5, "R": 0.3, "Q": 0.2},  # the most probable alternative wins
        ),
        threshold=0.2,
    )
    assert design.sequence == "MRR"
    assert design.ambiguous == 2
    assert [(s.site, s.ml, s.alt) for s in design.swaps] == [(2, "K", "R"), (3, "K", "R")]
    assert design.swaps[0].pp == 0.4 and design.swaps[0].ml_pp == 0.6


def test_alt_all_threshold_is_inclusive_and_changes_the_result() -> None:
    columns = dists({"M": 0.8, "K": 0.2})
    assert library.alt_all("M", columns, threshold=0.2).sequence == "K"
    assert library.alt_all("M", columns, threshold=0.21).sequence == "M"
    assert library.alt_all("M", columns, threshold=0.3).swaps == []


def test_alt_all_keeps_masked_columns_as_gaps() -> None:
    """A gap version's masking is a separate axis: AltAll never fills one in."""
    design = library.alt_all("M-K", dists({"M": 0.5, "R": 0.5}, {"M": 0.9, "R": 0.1},
                                          {"K": 0.5, "R": 0.5}), threshold=0.2)  # fmt: skip
    assert design.sequence == "R-R"
    assert design.unaligned == "RR"
    assert [s.site for s in design.swaps] == [1, 3]


def test_alt_all_breaks_ties_by_state_whatever_the_input_order() -> None:
    forwards = library.alt_all("M", [[("R", 0.3), ("Q", 0.3)]], threshold=0.2)
    backwards = library.alt_all("M", [[("Q", 0.3), ("R", 0.3)]], threshold=0.2)
    assert forwards.sequence == backwards.sequence == "Q"


def test_alt_all_without_alternatives_returns_the_ml_sequence() -> None:
    design = library.alt_all("MK", dists({"M": 1.0}, {"K": 1.0}), threshold=0.2)
    assert design.sequence == "MK" and design.swaps == []


def test_alt_all_refuses_a_mismatched_alignment() -> None:
    with pytest.raises(ValueError, match="columns differ"):
        library.alt_all("MKK", dists({"M": 1.0}))


def test_sites_tsv_lists_the_ambiguous_sites() -> None:
    design = library.alt_all("MK", dists({"M": 0.9, "K": 0.1}, {"K": 0.6, "R": 0.4}), 0.2)
    assert design.sites_tsv().splitlines() == [
        "site\tml\tml_pp\talt\talt_pp",
        "2\tK\t0.60000\tR\t0.40000",
    ]


def test_alternatives_are_threshold_free_and_counted_at_any_threshold() -> None:
    found = library.alternatives(
        "MKK", dists({"M": 0.9, "K": 0.1}, {"K": 0.6, "R": 0.4}, {"K": 0.7, "Q": 0.3})
    )
    assert [(s.site, s.alt, s.pp) for s in found] == [(1, "K", 0.1), (2, "R", 0.4), (3, "Q", 0.3)]
    assert library.count_at(found, 0.1) == 3
    assert library.count_at(found, 0.3) == 2
    assert library.count_at(found, 0.45) == 0
    # The same alternatives drive the sequence, so count and file always agree.
    assert library.alt_all("MKK", [], 0.3, found=found).sequence == "MRQ"


def test_alt_all_on_codons() -> None:
    """Codon reconstructions switch whole codons: three characters a column."""
    columns = dists({"ATG": 0.9, "TTG": 0.1}, {"AAA": 0.6, "AGA": 0.4})
    design = library.alt_all("ATGAAA", columns, threshold=0.2, width=3)
    assert design.sequence == "ATGAGA"
    assert [(s.site, s.ml, s.alt) for s in design.swaps] == [(2, "AAA", "AGA")]
    masked = library.alt_all("ATG---", columns, threshold=0.2, width=3)
    assert masked.sequence == "ATG---" and masked.unaligned == "ATG"


def test_split_columns_refuses_a_ragged_row() -> None:
    assert library.split_columns("ATGAAA", 3) == ["ATG", "AAA"]
    with pytest.raises(ValueError, match="whole columns"):
        library.split_columns("ATGAA", 3)


# --- The library proper (L3) ---------------------------------------------------------


def target(
    template: str,
    *columns: dict[str, float],
    combination: int = 1,
    node: str = "Node15",
    gap: str = "",
    weight: float = 0.0,
    width: int = 1,
) -> library.Target:
    source = library.Source(combination=combination, node=node, gap=gap, weight=weight)
    return library.Target(source=source, template=template, dists=dists(*columns), width=width)


def test_generate_makes_ml_altall_and_single_site_alternates() -> None:
    one = target("MKK", {"M": 0.95, "L": 0.05}, {"K": 0.6, "R": 0.4}, {"K": 0.7, "Q": 0.3})
    made = library.generate([one], library.Axes(threshold=0.2, max_substitutions=1))
    assert [(c.variant, c.kind, c.sequence) for c in made] == [
        ("ML", "ml", "MKK"),
        ("ALT", "altall", "MRQ"),
        ("K2R", "single", "MRK"),
        ("K3Q", "single", "MKQ"),
    ]
    # The site in a name is the alignment column, not the residue number.
    assert library.swap_label(library.Swap(site=300, ml="K", alt="R", pp=0.3, ml_pp=0.7)) == "K300R"


def test_generate_counts_what_it_will_make_before_making_it() -> None:
    one = target("MKK", {"M": 0.95, "L": 0.05}, {"K": 0.6, "R": 0.4}, {"K": 0.7, "Q": 0.3})
    for axes in (
        library.Axes(),
        library.Axes(alt_all=False),
        library.Axes(samples=5),
        library.Axes(max_substitutions=1),
        library.Axes(max_substitutions=2),
        library.Axes(alt_all=False, samples=3, max_substitutions=2, threshold=0.35),
    ):
        assert library.count_candidates([one], axes) == len(library.generate([one], axes))


def test_generate_refuses_more_than_the_candidate_limit() -> None:
    wide = target("K" * 30, *[{"K": 0.6, "R": 0.4}] * 30)
    # 30 ambiguous sites: 4,060 triples, but every subset up to 6 is far over the limit.
    with pytest.raises(ValueError, match="over the 20000 limit"):
        library.generate([wide], library.Axes(max_substitutions=6))


def test_posterior_samples_are_reproducible_and_stable_as_the_count_grows() -> None:
    one = target("MKK", {"M": 0.6, "L": 0.4}, {"K": 0.6, "R": 0.4}, {"K": 0.5, "Q": 0.5})
    few = library.generate([one], library.Axes(samples=3), seed=11)
    many = library.generate([one], library.Axes(samples=20), seed=11)
    other = library.generate([one], library.Axes(samples=3), seed=12)
    drawn = [c.sequence for c in few if c.kind == "sample"]
    assert drawn == [c.sequence for c in many if c.kind == "sample"][:3]
    assert drawn == [
        c.sequence
        for c in library.generate([one], library.Axes(samples=3), seed=11)
        if c.kind == "sample"
    ]
    assert drawn != [c.sequence for c in other if c.kind == "sample"]
    assert [c.variant for c in few if c.kind == "sample"] == ["S1", "S2", "S3"]


def test_posterior_samples_keep_masked_columns_as_gaps() -> None:
    masked = target("M-K", {"M": 0.6, "L": 0.4}, {"K": 0.6, "R": 0.4}, {"K": 0.5, "Q": 0.5})
    for candidate in library.generate([masked], library.Axes(samples=10), seed=1):
        assert candidate.sequence[1] == "-"


def test_samples_do_not_shift_when_another_target_joins_the_library() -> None:
    columns = ({"M": 0.6, "L": 0.4}, {"K": 0.6, "R": 0.4})
    first = target("MK", *columns, node="Node15")
    second = target("MK", *columns, node="Node20")
    axes = library.Axes(samples=4)
    alone = [c.sequence for c in library.generate([first], axes, seed=7) if c.kind == "sample"]
    together = [
        c.sequence
        for c in library.generate([first, second], axes, seed=7)
        if c.kind == "sample" and c.source.node == "Node15"
    ]
    assert alone == together


def test_reduce_collapses_the_same_protein_and_counts_the_choices() -> None:
    # Two combinations whose ML ancestors agree on MK, with MR as their AltAll, and a
    # third that reconstructs MR outright and has no plausible alternative.
    agree = [
        target("MK", {"M": 0.95, "L": 0.05}, {"K": 0.6, "R": 0.4}, combination=k, weight=0.4)
        for k in (1, 2)
    ]
    differs = target("MR", {"M": 0.95, "L": 0.05}, {"R": 0.9, "Q": 0.1}, combination=3, weight=0.2)
    members = library.reduce(library.generate([*agree, differs], library.Axes()))
    assert [m.unaligned for m in members] == ["MR", "MK"]  # support, then weight
    # MR is combination 3's ML ancestor and the other two's AltAll: three choices, and it
    # takes the best name among them.
    shared = members[0]
    assert shared.support == 3 and shared.weight == pytest.approx(1.0)
    assert shared.kind == "ml" and shared.variant == "ML" and shared.source.combination == 3
    assert sorted({c.variant for c in shared.candidates}) == ["ALT", "ML"]
    # Combination 3's own ML and AltAll agree, so that one choice is still counted once.
    assert len([c for c in shared.candidates if c.source.combination == 3]) == 2
    assert members[1].support == 2 and members[1].weight == pytest.approx(0.8)


def test_reduce_counts_one_source_once_however_many_variants_reach_it() -> None:
    # A two-state column: samples land on the ML sequence again and again.
    one = target("MK", {"M": 1.0}, {"K": 0.6, "R": 0.4}, weight=0.5)
    members = library.reduce(library.generate([one], library.Axes(samples=30), seed=3))
    ml = next(m for m in members if m.unaligned == "MK")
    assert ml.support == 1 and ml.weight == pytest.approx(0.5)
    assert len(ml.candidates) > 1  # many paths, one choice


def test_reduce_treats_the_same_protein_from_two_gap_versions_as_one_member() -> None:
    columns = ({"M": 0.95, "L": 0.05}, {"K": 0.6, "R": 0.4}, {"Q": 0.9, "E": 0.1})
    original = target("MKQ", *columns, gap="")
    masked = target("MKQ", *columns, gap="fitch")
    members = library.reduce(library.generate([original, masked], library.Axes(alt_all=False)))
    assert len(members) == 1
    assert members[0].support == 2 and members[0].gaps == ["", "fitch"]


def test_member_statistics() -> None:
    one = target("MKQ", {"M": 0.9, "L": 0.1}, {"K": 0.6, "R": 0.4}, {"Q": 0.8, "E": 0.2})
    members = library.reduce(library.generate([one], library.Axes()))
    ml = next(m for m in members if m.kind == "ml")
    altall = next(m for m in members if m.kind == "altall")
    assert ml.length == 3 and ml.differences(ml.sequence) == []
    assert altall.differences(ml.sequence) == ["K2R", "Q3E"]
    # The ML ancestor's states are the most probable by construction, AltAll's are not.
    assert ml.mean_posterior(one) == pytest.approx((0.9 + 0.6 + 0.8) / 3)
    assert altall.mean_posterior(one) == pytest.approx((0.9 + 0.4 + 0.2) / 3)
    assert altall.mean_posterior(one) < ml.mean_posterior(one)


def test_member_differences_refuse_a_mismatched_reference() -> None:
    one = target("MK", {"M": 0.9, "L": 0.1}, {"K": 0.6, "R": 0.4})
    member = library.reduce(library.generate([one], library.Axes(alt_all=False)))[0]
    with pytest.raises(ValueError, match="cannot be compared"):
        member.differences("MKQ")


def test_member_names_follow_the_agreed_scheme() -> None:
    source = library.Source(combination=1, node="Node15", gap="fitch")
    assert library.member_name("ime2", source, "ALT") == "IME2-C1-N15-ALT-F"
    assert library.member_name("ime2", source, "K300R") == "IME2-C1-N15-K300R-F"
    # The original reconstruction is the one segment that can be absent.
    plain = library.Source(combination=3, node="Node15")
    assert library.member_name("ime2", plain, "S42") == "IME2-C3-N15-S42"
    # Every gap method has its letter, and an unknown one still gives a usable name.
    for method, letter in [("parsimony", "P"), ("ml", "M"), ("prank", "K"), ("grasp", "G")]:
        named = library.member_name("x", library.Source(1, "Node1", method), "ML")
        assert named == f"X-C1-N1-ML-{letter}"
    # codeml reconstructions name their model, so two models never share a name.
    names = {
        model: library.member_name("ime2", library.Source(1, "Node15", f"codeml-{model}"), "ML")
        for model in ["M0", "M2a", "M3", "FR"]
    }
    assert names == {
        "M0": "IME2-C1-N15-ML-DM0",
        "M2a": "IME2-C1-N15-ML-DM2A",
        "M3": "IME2-C1-N15-ML-DM3",
        "FR": "IME2-C1-N15-ML-DFR",
    }
    # An unknown reconstruction is refused, never given a guessed letter.
    for version in ["newmethod", "codeml-", "codeml-!!"]:
        with pytest.raises(ValueError, match="no member-name segment"):
            library.member_name("x", library.Source(1, "Node1", version), "ML")
    # Codes are cleaned up and capped; an unusable one falls back rather than failing.
    assert library.member_name("my job 2024!", plain, "ML").startswith("MYJOB2-")
    assert library.member_name("...", plain, "ML") == "LIB-C3-N15-ML"
    assert library.member_name("x", library.Source(1, "AncKLSC"), "ML") == "X-C1-NANCKLSC-ML"


def test_name_fits_flags_the_cap_without_renaming() -> None:
    source = library.Source(combination=12, node="Node1500", gap="grasp")
    assert library.name_fits(library.member_name("IME2AB", source, "ALT"))
    # A wide alignment can push a single-site name over the cap: the page warns.
    long_name = library.member_name("IME2AB", source, "K12345R")
    assert len(long_name) > library.NAME_LIMIT and not library.name_fits(long_name)
    assert not library.name_fits("ime2-c1")  # lower case is not in the character set
    assert not library.name_fits("")


def members_for(*templates: str, support: int = 1) -> list[library.Member]:
    """Members straight from sequences, for the selector's own tests: one source each,
    all in one alignment, kinds as given by the first letter of the template's role."""
    made: list[library.Member] = []
    for index, template in enumerate(templates, start=1):
        source = library.Source(combination=index, node="Node15", weight=0.1 * index)
        candidate = library.Candidate(template, f"S{index}", "sample", source)
        made.append(
            library.Member(
                sequence=template,
                kind="sample",
                variant=f"S{index}",
                source=source,
                candidates=[candidate] * support,
            )
        )
    return made


def test_select_keeps_the_controls_and_respects_the_budget() -> None:
    one = target("MKQ", {"M": 0.9, "L": 0.1}, {"K": 0.6, "R": 0.4}, {"Q": 0.7, "E": 0.3})
    members = library.reduce(library.generate([one], library.Axes(samples=20), seed=5))
    picked = library.select(members, budget=4)
    assert len(picked.members) == 4
    assert [m.kind for m in picked.members[:2]] == ["ml", "altall"]
    assert picked.dropped == len(members) - 4
    # The budget is wells, so it is never exceeded, not even by the controls.
    assert len(library.select(members, budget=1).members) == 1
    # controls=False reserves no places, so the picks are just the top of reduce's
    # order; the ML sequence is still picked where it earns it (samples land on it).
    plain = library.select(members, budget=3, controls=False)
    assert [m.unaligned for m in plain.members] == [m.unaligned for m in members[:3]]


def test_select_orders_controls_by_fit_then_the_default_gap_method() -> None:
    """The first well holds the ancestor the user would have ordered: the best-fitting
    combination's, masked by the job's default gap method. The old order put the
    unmasked reconstruction (gap "") first, -0.188 identity on ime2/indels (L7)."""

    def ml(sequence: str, combination: int, gap: str, weight: float) -> library.Member:
        source = library.Source(combination=combination, node="Node2", gap=gap, weight=weight)
        candidate = library.Candidate(sequence, "ML", "ml", source)
        return library.Member(sequence, "ml", "ML", source, candidates=[candidate])

    members = [
        ml("MKQ-", 1, "", 0.3),
        ml("MK--", 1, "parsimony", 0.3),
        ml("M---", 1, "fitch", 0.3),
        ml("LKQ-", 2, "", 0.7),
        ml("LK--", 2, "parsimony", 0.7),
        ml("L---", 2, "fitch", 0.7),
    ]
    picked = library.select(members, budget=6, preferred_gaps=["fitch"]).members
    assert [m.unaligned for m in picked] == ["L", "LK", "LKQ", "M", "MK", "MKQ"]
    # Without an outgroup the default is shared presence; the original stays last.
    picked = library.select(members, budget=3, preferred_gaps=["parsimony"]).members
    assert [m.unaligned for m in picked] == ["LK", "L", "LKQ"]
    # With no preference the masks still go ahead of the unmasked original.
    assert library.select(members, budget=1).members[0].source.gap != ""


def test_select_most_supported_takes_reduce_s_order() -> None:
    picked = library.select(members_for("AAAA", "AAAB", "AAAC"), budget=2, controls=False)
    assert [m.sequence for m in picked.members] == ["AAAA", "AAAB"]
    assert picked.strategy == "supported"


def test_select_spanning_covers_the_disagreement_rather_than_near_identical_variants() -> None:
    # Three sequences one step apart and one far away: spanning must reach the far one.
    members = members_for("AAAA", "AAAB", "AAAC", "DDDD")
    spanning = library.select(members, budget=2, strategy="spanning", controls=False)
    assert [m.sequence for m in spanning.members] == ["AAAA", "DDDD"]
    supported = library.select(members, budget=2, controls=False)
    assert [m.sequence for m in supported.members] == ["AAAA", "AAAB"]


def test_select_spanning_spans_each_alignment_of_a_mixed_library() -> None:
    # Targets from two alignments have rows of different widths, which `distance`
    # refuses to compare: spanning used to raise ValueError from inside `select`.
    members = members_for("AAAA", "AAAB", "DDDD", "CCCCCC", "CCCCCD", "EEEEEE")
    spanning = library.select(members, budget=4, strategy="spanning", controls=False)
    # Both alignments are represented, each spanned on its own columns.
    assert [m.sequence for m in spanning.members] == ["AAAA", "CCCCCC", "DDDD", "EEEEEE"]
    # A budget that cannot reach every alignment still never compares across them.
    assert len(library.select(members, budget=2, strategy="spanning", controls=False).members) == 2


def test_select_per_node_gives_each_node_a_share() -> None:
    # Distinct ancestors, as nodes along a lineage are: identical ones would reduce to
    # one member and one well, which the per-node share then cannot split.
    targets = [
        target(template, *columns, node=node)
        for node, template, columns in [
            ("Node15", "MKQ", ({"M": 0.9, "L": 0.1}, {"K": 0.6, "R": 0.4}, {"Q": 0.7, "E": 0.3})),
            ("Node20", "WYF", ({"W": 0.9, "C": 0.1}, {"Y": 0.6, "H": 0.4}, {"F": 0.7, "V": 0.3})),
            ("Node25", "DEG", ({"D": 0.9, "N": 0.1}, {"E": 0.6, "Q": 0.4}, {"G": 0.7, "A": 0.3})),
        ]
    ]
    members = library.reduce(library.generate(targets, library.Axes(samples=10), seed=2))
    picked = library.select(members, budget=6, strategy="per_node", controls=False)
    nodes = [m.source.node for m in picked.members]
    assert len(picked.members) == 6
    # Round robin, so no node is served twice before another is served once.
    assert sorted(nodes) == ["Node15", "Node15", "Node20", "Node20", "Node25", "Node25"]


def test_select_refuses_an_unknown_strategy() -> None:
    with pytest.raises(ValueError, match="unknown strategy"):
        library.select(members_for("AAAA"), budget=2, strategy="whatever")


def test_distance_counts_columns_and_refuses_rows_of_different_widths() -> None:
    assert library.distance("MKQ", "MRQ") == 1
    assert library.distance("AAAAAA", "AAACCC", width=3) == 1  # one codon, three letters
    with pytest.raises(ValueError, match="cannot be compared"):
        library.distance("MKQ", "MK")


def test_coverage_states_what_the_set_achieves() -> None:
    one = target("MKQ", {"M": 0.95, "L": 0.05}, {"K": 0.6, "R": 0.4}, {"Q": 0.7, "E": 0.3})
    axes = library.Axes(threshold=0.2)
    members = library.reduce(library.generate([one], axes))
    # ML and AltAll together carry both states at both ambiguous sites.
    full = library.coverage(library.select(members, budget=2).members, [one])
    assert full.sites == 2 and full.total == 4 and full.complete and not full.missing
    assert "all 2 ambiguous sites" in full.sentence
    # The ML ancestor alone misses the alternative at each of them.
    alone = library.coverage([m for m in members if m.kind == "ml"], [one])
    assert alone.covered == 2 and alone.total == 4
    assert alone.missing == [("Node15", 2, "R"), ("Node15", 3, "E")]
    assert alone.fraction == pytest.approx(0.5) and "2 of 4 states" in alone.sentence


def test_coverage_skips_masked_columns_and_unambiguous_sites() -> None:
    masked = target("M-Q", {"M": 0.95, "L": 0.05}, {"K": 0.6, "R": 0.4}, {"Q": 0.99, "E": 0.01})
    members = library.reduce(library.generate([masked], library.Axes()))
    # Column 2 is a gap (no state to cover), columns 1 and 3 have one plausible state.
    nothing = library.coverage(members, [masked])
    assert nothing.sites == 0 and nothing.total == 0 and nothing.complete
    assert "nothing is ambiguous" in nothing.sentence


def test_coverage_ignores_members_from_another_alignment() -> None:
    one = target("MKQ", {"M": 0.95, "L": 0.05}, {"K": 0.6, "R": 0.4}, {"Q": 0.7, "E": 0.3})
    members = library.reduce(library.generate([one], library.Axes()))
    elsewhere = members_for("MRE-")  # a wider alignment: it cannot vouch for this one
    mixed = library.coverage([*members, *elsewhere], [one])
    assert mixed.covered == library.coverage(members, [one]).covered


def test_the_library_works_on_a_codon_reconstruction() -> None:
    codons = target(
        "ATGAAACAG",
        {"ATG": 0.95, "CTG": 0.05},
        {"AAA": 0.6, "AGA": 0.4},
        {"CAG": 0.7, "GAG": 0.3},
        width=3,
    )
    members = library.reduce(library.generate([codons], library.Axes(samples=5), seed=4), width=3)
    ml = next(m for m in members if m.kind == "ml")
    altall = next(m for m in members if m.kind == "altall")
    assert ml.sequence == "ATGAAACAG" and ml.length == 3  # codons, not nucleotides
    assert altall.sequence == "ATGAGAGAG"
    assert altall.differences(ml.sequence) == ["AAA2AGA", "CAG3GAG"]
    # Samples are whole codons, never a mixture of two.
    for member in members:
        assert all(
            codon in {"ATG", "CTG", "AAA", "AGA", "CAG", "GAG"}
            for codon in library.split_columns(member.sequence, 3)
        )
    cov = library.coverage(library.select(members, budget=2).members, [codons])
    assert cov.complete and cov.sites == 2


# --- The saved spec (L4) ---------------------------------------------------------------


def entry(node: str, version: str = "original") -> dict[str, object]:
    return library.target_entry("mafft/iqtree/LG", node, version, [["c", "d"], ["a"], ["b"]], 2)


def test_target_entries_store_a_sorted_partition() -> None:
    stored = library.target_entry("k", "Node2", "fitch", [{"d", "c"}, {"b"}, {"a"}], 2)
    assert stored == {
        "combination": "k",
        "node": "Node2",
        "version": "fitch",
        "partition": [["a"], ["b"], ["c", "d"]],  # JSON, and the same however it was built
        "below": 2,
    }


def test_add_targets_skips_what_is_already_there() -> None:
    spec = library.new_spec(seed=7)
    spec, added = library.add_targets(spec, [entry("Node1"), entry("Node2")])
    assert added == 2
    # A second click, or a lineage overlapping a node already added, changes nothing...
    spec, added = library.add_targets(spec, [entry("Node2"), entry("Node3")])
    assert added == 1
    assert [t["node"] for t in library.spec_targets(spec)] == ["Node1", "Node2", "Node3"]
    # ...but the same node in another reconstruction is another target.
    spec, added = library.add_targets(spec, [entry("Node2", "fitch")])
    assert added == 1 and len(library.spec_targets(spec)) == 4
    assert spec["seed"] == 7


def test_remove_target_goes_by_key_and_leaves_the_rest() -> None:
    spec, _ = library.add_targets(library.new_spec(1), [entry("Node1"), entry("Node2")])
    spec = library.remove_target(spec, ("mafft/iqtree/LG", "Node1", "original"))
    assert [t["node"] for t in library.spec_targets(spec)] == ["Node2"]
    # An out-of-date page naming a target that is gone removes nothing.
    assert library.remove_target(spec, ("mafft/iqtree/LG", "Node1", "original")) == spec


def test_spec_targets_tolerates_a_damaged_spec() -> None:
    assert library.spec_targets({}) == []
    assert library.spec_targets({"targets": "nonsense"}) == []
    assert library.spec_targets({"targets": [entry("Node1"), "junk", 3]}) == [entry("Node1")]


def test_default_code() -> None:
    assert library.default_code("Ime2 kinases", 1) == "IME2KI"
    assert library.default_code("", 3) == "LIB3"
    assert library.default_code("!!!", 12) == "LIB12"


def test_design_round_trips_through_the_spec() -> None:
    spec = library.new_spec(seed=7)
    assert library.design_from(spec).budget == library.DEFAULT_BUDGET
    # A library saved before L5a has no design: it gets the page's defaults, samples
    # included (10, the owner's decision), while `Axes` itself still defaults to none.
    assert library.design_from(spec).axes.samples == library.DEFAULT_SAMPLES
    assert library.Axes().samples == 0
    design = library.Design(
        axes=library.Axes(alt_all=False, threshold=0.3, samples=25, max_substitutions=1),
        budget=96,
        strategy="spanning",
    )
    saved = library.with_design(spec, design)
    assert saved["targets"] == [] and saved["seed"] == 7  # the targets are untouched
    back = library.design_from(saved)
    assert back == design
    # The page's threshold wins over the stored one: it is a View control.
    assert library.design_from(saved, threshold=0.1).axes.threshold == 0.1


def test_design_from_refuses_nonsense_rather_than_failing() -> None:
    spec = {
        "axes": {"samples": "lots", "max_substitutions": 99, "threshold": None},
        "budget": -4,
        "strategy": "whatever",
    }
    design = library.design_from(spec)
    assert design.axes.samples == 0  # unusable: the low end of the range
    assert design.axes.max_substitutions == library.MAX_SUBSTITUTIONS  # clamped
    assert design.axes.threshold == library.ALT_ALL_THRESHOLD
    assert design.budget == 1 and design.strategy == library.DEFAULT_STRATEGY
    assert library.clamp_int("12", 0, 100) == 12 and library.clamp_int(None, 2, 100) == 2


def test_select_spanning_is_the_same_set_however_it_is_computed() -> None:
    """Spanning keeps a running distance to what it has already chosen rather than
    recomputing it at every pick (a 629-member library took fifteen seconds). The set
    must not change: this is the slow definition, written out."""
    members = members_for("AAAA", "AAAB", "AAAC", "DDDD", "ADDA", "DADA", "AADD")

    def slowly(budget: int) -> list[str]:
        chosen: list[str] = []
        rest = [m.sequence for m in members]
        while rest and len(chosen) < budget:
            far = [min((library.distance(r, c) for c in chosen), default=inf) for r in rest]
            best = max(range(len(rest)), key=lambda i: (far[i], -i))
            chosen.append(rest.pop(best))
        return chosen

    for budget in range(1, len(members) + 1):
        picked = library.select(members, budget=budget, strategy="spanning", controls=False)
        assert [m.sequence for m in picked.members] == slowly(budget)


def test_coverage_curve_matches_coverage_at_every_budget() -> None:
    """The curve is what the page draws instead of running the selector once per budget,
    so it has to agree with `coverage` at each of them."""
    one = target("MKQ", {"M": 0.9, "L": 0.1}, {"K": 0.6, "R": 0.4}, {"Q": 0.5, "E": 0.5})
    two = target("MKQ", {"M": 0.5, "L": 0.5}, {"K": 0.6, "R": 0.4}, {"Q": 0.9, "E": 0.1})
    targets = [one, dataclasses.replace(two, source=library.Source(2, "Node2"))]
    axes = library.Axes(samples=6, max_substitutions=1)
    members = library.reduce(library.generate(targets, axes, seed=3))
    for strategy in library.STRATEGIES:
        picked = library.select(members, budget=len(members), strategy=strategy)
        curve = library.coverage_curve(picked.members, targets)
        assert curve[0] == 0 and len(curve) == len(picked.members) + 1
        assert curve == sorted(curve)  # adding a member never uncovers a state
        for budget in (1, 2, 5, len(picked.members)):
            at = library.select(members, budget=budget, strategy=strategy)
            assert curve[len(at.members)] == library.coverage(at.members, targets).covered
    # With nothing ambiguous there is nothing to cover, at any budget.
    certain = target("MK", {"M": 1.0}, {"K": 1.0})
    flat = library.reduce(library.generate([certain], library.Axes()))
    assert set(library.coverage_curve(flat, [certain])) == {0}


def test_strategy_choices_never_hide_the_strategy_in_use() -> None:
    """The page loads with `supported` and keeps the other two behind "Advanced" (owner
    decision 2026-10-13, after L7 found no accuracy case for them). A library saved with
    one of those must open the disclosure, or its own setting would be invisible."""
    keys = {s.key for s in library.STRATEGY_CHOICES}
    assert keys == set(library.STRATEGIES), "every strategy needs a presentation entry"
    plain, advanced, opened = library.strategy_choices(library.DEFAULT_STRATEGY)
    assert [s.key for s in plain] == [library.DEFAULT_STRATEGY]
    assert not opened, "the default must not open the disclosure"
    for key in (s.key for s in advanced):
        shown, hidden, open_now = library.strategy_choices(key)
        assert open_now, f"a library saved with {key} must open Advanced"
        assert key in {s.key for s in (*shown, *hidden)}
    # Every strategy is still selectable: Advanced is editorial, not a capability.
    assert {s.key for s in (*plain, *advanced)} == set(library.STRATEGIES)
