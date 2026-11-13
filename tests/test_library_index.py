"""The libraries landing page's stored-data summaries (phylobot.library_index, step L5d).

Everything here is pure: a spec, a job's params and a list of export receipts in, the
card's text out. Nothing reads storage or generates members, which is the decision this
page rests on (owner, 2026-10-07).
"""

from dataclasses import dataclass
from typing import Any

from phylobot import library, library_index


def test_state_follows_the_three_steps() -> None:
    empty = library.new_spec(seed=1)
    assert library_index.state(empty, exports=0) == library_index.DRAFT
    assert library_index.state_label(empty, exports=0) == "Draft"
    # Targets alone are still a draft: the design half has not been saved.
    with_targets, _ = library.add_targets(
        empty, [library.target_entry("mafft/iqtree/LG", "Node1", "original", [], 0)]
    )
    assert library_index.state(with_targets, exports=0) == library_index.DRAFT
    designed = dict(with_targets, budget=24)
    assert library_index.state(designed, exports=0) == library_index.DESIGNED
    assert library_index.state_label(designed, exports=0) == "Designed"
    # An export outranks the design, and several are counted.
    assert library_index.state(designed, exports=1) == library_index.EXPORTED
    assert library_index.state_label(designed, exports=1) == "Exported"
    assert library_index.state_label(designed, exports=3) == "Exported ×3"
    # A spec carrying only empty design values has not been designed.
    assert not library_index.designed(dict(with_targets, axes={}, filters=[]))


def test_model_count_counts_the_fixed_models_and_one_pick_per_finder() -> None:
    assert library_index.model_count({}) == 0
    assert library_index.model_count({"modelfinder": True, "candidates": ["LG+G4", "WAG"]}) == 1
    assert library_index.model_count({"fixed_models": True, "models": ["LG+G4", "WAG"]}) == 2
    # Models listed but the grid switched off do not count; codon models do.
    assert library_index.model_count({"fixed_models": False, "models": ["LG+G4"]}) == 0
    both = {
        "modelfinder": True,
        "fixed_models": True,
        "models": ["LG"],
        "codon_models": ["MG+F3X4"],
    }
    assert library_index.model_count(both) == 3


def test_axes_and_inventory_read_the_job_s_own_parameters() -> None:
    params = {
        "aligners": ["mafft", "muscle", "none"],  # "none" is "my own alignment", not an aligner
        "tree_builders": ["iqtree", "raxml"],
        "modelfinder": True,
        "gap_methods": ["fitch", "ml"],
    }
    assert [(a.key, a.count) for a in library_index.axes(params)] == [
        ("aligners", 2),
        ("trees", 2),
        ("models", 1),
        ("gaps", 2),
    ]
    assert (
        library_index.inventory(params) == "2 aligners · 2 tree searches · 1 model · 2 gap methods"
    )
    # An axis with nothing in it is left out rather than shown as zero.
    assert library_index.inventory({"aligners": ["mafft"]}) == "1 aligner"
    assert library_index.inventory({}) == ""
    # `thin` marks an axis with nothing to compare; no notice is drawn from it on this
    # page (the inventory already says "1 aligner"), but the axes controls will want it.
    assert [a.key for a in library_index.axes(params) if a.thin] == ["models"]


@dataclass
class FakeJob:
    params: dict[str, Any]


@dataclass
class FakeLibrary:
    spec: dict[str, Any]
    job: FakeJob

    @property
    def targets(self) -> list[dict[str, Any]]:
        return library.spec_targets(self.spec)


@dataclass
class FakeExport:
    member_count: int


def spec_with(*nodes: str) -> dict[str, Any]:
    spec = library.new_spec(seed=1)
    entries = [library.target_entry("mafft/iqtree/LG", node, "original", [], 0) for node in nodes]
    spec, _ = library.add_targets(spec, entries)
    return spec


def test_card_reports_what_is_stored_and_nothing_more() -> None:
    job = FakeJob({"aligners": ["mafft", "muscle"], "gap_methods": ["fitch"]})
    record = FakeLibrary(spec_with("Node1", "Node2"), job)
    card = library_index.card(record, [])
    assert card.targets == 2 and card.nodes == ["Node1", "Node2"] and card.more_nodes == 0
    assert card.state == library_index.DRAFT
    # The job is a sub-item of the library, not the other way round: one today, a list
    # when cross-job collections arrive.
    assert len(card.jobs) == 1
    assert card.jobs[0].inventory == "2 aligners · 1 gap method"
    # Without an export there is no honest member count: the page never generates.
    assert card.members is None
    exported = library_index.card(record, [FakeExport(24), FakeExport(12)])
    assert exported.members == 24 and exported.state_label == "Exported ×2"


def test_card_shows_the_first_nodes_and_counts_the_rest() -> None:
    job = FakeJob({"aligners": ["mafft"]})
    nodes = [f"Node{n}" for n in range(1, library_index.NODES_SHOWN + 4)]
    card = library_index.card(FakeLibrary(spec_with(*nodes), job), [])
    assert card.nodes == nodes[: library_index.NODES_SHOWN]
    assert card.more_nodes == 3
    # One ancestor per node here, so the target count is the node count.
    assert card.targets == len(nodes)
