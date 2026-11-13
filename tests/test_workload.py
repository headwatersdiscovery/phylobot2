from collections import Counter
from typing import Any

import pytest

from phylobot import catalog, stages, workload
from phylobot.models import Job
from phylobot.modes import local_user
from tests.test_stages import FULL


@pytest.mark.parametrize(
    "params",
    [
        {},
        {"alignment": ["a1", "a2"], "tree": ["t1"]},
        {"alignment_labels": ["a"], "tree_labels": ["t"], "aligners": "none"},
        {"models": ["LG+G4", "WAG"], "fixed_models": True},
        {"models": ["LG+G4"], "fixed_models": True, "modelfinder": False},
        {"tree_builders": ["iqtree", "raxml"], "tree": ["t1"]},
        {"tree_builders": ["raxml"], "stop_after": "trees"},
        FULL,
        {"gap_methods": ["parsimony", "ml", "prank", "grasp"], "outgroup": ["s1"]},
        {"stop_after": "alignments", "aligners": ["mafft", "muscle"], "alignment_labels": ["x"]},
        {"stop_after": "trees", "gap_methods": ["ml"], "outgroup": ["s1"]},
        {"aligners": "mafft,prank", "tree_builders": "iqtree,raxml", "models": "LG,WAG+R4",
         "fixed_models": True, "bootstrap": 1000, "gap_methods": "fitch,grasp"},
        {"seqtype": "codon", "genetic_code": 2, "aligners": ["mafft", "prank_codon"]},
        {"seqtype": "codon", "aligners": ["muscle"], "stop_after": "alignments"},
        {"seqtype": "codon", "codon_models": ["GY+F3X4+G4", "MG+F1X4"], "tree": ["t1"],
         "tree_builders": ["iqtree", "raxml"], "gap_methods": ["ml"]},
        {"seqtype": "codon", "codon_models": ["GY+F3X4"], "modelfinder": False,
         "stop_after": "trees"},
        {"seqtype": "codon", "codon_modelfinder": True, "codeml_models": ["M0", "M8"],
         "tree": ["t1"], "aligners": ["mafft", "muscle"]},
        {"seqtype": "codon", "codon_models": ["GY+F3X4", "MG+F1X4"], "codeml_trees": "all",
         "codeml_models": ["M3", "M2a"], "tree_builders": ["iqtree", "raxml"]},
        {"seqtype": "codon", "codon_models": ["GY+F3X4+G4"], "codon_trees": ["amino"],
         "tree_builders": ["iqtree", "raxml"], "tree": ["t1"], "codeml_models": ["M0"],
         "codeml_trees": "all"},
        {"seqtype": "codon", "codon_modelfinder": True, "codon_trees": "search,amino",
         "codon_amino_trees": "all", "models": ["LG", "WAG"], "fixed_models": True,
         "tree_builders": ["raxml"], "outgroup": ["s1"]},
        # An extension adding RAxML-NG and JTT+G4 to its inherited MAFFT alignment and a
        # new MUSCLE one: the parent's combinations and ModelFinder pick are not rerun.
        {"alignment": ["a1"], "alignment_keys": ["mafft"], "aligners": ["muscle"],
         "tree_builders": ["iqtree", "raxml"], "models": ["WAG+G4", "JTT+G4"],
         "fixed_models": True, "outgroup": ["s1"], "gap_methods": ["ml"],
         "skip_combinations": ["mafft/iqtree/modelfinder", "mafft/iqtree/WAG+G4"],
         "modelfinder_reused": ["mafft"]},
        {"alignment_labels": ["MAFFT", "mine"], "alignment_keys": ["mafft", "provided"],
         "aligners": "none", "tree_builders": ["iqtree", "raxml"], "stop_after": "trees",
         "skip_combinations": ["mafft/iqtree/modelfinder", "provided1/iqtree/modelfinder"]},
        {"simulate": True, "spec": {"simulation": {"seeds": [1, 2, 3]}}},
        {"reconcile": True, "tree_builders": ["iqtree", "raxml"], "outgroup": ["s1"]},
        {"reconcile": True, "stop_after": "trees", "tree": ["t1"]},
        {"reconcile": True, "stop_after": "alignments"},
        {"reconcile": True, "tree_builders": ["iqtree", "raxml", "generax"], "outgroup": ["s1"],
         "models": ["LG+G4", "WAG+G4"], "fixed_models": True},
        {"reconcile": True, "tree_builders": ["generax"], "tree": ["t1"], "stop_after": "trees",
         "seqtype": "codon", "codon_models": ["GY+F3X4+G4"]},
        {"uncertainty_only": True, "uncertainty_of": {"combinations": 9}},
    ],
)  # fmt: skip
def test_tasks_are_the_pipelines_tasks(params: dict[str, Any]) -> None:
    assert Counter(kind for kind, _ in workload.tasks(params)) == stages.expected_kinds(params)


def units(**params: Any) -> float:
    return workload.estimate(params).units


def test_codon_fits_on_amino_acid_trees_cost_less_than_codon_searches() -> None:
    codon = {"seqtype": "codon", "codon_models": ["GY+F3X4+G4"]}
    search = units(**codon)
    fit = units(**codon, codon_trees=["amino"])
    assert fit < search / 3
    assert units(**codon, codon_trees=["search", "amino"]) > search
    every = {"codon_trees": ["amino"], "codon_amino_trees": "all", "fixed_models": True,
             "models": ["LG+G4", "WAG+G4"]}  # fmt: skip
    assert units(**codon, **every) > units(**codon, **(every | {"codon_amino_trees": "best"}))


def test_options_change_the_work() -> None:
    base = units(models=["LG+G4"], fixed_models=True)
    assert units(models=["LG+G4", "WAG+G4"], fixed_models=True) > base
    assert units(models=["LG"], fixed_models=True) < base  # no rate heterogeneity
    # Branch support had no measurable cost in the benchmark (up to 40 sequences).
    assert units(models=["LG+G4"], fixed_models=True, bootstrap=2000) == base
    assert units(models=["LG+G4"], fixed_models=True, aligners=["mafft", "prank"]) > base
    assert units(models=["LG+G4"], fixed_models=True, stop_after="trees") < base
    # ModelFinder: more ticked candidates cost more; RAxML-NG widens the default set.
    assert units(models=["LG+G4", "WAG+R4"]) > units(models=["LG+G4"])
    assert units(tree_builders=["iqtree", "raxml"]) > units()
    # Every rate and frequency type ticked: IQ-TREE's shortcuts, so fewer candidates
    # (no +I types) cost more than all of them.
    every = [m for m in catalog.all_models() if not m.startswith("mtREV")]
    assert units(candidates=[m for m in every if "+I" not in m]) > units(candidates=every)


def test_estimate_groups_by_stage() -> None:
    estimate = workload.estimate({"gap_methods": ["ml"]})
    keys = [key for key, _, _ in estimate.groups]
    assert keys == ["setup", "align", "modelfinder", "trees", "ancestors"]
    assert estimate.units == pytest.approx(sum(u for _, _, u in estimate.groups))
    assert estimate.minutes is None
    timed = workload.estimate({"gap_methods": ["ml"]}, seconds_per_unit=30)
    assert timed.minutes == pytest.approx(estimate.units / 2)
    assert [key for key, _, _ in workload.estimate({"stop_after": "alignments"}).groups] == [
        "setup",
        "align",
    ]


@pytest.mark.django_db
def test_seconds_per_unit_learns_from_finished_jobs() -> None:
    assert workload.seconds_per_unit() is None
    reference = workload.reference_seconds({})
    for used in (1.0, 2.0, 3.0):  # this server is 2x the reference, by the median
        Job.objects.create(
            owner=local_user(), status=Job.Status.SUCCEEDED, params={}, cpu_seconds=used * reference
        )
    assert workload.seconds_per_unit() == pytest.approx(120)


def test_modelfinder_cost_follows_its_own_candidates() -> None:
    few = units(candidates=["LG+G4"], models=["JTT"], fixed_models=True)
    many = units(candidates=["LG+G4", "WAG+R4", "JTT+F"], models=["JTT"], fixed_models=True)
    assert many > few
    # Jobs from before the separate grid: the models were the candidates.
    assert units(models=["LG+G4"]) == units(candidates=["LG+G4"])


def test_size_of_a_fasta() -> None:
    size = workload.size_of(b">a\nMKT-AY\n>b\nMKTAYW\n>c\nMK\n")
    assert size == workload.Size(3, (5 + 6 + 2) / 3)  # gaps not counted
    assert workload.size_of(b"not fasta") is None and workload.size_of(b"") is None


def test_more_and_longer_sequences_cost_more() -> None:
    small, large = workload.Size(10, 300), workload.Size(80, 300)
    longer = workload.Size(10, 1200)
    params = {"aligners": ["mafft", "prank"], "gap_methods": ["grasp"], "outgroup": ["a"]}
    cost = {s: workload.estimate(params, s).units for s in (small, large, longer)}
    assert cost[large] > 4 * cost[small] and cost[longer] > 1.5 * cost[small]


@pytest.mark.django_db
def test_stored_estimates_take_over_from_guesses() -> None:
    owner = local_user()
    for _ in range(3):  # older jobs: no stored estimate, guessed at the Ime2 size
        Job.objects.create(owner=owner, status=Job.Status.SUCCEEDED, params={},
                           cpu_seconds=10 * workload.reference_seconds({}))  # fmt: skip
    value, jobs = workload.calibration()
    assert value == pytest.approx(600) and jobs == 3
    for _ in range(3):  # jobs with stored estimates: only these count once there are 3
        Job.objects.create(owner=owner, status=Job.Status.SUCCEEDED, params={},
                           work_units=1.0, cpu_seconds=90.0)  # fmt: skip
    value, jobs = workload.calibration()
    assert value == pytest.approx(90) and jobs == 3


def test_a_neural_property_is_priced_per_draw() -> None:
    """The predictors cost per draw, so ticking ESM-2 must show in the
    estimate, and scale with the draw count; the pure-Python three stay a flat step."""
    base = {"properties": ["charge", "moment", "hydropathy_regions"]}
    with_esm = {"properties": [*base["properties"], "esm_likelihood"]}
    cheap = dict(workload.tasks(base))["PROPERTIES"]
    dear = dict(workload.tasks(with_esm))["PROPERTIES"]
    fewer = dict(workload.tasks({**with_esm, "property_draws": 100}))["PROPERTIES"]
    assert dear > 100 * cheap
    assert cheap < fewer < dear
    with_disorder = dict(workload.tasks({"properties": ["charge", "disorder"]}))["PROPERTIES"]
    assert cheap < with_disorder < dear


def test_a_follow_up_run_prices_each_combinations_ancestors() -> None:
    params = {
        "properties_only": True,
        "properties": ["esm_likelihood"],
        "property_draws": 1000,
        "properties_of": {"combinations": {"mafft/iqtree/LG": ["Node1"], "a/b/c": ["N1", "N2"]}},
    }
    costs = [cost for kind, cost in workload.tasks(params) if kind == "PROPERTIES"]
    assert len(costs) == 2
    assert max(costs) > 1.9 * min(costs) - workload.STEPS["PROPERTIES"].at(workload.IME2_SIZE)


def test_a_mixture_costs_more_than_its_gamma_model() -> None:
    gamma = {"fixed_models": True, "models": ["LG+G4"], "modelfinder": False}
    mixture = {**gamma, "models": ["LG+C20+G4"]}
    plain, slow = workload.tasks(gamma), workload.tasks(mixture)
    assert sorted(kind for kind, _ in plain) == sorted(kind for kind, _ in slow)
    cost = lambda rows, kind: sum(c for k, c in rows if k == kind)  # noqa: E731
    assert cost(slow, "TREE_IQTREE") == pytest.approx(cost(plain, "TREE_IQTREE") * 20)
    assert cost(slow, "ANCESTRAL_IQTREE") == pytest.approx(cost(plain, "ANCESTRAL_IQTREE") * 12)
