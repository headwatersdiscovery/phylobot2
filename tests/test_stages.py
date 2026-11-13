from datetime import timedelta
from typing import Any

import pytest
from django.test import Client
from django.utils import timezone

from phylobot import jobs, stages
from phylobot.executors import expected_tasks
from phylobot.models import Job
from phylobot.modes import local_user
from tests.conftest import FakeStorage

FULL = {"aligners": ["mafft", "muscle"], "outgroup": ["s1"], "gap_methods": ["ml", "prank"]}


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
        {"seqtype": "codon", "aligners": ["mafft", "prank_codon"], "alignment": ["a1"]},
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
        {"uncertainty_only": True, "uncertainty_of": {"combinations": 9}},
        {"reconcile": True, "tree_builders": ["iqtree", "raxml"], "outgroup": ["s1"]},
        {"reconcile": True, "stop_after": "trees", "tree": ["t1"]},
        {"reconcile": True, "stop_after": "alignments"},
        {"reconcile": True, "tree_builders": ["iqtree", "raxml", "generax"], "outgroup": ["s1"],
         "models": ["LG+G4", "WAG+G4"], "fixed_models": True},
        {"reconcile": True, "tree_builders": ["generax"], "tree": ["t1"], "stop_after": "trees",
         "seqtype": "codon", "codon_models": ["GY+F3X4+G4"]},
    ],
)  # fmt: skip
def test_expected_kinds_add_up_to_expected_tasks(params: dict[str, Any]) -> None:
    assert sum(stages.expected_kinds(params).values()) == expected_tasks(params)


def test_codon_fits_on_amino_acid_trees_counted_as_the_pipeline_ran_them() -> None:
    # Lysin (2026-10-02): MAFFT, IQ-TREE + RAxML-NG, ModelFinder + LG+G4, GY+F3X4+G4 fitted
    # on the best amino-acid tree of each source; the trace's task counts.
    params = {"seqtype": "codon", "tree_builders": "iqtree,raxml", "fixed_models": True,
              "models": "LG+G4", "codon_models": "GY+F3X4+G4", "codon_trees": "amino",
              "gap_methods": "parsimony"}  # fmt: skip
    assert stages.expected_kinds(params) == {
        "VALIDATE": 1, "ALIGN_MAFFT": 1, "CODON_ALIGNMENT": 1, "EXPORT_ALIGNMENT": 1,
        "MODELFINDER": 1, "TREE_IQTREE": 2, "TREE_RAXML": 2, "FIT_IQTREE": 2,
        "FIT_CODON_AA": 2, "ANCESTRAL_IQTREE": 6, "EXPORT": 6, "SUMMARY": 1,
        "ANCESTOR_MAP": 1,  # added 2026-10-30, after that run
    }  # fmt: skip


def test_edges_drop_those_implied_by_a_longer_path() -> None:
    names = list(stages.expected_kinds(FULL))
    links = set(stages.edges(names))
    assert ("ALIGN_MAFFT", "MODELFINDER") in links and ("ROOT_TREE", "GAPS_PRANK") in links
    assert ("MODELFINDER", "TREE_IQTREE") in links
    assert ("ALIGN_MAFFT", "TREE_IQTREE") not in links  # via ModelFinder
    assert ("VALIDATE", "MODELFINDER") not in links  # via the aligners
    assert ("ANCESTRAL_IQTREE", "EXPORT") not in links  # via root / gaps
    assert ("EXPORT", "SUMMARY") in links and ("VALIDATE", "SUMMARY") not in links
    # without aligners, provided alignments come straight from VALIDATE
    provided = stages.edges(list(stages.expected_kinds({"aligners": "none", "alignment": ["a"]})))
    assert ("VALIDATE", "MODELFINDER") in provided
    fixed = {"aligners": "none", "alignment": ["a"], "modelfinder": False, "fixed_models": True,
             "models": ["LG+G4"]}  # fmt: skip
    assert ("VALIDATE", "TREE_IQTREE") in stages.edges(list(stages.expected_kinds(fixed)))
    # alignment-only: SUMMARY waits for the aligners directly
    alignments = stages.edges(list(stages.expected_kinds({"stop_after": "alignments"})))
    assert set(alignments) == {("VALIDATE", "ALIGN_MAFFT"), ("ALIGN_MAFFT", "EXPORT_ALIGNMENT"),
                               ("ALIGN_MAFFT", "SUMMARY")}  # fmt: skip


def iso(minutes: float) -> str:
    return (timezone.now() - timedelta(minutes=minutes)).isoformat()


def done(name: str, start: float, end: float) -> dict[str, str]:
    return {"name": name, "status": "COMPLETED", "realtime": "1s", "start": iso(start),
            "complete": iso(end)}  # fmt: skip


def test_graph_status_and_times() -> None:
    progress = {"tasks": [
        done("VALIDATE (1)", 10, 9.5),
        done("ALIGN_MAFFT (mafft)", 9.5, 8),
        done("ALIGN_MUSCLE (muscle)", 9.5, 7),
        done("TREE_IQTREE (mafft/modelfinder)", 7, 2),
        {"name": "TREE_IQTREE (muscle/modelfinder)", "status": "RUNNING", "start": iso(6)},
    ]}  # fmt: skip
    graph = stages.graph(FULL, progress, finished=False)
    by_name = {s.name: s for s in graph["stages"]}
    assert by_name["ALIGN_MAFFT"].status == "succeeded"
    assert by_name["ALIGN_MAFFT"].elapsed == pytest.approx(90, abs=1)
    tree = by_name["TREE_IQTREE"]
    assert (tree.status, tree.count, tree.completed) == ("running", "1/2", None)
    assert tree.elapsed == pytest.approx(7 * 60, abs=2)  # since its first task started
    assert by_name["EXPORT"].status == "waiting" and by_name["EXPORT"].elapsed is None
    assert {link.state for link in graph["links"]} == {"done", "flowing", ""}
    # nothing overlaps, and everything fits the drawing
    boxes = [(s.x, s.y) for s in graph["stages"]]
    assert len(set(boxes)) == len(boxes)
    assert max(x for x, _ in boxes) + stages.NODE_W <= graph["width"]
    assert max(y for _, y in boxes) + stages.NODE_H <= graph["height"]
    assert by_name["VALIDATE"].y == by_name["TREE_IQTREE"].y == by_name["SUMMARY"].y  # main line

    # once the job has ended: unfinished work is stopped, the rest never ran
    failed = stages.graph(FULL, progress | {"tasks": [*progress["tasks"], {
        "name": "ALIGN_MUSCLE (x)", "status": "FAILED"}]}, finished=True)  # fmt: skip
    by_name = {s.name: s for s in failed["stages"]}
    assert by_name["ALIGN_MUSCLE"].status == "failed"
    assert by_name["TREE_IQTREE"].status == "stopped" and by_name["EXPORT"].status == "not run"


def test_stage_groups_name_the_pipeline_stage() -> None:
    # Skins colour the steps by stage; rooting counts as tree work there.
    graph = stages.graph(FULL, {}, finished=False)
    groups = {stage.name: stage.group for stage in graph["stages"]}
    assert groups["VALIDATE"] == "setup" and groups["ALIGN_MAFFT"] == "align"
    assert groups["MODELFINDER"] == "modelfinder" and groups["TREE_IQTREE"] == "trees"
    assert groups["ROOT_TREE"] == "trees" and groups["ANCESTRAL_IQTREE"] == "ancestors"


def test_graph_for_jobs_without_timestamps() -> None:
    tasks = [{"name": "VALIDATE (1)", "status": "COMPLETED", "realtime": "771ms"},
             {"name": "NEW_STEP (1)", "status": "COMPLETED", "realtime": "2s"}]  # fmt: skip
    graph = stages.graph({}, {"tasks": tasks}, finished=True)
    by_name = {s.name: s for s in graph["stages"]}
    assert by_name["VALIDATE"].elapsed == 0.771 and by_name["VALIDATE"].completed is None
    assert by_name["NEW_STEP"].label == "New_Step"  # unknown stages are still drawn


@pytest.mark.django_db
def test_job_page_steps_panel_open_only_while_running(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = Job.objects.create(owner=local_user(), params={})
    assert jobs.mark(job.id, Job.Status.RUNNING)
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert 'id="job-stages" open' in page and page.count('class="stage waiting"') == 9
    assert 'data-scheme="simple"' in page and '<option value="safe">' in page  # colour menu
    assert jobs.mark(job.id, Job.Status.SUCCEEDED)
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert 'id="job-stages">' in page and "Check input" in page


def test_generax_tree_source_counts_match_a_real_run() -> None:
    """Ime2, MAFFT, IQ-TREE + RAxML-NG + GeneRax, LG+G4, outgroup: the hand run of
    2026-10-17 submitted 26 tasks, before ANCESTOR_MAP (one more since 2026-10-30)."""
    params = {
        "reconcile": True, "tree_builders": ["iqtree", "raxml", "generax"], "modelfinder": False,
        "models": ["LG+G4"], "fixed_models": True, "outgroup": ["Yarrowia.lipolytica.IME2"],
    }  # fmt: skip
    kinds = stages.expected_kinds(params)
    assert expected_tasks(params) == 27
    assert kinds["TREE_GENERAX"] == kinds["FIT_GENERAX"] == kinds["TOPOLOGY_TEST"] == 1
    assert kinds["RECONCILE"] == kinds["ANCESTRAL_IQTREE"] == 3


def test_depends_is_in_topological_order() -> None:
    """`ORDER = list(DEPENDS)` is relied on as a topological order, and `layout` reads a
    stage's parents out of a dict it fills as it goes — so a stage listed *before* one it
    depends on raises KeyError and the job page 500s. `PROPERTIES` shipped three lines
    above its own `EXPORT` on 2026-10-20, breaking every job that ticked the box."""
    seen: set[str] = set()
    for name, parents in stages.DEPENDS.items():
        missing = [p for p in parents if p in stages.DEPENDS and p not in seen]
        assert not missing, f"{name} is listed before its parent(s) {missing}"
        seen.add(name)


def test_the_steps_graph_draws_every_optional_stage() -> None:
    """One job with everything ticked: the graph must lay out, not raise. Each option
    adds stages and edges, and layout's assumptions are easy to break from DEPENDS."""
    params: dict[str, Any] = {
        "aligners": ["mafft", "prank"],
        "tree_builders": ["iqtree", "raxml"],
        "modelfinder": True,
        "fixed_models": True,
        "models": ["LG+G4"],
        "gap_methods": ["parsimony", "fitch", "ml", "prank", "grasp"],
        "outgroup": ["x"],
        "properties": ["charge", "moment", "hydropathy_regions"],
        "stop_after": "ancestors",
    }
    graph = stages.graph(params, {}, finished=False)
    drawn = {s.name: s for s in graph["stages"]}
    assert {"PROPERTIES", "EXPORT", "GAPS_ML", "GAPS_PRANK", "GAPS_GRASP"} <= set(drawn)
    # PROPERTIES reads EXPORT's output, so it is never drawn to the left of it.
    assert drawn["PROPERTIES"].x > drawn["EXPORT"].x


def test_an_extension_counts_only_what_it_runs() -> None:
    params = {"alignment": ["a1"], "alignment_keys": ["mafft"], "aligners": ["muscle"],
              "tree_builders": ["iqtree", "raxml"], "models": ["WAG+G4", "JTT+G4"],
              "fixed_models": True,
              "skip_combinations": ["mafft/iqtree/modelfinder", "mafft/iqtree/WAG+G4"],
              "modelfinder_reused": ["mafft"]}  # fmt: skip
    kinds = stages.expected_kinds(params)
    # MAFFT: iqtree/JTT+G4 and raxml x 3 models; MUSCLE: all 6. ModelFinder on MUSCLE only.
    assert kinds["MODELFINDER"] == 1
    assert kinds["TREE_IQTREE"] == 1 + 3 and kinds["TREE_RAXML"] == 3 + 3
    assert kinds["FIT_IQTREE"] == kinds["TREE_RAXML"]
    assert kinds["ANCESTRAL_IQTREE"] == 10 and kinds["EXPORT_ALIGNMENT"] == 2
    assert "ALIGN_MAFFT" not in kinds and kinds["ALIGN_MUSCLE"] == 1
