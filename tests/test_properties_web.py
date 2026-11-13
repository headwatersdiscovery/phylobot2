"""The New Job page's property-posterior checkbox and what it sends the pipeline
(tier 0 runs inline, with a box so it can be turned off).
"""

import re
from pathlib import Path
from typing import Any

import pytest
from django.test import Client

from phylobot import executors, forms, properties, results, stages, workload
from phylobot.models import Job

ROOT = Path(__file__).resolve().parents[1]
IME2 = (ROOT / "src" / "phylobot" / "examples" / "ime2.fasta").read_text()

pytestmark = pytest.mark.django_db

BASE: dict[str, Any] = {
    "aligners": ["mafft"],
    "tree_builders": ["iqtree"],
    "modelfinder": "on",
    "run_trees": "on",
    "run_ancestors": "on",
    "bootstrap": 1000,
    "seed": 1,
    "sequences_text": IME2,
}


def post(client: Client, **fields: Any) -> Any:
    return client.post("/jobs/new/", {**BASE, **fields})


def box_of(page: str) -> str:
    tags: list[str] = [
        tag for tag in re.findall(r"<input[^>]*>", page) if 'name="properties"' in tag
    ]
    assert len(tags) == 1, tags
    return tags[0]


class TestThePipelineDefault:
    """The pipeline measures properties by default (owner, 2026-10-20), so a hand-run
    `nextflow run` matches a web job. Its default and `properties.DEFAULTS` are two
    lists of the same keys in different languages — exactly the drift that bit
    `gap_methods`, so they are pinned together here. The predictors are
    opt-in (owner, 2026-10-21), so they are registered but not in the default."""

    def default(self) -> list[str]:
        config = (ROOT / "pipeline" / "nextflow.config").read_text()
        line = next(line for line in config.splitlines() if line.strip().startswith("properties "))
        return [k for k in line.split("=", 1)[1].strip().strip("'\"").split(",") if k]

    def test_it_is_every_default_property(self) -> None:
        assert self.default() == list(properties.DEFAULTS)

    def test_main_nf_allows_exactly_those_keys(self) -> None:
        """main.nf repeats the keys, because they reach a command line. It lists them
        once (parsePropertyKeys), for the inline run and for PROPERTIES_ONLY alike."""
        main = (ROOT / "pipeline" / "main.nf").read_text()
        lines = [line for line in main.splitlines() if "'properties', [" in line]
        assert len(lines) == 1, "the keys are listed in one place"
        allowed = re.findall(r"'([a-z_]+)'", lines[0].split("'properties',", 1)[1])
        assert allowed == list(properties.REGISTRY)

    def test_main_nf_knows_which_need_phylobot_ml(self) -> None:
        """The PROPERTIES task picks its container from this list, so a predictor
        missing from it would run in phylobot-tools and fail to import torch."""
        main = (ROOT / "pipeline" / "main.nf").read_text()
        body = main.split("def mlPropertyKeys()", 1)[1].split("}", 1)[0]
        listed = re.findall(r"'([a-z_]+)'", body)
        assert listed == [k for k, p in properties.REGISTRY.items() if p.image == properties.ML]

    def test_the_ml_image_is_configured(self) -> None:
        config = (ROOT / "pipeline" / "nextflow.config").read_text()
        assert re.search(r"^\s*ml\s*=\s*'phylobot-ml:local'", config, re.M)
        module = (ROOT / "pipeline" / "modules" / "properties.nf").read_text()
        assert "params.containers.ml" in module


class TestTheBox:
    def test_it_is_on_the_page_and_ticked(self, client: Client) -> None:
        assert "checked" in box_of(client.get("/jobs/new/").content.decode())

    def test_it_sits_in_the_ancestors_step(self, client: Client) -> None:
        """So the step's switch hides it and its "Use defaults" resets it."""
        page = client.get("/jobs/new/").content.decode()
        step = page.split('data-needs="ancestors"')[1]
        assert 'name="properties"' in step.split("</section>")[0]

    def test_use_defaults_restores_the_tick(self) -> None:
        assert forms.section_defaults()["ancestors"]["properties"] is True


@pytest.mark.usefixtures("fake_storage")
class TestWhatItSends:
    def test_ticked_sends_every_default_property_key(self, client: Client) -> None:
        assert post(client, properties="on").status_code == 302
        job = Job.objects.latest("created_at")
        assert job.params["properties"] == list(properties.DEFAULTS)

    def test_a_ticked_predictor_is_sent_too(self, client: Client) -> None:
        assert post(client, properties="on", extra_properties="disorder").status_code == 302
        job = Job.objects.latest("created_at")
        assert job.params["properties"] == [*properties.DEFAULTS, "disorder"]

    def test_a_predictor_without_the_main_box_sends_nothing(self, client: Client) -> None:
        assert post(client, extra_properties="esm_likelihood").status_code == 302
        assert "properties" not in Job.objects.latest("created_at").params

    def test_the_predictors_start_unticked(self, client: Client) -> None:
        page = client.get("/jobs/new/").content.decode()
        for key in ("disorder", "esm_likelihood"):
            box = re.search(rf'<input[^>]*name="extra_properties"[^>]*value="{key}"[^>]*>', page)
            assert box and "checked" not in box.group(0)
        assert forms.section_defaults()["ancestors"]["extra_properties"] == []

    def test_unticked_sends_nothing(self, client: Client) -> None:
        """The pipeline's own default is none, so an absent param means no properties and
        a hand-run pipeline or golden scenario is unaffected."""
        assert post(client).status_code == 302
        assert "properties" not in Job.objects.latest("created_at").params

    def test_a_trees_only_job_sends_nothing(self, client: Client) -> None:
        assert post(client, properties="on", run_ancestors="").status_code == 302
        assert "properties" not in Job.objects.latest("created_at").params

    def test_the_worker_forwards_them_to_nextflow(self, client: Client, tmp_path: Path) -> None:
        """Storing the keys on the Job is not enough: `executors.PARAM_ALLOWLIST` drops
        anything not listed, and it drops it *silently* — the job still succeeds, just
        without properties. That is how this shipped broken, found on a real run."""
        assert post(client, properties="on").status_code == 302
        job = Job.objects.latest("created_at")
        sent = executors.build_params(job, tmp_path, tmp_path / "out")
        assert sent["properties"] == list(properties.DEFAULTS)

    def test_no_pipeline_param_is_dropped_on_the_way(self, client: Client, tmp_path: Path) -> None:
        """Every key the form puts on a Job either reaches Nextflow or is one of the few
        kept for display only. A new pipeline option that forgets the allowlist fails
        here rather than in a finished job with files missing."""
        display_only = {
            "sequence_count",  # My Jobs
            "simulated",  # badges; never a pipeline option
            "simulated_from",
            "left_out",  # the homolog review's record of what was not used
            "left_out_count",
        }
        assert post(client, properties="on").status_code == 302
        job = Job.objects.latest("created_at")
        sent = executors.build_params(job, tmp_path, tmp_path / "out")
        missing = set(job.params) - set(sent) - display_only
        assert not missing, f"dropped by PARAM_ALLOWLIST: {sorted(missing)}"

    def test_the_draw_count_is_the_users(self, client: Client, tmp_path: Path) -> None:
        """More draws narrow the interval and cost proportionally, so the number is the
        user's; 1000 stays the default (owner, 2026-10-21)."""
        assert post(client, properties="on", property_draws=2500).status_code == 302
        job = Job.objects.latest("created_at")
        assert job.params["property_draws"] == 2500
        sent = executors.build_params(job, tmp_path, tmp_path / "out")
        assert sent["property_draws"] == 2500

    def test_the_default_is_still_a_thousand(self, client: Client) -> None:
        assert post(client, properties="on").status_code == 302
        assert Job.objects.latest("created_at").params["property_draws"] == 1000
        assert forms.SubmitForm.base_fields["property_draws"].initial == 1000

    def test_out_of_range_is_refused_on_the_form(self, client: Client) -> None:
        """Unlike the viewer's field, this one can report: the page is right there."""
        assert post(client, properties="on", property_draws=10).status_code == 200
        assert post(client, properties="on", property_draws=10**9).status_code == 200
        assert not Job.objects.exists()

    def test_none_when_the_box_is_unticked(self, client: Client) -> None:
        assert post(client, property_draws=2500).status_code == 302
        assert "property_draws" not in Job.objects.latest("created_at").params

    def test_the_sentence_says_how_many(self, client: Client) -> None:
        assert "2500 draws each" in " ".join(
            forms.describe({"properties": ["charge"], "property_draws": 2500}, 0, 0)
        )

    def test_use_defaults_restores_the_count(self) -> None:
        assert forms.section_defaults()["ancestors"]["property_draws"] == 1000

    def test_duplicate_and_run_again_keeps_the_draw_count(self) -> None:
        asked = {"properties": ["charge"], "property_draws": 250}
        on, _ticks = forms.initial_from_job(asked, "", "")
        assert on["property_draws"] == 250
        off, _ = forms.initial_from_job({}, "", "")
        assert off["property_draws"] == 1000

    def test_duplicate_and_run_again_keeps_the_choice(self) -> None:
        on, _ticks = forms.initial_from_job({"properties": ["charge"]}, "", "")
        off, _ = forms.initial_from_job({}, "", "")
        assert on["properties"] is True
        assert off["properties"] is False


class TestDescribe:
    def test_the_sentence_mentions_properties_when_asked(self) -> None:
        steps = forms.describe({"aligners": ["mafft"], "properties": ["charge"]}, 0, 0)
        assert any("sequence properties" in step for step in steps)

    def test_and_not_otherwise(self) -> None:
        steps = forms.describe({"aligners": ["mafft"]}, 0, 0)
        assert not any("sequence properties" in step for step in steps)


class TestAccounting:
    """A new process has to appear in the steps graph and the work estimate together, or
    the job page's progress and its estimate disagree."""

    params: dict[str, Any] = {
        "aligners": ["mafft"],
        "tree_builders": ["iqtree"],
        "modelfinder": True,
        "properties": ["charge"],
    }

    def test_one_task_per_combination(self) -> None:
        kinds = stages.expected_kinds(self.params)
        assert kinds["PROPERTIES"] == kinds["EXPORT"]

    def test_none_without_the_param(self) -> None:
        without = {k: v for k, v in self.params.items() if k != "properties"}
        assert "PROPERTIES" not in stages.expected_kinds(without)

    def test_none_for_a_trees_only_job(self) -> None:
        assert "PROPERTIES" not in stages.expected_kinds({**self.params, "stop_after": "trees"})

    def test_it_runs_after_export(self) -> None:
        assert stages.DEPENDS["PROPERTIES"] == ["EXPORT"]

    def test_it_has_a_label_and_a_cost(self) -> None:
        assert stages.LABELS["PROPERTIES"]
        assert workload.STEPS["PROPERTIES"].seconds > 0

    def test_the_work_estimate_counts_the_same_tasks(self) -> None:
        counted = [kind for kind, _cost in workload.tasks(self.params)]
        assert counted.count("PROPERTIES") == stages.expected_kinds(self.params)["PROPERTIES"]


class TestJobPageLink:
    """The Ancestors cell links the properties index, and only when it was written."""

    key = "mafft/iqtree/modelfinder"
    summary = {
        "combinations": [
            {
                "key": key,
                "aligner": "mafft",
                "tree": "iqtree",
                "model": "modelfinder",
                "path": key,
            }
        ]
    }

    def links_with(self, available: set[str]) -> list[str]:
        [result] = results.alignment_results(
            self.summary, None, {}, lambda path: path, available=available
        )
        # The row's ZIP lists what it holds in its tooltip.
        return [cell.title for row in result.table.rows for cell in row if cell.title]

    def test_linked_when_written(self) -> None:
        available = {f"{self.key}/ancestors.fasta", f"{self.key}/properties/index.json"}
        assert any("properties" in title for title in self.links_with(available))

    def test_absent_otherwise(self) -> None:
        assert not any("properties" in t for t in self.links_with({f"{self.key}/ancestors.fasta"}))


class TestAfterTheFactJob:
    """A properties-only job: one PROPERTIES per combination and nothing else.
    Its plan reaches the job page as `properties_of` on the Job and the worker as a staged
    TSV, so both shapes have to count the same."""

    plan: dict[str, Any] = {
        "properties_only": True,
        "properties": ["charge"],
        "properties_of": {
            "parent": "00000000-0000-0000-0000-000000000000",
            "combinations": {"mafft/iqtree/LG+G4": ["Node1"], "muscle/iqtree/LG+G4": ["Node3"]},
        },
    }

    def test_the_steps_graph_is_one_kind(self) -> None:
        kinds = stages.expected_kinds(self.plan)
        assert dict(kinds) == {"PROPERTIES": 2}

    def test_the_task_count_agrees(self) -> None:
        assert executors.expected_tasks(self.plan) == sum(stages.expected_kinds(self.plan).values())

    def test_the_work_estimate_counts_the_same_tasks(self) -> None:
        assert [kind for kind, _cost in workload.tasks(self.plan)] == ["PROPERTIES"] * 2

    def test_the_worker_counts_from_the_staged_plan(self, tmp_path: Path) -> None:
        """The worker holds only the filtered pipeline params, where the plan is a file."""
        targets = tmp_path / executors.PROPERTY_TARGETS
        targets.write_text(
            "combination\tnode\n"
            "mafft/iqtree/LG+G4\tNode1\n"
            "mafft/iqtree/LG+G4\tNode4\n"
            "muscle/iqtree/LG+G4\tNode3\n"
        )
        params = {"properties_only": True, "property_targets": str(targets)}
        assert executors.expected_tasks(params) == 2  # two combinations, three ancestors

    def test_an_unreadable_plan_counts_as_one(self) -> None:
        """Never zero: progress would read "done 1 of 0" and the page would show -100%."""
        params = {"properties_only": True, "property_targets": "/nowhere/targets.tsv"}
        assert executors.expected_tasks(params) == 1

    def test_the_staged_inputs_become_pipeline_params(self, tmp_path: Path) -> None:
        (tmp_path / executors.PROPERTY_DIR / "mafft" / "iqtree" / "LG+G4").mkdir(parents=True)
        (tmp_path / executors.PROPERTY_TARGETS).write_text("combination\tnode\n")
        sent = executors.staged_inputs(tmp_path)
        assert sent["properties_dir"] == str(tmp_path / executors.PROPERTY_DIR)
        assert sent["property_targets"] == str(tmp_path / executors.PROPERTY_TARGETS)

    def test_the_flag_reaches_nextflow(self) -> None:
        assert "properties_only" in executors.PARAM_ALLOWLIST
        # The record of what was asked for is the Job's, not the pipeline's: the
        # pipeline is handed the staged TSV instead.
        assert "properties_of" not in executors.PARAM_ALLOWLIST

    def test_main_nf_has_the_entry_workflow(self) -> None:
        main = (ROOT / "pipeline" / "main.nf").read_text()
        assert "workflow PROPERTIES_ONLY" in main
        assert "params.properties_only" in main
        config = (ROOT / "pipeline" / "nextflow.config").read_text()
        for param in ("properties_only", "properties_dir", "property_targets"):
            assert re.search(rf"^\s*{param}\s*=", config, re.M), param
