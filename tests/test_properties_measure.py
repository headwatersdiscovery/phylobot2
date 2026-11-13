"""Measuring properties after a job has finished (`phylobot.properties_submit`):
the viewer's button, the follow-up job it queues, what it stages, and the panel reading
that run's results back.
"""

from __future__ import annotations

from typing import Any

import pytest

from phylobot import ancestors, properties, properties_run, property_views
from phylobot.executors import PROPERTY_DIR, PROPERTY_TARGETS
from phylobot.models import Job
from tests.conftest import FakeStorage
from tests.test_views import store_properties, viewer_job

pytestmark = pytest.mark.django_db

MEASURE: dict[str, Any] = {
    "combination": "mafft/iqtree/modelfinder",
    "node": "Node1",
    "gaps": "fitch",
    properties_run.FIELD: ["charge", "moment"],
    properties_run.SCOPE_FIELD: "node",
}


def measure(client: Any, job: Job, **changed: Any) -> Any:
    return client.post(f"/jobs/{job.id}/properties/measure/", {**MEASURE, **changed})


def queued() -> Job:
    return Job.objects.filter(params__properties_only=True).latest("created_at")


class TestTheButtonIsOnThePage:
    def test_it_posts_the_ticked_properties(self, client: Any, fake_storage: FakeStorage) -> None:
        job = viewer_job(fake_storage, gap_methods=["fitch"])
        page = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
        assert f'action="/jobs/{job.id}/properties/measure/"' in page
        assert 'form="properties-measure"' in page and "Launch" in page
        # The viewer is a GET form, so the button's values are tied to the page-level
        # POST form by `form=` and never reach the viewer's URL.
        assert 'value="charge" form="properties-measure"' in page

    def test_no_button_with_nothing_ticked(self, client: Any, fake_storage: FakeStorage) -> None:
        job = viewer_job(fake_storage, gap_methods=["fitch"])
        url = f"/jobs/{job.id}/ancestors/"
        page = client.get(url, {properties_run.FIELD: ""}).content.decode()
        # The dialog's Launch is always in the markup (script opens it per property);
        # what must be absent with nothing ticked is the plain form's run.
        assert "offer-actions" not in page


class TestTheFollowUpJob:
    def test_it_is_its_own_job_with_its_own_plan(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        response = measure(client, parent)
        assert response.status_code == 302
        job = queued()
        assert job.id != parent.id and job.status == Job.Status.QUEUED
        of = job.params["properties_of"]
        assert of["parent"] == str(parent.id) and of["gaps"] == "fitch"
        assert of["combinations"] == {"mafft/iqtree/modelfinder": ["Node1"]}
        assert job.params["properties"] == ["charge", "moment"]
        assert job.params["properties_only"] is True
        # The masked version the viewer was showing, so the run measures that one.
        assert job.params["gap_methods"] == ["fitch"]
        assert response["Location"] == f"/jobs/{job.id}/"

    def test_its_results_are_its_own(self, client: Any, fake_storage: FakeStorage) -> None:
        """§18.3: nothing is written into the parent, so deleting either job cannot
        strip the other."""
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        before = set(fake_storage.list_keys(f"jobs/{parent.id}/"))
        measure(client, parent)
        assert set(fake_storage.list_keys(f"jobs/{parent.id}/")) == before

    def test_it_stages_the_ancestors_and_the_plan(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent)
        job = queued()
        staged = {
            key[len(f"jobs/{job.id}/inputs/") :]
            for key in fake_storage.list_keys(f"jobs/{job.id}/inputs/")
        }
        combination = MEASURE["combination"]
        assert f"{PROPERTY_DIR}/{combination}/ancestors.gaps-fitch.fasta" in staged
        assert f"{PROPERTY_DIR}/{combination}/ancestors_posteriors.tsv" in staged
        # Not the unmasked reconstruction: the run measures the version asked for.
        assert f"{PROPERTY_DIR}/{combination}/ancestors.fasta" not in staged
        plan = fake_storage.read_bytes(f"jobs/{job.id}/inputs/{PROPERTY_TARGETS}").decode()
        assert plan == f"combination\tnode\n{combination}\tNode1\n"

    def test_with_an_index_only_the_measured_nodes_rows_are_staged(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        """Read by ranges through the parent's byte index: one ancestor's rows, not the
        whole table."""
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        combination = MEASURE["combination"]
        table = f"jobs/{parent.id}/results/{combination}/ancestors_posteriors.tsv"
        text = fake_storage.objects[table].decode()
        header = len(text.splitlines(keepends=True)[0])
        runs: dict[str, list[tuple[int, int]]] = {}
        offset = header
        for line in text.splitlines(keepends=True)[1:]:
            runs.setdefault(line.split("\t")[0], []).append((offset, len(line)))
            offset += len(line)
        index = ancestors.ByteIndex(header, runs)
        fake_storage.objects[table.replace(".tsv", ".index.tsv")] = ancestors.write_index(
            index
        ).encode()
        measure(client, parent)
        job = queued()
        staged = fake_storage.read_bytes(
            f"jobs/{job.id}/inputs/{PROPERTY_DIR}/{combination}/ancestors_posteriors.tsv"
        ).decode()
        lines = staged.splitlines()
        assert lines[0] == text.splitlines()[0] and len(lines) == 4  # header, Node1's 3 sites
        assert {line.split("\t")[0] for line in lines[1:]} == {"Node1"}

    def test_the_scope_names_each_combination_its_own_node(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        """MUSCLE numbers the same clade Node2, so the plan must not repeat Node1."""
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent, **{properties_run.SCOPE_FIELD: "analogs"})
        of = queued().params["properties_of"]
        assert of["combinations"] == {
            "mafft/iqtree/WAG": ["Node1"],
            "mafft/iqtree/modelfinder": ["Node1"],
            "muscle/iqtree/modelfinder": ["Node2"],
        }

    def test_asking_twice_leads_to_the_one_run(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        first = measure(client, parent)
        again = measure(client, parent)
        assert again["Location"] == first["Location"]
        assert Job.objects.filter(params__properties_only=True).count() == 1

    def test_a_different_plan_is_a_different_run(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent)
        measure(client, parent, **{properties_run.FIELD: ["charge"]})
        assert Job.objects.filter(params__properties_only=True).count() == 2

    def test_a_failed_run_can_be_asked_for_again(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent)
        Job.objects.filter(params__properties_only=True).update(status=Job.Status.FAILED)
        measure(client, parent)
        assert Job.objects.filter(params__properties_only=True).count() == 2


class TestNothingFromTheRequestIsTrusted:
    """Every value is a choice among the parent job's own: the viewer is public input."""

    def setup_parent(self, fake_storage: FakeStorage) -> Job:
        return viewer_job(fake_storage, gap_methods=["fitch"])

    @pytest.mark.parametrize(
        "changed",
        [
            {"combination": "mafft/iqtree/../../etc"},
            {"combination": "muscle/iqtree/WAG"},  # not a combination of this job
            {"node": "Node9"},
            {"node": "../../etc/passwd"},
            {"gaps": "prank"},  # no posteriors to draw from
            {"gaps": "../.."},
            {properties_run.FIELD: ["solubility"]},  # not in the registry
            {properties_run.FIELD: []},
        ],
    )
    def test_refused(self, client: Any, fake_storage: FakeStorage, changed: dict[str, Any]) -> None:
        parent = self.setup_parent(fake_storage)
        assert measure(client, parent, **changed).status_code == 404
        assert not Job.objects.filter(params__properties_only=True).exists()

    def test_another_users_job(self, client: Any, fake_storage: FakeStorage) -> None:
        from tests.test_views import other_users_job

        assert measure(client, other_users_job()).status_code == 404

    def test_an_unfinished_job(self, client: Any, fake_storage: FakeStorage) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        Job.objects.filter(id=parent.id).update(status=Job.Status.RUNNING)
        assert measure(client, parent).status_code == 404

    def test_an_unknown_scope_falls_back_rather_than_failing(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        assert measure(client, parent, **{properties_run.SCOPE_FIELD: "everything!"}).status_code
        assert queued().params["properties_of"]["scope"] == properties_run.DEFAULT_SCOPE


def finish_run(
    fake_storage: FakeStorage,
    run: Job,
    combination: str,
    node: str,
    mean: float = 0.31,
    draws: int = 1000,
) -> None:
    """A finished properties run's own results: one property, under its own prefix."""
    Job.objects.filter(id=run.id).update(status=Job.Status.SUCCEEDED)
    body = {
        "node": node, "gaps": "fitch", "draws": draws, "seed": "x",
        "sampling": "per-column marginal posteriors, drawn independently", "weighted": False,
        "properties": {
            "moment": {"label": "Hydrophobic moment", "units": "", "level": 0.95,
                       "mean": mean, "sd": 0.01, "low": 0.29, "high": 0.33, "ml": 0.3,
                       "standardised_shift": -1.0},
        },
    }  # fmt: skip
    fake_storage.store_properties(
        f"jobs/{run.id}/results/", combination, node, "fitch", body, {"moment": [0.3] * 20}
    )


class TestTheViewerReadsTheRunBack:
    """§18.3: the run owns its results, so the panel is assembled from both prefixes."""

    def test_a_later_measurement_joins_the_panel(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        store_properties(fake_storage, parent)  # the job's own: charge
        measure(client, parent, **{properties_run.FIELD: ["moment"]})
        run = queued()
        finish_run(fake_storage, run, MEASURE["combination"], "Node1")

        page = client.get(f"/jobs/{parent.id}/ancestors/").content.decode()
        assert "Net charge" in page  # the job's own measurement, still there
        assert "Hydrophobic moment" in page  # the later run's
        assert "measured later" in page  # and said so
        # Both now count as measured on the checklist; only the third is still offered.
        assert page.count("already measured") == 2

    def test_a_queued_run_is_visible_before_it_finishes(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent)
        page = client.get(f"/jobs/{parent.id}/ancestors/").content.decode()
        assert "A properties run is queued" in page

    def test_the_comparison_page_reads_them_too(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        """§18.6 step 6: the scope "every combination that recovered it" exists to feed
        the comparison page, so that page has to read the run as well."""
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent, **{properties_run.SCOPE_FIELD: "analogs"})
        run = queued()
        for combination, node in run.params["properties_of"]["combinations"].items():
            finish_run(fake_storage, run, combination, node[0])
        page = client.get(f"/jobs/{parent.id}/properties/", {"node": "Node1"}).content.decode()
        assert "Hydrophobic moment" in page
        # Every combination that recovered the ancestor, from a job that measured none.
        for label in ("MAFFT · ModelFinder", "MAFFT · WAG", "MUSCLE · ModelFinder"):
            assert label in page, label

    def test_the_run_of_another_node_is_not_read_here(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent, node="Node2")
        run = queued()
        finish_run(fake_storage, run, MEASURE["combination"], "Node2")
        page = client.get(f"/jobs/{parent.id}/ancestors/", {"node": "Node1"}).content.decode()
        assert "Hydrophobic moment (" not in page
        assert "Measure or recompute properties" in page  # Node1 still has nothing


class TestItIsItsOwnKindOfJob:
    """My Jobs must not file a properties run under Alignment / Trees / Ancestors: it
    runs none of those, and its sentence should say what it does."""

    def test_the_type_and_its_filter(self, client: Any, fake_storage: FakeStorage) -> None:
        from phylobot import job_types
        from phylobot.models import Job as JobModel

        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent)
        run = queued()
        assert [t.key for t in job_types.types_of(run.params)] == ["properties"]
        assert job_types.BY_KEY["properties"].label == "Ancestral properties"
        assert [t.key for t in job_types.types_of(parent.params)] != ["properties"]
        listed = JobModel.objects.filter(job_types.type_filter("properties"))
        assert list(listed.values_list("id", flat=True)) == [run.id]
        # ...and it is not swept up by the analysis filters either.
        for key in ("alignments", "trees", "ancestors"):
            assert run.id not in set(
                JobModel.objects.filter(job_types.type_filter(key)).values_list("id", flat=True)
            ), key

    def test_the_sentence_says_what_it_measures(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        from phylobot.forms import describe

        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent)
        [sentence] = describe(queued().params, 0, 0)
        assert "Net charge" in sentence and "1 ancestor" in sentence and "1000 draws" in sentence
        page = client.get(f"/jobs/{queued().id}/").content.decode()
        assert "Net charge" in page
        # and the way back to the ancestor it measured
        assert f"/jobs/{parent.id}/ancestors/?combination=" in page
        assert "open it in the ancestor viewer" in page
        assert "alignments only" not in page  # it aligns nothing


class TestDrawingAgain:
    """ "Recompute": the same plan, a new seed. A property posterior is a sample, so the
    numbers move; what must not move is which run the viewer shows."""

    def test_the_control_offers_the_measured_properties(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        job = viewer_job(fake_storage, gap_methods=["fitch"])
        # Nothing measured: the list offers each property, and none of them a recompute.
        bare = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
        assert "Recompute" not in bare
        store_properties(fake_storage, job)
        page = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
        # Measuring and recomputing are one list and one dialog (owner, 2026-10-21):
        # a measured property's row opens it in recompute mode.
        assert '<dialog id="measure-dialog"' in page
        assert 'data-offer-key="charge" data-offer-label="Net charge (pH 7)"' in page
        assert 'data-offer-redraw="on"' in page  # charge is measured
        assert 'data-offer-key="moment"' in page and 'data-offer-redraw=""' in page
        assert "Recompute" in page and "Launch" in page
        # The dialog is useless without its script, so the two ship together.
        assert "property_measure.js" in page

    def test_the_dialog_is_not_shown_without_its_script(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        """Pins the markup to the script that opens it: a <dialog> no one can open is a
        dead end, and this is how that would be caught."""
        from pathlib import Path

        template = Path("src/phylobot/templates/phylobot/ancestors.html").read_text()
        assert "property_measure.js" in template

    def test_it_queues_a_new_run_with_the_next_seed(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent, **{properties_run.FIELD: ["charge"]})
        first = queued()
        assert first.params["seed"] == 1 and first.params["properties_of"]["attempt"] == 0
        finish_run(fake_storage, first, MEASURE["combination"], "Node1")

        measure(client, parent, **{properties_run.FIELD: ["charge"]}, redraw="on")
        again = queued()
        assert again.id != first.id
        assert again.params["properties_of"]["attempt"] == 1
        # A seed past every measurement this ancestor has: PROPERTIES seeds per node, so
        # the parent's seed would have reproduced the first run's numbers exactly.
        assert again.params["seed"] == 2
        # The same work, told apart by the attempt: the base is what dedup compares.
        assert again.params["properties_of"]["base"] == first.params["properties_of"]["base"]

    def test_a_narrower_redraw_still_moves_the_seed(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        """The trap: recomputing one ancestor after a four-ancestor run shares no plan
        with it, so an attempt count alone would hand back the same draw under a new
        name. The seed counts this ancestor's measurements, not the plan's attempts."""
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent, **{properties_run.SCOPE_FIELD: "analogs"})
        wide = queued()
        for combination, nodes in wide.params["properties_of"]["combinations"].items():
            finish_run(fake_storage, wide, combination, nodes[0])
        measure(client, parent, redraw="on")  # this ancestor only
        assert queued().params["seed"] == wide.params["seed"] + 1

    def test_a_second_click_while_it_runs_leads_to_the_run(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        first = measure(client, parent, **{properties_run.FIELD: ["charge"]})
        again = measure(client, parent, **{properties_run.FIELD: ["charge"]}, redraw="on")
        assert again["Location"] == first["Location"]
        assert Job.objects.filter(params__properties_only=True).count() == 1

    def test_without_the_flag_it_is_still_one_run(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent)
        measure(client, parent)
        assert Job.objects.filter(params__properties_only=True).count() == 1

    def test_the_viewer_shows_the_newer_run(self, client: Any, fake_storage: FakeStorage) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent, **{properties_run.FIELD: ["moment"]})
        finish_run(fake_storage, queued(), MEASURE["combination"], "Node1")
        measure(client, parent, **{properties_run.FIELD: ["moment"]}, redraw="on")
        newer = queued()
        finish_run(fake_storage, newer, MEASURE["combination"], "Node1", mean=0.99)
        page = client.get(f"/jobs/{parent.id}/ancestors/").content.decode()
        assert "0.99" in page and "<strong>0.31</strong>" not in page


class TestHowManyDraws:
    def test_the_viewer_asks_for_the_number_shown(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        page = client.get(f"/jobs/{parent.id}/ancestors/", {"draws": "2500"}).content.decode()
        assert 'name="draws" value="2500"' in page
        measure(client, parent, draws="2500")
        job = queued()
        assert job.params["property_draws"] == 2500
        assert job.params["properties_of"]["draws"] == 2500

    def test_any_whole_number_in_range_is_allowed(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        """step="100" on the viewer's fields made the browser refuse 1500 or 250 — the
        bound is a range, not a grid (owner, 2026-10-21)."""
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        store_properties(fake_storage, parent)
        page = client.get(f"/jobs/{parent.id}/ancestors/", {"node": "Node2"}).content.decode()
        assert 'step="100"' not in page
        measure(client, parent, draws="1234")
        assert queued().params["property_draws"] == 1234

    def test_a_silly_number_is_clamped_not_refused(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        """The viewer's own number field: a 404 for a typo would be unkind, and the
        bounds are what keeps a tick box from queueing an hour of work."""
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        measure(client, parent, draws="99999999")
        assert queued().params["property_draws"] == properties.MAX_DRAWS
        measure(client, parent, draws="three", **{properties_run.FIELD: ["charge"]})
        assert queued().params["property_draws"] == properties.DEFAULT_DRAWS

    def test_the_cost_follows_the_number(self, client: Any, fake_storage: FakeStorage) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        url = f"/jobs/{parent.id}/ancestors/"
        few = client.get(url, {"draws": "100"}).content.decode()
        many = client.get(url, {"draws": "10000"}).content.decode()
        assert "at 100 draws" in few and "at 10000 draws" in many

    def test_the_caption_says_what_was_actually_used(
        self, client: Any, fake_storage: FakeStorage
    ) -> None:
        """The caption claimed 1000 whatever the run did; with runs at different counts
        it now says so per property."""
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        store_properties(fake_storage, parent)  # charge, at 1000
        measure(client, parent, **{properties_run.FIELD: ["moment"]}, draws="2500")
        finish_run(fake_storage, queued(), MEASURE["combination"], "Node1", draws=2500)
        page = client.get(f"/jobs/{parent.id}/ancestors/").content.decode()
        assert "the number of sequences shown with it" in page
        assert "1000 draws" in page and "2500 draws" in page


class TestMergingSummaries:
    """Pure: which number wins when an ancestor was measured more than once."""

    def summary(self, **entries: Any) -> dict[str, Any]:
        return {"draws": 1000, "properties": {k: {"mean": v} for k, v in entries.items()}}

    def test_later_wins_per_property(self) -> None:
        merged = property_views.merge([
            self.summary(charge=1.0, moment=2.0), self.summary(moment=9.0)
        ])  # fmt: skip
        assert merged is not None
        assert merged["properties"]["charge"]["mean"] == 1.0
        assert merged["properties"]["moment"]["mean"] == 9.0

    def test_nothing_measured_is_none(self) -> None:
        assert property_views.merge([None, {}, {"properties": {}}]) is None

    def test_one_summary_is_itself(self) -> None:
        one = self.summary(charge=1.0)
        merged = property_views.merge([None, one])
        assert merged is not None
        assert merged["properties"]["charge"]["mean"] == 1.0

    def test_each_property_keeps_its_own_draw_count(self) -> None:
        """Runs can be asked for different numbers, so the count belongs to the
        property, not to the panel (it is what the caption claims)."""
        merged = property_views.merge([
            {"draws": 1000, "properties": {"charge": {"mean": 1.0}}},
            {"draws": 5000, "properties": {"moment": {"mean": 2.0}}},
        ])  # fmt: skip
        assert merged is not None
        assert merged["properties"]["charge"]["draws"] == 1000
        assert merged["properties"]["moment"]["draws"] == 5000


class TestWhatIsOfferedAfterwards:
    def test_the_registry_minus_what_every_run_measured(self) -> None:
        both = property_views.merge([
            {"properties": {"charge": {"mean": 1}}}, {"properties": {"moment": {"mean": 2}}}
        ])  # fmt: skip
        assert properties_run.missing(both) == [
            k for k in properties.REGISTRY if k not in ("charge", "moment")
        ]

    def test_only_versions_with_posteriors_can_be_measured(self) -> None:
        from phylobot import properties_submit

        assert properties_submit.measurable(ancestors.ORIGINAL)
        assert properties_submit.measurable("fitch")
        assert not properties_submit.measurable("prank")
        assert not properties_submit.measurable("grasp")
        # codeml's reconstructions have posteriors, but under their own folder, which
        # PROPERTIES_ONLY does not stage: not offered rather than offered and refused.
        assert not properties_submit.measurable("codeml:M0")
