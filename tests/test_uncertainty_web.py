"""The Uncertainty page and its follow-up job (phylobot.uncertainty_submit,
phylobot.uncertainty_views): what the button queues and stages,
that asking twice queues once, that only the job's own values are accepted, and the page
reading a run's files back."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from django.test import Client

from phylobot import executors, job_types, jobs, usage
from phylobot import uncertainty_run as r
from phylobot.executors import UNCERTAINTY_DIR, UNCERTAINTY_SETTINGS
from phylobot.models import Job
from tests.conftest import FakeStorage
from tests.test_views import viewer_job

pytestmark = pytest.mark.django_db


def queued() -> Job:
    return Job.objects.filter(params__uncertainty_only=True).latest("created_at")


def compute(client: Client, job: Job, **fields: Any) -> Any:
    return client.post(f"/jobs/{job.id}/uncertainty/run/", {"weighting": "flat", **fields})


def finish(run: Job, parent: Job, fake_storage: FakeStorage) -> None:
    """Run the pipeline step on the parent's files, as the worker would, and store them."""
    prefix = f"jobs/{parent.id}/results/"

    def read(relative: str) -> str | None:
        data = fake_storage.objects.get(prefix + relative)
        return data.decode() if data is not None else None

    settings = r.Settings.parse(run.params["uncertainty_of"]["settings"])
    for name, text in r.run(read, settings).items():
        fake_storage.objects[f"jobs/{run.id}/results/uncertainty/{name}"] = text.encode()
    assert jobs.mark(run.id, Job.Status.RUNNING) and jobs.mark(run.id, Job.Status.SUCCEEDED)


class TestThePageBeforeARun:
    def test_it_offers_the_settings(self, client: Client, fake_storage: FakeStorage) -> None:
        job = viewer_job(fake_storage, gap_methods=["fitch"])
        page = client.get(f"/jobs/{job.id}/uncertainty/").content.decode()
        assert "No uncertainty run on this job yet" in page
        assert f'action="/jobs/{job.id}/uncertainty/run/"' in page
        assert 'name="weighting" value="bic"' in page
        assert 'name="weight-alignment-1"' in page  # MUSCLE, by position, never by name
        assert 'value="fitch" checked' in page  # the job's mask, ticked

    def test_the_job_page_links_to_it(self, client: Client, fake_storage: FakeStorage) -> None:
        job = viewer_job(fake_storage)
        page = client.get(f"/jobs/{job.id}/").content.decode()
        assert f'href="/jobs/{job.id}/uncertainty/"' in page


class TestTheFollowUpJob:
    def test_it_queues_its_own_job_and_stages_the_parents_files(
        self, client: Client, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        response = compute(client, parent)
        run = queued()
        assert response.status_code == 302
        assert response["Location"] == f"/jobs/{parent.id}/uncertainty/?run={run.id}"
        assert run.id != parent.id and run.status == Job.Status.QUEUED
        assert run.params["uncertainty_only"] is True
        of = run.params["uncertainty_of"]
        assert of["parent"] == str(parent.id) and of["combinations"] == 3
        staged = f"jobs/{run.id}/inputs/{UNCERTAINTY_DIR}/"
        names = {k[len(staged) :] for k in fake_storage.objects if k.startswith(staged)}
        assert "summary.json" in names
        assert "mafft/iqtree/WAG/ancestors_posteriors.tsv" in names
        assert "muscle/iqtree/modelfinder/ancestors.gaps-fitch.fasta" in names
        assert "mafft/iqtree/WAG/supports.tsv" not in names  # only what UNCERTAINTY reads
        settings = json.loads(fake_storage.objects[f"jobs/{run.id}/inputs/{UNCERTAINTY_SETTINGS}"])
        assert r.Settings.parse(settings) == r.Settings()

    def test_asking_twice_queues_once(self, client: Client, fake_storage: FakeStorage) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        compute(client, parent)
        compute(client, parent)
        assert Job.objects.filter(params__uncertainty_only=True).count() == 1
        compute(client, parent, weighting="bic")
        assert Job.objects.filter(params__uncertainty_only=True).count() == 2

    def test_multipliers_name_levels_by_position(
        self, client: Client, fake_storage: FakeStorage
    ) -> None:
        parent = viewer_job(fake_storage)
        compute(client, parent, **{"weight-alignment-1": "0", "weight-model-0": "2"})
        settings = queued().params["uncertainty_of"]["settings"]
        assert settings["multipliers"] == {
            "alignment": {"muscle": 0.0},
            "model": {"modelfinder": 2.0},
        }

    @pytest.mark.parametrize(
        "fields",
        [
            {"weighting": "aic"},
            {"reference": "x/y/z"},
            {"versions": "grasp"},  # the job has no GRASP reconstruction
            {"weight-tree-0": "many"},
            {"weight-tree-0": "-1"},
        ],
    )
    def test_only_the_jobs_own_values(
        self, client: Client, fake_storage: FakeStorage, fields: dict[str, str]
    ) -> None:
        parent = viewer_job(fake_storage)
        assert compute(client, parent, **fields).status_code == 404
        assert not Job.objects.filter(params__uncertainty_only=True).exists()

    def test_reaches_the_pipeline(self, tmp_path: Path) -> None:
        assert "uncertainty_only" in executors.PARAM_ALLOWLIST
        assert "uncertainty_of" not in executors.PARAM_ALLOWLIST  # display only
        (tmp_path / UNCERTAINTY_DIR).mkdir()
        (tmp_path / UNCERTAINTY_SETTINGS).write_text("{}")
        sent = executors.staged_inputs(tmp_path)
        assert sent["uncertainty_dir"] == str(tmp_path / UNCERTAINTY_DIR)
        assert sent["uncertainty_settings"] == str(tmp_path / UNCERTAINTY_SETTINGS)
        main = (Path(__file__).parents[1] / "pipeline" / "main.nf").read_text()
        assert "workflow UNCERTAINTY_ONLY" in main and "params.uncertainty_only" in main
        config = (Path(__file__).parents[1] / "pipeline" / "nextflow.config").read_text()
        for param in ("uncertainty_only", "uncertainty_dir", "uncertainty_settings"):
            assert f"    {param} =" in config

    def test_it_is_a_kind_of_its_own(self) -> None:
        params = {"uncertainty_only": True}
        assert [t.key for t in job_types.types_of(params)] == ["uncertainty"]
        assert usage.kind_of(params) == "uncertainty"
        assert executors.expected_tasks(params) == 1


class TestThePageWithARun:
    def test_it_reads_the_run(self, client: Client, fake_storage: FakeStorage) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        compute(client, parent)
        run = queued()
        finish(run, parent, fake_storage)
        page = client.get(f"/jobs/{parent.id}/uncertainty/").content.decode()
        assert (
            "On average over the ancestors, a sample from each of two reconstructions picked"
            in page
        )
        assert "How to Read This Page" in page and "3 positions differ" in page
        assert "A worked example: two columns of one ancestor" in page
        assert "Each Ancestor" in page and ">Node1</text>" in page
        assert "Map of the Reconstructions" in page and 'name="colour_by"' in page
        assert "Show every ancestor on the tree" in page and 'class="pie"' in page
        assert 'href="?node=Node1#node"' in page or "node=Node1" in page
        assert 'data-tip="Node1 · Posterior, residue identity' in page
        assert 'fill="url(#nd-hatch-' in page  # presence is hatched
        assert "Not varied in this job" in page and "tree" in page

    def test_a_node_site_by_site(self, client: Client, fake_storage: FakeStorage) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        compute(client, parent)
        finish(queued(), parent, fake_storage)
        page = client.get(f"/jobs/{parent.id}/uncertainty/", {"node": "Node1"}).content.decode()
        assert "Node1, Site by Site" in page
        assert "How far apart the reconstructions are" in page
        page = client.get(
            f"/jobs/{parent.id}/uncertainty/", {"node": "Node1", "mode": "map", "scale": "share"}
        ).content.decode()
        assert '<option value="map" selected>' in page and "100%" in page

    def test_a_run_still_queued(self, client: Client, fake_storage: FakeStorage) -> None:
        parent = viewer_job(fake_storage)
        compute(client, parent)
        page = client.get(f"/jobs/{parent.id}/uncertainty/").content.decode()
        assert "This run is" in page and "queued" in page

    def test_the_runs_job_page_links_back(self, client: Client, fake_storage: FakeStorage) -> None:
        parent = viewer_job(fake_storage)
        compute(client, parent)
        run = queued()
        page = client.get(f"/jobs/{run.id}/").content.decode()
        assert f"/jobs/{parent.id}/uncertainty/?run={run.id}" in page

    def test_not_on_a_run_itself(self, client: Client, fake_storage: FakeStorage) -> None:
        parent = viewer_job(fake_storage)
        compute(client, parent)
        run = queued()
        finish(run, parent, fake_storage)
        assert client.get(f"/jobs/{run.id}/uncertainty/").status_code == 404


class TestDownloads:
    def page_with_run(self, client: Client, fake_storage: FakeStorage) -> tuple[Job, Job]:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        compute(client, parent)
        run = queued()
        finish(run, parent, fake_storage)
        return parent, run

    def test_every_site_as_csv(self, client: Client, fake_storage: FakeStorage) -> None:
        parent, run = self.page_with_run(client, fake_storage)
        page = client.get(f"/jobs/{parent.id}/uncertainty/", {"node": "Node1"}).content.decode()
        link = f"/jobs/{parent.id}/uncertainty/download/?run={run.id}&amp;mode=posterior"
        assert f"{link}&amp;table=sites&amp;format=csv&amp;node=Node1" in page
        response = client.get(
            f"/jobs/{parent.id}/uncertainty/download/",
            {"run": str(run.id), "table": "sites", "node": "Node1", "format": "csv"},
        )
        assert response["Content-Disposition"] == (
            'attachment; filename="uncertainty.posterior.Node1.sites.csv"'
        )
        header, *rows = response.content.decode().splitlines()
        assert header.startswith("Column,Alignment agreement,Variance,Variance: Posterior")
        assert "Reconstruction, Alignment: MAFFT" in header
        assert "Reconstruction, Gap placement: Parsimony (Fitch)" in header
        assert [r.split(",")[0] for r in rows] == ["1", "2", "3"]  # every column, in order

    def test_ancestors_as_excel(self, client: Client, fake_storage: FakeStorage) -> None:
        from io import BytesIO

        from openpyxl import load_workbook

        parent, run = self.page_with_run(client, fake_storage)
        page = client.get(f"/jobs/{parent.id}/uncertainty/").content.decode()
        assert "table=ancestors&amp;format=xlsx" in page
        response = client.get(
            f"/jobs/{parent.id}/uncertainty/download/",
            {"run": str(run.id), "table": "ancestors", "format": "xlsx", "mode": "map"},
        )
        book = load_workbook(BytesIO(response.content))
        assert book.sheetnames == ["Ancestors", "About"]
        rows = list(book["Ancestors"].iter_rows(values_only=True))
        assert rows[0][:5] == (
            "Ancestor",
            "Present in",
            "Of combinations",
            "Sequences below",
            "Total",
        )
        assert {r[0] for r in rows[1:]} == {"Node1", "Node2"}
        assert isinstance(rows[1][4], float)  # numbers stay numbers
        about = {str(k): v for k, v in book["About"].iter_rows(values_only=True)}
        assert about["Compared"] == "Most probable sequences"

    @pytest.mark.parametrize(
        "fields",
        [{"table": "nodes"}, {"format": "pdf"}, {"table": "sites", "node": "Nope"}],
    )
    def test_refused(
        self, client: Client, fake_storage: FakeStorage, fields: dict[str, str]
    ) -> None:
        parent, run = self.page_with_run(client, fake_storage)
        query = {"run": str(run.id), "table": "ancestors", "format": "csv"} | fields
        assert client.get(f"/jobs/{parent.id}/uncertainty/download/", query).status_code == 404


class TestHeatmapColours:
    def test_a_palette_and_its_key(self, client: Client, fake_storage: FakeStorage) -> None:
        parent = viewer_job(fake_storage, gap_methods=["fitch"])
        compute(client, parent)
        finish(queued(), parent, fake_storage)
        url = f"/jobs/{parent.id}/uncertainty/"
        page = client.get(url, {"node": "Node1", "heat": "viridis"}).content.decode()
        assert '<option value="viridis" selected>' in page
        assert 'name="heat" form="uncertainty"' in page
        assert "Expected differing positions" in page and 'fill="url(#heat-key)"' in page
        assert 'stop-color="#fde725"' in page  # viridis's near end
        page = client.get(url, {"node": "Node1", "heat": "rainbow"}).content.decode()
        assert '<option value="blues" selected>' in page  # unknown: the default


def test_the_track_keeps_its_text_size_on_any_frame() -> None:
    """The figure scales to the card, so a short frame must not be drawn small."""
    from phylobot import uncertainty_views as v

    def rows(n: int) -> list[dict[str, str]]:
        values = ("total", "within", "alignment", "tree", "model", "gaps", "other")
        parts = [f"{p}_{k}" for p in ("presence", "residue") for k in values]
        names = ["mode", "column", "agreement", *parts]
        return [
            dict(zip(names, ["posterior", str(c), "1"] + ["0.1"] * 14, strict=True))
            for c in range(1, n + 1)
        ]

    for n in (37, 181, 956, 3000):
        drawn = v.track(rows(n), "posterior")
        assert drawn.width == v.TRACK_WIDTH
        assert 3 <= len(drawn.ticks) <= 10
    assert v.column_step(181) == 20 and v.column_step(956) == 100 and v.column_step(37) == 5
