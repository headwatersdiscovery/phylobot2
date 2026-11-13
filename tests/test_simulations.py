"""Users' simulations (phylobot.simulations, simulation_views): the Simulations page's
form, the jobs it queues, the replicates' registration as simulated alignments, the job
page, and "use in a new job" marking jobs on simulated data."""

import json
import re
from pathlib import Path
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.test import Client

from phylobot import alignments, filelist, jobs, simulations
from phylobot.executors import build_params, expected_tasks
from phylobot.models import Alignment, Job
from phylobot.modes import local_user
from phylobot.simulate import parse
from tests.conftest import FakeStorage

ROOTED = "((a:0.1,b:0.2)Node2:0.15,(c:0.3,d:0.4)Node3:0.15)Node1;"
COMBO = "mafft/iqtree/modelfinder"


def finished_job(fake_storage: FakeStorage, owner: User | None = None) -> Job:
    """A finished analysis: one combination with a rooted tree and masked ancestors."""
    job = Job.objects.create(
        owner=owner or local_user(), name="ime2", params={"gap_methods": ["fitch"]}
    )
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combo = {"key": COMBO, "aligner": "mafft", "tree": "iqtree", "model": "modelfinder",
             "label": "MAFFT · IQ-TREE · ModelFinder"}  # fmt: skip
    files = {
        "summary.json": json.dumps({"combinations": [combo]}),
        f"{COMBO}/tree.nwk": "(a:0.1,b:0.2,(c:0.3,d:0.4)Node3:0.3)Node1;",
        f"{COMBO}/tree.rooted.nwk": ROOTED,
        f"{COMBO}/ancestors.fasta": ">Node1\nMKLV\n>Node2\nMKLA\n>Node3\nMKLW\n",
        f"{COMBO}/ancestors.gaps-fitch.fasta": ">Node1\nMK-V\n>Node2\nMKLA\n>Node3\nM-LW\n",
    }
    for name, text in files.items():
        fake_storage.objects[prefix + name] = text.encode()
    return job


def form(**changes: Any) -> dict[str, Any]:
    data = {
        "name": "test", "seed": 7, "tree_source": "random", "shape": "yh", "taxa": 10,
        "length_min": 0.02, "length_mean": 0.1, "length_max": 0.3, "birth": 0.1, "death": 0.05,
        "root_source": "random", "root_length": 200, "matrix": "LG", "rate": "G4",
        "alpha": 0.5, "invariant": 0.1, "insertion_rate": 0.03, "deletion_rate": 0.03,
        "indel_mean": 3, "replicates": 3,
    }  # fmt: skip
    return data | changes


@pytest.mark.django_db
def test_page_and_a_random_simulation(client: Client, fake_storage: FakeStorage) -> None:
    page = client.get("/simulations/").content.decode()
    assert "<h1>Simulations</h1>" in page
    assert 'aria-current="page" title="Evolve sequences on a known tree' in page
    # AliSim is a link to its paper; the long blurb is gone; an "On this page" sidebar.
    assert 'href="https://doi.org/10.1093/molbev/msac092"' in page and ">AliSim</a>" in page
    assert "keep the truth" not in page
    assert 'data-toc aria-label="On this page"' in page and 'href="#my-simulations"' in page
    response = client.post("/simulations/", form(indels="on"))
    job = Job.objects.get()
    assert response.status_code == 302 and response["Location"] == f"/jobs/{job.id}/"
    spec = parse(job.params["spec"])  # the stored spec is a valid spec
    assert job.params["simulate"] and spec.seeds == (7, 8, 9) and spec.names == "neutral"
    assert job.params["spec"]["tree"] == {"random": "yh", "taxa": 10, "lengths": [0.02, 0.1, 0.3]}
    assert job.params["spec"]["model"] == {"substitution": "LG+G4{0.5}"}
    assert job.params["spec"]["indels"] == {"rates": "0.03,0.03", "sizes": "GEO{3},GEO{3}"}
    inputs = f"jobs/{job.id}/inputs/"
    assert fake_storage.list_keys(inputs) == [inputs + "spec.json"]
    assert json.loads(fake_storage.objects[inputs + "spec.json"]) == job.params["spec"]
    assert job.work_units is not None and expected_tasks(job.params) == 9


@pytest.mark.django_db
def test_simulation_on_a_job_tree(client: Client, fake_storage: FakeStorage) -> None:
    source = finished_job(fake_storage)
    options = client.get("/simulations/trees/", {"job": str(source.id)}).content.decode()
    assert f'value="{COMBO}" selected' in options and '<option value="Node2"' in options
    data = form(tree_source="job", job=source.id, combination=COMBO, root_source="ancestor",
                node="Node3", rate="I+G4", invariant=0.2)  # fmt: skip
    assert client.post("/simulations/", data).status_code == 302
    job = Job.objects.exclude(id=source.id).get()
    inputs = f"jobs/{job.id}/inputs/"
    # Leaves renamed sim001...; the job's node names kept; the real names beside it.
    tree = fake_storage.objects[inputs + "sim_tree.nwk"].decode()
    assert "a:" not in tree and "sim001" in tree and "Node2" in tree
    names = fake_storage.objects[inputs + "names.tsv"].decode()
    assert names == "simulated\toriginal\nsim001\ta\nsim002\tb\nsim003\tc\nsim004\td\n"
    # The root: Node3 with the job's own gap placement (Fitch), gaps removed.
    assert fake_storage.objects[inputs + "sim_root.fasta"] == b">root\nMLW\n"
    assert job.params["spec"]["model"] == {"substitution": "LG+I{0.2}+G4{0.5}"}
    assert job.params["spec"]["simulation"]["names"] == "keep"
    assert job.params["simulation_of"] == {
        "kind": "job", "job": str(source.id), "combination": COMBO,
        "label": "ime2 · MAFFT · IQ-TREE · ModelFinder", "node": "Node3", "version": "fitch",
    }  # fmt: skip
    plan = simulations.describe(job.params)
    assert plan[0].startswith("simulate 3 replicates on the tree of ime2") and "Node3" in plan[1]


@pytest.mark.django_db
def test_build_params_passes_the_staged_files_only(tmp_path: Path) -> None:
    for name in ("spec.json", "sim_tree.nwk", "sim_root.fasta"):
        (tmp_path / name).write_text("x")
    job = Job.objects.create(owner=local_user(), params={
        "simulate": True, "spec": {"simulation": {"seeds": [1]}}, "simulation_of": {"x": 1},
    })  # fmt: skip
    params = build_params(job, tmp_path, tmp_path / "out")
    assert params["simulate"] == str(tmp_path / "spec.json")  # a path, never the stored True
    assert params["sim_tree"] == str(tmp_path / "sim_tree.nwk")
    assert params["sim_root"] == str(tmp_path / "sim_root.fasta")
    assert "spec" not in params and "simulation_of" not in params


@pytest.mark.django_db
@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"root_source": "ancestor"}, "A random tree has no ancestors"),
        ({"length_min": 0.1, "length_max": 0.1, "length_mean": 0.1}, "minimum below"),
        ({"shape": "bd", "birth": 0.1, "death": 0.2}, "death rate must be below"),
        ({"root_source": "paste", "root_text": ">x\nMKLB"}, "20 amino acids"),
        ({"replicates": 101}, "100"),
        ({"taxa": 400, "root_length": 5000, "replicates": 100}, "the limit is"),
        ({"tree_source": "job"}, "Choose one of your finished jobs"),
    ],
)
def test_bad_simulations_are_refused(
    client: Client, fake_storage: FakeStorage, changes: dict[str, Any], message: str
) -> None:
    page = client.post("/simulations/", form(**changes)).content.decode()
    assert message in page and not Job.objects.exists()


@pytest.mark.django_db
def test_another_users_job_tree_is_refused(client: Client, fake_storage: FakeStorage) -> None:
    other = finished_job(fake_storage, User.objects.create(username="other"))
    data = form(tree_source="job", job=other.id, combination=COMBO, root_source="random")
    assert "Choose one of your finished jobs" in client.post("/simulations/", data).content.decode()
    options = client.get("/simulations/trees/", {"job": str(other.id)}).content.decode()
    assert COMBO not in options


def finished_simulation(fake_storage: FakeStorage) -> Job:
    spec = simulations.build_spec(
        description="", tree={"random": "yh", "taxa": 3, "lengths": [0.1, 0.2, 0.3]},
        root={"length": 4}, model="LG", indels=None, seed=5, replicates=2,
    )  # fmt: skip
    job = Job.objects.create(owner=local_user(), name="sims", params={"simulate": True,
                             "spec": spec, "simulation_of": {"kind": "random"}})  # fmt: skip
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    for seed in (5, 6):
        folder = f"jobs/{job.id}/results/replicates/{seed}/"
        unaligned = b">sim001\nMKL\n>sim002\nMQL\n>sim003\nMKQL\n"
        fake_storage.objects[folder + "input.fasta"] = unaligned
        fake_storage.objects[folder + "true_alignment.fasta"] = (
            b">sim001\nMK-L\n>sim002\nM-QL\n>sim003\nMKQL\n"
        )
        fake_storage.objects[folder + "true_tree.nwk"] = b"(sim001:0.1,sim002:0.2,sim003:0.3)4;\n"
        fake_storage.objects[folder + "provenance.json"] = b"{}"
    job.refresh_from_db()
    return job


@pytest.mark.django_db
def test_replicates_join_my_alignments_as_simulated(fake_storage: FakeStorage) -> None:
    job = finished_simulation(fake_storage)
    records = alignments.register_job(job)
    assert [(r.job_key, r.simulated, r.sequences, r.columns) for r in records] == [
        ("replicate5", True, 3, 4), ("replicate6", True, 3, 4),
    ]  # fmt: skip
    assert alignments.register_job(job) == records  # idempotent
    rows = [alignments.Row("sim001", "", "MKL")]
    copy = alignments.store(local_user(), "edited", rows, "protein", 1,
                            origin=Alignment.Origin.EDIT, derived_from=records[0])  # fmt: skip
    assert copy.simulated  # a copy of simulated data is simulated too


@pytest.mark.django_db
def test_job_page_and_use_in_a_new_job(client: Client, fake_storage: FakeStorage) -> None:
    job = finished_simulation(fake_storage)
    alignments.register_job(job)
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Simulated Replicates" in page and "Duplicate and run again" not in page
    assert client.get("/jobs/new/", {"from": str(job.id)}).status_code == 404
    assert f"/jobs/new/?simulation={job.id}&amp;replicate=6" in page
    assert "true_alignment.fasta" in page and "View alignment" in page

    new = client.get("/jobs/new/", {"simulation": str(job.id), "replicate": "6"}).content.decode()
    assert "Step 2 holds replicate 6" in new and "&gt;sim001" in new
    assert f'name="simulated_from" value="{job.id}:6"' in new
    response = client.post("/jobs/new/", {
        "sequences_text": ">sim001\nMKLAAAAAAAA\n>sim002\nMQLAAAAAAAA\n>sim003\nMKQLAAAAAAA\n"
                          ">sim004\nMKQLAAAAAAW\n",
        "simulated_from": f"{job.id}:6", "aligners": ["mafft"], "seed": 1, "seqtype": "auto",
        "analyses": "both", "stop_after": "alignments",
    })  # fmt: skip
    assert response.status_code == 302, response.content.decode()[:2000]
    analysis = Job.objects.exclude(id=job.id).get()
    assert analysis.params["simulated"] is True
    assert analysis.params["simulated_from"] == {"job": str(job.id), "seed": 6, "name": "sims"}
    banner = client.get(f"/jobs/{analysis.id}/").content.decode()
    assert "This job ran on simulated sequences" in banner and "replicate\n    6" in banner


@pytest.mark.django_db
def test_simulated_from_must_be_the_users_simulation(
    client: Client, fake_storage: FakeStorage
) -> None:
    owner = User.objects.create(username="other")
    other = Job.objects.create(owner=owner, params={"simulate": True})
    query = {"simulation": str(other.id), "replicate": "1"}
    assert client.get("/jobs/new/", query).status_code == 404
    page = client.post("/jobs/new/", {
        "sequences_text": ">a\nMKL\n>b\nMKL\n>c\nMKL\n>d\nMKL\n", "simulated_from": f"{other.id}:1",
        "aligners": ["mafft"], "seed": 1,
    }).content.decode()  # fmt: skip
    assert "That simulation is not one of yours." in page and Job.objects.count() == 1


def test_replicate_files_are_described() -> None:
    keys = [f"replicates/5/{name}" for name in filelist.SIMULATION]
    tables = filelist.file_tables(keys)
    assert [t.title for t in tables] == ["Replicate (seed 5)"]
    assert all(description for _, _, description in tables[0].rows)


def test_model_strings() -> None:
    assert simulations.model_string("WAG", "", 0.5, 0.1) == "WAG"
    assert simulations.model_string("LG", "G4", 0.75, 0.1) == "LG+G4{0.75}"
    assert simulations.model_string("LG", "I+G4", 1, 0.25) == "LG+I{0.25}+G4{1}"
    with pytest.raises(ValueError):
        simulations.model_string("LG; rm", "G4", 1, 0)


@pytest.mark.django_db
def test_job_picker_lists_finished_jobs_with_trees(
    client: Client, monkeypatch: pytest.MonkeyPatch
) -> None:
    owner = local_user()
    made = {}
    for name, params, status in (
        ("ime2 full", {}, Job.Status.SUCCEEDED),
        ("ime2 aligned", {"stop_after": "alignments"}, Job.Status.SUCCEEDED),
        ("ime2 failed", {}, Job.Status.FAILED),
        ("a simulation", {"simulate": True}, Job.Status.SUCCEEDED),
        ("lysin trees", {"stop_after": "trees"}, Job.Status.SUCCEEDED),
    ):
        made[name] = Job.objects.create(owner=owner, name=name, params=params, status=status)
    Job.objects.create(owner=User.objects.create_user("eve"), name="ime2 theirs",
                       status=Job.Status.SUCCEEDED)  # fmt: skip

    def picked(query: str = "") -> list[str]:
        page = client.get(f"/simulations/jobs/?{query}").content.decode()
        return re.findall(r'data-label="([^"·]*) ·', page)

    assert picked() == ["lysin trees", "ime2 full"]  # newest first
    assert picked("q=IME2") == ["ime2 full"]
    assert picked(f"q={str(made['lysin trees'].id)[:6]}") == ["lysin trees"]
    assert (
        "No finished job with trees matches"
        in client.get("/simulations/jobs/?q=zz").content.decode()
    )

    page = client.get("/simulations/").content.decode()
    assert 'name="job" id="sim-job" value=""' in page and "Choose a job…" in page
    assert '<select name="job"' not in page  # no long menu

    monkeypatch.setattr("phylobot.simulation_views.PICKER_PAGE", 1)
    first = client.get("/simulations/jobs/").content.decode()
    assert first.count("data-pick=") == 1 and "page=1" in first
    assert picked("page=1") == ["ime2 full"]


@pytest.mark.django_db
def test_my_simulations_sort_and_search(client: Client) -> None:
    Job.objects.create(owner=local_user(), name="sim A", params={"simulate": True})
    page = client.get("/simulations/").content.decode()
    assert 'id="simulations-table"' in page and 'data-table-search="simulations-table"' in page
    assert '<button type="button" class="sort">Submitted</button>' in page
