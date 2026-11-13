"""Extend (phylobot.extend): a new job that builds on a finished job's alignments."""

import json
import re
from pathlib import Path
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import Client

from phylobot import (
    analogs,
    bundle,
    executors,
    extend,
    jobs,
    joined,
    properties_run,
    properties_submit,
    stages,
)
from phylobot.executors import PROPERTY_DIR
from phylobot.models import Job
from phylobot.modes import local_user
from phylobot.pipeline.summary import COMPARISON_FIELDS
from tests.conftest import FakeStorage
from tests.test_views import viewer_job

pytestmark = pytest.mark.django_db
FOUR = "".join(f">seq{i}\nMKTAYIAK{'LV'[i % 2]}Q\n" for i in range(4))
MAFFT = "".join(f">seq{i}\nMKTA-YIAK{'LV'[i % 2]}Q\n" for i in range(4))
MINE = "".join(f">seq{i}\nMKTAYIA-K{'LV'[i % 2]}Q\n" for i in range(4))


def finished(fake_storage: FakeStorage, owner: Any = None, **params: Any) -> Job:
    """An alignments-only job that ran MAFFT on FOUR, with one alignment of the user's."""
    params = {"aligners": ["mafft"], "stop_after": "alignments", "seqtype": "protein",
              "alignment_labels": ["mine"], "outgroup": ["seq3"], "seed": 7, **params}  # fmt: skip
    job = Job.objects.create(owner=owner or local_user(), name="parent", params=params)
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    job.refresh_from_db()
    stored = fake_storage.objects
    stored[f"jobs/{job.id}/inputs/sequences.fasta"] = FOUR.encode()
    stored[f"jobs/{job.id}/inputs/alignment1.fasta"] = MINE.encode()
    results = f"jobs/{job.id}/results/"
    stored[results + "summary.json"] = json.dumps({"alignments": [
        {"key": "provided1", "label": "mine", "path": "provided1/alignment.fasta"},
        {"key": "mafft", "label": "MAFFT", "path": "mafft/alignment.fasta"},
    ]}).encode()  # fmt: skip
    stored[results + "mafft/alignment.fasta"] = MAFFT.encode()
    stored[results + "provided1/alignment.fasta"] = MINE.encode()
    return job


def submit(client: Client, parent: Job, **fields: Any) -> Any:
    data = {"extend_job": str(parent.id), "aligners": ["mafft", "muscle"],
            "tree_builders": ["iqtree"], "modelfinder": "on", "run_trees": "on",
            "run_ancestors": "on", "bootstrap": 1000, "seed": 99,
            "sequences_text": ">x\nMKV\n", **fields}  # fmt: skip
    return client.post("/jobs/new/", data)


def test_the_page_locks_what_the_parent_gave(client: Client, fake_storage: FakeStorage) -> None:
    parent = finished(fake_storage)
    page = client.get(f"/jobs/new/?extend={parent.id}").content.decode()
    text = re.sub(r"\s+", " ", page)
    assert "Extending “parent”" in text and f'name="extend_job" value="{parent.id}"' in page
    assert "Kept from “parent”: <strong>MAFFT</strong>, <strong>mine</strong>." in text
    assert not re.search(r'name="aligners" value="mafft"', page)  # shown ticked, locked
    assert re.search(r'name="aligners" value="muscle"', page)  # free to add
    assert 'name="seed" value="7"' in page and 'type="hidden" name="seed"' in page
    assert 'data-defaults="sequences"' not in page  # "Use defaults" would unlock them
    assert re.search(r'name="outgroup" value="seq3"\s*checked', page)
    assert not Job.objects.exclude(id=parent.id).exists()  # nothing runs until Submit


def test_submit_queues_a_child_with_the_parents_alignments(
    client: Client, fake_storage: FakeStorage, tmp_path: Path
) -> None:
    parent = finished(fake_storage)
    parent.add_tags(["gba1"])
    response = submit(client, parent)
    child = Job.objects.exclude(id=parent.id).get()
    assert response["Location"] == f"/jobs/{child.id}/"
    assert child.extends_id == parent.id and child.rerun_of_id is None
    params = child.params
    # MAFFT kept under its key (and first: the primary), the user's own as provided1;
    # MUSCLE is the only aligner that runs. Posted sequences and seed are ignored.
    assert params["alignment_keys"] == ["mafft", "provided"]
    assert params["alignment_labels"] == ["MAFFT", "mine"]
    assert params["aligners"] == ["muscle"] and params["seed"] == 7
    assert params["stop_after"] == "ancestors"
    assert params["extends"] == {"job": str(parent.id), "name": "parent", "alignments": 2}
    stored = fake_storage.objects
    assert stored[f"jobs/{child.id}/inputs/sequences.fasta"] == FOUR.encode()
    assert stored[f"jobs/{child.id}/inputs/alignment1.fasta"] == MAFFT.encode()
    assert stored[f"jobs/{child.id}/inputs/alignment2.fasta"] == MINE.encode()
    assert list(child.tags.values_list("name", flat=True)) == ["gba1"]
    # The keys reach Nextflow (PARAM_ALLOWLIST), in the staged files' order.
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for n in (1, 2):
        (inputs / f"alignment{n}.fasta").write_bytes(
            stored[f"jobs/{child.id}/inputs/alignment{n}.fasta"]
        )
    sent = executors.build_params(child, inputs, tmp_path / "out")
    assert sent["alignment_keys"] == ["mafft", "provided"]
    assert [Path(p).name for p in sent["alignment"]] == ["alignment1.fasta", "alignment2.fasta"]
    assert "extends" not in sent  # display only
    # Two alignments of the parent's, one new: what the job page says it will do.
    page = re.sub(r"\s+", " ", client.get(f"/jobs/{child.id}/").content.decode())
    assert "Use the 2 alignments of “parent” → also align with MUSCLE" in page
    assert f'Extends <a href="/jobs/{parent.id}/">parent</a>' in page
    parent_page = re.sub(r"\s+", " ", client.get(f"/jobs/{parent.id}/").content.decode())
    assert f'Extended in <a href="/jobs/{child.id}/">' in parent_page
    parent.refresh_from_db()
    assert parent.status == Job.Status.SUCCEEDED and parent.params["stop_after"] == "alignments"


def test_extend_sits_in_the_workflow(client: Client, fake_storage: FakeStorage) -> None:
    """Below the steps graph in its panel (not as a step of its own: it is not one of
    this job's), and only there; the plan line sits in that panel too, below it."""
    parent = finished(fake_storage)
    page = client.get(f"/jobs/{parent.id}/").content.decode()
    assert page.count(f"?extend={parent.id}") == 1
    panel = page[
        page.index('id="job-stages"') : page.index("</details>", page.index('id="job-stages"'))
    ]
    assert f'class="button" href="/jobs/new/?extend={parent.id}"' in panel
    assert panel.index("</svg>") < panel.index("?extend=")  # below the graph
    assert panel.index("?extend=") < panel.index('class="job-plan"')
    assert "stage-extend" not in page and "more aligners, trees and ancestors" in panel
    assert 'class="head-action" href="/jobs/new/?extend=' not in page  # moved from the header


def test_nothing_to_add_is_refused(client: Client, fake_storage: FakeStorage) -> None:
    parent = finished(fake_storage)
    response = submit(client, parent, aligners=["mafft"], run_trees="")
    assert response.status_code == 200 and "Nothing to add" in response.content.decode()
    assert not Job.objects.exclude(id=parent.id).exists()


def test_more_aligners_alone_is_an_addition(client: Client, fake_storage: FakeStorage) -> None:
    parent = finished(fake_storage)
    assert submit(client, parent, run_trees="").status_code == 302
    child = Job.objects.exclude(id=parent.id).get()
    assert child.params["stop_after"] == "alignments" and child.params["aligners"] == ["muscle"]
    # ...and that child can be extended in turn: its inherited aligners come first.
    assert extend.counts(child) == (["mafft", "muscle"], 1)


def test_only_finished_alignment_jobs_of_your_own(
    client: Client, fake_storage: FakeStorage
) -> None:
    other = User.objects.create(username="other")
    cases = [
        finished(fake_storage, owner=other),  # someone else's: 404 like any job of theirs
        finished(fake_storage, stop_after="ancestors", reconcile=True),  # reconciled
    ]
    running = Job.objects.create(owner=local_user(), params={"stop_after": "alignments"})
    cases.append(running)
    for job in cases:
        assert client.get(f"/jobs/new/?extend={job.id}").status_code == 404
        assert submit(client, job).status_code == 404
    assert client.get("/jobs/new/?extend=not-a-uuid").status_code == 404
    page = client.get(f"/jobs/{cases[1].id}/").content.decode()
    assert "?extend=" not in page  # no button where it cannot work
    assert (
        f"?extend={finished(fake_storage).id}"
        in client.get(f"/jobs/{Job.objects.latest('created_at').id}/").content.decode()
    )


def test_refusals_say_why() -> None:
    job = Job(params={"stop_after": "alignments"}, status=Job.Status.SUCCEEDED)
    assert extend.refusal(job) == ""
    assert "finished" in extend.refusal(Job(params=job.params, status=Job.Status.FAILED))
    assert "imported" in extend.refusal(
        Job(params=job.params, status=Job.Status.SUCCEEDED, imported_at="2026-10-22T00:00Z")
    )
    assert "analysis" in extend.refusal(Job(params={"simulate": {"replicates": 1}}))


def test_coding_dna_inherits_the_codon_alignments(fake_storage: FakeStorage) -> None:
    parent = finished(fake_storage, seqtype="codon", alignment_labels=[])
    results = f"jobs/{parent.id}/results/"
    fake_storage.objects[results + "summary.json"] = json.dumps(
        {"alignments": [{"key": "mafft", "label": "MAFFT", "path": "mafft/alignment.fasta"}]}
    ).encode()
    fake_storage.objects[results + "mafft/alignment.cds.fasta"] = b">seq0\nATG\n"
    inherited = extend.inherit(parent)
    assert [(a.key, a.content) for a in inherited.alignments] == [("mafft", b">seq0\nATG\n")]
    del fake_storage.objects[results + "mafft/alignment.cds.fasta"]
    with pytest.raises(LookupError):
        extend.inherit(parent)


def test_the_preview_counts_the_inherited_alignments(
    client: Client, fake_storage: FakeStorage
) -> None:
    parent = finished(fake_storage)
    fields = {"sequences_text": FOUR, "aligners": ["muscle"], "run_trees": "on",
              "modelfinder": "on", "tree_builders": ["iqtree"],
              "extend_job": str(parent.id)}  # fmt: skip
    page = re.sub(r"\s+", " ", client.post("/preview/", fields).content.decode())
    assert "use the 2 alignments of “parent”" in page.lower()
    assert "<strong>3</strong> combinations" in page  # MAFFT, mine and MUSCLE


# Extending a job that has trees (results mode): the parent's results are read in place.


def comparison_row(aligner: str, model: str) -> str:
    values = {"aligner": aligner, "tree": "iqtree", "model_key": model, "model": model,
              "data": "amino", "rank_bic": "1"}  # fmt: skip
    return "\t".join(values.get(f, "0") for f in COMPARISON_FIELDS) + "\n"


def split(fake_storage: FakeStorage, gap_methods: list[str] | None = None) -> tuple[Job, Job]:
    """`viewer_job` (MAFFT and MUSCLE) as a MAFFT job and an extension that ran MUSCLE:
    the MUSCLE files move to the extension, each summary keeps its own combinations."""
    parent = viewer_job(fake_storage, gap_methods)
    parent.refresh_from_db()  # marked succeeded in the database
    parent.name = "parent"
    parent.params |= {"stop_after": "ancestors", "aligners": ["mafft"], "seed": 7,
                      "tree_builders": ["iqtree"], "modelfinder": True}  # fmt: skip
    parent.save()
    child = Job.objects.create(
        owner=local_user(), name="child", extends=parent,
        params={**parent.params, "aligners": ["muscle"],
                "inherits": {"job": str(parent.id), "alignments": ["mafft"]}},
    )  # fmt: skip
    assert jobs.mark(child.id, Job.Status.RUNNING) and jobs.mark(child.id, Job.Status.SUCCEEDED)
    child.refresh_from_db()
    above, below = f"jobs/{parent.id}/results/", f"jobs/{child.id}/results/"
    stored = fake_storage.objects
    whole = json.loads(stored[above + "summary.json"])
    for key in [k for k in stored if k.startswith(above + "muscle/")]:
        stored[below + key[len(above) :]] = stored.pop(key)
    for prefix, aligner in ((above, "mafft"), (below, "muscle")):
        mine = [c for c in whole["combinations"] if c["aligner"] == aligner]
        stored[prefix + "summary.json"] = json.dumps({
            "primary": mine[0]["key"], "combinations": mine,
            "alignments": [{"key": aligner, "label": aligner.upper(),
                            "path": f"{aligner}/alignment.fasta"}],
        }).encode()  # fmt: skip
        stored[prefix + f"{aligner}/alignment.fasta"] = b">a\nMK-\n>b\nM-K\n>c<i>\nMKK\n"
        stored[prefix + "model_comparison.tsv"] = (
            "\t".join(COMPARISON_FIELDS)
            + "\n"
            + "".join(comparison_row(c["aligner"], c["model"]) for c in mine)
        ).encode()
        stored[prefix + "versions.tsv"] = f"tool\tversion\n{aligner}\tv1\niqtree\t3.1.3\n".encode()
    stored[f"jobs/{parent.id}/inputs/sequences.fasta"] = FOUR.replace("seq0", "a").encode()
    fake_storage.map_ancestors(parent)  # each job's map covers what it ran
    fake_storage.map_ancestors(child)
    return parent, child


def test_an_extension_reads_as_one_job(fake_storage: FakeStorage) -> None:
    parent, child = split(fake_storage)
    found = joined.results(child)
    assert [c["key"] for c in found.summary["combinations"]] == [
        "mafft/iqtree/WAG", "mafft/iqtree/modelfinder", "muscle/iqtree/modelfinder",
    ]  # fmt: skip
    assert found.summary["primary"] == "mafft/iqtree/modelfinder"  # the parent's
    assert [a["key"] for a in found.summary["alignments"]] == ["mafft", "muscle"]
    # Files live where they were written; nothing was copied into the extension.
    assert (
        found.key("mafft/iqtree/WAG/tree.nwk")
        == f"jobs/{parent.id}/results/mafft/iqtree/WAG/tree.nwk"
    )
    assert found.key("muscle/iqtree/modelfinder/tree.nwk").startswith(f"jobs/{child.id}/")
    assert not fake_storage.list_keys(f"jobs/{child.id}/results/mafft/")
    # The shared tables hold both jobs' rows, the parent's first, each row once.
    assert [r.split("\t")[:3] for r in found.read("model_comparison.tsv").splitlines()] == [
        ["aligner", "tree", "model_key"], ["mafft", "iqtree", "modelfinder"],
        ["mafft", "iqtree", "WAG"], ["muscle", "iqtree", "modelfinder"],
    ]  # fmt: skip
    assert found.read("versions.tsv").splitlines() == [
        "tool\tversion",
        "mafft\tv1",
        "iqtree\t3.1.3",
        "muscle\tv1",
    ]  # the inherited aligner's version too (B1's gap, closed here)
    # The parent itself is untouched, and reads as it always did.
    assert [c["aligner"] for c in joined.results(parent).summary["combinations"]] == ["mafft"] * 2
    # The ancestor map: the parent's rows for the inherited alignment, then the child's.
    mapped = {r.combination for r in analogs.read_map(found.read(analogs.MAP))}
    assert mapped == {"mafft/iqtree/WAG", "mafft/iqtree/modelfinder", "muscle/iqtree/modelfinder"}
    # Sizes come from the manifest of whichever job holds each file.
    sizes = found.sizes("mafft/iqtree/WAG/")
    tree = f"jobs/{parent.id}/results/mafft/iqtree/WAG/tree.nwk"
    assert sizes["mafft/iqtree/WAG/tree.nwk"] == len(fake_storage.objects[tree])
    assert set(found.sizes()) == found.keys


def test_its_pages_show_both_alignments(client: Client, fake_storage: FakeStorage) -> None:
    parent, child = split(fake_storage, gap_methods=["fitch"])
    page = re.sub(r"\s+", " ", client.get(f"/jobs/{child.id}/").content.decode())
    assert "3 combinations" in page
    assert f"/jobs/{child.id}/ancestors/?combination=mafft%2Fiqtree%2FWAG" in page
    assert f"/jobs/{child.id}/ancestors/?combination=muscle%2Fiqtree%2Fmodelfinder" in page
    assert "the results of its mafft alignment are read from it" in page
    viewer = client.get(f"/jobs/{child.id}/ancestors/?combination=mafft/iqtree/WAG").content
    assert b'<option value="mafft/iqtree/WAG" selected>' in viewer
    assert client.get(f"/jobs/{child.id}/tree/?combination=mafft/iqtree/WAG").status_code == 200
    # Files: an inherited one streams from the parent, the summary is the joined one.
    file = client.get(f"/jobs/{child.id}/files/mafft/iqtree/WAG/tree.nwk")
    assert b"".join(file.streaming_content) == b"((a,b)Node2,c<i>)Node1;"  # type: ignore[attr-defined]
    summary = client.get(f"/jobs/{child.id}/files/summary.json")
    body = b"".join(summary.streaming_content)  # type: ignore[attr-defined]
    assert len(json.loads(body)["combinations"]) == 3
    # The ZIP stands alone: the inherited files are in it, read from the parent.
    zipped = {path: key for key, path in bundle.keys(child)}
    name = next(p for p in zipped if p.endswith("results/mafft/iqtree/WAG/tree.nwk"))
    assert zipped[name].startswith(f"jobs/{parent.id}/")


def test_properties_measured_from_the_extension_belong_to_the_files(
    client: Client, fake_storage: FakeStorage
) -> None:
    """Measuring MAFFT's and MUSCLE's Node1 from the extension: each is staged from the
    job that holds it, and both jobs show the numbers for the files they read."""
    parent, child = split(fake_storage, gap_methods=["fitch"])
    response = client.post(f"/jobs/{child.id}/properties/measure/", {
        "combination": "mafft/iqtree/modelfinder", "node": "Node1", "gaps": "fitch",
        properties_run.FIELD: ["charge"], properties_run.SCOPE_FIELD: "analogs",
    })  # fmt: skip
    assert response.status_code == 302
    run = Job.objects.filter(params__properties_only=True).get()
    of = run.params["properties_of"]
    assert of["parent"] == str(child.id)
    assert of["owners"] == {
        "mafft/iqtree/WAG": str(parent.id),
        "mafft/iqtree/modelfinder": str(parent.id),
        "muscle/iqtree/modelfinder": str(child.id),
    }
    staged = f"jobs/{run.id}/inputs/{PROPERTY_DIR}/"
    assert fake_storage.objects[staged + "mafft/iqtree/WAG/ancestors.gaps-fitch.fasta"]
    assert fake_storage.objects[staged + "muscle/iqtree/modelfinder/ancestors.gaps-fitch.fasta"]
    # The parent's pages read the run too, for the files it holds; not for MUSCLE's.
    assert run in properties_submit.runs_for(parent)
    assert properties_submit.measured_in(run, "mafft/iqtree/WAG") == str(parent.id)
    assert properties_submit.measured_in(run, "muscle/iqtree/modelfinder") == str(child.id)


def test_deleting_the_parent_warns_and_the_extension_says_what_it_lost(
    client: Client, fake_storage: FakeStorage
) -> None:
    parent, child = split(fake_storage)
    page = client.get(f"/jobs/{parent.id}/").content.decode()
    assert 'data-dependents="child"' in page
    listing = client.get("/jobs/").content.decode()
    assert re.search(
        rf'value="{parent.id}"\s+data-name="parent"[^>]*data-dependents="child"', listing
    )
    assert 'data-dependents=""' in client.get(f"/jobs/{child.id}/").content.decode()
    jobs.delete([parent])
    child.refresh_from_db()
    found = joined.results(child)
    assert found.lost and [c["aligner"] for c in found.summary["combinations"]] == ["muscle"]
    page = re.sub(r"\s+", " ", client.get(f"/jobs/{child.id}/").content.decode())
    assert "The job this one extends was deleted" in page
    assert "can't be extended" not in page and "?extend=" not in page
    assert "incomplete" in extend.refusal(child)


def full_parent(fake_storage: FakeStorage, **params: Any) -> Job:
    """A finished MAFFT job with trees and ancestors (split's parent, without its child)."""
    parent, child = split(fake_storage)
    parent.params |= params
    parent.save()
    jobs.delete([child])
    return parent


def test_submit_takes_every_later_setting_from_the_parent(
    client: Client, fake_storage: FakeStorage
) -> None:
    parent = full_parent(
        fake_storage, models=["WAG+G4"], fixed_models=True, candidates=["LG+G4", "JTT+G4"],
        gap_methods=["fitch"], outgroup=["a"], alrt=1000, bootstrap=0, seqtype="protein",
    )  # fmt: skip
    page = re.sub(r"\s+", " ", client.get(f"/jobs/new/?extend={parent.id}").content.decode())
    assert "Their trees and ancestors are kept too" in page
    assert "Trees and Ancestors</h2>" in page and "Add to every alignment" in page
    # Posted later-step settings are ignored: a different setting would make the new
    # alignment's results incomparable with the kept ones (models and tree sources are
    # additions instead: test_adding_a_tree_source_and_a_model_runs_only_what_is_missing).
    response = submit(
        client, parent, aligners=["muscle"], gap_methods=[], sh_alrt="", outgroup=["seq1"]
    )
    assert response.status_code == 302
    child = Job.objects.filter(extends=parent, deleted_at__isnull=True).get()
    assert response["Location"] == f"/jobs/{child.id}/"
    params = child.params
    assert params["aligners"] == ["muscle"]
    assert params["inherits"] == {"job": str(parent.id), "alignments": ["mafft"]}
    for key in ("models", "candidates", "gap_methods", "outgroup", "alrt", "bootstrap", "seed",
                "stop_after", "tree_builders", "modelfinder", "fixed_models"):  # fmt: skip
        assert params[key] == parent.params[key], key
    assert "alignment_labels" not in params and "alignment_keys" not in params
    inputs = {k.rsplit("/", 1)[1] for k in fake_storage.list_keys(f"jobs/{child.id}/inputs/")}
    assert inputs == {"sequences.fasta"}  # the alignments are read, not given
    sent = executors.build_params(child, Path("/nonexistent"), Path("/out"))
    assert "inherits" not in sent and sent["aligners"] == ["muscle"]


def test_what_cannot_be_added_to_a_job_with_trees(
    client: Client, fake_storage: FakeStorage
) -> None:
    parent = full_parent(fake_storage)
    response = submit(client, parent, aligners=["mafft"])
    assert "Nothing to add" in response.content.decode()
    mine = SimpleUploadedFile("mine.fasta", MINE.encode())
    response = submit(client, parent, aligners=["muscle"], alignment_files=[mine])
    assert "Only aligners can be added" in response.content.decode()
    reconciled = full_parent(fake_storage, reconcile=True)
    assert "reconciled" in extend.refusal(reconciled)
    assert client.get(f"/jobs/new/?extend={reconciled.id}").status_code == 404


def test_the_preview_counts_only_the_new_combinations(
    client: Client, fake_storage: FakeStorage
) -> None:
    parent = full_parent(fake_storage, fixed_models=True, models=["WAG+G4"], modelfinder=True)
    fields = {"extend_job": str(parent.id), "aligners": ["muscle"], "sequences_text": FOUR}
    page = re.sub(r"\s+", " ", client.post("/preview/", fields).content.decode())
    assert "keep the results of the 1 alignment of “parent”" in page.lower()
    assert "Also align with MUSCLE" in page
    assert "<strong>2</strong> combinations" in page  # MUSCLE x (ModelFinder, WAG+G4)


# Adding tree sources and models to every alignment, the inherited ones too.


def with_picks(fake_storage: FakeStorage, parent: Job) -> None:
    """The parent's MAFFT alignment of its own sequences, and its ModelFinder pick."""
    results = f"jobs/{parent.id}/results/"
    fake_storage.objects[results + "mafft/alignment.fasta"] = MAFFT.replace("seq0", "a").encode()
    combo = results + "mafft/iqtree/modelfinder/"
    fake_storage.objects[combo + "model.txt"] = b"LG+G4\n"
    fake_storage.objects[combo + "modelfinder.tsv"] = b"model\tbic\nLG+G4\t100\n"


def test_adding_a_tree_source_and_a_model_runs_only_what_is_missing(
    client: Client, fake_storage: FakeStorage, tmp_path: Path
) -> None:
    parent = full_parent(fake_storage, models=["WAG"], fixed_models=True, modelfinder=True)
    with_picks(fake_storage, parent)
    page = client.get(f"/jobs/new/?extend={parent.id}").content.decode()
    assert re.search(r'name="models" value="WAG\+G4"[^>]*>', page)  # free to add
    response = submit(
        client, parent, aligners=[], tree_builders=["raxml"], models=["JTT+G4"],
        modelfinder="on",
    )  # fmt: skip
    assert response.status_code == 302, re.findall(r"<p>(.*?)</p>", response.content.decode())[:3]
    child = Job.objects.filter(extends=parent, deleted_at__isnull=True).get()
    params = child.params
    assert params["aligners"] == ["none"]
    assert params["tree_builders"] == ["iqtree", "raxml"]
    assert params["models"] == ["WAG", "JTT+G4"] and params["fixed_models"] is True
    # The parent's alignment comes again, under its key; what it has is not run again.
    assert params["alignment_keys"] == ["mafft"] and params["alignment_labels"] == ["MAFFT"]
    assert params["skip_combinations"] == ["mafft/iqtree/WAG", "mafft/iqtree/modelfinder"]
    assert params["modelfinder_reused"] == ["mafft"]  # RAxML-NG's ModelFinder combination
    stored = fake_storage.objects
    inputs = f"jobs/{child.id}/inputs/"
    assert stored[inputs + "alignment1.fasta"] == MAFFT.replace("seq0", "a").encode()
    assert stored[inputs + "picks/mafft.model.txt"] == b"LG+G4\n"
    assert stored[inputs + "picks/mafft.modelfinder.tsv"].startswith(b"model\tbic")
    # The worker hands the pipeline the skip list, the reused picks and their folder.
    staged = tmp_path / "inputs"
    (staged / "picks").mkdir(parents=True)
    (staged / "alignment1.fasta").write_bytes(stored[inputs + "alignment1.fasta"])
    sent = executors.build_params(child, staged, tmp_path / "out")
    assert sent["skip_combinations"] == params["skip_combinations"]
    assert sent["modelfinder_reused"] == ["mafft"]
    assert sent["modelfinder_picks"] == str(staged / "picks")
    # Four combinations run: IQ-TREE under JTT+G4, RAxML-NG under all three; no MODELFINDER.
    assert executors.combination_count(params) == 4
    assert "MODELFINDER" not in stages.expected_kinds(params)


def test_the_preview_counts_the_additions(client: Client, fake_storage: FakeStorage) -> None:
    parent = full_parent(fake_storage, models=["WAG"], fixed_models=True, modelfinder=True)
    with_picks(fake_storage, parent)
    fields = {"extend_job": str(parent.id), "sequences_text": FOUR, "aligners": ["muscle"],
              "tree_builders": ["raxml"], "models": ["JTT+G4"]}  # fmt: skip
    page = re.sub(r"\s+", " ", client.post("/preview/", fields).content.decode())
    # MUSCLE x (iqtree, raxml) x (modelfinder, WAG, JTT+G4) + MAFFT's four new ones.
    assert "<strong>10</strong> combinations" in page
    assert "only what is new (2 combinations read, not rerun; 1 ModelFinder pick reused)" in page


def test_the_joined_ranking_spans_both_jobs() -> None:
    """A model added to an inherited alignment changes that alignment's ranking."""
    fit = {"model": "x", "log_likelihood": "-1", "free_parameters": "1", "aic": "1",
           "aicc": "1", "sites": "10", "tree_length": "1"}  # fmt: skip

    def combo(key: str, bic: float) -> dict[str, Any]:
        aligner, tree, model = key.split("/")
        return {"key": key, "aligner": aligner, "tree": tree, "model": model,
                "fit": {**fit, "model": model, "bic": str(bic)}, "rank_bic": 1}  # fmt: skip

    parent: dict[str, Any] = {
        "primary": "mafft/iqtree/LG",
        "alignments": [{"key": "mafft"}],
        "combinations": [combo("mafft/iqtree/LG", 100), combo("mafft/iqtree/WAG", 110)],
    }
    own = {"primary": "mafft/raxml/LG", "combinations": [combo("mafft/raxml/LG", 90)]}
    summary = joined.join_summaries(parent, own, {"mafft"})
    rows = joined.rerank(summary)
    ranks = {c["key"]: c["rank_bic"] for c in summary["combinations"]}
    assert ranks == {"mafft/raxml/LG": 1, "mafft/iqtree/LG": 2, "mafft/iqtree/WAG": 3}
    assert summary["primary"] == "mafft/iqtree/LG" and len(rows) == 3
    assert parent["combinations"][0]["rank_bic"] == 1  # the parent's own summary unchanged
