import csv
import html
import io
import json
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.parse import quote

import pytest
from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import Client, override_settings
from django.utils import timezone

from phylobot import (
    analogs,
    eta,
    executors,
    jobs,
    properties,
    properties_run,
    property_compare,
    views,
)
from phylobot.durations import seconds
from phylobot.models import Job, JobAttempt, SiteSettings, UserLimits
from phylobot.modes import local_user
from phylobot.settings import middleware_for
from tests.conftest import FakeStorage

pytestmark = pytest.mark.django_db
ROOT = Path(__file__).resolve().parent.parent
IME2 = (ROOT / "src" / "phylobot" / "examples" / "ime2.fasta").read_text()
FOUR = "".join(f">seq{i}\nMKTAYIAK{'LV'[i % 2]}Q\n" for i in range(4))


def post(client: Client, **fields: Any) -> Any:
    data = {"aligners": ["mafft"], "tree_builders": ["iqtree"], "modelfinder": "on",
            "run_trees": "on", "run_ancestors": "on", "bootstrap": 1000, "seed": 1,
            "codon_search": "on", **fields}  # fmt: skip
    return client.post("/jobs/new/", data)


def test_submit_creates_job_and_stores_input(client: Client, fake_storage: FakeStorage) -> None:
    response = post(
        client, name="My run", sequences_text=FOUR, candidates=["LG+G4", "WAG+R4"],
        models=["JTT"], gap_methods=["parsimony", "ml"],
    )  # fmt: skip
    job = Job.objects.get()
    assert response.status_code == 302 and response["Location"] == f"/jobs/{job.id}/"
    assert job.owner == local_user() and job.status == Job.Status.QUEUED and job.name == "My run"
    assert job.params["candidates"] == ["LG+G4", "WAG+R4"] and job.params["modelfinder"] is True
    assert job.params["models"] is None  # full searches unticked: that grid is ignored
    assert job.params["gap_methods"] == ["parsimony", "ml"]
    assert job.work_units and job.work_units > 0  # the estimate, kept beside cpu_seconds
    assert job.params["sequence_count"] == 4  # for My Jobs
    assert fake_storage.objects[f"jobs/{job.id}/inputs/sequences.fasta"] == FOUR.encode()


def test_submit_with_raxml(client: Client, fake_storage: FakeStorage) -> None:
    post(client, sequences_text=FOUR, tree_builders=["iqtree", "raxml"], sh_alrt="on")
    job = Job.objects.get()
    assert job.params["tree_builders"] == ["iqtree", "raxml"]
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "IQ-TREE and RAxML-NG (best model by ModelFinder) with SH-aLRT support" in page
    assert "2 combinations" in page


def test_branch_support_choices(client: Client, fake_storage: FakeStorage) -> None:
    page = client.get("/jobs/new/").content.decode()
    boxes = {name: tag for tag in re.findall(r"<input[^>]*>", page)
             for name in re.findall(r'name="(sh_alrt|ufboot)"', tag)}  # fmt: skip
    assert "checked" in boxes["sh_alrt"] and "checked" not in boxes["ufboot"]  # owner decision
    assert re.search(r'name="alrt" value="1000"', page)

    post(client, sequences_text=FOUR, sh_alrt="on", alrt=2000, ufboot="on")
    assert {k: Job.objects.get().params[k] for k in ("alrt", "bootstrap")} == {
        "alrt": 2000,
        "bootstrap": 1000,  # the helper posts 1000
    }
    Job.objects.all().delete()
    post(client, sequences_text=FOUR, sh_alrt="on", alrt="")  # blank count: the default
    assert {k: Job.objects.get().params[k] for k in ("alrt", "bootstrap")} == {
        "alrt": 1000,
        "bootstrap": 0,  # unticked: its count is ignored
    }
    # IQ-TREE refuses fewer than 1000 replicates
    response = post(client, sequences_text=FOUR, sh_alrt="on", alrt=500)
    assert response.status_code == 200 and Job.objects.count() == 1


def test_submit_with_file_upload(client: Client, fake_storage: FakeStorage) -> None:
    upload = SimpleUploadedFile("mine.fasta", FOUR.encode())
    assert post(client, sequences_file=upload).status_code == 302


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        ({"sequences_text": ">a\nMK\n>b\nMK\n"}, "Input error: need at least 4 sequences"),
        ({"sequences_text": FOUR.replace("MKTAY", "MK1AY", 1)}, "invalid characters: 1"),
        ({}, "Paste or upload sequences, or upload or choose your own alignment"),
        ({"sequences_text": FOUR, "modelfinder": ""}, "Choose ModelFinder, a full search"),
        ({"sequences_text": FOUR, "fixed_models": "on"}, "Tick the models to search in full"),
        ({"sequences_text": FOUR, "candidates": ["LG+G4;rm -rf /"]}, "Select a valid choice"),
        ({"sequences_text": FOUR, "models": ["LG+G4;rm -rf /"]}, "Select a valid choice"),
        ({"sequences_text": FOUR, "aligners": ["clustalw"]}, "Select a valid choice"),
    ],
)
def test_submit_rejects_bad_input(
    client: Client, fake_storage: FakeStorage, fields: dict[str, Any], message: str
) -> None:
    response = post(client, **fields)
    assert response.status_code == 200
    assert message in response.content.decode()
    assert not Job.objects.exists() and not fake_storage.objects


def test_per_user_upload_and_sequence_limits(client: Client, fake_storage: FakeStorage) -> None:
    UserLimits.objects.create(user=local_user(), max_upload_bytes=50, max_sequences=3)
    assert "your limit is 50" in post(client, sequences_text=FOUR).content.decode()
    UserLimits.objects.filter(user=local_user()).update(max_upload_bytes=None)
    assert "at most 3 sequences" in post(client, sequences_text=FOUR).content.decode()


def test_example_loads_ime2(client: Client) -> None:
    response = client.get("/examples/ime2/")
    assert response.status_code == 200
    assert "Saccharomyces.cerevisiae.IME2" in response.content.decode()
    assert client.get("/examples/nope/").status_code == 404


def test_preview_explains_combinations_and_candidates(client: Client) -> None:
    data = {"aligners": ["mafft"], "tree_builders": ["iqtree"], "modelfinder": "on",
            "fixed_models": "on", "models": ["LG+G4", "WAG+R4"], "run_trees": "on",
            "candidates": ["LG+G4", "WAG+R4"], "run_ancestors": "on", "gap_methods": ["fitch"],
            "sh_alrt": "on"}  # fmt: skip
    text = client.post("/preview/", data).content.decode()
    assert "<strong>3</strong> combinations" in text
    assert "LG+G4, LG+R4, WAG+G4, WAG+R4" in text and "more than the cells you ticked" in text
    assert (
        "<li>Align with MAFFT</li><li>Infer trees with IQ-TREE (best model by ModelFinder and "
        "each of 2 ticked models) "
        "with SH-aLRT support</li>"
        "<li>Reconstruct ancestral sequences with gap placement by Parsimony (Fitch)</li>"
    ) in text
    assert "work units" in text and '<span class="seg modelfinder"' in text
    data = {"aligners": ["mafft"], "tree_builders": ["iqtree"], "modelfinder": "on",
            "run_trees": "on"}  # fmt: skip
    text = client.post("/preview/", data).content.decode()
    assert "IQ-TREE's default candidate models" in text
    assert ", no branch support</li><li>Stop (no ancestral reconstruction)</li>" in text
    assert '<span class="seg ancestors"' not in text


def test_preview_counts_your_own_files_and_early_stops(client: Client) -> None:
    data: dict[str, Any] = {
        "aligners": ["mafft", "muscle"],
        "alignment_files": [SimpleUploadedFile("mine.fasta", b">a\nM-K\n"),
                            SimpleUploadedFile("yours.fasta", b">a\nMK-\n")],
    }  # fmt: skip
    text = client.post("/preview/", data).content.decode()
    assert (
        "<li>Align with MAFFT, MUSCLE and your 2 alignments</li><li>Stop (alignments only)</li>"
    ) in text
    assert "<strong>4</strong> alignments" in text
    data = {"tree_builders": [], "modelfinder": "on", "run_trees": "on", "run_ancestors": "on",
            "alignment_files": [SimpleUploadedFile("mine.fasta", b">a\nM-K\n")],
            "tree_files": [SimpleUploadedFile("t.nwk", b"(a,b,c);")]}  # fmt: skip
    text = client.post("/preview/", data).content.decode()
    assert "<li>Use your alignment</li><li>Fit your tree (best model by ModelFinder) (" in text
    assert "<strong>1</strong> combination\n" in text


def test_several_sequence_files_are_combined(client: Client, fake_storage: FakeStorage) -> None:
    first = ">seq0\nMKTAYIAKQR\n>seq1\nMKTAYIAKQL\n>seq2\nMKTAYLAKQR\n"
    second = ">seq2\nMKTAYLAKQR\n>seq3\nMKTAHIAKQR\n>seq1\nMKTAYIAKQV\n"  # seq2 again, seq1 not
    files = [upload("a.fasta", first), upload("b.fasta", second)]
    response = post(client, sequences_file=files)
    assert response.status_code == 200 and not Job.objects.exists()
    assert "1 sequence name is given to different sequences" in response.content.decode()
    files = [upload("a.fasta", first), upload("b.fasta", second)]
    assert post(client, sequences_file=files, duplicate_names="rename").status_code == 302
    stored = fake_storage.objects[f"jobs/{Job.objects.get().id}/inputs/sequences.fasta"].decode()
    assert stored.split()[::2] == [">seq0", ">seq1", ">seq2", ">seq3", ">seq1.2"]


def test_unticked_files_are_left_out(client: Client, fake_storage: FakeStorage) -> None:
    # b.fasta is listed but unticked: only a.fasta is used, so there are too few sequences.
    files = [upload("a.fasta", FOUR), upload("b.fasta", ">x\nMKTAYIAKQR\n")]
    response = post(client, sequences_file=files, sequences_listed=["a.fasta", "b.fasta"],
                    sequences_use=["a.fasta"])  # fmt: skip
    assert response.status_code == 302
    stored = fake_storage.objects[f"jobs/{Job.objects.get().id}/inputs/sequences.fasta"]
    assert b">x" not in stored


def test_sequence_check_lists_the_files(client: Client) -> None:
    files = [upload("a.fasta", FOUR), upload("b.fasta", FOUR)]
    page = check(client, sequences_file=files, sequences_listed=["a.fasta", "b.fasta"],
                 sequences_use=["a.fasta", "b.fasta"])  # fmt: skip
    assert 'id="sequence-files"' in page and 'hx-swap-oob="true"' in page
    assert 'a.fasta</span></label> <span class="hint">4 protein sequences, 10 residues' in page
    assert "4 sequences appear again with identical residues" in page
    assert "✓ 4 protein sequences" in page


# An aligned set that names two different sequences MonoDLG1- (as the MAGUK alignment of
# Anderson et al. 2016 does), invented so no real dataset is needed.
MAGUK = (
    ">MonoDLG1-\nMKT--AYIAKQRQISFVKSHFSRQ\n"
    ">ApisMPP1-\nMKTLLAYIAKQRQ--FVKSHWSRQ\n"
    ">MonoDLG1-\nMRT--AYLAKQRQISFVRSHFSKQ\n"
    ">HomoDLG4-\nMKS--GYIAKERQISYVKSHFTRQ\n"
    ">DrosDLG1-\nMKT--AFIGKQRNISFVKAHFSRE\n"
)


def test_a_name_used_twice_is_resolved_not_refused(
    client: Client, fake_storage: FakeStorage
) -> None:
    # The MAGUK alignment names two different sequences MonoDLG1-.
    page = check(client, sequences_text=MAGUK)
    assert "1 sequence name is given to different sequences" in page
    assert "<strong>MonoDLG1-</strong> twice in your pasted sequences" in page
    assert "<code>MonoDLG1-.2</code>" in page and "✓ 4 protein sequences" in page
    assert "✓ 5 protein sequences" in check(client, sequences_text=MAGUK,
                                              duplicate_names="rename")  # fmt: skip
    response = post(client, sequences_text=MAGUK)
    assert response.status_code == 200 and not Job.objects.exists()
    assert "given to different sequences: choose" in response.content.decode()
    # The same sequence twice is used once, without asking.
    twice = FOUR + FOUR.split("\n", 2)[0] + "\n" + FOUR.split("\n", 2)[1] + "\n"
    assert "1 sequence appears again with identical residues" in check(client, sequences_text=twice)
    assert post(client, sequences_text=twice).status_code == 302
    # As an uploaded alignment: the same choice, in step 3.
    files = {"alignment_files": [upload("maguk.fasta", MAGUK)]}
    page = re.sub(r"\s+", " ", client.post("/alignment-files/", files).content.decode())
    assert "1 sequence name is used for different rows" in page
    assert "<strong>MonoDLG1-</strong> twice in maguk.fasta" in page
    data = {"aligners": [], "run_trees": "", "run_ancestors": ""}
    response = post(client, **data, alignment_files=[upload("maguk.fasta", MAGUK)])
    assert "uses a sequence name more than once" in response.content.decode()
    response = post(client, **data, alignment_duplicates="rename",
                    alignment_files=[upload("maguk.fasta", MAGUK)])  # fmt: skip
    assert response.status_code == 302, response.content.decode()[-1500:]
    job = Job.objects.latest("created_at")
    stored = fake_storage.objects[f"jobs/{job.id}/inputs/alignment1.fasta"].decode()
    assert ">MonoDLG1-.2" in stored and stored.count(">") == 5


def test_alignments_with_different_sequences(client: Client, fake_storage: FakeStorage) -> None:
    fewer = "".join(f">seq{i}\nMKTAYIAK{'LV'[i % 2]}Q-\n" for i in range(1, 5))  # no seq0
    five = ALIGNED_FOUR + ">seq4\nMKTAYIAKLQ-\n"
    page = re.sub(r"\s+", " ", client.post("/alignment-files/", {
        "alignment_files": [upload("all.fasta", five), upload("fewer.fasta", fewer)],
    }).content.decode())  # fmt: skip
    assert "These don't all hold the same sequences" in page
    assert "fewer.fasta: 4 sequences, without seq0" in page
    assert "Use the 4 sequences they all share (leaves out 1)" in page
    data = {"aligners": [], "run_trees": "", "run_ancestors": ""}
    response = post(client, **data, alignment_files=[upload("all.fasta", five),
                                                     upload("fewer.fasta", fewer)])  # fmt: skip
    assert response.status_code == 200 and "names do not match" in response.content.decode()
    response = post(
        client,
        **data,
        shared_names="on",
        alignment_files=[upload("all.fasta", five), upload("fewer.fasta", fewer)],
    )
    assert response.status_code == 302, response.content.decode()[-2000:]
    job = Job.objects.get()
    assert job.params["left_out"] == ["seq0"] and job.params["left_out_count"] == 1
    stored = fake_storage.objects[f"jobs/{job.id}/inputs/alignment1.fasta"].decode()
    assert stored.split()[::2] == [">seq1", ">seq2", ">seq3", ">seq4"]


def test_failed_job_page_says_why(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(
        owner=local_user(), name="x", status=Job.Status.FAILED,
        error="The pipeline failed. The log below may help; if not, please report it.",
        progress={"tasks": [{"name": "PHYLOBOT:VALIDATE (1)", "status": "COMPLETED"},
                            {"name": "PHYLOBOT:TREE_IQTREE (mafft/LG+G4)", "status": "FAILED",
                             "exit": "137"}]},
    )  # fmt: skip
    fake_storage.objects[f"{job.storage_prefix}results/logs/tasks/TREE_IQTREE_mafft_LG+G4.err"] = (
        b"std::bad_alloc\n"
    )
    page = re.sub(r"\s+", " ", client.get(f"/jobs/{job.id}/").content.decode())
    assert "Why It Failed: The step ran out of memory" in page
    assert "Step: <strong>IQ-TREE search (mafft/LG+G4)</strong> · your input had passed" in page
    assert "<pre>std::bad_alloc" in page and "own message" not in page


def test_failed_job_can_be_resumed(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(
        owner=local_user(), status=Job.Status.FAILED, error=executors.INTERRUPTED,
        started_at=timezone.now(), progress={"done": 3, "total": 8}, resume_cache=True,
    )  # fmt: skip
    JobAttempt.objects.create(
        job=job, number=1, started_at=job.started_at, status="failed", tasks_run=4,
        completed_cpu_seconds=50, failed_cpu_seconds=10, failed_rate_percent=50,
        billable_cpu_seconds=55,
    )  # fmt: skip
    page = re.sub(r"\s+", " ", client.get(f"/jobs/{job.id}/").content.decode())
    assert "Resume this job" in page and "(3 of 8) are reused" in page
    failure_box = page[page.index('<section class="job-failure') :]
    assert "Resume this job" in failure_box[: failure_box.index("</section>")]
    assert page.count("Resume this job") == 1  # inside "Why It Failed", not also beneath
    assert "Attempts" not in page  # only once there is more than one
    assert client.get(f"/jobs/{job.id}/resume/").status_code == 405
    response = client.post(f"/jobs/{job.id}/resume/")
    assert response.status_code == 302
    job.refresh_from_db()
    assert job.status == Job.Status.QUEUED and job.attempt == 2
    page = re.sub(r"\s+", " ", client.get(f"/jobs/{job.id}/").content.decode())
    assert "attempt 2" in page and "Resume this job" not in page
    assert "<h2>Attempts</h2>" in page and "(failed steps at 50%)" in page
    Job.objects.filter(id=job.id).update(
        status=Job.Status.FAILED, resume_cache=False, error=executors.INTERRUPTED
    )
    page = re.sub(r"\s+", " ", client.get(f"/jobs/{job.id}/").content.decode())
    assert "starts again from the beginning, as attempt 3" in page  # an older job
    Job.objects.filter(id=job.id).update(work_expired_at=datetime(2026, 10, 1, tzinfo=UTC))
    page = re.sub(r"\s+", " ", client.get(f"/jobs/{job.id}/").content.decode())
    assert "removed from the worker on 2026-10-01" in page
    theirs = other_users_job()
    Job.objects.filter(id=theirs.id).update(status=Job.Status.FAILED)
    assert client.post(f"/jobs/{theirs.id}/resume/").status_code == 404
    theirs.refresh_from_db()
    assert theirs.status == Job.Status.FAILED and theirs.attempt == 1


def other_users_job() -> Job:
    return Job.objects.create(owner=User.objects.create_user("bob"))


def test_jobs_are_private(client: Client, fake_storage: FakeStorage) -> None:
    mine = Job.objects.create(owner=local_user(), name="mine")
    theirs = other_users_job()
    assert client.get(f"/jobs/{mine.id}/").status_code == 200
    for url in (f"/jobs/{theirs.id}/", f"/jobs/{theirs.id}/status/",
                f"/jobs/{theirs.id}/files/summary.json"):  # fmt: skip
        assert client.get(url).status_code == 404
    listing = client.get("/jobs/").content.decode()
    assert "mine" in listing and str(theirs.id)[:8] not in listing


def test_status_polls_until_finished(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user())
    assert jobs.mark(job.id, Job.Status.RUNNING)
    jobs.heartbeat(job.id, {"done": 2, "total": 6, "last": "ALIGN_MAFFT (mafft)", "tasks": []})
    prefix = f"jobs/{job.id}/results/"
    fake_storage.objects[prefix + "logs/pipeline.log"] = b"N E X T F L O W\nstep two\n"
    text = client.get(f"/jobs/{job.id}/status/").content.decode()
    assert 'hx-trigger="every 5s"' in text and "2 of 6 steps" in text and "step two" in text

    combo = {"key": "mafft/iqtree/modelfinder", "aligner": "mafft", "tree": "iqtree",
             "model": "modelfinder"}  # fmt: skip
    summary = {"primary": combo["key"], "combinations": [combo]}
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    fake_storage.objects[prefix + "mafft/iqtree/modelfinder/ancestors.fasta"] = b">Node1\nMKT\n"
    assert jobs.mark(job.id, Job.Status.SUCCEEDED)
    text = client.get(f"/jobs/{job.id}/status/").content.decode()
    assert "hx-trigger" not in text and "Alignment: mafft" in text
    assert f"/jobs/{job.id}/files/mafft/iqtree/modelfinder/ancestors.fasta" in text


def test_failed_job_shows_only_the_user_message(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user())
    assert jobs.mark(job.id, Job.Status.FAILED, "Input error: need at least 4 sequences")
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Input error: need at least 4 sequences" in page


def test_file_downloads(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user())
    prefix = f"jobs/{job.id}/results/"
    fake_storage.objects[prefix + "summary.json"] = b"{}"
    fake_storage.objects[prefix + "pipeline_info/report.html"] = b"<script>alert(1)</script>"
    fake_storage.objects[f"jobs/{job.id}/inputs/sequences.fasta"] = b">x\nMK\n"
    text = client.get(f"/jobs/{job.id}/files/summary.json")
    assert text["Content-Type"].startswith("text/plain")
    assert text["Content-Disposition"].startswith("inline")
    html = client.get(f"/jobs/{job.id}/files/pipeline_info/report.html")
    assert html["Content-Type"] == "application/octet-stream"
    assert html["Content-Disposition"].startswith("attachment")
    assert html["Content-Security-Policy"] == "sandbox"
    for bad in ("../inputs/sequences.fasta", "missing.txt", "%2e%2e/inputs/sequences.fasta"):
        assert client.get(f"/jobs/{job.id}/files/{bad}").status_code == 404


def test_file_tables_describe_every_file() -> None:
    from phylobot import filelist

    combo = "mafft/iqtree/modelfinder_codon"
    keys = [  # every kind of file a job writes (protein, codon, gaps, codeml, logs)
        "summary.json", "model_comparison.tsv", "validation.json", "versions.tsv",
        "sources.tsv", "mafft/alignment.fasta", "mafft/alignment.cds.fasta",
        *(f"{combo}/{name}" for name in (
            "alignment.fasta", "alignment.cds.fasta", "tree.nwk", "tree.rooted.nwk",
            "rooting.json", "supports.tsv", "site_rates.tsv", "site_lh.tsv", "model.txt",
            "fit.tsv", "modelfinder.tsv",
            "ancestors.fasta", "ancestors_pp.tsv", "ancestors_posteriors.tsv",
            "ancestors.cds.fasta", "ancestors_codon_pp.tsv", "ancestors_codon_posteriors.tsv",
            "gaps-ml.tsv", "grasp_model.txt", "ancestors.gaps-parsimony.fasta",
            "ancestors.gaps-fitch.cds.fasta", "ancestors.gaps-ml.fasta",
            "ancestors.gaps-prank.fasta", "ancestors.gaps-grasp.fasta",
            "codeml/codeml.tsv", "codeml/selection.tsv", "codeml/beb.tsv",
            "codeml/M8/ancestors.fasta", "codeml/M8/ancestors_codon_posteriors.tsv",
        )),
        "logs/nextflow.log", "logs/pipeline.log", "logs/params.json",
        "logs/tasks/CODEML_mafft_iqtree_modelfinder_codon_M0.start1.log",
        "logs/tasks/ALIGN_MAFFT_mafft.err", "pipeline_info/report.html",
        "pipeline_info/trace.tsv",
    ]  # fmt: skip
    label = "MAFFT × IQ-TREE × ModelFinder (codon)"
    summary = {"alignments": [{"key": "mafft", "label": "MAFFT"}],
               "combinations": [{"key": combo, "label": label}]}  # fmt: skip
    tables = filelist.file_tables(keys, summary)
    assert [t.title for t in tables] == [
        "Job", "Alignment: MAFFT", "MAFFT × IQ-TREE × ModelFinder (codon)", "Logs",
        "Pipeline report",
    ]  # fmt: skip
    undescribed = [key for t in tables for key, _, text in t.rows if not text]
    assert undescribed == []
    rows = {name: text for t in tables for _, name, text in t.rows}
    assert rows["codeml/M8/ancestors.fasta"].startswith("codeml M8: ancestral sequences")
    assert rows["ancestors.gaps-fitch.cds.fasta"].startswith("Coding ancestors with columns")
    assert rows["tasks/CODEML_mafft_iqtree_modelfinder_codon_M0.start1.log"] == (
        "Step “codeml” (mafft_iqtree_modelfinder_codon_M0.start1): its output."
    )
    assert (
        rows["tasks/ALIGN_MAFFT_mafft.err"] == "Step “Align: MAFFT” (mafft): its errors (stderr)."
    )
    assert sum(len(t.rows) for t in tables) == len(keys)


def test_eta_from_history() -> None:
    assert seconds("1m 3s") == 63 and seconds("45ms") == 0.045 and seconds("") is None
    owner = local_user()
    done = Job.objects.create(owner=owner, params={})
    tasks = [{"name": f"{kind} (x)", "status": "COMPLETED", "realtime": t} for kind, t in (
        ("VALIDATE", "1s"), ("ALIGN_MAFFT", "2s"), ("EXPORT_ALIGNMENT", "1s"),
        ("MODELFINDER", "10s"), ("TREE_IQTREE", "1m"), ("ANCESTRAL_IQTREE", "5s"),
        ("EXPORT", "1s"), ("SUMMARY", "1s"), ("ANCESTOR_MAP", "1s"))]  # fmt: skip
    Job.objects.filter(id=done.id).update(status=Job.Status.SUCCEEDED, progress={"tasks": tasks})
    running = Job.objects.create(owner=owner, params={}, progress={"tasks": tasks[:2]})
    assert eta.remaining_seconds(running) == 1 + 10 + 60 + 5 + 1 + 1 + 1
    fixed = Job.objects.create(owner=owner, params={"fixed_models": True, "models": ["LG"]})
    assert eta.remaining_seconds(fixed) == 1 + 2 + 1 + 10 + 2 * (60 + 5 + 1) + 1 + 1
    only_fixed = {"fixed_models": True, "models": ["LG"], "modelfinder": False}
    no_mf = Job.objects.create(owner=owner, params=only_fixed)
    assert eta.remaining_seconds(no_mf) == 1 + 2 + 1 + (60 + 5 + 1) + 1 + 1
    assert eta.describe(None).startswith("Estimate available")
    assert eta.describe(67) == "about 1 min"


@override_settings(PHYLOBOT_MODE="cloud", MIDDLEWARE=middleware_for("cloud"))
def test_cloud_mode_pages_need_login(client: Client) -> None:
    for url in ("/", "/jobs/new/", "/jobs/", "/preview/"):
        assert client.get(url).status_code == 302


def test_tool_lists_match_pipeline() -> None:
    from phylobot.forms import ALIGNERS, GAP_METHODS, RECONCILED_TREE, TREE_BUILDERS

    main = (ROOT / "pipeline" / "main.nf").read_text()
    from phylobot.forms import NEED_OUTGROUP
    from phylobot.gaps import NEED_OUTGROUP as PIPELINE_NEED_OUTGROUP

    assert NEED_OUTGROUP == set(PIPELINE_NEED_OUTGROUP)
    gap_list = re.search(r"'gap_methods', \[(.*?)\]\)", main)
    assert gap_list and [m for m, _ in GAP_METHODS] == re.findall(r"'(\w+)'", gap_list.group(1))
    aligners = re.search(r"parseTools\(params.aligners, 'aligners', \[(.*?)\]\)", main)
    builders = re.search(r"parseTools\(params.tree_builders, 'tree_builders', \[(.*?)\]\)", main)
    assert aligners and builders
    assert [a for a, _ in ALIGNERS] == re.findall(r"'(\w+)'", aligners.group(1))
    # An inherited alignment's key (phylobot.extend) is an aligner's, in both languages.
    from phylobot.pipeline.validate import ALIGNMENT_KEYS

    keys = re.search(r"key in \['provided', (.*?)\]", main)
    assert keys and list(ALIGNMENT_KEYS) == re.findall(r"'(\w+)'", keys.group(1))
    assert list(ALIGNMENT_KEYS) == [a for a, _ in ALIGNERS]
    # GeneRax's reconciled tree is offered in the reconciliation block, not beside the others.
    offered = [t for t, _ in TREE_BUILDERS] + [RECONCILED_TREE]
    assert offered == re.findall(r"'(\w+)'", builders.group(1))


MODELFINDER_TSV = "\n".join(
    "\t".join(row.split())
    for row in (
        "model log_likelihood aic aic_in_95 w_aic aicc aicc_in_95 w_aicc bic bic_in_95 w_bic",
        "LG+F+G4 -9414.478 18890.956 yes 0.328 18892.901 yes 0.382 19044.668 yes 0.999",
        "LG+F+R3 -9411.097 18890.195 yes 0.48 18892.535 yes 0.459 19058.782 no 0.000861",
        "WAG+G4 -9516.0 19120.0 no 0 19121.0 no 0 19200.0 no 0",
    )
)


def comparison_tsv() -> str:
    fields = (
        "aligner tree model_key model log_likelihood free_parameters "
        "aic aicc bic delta_bic w_bic w_aic w_aicc rank_bic"
    ).split()
    rows = [
        ["mafft", "iqtree", "modelfinder", "LG+F+G4", "-9414.4707", "31", "18890.9", "18892.9",
         "19044.6533", "0.0000", "1", "1", "1", "1"],
        ["mafft", "iqtree", "WAG+G4", "WAG+G4", "-9516.0431", "12", "19056.1", "19056.4",
         "19115.5877", "70.9344", "3.952e-16", "0", "0", "2"],
    ]  # fmt: skip
    return "\n".join("\t".join(r) for r in [fields, *rows]) + "\n"


def comparison_job(fake_storage: FakeStorage) -> Job:
    """A finished job: MAFFT with ModelFinder and WAG+G4, ranked in model_comparison.tsv."""
    job = Job.objects.create(owner=local_user(), name="Two models")
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combos = [{"key": f"mafft/iqtree/{m}", "aligner": "mafft", "tree": "iqtree", "model": m,
               "aligner_label": "MAFFT"} for m in ("modelfinder", "WAG+G4")]  # fmt: skip
    summary = {"primary": combos[0]["key"], "combinations": combos}
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    fake_storage.objects[prefix + "model_comparison.tsv"] = comparison_tsv().encode()
    return job


def test_results_table_downloads_as_csv_and_excel(
    client: Client, fake_storage: FakeStorage
) -> None:
    """Above each table: CSV and Excel of it, with the statistics unrounded and each
    row's combination, ranked as on the page."""
    import openpyxl

    job = comparison_job(fake_storage)
    page = client.get(f"/jobs/{job.id}/").content.decode()
    csv_url = f"/jobs/{job.id}/results-tables/alignment-mafft.csv"
    assert f'href="{csv_url}"' in page and f'href="{csv_url[:-4]}.xlsx"' in page
    response = client.get(csv_url)
    assert response.status_code == 200 and response["Content-Type"].startswith("text/csv")
    assert 'filename="Two_models-' in response["Content-Disposition"]
    assert response["Content-Disposition"].endswith('_mafft_models.csv"')
    rows = list(csv.reader(io.StringIO(response.content.decode())))
    assert rows[0][:5] == ["Rank", "Combination", "Model", "Fitted model", "log L"]
    assert rows[1][:5] == ["1", "mafft/iqtree/modelfinder", "ModelFinder", "LG+F+G4", "-9414.4707"]
    assert rows[2][1] == "mafft/iqtree/WAG+G4" and rows[2][-1] == "3.952e-16"  # BIC weight
    workbook = openpyxl.load_workbook(io.BytesIO(client.get(csv_url[:-4] + ".xlsx").content))
    sheet = workbook["Models"]
    values = [[c.value for c in row] for row in sheet.iter_rows()]
    assert values[0] == rows[0] and values[1][0] == 1 and values[1][4] == -9414.4707  # numbers
    assert workbook["About"]["A1"].value == "Job"


@pytest.mark.parametrize("path", ["alignment-muscle.csv", "alignment-mafft.json", "x.csv"])
def test_results_table_download_refuses_unknown(
    client: Client, fake_storage: FakeStorage, path: str
) -> None:
    job = comparison_job(fake_storage)
    assert client.get(f"/jobs/{job.id}/results-tables/{path}").status_code == 404


def test_results_table_per_alignment(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user())
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combos = [
        {"key": "mafft/iqtree/modelfinder", "aligner": "mafft", "tree": "iqtree",
         "model": "modelfinder", "aligner_label": "MAFFT", "tree_label": "IQ-TREE",
         "model_label": "ModelFinder", "fit": {"sites": "1052"}},
        {"key": "mafft/iqtree/WAG+G4", "aligner": "mafft", "tree": "iqtree", "model": "WAG+G4",
         "aligner_label": "MAFFT", "tree_label": "IQ-TREE", "model_label": "WAG+G4"},
    ]  # fmt: skip
    summary = {"primary": "mafft/iqtree/modelfinder", "combinations": combos}
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    fake_storage.objects[prefix + "model_comparison.tsv"] = comparison_tsv().encode()
    fake_storage.objects[prefix + "mafft/iqtree/modelfinder/modelfinder.tsv"] = (
        MODELFINDER_TSV.encode()
    )
    fake_storage.objects[prefix + "mafft/iqtree/modelfinder/alignment.fasta"] = (
        b">a\nMK-T\n>b\nM-TT\n"  # 2 gaps in 8 cells = 25.0%
    )
    page = client.get(f"/jobs/{job.id}/").content.decode()
    results = page[page.index(">Results</button></h2>") :]
    assert results.count(">Alignment: MAFFT</button></h3>") == 1  # one table for the alignment
    assert "mafft/iqtree/modelfinder/alignment.fasta" in results
    assert "4 columns · 2 sequences · 25.0% gaps" in re.sub(r"\s+", " ", results)
    start = results.index('id="alignments-overview">Alignments</h3>')
    overview = results[start : results.index("</table>")]
    cells = re.findall(r"<td[^>]*>(.*?)</td>", overview, re.S)
    assert [re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", c)).strip() for c in cells] == [
        "MAFFT", "4", "2", "25.0%", "view · alignment.fasta"
    ]  # fmt: skip
    # one row per combination, ranked; ModelFinder shows its pick; no "primary" badge
    table = results[results.index(">Alignment: MAFFT</button></h3>") :]
    table = table[: table.index("</table>")]
    assert table.count("<tr>") == 3  # header + 2 combinations
    assert "ModelFinder\u00a0→ LG+F+G4" in table and "primary" not in table
    assert table.index("ModelFinder\u00a0→") < table.index("WAG+G4")
    for link in ("mafft/iqtree/WAG%2BG4/tree.nwk", "mafft/iqtree/WAG+G4/tree.nwk"):
        if link in table:
            break
    else:
        raise AssertionError("tree.nwk link missing")
    assert "ancestors.zip?combination=mafft%2Fiqtree%2Fmodelfinder" in table
    assert "Tree source" not in table  # only one tree source
    best, worst = "rgb(161, 201, 225)", "rgb(234, 176, 166)"
    assert best in table and worst in table
    # ModelFinder's candidates: a collapsed secondary table under the alignment
    assert '<details class="candidates">' in results
    assert "ModelFinder candidates (3 models)" in results
    assert "LG+F+G4 ✓" in results
    # the old sections are gone
    assert "Model Comparisons" not in page and "Ancestral sequences" not in page


def test_single_combination_still_gets_a_table(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user())
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combo = {"key": "mafft/iqtree/modelfinder", "aligner": "mafft", "tree": "iqtree",
             "model": "modelfinder"}  # fmt: skip
    # (an older summary.json, without label fields)
    summary = {"primary": combo["key"], "combinations": [combo]}
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    fake_storage.objects[prefix + "model_comparison.tsv"] = comparison_tsv().encode()
    page = client.get(f"/jobs/{job.id}/").content.decode()
    table = page[page.index(">Alignment: mafft</button></h3>") :]
    table = table[: table.index("</table>")]
    assert table.count("<tr>") == 2 and "ModelFinder\u00a0→ LG+F+G4" in table
    assert "background-color" not in table  # nothing to compare a single row against


def test_delete_jobs_keeps_records_and_removes_files(
    client: Client, fake_storage: FakeStorage
) -> None:
    owner = local_user()
    done, failed, running = (Job.objects.create(owner=owner, name=n) for n in "abc")
    assert jobs.mark(done.id, Job.Status.RUNNING) and jobs.mark(done.id, Job.Status.SUCCEEDED)
    assert jobs.mark(failed.id, Job.Status.FAILED, "Input error: x")
    assert jobs.mark(running.id, Job.Status.RUNNING)
    theirs = other_users_job()
    for job in (done, failed, running, theirs):
        fake_storage.objects[f"jobs/{job.id}/results/summary.json"] = b"{}"

    response = client.post(
        "/jobs/delete/",
        {"job": [str(done.id), str(failed.id), str(running.id), str(theirs.id), "nonsense"]},
        follow=True,
    )
    assert "Deleted 3 jobs." in response.content.decode()
    for job in (done, failed, running, theirs):
        job.refresh_from_db()
    # Finished jobs: files gone now; records kept with their history.
    assert done.deleted_at and done.purged_at and done.status == Job.Status.SUCCEEDED
    assert failed.deleted_at and failed.purged_at and failed.error == "Input error: x"
    assert not fake_storage.list_keys(f"jobs/{done.id}/")
    # Running job: cancelled (its worker stops at the next heartbeat) and purged too.
    assert running.status == Job.Status.CANCELLED and running.deleted_at and running.purged_at
    assert not fake_storage.list_keys(f"jobs/{running.id}/")
    # Another user's job is untouched.
    assert theirs.deleted_at is None and fake_storage.list_keys(f"jobs/{theirs.id}/")
    # Deleted Jobs vanish from the user's pages but stay in the database.
    listing = client.get("/jobs/").content.decode()
    assert "No jobs yet" in listing
    assert client.get(f"/jobs/{done.id}/").status_code == 404
    assert Job.objects.filter(deleted_at__isnull=False).count() == 3
    # Nothing left for the worker's maintenance sweep to retry.
    assert jobs.purge_deleted() == 0


def test_delete_requires_post(client: Client) -> None:
    assert client.get("/jobs/delete/").status_code == 405


def test_my_jobs_has_selection_and_confirmation(client: Client) -> None:
    Job.objects.create(owner=local_user(), name="one")
    page = client.get("/jobs/").content.decode()
    assert 'id="select-all"' in page and 'name="job"' in page
    assert '<dialog id="confirm-delete"' in page and "phylobot/delete.js" in page


def test_parse_tags() -> None:
    from phylobot.models import parse_tags

    assert parse_tags("#Ime2, review  #draft #ime2") == ["ime2", "review", "draft"]
    assert parse_tags("  ") == []
    for bad in ("#has space?", "#" + "x" * 41, "#-leading", "#émoji"):
        with pytest.raises(ValueError, match="not a valid tag"):
            parse_tags(bad)


def test_bulk_tags_on_own_jobs_only(client: Client) -> None:
    a, b = (Job.objects.create(owner=local_user(), name=n) for n in "ab")
    theirs = other_users_job()
    ids = [str(a.id), str(b.id), str(theirs.id)]
    response = client.post(
        "/jobs/tags/", {"action": "add", "tags": "#Ime2 review", "job": ids}, follow=True
    )
    assert "Added #ime2 #review to 2 jobs." in response.content.decode()
    assert sorted(a.tags.values_list("name", flat=True)) == ["ime2", "review"]
    assert not theirs.tags.exists()
    client.post("/jobs/tags/", {"action": "remove", "tags": "review", "job": [str(a.id)]})
    assert list(a.tags.values_list("name", flat=True)) == ["ime2"]
    assert list(b.tags.values_list("name", flat=True)) == ["ime2", "review"]


def test_bulk_tags_errors_and_safe_redirect(client: Client) -> None:
    job = Job.objects.create(owner=local_user())
    bad = client.post("/jobs/tags/", {"action": "add", "tags": "#no way!", "job": [str(job.id)]},
                      follow=True).content.decode()  # fmt: skip
    assert "not a valid tag" in bad
    none = client.post("/jobs/tags/", {"action": "add", "tags": "x"}, follow=True)
    assert "Select one or more jobs" in none.content.decode()
    outside = client.post("/jobs/tags/", {"action": "add", "tags": "x", "job": [str(job.id)],
                                          "next": "https://evil.example/"})  # fmt: skip
    assert outside["Location"] == "/jobs/"
    kept = client.post("/jobs/tags/", {"action": "add", "tags": "x", "job": [str(job.id)],
                                       "next": "/jobs/?tag=x"})  # fmt: skip
    assert kept["Location"] == "/jobs/?tag=x"


def test_tag_filter_and_chips(client: Client) -> None:
    tagged, plain = (Job.objects.create(owner=local_user(), name=n) for n in ("tagged", "plain"))
    tagged.add_tags(["ime2"])
    page = client.get("/jobs/?tag=%23IME2").content.decode()
    assert ">tagged<" in page and ">plain<" not in page and "Showing 1 of 2 jobs" in page
    assert '<option value="ime2" selected>#ime2</option>' in page
    assert 'href="/jobs/?tag=ime2"' in client.get("/jobs/").content.decode()
    assert "No jobs match these filters" in client.get("/jobs/?tag=nothing").content.decode()


def test_tags_are_per_user(client: Client) -> None:
    mine = Job.objects.create(owner=local_user())
    theirs = other_users_job()
    mine.add_tags(["shared-name"])
    theirs.add_tags(["shared-name"])
    assert mine.tags.get().pk != theirs.tags.get().pk
    assert '<option value="#shared-name">' in client.get("/jobs/").content.decode()
    assert client.get("/jobs/?tag=shared-name").content.decode().count("data-delete=") == 1


def test_job_page_tag_editing(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user())
    htmx = {"HX-Request": "true"}
    fragment = client.post(f"/jobs/{job.id}/tags/", {"action": "add", "tags": "#A b"}, headers=htmx)
    assert fragment.status_code == 200 and 'id="job-tags"' in fragment.content.decode()
    assert "#a" in fragment.content.decode() and "<html" not in fragment.content.decode()
    client.post(f"/jobs/{job.id}/tags/", {"action": "remove", "tags": "a"}, headers=htmx)
    assert list(job.tags.values_list("name", flat=True)) == ["b"]
    bad = {"action": "add", "tags": "#bad tag!"}
    error = client.post(f"/jobs/{job.id}/tags/", bad, headers=htmx)
    assert "not a valid tag" in error.content.decode()
    theirs = other_users_job()
    assert client.post(f"/jobs/{theirs.id}/tags/", {"tags": "x"}).status_code == 404
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert 'id="job-tags"' in page and "#b" in page and 'id="delete-job"' in page


def test_delete_from_job_page(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user(), name="solo")
    assert jobs.mark(job.id, Job.Status.FAILED)
    response = client.post("/jobs/delete/", {"job": [str(job.id)]}, follow=True)
    assert response.redirect_chain[-1][0] == "/jobs/"
    assert "Deleted 1 job." in response.content.decode()


def test_deleted_jobs_table(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user(), name="old run")
    assert jobs.mark(job.id, Job.Status.FAILED)
    job.add_tags(["archive"])
    Job.objects.filter(id=job.id).update(cpu_seconds=125)
    assert "Deleted Jobs" not in client.get("/jobs/").content.decode()
    jobs.delete([job])
    other_deleted = other_users_job()
    jobs.delete([other_deleted])
    page = client.get("/jobs/").content.decode()
    deleted = page[page.index("Deleted Jobs") :]
    assert "old run" in deleted and "#archive" in deleted and "2 min 5 s" in deleted
    assert f'href="/jobs/{job.id}/"' not in page  # records only: no link to a gone page
    assert str(other_deleted.id)[:8] not in page


def test_my_jobs_shows_eta_for_running_jobs(client: Client) -> None:
    owner = local_user()
    tasks = [{"name": f"{kind} (x)", "status": "COMPLETED", "realtime": t} for kind, t in (
        ("VALIDATE", "1s"), ("ALIGN_MAFFT", "2s"), ("EXPORT_ALIGNMENT", "1s"),
        ("MODELFINDER", "10s"), ("TREE_IQTREE", "3m"), ("ANCESTRAL_IQTREE", "5s"),
        ("EXPORT", "1s"), ("SUMMARY", "1s"), ("ANCESTOR_MAP", "1s"))]  # fmt: skip
    history = Job.objects.create(owner=owner, name="history")
    Job.objects.filter(id=history.id).update(status=Job.Status.SUCCEEDED, progress={"tasks": tasks})
    running = Job.objects.create(owner=owner, name="running", progress={"tasks": tasks[:2]})
    assert jobs.mark(running.id, Job.Status.RUNNING)
    queued = Job.objects.create(owner=owner, name="queued")
    page = client.get("/jobs/").content.decode()
    assert '<em class="eta">ETA: about 3 min</em>' in page
    assert page.count('class="eta"') == 1  # not for queued or finished jobs
    assert queued.status == Job.Status.QUEUED


def test_my_jobs_eta_before_any_history(client: Client) -> None:
    running = Job.objects.create(owner=local_user(), name="first ever")
    assert jobs.mark(running.id, Job.Status.RUNNING)
    assert '<em class="eta">ETA: estimating…</em>' in client.get("/jobs/").content.decode()


def test_job_page_shows_input_as_submitted(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user())
    raw = "> seq0  first\nmktayiak\n>seq1\nMKTAYIAR\n"  # odd spacing, lowercase: kept
    fake_storage.objects[f"jobs/{job.id}/inputs/sequences.fasta"] = raw.encode()
    report = {"sequences": 2, "min_length": 8, "max_length": 8}
    fake_storage.objects[f"jobs/{job.id}/results/validation.json"] = json.dumps(report).encode()
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Input Sequences" in page and "2 sequences" in page and "8 residues" in page
    assert "&gt; seq0  first\nmktayiak" in page  # verbatim (HTML-escaped)
    download = client.get(f"/jobs/{job.id}/inputs/sequences.fasta")
    assert b"".join(download.streaming_content) == raw.encode()  # type: ignore[attr-defined]
    assert download["Content-Type"].startswith("text/plain")
    for bad in ("../results/validation.json", "missing.fasta"):
        assert client.get(f"/jobs/{job.id}/inputs/{bad}").status_code == 404
    theirs = other_users_job()
    fake_storage.objects[f"jobs/{theirs.id}/inputs/sequences.fasta"] = b">x\nMK\n"
    assert client.get(f"/jobs/{theirs.id}/inputs/sequences.fasta").status_code == 404


def test_job_page_lists_the_outgroup(client: Client, fake_storage: FakeStorage) -> None:
    rooted = Job.objects.create(owner=local_user(), params={"outgroup": ["seq0", "<b>x</b>"]})
    unrooted = Job.objects.create(owner=local_user())
    for job in (rooted, unrooted):
        fake_storage.objects[f"jobs/{job.id}/inputs/sequences.fasta"] = b">seq0\nMK\n"
    page = client.get(f"/jobs/{rooted.id}/").content.decode()
    assert "2 in the outgroup" in page
    assert "<strong>Outgroup:</strong>\n      seq0, &lt;b&gt;x&lt;/b&gt;" in page  # escaped
    page = client.get(f"/jobs/{unrooted.id}/").content.decode()
    assert "in the outgroup" not in page and "none given" in page


def test_large_input_offered_as_download_only(client: Client, fake_storage: FakeStorage) -> None:
    from phylobot.views import INPUT_PREVIEW_BYTES

    job = Job.objects.create(owner=local_user())
    big = (">s\n" + "M" * INPUT_PREVIEW_BYTES).encode()
    fake_storage.objects[f"jobs/{job.id}/inputs/sequences.fasta"] = big
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Too large to show here" in page and "M" * 1000 not in page


def test_every_alignment_and_combination_listed(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user())
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    labels = {"mafft": "MAFFT", "prank": "PRANK"}
    combos = [
        {"key": f"{a}/{t}/{m}", "aligner": a, "tree": t, "model": m, "aligner_label": labels[a],
         "tree_label": {"iqtree": "IQ-TREE", "provided1": "My tree"}[t],
         "model_label": "ModelFinder" if m == "modelfinder" else m}
        for a in ("prank", "mafft") for t in ("iqtree", "provided1") for m in ("modelfinder", "JTT")
    ]  # fmt: skip
    summary = {"primary": "mafft/iqtree/modelfinder", "combinations": combos}
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    page = client.get(f"/jobs/{job.id}/").content.decode()
    results = page[page.index(">Results</button></h2>") :]
    # an overview table first, then one table per alignment, the primary's first
    start = results.index('id="alignments-overview">Alignments</h3>')
    overview = results[start : results.index("</table>")]
    assert overview.index(">MAFFT</a>") < overview.index(">PRANK</a>")
    assert 'href="#alignment-mafft"' in overview and 'id="alignment-mafft"' in results
    for aligner in ("mafft", "prank"):
        assert f"{aligner}/iqtree/modelfinder/alignment.fasta" in overview
    assert results.index("Alignment: MAFFT") < results.index("Alignment: PRANK")
    assert results[: results.index("All Files")].count("<table") == 1 + 2
    for combo in combos:  # every combination has its tree and ancestors linked
        assert f"ancestors.zip?combination={quote(combo['key'], safe='')}" in results
    # rows differ by tree too, so a Tree source column appears
    assert "<th>Tree source</th>" in results and "My tree" in results
    # no statistics yet (older job): rows still listed, with dashes
    assert "<td>–</td>" in results


def test_job_page_contents_and_folding(client: Client, fake_storage: FakeStorage) -> None:
    """An "On this page" nav lists each alignment's table; Results and each alignment fold;
    a table of more than SCROLL_ROWS rows scrolls in its own box."""
    from phylobot.results import SCROLL_ROWS

    job = Job.objects.create(owner=local_user())
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    many = [f"M{i}" for i in range(SCROLL_ROWS + 1)]
    combos = [
        {"key": f"{a}/iqtree/{m}", "aligner": a, "tree": "iqtree", "model": m}
        for a, models in (("mafft", many), ("muscle", ["LG"])) for m in models
    ]  # fmt: skip
    summary = {"primary": combos[0]["key"], "combinations": combos}
    fake_storage.objects[f"jobs/{job.id}/results/summary.json"] = json.dumps(summary).encode()
    page = client.get(f"/jobs/{job.id}/").content.decode()
    start = page.index('<nav class="viewer-toc page-toc job-toc"')
    nav = page[start : page.index("</nav>", start)]
    for target in ("#job-stages", "#results", "#alignments-overview", "#alignment-mafft",
                   "#alignment-muscle", "#files"):  # fmt: skip
        assert f'href="{target}"' in nav and f'id="{target[1:]}"' in page
    assert 'aria-controls="results-body"' in page and 'id="results-body"' in page
    assert 'aria-controls="alignment-mafft-body"' in page and "data-fold-all" not in page
    assert page.count('class="results-scroll"') == 1  # mafft's table only
    assert f"{SCROLL_ROWS + 1} rows: scroll the table for more." in page
    # the refreshed status block carries the nav's list out of band; the page itself not
    assert 'hx-swap-oob="true"' not in page
    status = client.get(f"/jobs/{job.id}/status/").content.decode()
    assert '<ol id="job-toc-live" hx-swap-oob="true">' in status and "#alignment-muscle" in status


def test_results_offer_each_rows_ancestors_as_one_zip(
    client: Client, fake_storage: FakeStorage
) -> None:
    """The Ancestors column is "view" and "download" (a ZIP of the row's ancestor files),
    with the ZIP's contents in its tooltip; the gap versions are those that row has."""
    job = Job.objects.create(owner=local_user())
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combos = [{"key": f"mafft/iqtree/{m}", "aligner": "mafft", "tree": "iqtree", "model": m}
              for m in ("modelfinder", "LG")]  # fmt: skip
    fake_storage.objects[prefix + "summary.json"] = json.dumps(
        {"primary": combos[0]["key"], "combinations": combos}
    ).encode()
    for name in ("ancestors.gaps-parsimony.fasta", "ancestors.gaps-ml.fasta", "gaps-ml.tsv"):
        fake_storage.objects[f"{prefix}mafft/iqtree/modelfinder/{name}"] = b">N\nM\n"
    page = client.get(f"/jobs/{job.id}/").content.decode()
    table = page[page.index(">Alignment: mafft</button></h3>") :]
    first, second = table.split("</tr>")[1], table.split("</tr>")[2]
    zipped = f"/jobs/{job.id}/ancestors.zip?combination=mafft%2Fiqtree%2Fmodelfinder"
    # as icons, their words the tooltips
    assert f'<a class="icon-link" href="{zipped}" title="Download" aria-label="Download">' in first
    assert 'title="View" aria-label="View"' in first
    assert (
        'title="ZIP of original, shared presence, ML, site probabilities, gap probabilities"'
        in first
    )
    assert 'title="ZIP of original, site probabilities"' in second  # no gap versions run
    assert ">shared presence</a>" not in table  # one link per file no longer


def test_ancestors_zip_holds_the_rows_ancestor_files(
    client: Client, fake_storage: FakeStorage
) -> None:
    import zipfile

    from phylobot import ancestor_names

    job = viewer_job(fake_storage)
    key = "mafft/iqtree/WAG"
    prefix = f"jobs/{job.id}/results/{key}/"
    fake_storage.objects[prefix + "ancestors_posteriors.index.tsv"] = b"node\tstart\tlength\n"
    fake_storage.objects[prefix + "properties/index.json"] = b"{}"
    fake_storage.objects[prefix + "properties/summaries.json"] = b"{}"
    response = client.get(f"/jobs/{job.id}/ancestors.zip", {"combination": key})
    assert response.status_code == 200 and response["Content-Type"] == "application/zip"
    assert "_mafft_iqtree_WAG_ancestors.zip" in response["Content-Disposition"]
    archive = zipfile.ZipFile(io.BytesIO(response.getvalue()))
    [folder] = {name.split("/")[0] for name in archive.namelist()}
    assert sorted(name[len(folder) + 1 :] for name in archive.namelist()) == [
        "ancestors.fasta", "ancestors.gaps-fitch.fasta", "ancestors.gaps-ml.fasta",
        "ancestors_posteriors.tsv", "ancestors_pp.tsv", "gaps-ml.tsv",
        "properties/index.json", "properties/summaries.json",
    ]  # fmt: skip
    assert archive.read(f"{folder}/ancestors.fasta") == b">Node1\nMKK\n>Node2\nMRK\n"
    # Names live in PhyloBot, not in the files: the ZIP carries them beside the files.
    mapped = fake_storage.objects[f"jobs/{job.id}/results/{analogs.MAP}"].decode()
    [row] = [r for r in analogs.read_map(mapped) if (r.combination, r.node) == (key, "Node2")]
    parts = str(row.place.parts)  # no outgroup: names are keyed by the subtrees
    assert ancestor_names.save(job, parts, parts, "Anc.AB") == ""
    response = client.get(f"/jobs/{job.id}/ancestors.zip", {"combination": key})
    archive = zipfile.ZipFile(io.BytesIO(response.getvalue()))
    names = archive.read(f"{folder}/ancestor_names.tsv").decode()
    assert names == "node\tname\tsame_branching\nNode2\tAnc.AB\tyes\n"


@pytest.mark.parametrize("combination", ["", "mafft/iqtree", "../x", "mafft/iqtree/nothere"])
def test_ancestors_zip_refuses_what_is_not_a_combination(
    client: Client, fake_storage: FakeStorage, combination: str
) -> None:
    job = viewer_job(fake_storage)
    assert (
        client.get(f"/jobs/{job.id}/ancestors.zip", {"combination": combination}).status_code == 404
    )


def test_ancestors_zip_is_the_owners_only(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    job.owner = User.objects.create(username="someone-else")
    job.save()
    query = {"combination": "mafft/iqtree/WAG"}
    assert client.get(f"/jobs/{job.id}/ancestors.zip", query).status_code == 404


def test_parsimony_gap_placement_is_ticked_by_default(client: Client) -> None:
    page = client.get("/jobs/new/").content.decode()
    box = re.search(r'<input[^>]*name="gap_methods"[^>]*value="parsimony"[^>]*>', page)
    assert box and "checked" in box.group(0)
    ml = re.search(r'<input[^>]*name="gap_methods"[^>]*value="ml"[^>]*>', page)
    assert ml and "checked" not in ml.group(0)


SIX = "".join(f">seq{i}\nMKTAYIAK{'LV'[i % 2]}Q\n" for i in range(6))


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        ({"gap_methods": ["prank"]}, "Gap placement by PRANK needs an outgroup."),
        ({"gap_methods": ["grasp", "prank"]}, "by GRASP and PRANK needs an outgroup."),
        ({"outgroup": ["nope"]}, "outgroup sequences not among the sequences: &#x27;nope&#x27;"),
        ({"outgroup": ["seq0", "seq1", "seq2", "seq3"]}, "at least 3 ingroup sequences"),
    ],
)
def test_outgroup_rules(
    client: Client, fake_storage: FakeStorage, fields: dict[str, Any], message: str
) -> None:
    page = post(client, sequences_text=SIX, **fields).content.decode()
    assert message in page and not Job.objects.exists()


def test_outgroup_saved_with_prank(client: Client, fake_storage: FakeStorage) -> None:
    response = post(client, sequences_text=SIX, outgroup=["seq0"], gap_methods=["prank"])
    assert response.status_code == 302
    job = Job.objects.get()
    assert job.params["outgroup"] == ["seq0"] and job.params["gap_methods"] == ["prank"]


def test_outgroup_checklist_follows_the_sequences(client: Client) -> None:
    fragment = client.post(
        "/outgroup-options/", {"sequences_text": SIX, "outgroup": ["seq2", "gone"]}
    ).content.decode()
    assert fragment.count('name="outgroup"') == 6
    assert re.search(r'value="seq2"\s+checked', fragment)
    assert "gone" not in fragment  # stale choices are dropped
    empty = client.post("/outgroup-options/", {"sequences_text": ""}).content.decode()
    assert "Enter or upload sequences above" in empty


def test_example_ticks_its_outgroup(client: Client) -> None:
    fragment = client.get("/examples/ime2/").content.decode()
    assert 'id="outgroup-section" hx-swap-oob="true"' in fragment
    assert re.search(r'value="Yarrowia.lipolytica.IME2"\s+checked', fragment)
    assert fragment.count('name="outgroup"') == 7


def test_rooted_tree_link_and_warning(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user(), params={"outgroup": ["a", "b"]})
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combos = [{"key": f"mafft/iqtree/{m}", "aligner": "mafft", "tree": "iqtree", "model": m}
              for m in ("modelfinder", "LG")]  # fmt: skip
    fake_storage.objects[prefix + "summary.json"] = json.dumps(
        {"primary": combos[0]["key"], "combinations": combos}
    ).encode()
    for key, monophyletic in ((combos[0]["key"], True), (combos[1]["key"], False)):
        fake_storage.objects[f"{prefix}{key}/tree.rooted.nwk"] = b"(a,b);"
        report = {"outgroup": ["a", "b"], "monophyletic": monophyletic,
                  "intruders": [] if monophyletic else ["c"], "method": "outgroup"}  # fmt: skip
        fake_storage.objects[f"{prefix}{key}/rooting.json"] = json.dumps(report).encode()
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Trees are rooted on the outgroup" not in page
    table = page[page.index(">Alignment: mafft</button></h3>") :]
    first, second = table.split("</tr>")[1], table.split("</tr>")[2]
    assert ">rooted</a>" in first and "⚠" not in first
    assert "⚠" in second and "also inside: c" in second


def viewer_job(fake_storage: FakeStorage, gap_methods: list[str] | None = None) -> Job:
    """A finished job with two alignments (mafft: two models) and ancestors Node1-2."""
    job = Job.objects.create(owner=local_user(), params={"gap_methods": gap_methods or []})
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combos = [
        {"key": f"{a}/iqtree/{m}", "aligner": a, "tree": "iqtree", "model": m,
         "aligner_label": a.upper(), "fit": {"model": "LG+G4" if m == "modelfinder" else m},
         "rank_bic": rank}
        for a, m, rank in (("mafft", "modelfinder", 1), ("mafft", "WAG", 2),
                           ("muscle", "modelfinder", 1))
    ]  # fmt: skip
    fake_storage.objects[prefix + "summary.json"] = json.dumps(
        {"primary": combos[0]["key"], "combinations": combos}
    ).encode()
    for c in combos:
        files = {
            "alignment.fasta": ">a\nMK-\n>b\nM-K\n>c<i>\nMKK\n",
            "tree.nwk": "((a,b)Node2,c<i>)Node1;",
            "ancestors.fasta": ">Node1\nMKK\n>Node2\nMRK\n",
            "ancestors.gaps-fitch.fasta": ">Node1\nMK-\n>Node2\nM-K\n",
            "ancestors_posteriors.tsv": "node\tsite\tM\tK\tR\n"
            + "".join(f"Node1\t{s}\t0.9\t0.1\t0\nNode2\t{s}\t0.2\t0.3\t0.5\n" for s in (1, 2, 3)),
            "ancestors_pp.tsv": "node\tsite\tstate\tpp\nNode1\t1\tM\t0.9\nNode2\t1\tR\t0.5\n",
            "supports.tsv": "node\tufboot\nNode2\t88\n",
            "ancestors.gaps-ml.fasta": ">Node1\nMK-\n>Node2\nMRK\n",
            "gaps-ml.tsv": "node\tsite\tp_gap\n"
            + "".join(
                f"Node1\t{s}\t{p}\nNode2\t{s}\t0.0\n" for s, p in ((1, 0.0), (2, 0.3), (3, 0.8))
            ),
        }
        if c["aligner"] == "muscle":  # the same clades, numbered the other way round
            files["tree.nwk"] = "((a,b)Node1,c<i>)Node2;"
        for name, text in files.items():
            fake_storage.objects[f"{prefix}{c['key']}/{name}"] = text.encode()
    fake_storage.map_ancestors(job)
    return job


def test_viewers_open_on_the_best_ranked_combination(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = viewer_job(fake_storage)
    key = f"jobs/{job.id}/results/summary.json"
    summary = json.loads(fake_storage.objects[key])
    ranks = {"mafft/iqtree/modelfinder": 2, "mafft/iqtree/WAG": 1, "muscle/iqtree/modelfinder": 1}
    for combo in summary["combinations"]:
        combo["rank_bic"] = ranks[combo["key"]]
    fake_storage.objects[key] = json.dumps(summary).encode()
    for viewer in ("tree", "ancestors"):
        page = client.get(f"/jobs/{job.id}/{viewer}/").content.decode()
        assert '<option value="mafft/iqtree/WAG" selected>' in page


def test_ancestor_viewer_defaults_and_tiles(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    page = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
    assert '<option value="mafft/iqtree/modelfinder" selected>ModelFinder → LG+G4' in page
    assert '<option value="fitch" selected>Parsimony (Fitch)' in page  # the job's method
    assert '<option value="Node1" selected>Node1 · 3 sequences' in page
    assert "Node2 · 2 sequences · support 88" in page
    assert "Site in MAFFT" in page and '<option value="a" selected>a</option>' in page
    assert "c&lt;i&gt;" in page and "c<i>" not in page  # extant names are escaped
    assert page.count('class="tile gap"') == 1  # Node1 masked at column 3
    assert "M<br>(0.90)" in page and "K<br>(0.10)" in page
    assert "<strong>Node1</strong>" in page and "2 residues · " in page
    assert "sites table (Excel)" in page  # in Downloads too
    assert "Download (Node1)" in page
    # support plots over Node1's two unmasked sites (ML pp 0.9 each)
    assert "mean posterior probability 0.900" in page and "over 2 sites" in page
    assert 'class="pp-trace"' in page and 'id="pp-sites"' in page
    assert "Site in ancestor Node1" in page
    # the job page links here, per combination and overall
    job_page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Ancestor Viewer →" in job_page
    viewer_link = f'/jobs/{job.id}/ancestors/?combination=muscle%2Fiqtree%2Fmodelfinder"'
    assert f'{viewer_link} title="View"' in job_page


def test_ancestor_viewer_raw_by_default_calibrated_on_request(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/ancestors/"
    page = client.get(url).content.decode()
    assert '<option value="raw" selected>Raw (IQ-TREE\'s)</option>' in page
    assert "data-calibration" not in page and 'class="temperature"' not in page
    # The built-in coefficients ship the fitted global T = 1.24: calibrated flattens raw.
    page = client.get(url, {"pp": "calibrated"}).content.decode()
    assert '<option value="calibrated" selected>Calibrated</option>' in page
    assert "data-calibration" in page and "T is 1 everywhere" not in page
    assert "Fitted on a benchmark" in page  # the fitted provenance, not the placeholder note
    assert "M<br>(0.85)" in page and "K<br>(0.15)" in page  # 0.90/0.10 raw, flattened at T=1.24
    # Coefficients set on the Settings page: T = 1.5 at every site here.
    SiteSettings.objects.create(pk=1, calibration_intercept=0.4054651081)
    page = client.get(url, {"pp": "calibrated"}).content.decode()
    assert "M<br>(0.81)" in page and "K<br>(0.19)" in page  # flatter, M still first
    assert "this site&#x27;s defaults" in page and "T is 1 everywhere" not in page
    assert '<td class="temperature">1.50</td>' in page
    assert "pp=calibrated" in page  # the downloads follow the view
    # Downloads: the page's probabilities, with T per site; sequences unchanged.
    query = {"combination": "mafft/iqtree/modelfinder", "node": "Node1", "pp": "calibrated"}
    table = client.get(f"{url}download/", query | {"what": "posteriors"}).content.decode()
    header, first, *_ = table.splitlines()
    assert header == "node\tsite\tM\tK\tR\ttemperature"
    assert first.startswith("Node1\t1\t0.81") and first.endswith("\t1.5000")
    states = client.get(f"{url}download/", query | {"what": "states"}).content.decode()
    assert states.splitlines()[1].startswith("Node1\t1\tM\t0.81")
    raw = client.get(f"{url}download/", query | {"what": "posteriors", "pp": "raw"})
    assert raw.content.decode().splitlines()[1] == "Node1\t1\t0.9\t0.1\t0"
    fasta = client.get(f"{url}download/", query | {"what": "sequence"}).content.decode()
    assert fasta.splitlines()[1] == "MKK"
    # The user's own coefficients beat the site's: here no effect at all.
    client.post("/calibration/", {"action": "none"})
    page = client.get(url, {"pp": "calibrated"}).content.decode()
    assert "M<br>(0.90)" in page and "your own" in page and "T is 1 everywhere" in page


def test_ancestor_viewer_offers_calibration_only_with_posteriors(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = viewer_job(fake_storage)
    key = f"jobs/{job.id}/results/mafft/iqtree/modelfinder/ancestors.gaps-prank.fasta"
    fake_storage.objects[key] = b">Node1\nMK-\n>Node2\nM-K\n"
    page = client.get(f"/jobs/{job.id}/ancestors/", {"gaps": "prank", "pp": "calibrated"})
    text = page.content.decode()
    assert 'name="pp"' not in text and "data-calibration" not in text


def test_ancestor_viewer_sections_and_sidebar(client: Client, fake_storage: FakeStorage) -> None:
    """The viewer is laid out as discrete sections with an "On this page" nav (owner,
    2026-10-20). The nav lives inside the swapped form, so it can never offer a section
    the view does not have — which is the failure worth testing."""
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    page = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
    nav_start = page.index("viewer-toc")
    nav = page[nav_start : page.index("</nav>", nav_start)]
    targets = re.findall(r'href="#([\w-]+)"', nav)
    # Sections are collapsible, so they are <details class="viewer-section">.
    sections = re.findall(SECTION_RE, page)
    assert "choose" in targets and "downloads" in targets and "sites" in targets
    # Every link points at a section that is on the page, exactly once...
    assert set(targets) <= set(sections), set(targets) - set(sections)
    assert all(page.count(f'id="{t}"') == 1 for t in targets)
    # ...and every section is reachable from the nav.
    assert set(sections) <= set(targets), set(sections) - set(targets)


def add_prank(fake_storage: FakeStorage, job: Job) -> None:
    """PRANK's own reconstruction for every combination: residues, no posteriors."""
    prefix = f"jobs/{job.id}/results/"
    for key in {k.rsplit("/", 1)[0] for k in fake_storage.objects if k.startswith(prefix)}:
        if key.count("/") == 5:  # jobs/<id>/results/<aligner>/<tree>/<model>
            prank = b">Node1\nMK-\n>Node2\nM-K\n"
            fake_storage.objects[f"{key}/ancestors.gaps-prank.fasta"] = prank


def test_ancestor_viewer_sidebar_drops_absent_sections(
    client: Client, fake_storage: FakeStorage
) -> None:
    """PRANK's reconstruction reports no posteriors, so nothing can be drawn from it:
    neither the support plots nor the properties — measured or offered — apply."""
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    add_prank(fake_storage, job)
    page = client.get(f"/jobs/{job.id}/ancestors/", {"gaps": "prank"}).content.decode()
    assert '<option value="prank" selected>' in page  # it really is PRANK's version
    nav_start = page.index("viewer-toc")
    nav = page[nav_start : page.index("</nav>", nav_start)]
    assert "#properties" not in nav
    sections = re.findall(SECTION_RE, page)
    assert "properties" not in sections


#: The viewer's sections: collapsible, so <details class="viewer-section" id="...">.
SECTION_RE = r'<(?:details|section) class="[^"]*viewer-section[^"]*" id="([\w-]+)"'


def store_properties(fake_storage: FakeStorage, job: Job) -> None:
    """PROPERTIES' output for Node1 and Node2 of the first combination, fitch version.
    "charge" is a real registry key, so the panel finds its method note and its draws."""
    prefix = f"jobs/{job.id}/results/"
    summary: dict[str, Any] = {
        "node": "Node1", "gaps": "fitch", "draws": 1000, "seed": "1:Node1:fitch",
        "sampling": "per-column marginal posteriors, drawn independently", "weighted": False,
        "properties": {
            "charge": {"label": "Net charge", "units": "e", "level": 0.95, "mean": -3.2,
                       "sd": 4.0, "low": -11.0, "high": 4.6, "ml": -1.5,
                       "standardised_shift": 0.425},
        },
    }  # fmt: skip
    charge: dict[str, Any] = summary["properties"]["charge"]
    for node in ("Node1", "Node2"):
        body: dict[str, Any] = dict(summary, node=node)
        if node == "Node2":  # a different value, so the panel is visibly per node
            body["properties"] = {"charge": dict(charge, mean=7.5)}
        draws = [float(-6 + i % 9) for i in range(1, 61)]
        fake_storage.store_properties(
            prefix, "mafft/iqtree/modelfinder", node, "fitch", body, {"charge": draws}
        )


def test_ancestor_viewer_downloads_table(client: Client, fake_storage: FakeStorage) -> None:
    """Each download is a link, and what the file is for sits behind the list's "?"."""
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    page = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
    start = page.index('class="download-links"')
    links = page[start : page.index("</ul>", start)]
    # Every download this view offers is a link: all but the ML-gap-probabilities one,
    # which needs ml gap placement.
    # Plus the Sites table as Excel, the same link as under Sites.
    assert links.count("<li><a ") == len(views.ANCESTOR_DOWNLOADS) - 1 + 1
    assert "sites table (Excel)" in links
    # The notes: each link's tooltip, and all of them in the popover the "?" opens.
    notes = links[links.index('id="help-downloads"') :]
    assert "gaps removed" in notes  # the ML sequence's note
    assert "does not rest on the uncertain sites" in notes  # AltAll's


def test_every_download_explains_itself() -> None:
    """A download with no note renders a blank cell, which tells the reader nothing."""
    assert set(views.ANCESTOR_DOWNLOADS) == set(views.DOWNLOAD_NOTES)
    assert all(len(note) > 40 for note in views.DOWNLOAD_NOTES.values())


def test_ancestor_properties_csv(client: Client, fake_storage: FakeStorage) -> None:
    """The panel as shown, as a file: same numbers, plus the provenance to read it later."""
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    store_properties(fake_storage, job)
    url = f"/jobs/{job.id}/ancestors/properties.csv"
    response = client.get(url, {"node": "Node1", "gaps": "fitch"})
    assert response.status_code == 200
    assert response["Content-Type"].startswith("text/csv")
    assert "attachment" in response["Content-Disposition"]
    rows = list(csv.reader(io.StringIO(response.content.decode())))
    assert rows[0][:5] == ["property", "label", "units", "mean", "sd"]
    body = {row[0]: row for row in rows[1:]}
    assert body["charge"][3] == "-3.2"  # the mean the page shows
    assert body["charge"][rows[0].index("node")] == "Node1"
    assert body["charge"][rows[0].index("gaps")] == "fitch"
    assert body["charge"][rows[0].index("sampling")].startswith("per-column marginal")
    # A node with nothing measured is a 404, not an empty file that looks like a result.
    assert client.get(url, {"node": "Node1", "gaps": "original"}).status_code == 404
    # Another user's job stays private.
    assert client.get(f"/jobs/{other_users_job().id}/ancestors/properties.csv").status_code == 404


def test_properties_across_combinations(client: Client, fake_storage: FakeStorage) -> None:
    """The same ancestor under every combination that recovered it."""
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    store_properties(fake_storage, job)
    # A second combination measuring the same clade, with a visibly different centre.
    prefix = f"jobs/{job.id}/results/"
    first = f"{prefix}mafft/iqtree/modelfinder/properties/summaries.json"
    summary = json.loads(fake_storage.objects[first])["nodes"]["fitch"]["Node1"]
    summary["properties"]["charge"] = dict(summary["properties"]["charge"], mean=5.5)
    draws = [float(2 + i % 9) for i in range(1, 61)]
    fake_storage.store_properties(
        prefix, "mafft/iqtree/WAG", "Node1", "fitch", summary, {"charge": draws}
    )

    url = f"/jobs/{job.id}/properties/"
    page = client.get(url, {"node": "Node1", "gaps": "fitch"}).content.decode()
    assert "Properties Across Combinations" in page
    assert "Net charge" in page
    # Both combinations are offered in the sidebar, ticked.
    assert page.count('name="show"') == 2
    assert page.count("checked") >= 2
    assert "compare-overlay" in page and "viz-area" in page
    # The ancestor viewer links here.
    viewer = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
    assert f"/jobs/{job.id}/properties/" in viewer

    # Unticking one leaves it in the table but off the plot.
    one = client.get(
        url, {"node": "Node1", "gaps": "fitch", "show": "mafft/iqtree/WAG"}
    ).content.decode()
    assert one.count("<polyline") == 1  # one series left on the overlay
    assert 'class="off"' in one  # the hidden combination, still in the table

    # A combination whose tree lacks the clade is named, not silently dropped.
    assert "Not recovered under" in page or "muscle" in page


def test_properties_page_compares_the_best_few_per_alignment(
    client: Client, fake_storage: FakeStorage
) -> None:
    """Owner, 2026-10-30: the best N per alignment by BIC (3 by default, a menu to change
    it, a warning above 3), and the combination the viewer came from whatever its rank."""
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    store_properties(fake_storage, job)
    prefix = f"jobs/{job.id}/results/"
    stored = f"{prefix}mafft/iqtree/modelfinder/properties/summaries.json"
    summary = json.loads(fake_storage.objects[stored])
    body = summary["nodes"]["fitch"]["Node1"]
    fake_storage.store_properties(
        prefix, "mafft/iqtree/WAG", "Node1", "fitch", body, {"charge": [1.0, 2.0, 3.0]}
    )
    url = f"/jobs/{job.id}/properties/"
    base = {"node": "Node1", "gaps": "fitch", "combination": "mafft/iqtree/modelfinder"}
    # The default (3) covers both MAFFT combinations; the menu offers only what changes it.
    page = client.get(url, base).content.decode()
    assert page.count('name="show"') == 2
    assert '<option value="1">Best 1 by BIC</option>' in page
    assert '<option value="all" selected>All (2)</option>' in page  # 3 covers both
    assert 'id="top-warning"' in page and 'data-top-default="3"' in page
    # Best 1: MAFFT's rank-2 WAG is left out, and is neither read nor listed as missing.
    one = client.get(url, base | {"top": "1"}).content.decode()
    assert one.count('name="show"') == 1 and "mafft/iqtree/WAG" not in one.split("compare-main")[0]
    assert "Comparing 2 of this job's 3 combinations." in one  # MAFFT's best and MUSCLE's
    # The combination the viewer came from stays, even when it is not among the best.
    came = client.get(url, base | {"top": "1", "combination": "mafft/iqtree/WAG"}).content.decode()
    assert 'value="mafft/iqtree/WAG"' in came
    # Nonsense is the default, never an error; the choice travels in the shareable URL.
    assert client.get(url, base | {"top": "999"}).status_code == 200
    assert property_compare.top_of("999") == 3 and property_compare.top_of("all") is None
    assert '<option value="1" selected>' in one  # a GET form: the choice is in the URL


def test_properties_page_display_options(client: Client, fake_storage: FakeStorage) -> None:
    """The mean lines and the palette are in the URL, like every other viewer choice."""
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    store_properties(fake_storage, job)
    url = f"/jobs/{job.id}/properties/"
    base = {"node": "Node1", "gaps": "fitch"}
    page = client.get(url, base).content.decode()
    assert 'class="viz-mean' in page  # means on by default
    assert 'data-palette="bright"' in page
    off = client.get(url, base | {"means": "off"}).content.decode()
    assert "viz-mean" not in off
    okabe = client.get(url, base | {"palette": "okabe"}).content.decode()
    assert 'data-palette="okabe"' in okabe
    junk = client.get(url, base | {"palette": "../../etc"}).content.decode()
    assert 'data-palette="bright"' in junk


def test_property_help_buttons_have_their_script(client: Client, fake_storage: FakeStorage) -> None:
    """The (?) markup is shared by two pages but its behaviour is one script, so a page
    can draw the buttons and leave them dead — which is how the ancestor viewer's broke
    when the handler moved out of ancestors.js."""
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    store_properties(fake_storage, job)
    for url in (f"/jobs/{job.id}/ancestors/", f"/jobs/{job.id}/properties/"):
        page = client.get(url).content.decode()
        assert 'class="property-help"' in page, url
        assert "property_help.js" in page, url


def test_properties_page_is_private(client: Client, fake_storage: FakeStorage) -> None:
    assert client.get(f"/jobs/{other_users_job().id}/properties/").status_code == 404


def test_ancestor_viewer_property_panel(client: Client, fake_storage: FakeStorage) -> None:
    """The property panel follows the viewer's own node and gap menus, and is
    absent (not broken) for a job that did not measure properties."""
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    url = f"/jobs/{job.id}/ancestors/"
    # Nothing measured: the section is the offer to measure it, not the panel.
    empty = client.get(url).content.decode()
    assert "No properties have been measured" in empty
    assert 'class="property-plot"' not in empty

    store_properties(fake_storage, job)

    page = client.get(url).content.decode()
    assert "Properties of This Ancestor" in page
    assert "1000 sequences drawn" in page and "95% credible interval" in page
    # the (?) affordance and its note; charge is our own model, so it cites nothing
    assert 'class="property-help"' in page and "pKa" in page
    assert "no method to cite" in page
    assert "Net charge" in page and "<strong>-3.2</strong>" in page
    assert "+0.42 sd" in page  # signed: the offset's direction is random
    assert 'class="property-plot"' in page  # the histogram of the draws
    assert "per-column marginal posteriors, drawn independently" in page  # the caveat
    # It follows the node menu...
    other = client.get(url, {"node": "Node2"}).content.decode()
    assert "<strong>7.5</strong>" in other
    # ...and the gap-version menu: nothing was measured for the unmasked version.
    unmasked = client.get(url, {"node": "Node1", "gaps": "original"}).content.decode()
    assert "<strong>-3.2</strong>" not in unmasked
    assert "No properties have been measured" in unmasked


def test_ancestor_viewer_property_offer(client: Client, fake_storage: FakeStorage) -> None:
    """The offer: what the registry can measure and this ancestor lacks.

    The offer is the *difference* between the registry and the stored summary, so a job
    that measured one property still offers the others — which is how a job finished
    before a property existed ever gets it.
    """
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    url = f"/jobs/{job.id}/ancestors/"
    store_properties(fake_storage, job)  # charge only
    page = client.get(url).content.decode()
    assert "Measure or recompute properties" in page
    # Measured: ticked and disabled, since the list doubles as "what do I have?".
    assert re.search(r'value="charge"[^>]*checked disabled', page)
    assert "already measured" in page
    # Missing: ticked by default on a first load, so the cost shown is the whole offer.
    missing = set(properties.REGISTRY) - {"charge"}
    for key in missing:
        assert f'value="{key}"' in page
    # The hidden sentinel: without it an all-unticked checklist would read as a first
    # load and the cost would be of properties nobody asked for (an unticked checkbox
    # submits nothing).
    assert f'<input type="hidden" name="{properties_run.FIELD}" value="">' in page
    none_ticked = client.get(url, {properties_run.FIELD: ""}).content.decode()
    assert "offer-actions" not in none_ticked  # nothing ticked: the plain form offers no run

    # Scope: the counts come from the trees, and grow with the scope.
    counts = re.findall(r'<option value="(\w+)"[^>]*>[^<]*\((\d+)\)</option>', page)
    scopes = {key: int(count) for key, count in counts if key in properties_run.SCOPES}
    # Scope counts are not one ladder — three combinations of two ancestors each means
    # "this ancestor everywhere" (3) is dearer than "every ancestor here" (2) — so the
    # reader is shown each count rather than told the order.
    assert scopes["node"] == 1
    assert scopes["analogs"] == 3 and scopes["combination"] == 2 and scopes["job"] == 6
    assert '<option value="analogs" selected>' in page  # the default: across combinations

    # One property over more ancestors costs more, and the figure is in the page.
    one = client.get(url, {properties_run.FIELD: "moment"}).content.decode()
    many = client.get(
        url, {properties_run.FIELD: "moment", properties_run.SCOPE_FIELD: "job"}
    ).content.decode()
    assert "<strong>1</strong>" in one
    # The estimate reads as a duration, never as raw CPU-seconds (owner, 2026-10-21).
    assert "CPU-s" not in page
    assert (
        properties_run.cost_text(
            properties_run.total_cost(["moment"], scopes["job"], properties.DEFAULT_DRAWS, 2)
        )
        in many
    )


def test_ancestor_viewer_histogram_counts_its_own_draws(
    client: Client, fake_storage: FakeStorage
) -> None:
    """Each property's histogram reports the draws *it* was measured on.

    The panel-wide count is 0 by design when the rows disagree, which is the normal case
    once a property can be measured later at its own draw count — every bar's tooltip
    read "of 0 draws" until each row used its own.
    """
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    prefix = f"jobs/{job.id}/results/"
    entry = {"label": "x", "units": "", "level": 0.95, "mean": 1.0, "sd": 1.0,
             "low": 0.0, "high": 2.0, "ml": 1.0, "standardised_shift": 0.0}  # fmt: skip
    body = {
        "node": "Node1", "gaps": "fitch", "draws": 200, "seed": "s",
        "sampling": "per-column marginal posteriors, drawn independently",
        "weighted": False,
        # Measured in two runs, so the two disagree on how many draws they used.
        "properties": {"charge": dict(entry, draws=200), "moment": dict(entry, draws=50)},
    }  # fmt: skip
    counts = (("charge", 200), ("moment", 50))
    draws = {key: [float(i % 7) for i in range(1, n + 1)] for key, n in counts}
    fake_storage.store_properties(prefix, "mafft/iqtree/modelfinder", "Node1", "fitch", body, draws)
    page = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
    assert "of 0 draws" not in page
    assert "of 200 draws" in page and "of 50 draws" in page


def test_ancestor_viewer_offer_needs_posteriors(client: Client, fake_storage: FakeStorage) -> None:
    """A reconstruction with no posteriors cannot be drawn from, so nothing is offered on
    it — the panel would otherwise promise a measurement the pipeline cannot make."""
    job = viewer_job(fake_storage, gap_methods=["fitch"])
    add_prank(fake_storage, job)
    page = client.get(f"/jobs/{job.id}/ancestors/", {"gaps": "prank"}).content.decode()
    assert '<option value="prank" selected>' in page
    assert "Measure or recompute properties" not in page


def test_ancestor_viewer_choices_and_fallbacks(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/ancestors/"
    page = client.get(url, {"combination": "mafft/iqtree/WAG", "node": "Node2",
                            "gaps": "original", "extant": "b", "palette": "greys"},
                      headers={"HX-Request": "true"}).content.decode()  # fmt: skip
    assert "<header>" not in page and page.lstrip().startswith("<form")  # the fragment only
    assert '<option value="mafft/iqtree/WAG" selected>' in page
    assert '<option value="b" selected>' in page and '<option value="greys" selected>' in page
    assert "R<br>(0.50)" in page and 'class="tile gap"' not in page
    # extant b (M-K) has a gap at column 2: residue and number, or a dash
    extant = r'<td class="column">(\d+)</td><td class="extant">([^<]+)(?:<br>\((\d+)\))?</td>'
    cells = re.findall(extant, page)
    assert cells == [("1", "M", "1"), ("2", "–", ""), ("3", "K", "2")]
    # switching alignment drops a model from the other alignment; junk falls back
    junk = {"aligner": "muscle", "combination": "mafft/iqtree/WAG", "node": "<script>",
            "gaps": "grasp", "palette": "x"}  # fmt: skip
    page = client.get(url, junk).content.decode()
    assert '<option value="muscle/iqtree/modelfinder" selected>' in page
    assert '<option value="Node2" selected>Node2 · 3' in page  # the root, listed first
    assert "<script>" not in page
    assert '<option value="original" selected>' in page
    assert '<option value="coolwarm" selected>' in page


def test_ancestor_downloads(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/ancestors/download/"
    base = {"combination": "mafft/iqtree/WAG", "node": "Node2", "gaps": "fitch"}
    sequence = client.get(url, base | {"what": "sequence"})
    assert sequence.content == b">Node2 mafft/iqtree/WAG gaps=fitch\nMK\n"
    assert 'filename="Node2.mafft_iqtree_WAG.fitch.fasta"' in sequence["Content-Disposition"]
    assert client.get(url, base | {"what": "aligned"}).content.endswith(b"\nM-K\n")
    posteriors = client.get(url, base | {"what": "posteriors"}).content.decode()
    assert posteriors.splitlines()[0] == "node\tsite\tM\tK\tR"
    assert len(posteriors.splitlines()) == 4 and "Node1" not in posteriors
    assert client.get(url, base | {"what": "states"}).content == (
        b"node\tsite\tstate\tpp\nNode2\t1\tR\t0.5\n"
    )
    for bad in ({"what": "secrets"}, {"what": "sequence", "node": "Node9"}):
        assert client.get(url, base | bad).status_code == 404
    theirs = other_users_job()
    assert client.get(f"/jobs/{theirs.id}/ancestors/").status_code == 404
    assert client.get(f"/jobs/{theirs.id}/ancestors/download/", base).status_code == 404


def test_ancestor_altall_downloads(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/ancestors/download/"
    base = {"combination": "mafft/iqtree/WAG", "node": "Node2", "gaps": "fitch"}
    # Node2's states are M 0.2, K 0.3, R 0.5 per site; masked to M-K by fitch.
    fasta = client.get(url, base | {"what": "altall"}).content.decode()
    assert fasta.splitlines()[0] == (
        ">Node2_AltAll mafft/iqtree/WAG gaps=fitch altall threshold=0.2 switched=2"
    )
    assert fasta.splitlines()[1] == "RR"
    aligned = client.get(url, base | {"what": "altall_aligned"})
    assert aligned.content.decode().endswith("\nR-R\n")  # the mask is kept
    assert (
        'filename="Node2.mafft_iqtree_WAG.fitch.altall0.2.fasta"' in aligned["Content-Disposition"]
    )
    sites = client.get(url, base | {"what": "altall_sites"}).content.decode()
    assert sites.splitlines() == [
        "site\tml\tml_pp\talt\talt_pp",
        "1\tM\t0.20000\tR\t0.50000",
        "3\tK\t0.30000\tR\t0.50000",
    ]
    # The threshold is a control: Node1 is M 0.9 / K 0.1 / R 0 per site, reconstructed MKK,
    # so lowering it switches site 1 as well.
    node1 = {"combination": "mafft/iqtree/WAG", "node": "Node1", "gaps": "original"}
    default = client.get(url, node1 | {"what": "altall"}).content.decode()
    assert default.splitlines()[1] == "MMM" and "switched=2" in default
    low = client.get(url, node1 | {"what": "altall", "altall": "0.05"}).content.decode()
    assert low.splitlines()[1] == "KMM" and "switched=3" in low
    assert "threshold=0.05" in low
    # Junk falls back to the default rather than erroring.
    junk = client.get(url, base | {"what": "altall", "altall": "../etc"}).content.decode()
    assert junk.splitlines()[1] == "RR" and "threshold=0.2" in junk
    # Reconstructions without posteriors (PRANK, GRASP) have no alternatives to offer.
    prefix = f"jobs/{job.id}/results/mafft/iqtree/WAG/"
    fake_storage.objects[prefix + "ancestors.gaps-prank.fasta"] = b">Node1\nMK-\n>Node2\nM-K\n"
    assert client.get(url, base | {"what": "altall", "gaps": "prank"}).status_code == 404
    page = client.get(f"/jobs/{job.id}/ancestors/", base | {"gaps": "prank"}).content.decode()
    assert "AltAll threshold" not in page and "AltAll sequence" not in page


def test_ancestor_viewer_offers_altall_with_a_threshold(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = viewer_job(fake_storage)
    page = client.get(f"/jobs/{job.id}/ancestors/", {"node": "Node2"}).content.decode()
    assert "AltAll sequence (FASTA)" in page and "AltAll threshold" in page
    assert 'name="altall" value="0.2"' in page
    # The slider carries no name, so the number field stays what the form submits.
    assert 'type="range" class="altall-slider" value="0.2"' in page
    # Node2 is MRK with M 0.2 / K 0.3 / R 0.5 at every site: all three are ambiguous at 0.2.
    assert "<strong>3</strong>" in page and "ambiguous sites of 3" in page
    assert '"altall-pps"' in page and "[0.5, 0.5, 0.3]" in page
    page = client.get(
        f"/jobs/{job.id}/ancestors/", {"node": "Node2", "altall": "0.35"}
    ).content.decode()
    assert 'name="altall" value="0.35"' in page
    assert "altall=0.35" in page  # the download links carry it


def test_ancestor_viewer_before_results(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user())
    page = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
    assert "No ancestors to show" in page


def test_ancestor_viewer_keeps_the_clade(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/ancestors/"
    move = {"aligner": "muscle", "combination": "mafft/iqtree/WAG", "node": "Node2"}
    page = client.get(url, move | {"shown": "mafft/iqtree/WAG"}).content.decode()
    assert '<option value="Node1" selected>Node1 · 2 sequences' in page
    # without "shown" (e.g. a bookmarked URL), the name is taken as given
    page = client.get(url, move).content.decode()
    assert '<option value="Node2" selected>' in page
    # an unknown "shown" is ignored
    page = client.get(url, move | {"shown": "../../etc"}).content.decode()
    assert '<option value="Node2" selected>' in page


def test_ancestor_viewer_shows_ml_gap_probabilities(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/ancestors/"
    page = client.get(url, {"node": "Node1", "gaps": "ml"}).content.decode()
    assert "tiles give the probability of a" in page
    assert re.search(r'class="tile gap" style="[^"]+">–<br>\(0\.80\)', page)  # masked
    assert re.search(r'class="tile gap" style="[^"]+">–<br>\(0\.30\)', page)  # beside K
    assert "gap probabilities (TSV)" in page
    page = client.get(url, {"node": "Node1", "gaps": "fitch"}).content.decode()
    assert "gap probabilities (TSV)" not in page and "(0.80)" not in page
    download = f"/jobs/{job.id}/ancestors/download/"
    tsv = client.get(download, {"node": "Node1", "gaps": "ml", "what": "gaps"}).content
    assert tsv.decode().splitlines() == ["node\tsite\tp_gap", "Node1\t1\t0.0",
                                         "Node1\t2\t0.3", "Node1\t3\t0.8"]  # fmt: skip


def test_ancestor_viewer_extant_columns(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/ancestors/"
    hx = {"HX-Request": "true"}
    # + adds the next sequence not yet shown; the address bar gets the view, not the button
    response = client.get(url, {"extant": "a", "add_extant": "1"}, headers=hx)
    page = response.content.decode()
    assert page.count('<th class="extant">') == 2
    assert '<option value="b" selected>' in page
    assert response["HX-Push-Url"].endswith(
        "&extant=a&extant=b&palette=coolwarm&lines=standard&altall=0.2"
    )
    assert "add_extant" not in response["HX-Push-Url"]
    extant = r'<td class="extant">([^<]+)(?:<br>\((\d+)\))?</td>'
    assert re.findall(extant, page)[:2] == [("M", "1"), ("M", "1")]  # column 1: a, b
    assert page.count('name="remove_extant"') == 1  # the first column stays
    # × removes that column; the first cannot be removed; unknown names are dropped
    shown = {"extant": ["a", "b", "c<i>", "nope"]}
    page = client.get(url, shown | {"remove_extant": "1"}, headers=hx).content.decode()
    assert page.count('<th class="extant">') == 2 and '<option value="b" selected>' not in page
    page = client.get(url, shown | {"remove_extant": "0"}, headers=hx).content.decode()
    assert page.count('<th class="extant">') == 3
    # at the cap, no + button; without HTMX the buttons redirect to the resulting view
    many = {"extant": ["a"] * 10}
    assert 'name="add_extant"' not in client.get(url, many).content.decode()
    response = client.get(url, {"extant": "a", "add_extant": "1"})
    assert response.status_code == 302 and "add_extant" not in response["Location"]
    assert response["Location"].endswith(
        "&extant=a&extant=b&palette=coolwarm&lines=standard&altall=0.2"
    )


def test_site_table_as_excel(client: Client, fake_storage: FakeStorage) -> None:
    from io import BytesIO

    from openpyxl import load_workbook

    from phylobot import ancestors

    job = viewer_job(fake_storage)
    page = client.get(f"/jobs/{job.id}/ancestors/", {"extant": ["a", "b"], "palette": "greys"})
    link = re.search(r'<a href="([^"]+)">Download Excel</a>', page.content.decode())
    assert link
    response = client.get(link.group(1).replace("&amp;", "&"))
    assert response["Content-Type"].startswith("application/vnd.openxmlformats")
    assert (
        'filename="Node1.mafft_iqtree_modelfinder.original.sites.xlsx"'
        in (response["Content-Disposition"])
    )
    book = load_workbook(BytesIO(response.content))
    sites = book["Sites"]
    header = [c.value for c in sites[1]]
    assert header[:5] == ["Site in MAFFT", "a residue", "a site", "b residue", "b site"]
    assert header[5:9] == ["State 1", "PP 1", "State 2", "PP 2"]
    # column 2: a has K (residue 2), b a gap; Node1's states M 0.9 then K 0.1
    assert [c.value for c in sites[3]][:9] == [2, "K", 2, "–", None, "M", 0.9, "K", 0.1]
    expected = "{:02X}{:02X}{:02X}".format(*ancestors.rgb(0.9, "greys"))
    assert sites["F3"].fill.start_color.rgb.endswith(expected)
    assert sites["G3"].fill.start_color.rgb.endswith(expected)
    about = {row[0].value: str(row[1].value) for row in book["About"].iter_rows() if row[0].value}
    assert about["Ancestral node"].startswith("Node1 · 3 sequences")
    assert about["Extant sequences"] == "a, b" and about["Palette"].startswith("Greys")
    assert about["ML residues"].startswith("mean posterior probability 0.900")


def test_ancestor_plots_line_styles(client: Client, fake_storage: FakeStorage) -> None:
    from phylobot import ancestors

    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/ancestors/"
    page = client.get(url, {"lines": "dashes"}).content.decode()
    assert '<option value="dashes" selected>' in page
    assert 'stroke-dasharray="6 4"' in page and "data-plot-download" in page
    assert len(ancestors.LINE_STYLES) >= 3
    page = client.get(url, {"lines": "<x>"}).content.decode()
    assert '<option value="standard" selected>' in page


def test_tree_viewer(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/tree/"
    page = client.get(url, {"node": "Node2"}).content.decode()
    assert "<h1>Tree Viewer</h1>" in page and '<option value="adjusted" selected>' in page
    # ancestral nodes link to the ancestor viewer with this combination and node
    link = "/ancestors/?aligner=mafft&amp;combination=mafft%2Fiqtree%2Fmodelfinder&amp;node=Node2"
    assert link in page and 'class="tree-node"' in page
    assert 'stroke="#e08214"' in page  # Node2 highlighted
    assert "c&lt;i&gt;" in page and "c<i>" not in page  # names escaped
    assert "Not rooted on an" in page
    query = {"aligner": "muscle", "lengths": "true"}
    response = client.get(url, query, headers={"HX-Request": "1"})
    fragment = response.content.decode()
    assert '<option value="muscle/iqtree/modelfinder" selected>' in fragment
    assert "substitutions per site" in fragment
    pushed = response["HX-Push-Url"]
    assert "combination=muscle%2Fiqtree%2Fmodelfinder&lengths=true" in pushed
    theirs = other_users_job()
    assert client.get(f"/jobs/{theirs.id}/tree/").status_code == 404
    # the job page and the ancestor viewer link here
    assert f'/jobs/{job.id}/tree/?combination=mafft%2Fiqtree%2FWAG" title="View"' in (
        client.get(f"/jobs/{job.id}/").content.decode()
    )
    ancestors_page = client.get(f"/jobs/{job.id}/ancestors/", {"node": "Node2"}).content.decode()
    assert "tree/?aligner=mafft&amp;combination=mafft%2Fiqtree%2Fmodelfinder&amp;node=Node2" in (
        ancestors_page
    )


def test_tree_viewer_branch_display(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)  # tree ((a,b)Node2,c<i>)Node1; UFBoot 88 at Node2
    url = f"/jobs/{job.id}/tree/"
    page = client.get(url).content.decode()
    menu = page[page.index('name="colour"') : page.index("</select>", page.index('name="colour"'))]
    assert 'value="ufboot"' in menu and 'value="dnds"' not in menu  # no free ratios
    assert 'name="dnds"' not in page and "Node2 <tspan" in page  # labels on by default
    # Unknown or unavailable choices fall back; the URL keeps only resolved ones.
    response = client.get(url, {"colour": "dnds", "width": "<x>", "labels": "nope"},
                          headers={"HX-Request": "1"})  # fmt: skip
    assert "colour=&" not in response["HX-Push-Url"] and "labels=all" in response["HX-Push-Url"]
    assert 'name="branch_palette"' not in page  # palettes only with a branch colour
    query = {"colour": "ufboot", "marks": "stars", "branch_palette": "viridis"}
    page = client.get(url, query).content.decode()
    assert "Branch colour: UFBoot support (%)" in page and "UFBoot * &ge; 80" in page
    assert '<option value="viridis" selected>' in page and 'class="light"' not in page
    assert '<tspan fill="#5d6675">*</tspan>' in page  # 88: one star
    page = client.get(url, {"labels": "none"}).content.decode()
    assert "Node2 <tspan" not in page
    # With codeml's free ratios: dN/dS for colour, thickness and labels.
    results = f"jobs/{job.id}/results/mafft/iqtree/modelfinder/"
    fake_storage.objects[results + "codeml/FR/branches.tsv"] = (
        b"parent\tchild\tt\tN\tS\tdN_dS\tdN\tdS\tN_dN\tS_dS\n"
        b"Node1\tNode2\t0.1\t300\t100\t2.5\t0.05\t0.02\t15\t2\n"
    )
    page = client.get(url, {"width": "dnds", "dnds": "show"}).content.decode()
    assert "Branch thickness: dN/dS (codeml free ratios)" in page
    assert ">2.50</text>" in page and "dN/dS 2.5;" in page


def test_tip_spacing(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/tree/"

    def tips(page: str) -> list[float]:
        return [
            float(y)
            for y in re.findall(r'class="tree-leaf" data-leaf="[^"]*" x="[^"]*" y="([^"]*)"', page)
        ]

    page = client.get(url).content.decode()
    assert 'name="spacing"' in page and 'value="1"' in page
    assert (
        'aria-controls="display-options-body">' in page
        and 'id="display-options-body" hidden' in page
    )
    one = tips(page)
    response = client.get(url, {"spacing": "2"}, headers={"HX-Request": "1"})
    two = tips(response.content.decode())
    assert abs((two[1] - two[0]) - 2 * (one[1] - one[0])) < 0.05
    assert "spacing=2" in response["HX-Push-Url"]
    assert (
        "spacing"
        not in client.get(url, {"spacing": "1"}, headers={"HX-Request": "1"})["HX-Push-Url"]
    )
    assert tips(client.get(url, {"spacing": "99"}).content.decode()) == tips(
        client.get(url, {"spacing": "2.5"}).content.decode()
    )  # clamped


ALIGNED_FOUR = "".join(f">seq{i}\nMKTAYIAK{'LV'[i % 2]}Q-\n" for i in range(4))


def upload(name: str, text: str) -> SimpleUploadedFile:
    return SimpleUploadedFile(name, text.encode())


def test_submit_alignments_only_with_your_own_alignment(
    client: Client, fake_storage: FakeStorage
) -> None:
    response = post(client, aligners=[], run_trees="", run_ancestors="",
                    alignment_files=[upload("curated v2.fasta", ALIGNED_FOUR)])  # fmt: skip
    assert response.status_code == 302, response.content.decode()[-2000:]
    job = Job.objects.get()
    assert job.params["stop_after"] == "alignments" and job.params["aligners"] == ["none"]
    assert job.params["alignment_labels"] == ["curated v2"]
    assert "tree_builders" not in job.params and "gap_methods" not in job.params
    keys = sorted(k.rsplit("/", 1)[1] for k in fake_storage.objects if "/inputs/" in k)
    assert keys == ["alignment1.fasta"]  # no sequences: they come from the alignment


def test_submit_a_phylip_alignment_with_missing_data(
    client: Client, fake_storage: FakeStorage
) -> None:
    # Interleaved PHYLIP with "?" (missing data): stored as FASTA with gaps.
    phylip = (
        "4 11\n"
        + "".join(f"seq{i} MKTAYI\n" for i in range(4))
        + "\n"
        + "".join(f"AK{'LV'[i % 2]}{'Q?' if i == 3 else 'Q-'}\n" for i in range(4))
    )
    page = client.post("/alignment-files/", {"alignment_files": [upload("mine.phy", phylip)]})
    assert "error" not in page.content.decode().lower()
    response = post(client, aligners=[], run_trees="", run_ancestors="",
                    alignment_files=[upload("mine.phy", phylip)])  # fmt: skip
    assert response.status_code == 302, response.content.decode()[-2000:]
    job = Job.objects.get()
    assert job.params["alignment_labels"] == ["mine"]
    stored = next(v for k, v in fake_storage.objects.items() if k.endswith("alignment1.fasta"))
    assert stored.decode() == ALIGNED_FOUR


def test_a_phylip_alignment_that_does_not_fit_its_header(
    client: Client, fake_storage: FakeStorage
) -> None:
    bad = "4 12\n" + "".join(f"seq{i} MKTAYIAKLQ-\n" for i in range(4))
    page = client.post("/alignment-files/", {"alignment_files": [upload("bad.phy", bad)]})
    assert "PHYLIP header says 4 sequences × 12 columns" in page.content.decode()
    page = post(client, aligners=[], run_trees="", run_ancestors="",
                alignment_files=[upload("bad.phy", bad)]).content.decode()  # fmt: skip
    assert "bad.phy: the PHYLIP header says 4 sequences" in page and not Job.objects.exists()


def test_submit_trees_only_with_your_own_tree(client: Client, fake_storage: FakeStorage) -> None:
    tree = "((seq0,seq1),(seq2,seq3));"
    response = post(client, sequences_text=FOUR, run_ancestors="", tree_builders=[],
                    tree_files=[upload("species.nwk", tree)])  # fmt: skip
    assert response.status_code == 302, response.content.decode()[-2000:]
    job = Job.objects.get()
    assert job.params["stop_after"] == "trees" and job.params["tree_labels"] == ["species"]
    assert job.params["tree_builders"] == ["none"] and "gap_methods" not in job.params
    keys = sorted(k.rsplit("/", 1)[1] for k in fake_storage.objects if "/inputs/" in k)
    assert keys == ["sequences.fasta", "tree1.nwk"]


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        ({"sequences_text": FOUR, "aligners": []},
         "Tick an aligner, or upload or choose your own alignment"),
        ({"sequences_text": FOUR, "tree_builders": []},
         "Tick IQ-TREE or RAxML-NG, or upload your own tree"),
        ({"sequences_text": FOUR, "tree_builders": ["raxml"], "fixed_models": "on",
          "models": ["JTTDCMut+G4", "LG+G4"]},
         "RAxML-NG does not have the JTTDCMut matrix (JTTDCMut+G4)"),
        ({"sequences_text": FOUR, "tree_builders": ["raxml"], "fixed_models": "on",
          "models": ["LG+C20+G4", "LG+G4"]},
         "RAxML-NG has no site-profile mixtures (LG+C20+G4)"),
        ({"sequences_text": FOUR, "modelfinder": "on", "candidates": ["LG+C20+G4"]},
         "Select a valid choice"),
        ({"sequences_text": FOUR, "run_trees": "", "tree_files": [upload("t.nwk", "(a,b,c);")]},
         "Your own trees need “Infer trees” ticked"),
        ({"alignment_files": [upload("bad.fasta", ALIGNED_FOUR.replace("seq3", "seqX"))],
          "sequences_text": FOUR}, "Input error: alignment 1"),
        ({"sequences_text": FOUR, "tree_files": [upload("t.nwk", "((seq0,seq1),(seq2,nope));")]},
         "Input error: tree 1"),
    ],
)  # fmt: skip
def test_submit_rejects_bad_compositions(
    client: Client, fake_storage: FakeStorage, fields: dict[str, Any], message: str
) -> None:
    page = post(client, **fields).content.decode()
    assert message in page and not Job.objects.exists()


def test_a_mixture_is_offered_for_a_full_search_only(
    client: Client, fake_storage: FakeStorage
) -> None:
    page = client.get("/jobs/new/").content.decode()
    assert page.count("data-mixture") >= 28 * 2
    cells = re.findall(r'name="(\w+)" value="LG\+C20\+G4"', page)
    assert cells == ["models"]
    response = post(client, sequences_text=FOUR, fixed_models="on", models=["LG+C20+F+G4"],
                    candidates=["LG+G4"])  # fmt: skip
    assert response.status_code == 302
    assert Job.objects.get().params["models"] == ["LG+C20+F+G4"]


def test_ancestor_options_ignored_when_stopping_early(
    client: Client, fake_storage: FakeStorage
) -> None:
    # PRANK gap placement needs an outgroup, but a trees-only job places no gaps
    response = post(client, sequences_text=FOUR, run_ancestors="", gap_methods=["prank"])
    assert response.status_code == 302


def finished(params: dict[str, Any]) -> Job:
    job = Job.objects.create(owner=local_user(), params=params)
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    return job


def test_job_page_for_an_alignment_only_job(client: Client, fake_storage: FakeStorage) -> None:
    job = finished({"stop_after": "alignments", "aligners": ["mafft"],
                    "alignment_labels": ["Curated"]})  # fmt: skip
    prefix = f"jobs/{job.id}/results/"
    summary = {"stop_after": "alignments", "primary": "", "combinations": [], "alignments": [
        {"key": "mafft", "label": "MAFFT", "path": "mafft/alignment.fasta"},
        {"key": "provided1", "label": "Curated", "path": "provided1/alignment.fasta"},
    ]}  # fmt: skip
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    for key in ("mafft", "provided1"):
        fake_storage.objects[f"{prefix}{key}/alignment.fasta"] = ALIGNED_FOUR.encode()
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "alignments only" in page
    assert "Align with MAFFT and your alignment → stop (alignments only)." in page
    assert "This job stopped after the alignments" in page
    assert page.count('<a href="/jobs/') >= 2 and "Curated" in page
    assert "<td>11</td>" in page  # columns of the uploaded alignment
    # Skipped stages' viewers are greyed tiles, not links.
    assert f"/jobs/{job.id}/tree/" not in page and f"/jobs/{job.id}/ancestors/" not in page
    assert page.count('class="viewer-card disabled"') == 2
    assert client.get(f"/jobs/{job.id}/tree/").content.decode().count("No tree to show") == 1


def test_job_page_and_viewers_for_a_trees_only_job(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished({"stop_after": "trees", "aligners": ["mafft"], "tree_builders": ["iqtree"],
                    "modelfinder": True})  # fmt: skip
    prefix = f"jobs/{job.id}/results/"
    key = "mafft/iqtree/modelfinder"
    combo = {"key": key, "aligner": "mafft", "tree": "iqtree", "model": "modelfinder",
             "path": key, "fit": {"model": "LG+G4"}}  # fmt: skip
    mafft = {"key": "mafft", "label": "MAFFT", "path": "mafft/alignment.fasta"}
    summary = {"stop_after": "trees", "primary": key, "combinations": [combo],
               "alignments": [mafft]}  # fmt: skip
    files = {
        "summary.json": json.dumps(summary),
        "mafft/alignment.fasta": ALIGNED_FOUR,
        f"{key}/alignment.fasta": ALIGNED_FOUR,
        f"{key}/tree.nwk": "((seq0,seq1)Node2,seq2,seq3)Node1;",
        f"{key}/supports.tsv": "node\tufboot\nNode2\t97\n",
        f"{key}/fit.tsv": "model\nLG+G4\n",
    }
    for name, text in files.items():
        fake_storage.objects[prefix + name] = text.encode()
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Tree Viewer →" in page and f"/jobs/{job.id}/ancestors/" not in page
    assert page.count('class="viewer-card disabled"') == 1
    assert "stopped after the trees" in page and ">site probabilities<" not in page
    assert "<th>Ancestors</th>" not in page
    tree = client.get(f"/jobs/{job.id}/tree/").content.decode()
    assert "Node2 97" in tree and 'class="tree-node"' not in tree  # labelled, not linked
    assert "Numbers beside nodes are UFBoot support (%)." in tree  # jobs before SH-aLRT
    both = "node\tufboot\tsh_alrt\nNode1\t\t\nNode2\t97\t85.5\n"
    fake_storage.objects[prefix + f"{key}/supports.tsv"] = both.encode()
    tree = client.get(f"/jobs/{job.id}/tree/").content.decode()
    assert "Node2 85.5/97" in tree and "are SH-aLRT / UFBoot support (%)." in tree
    none = "node\tufboot\tsh_alrt\nNode1\t\t\nNode2\t\t\n"
    fake_storage.objects[prefix + f"{key}/supports.tsv"] = none.encode()
    tree = client.get(f"/jobs/{job.id}/tree/").content.decode()
    assert ">Node2</text>" in tree and "No branch support was computed" in tree
    assert "Mutational paths need ancestral sequences" in tree
    viewer = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
    assert "this job did not" in viewer


def test_homologs_page_and_saved_searches(client: Client) -> None:
    page = client.get("/homologs/").content.decode()
    assert "Find Homologs" in page and "NCBI · ClusteredNR" in page
    assert 'id="homolog-filters"' not in page  # nothing to review until one is opened
    assert "Load from Homologs" in client.get("/jobs/new/").content.decode()

    page = client.get("/homologs/").content.decode()
    # The Ime2 examples are no longer featured on the page (kept as fixtures, reachable by
    # URL); nothing is featured yet, so the "Try a Saved Search" section is hidden.
    assert 'href="/homologs/examples/ime2-ebi/"' not in page and "Try a Saved Search" not in page
    old_link = client.get("/homologs/?example=ime2-ebi&min_coverage=55")
    assert old_link.status_code == 302
    assert old_link["Location"] == "/homologs/examples/ime2-ebi/?min_coverage=55"

    page = client.get("/homologs/examples/ime2-ebi/").content.decode()
    assert re.findall(r'class="count">([\d,]+)<', page) == ["1,000", "1,000", "74", "62", "63"]
    assert "stopped at its limit (1,000 hits)" in " ".join(page.split())
    assert "Search results: 63 selected" in page and "937 not selected" in page
    assert "fails identity 30-90%" in page  # the query itself: 100% identical

    # The filters re-request the page; HTMX gets only the review panel.
    lower = "/homologs/examples/ime2-ebi/?max_evalue=1e-05&min_coverage=55"
    panel = client.get(lower + "&min_identity=30&max_identity=90", HTTP_HX_REQUEST="true")
    text = panel.content.decode()
    assert text.lstrip().startswith('<section id="homolog-filters"')
    assert re.findall(r'class="count">([\d,]+)<', text)[2:] == ["638", "626", "627"]
    assert "Yarrowia lipolytica" in text  # the ortholog 70% removed is back

    bad = client.get(lower.replace("55", "150") + "&min_identity=30&max_identity=90")
    assert "less than or equal to 100" in bad.content.decode()
    assert client.get("/homologs/?example=../../etc/passwd").status_code == 404
    assert client.get("/homologs/examples/nope/").status_code == 404


@override_settings(PHYLOBOT_MODE="cloud", MIDDLEWARE=middleware_for("cloud"))
def test_homologs_hides_ncbi_in_cloud_mode(client: Client) -> None:
    assert client.get("/homologs/").status_code == 302  # login required
    client.force_login(User.objects.create_user("someone"))
    page = client.get("/homologs/").content.decode()
    assert 'value="ebi:uniref90"' in page and 'value="ncbi:' not in page
    assert "NCBI" not in page.split("<legend>")[1].split("</fieldset>")[0]


def test_rename_job(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user(), name="old")
    response = client.post(
        f"/jobs/{job.id}/rename/", {"name": "  Ime2   fungi <b>x</b> "}, HTTP_HX_REQUEST="true"
    )
    job.refresh_from_db()
    assert job.name == "Ime2 fungi <b>x</b>"
    text = response.content.decode()
    assert text.lstrip().startswith('<div class="editable-title" id="job-title"')
    assert "&lt;b&gt;" in text and "<html" not in text
    client.post(f"/jobs/{job.id}/rename/", {"name": "x" * 500})
    job.refresh_from_db()
    assert len(job.name) == jobs.MAX_NAME
    theirs = other_users_job()
    assert client.post(f"/jobs/{theirs.id}/rename/", {"name": "mine"}).status_code == 404
    assert client.get(f"/jobs/{job.id}/rename/").status_code == 405
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert 'id="job-title"' in page and "data-rename" in page and "phylobot/rename.js" in page


def test_download_zip(client: Client, fake_storage: FakeStorage) -> None:
    import io
    import zipfile

    job = Job.objects.create(owner=local_user(), name="Ime2 / fungi")
    prefix = f"jobs/{job.id}/"
    fake_storage.objects[prefix + "inputs/sequences.fasta"] = b">a\nMK\n"
    fake_storage.objects[prefix + "results/mafft/iqtree/LG/tree.nwk"] = b"(a,b);\n"
    fake_storage.objects[prefix + "results/summary.json"] = b"{}" * 10_000
    fake_storage.objects[f"jobs/{job.id}x/inputs/other"] = b"not this job's"
    response = client.get(f"/jobs/{job.id}/download/")
    folder = f"Ime2_fungi-{str(job.id)[:8]}"
    assert response["Content-Type"] == "application/zip"
    assert response["Content-Disposition"] == f'attachment; filename="{folder}.zip"'
    content = b"".join(response.streaming_content)  # type: ignore[attr-defined]
    archive = zipfile.ZipFile(io.BytesIO(content))
    assert archive.testzip() is None
    assert sorted(archive.namelist()) == [
        f"{folder}/inputs/sequences.fasta",
        f"{folder}/phylobot-job.json",  # the record an import reads (phylobot.job_import)
        f"{folder}/results/mafft/iqtree/LG/tree.nwk",
        f"{folder}/results/summary.json",
    ]
    assert archive.read(f"{folder}/results/summary.json") == b"{}" * 10_000
    unnamed = Job.objects.create(owner=local_user())
    assert (
        f"job-{str(unnamed.id)[:8]}.zip"
        in client.get(f"/jobs/{unnamed.id}/download/")["Content-Disposition"]
    )
    assert client.get(f"/jobs/{other_users_job().id}/download/").status_code == 404


def test_duplicate_opens_a_filled_new_job_page(client: Client, fake_storage: FakeStorage) -> None:
    params = {"aligners": ["mafft", "muscle"], "tree_builders": ["iqtree"], "modelfinder": True,
              "fixed_models": True, "models": ["WAG+G4"], "candidates": ["LG+G4", "JTT+G4"],
              "alrt": 0, "bootstrap": 1000, "seed": 7, "gap_methods": ["fitch"],
              "outgroup": ["seq3"], "stop_after": "ancestors", "seqtype": "protein"}  # fmt: skip
    old = Job.objects.create(owner=local_user(), name="first", params=params)
    old.add_tags(["ime2"])
    assert jobs.mark(old.id, Job.Status.RUNNING) and jobs.mark(old.id, Job.Status.FAILED, "x")
    fake_storage.objects[f"jobs/{old.id}/inputs/sequences.fasta"] = FOUR.encode()
    page = client.get(f"/jobs/new/?from={old.id}").content.decode()
    assert "Duplicating “first”" in page and f'value="{old.id}"' in page  # hidden from_job
    assert "seq0" in page and 'name="seed" value="7"' in page  # sequences and settings
    for value in ("mafft", "muscle"):
        assert re.search(rf'name="aligners" value="{value}"[^>]*checked', page)
    assert re.search(r'name="candidates" value="LG\+G4"[^>]*checked', page)
    assert not re.search(r'name="candidates" value="WAG\+G4"[^>]*checked', page)
    assert re.search(r'name="models" value="WAG\+G4"[^>]*checked', page)
    assert re.search(r'name="outgroup" value="seq3"\s*checked', page)
    assert not Job.objects.exclude(id=old.id).exists()  # nothing runs until Submit
    # Submitting (here unchanged but for the seed) queues a new job linked to the old one.
    response = post(client, sequences_text=FOUR, from_job=str(old.id), seed=8)
    new = Job.objects.exclude(id=old.id).get()
    assert response["Location"] == f"/jobs/{new.id}/"
    assert new.rerun_of_id == old.id and new.params["seed"] == 8
    assert list(new.tags.values_list("name", flat=True)) == ["ime2"]
    old.refresh_from_db()
    assert (old.status, old.error, old.deleted_at) == (Job.Status.FAILED, "x", None)
    # Deleting the original later links it to its copy in the Deleted Jobs list.
    jobs.delete([old])
    listing = client.get("/jobs/").content.decode()
    deleted = listing[listing.index("Deleted Jobs") :]
    assert f'duplicated as\n                <a href="/jobs/{new.id}/">' in deleted


def test_duplicate_keeps_the_codon_tree_choices(client: Client, fake_storage: FakeStorage) -> None:
    params = {"seqtype": "codon", "aligners": ["mafft"], "tree_builders": ["iqtree"],
              "modelfinder": True, "codon_modelfinder": True, "codon_trees": ["amino"],
              "codon_amino_trees": "all", "stop_after": "trees"}  # fmt: skip
    old = Job.objects.create(owner=local_user(), name="codon", params=params)
    fake_storage.objects[f"jobs/{old.id}/inputs/sequences.fasta"] = LYSIN.encode()
    page = client.get(f"/jobs/new/?from={old.id}").content.decode()
    assert re.search(r'name="codon_on_amino"[^>]*checked', page)
    assert not re.search(r'name="codon_search"[^>]*checked', page)
    assert re.search(r'name="codon_amino_trees" value="all"[^>]*checked', page)


def test_duplicate_reuses_uploaded_alignments_unless_replaced(
    client: Client, fake_storage: FakeStorage
) -> None:
    params = {"aligners": ["none"], "stop_after": "alignments", "seqtype": "protein",
              "alignment_labels": ["mine"], "outgroup": [], "seed": 1}  # fmt: skip
    old = Job.objects.create(owner=local_user(), name="aln", params=params)
    rows = "".join(f">seq{i}\nMKTA-YIAK{'LV'[i % 2]}Q\n" for i in range(4))
    fake_storage.objects[f"jobs/{old.id}/inputs/alignment1.fasta"] = rows.encode()
    page = client.get(f"/jobs/new/?from={old.id}").content.decode()
    assert "Its uploaded alignment (mine) is reused" in re.sub(r"\s+", " ", page)
    assert re.search(r'name="outgroup" value="seq0"', page)  # names from the alignment
    post(client, from_job=str(old.id), aligners=[], run_trees="")
    new = Job.objects.exclude(id=old.id).get()
    assert fake_storage.objects[f"jobs/{new.id}/inputs/alignment1.fasta"] == rows.encode()
    assert new.params["alignment_labels"] == ["mine"]


def test_duplicate_only_your_own_jobs(client: Client, fake_storage: FakeStorage) -> None:
    theirs = other_users_job()
    assert client.get(f"/jobs/new/?from={theirs.id}").status_code == 404
    assert client.get("/jobs/new/?from=not-a-uuid").status_code == 404
    response = post(client, sequences_text=FOUR, from_job=str(theirs.id))
    assert response.status_code == 404 and not Job.objects.filter(owner=local_user()).exists()


def ticked(page: str, name: str) -> set[str]:
    return set(re.findall(rf'name="{name}" value="([^"]+)"[^>]*checked', page))


def test_new_job_defaults_give_a_library_axes(client: Client) -> None:
    page = client.get("/jobs/new/").content.decode()
    assert ticked(page, "aligners") == {"mafft", "muscle", "clustalo"}
    assert ticked(page, "models") == {"LG+G4", "WAG+G4"}
    assert ticked(page, "gap_methods") == {"parsimony", "fitch"}
    assert re.search(r'name="fixed_models"[^>]*checked', page)
    assert re.search(r'name="modelfinder"[^>]*checked', page)
    # Each step's "Use defaults" button, and the values it sets (the fields' initials).
    for stage in ("sequences", "alignments", "trees", "ancestors"):
        assert f'data-defaults="{stage}">Use defaults</button>' in page
    defaults = json.loads(page.split('id="section-defaults" type="application/json">')[1]
                          .split("</script>")[0])  # fmt: skip
    assert defaults["alignments"] == {"aligners": ["mafft", "muscle", "clustalo"]}
    assert defaults["trees"]["models"] == ["LG+G4", "WAG+G4"] and defaults["trees"]["sh_alrt"]
    assert defaults["ancestors"]["gap_methods"] == ["parsimony", "fitch"]
    assert "sequences_text" not in defaults["sequences"] and "outgroup" not in defaults["trees"]


def test_duplicate_keeps_its_own_settings_not_the_defaults(
    client: Client, fake_storage: FakeStorage
) -> None:
    params = {"aligners": ["prank"], "tree_builders": ["iqtree"], "modelfinder": True,
              "gap_methods": ["ml"], "stop_after": "ancestors", "seqtype": "protein"}  # fmt: skip
    old = Job.objects.create(owner=local_user(), name="thin", params=params)
    fake_storage.objects[f"jobs/{old.id}/inputs/sequences.fasta"] = FOUR.encode()
    page = client.get(f"/jobs/new/?from={old.id}").content.decode()
    assert ticked(page, "aligners") == {"prank"} and not ticked(page, "models")
    assert ticked(page, "gap_methods") == {"ml"}


def test_job_actions_on_both_pages(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user(), name="one")
    listing = client.get("/jobs/").content.decode()
    assert f'href="/jobs/{job.id}/download/"' in listing
    assert f'href="/jobs/new/?from={job.id}"' in listing
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert f'href="/jobs/{job.id}/download/"' in page and f'href="/jobs/new/?from={job.id}"' in page
    assert 'id="confirm-action"' in page


def test_new_job_page_has_two_model_grids(client: Client) -> None:
    page = client.get("/jobs/new/").content.decode()
    ticked = {name: len(re.findall(rf'name="{name}" value="[^"]+"\s+title="[^"]+" '
                                   r'data-col="\d+"\s+checked', page))
              for name in ("candidates", "models")}  # fmt: skip
    # Every candidate ticked (IQ-TREE's default ModelFinder set); full searches under
    # the two default models (forms.DEFAULT_MODELS).
    assert ticked == {"candidates": 252, "models": 2}
    assert page.count('type="checkbox" name="models"') == 252 + 56


def test_template_comments_are_single_line() -> None:
    """Django's {# #} comments end at the line's end: a multi-line one renders as text."""
    templates = Path(__file__).parent.parent / "src/phylobot/templates"
    assert len(list(templates.rglob("*.html"))) > 20
    bad = [f"{path.name}:{n}" for path in templates.rglob("*.html")
           for n, line in enumerate(path.read_text().splitlines(), 1)
           if "{#" in line and "#}" not in line]  # fmt: skip
    assert not bad


def test_work_estimate_page(client: Client) -> None:
    mine = Job.objects.create(owner=local_user(), name="mine", status=Job.Status.SUCCEEDED,
                              work_units=2.0, cpu_seconds=600.0)  # fmt: skip
    Job.objects.create(owner=User.objects.create_user("bob"), name="theirs",
                       status=Job.Status.SUCCEEDED, work_units=1.0, cpu_seconds=60.0)  # fmt: skip
    page = client.get("/settings/work-estimate/").content.decode()
    assert "Reference Costs" in page and "IQ-TREE search" in page and "+G4 or +R4" in page
    # Only the viewer's own jobs; 600 s against 2 units x 60 s is 5x: flagged.
    assert f"/jobs/{mine.id}/" in page and "theirs" not in page
    assert "5.00 ⚑" in page
    # Two jobs are not enough to calibrate (MIN_HISTORY = 3).
    assert "Not calibrated yet" in page
    # Linked from the Settings page and from the New Job page's estimate.
    assert "/settings/work-estimate/" in client.get("/settings/").content.decode()
    preview = client.post("/preview/", {"aligners": ["mafft"], "run_trees": "on"})
    assert "/settings/work-estimate/" in preview.content.decode()


# Coding DNA

LYSIN = (ROOT / "src" / "phylobot" / "examples" / "lysin.fasta").read_text()


def test_coding_dna_detected_and_submitted(client: Client, fake_storage: FakeStorage) -> None:
    response = post(client, sequences_text=LYSIN, seqtype="auto", genetic_code="1",
                    aligners=["mafft", "prank_codon"], analyses="amino")  # fmt: skip
    assert response.status_code == 302
    job = Job.objects.get()
    assert job.params["seqtype"] == "codon" and job.params["genetic_code"] == 1
    assert job.params["aligners"] == ["mafft", "prank_codon"]
    assert "codon_models" not in job.params  # amino-acid models only
    # The estimate counts codons (lysin: about 130), not nucleotides.
    assert job.work_units and job.work_units > 0


def test_protein_jobs_record_their_type_without_a_genetic_code(
    client: Client, fake_storage: FakeStorage
) -> None:
    post(client, sequences_text=FOUR, genetic_code="2")
    job = Job.objects.get()
    assert job.params["seqtype"] == "protein" and "genetic_code" not in job.params


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        ({"sequences_text": FOUR, "aligners": ["prank_codon"]}, "codon model needs coding DNA"),
        ({"sequences_text": LYSIN, "seqtype": "protein"}, "--seqtype codon"),
        ({"sequences_text": FOUR, "seqtype": "codon", "analyses": "amino"}, "look like protein"),
        ({"sequences_text": LYSIN, "analyses": "codon"}, "Choose codon ModelFinder"),
        (
            {"sequences_text": LYSIN, "analyses": "codon", "codon_fixed": "on"},
            "Tick the codon models to search in full",
        ),
        (
            {
                "sequences_text": LYSIN,
                "analyses": "codon",
                "codon_modelfinder": "on",
                "tree_builders": ["raxml"],
            },
            "RAxML-NG has no codon models",
        ),
        (
            {"sequences_text": LYSIN, "codon_modelfinder": "on", "codon_search": ""},
            "Codon models need a tree",
        ),
        (
            {
                "sequences_text": LYSIN,
                "analyses": "codon",
                "codon_modelfinder": "on",
                "codon_on_amino": "on",
            },
            "needs amino-acid models",
        ),  # fmt: skip
        ({"sequences_text": LYSIN, "codon_models": ["MG+F"]}, "Select a valid choice"),
        ({"sequences_text": LYSIN, "genetic_code": "7"}, "Select a valid choice"),
    ],
)
def test_coding_dna_choices_checked(
    client: Client, fake_storage: FakeStorage, fields: dict[str, Any], message: str
) -> None:
    response = post(client, **fields)
    assert response.status_code == 200 and not Job.objects.exists()
    assert message in response.content.decode()


def check(client: Client, **fields: Any) -> str:
    response = client.post("/sequence-check/", fields)
    assert response.status_code == 200
    return re.sub(r"\s+", " ", response.content.decode())


def test_sequence_check_reports_type_size_and_problems(client: Client) -> None:
    page = check(client, sequences_text=LYSIN)
    assert 'data-seqtype="codon"' in page
    assert "✓ 25 coding DNA sequences, 125–135 codons each" in page and "(detected" in page
    assert "✓ 4 protein sequences, 10 residues each" in check(client, sequences_text=FOUR)
    frameshift = LYSIN.replace("ATG", "AT", 1)
    assert "is not a multiple of 3 (a frameshift" in check(client, sequences_text=frameshift)
    stop = LYSIN.split("\n")
    stop[1] = stop[1][:30] + "TAA" + stop[1][33:]
    assert "internal stop codon TAA at codon 11" in check(client, sequences_text="\n".join(stop))
    # The chosen type wins over detection, and the genetic code is applied.
    assert "look like protein" in check(client, sequences_text=FOUR, seqtype="codon")
    mito = ">a\nATGTGA\n>b\nATGTGG\n>c\nATGTGG\n>d\nATGTGA\n"
    assert "internal stop codon" not in check(client, sequences_text=mito, genetic_code="2")
    assert 'data-seqtype=""' in check(client, sequences_text="")


def ambiguous_lysin() -> str:
    """The lysin example with two ambiguous codons: AMC and NNN."""
    lines = LYSIN.split("\n")
    lines[1] = lines[1][:3] + "AMC" + lines[1][6:]  # AAC or ACC: Asn or Thr
    lines[3] = lines[3][:3] + "NNN" + lines[3][6:]
    return "\n".join(lines)


def test_sequence_check_asks_to_confirm_ambiguous_codons(client: Client) -> None:
    page = check(client, sequences_text=ambiguous_lysin())
    assert "✓ 25 coding DNA sequences" in page
    assert "2 codons in 2 sequences have ambiguous bases" in page
    assert 'id="ambiguity-dialog"' in page and "data-fingerprint=" in page
    assert "<code>AMC</code></td><td>Asn or Thr</td><td>X (unknown)</td>" in page
    assert "<code>NNN</code></td><td>any of 20 amino acids or a stop</td>" in page
    assert "codon models</strong> treat the whole codon as unknown" in page
    assert "ambiguous bases: 2 codons" not in page  # the pipeline's warning, shown as the line
    assert "ambiguity-dialog" not in check(client, sequences_text=LYSIN)


def test_ambiguous_codons_must_be_confirmed(client: Client, fake_storage: FakeStorage) -> None:
    text = ambiguous_lysin()
    response = post(client, sequences_text=text, analyses="amino")
    assert response.status_code == 200 and not Job.objects.exists()
    assert "codons with ambiguous bases: review them" in response.content.decode()
    fingerprint = re.search(r'data-fingerprint="(\w+)"', check(client, sequences_text=text))
    assert fingerprint
    stale = post(client, sequences_text=text, analyses="amino", ambiguity_ack="0" * 16)
    assert stale.status_code == 200 and not Job.objects.exists()
    response = post(client, sequences_text=text, analyses="amino",
                    ambiguity_ack=fingerprint.group(1))  # fmt: skip
    assert response.status_code == 302 and Job.objects.count() == 1


def test_lysin_example_refreshes_the_check(client: Client) -> None:
    page = client.get("/examples/lysin/").content.decode()
    assert "Haliotis_rufescens" in page and 'id="sequence-check"' in page
    assert 'data-seqtype="codon"' in page and 'hx-swap-oob="true"' in page


def test_preview_counts_codon_steps(client: Client) -> None:
    page = client.post("/preview/", {
        "sequences_text": LYSIN, "aligners": ["mafft", "prank_codon"], "run_trees": "on",
        "run_ancestors": "on", "tree_builders": ["iqtree"], "modelfinder": "on",
    }).content.decode()  # fmt: skip
    assert "Check and translate the coding DNA (the standard code)" in page
    assert "Align with MAFFT and PRANK (codon model), each also as a codon alignment" in page
    protein = client.post("/preview/", {"sequences_text": FOUR, "aligners": ["prank_codon"]})
    assert "PRANK (codon model)" not in protein.content.decode()  # coding DNA only


def test_codon_models_both_and_codon_only(client: Client, fake_storage: FakeStorage) -> None:
    # Both (the default for coding DNA): protein ModelFinder plus codon ModelFinder among
    # the ticked candidates, plus a full search under each ticked codon model.
    post(client, sequences_text=LYSIN, codon_modelfinder="on",
         codon_candidates=["MG+F1X4+G4", "GY+F3X4"], codon_fixed="on",
         codon_models=["GY+F3X4+G4"])  # fmt: skip
    both = Job.objects.get()
    assert both.params["modelfinder"] is True and both.params["codon_modelfinder"] is True
    assert both.params["codon_candidates"] == ["MG+F1X4+G4", "GY+F3X4"]
    assert both.params["codon_models"] == ["GY+F3X4+G4"]
    Job.objects.all().delete()
    # Codon models only: the amino-acid options are ignored, whatever they say; a codon
    # grid whose option is unticked is ignored too.
    post(client, sequences_text=LYSIN, analyses="codon", codon_modelfinder="on",
         codon_models=["GY+F3X4+G4"], fixed_models="on", models=["LG"])  # fmt: skip
    codon = Job.objects.get()
    assert codon.params["modelfinder"] is False and codon.params["fixed_models"] is False
    assert codon.params["models"] is None and codon.params["candidates"] is None
    assert codon.params["codon_modelfinder"] is True and "codon_models" not in codon.params
    Job.objects.all().delete()
    # Protein jobs never send codon options, even if the (hidden) controls were ticked.
    post(client, sequences_text=FOUR, codon_modelfinder="on", codon_fixed="on",
         codon_models=["GY+F3X4+G4"])  # fmt: skip
    params = Job.objects.get().params
    assert "codon_models" not in params and "codon_modelfinder" not in params


def test_preview_describes_codon_models(client: Client) -> None:
    fields = {"sequences_text": LYSIN, "aligners": ["mafft"], "run_trees": "on",
              "tree_builders": ["iqtree"], "modelfinder": "on", "codon_modelfinder": "on",
              "codon_fixed": "on", "codon_models": ["GY+F3X4+G4"],
              "codon_search": "on"}  # fmt: skip
    both = client.post("/preview/", fields).content.decode()
    assert ("best model by ModelFinder, best codon model by ModelFinder and codon model "
            "GY+F3X4+G4") in both  # fmt: skip
    assert "<strong>3</strong> combination" in re.sub(r"\s+", " ", both)
    codon = client.post("/preview/", fields | {"analyses": "codon"}).content.decode()
    plan = codon.split("This job will")[1].split("</ol>")[0]
    assert "best model by ModelFinder" not in plan and "best codon model" in plan


def test_codon_models_on_the_amino_acid_trees(client: Client, fake_storage: FakeStorage) -> None:
    # Owner decision 2026-10-02: the best amino-acid tree by default, every one on request;
    # a codon search too unless unticked.
    post(client, sequences_text=LYSIN, codon_modelfinder="on", codon_on_amino="on",
         codon_search="", tree_builders=["iqtree", "raxml"])  # fmt: skip
    params = Job.objects.get().params
    assert params["codon_trees"] == ["amino"] and params["codon_amino_trees"] == "best"
    Job.objects.all().delete()
    post(client, sequences_text=LYSIN, codon_modelfinder="on", codon_on_amino="on",
         codon_amino_trees="all")  # fmt: skip
    params = Job.objects.get().params
    assert params["codon_trees"] == ["search", "amino"] and params["codon_amino_trees"] == "all"
    Job.objects.all().delete()
    post(client, sequences_text=LYSIN, codon_modelfinder="on")
    params = Job.objects.get().params
    assert params["codon_trees"] == ["search"] and "codon_amino_trees" not in params
    Job.objects.all().delete()
    # Only the user's own tree: no IQ-TREE or RAxML-NG trees for the codon models.
    post(client, sequences_text=FOUR, codon_modelfinder="on", codon_on_amino="on")
    assert "codon_trees" not in Job.objects.get().params  # protein: no codon options


def test_preview_describes_codon_fits_on_amino_acid_trees(client: Client) -> None:
    fields = {"sequences_text": LYSIN, "aligners": ["mafft"], "run_trees": "on",
              "tree_builders": ["iqtree", "raxml"], "modelfinder": "on",
              "codon_modelfinder": "on", "codon_on_amino": "on"}  # fmt: skip
    page = re.sub(r"\s+", " ", client.post("/preview/", fields).content.decode())
    assert "fit the codon models on the best amino-acid tree (no codon tree search)" in page
    # IQ-TREE and RAxML-NG: 2 amino-acid combinations, 2 codon ones (one per tree source).
    assert "<strong>4</strong> combination" in page
    every = fields | {"codon_amino_trees": "all", "codon_search": "on", "fixed_models": "on",
                      "models": ["LG+G4"]}  # fmt: skip
    page = re.sub(r"\s+", " ", client.post("/preview/", every).content.decode())
    assert "fit the codon models on every amino-acid tree" in page
    assert "no codon tree search" not in page
    # 4 amino-acid; codon: 1 search + 2 sources x 2 amino-acid models.
    assert "<strong>9</strong> combination" in page


def test_new_job_page_has_the_codon_grids(client: Client) -> None:
    page = client.get("/jobs/new/").content.decode()
    assert "data-codon-models" in page and 'name="codon_modelfinder"' in page
    candidates = re.findall(r'name="codon_candidates" value="([^"]+)"[^>]*checked', page)
    assert len(candidates) == 76  # all ticked: IQ-TREE's default codon set
    assert not re.findall(r'name="codon_models" value="[^"]+"[^>]*checked', page)
    # Codon trees: both a codon search and a fit on the amino-acid tree (owner decision).
    assert re.search(r'name="codon_search"[^>]*checked', page)
    assert re.search(r'name="codon_on_amino"[^>]*checked', page)
    assert 'value="MG+F"' not in page and "IQ-TREE does not offer MG+F" in page


def test_codeml_options(client: Client, fake_storage: FakeStorage) -> None:
    codon = {"sequences_text": LYSIN, "codon_modelfinder": "on"}
    post(client, **codon, codeml="on", codeml_models=["M0", "M3", "M7", "M8"])
    params = Job.objects.get().params
    assert params["codeml_models"] == ["M0", "M3", "M7", "M8"]
    assert params["codeml_trees"] == "best"  # owner decision: each alignment's best
    Job.objects.all().delete()
    # Ignored without codon models (amino-acid models only) and for protein input.
    post(client, sequences_text=LYSIN, analyses="amino", codeml="on", codeml_models=["M0"])
    assert "codeml_models" not in Job.objects.get().params
    Job.objects.all().delete()
    post(client, sequences_text=FOUR, codeml="on", codeml_models=["M0"])
    assert "codeml_models" not in Job.objects.get().params
    Job.objects.all().delete()
    refused = post(client, **codon, codeml="on").content.decode()
    assert "Tick at least one codeml model" in refused
    mito_sr1 = post(client, **codon, codeml="on", codeml_models=["M0"], genetic_code="25")
    assert "PAML doesn&#x27;t have genetic code 25" in mito_sr1.content.decode()
    assert not Job.objects.exists()


def test_preview_names_codeml(client: Client) -> None:
    fields = {"sequences_text": LYSIN, "aligners": ["mafft"], "run_trees": "on",
              "run_ancestors": "on", "tree_builders": ["iqtree"], "codon_modelfinder": "on",
              "codeml": "on", "codeml_models": ["M1a", "M2a"], "codeml_trees": "all"}  # fmt: skip
    page = client.post("/preview/", fields).content.decode()
    assert "Run PAML codeml (M1a, M2a) on every codon combination" in page
    page = client.get("/jobs/new/").content.decode()
    ticked = re.findall(r'name="codeml_models" value="(M\w+)"[^>]*checked', page)
    assert ticked == ["M0", "M3", "M7", "M8"]
    assert "M1a vs M2a positive-selection test" in page


def test_preview_counts_a_duplicated_jobs_reused_files(client: Client) -> None:
    params = {"aligners": ["mafft"], "tree_labels": ["mine"], "tree_builders": ["none"]}
    old = Job.objects.create(owner=local_user(), params=params)
    fields = {"sequences_text": FOUR, "aligners": ["mafft"], "run_trees": "on",
              "modelfinder": "on", "from_job": str(old.id)}  # fmt: skip
    page = re.sub(r"\s+", " ", client.post("/preview/", fields).content.decode())
    assert "Fit your tree" in page and "<strong>1</strong> combination" in page


def test_ancestor_altall_follows_the_states_menu_into_codons(
    client: Client, fake_storage: FakeStorage
) -> None:
    """A codon reconstruction gets codon AltAll: whole codons switched, the CDS out."""
    job = viewer_job(fake_storage)
    results = f"jobs/{job.id}/results/mafft/iqtree/WAG/"
    fake_storage.objects[results + "ancestors.cds.fasta"] = (
        b">Node1\nATGAAAAAA\n>Node2\nATGAGAAAG\n"
    )
    fake_storage.objects[results + "ancestors_codon_posteriors.tsv"] = (
        b"node\tsite\tcodon\tpp\n"
        b"Node2\t1\tATG\t1.0\nNode2\t2\tAGA\t0.7\nNode2\t2\tCGA\t0.3\n"
        b"Node2\t3\tAAG\t0.9\nNode2\t3\tAAA\t0.1\n"
    )
    url = f"/jobs/{job.id}/ancestors/"
    base = {"combination": "mafft/iqtree/WAG", "node": "Node2", "gaps": "original"}
    page = client.get(url, base | {"states": "codons"}).content.decode()
    assert "<strong>1</strong>" in page and "ambiguous codon of 3" in page  # singular
    download = f"/jobs/{job.id}/ancestors/download/"
    fasta = client.get(
        download, base | {"states": "codons", "what": "altall_aligned"}
    ).content.decode()
    assert fasta.splitlines()[1] == "ATGCGAAAG"  # only codon 2 reaches 0.2
    assert "coding DNA" in fasta.splitlines()[0]
    sites = client.get(download, base | {"states": "codons", "what": "altall_sites"}).content
    assert sites.decode().splitlines()[1].startswith("2\tAGA\t0.70000\tCGA")
    # In amino-acid mode the same node gets the protein AltAll, as before.
    page = client.get(url, base).content.decode()
    assert "ambiguous sites of 3" in page


def test_viewer_reads_one_node_through_the_byte_index(
    client: Client, fake_storage: FakeStorage, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """EXPORT indexes ancestors_posteriors.tsv by node; the viewer reads one node's rows
    by ranged reads and never the whole table, falling back to it without an index."""
    from phylobot import storage
    from phylobot.pipeline import export

    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/ancestors/"
    base = {"combination": "mafft/iqtree/WAG", "node": "Node2", "gaps": "original"}
    whole = f"jobs/{job.id}/results/mafft/iqtree/WAG/ancestors_posteriors.tsv"
    # Without an index: the whole table, as for an older job. Node2 has 3 ambiguous sites.
    assert "<strong>3</strong>" in client.get(url, base).content.decode()
    local = tmp_path / "ancestors_posteriors.tsv"
    local.write_bytes(fake_storage.objects[whole])
    export.index_by_node(local)
    fake_storage.objects[whole.replace(".tsv", ".index.tsv")] = (
        tmp_path / "ancestors_posteriors.index.tsv"
    ).read_bytes()
    read = fake_storage.read_bytes
    monkeypatch.setattr(
        storage, "read_bytes", lambda key: pytest.fail("read whole") if key == whole else read(key)
    )
    assert "<strong>3</strong>" in client.get(url, base).content.decode()  # the same answer
    assert "<strong>2</strong>" in client.get(url, base | {"node": "Node1"}).content.decode()
    rows = client.get(
        f"/jobs/{job.id}/ancestors/download/", base | {"what": "posteriors"}
    ).content.decode()
    assert {r.split("\t")[0] for r in rows.splitlines()[1:]} == {"Node2"}


def test_the_posterior_index_is_described_and_travels_in_the_zip(
    client: Client, fake_storage: FakeStorage
) -> None:
    from phylobot import bundle, filelist

    job = viewer_job(fake_storage)
    folder = "mafft/iqtree/WAG"
    keys = [f"{folder}/ancestors_posteriors.tsv", f"{folder}/ancestors_posteriors.index.tsv"]
    rows = [row for table in filelist.file_tables(keys) for row in table.rows]
    index = next(text for _k, name, text in rows if name.endswith("index.tsv"))
    assert "byte offset" in index and "ancestors_posteriors.tsv" in index
    fake_storage.objects[f"jobs/{job.id}/results/{keys[1]}"] = b"node\tstart\tlength\n\t0\t5\n"
    zipped = [path for _key, path in bundle.keys(job)]
    assert [p for p in zipped if p.endswith("WAG/ancestors_posteriors.index.tsv")]


def test_ancestor_note_when_the_two_trees_disagree(
    client: Client, fake_storage: FakeStorage
) -> None:
    """IQ-TREE groups a with b; GeneRax's reconciled tree puts them on both sides of a
    duplication at its root. The note reports the disagreement and the AU test's view,
    from either side, and does not claim which tree is right."""
    job = viewer_job(fake_storage)
    prefix = f"jobs/{job.id}/results/"
    header = "node\tevent\tspecies\tlosses\tlost_in\n"
    leaves = "a\tleaf\tA\t0\t\nb\tleaf\tB\t0\t\nc<i>\tleaf\tC\t0\t\n"
    gx = "mafft/generax/WAG"
    summary = json.loads(fake_storage.objects[prefix + "summary.json"])
    summary["combinations"].append({**summary["combinations"][1], "key": gx, "tree": "generax"})
    fake_storage.objects[prefix + "summary.json"] = json.dumps(summary).encode()
    for name in ("alignment.fasta", "ancestors.fasta", "ancestors_posteriors.tsv",
                 "ancestors_pp.tsv", "supports.tsv"):  # fmt: skip
        fake_storage.objects[f"{prefix}{gx}/{name}"] = fake_storage.objects[
            f"{prefix}mafft/iqtree/WAG/{name}"
        ]
    files = {
        "mafft/iqtree/WAG/reconciliation.tsv": header + "Node2\tspeciation\tS\t0\t\n"
        "(root)\tspeciation\tS\t0\t\n" + leaves,
        "mafft/iqtree/WAG/reconciliation.nwk": "((a,b)Node2,c<i>);",
        f"{gx}/tree.nwk": "((a,c<i>)Node2,b)Node1;",
        f"{gx}/reconciliation.tsv": header + "Node2\tspeciation\tS\t0\t\n"
        "(root)\tduplication\tS\t0\t\n" + leaves,
        f"{gx}/reconciliation.nwk": "((a,c<i>)Node2,b);",
        "mafft/topology_tests/WAG.tsv": "tree\tlogL\tdeltaL\tbp_rell\tp_kh\tp_sh\tc_elw\tp_au\n"
        "iqtree\t-10\t0\t0.9\t1\t1\t0.9\t1\ngenerax\t-12\t2\t0.1\t0.2\t0.2\t0.1\t0.21\n",
    }  # fmt: skip
    for name, text in files.items():
        fake_storage.objects[prefix + name] = text.encode()
    page = client.get(
        f"/jobs/{job.id}/ancestors/", {"combination": "mafft/iqtree/WAG", "node": "Node2"}
    ).content.decode()
    text = " ".join(html.unescape(page).split())
    assert "The two trees disagree about this ancestor." in text
    assert "form a group in IQ-TREE's tree, but not in GeneRax's reconciled tree" in text
    assert "cannot tell the two trees apart (AU test: IQ-TREE p = 1, GeneRax p = 0.21)" in text
    assert "if that tree is right, this ancestor mixes paralogs" in text
    # and from GeneRax's side, against IQ-TREE's tree
    other = client.get(f"/jobs/{job.id}/ancestors/", {"combination": gx, "node": "Node2"})
    text = " ".join(html.unescape(other.content.decode()).split())
    assert "form a group in GeneRax's reconciled tree, but not in IQ-TREE's tree" in text


def test_tree_viewer_agreement_across_combinations(
    client: Client, fake_storage: FakeStorage
) -> None:
    """Each branch: how many of the other alignment x model combinations (same tree
    source) have it. muscle's tree has the same clades, numbered the other way round."""
    job = viewer_job(fake_storage)
    prefix = f"jobs/{job.id}/results/"
    trees = {
        "mafft/iqtree/modelfinder": "((a,b)Node2,(c<i>,d)Node3)Node1;",
        "mafft/iqtree/WAG": "((a,c<i>)Node2,(b,d)Node3)Node1;",
        "muscle/iqtree/modelfinder": "((a,b)Node1,(c<i>,d)Node3)Node2;",
    }
    for name, newick in trees.items():
        fake_storage.objects[f"{prefix}{name}/tree.nwk"] = newick.encode()
    url = f"/jobs/{job.id}/tree/"
    query = {"combination": "mafft/iqtree/modelfinder", "labels": "agree", "width": "agree"}
    page = html.unescape(client.get(url, query).content.decode())
    assert "Agreement across combos" in page and "Agreement (k of n)" in page
    assert re.search(r">Node2 <tspan[^>]*>1/2<", page)  # muscle has (a,b); mafft WAG not
    assert not re.search(r">Node1 <tspan", page)  # the root branch: no split to compare
    assert "k of the 2 other alignment × model combinations" in " ".join(page.split())
    assert "share of the other combinations with this branch" in page


def test_stop_job_from_the_job_page(client: Client, fake_storage: FakeStorage) -> None:
    running = Job.objects.create(owner=local_user(), name="long one")
    assert jobs.mark(running.id, Job.Status.RUNNING)

    # The affordance is offered while it runs, and says it cannot be resumed afterwards.
    page = client.get(f"/jobs/{running.id}/").content.decode()
    assert f'action="/jobs/{running.id}/stop/"' in page
    assert "Stop this job" in page and "cannot be resumed" in page

    response = client.post(f"/jobs/{running.id}/stop/", follow=True)
    assert "Stopping the job" in response.content.decode()
    running.refresh_from_db()
    assert running.status == Job.Status.CANCELLED
    assert running.error == jobs.STOPPED_BY_OWNER
    assert running.deleted_at is None  # stopping is not deleting: the files stay

    # Gone once the job has finished, and a stale click does nothing.
    page = client.get(f"/jobs/{running.id}/").content.decode()
    assert "Stop this job" not in page
    assert client.post(f"/jobs/{running.id}/stop/").status_code == 302
    running.refresh_from_db()
    assert running.status == Job.Status.CANCELLED


def test_stop_job_is_post_only_and_private(client: Client, fake_storage: FakeStorage) -> None:
    theirs = other_users_job()
    assert jobs.mark(theirs.id, Job.Status.RUNNING)
    assert client.post(f"/jobs/{theirs.id}/stop/").status_code == 404
    theirs.refresh_from_db()
    assert theirs.status == Job.Status.RUNNING  # another user's job is untouched

    mine = Job.objects.create(owner=local_user())
    assert jobs.mark(mine.id, Job.Status.RUNNING)
    assert client.get(f"/jobs/{mine.id}/stop/").status_code == 405  # never from a GET
    mine.refresh_from_db()
    assert mine.status == Job.Status.RUNNING


def test_my_jobs_offers_stop_only_for_unfinished_jobs(
    client: Client, fake_storage: FakeStorage
) -> None:
    owner = local_user()
    queued = Job.objects.create(owner=owner, name="waiting")
    running = Job.objects.create(owner=owner, name="going")
    done = Job.objects.create(owner=owner, name="finished")
    assert jobs.mark(running.id, Job.Status.RUNNING)
    assert jobs.mark(done.id, Job.Status.RUNNING) and jobs.mark(done.id, Job.Status.SUCCEEDED)

    listing = client.get("/jobs/").content.decode()
    for job in (queued, running):
        assert f'id="stop-{job.id}"' in listing and f'data-stop="{job.id}"' in listing
        assert f'action="/jobs/{job.id}/stop/"' in listing
    assert f'data-stop="{done.id}"' not in listing
    assert f'id="stop-{done.id}"' not in listing

    # It posts back to My jobs rather than the job page (safe_next).
    response = client.post(f"/jobs/{running.id}/stop/", {"next": "/jobs/"}, follow=True)
    assert response.redirect_chain[-1][0] == "/jobs/"
    running.refresh_from_db()
    assert running.status == Job.Status.CANCELLED
