"""The ESR page and its jobs (phylobot.esr_views): the form, from a job's combination and
from uploads, what it stages and queues, the job type's place in the step graph, the
work estimate and My usage, and the job page's ESR section."""

import json
import math
import re
from pathlib import Path
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import Client

from phylobot import (
    esr,
    esr_lengths,
    esr_views,
    extend,
    filelist,
    job_types,
    jobs,
    stages,
    usage,
    workload,
)
from phylobot.executors import PARAM_ALLOWLIST, build_params, expected_tasks
from phylobot.models import Job
from phylobot.modes import local_user
from tests.conftest import FakeStorage

ROOT = Path(__file__).resolve().parent.parent
COMBO = "mafft/iqtree/LG+G4"
CODON = "mafft/iqtree/GY+F3X4+G4"
ALIGNMENT = ">a\nMKLV-\n>b\nMKLAA\n>c\nMRLWA\n>d\nMRLW-\n>e\nMKIW-\n"
TREE = "(a:0.1,b:0.2,((c:0.3,d:0.4)Node3:0.3,e:0.1)Node4:0.2)Node1;"


def finished_job(fake_storage: FakeStorage, owner: User | None = None, **params: Any) -> Job:
    """A finished analysis: an amino-acid combination and a codon one."""
    job = Job.objects.create(owner=owner or local_user(), name="ime2", params=params)
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combos = [
        {"key": COMBO, "label": "MAFFT · IQ-TREE · LG+G4", "data": "amino",
         "fit": {"model": "LG+G4"}},
        {"key": CODON, "label": "MAFFT · IQ-TREE · GY", "data": "codon",
         "fit": {"model": "GY+F3X4+G4"}},
    ]  # fmt: skip
    files = {"summary.json": json.dumps({"combinations": combos})}
    for key, model in ((COMBO, "LG+G4"), (CODON, "GY+F3X4+G4")):
        files |= {f"{key}/alignment.fasta": ALIGNMENT, f"{key}/tree.nwk": TREE,
                  f"{key}/model.txt": model + "\n"}  # fmt: skip
    for name, text in files.items():
        fake_storage.objects[prefix + name] = text.encode()
    return job


def form(**changes: Any) -> dict[str, Any]:
    data = {"name": "esr test", "seed": 7, "source": "job", "tips": "all", "sample": 20,
            "matrix": "LG", "variant": "+G4", "e0": "on", "e1": "on"}  # fmt: skip
    return data | changes


def staged(fake_storage: FakeStorage, job: Job) -> dict[str, str]:
    inputs = f"jobs/{job.id}/inputs/esr/"
    return {k[len(inputs) :]: v.decode() for k, v in fake_storage.objects.items()
            if k.startswith(inputs)}  # fmt: skip


def esr_job() -> Job:
    return Job.objects.filter(params__esr=True).latest("created_at")


@pytest.mark.django_db
def test_page_and_nav(client: Client) -> None:
    page = client.get("/esr/").content.decode()
    assert "<h1>Extant Sequence Reconstruction</h1>" in page
    # In the bar it is under Validation, by its full name, and that menu is marked current.
    assert '<summary class="current">Validation</summary>' in page
    assert 'aria-current="page" title="How well do ancestors come out' in page
    assert ">Extant Sequence Reconstruction</a>" in page
    assert "My ESR Jobs" in page and "None yet." in page
    # The reference is a link, and the page has an "On this page" sidebar.
    assert 'href="https://doi.org/10.1007/s00239-024-10162-3"' in page
    assert 'data-toc aria-label="On this page"' in page and 'href="#my-esr"' in page


@pytest.mark.django_db
def test_on_a_combination(client: Client, fake_storage: FakeStorage) -> None:
    source = finished_job(fake_storage)
    options = client.get("/esr/combinations/", {"job": str(source.id)}).content.decode()
    assert f'value="{COMBO}" selected' in options
    assert CODON not in options  # codon models are not handled yet

    response = client.post("/esr/", form(job=source.id, combination=COMBO))
    job = esr_job()
    assert response.status_code == 302 and response["Location"] == f"/jobs/{job.id}/"
    files = staged(fake_storage, job)
    assert set(files) == {"alignment.fasta", "tree.nwk", "model.txt", "settings.json"}
    assert files["model.txt"] == "LG+G4\n" and files["tree.nwk"] == TREE + "\n"
    settings = json.loads(files["settings.json"])
    assert settings == {"levels": ["e0", "e1"], "tips": ["a", "b", "c", "d", "e"], "seed": 7,
                        "fit_lengths": False}  # fmt: skip
    assert job.params["esr"] is True and job.params["esr_of"] == {
        "kind": "job", "job": str(source.id), "combination": COMBO,
        "label": "ime2 · MAFFT · IQ-TREE · LG+G4", "model": "LG+G4", "levels": ["e0", "e1"],
        "tips": 5, "tip_choice": "all", "fit_lengths": False, "sequences": 5, "columns": 5,
    }  # fmt: skip
    assert "simulated" not in job.params
    assert job.work_units and job.work_units > 0
    assert [t.key for t in job_types.types_of(job.params)] == ["esr"]
    assert usage.kind_of(job.params) == "esr" and extend.refusal(job)
    plan = esr_views.describe(job.params)
    assert plan[0].startswith("fit the alignment, tree and model of ime2")
    assert "hold out 5 sequences" in plan[1] and "(E0 and E1)" in plan[1]


@pytest.mark.django_db
def test_simulated_data_stays_simulated(client: Client, fake_storage: FakeStorage) -> None:
    source = finished_job(fake_storage, simulated=True)
    assert client.post("/esr/", form(job=source.id, combination=COMBO)).status_code == 302
    assert esr_job().params["simulated"] is True


@pytest.mark.django_db
def test_someone_elses_job_is_refused(client: Client, fake_storage: FakeStorage) -> None:
    other = User.objects.create(username="other")
    source = finished_job(fake_storage, owner=other)
    page = client.post("/esr/", form(job=source.id, combination=COMBO)).content.decode()
    assert "Choose one of your finished jobs." in page
    assert not Job.objects.filter(params__esr=True).exists()


@pytest.mark.django_db
def test_a_codon_combination_is_refused(client: Client, fake_storage: FakeStorage) -> None:
    source = finished_job(fake_storage)
    page = client.post("/esr/", form(job=source.id, combination=CODON)).content.decode()
    assert "Choose one of that job&#x27;s combinations." in page


@pytest.mark.django_db
def test_which_sequences(client: Client, fake_storage: FakeStorage) -> None:
    source = finished_job(fake_storage)
    base = {"job": source.id, "combination": COMBO}
    client.post("/esr/", form(**base, tips="sample", sample=2))
    first = json.loads(staged(fake_storage, esr_job())["settings.json"])["tips"]
    client.post("/esr/", form(**base, tips="sample", sample=2))
    again = json.loads(staged(fake_storage, esr_job())["settings.json"])["tips"]
    assert len(first) == 2 and first == again  # drawn with the seed
    client.post("/esr/", form(**base, tips="list", tip_names="c\n\nd\nc\n", e1=""))
    job = esr_job()
    assert json.loads(staged(fake_storage, job)["settings.json"])["tips"] == ["c", "d"]
    assert job.params["esr_of"]["levels"] == ["e0"] and job.params["esr_of"]["tips"] == 2
    page = client.post("/esr/", form(**base, tips="list", tip_names="c\nzz\n")).content.decode()
    assert "Not in the alignment: zz." in page
    page = client.post("/esr/", form(**base, e0="", e1="")).content.decode()
    assert "Choose at least one level" in page


def upload(**changes: Any) -> dict[str, Any]:
    files = {
        "alignment_file": SimpleUploadedFile("mine.fasta", ALIGNMENT.encode()),
        "tree_file": SimpleUploadedFile("mine.nwk", b"((a,b),(c,d),e);"),
    }
    return form(source="upload", matrix="WAG", variant="+I+G4") | files | changes


@pytest.mark.django_db
def test_uploads(client: Client, fake_storage: FakeStorage) -> None:
    assert client.post("/esr/", upload()).status_code == 302
    job = esr_job()
    files = staged(fake_storage, job)
    assert files["model.txt"] == "WAG+I+G4\n"
    assert json.loads(files["settings.json"])["fit_lengths"] is True
    assert files["alignment.fasta"] == ALIGNMENT
    assert job.params["esr_of"]["kind"] == "upload"
    assert "uploaded alignment and tree under WAG+I+G4" in esr_views.describe(job.params)[0]


DNA = b">a\nACGT\n>b\nACGT\n>c\nACGA\n>d\nACGA\n>e\nACGG\n"
RAGGED = b">a\nMK\n>b\nMKL\n>c\nM\n>d\nM\n>e\nM\n"


def tree_file(text: bytes) -> dict[str, Any]:
    return {"tree_file": SimpleUploadedFile("t.nwk", text)}


def alignment_file(text: bytes) -> dict[str, Any]:
    return {"alignment_file": SimpleUploadedFile("a.fasta", text)}


@pytest.mark.django_db
@pytest.mark.parametrize(
    ("changes", "message"),
    [
        (tree_file(b"((a,b),(c,x),e);"), "the tree and the alignment name different sequences; "
         "not in the tree: d; not in the alignment: x"),
        (tree_file(b"((a,b,c),d,e);"), "needs a binary tree"),
        (tree_file(b"((a,b),(c"), "not a Newick tree"),
        (alignment_file(RAGGED), "not aligned"),
        (alignment_file(b">a\nMK\n>a\nMK\n"), "given twice"),
        (alignment_file(DNA), "protein alignments only"),
        (alignment_file(b"\xff\xfe"), "not a text file"),
        ({"variant": "+C20"}, "Select a valid choice"),
    ],
)  # fmt: skip
def test_bad_uploads_are_refused(
    client: Client, fake_storage: FakeStorage, changes: dict[str, Any], message: str
) -> None:
    page = client.post("/esr/", upload(**changes)).content.decode().replace("&#x27;", "'")
    assert message in page
    assert not Job.objects.filter(params__esr=True).exists()


@pytest.mark.django_db
def test_upload_limits(client: Client, fake_storage: FakeStorage, settings: Any) -> None:
    from phylobot.models import UserLimits

    UserLimits.objects.create(user=local_user(), max_sequences=4)
    page = client.post("/esr/", upload()).content.decode()
    assert "The alignment has 5 sequences; your limit is 4." in page


@pytest.mark.django_db
def test_counts_and_costs_agree(tmp_path: Path) -> None:
    """Tasks: ESR_FIT, an ESR_TIP per held-out sequence, ESR_COLLECT; the step graph,
    the work estimate and the worker's own count (from the staged settings, since only
    `esr` passes the allowlist) agree."""
    params = {"esr": True, "seed": 3, "esr_of": {"tips": 4, "levels": ["e0", "e1"],
                                                  "model": "LG+G4"}}  # fmt: skip
    assert expected_tasks(params) == 6 == sum(stages.expected_kinds(params).values())
    assert [k for k, _ in workload.tasks(params)] == ["ESR_FIT", *["ESR_TIP"] * 4, "ESR_COLLECT"]
    e0_only = workload.estimate(params | {"esr_of": {"tips": 4, "levels": ["e0"]}}).units
    assert 0 < e0_only < workload.estimate(params).units
    (tmp_path / "esr").mkdir()
    (tmp_path / "esr" / "settings.json").write_text(json.dumps({"tips": ["a", "b", "c"]}))
    job = Job.objects.create(owner=local_user(), params=params)
    built = build_params(job, tmp_path, tmp_path / "out")
    assert built["esr"] is True and built["esr_dir"] == str(tmp_path / "esr")
    assert "esr_of" not in built and expected_tasks(built) == 5
    assert "esr" in PARAM_ALLOWLIST


def test_job_type_filters_keep_esr_apart() -> None:
    assert job_types.BY_KEY["esr"].label == "ESR"
    assert usage.KINDS["esr"][1] not in {slot for k, (_, slot) in usage.KINDS.items() if k != "esr"}
    kinds = list(usage.KINDS)  # slot 8 is too close to slot 5 in the dark skin (site.css)
    assert abs(kinds.index("esr") - kinds.index("simulation")) > 1


def test_pipeline_wiring() -> None:
    main = (ROOT / "pipeline/main.nf").read_text()
    config = (ROOT / "pipeline/nextflow.config").read_text()
    assert "workflow ESR {" in main and "} else if (params.esr) {" in main
    for name in ("esr = false", "esr_dir = null"):
        assert name in config
    # ESR's IQ-TREE is the pipeline's own pinned one, copied into phylobot-tools.
    iqtree = config.split("iqtree = '", 1)[1].split("'", 1)[0]
    assert f"FROM {iqtree} AS iqtree" in (ROOT / "docker/tools.Dockerfile").read_text()


def test_files_are_described() -> None:
    keys = ["esr/esr_tips.tsv", "esr/esr_sites.tsv", "esr/esr_summary.json", "esr/tree.nwk"]
    (table,) = filelist.file_tables(keys)
    assert table.title == "ESR" and all(description for _, _, description in table.rows)


SUMMARY = {
    "levels": {
        "e0": {"tips": 2, "sites": 10, "accuracy": 0.6, "mean_pp": 0.7, "log_score": -1.2,
               "bins": [{"low": i / 10, "high": (i + 1) / 10, "sites": 5 if i in (5, 8) else 0,
                         "correct": 3, "pp": 0.55 if i == 5 else 0.85,
                         "accuracy": 0.6 if i in (5, 8) else None} for i in range(10)]},
    },
    "fit": {"model": "LG+G4", "fixed": "LG+G4{0.5}", "log_likelihood": -123.456},
}  # fmt: skip
TIPS = (
    "tip\tlevel\tpendant\tsites\taccuracy\tmean_pp\tlog_score\tlog_likelihood\tseconds\n"
    "a\te0\t0.1\t5\t0.6000\t0.7000\t-1.2\t-100\t0.1\n"
    "b\te0\t0.2\t5\t0.6000\t0.7000\t-1.2\t-100\t0.1\n"
)


@pytest.mark.django_db
def test_job_page(client: Client, fake_storage: FakeStorage) -> None:
    source = finished_job(fake_storage)
    client.post("/esr/", form(job=source.id, combination=COMBO, e1=""))
    job = esr_job()
    running = client.get(f"/jobs/{job.id}/").content.decode()
    assert "ESR on 5 sequences" in running and "Duplicate and run again" not in running
    assert f'href="/jobs/{source.id}/"' in running and "ESR Results" in running
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/esr/"
    fake_storage.put_bytes(prefix + "esr_summary.json", json.dumps(SUMMARY).encode())
    fake_storage.put_bytes(prefix + "esr_tips.tsv", TIPS.encode())
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert 'href="#esr-result"' in page and "60%" in page and "+0.100" in page
    # Colours are attributes, so a downloaded figure keeps them: two non-empty bins.
    assert page.count('r="5" fill="#2a78d6" stroke="#ffffff"') == 2
    assert 'data-plot-download="png"' in page and 'data-plot-download="svg"' in page
    assert 'data-filename="esr.' in page and "phylobot/plots.js" in page
    assert 'popovertarget="help-esr-levels"' in page and "sequence-wise ESR" in page
    for heading in ("Reliability", "Each sequence held out", "Downloads"):
        assert f"<h3>{heading}</h3>" in page
    assert (
        f"/jobs/{job.id}/esr/sequences.csv" in page and f"/jobs/{job.id}/esr/sequences.xlsx" in page
    )
    assert "right 60% of 5 sites" in page and "-123.46" in page
    row = '<td class="label">a</td><td>0.1</td><td>5</td><td>0.6000</td><td>0.7000</td>'
    assert row in re.sub(r">\s+<", "><", page)
    assert f"/jobs/{job.id}/files/esr/esr_tips.tsv" in page


@pytest.mark.django_db
def test_an_imported_esr_job_without_its_source(client: Client) -> None:
    """job_import moves an unknown source's id aside (LINKS); the page names it, unlinked."""
    from phylobot import job_import

    params = job_import._relink(local_user(), {"esr": True, "esr_of": {
        "kind": "job", "job": "0" * 32, "label": "gone · MAFFT", "tips": 3}})  # fmt: skip
    assert "job" not in params["esr_of"]
    job = Job.objects.create(owner=local_user(), params=params)
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "gone · MAFFT (not on this server)" in page


@pytest.mark.django_db
def test_sequences_export(client: Client, fake_storage: FakeStorage) -> None:
    import io

    from openpyxl import load_workbook

    job = Job.objects.create(owner=local_user(), name="gfp", params={"esr": True, "esr_of": {
        "label": "gfp · MAFFT", "model": "LG+G4"}})  # fmt: skip
    url = f"/jobs/{job.id}/esr/sequences"
    assert client.get(url + ".csv").status_code == 404  # not finished
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    two_levels = TIPS + "a\te1\t0.1\t5\t0.4000\t0.6000\t-1.5\t-101\t0.2\n"
    fake_storage.put_bytes(f"jobs/{job.id}/results/esr/esr_tips.tsv", two_levels.encode())
    response = client.get(url + ".csv")
    assert response.status_code == 200 and "_esr_sequences.csv" in response["Content-Disposition"]
    assert response.content.decode().splitlines() == [
        "Sequence,Pendant branch length,Sites,E0 accuracy,E0 mean posterior,E0 log score,"
        "E1 accuracy,E1 mean posterior,E1 log score",
        "a,0.1,5,0.6000,0.7000,-1.2,0.4000,0.6000,-1.5",
        "b,0.2,5,0.6000,0.7000,-1.2,,,",
    ]
    book = load_workbook(io.BytesIO(client.get(url + ".xlsx").content))
    rows = list(book["Sequences"].iter_rows(values_only=True))
    assert rows[1][:4] == ("a", 0.1, 5, 0.6) and book["About"]["B4"].value == "gfp · MAFFT"
    assert client.get(url + ".tsv").status_code == 404
    other = Job.objects.create(owner=local_user(), params={})
    assert jobs.mark(other.id, Job.Status.RUNNING) and jobs.mark(other.id, Job.Status.SUCCEEDED)
    assert client.get(f"/jobs/{other.id}/esr/sequences.csv").status_code == 404


# --- by branch length (esr_lengths) ----------------------------------------------------

SITE_HEADER = "tip\tlevel\tsite\ttrue\tml\tpp_ml\tpp_true\t" + "\t".join(esr.AMINO_ACIDS)


def site_line(tip: str, level: str, right: bool, pp: float, pp_true: float) -> str:
    cells = [tip, level, "1", "A", "A" if right else "K", pp, pp if right else pp_true]
    return "\t".join(str(c) for c in cells + [0] * 20)


def eight_tips() -> tuple[str, str]:
    """Eight sequences, pendants 0.01 ... 0.08 (h the longest, 0 for g's zero length); each
    with four sites per level: E0 three right at 0.95 and one wrong (at 0.95 for the
    shorter half, 0.5 with the truth at 0.3 for the longer half); E1 one more wrong."""
    sites, tips = [SITE_HEADER], [TIPS.splitlines()[0]]
    for i, tip in enumerate("abcdefgh"):
        pendant = 0 if tip == "g" else (i + 1) / 100
        for level in ("e0", "e1"):
            tips.append(f"{tip}\t{level}\t{pendant}\t4\t0\t0\t0\t0\t0")
            for k in range(4):
                wrong = k == 3 or (level == "e1" and k == 2)
                if wrong:
                    sites.append(site_line(tip, level, False, 0.95 if i < 4 else 0.5, 0.3))
                else:
                    sites.append(site_line(tip, level, True, 0.95, 0.95))
    return "\n".join(sites) + "\n", "\n".join(tips) + "\n"


def test_tallies_and_quarters() -> None:
    sites, tips = eight_tips()
    tallies = esr_lengths.tally_sites(sites.splitlines())
    a = tallies[("a", "e0")]
    assert (a.sites, a.correct, a.confident_wrong, a.caught) == (4, 3, 1, 1)
    assert a.ece == pytest.approx(abs(3 - 4 * 0.95) / 4)
    # A printed 0 for the true residue is "below 0.000005", not minus infinity.
    zero = esr_lengths.tally_sites([SITE_HEADER, site_line("z", "e0", False, 1.0, 0.0)])
    assert zero[("z", "e0")].log_score == pytest.approx(math.log(esr.LOG_FLOOR))
    lengths = esr_lengths.pendants(tips)
    assert lengths["g"] == 0 and lengths["h"] == 0.08
    rows = esr_lengths.quarters(tallies, lengths, ["e0"])
    assert [r.label for r in rows] == [*esr_lengths.QUARTERS, "All"]
    assert [r.t.tips for r in rows] == [2, 2, 2, 2, 8]
    assert (rows[0].low, rows[0].high) == (0, 0.01)  # g's zero length sorts first
    assert rows[-1].t.confident_wrong_per_tip == 0.5 and rows[-1].t.caught_share == 1.0
    # Too few sequences for quarters: only All.
    few = {k: v for k, v in tallies.items() if k[0] in "abc"}
    assert [r.label for r in esr_lengths.quarters(few, lengths, ["e0"])] == ["All"]
    # A file without the columns gives nothing, never an error.
    assert esr_lengths.tally_sites(["tip\tlevel", "a\te0"]) == {}


def test_axes() -> None:
    lo, hi, ticks = esr_lengths._log_axis([0.0, 0.002, 2.08])
    assert 10**lo == pytest.approx(0.001) and ticks[-1] == 10
    assert esr_lengths._symmetric_axis([0.012, -0.03], 0.05) == (0.05, [-0.05, 0, 0.05])
    m, ticks = esr_lengths._symmetric_axis([0.017], 0.005)
    assert m == pytest.approx(0.02) and len(ticks) == 5


@pytest.mark.django_db
def test_by_branch_length_on_the_job_page(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user(), name="esr", params={"esr": True, "esr_of": {}})
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/esr/"
    summary = json.loads(json.dumps(SUMMARY))
    summary["levels"]["e1"] = summary["levels"]["e0"]
    sites, tips = eight_tips()
    fake_storage.put_bytes(prefix + "esr_summary.json", json.dumps(summary).encode())
    fake_storage.put_bytes(prefix + "esr_tips.tsv", tips.encode())
    without = client.get(f"/jobs/{job.id}/").content.decode()
    assert "By branch length" not in without and "AltAll catches" not in without
    fake_storage.put_bytes(prefix + "esr_sites.tsv", sites.encode())
    page = re.sub(r">\s+<", "><", client.get(f"/jobs/{job.id}/").content.decode())
    assert "<h3>By branch length</h3>" in page and "<h4>Gap per sequence</h4>" in page
    assert "<h4>E1 against E0</h4>" in page and "Longest quarter" in page
    assert "<dt" in page and ">ECE</dt><dd>" in page and "<dd>100%</dd>" in page
    # Eight dots per level on the gap scatter, eight on E1 - E0; the colours as attributes.
    assert page.count('r="3.5" fill="#2a78d6"') == 8 and page.count('r="3.5" fill="#eb6834"') == 16
    assert 'data-filename="esr.' in page and '.gap"' in page and '.difference"' in page
    assert "E1 − E0 log score" in page


@pytest.mark.django_db
def test_my_esr_jobs_sort_and_search(client: Client) -> None:
    params = {"esr": True, "esr_of": {"tips": 3}}
    Job.objects.create(owner=local_user(), name="esr A", params=params)
    page = client.get("/esr/").content.decode()
    assert 'id="esr-table"' in page and 'data-table-search="esr-table"' in page
    assert '<button type="button" class="sort">Sequences held out</button>' in page
