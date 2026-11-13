"""The branch-site test (phylobot.branch_site): foreground marks, codeml's tree and
control files, output parsing, the export, validation, and the follow-up job's pages."""

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest
from django.test import Client

from phylobot import branch_site, jobs, stages, workload
from phylobot.executors import build_params, expected_tasks
from phylobot.models import Job
from phylobot.modes import local_user
from phylobot.pipeline.branch_site_export import export
from tests.conftest import FakeStorage

FIXTURES = Path(__file__).parent / "fixtures" / "codeml"
UNROOTED = "((a:0.1,b:0.2)Node2:0.3,c:0.4,d:0.5)Node1;"
ROOTED = "((a:0.1,b:0.2)Node2:0.15,(c:0.4,d:0.5)Node3:0.15);"  # the outgroup root


def test_marks_fold_into_labels_and_come_back_out() -> None:
    text = branch_site.encode_marks("((a #1,'b x'#1)Node2$1,c,d);")
    assert "#" not in text and "$" not in text
    assert "'b x__phylobot_fg_branch__'" in text  # inside the quotes
    assert branch_site.strip_mark("a__phylobot_fg_branch__") == ("a", "branch")
    assert branch_site.strip_mark("Node2__phylobot_fg_clade__") == ("Node2", "clade")
    assert branch_site.strip_mark("__phylobot_fg_branch__") == (None, "branch")
    assert branch_site.strip_mark("c") == ("c", None)
    for bad, message in (
        ("((a#2,b),c,d);", r"only #1 and \$1"),
        ("((a#,b),c,d);", "only allowed in codeml's"),
        ("((a__phylobot_fg_clade__,b),c,d);", "reserved name"),
    ):
        with pytest.raises(branch_site.MarkError, match=message):
            branch_site.encode_marks(bad)


def test_foreground_from_the_drawn_tree() -> None:
    clade = branch_site.foreground_splits(UNROOTED, "Node2", "clade")
    assert clade == {frozenset("cd"), frozenset("bcd"), frozenset("b")}  # Node2, a, b
    assert branch_site.foreground_splits(UNROOTED, "c", "branch") == {frozenset("c")}
    # On the outgroup-rooted tree the root's two branches are one unrooted branch.
    assert branch_site.foreground_splits(ROOTED, "Node3", "branch") == {frozenset("cd")}
    assert branch_site.foreground_splits(ROOTED, "Node2", "branch") == {frozenset("cd")}
    with pytest.raises(branch_site.MarkError, match="no branch leading to it"):
        branch_site.foreground_splits(UNROOTED, "Node1", "branch")
    with pytest.raises(branch_site.MarkError, match="every branch is foreground"):
        branch_site.foreground_splits(UNROOTED, "Node1", "clade")
    with pytest.raises(branch_site.MarkError, match="no node 'x'"):
        branch_site.foreground_splits(UNROOTED, "x", "branch")


def test_codeml_tree_is_unrooted_and_marked() -> None:
    clade = branch_site.foreground_splits(UNROOTED, "Node2", "clade")
    assert branch_site.mark_tree(UNROOTED, clade) == "((a #1:0.1,b #1:0.2) #1:0.3,c:0.4,d:0.5);"
    # A rooted tree loses its root; the merged root branch carries the mark.
    marked = branch_site.mark_tree(ROOTED, {frozenset("cd")})
    assert marked.count("#1") == 1 and marked.count("(") == 2
    assert branch_site.mark_tree("((a,b)N,'c d',e);", {frozenset({"c d"})}) == "((a,b),'c d' #1,e);"
    with pytest.raises(branch_site.MarkError, match="not in the tree"):
        branch_site.mark_tree(UNROOTED, {frozenset("ac")})


def test_keys_ignore_lengths_names_and_rooting() -> None:
    codons = b">a\nATG\n"
    tree = branch_site.tree_key(codons, UNROOTED, 1)
    assert tree == branch_site.tree_key(codons, "((b:9,a:1)X:5,d,c);", 1)  # same topology
    assert tree == branch_site.tree_key(codons, ROOTED, 1)  # rooted on the same split
    assert tree != branch_site.tree_key(codons, "((a,c),b,d);", 1)
    assert tree != branch_site.tree_key(b">a\nATA\n", UNROOTED, 1)
    assert tree != branch_site.tree_key(codons, UNROOTED, 2)
    one = branch_site.test_key(tree, {frozenset("cd")})
    assert one == branch_site.test_key(tree, {frozenset("cd")})
    assert one != branch_site.test_key(tree, {frozenset("b")})


def test_control_files(tmp_path: Path) -> None:
    branch_site.write_controls(tmp_path, icode=0)
    assert sorted(str(p.relative_to(tmp_path)) for p in tmp_path.rglob("*.ctl")) == [
        "A/start1.ctl", "A/start2.ctl", "A_null/start1.ctl",
    ]  # fmt: skip
    alternative = (tmp_path / "A" / "start2.ctl").read_text()
    null = (tmp_path / "A_null" / "start1.ctl").read_text()
    for text in (alternative, null):
        assert "model = 2" in text and "NSsites = 2" in text and "RateAncestor = 0" in text
    assert "fix_omega = 0" in alternative and "omega = 4.0" in alternative
    assert "fix_omega = 1" in null and "omega = 1.0" in null


def test_parse_mlc_and_test_on_lysin() -> None:
    # PAML 4.10.10, lysin, foreground the branch to Node10 (2026-10-02).
    a = branch_site.parse_mlc((FIXTURES / "lysin_branch_site_A.mlc").read_text())
    null = branch_site.parse_mlc((FIXTURES / "lysin_branch_site_A_null.mlc").read_text())
    assert (a.lnl, a.parameters, a.kappa) == (-4565.453719, 52, 1.45222)
    assert a.classes.proportions == [0.46056, 0.48391, 0.02708, 0.02845]
    assert a.classes.foreground == [0.12027, 1.0, 998.94025, 998.94025]
    assert a.beb[0] == (97, "D", 0.915) and len(a.beb) == 4  # BEB, not the NEB list
    assert (null.lnl, null.parameters, null.classes.foreground[2]) == (-4567.588984, 51, 1.0)
    result = branch_site.test(a.lnl, null.lnl)
    assert result.statistic == pytest.approx(4.27053, abs=1e-4)
    assert result.p_value == pytest.approx(0.0388, abs=1e-4)  # chi2(1)
    assert result.p_mixture == pytest.approx(result.p_value / 2)
    assert branch_site.test(-10.0, -9.0) == branch_site.Test(0.0, 1.0, 1.0)
    assert math.erfc(math.sqrt(3.841 / 2)) == pytest.approx(0.05, abs=1e-4)
    with pytest.raises(ValueError, match="not from the branch-site model"):
        branch_site.parse_mlc((FIXTURES / "lysin_M0.mlc").read_text())


def test_site_classes_from_rst() -> None:
    # Every site, BEB (not the NEB table before it), from Model A's rst on lysin.
    sites = branch_site.parse_site_classes((FIXTURES / "lysin_branch_site_A.rst").read_text())
    assert len(sites) == 135 and sites[0] == (1, "R", (0.9832, 0.00697, 0.00973, 0.0001))
    assert sites[96] == (97, "D", (0.08215, 0.00316, 0.91362, 0.00106))
    assert sites[-1][:2] == (135, "-")
    with pytest.raises(ValueError, match="no BEB"):
        branch_site.parse_site_classes("Naive Empirical Bayes (NEB) probabilities for 4 classes")


def test_p_values_for_the_page() -> None:
    from phylobot.branch_site_views import p_text

    assert [p_text(p) for p in (0.03878, 0.5, 1.0, 0.0001)] == ["0.0388", "0.5", "1", "0.0001"]
    assert p_text(2.6e-18) == "2.6 × 10⁻¹⁸" and p_text(0.0) == "< 10⁻³⁰⁰"


def test_export(tmp_path: Path) -> None:
    runs = tmp_path / "runs"
    runs.mkdir()
    for name in ("A", "A_null"):
        (runs / f"{name}.start1.mlc").write_text(
            (FIXTURES / f"lysin_branch_site_{name}.mlc").read_text()
        )
        rst = FIXTURES / "lysin_branch_site_A.rst" if name == "A" else None
        (runs / f"{name}.start1.rst").write_text(rst.read_text() if rst else "")
        (runs / f"{name}.start1.ctl").write_text(f"ctl {name}\n")
    tree = tmp_path / "tree.nwk"
    tree.write_text("((s0001 #1:0.1,s0002:0.2) #1:0.3,s0003,s0004);\n")
    alignment, names_file = tmp_path / "codons.fasta", tmp_path / "names.tsv"
    alignment.write_text(">s0001\nATG\n")
    names_file.write_text("id\tname\tdescription\ns0001\tseq one\t\n")
    out = tmp_path / "out"
    names = {"s0001": "seq one", "s0002": "b", "s0003": "c", "s0004": "d"}
    export(tree, runs, out, names, alignment, names_file)
    test_row = (out / "test.tsv").read_text().splitlines()
    assert test_row == ["statistic\tdf\tp_value\tp_mixture", "4.2705\t1\t0.03878\t0.01939"]
    models = (out / "models.tsv").read_text().splitlines()
    assert models[1].startswith("A\t-4565.453719\t52\t1.4522\t0.46056\t")
    assert "\t0.12027\t998.94\t1.5: -4565.4537" in models[1]
    assert models[2].split("\t")[0] == "A_null"
    assert (out / "beb.tsv").read_text().splitlines()[1] == "97\tD\t0.915"
    sites = (out / "sites.tsv").read_text().splitlines()
    assert sites[0] == "site\tresidue\tp0\tp1\tp2a\tp2b\tp_positive" and len(sites) == 136
    assert sites[97] == "97\tD\t0.08215\t0.00316\t0.91362\t0.00106\t0.91468"
    assert (out / "foreground.nwk").read_text() == "(('seq one' #1:0.1,b:0.2) #1:0.3,c,d);\n"
    assert json.loads((out / "summary.json").read_text())["branch_site"] is True
    # codeml's own files, to rerun as they ran.
    codeml = out / "codeml"
    assert sorted(p.name for p in codeml.iterdir()) == [
        "A.start1.ctl", "A.start1.mlc", "A.start1.rst", "A_null.start1.ctl", "A_null.start1.mlc",
        "A_null.start1.rst",
        "README.txt", "aln.fasta", "names.tsv", "tree.nwk",
    ]  # fmt: skip
    assert (codeml / "tree.nwk").read_text().startswith("4 1\n((s0001 #1:0.1,")
    assert (codeml / "A.start1.ctl").read_text() == "ctl A\n"
    (runs / "A_null.start1.mlc").unlink()
    with pytest.raises(ValueError, match="missing"):
        export(tree, runs, out, {})


CODONS = {"seq0": "ATGAAATGGCTG", "seq1": "ATGAAGTGGGTG", "seq2": "ATGCGTTGGCTT",
          "seq3": "ATGAAATTTCTG"}  # fmt: skip


def run_validate(tmp_path: Path, tree: str, **manifest: object) -> str:
    (tmp_path / "alignment1.fasta").write_text("".join(f">{n}\n{s}\n" for n, s in CODONS.items()))
    (tmp_path / "tree1.nwk").write_text(tree)
    content = {"sequences": False, "max_sequences": 500, "max_length": 5000, "max_provided": 5,
               "alignments": [{"filename": "a.fasta"}], "trees": [{"filename": "t.nwk"}],
               "seqtype": "codon", "branch_site": True} | manifest  # fmt: skip
    (tmp_path / "manifest.json").write_text(json.dumps(content))
    proc = subprocess.run(
        [sys.executable, "-m", "phylobot.pipeline.validate", "manifest.json"],
        cwd=tmp_path, capture_output=True, text=True,
    )  # fmt: skip
    return proc.stderr


def test_validate_reads_the_marks(tmp_path: Path) -> None:
    assert run_validate(tmp_path, "((seq0:1,seq1:2)Node2 $1:3,seq2,seq3);") == ""
    tree = (tmp_path / "provided" / "trees" / "provided1.nwk").read_text()
    assert tree == "((s0001 #1:1,s0002 #1:2) #1:3,s0003,s0004);\n"
    assert sorted(p.name for p in (tmp_path / "codeml").rglob("*.ctl")) == [
        "start1.ctl", "start1.ctl", "start2.ctl",
    ]  # fmt: skip
    assert (tmp_path / "provided" / "codon_alignments" / "provided1.fasta").exists()
    # A rooted tree: unrooted for codeml, with a warning when a root branch is marked.
    assert run_validate(tmp_path, "((seq0,seq1)#1,(seq2,seq3));") == ""
    report = json.loads((tmp_path / "validation.json").read_text())
    assert any("rooted" in w for w in report["warnings"])
    for tree, message in (
        ("((seq0,seq1),seq2,seq3);", "no foreground branch is marked"),
        ("((seq0,seq1),seq2,seq3)#1;", "top node has no branch"),
        ("((seq0,seq1),seq2,seq3)$1;", "every branch is foreground"),
        ("((seq0#1,seq1),seq2,seq9);", "names do not match"),
        ("((seq0#2,seq1),seq2,seq3);", "only #1 and $1"),
    ):
        assert message in run_validate(tmp_path, tree), tree
    assert "needs one codon alignment" in run_validate(
        tmp_path, "((seq0#1,seq1),seq2,seq3);", trees=[]
    )
    assert "codon models need coding DNA" in run_validate(
        tmp_path, "((seq0#1,seq1),seq2,seq3);", seqtype="protein"
    )


PARAMS = {"branch_site": True, "seqtype": "codon", "aligners": ["none"],
          "tree_builders": ["none"], "alignment_labels": ["x"], "tree_labels": ["y"]}  # fmt: skip


def test_tasks_and_estimate() -> None:
    assert expected_tasks(PARAMS) == 5  # VALIDATE, 3 x CODEML, BRANCH_SITE_EXPORT
    assert stages.expected_kinds(PARAMS) == {"VALIDATE": 1, "CODEML": 3, "BRANCH_SITE_EXPORT": 1}
    kinds = [kind for kind, _ in workload.tasks(PARAMS)]
    assert sorted(kinds) == sorted(stages.expected_kinds(PARAMS).elements())
    lysin = workload.Size(25, 135)
    seconds = workload.reference_seconds(PARAMS, lysin)
    assert seconds == pytest.approx(2 * 389 + 186 + 0.7 + 0.5, rel=0.01)
    assert workload.estimate(PARAMS, lysin).groups[-1][0] == "codeml"


@pytest.mark.django_db
def test_only_branch_site_reaches_the_pipeline(tmp_path: Path) -> None:
    job = Job(owner=local_user(), params=PARAMS | {"branch_site_of": {"node": "x"}})
    params = build_params(job, tmp_path, tmp_path / "out")
    assert params["branch_site"] is True and "branch_site_of" not in params


def codon_job(fake_storage: FakeStorage, outgroup: bool = False) -> Job:
    """A finished coding-DNA job: a codon combination, tree ((a,b)Node2,c,d)Node1, and an
    amino-acid one with the same topology (other lengths, the clade named Node3)."""
    job = Job.objects.create(
        owner=local_user(), name="lysin", params={"seqtype": "codon", "genetic_code": 1}
    )
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combo = {"key": "mafft/iqtree/modelfinder_codon", "aligner": "mafft", "tree": "iqtree",
             "model": "modelfinder_codon",
             "label": "MAFFT · IQ-TREE · ModelFinder (codon)"}  # fmt: skip
    amino = {"key": "mafft/iqtree/LG", "aligner": "mafft", "tree": "iqtree", "model": "LG",
             "label": "MAFFT · IQ-TREE · LG"}  # fmt: skip
    files = {
        "summary.json": json.dumps({"primary": combo["key"], "combinations": [combo, amino]}),
        "mafft/iqtree/LG/tree.nwk": "(c:0.9,d:0.1,(b:0.5,a:0.5)Node3:0.2)Node1;",
        "mafft/iqtree/LG/ancestors.fasta": ">Node1\nMK\n>Node3\nMK\n",
        "mafft/iqtree/LG/alignment.fasta": "".join(f">{n}\nMK\n" for n in "abcd"),
        "mafft/alignment.cds.fasta": "".join(f">{n}\nATGAAA\n" for n in "abcd"),
        f"{combo['key']}/tree.nwk": UNROOTED,
        f"{combo['key']}/alignment.fasta": "".join(f">{n}\nMK\n" for n in "abcd"),
        f"{combo['key']}/ancestors.fasta": ">Node1\nMK\n>Node2\nMK\n",
    }
    if outgroup:
        files[f"{combo['key']}/tree.rooted.nwk"] = "((a:0.1,b:0.2)Node2:0.15,(c,d)Node1:0.15);"
    for name, text in files.items():
        fake_storage.objects[prefix + name] = text.encode()
    return job


@pytest.mark.django_db
def test_tree_viewer_panel_and_submit(client: Client, fake_storage: FakeStorage) -> None:
    job = codon_job(fake_storage)
    url = f"/jobs/{job.id}/tree/"
    page = client.get(url).content.decode()
    assert "Test for Positive Selection (Branch-Site)" in page
    assert "Choose the\n    foreground from your hypothesis" in page
    assert "so the model above does not\n    affect it" in page
    assert '<option value="Node2">Node2 · 2 sequences</option>' in page
    panel = page.split("branch-site-section", 1)[1]
    assert '<option value="Node1"' not in panel  # the top node: no branch, clade = all
    assert 'id="branch-site-form"' in page and 'stroke="#8e44ad"' not in page
    page = client.get(url, {"fg": "Node2", "fg_scope": "clade"}).content.decode()
    assert page.count('stroke="#8e44ad"') == 3 and "3 branches as drawn" in page
    assert '<input type="hidden" name="fg" value="Node2" form="branch-site-form">' in page

    response = client.post(
        f"/jobs/{job.id}/branch-site/",
        {"combination": "mafft/iqtree/modelfinder_codon", "fg": "Node2", "fg_scope": "clade"},
    )
    test = Job.objects.exclude(id=job.id).get()
    assert response.status_code == 302 and response["Location"] == f"/jobs/{test.id}/"
    assert test.params["branch_site"] and test.params["genetic_code"] == 1
    of = dict(test.params["branch_site_of"])
    tree_key, key = of.pop("tree_key"), of.pop("key")
    assert of == {
        "job": str(job.id), "combination": "mafft/iqtree/modelfinder_codon", "aligner": "mafft",
        "label": "MAFFT · IQ-TREE · ModelFinder (codon)", "node": "Node2", "scope": "clade",
        "branches": 3,
    }  # fmt: skip
    assert len(tree_key) == len(key) == 64 and tree_key != key
    assert test.status == Job.Status.QUEUED and test.work_units
    inputs = f"jobs/{test.id}/inputs/"
    assert fake_storage.objects[inputs + "tree1.nwk"] == (
        b"((a #1:0.1,b #1:0.2) #1:0.3,c:0.4,d:0.5);\n"
    )
    assert fake_storage.objects[inputs + "alignment1.fasta"].startswith(b">a\nATGAAA")
    # The parent's viewer lists it.
    page = client.get(url, {"fg": "Node2"}).content.decode()
    assert "Tests on this tree (from any model" in page and str(test.id) in page

    # The same test again: the panel links to it, and submitting leads there too, also
    # from another model's tree with the same branching order.
    page = client.get(url, {"fg": "Node2", "fg_scope": "clade"}).content.decode()
    assert f'<a class="button primary" href="/jobs/{test.id}/">See the result</a>' in page
    assert "Test for positive selection here</button>" not in page
    assert "this test queued" in page
    lg = {"combination": "mafft/iqtree/LG", "fg": "Node3", "fg_scope": "clade"}
    page = client.get(url, lg).content.decode()
    panel = page.split("branch-site-section", 1)[1]
    assert "See the result" in panel
    assert "another model's tree" in panel
    assert str(test.id) in page.split("Tests on this tree")[1]
    response = client.post(f"/jobs/{job.id}/branch-site/", lg)
    assert response["Location"] == f"/jobs/{test.id}/" and Job.objects.count() == 2
    # A different foreground is a different test; a failed test can be run again.
    page = client.get(url, {"fg": "Node2", "fg_scope": "branch"}).content.decode()
    assert "Test for positive selection here</button>" in page
    assert jobs.mark(test.id, Job.Status.RUNNING) and jobs.mark(test.id, Job.Status.FAILED)
    response = client.post(f"/jobs/{job.id}/branch-site/", lg)
    assert Job.objects.count() == 3 and response["Location"] != f"/jobs/{test.id}/"

    for bad in (
        {"combination": "nope", "fg": "Node2", "fg_scope": "clade"},
        {"combination": "mafft/iqtree/modelfinder_codon", "fg": "Node1", "fg_scope": "clade"},
        {"combination": "mafft/iqtree/modelfinder_codon", "fg": "Node2", "fg_scope": "all"},
    ):
        assert client.post(f"/jobs/{job.id}/branch-site/", bad).status_code == 404
    assert client.get(f"/jobs/{job.id}/branch-site/").status_code == 405


@pytest.mark.django_db
def test_rooted_tree_and_protein_jobs(client: Client, fake_storage: FakeStorage) -> None:
    job = codon_job(fake_storage, outgroup=True)
    page = client.get(f"/jobs/{job.id}/tree/", {"fg": "Node2"}).content.decode()
    # Both of the root's branches: they are one branch for codeml.
    assert page.count('stroke="#8e44ad"') == 2 and "the outgroup root" in page
    job.params = {}  # protein: no codon alignment, no panel
    job.save()
    page = client.get(f"/jobs/{job.id}/tree/").content.decode()
    assert "Test for positive selection" not in page


@pytest.mark.django_db
def test_job_page(client: Client, fake_storage: FakeStorage) -> None:
    parent = codon_job(fake_storage)
    of = {"job": str(parent.id), "combination": "mafft/iqtree/modelfinder_codon",
          "label": "MAFFT · IQ-TREE", "node": "Node2", "scope": "branch",
          "branches": 1}  # fmt: skip
    job = Job.objects.create(owner=local_user(), name="bs", params=PARAMS | {"branch_site_of": of})
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Test for positive selection with codeml&#x27;s branch-site test" in page
    assert "Foreground: the branch leading to Node2 (1 branch)" in page
    assert "Duplicate and run again" not in page and "branch-site test</p>" in page
    assert "fg=Node2&amp;fg_scope=branch" in page
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    files = {
        "summary.json": '{"branch_site": true, "combinations": []}',
        "branch_site/test.tsv": "statistic\tdf\tp_value\tp_mixture\n4.2705\t1\t0.03878\t0.01939\n",
        "branch_site/models.tsv": "model\tlog_likelihood\tfree_parameters\tkappa\tp0\tp1\tp2a\t"
        "p2b\tomega0\tomega2\tstarts\nA\t-4565.45\t52\t1.45\t0.46\t0.48\t0.027\t0.028\t0.12\t"
        "998.94\t1.5: -4565.45\nA_null\t-4567.59\t51\t1.45\t0.43\t0.47\t0.044\t0.048\t0.117\t1\t"
        "1: -4567.59\n",
        "branch_site/beb.tsv": "site\tresidue\tp_positive\n97\tD\t0.915\n12\tK\t0.991\n",
        "branch_site/sites.tsv": "site\tresidue\tp0\tp1\tp2a\tp2b\tp_positive\n"
        + "".join(f"{i}\tK\t0.9\t0.1\t0\t0\t0\n" for i in range(1, 12))
        + "12\tK\t0.005\t0.004\t0.9\t0.091\t0.991\n97\tD\t0.08\t0.005\t0.91\t0.005\t0.915\n",
    }
    for name, text in files.items():
        fake_storage.objects[prefix + name] = text.encode()
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Evidence of positive selection on the foreground" in page
    assert "p = 0.0388" in page and "upper limit (999)" in page
    assert "<li>97 D: 0.915</li>" in page and "<li>12 K: 0.991 **</li>" in page
    assert "<strong>Above 0.95</strong> (1 site; ** above 0.99)" in page
    assert "<strong>12 K</strong> (0.991 **)" in page
    assert 'class="pp-trace"' in page and page.count("<circle cx=") == 1  # one site above 0.95
    assert 'id="bs-sites"' in page and "sites.tsv</a>" in page
    assert "phylobot/ancestors.js" in page
    assert "Branch-site test</h3>" in page  # the file table's group
    fake_storage.objects[prefix + "branch_site/test.tsv"] = (
        b"statistic\tdf\tp_value\tp_mixture\n1.0\t1\t0.3173\t0.1587\n"
    )
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "No significant evidence" in page and "so these are not evidence of" in page
