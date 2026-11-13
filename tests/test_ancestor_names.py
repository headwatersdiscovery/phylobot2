"""Names for ancestors (phylobot.ancestor_names): given once, shown on every combination
that has the ancestor, held to a safe character set, exported and imported with the job."""

from typing import Any

import pytest
from django.contrib.auth.models import User
from django.test import Client

from phylobot import analogs, ancestor_names
from phylobot.models import AncestorName, AncestorNote, Job
from tests.conftest import FakeStorage
from tests.test_extend import split
from tests.test_job_import import download, edit_manifest, refused, upload
from tests.test_views import viewer_job

pytestmark = pytest.mark.django_db

WAG, MUSCLE = "mafft/iqtree/WAG", "muscle/iqtree/modelfinder"
# viewer_job's trees: MAFFT's ((a,b)Node2,c<i>)Node1, MUSCLE's the same numbered the
# other way round, so MAFFT's Node2 is MUSCLE's Node1.


def name(
    client: Client, job: Job, node: str, text: str, combination: str = WAG, **more: str
) -> Any:
    aligner = combination.split("/")[0]
    data = {"combination": combination, "aligner": aligner, "gaps": "original", "node": node,
            "ancestor_name": text, **more}  # fmt: skip
    return client.post(f"/jobs/{job.id}/ancestors/name/", data)


def viewer(client: Client, job: Job, node: str, combination: str = MUSCLE) -> str:
    query = {"combination": combination, "aligner": combination.split("/")[0], "node": node}
    return client.get(f"/jobs/{job.id}/ancestors/", query).content.decode()


@pytest.mark.parametrize(
    ("text", "ok"),
    [
        ("Anc.GFP", True), ("AncMRCA_1-b", True), ("x" * 40, True),
        ("x" * 41, False), ("", False), ("1abc", False), ("Anc GFP", False),
        ("<script>", False), ("Anc|GFP", False), ("Ancé", False), ("Node12", False),
        ("node3", False), ("A", True),
    ],
)  # fmt: skip
def test_names_are_safe_as_fasta_ids_and_file_names(text: str, ok: bool) -> None:
    assert (ancestor_names.problem(text, set()) == "") is ok


def test_a_name_may_not_repeat_a_sequence() -> None:
    assert "sequence names" in ancestor_names.problem("Yarrowia", {"yarrowia"})


def test_named_once_shown_on_every_combination(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    response = name(client, job, "Node2", "Anc.AB")
    assert response.status_code == 302 and "node=Node2" in response["Location"]
    # MUSCLE numbers the same clade Node1: it carries the name there.
    page = viewer(client, job, "Node1")
    assert "<strong>Anc.AB</strong> (Node1)" in page and "Name and Description" in page
    assert ">Anc.AB (Node1) · 2 sequences<" in page  # the node menu
    tree = client.get(
        f"/jobs/{job.id}/tree/", {"combination": MUSCLE, "aligner": "muscle"}
    ).content.decode()
    assert '<tspan font-weight="bold">Anc.AB</tspan> Node1' in tree
    # Downloads carry the name as the FASTA id, the node after it, and in the file name.
    fasta = client.get(
        f"/jobs/{job.id}/ancestors/download/",
        {"combination": MUSCLE, "aligner": "muscle", "node": "Node1", "what": "sequence"},
    )
    assert fasta.content.decode().startswith(">Anc.AB Node1 muscle/iqtree/modelfinder")
    filename = 'filename="Anc.AB.muscle_iqtree_modelfinder.original.fasta"'
    assert filename in fasta["Content-Disposition"]
    # Renamed, then removed: one name per clade, stored once.
    name(client, job, "Node1", "Anc.Renamed", combination=MUSCLE)
    assert list(AncestorName.objects.values_list("name", flat=True)) == ["Anc.Renamed"]
    name(client, job, "Node2", "", remove="1")
    assert not AncestorName.objects.exists()
    assert (
        "<strong>Anc.AB" not in (p := viewer(client, job, "Node1"))
        and "<strong>Node1</strong>" in p
    )


def test_bad_names_are_refused_and_said_why(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    for text, says in (("<b>x</b>", "letters, digits"), ("Node7", "label nodes"),
                       ("A", "sequence names"), ("Anc", None)):  # fmt: skip
        response = name(client, job, "Node2", text)
        page = client.get(response["Location"]).content.decode()
        if says is None:
            assert "is now Anc" in page
        else:
            assert "Not named" in page and says in page
    assert list(AncestorName.objects.values_list("name", flat=True)) == ["Anc"]
    # Unique per job, whatever the case: the root (Node1) cannot also be "ANC".
    response = name(client, job, "Node1", "ANC")
    page = client.get(response["Location"]).content.decode()
    assert "already has that name" in page and "<b>x</b>" not in page


def test_only_the_owner_names_and_only_real_nodes(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = viewer_job(fake_storage)
    assert name(client, job, "Node99", "Anc").status_code == 404
    Job.objects.filter(id=job.id).update(owner=User.objects.create(username="other"))
    assert name(client, job, "Node2", "Anc").status_code == 404
    assert client.get(f"/jobs/{job.id}/ancestors/name/").status_code in (404, 405)
    assert not AncestorName.objects.exists()


def test_an_extension_shares_its_parents_names(client: Client, fake_storage: FakeStorage) -> None:
    parent, child = split(fake_storage)
    name(client, child, "Node2", "Anc.AB")  # named from the extension's viewer
    assert AncestorName.objects.get().job_id == parent.id  # held by the parent
    # MUSCLE is the child's own combination.
    assert "<strong>Anc.AB</strong> (Node1)" in viewer(client, child, "Node1")
    assert "<strong>Anc.AB</strong> (Node2)" in viewer(client, parent, "Node2", WAG)


def test_approximate_where_the_branching_inside_differs() -> None:
    named = analogs.Located("c" * 16, "p" * 16, 4)
    other = analogs.Located("c" * 16, "q" * 16, 4)
    job = Job.objects.create(owner=User.objects.create(username="u"))
    assert ancestor_names.save(job, "c" * 16, "p" * 16, "Anc") == ""
    found = ancestor_names.for_nodes(job, {"Node1": named, "Node9": other}, by_clade=True)
    assert found == {"Node1": ancestor_names.Named("Anc", True),
                     "Node9": ancestor_names.Named("Anc", False)}  # fmt: skip


def test_names_travel_in_the_zip(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    name(client, job, "Node2", "Anc.AB")
    data = download(client, job)
    Job.objects.filter(id=job.id).update(name="original")  # no name clash on import
    assert upload(client, data).status_code == 302
    imported = Job.objects.get(imported_at__isnull=False)
    assert list(imported.ancestor_names.values_list("name", flat=True)) == ["Anc.AB"]


@pytest.mark.parametrize(
    "names",
    [
        "x",
        [{"clade": "0" * 16, "parts": "0" * 16, "name": "<script>"}],
        [{"clade": "zz", "parts": "0" * 16, "name": "Anc"}],
        [{"clade": "0" * 16, "parts": "0" * 16, "name": "Anc"},
         {"clade": "1" * 16, "parts": "0" * 16, "name": "anc"}],
    ],
)  # fmt: skip
def test_bad_names_in_a_zip_are_refused(
    client: Client, fake_storage: FakeStorage, names: Any
) -> None:
    job = viewer_job(fake_storage)
    data = edit_manifest(download(client, job), lambda m: m["job"].update(ancestor_names=names))
    Job.objects.filter(id=job.id).update(name="original")
    assert "ancestor names" in refused(client, data)
    assert not AncestorName.objects.exists()


def test_the_tree_popover_views_and_names(client: Client, fake_storage: FakeStorage) -> None:
    """Owner, 2026-10-30: clicking a node opens its popover, which views the ancestor and
    names it, so naming happens on the tree."""
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/tree/"
    page = client.get(url, {"combination": WAG, "aligner": "mafft"}).content.decode()
    popover = page.split('class="node-popover" data-node="Node2"', 1)[1].split("</div>", 1)[0]
    assert f'data-combination="{WAG}"' in popover and ">View ancestor</a>" in popover
    assert "data-pop-name" in popover and "data-pop-save" in popover
    assert 'name="ancestor_name"' not in popover  # never sent by the tree's GET form
    assert 'id="tree-ancestor-name"' in page  # the POST form the popover goes through
    # Saved from the popover: back to the tree as it was, with the name drawn.
    here = f"{url}?combination={WAG}&aligner=mafft&lengths=cladogram"
    response = name(client, job, "Node2", "Anc.AB", return_to=here)
    assert response["Location"] == here
    page = client.get(response["Location"]).content.decode()
    assert "<strong>Anc.AB</strong> (Node2)" in page and "Node2 is now Anc.AB" in page


@pytest.mark.parametrize(
    "where",
    ["https://evil.example/", "//evil.example/jobs/", "/jobs/", "/jobs/{other}/tree/",
     "/jobs/{job}/ancestors/?x=1", "/jobs/{job}/treeX"],
)  # fmt: skip
def test_it_returns_only_to_this_jobs_tree(
    client: Client, fake_storage: FakeStorage, where: str
) -> None:
    job = viewer_job(fake_storage)
    asked = where.format(job=job.id, other="0" * 8 + "-0000-0000-0000-" + "0" * 12)
    response = name(client, job, "Node2", "Anc.AB", return_to=asked)
    assert response["Location"].startswith(f"/jobs/{job.id}/ancestors/?")


def test_every_ancestral_node_has_a_popover(client: Client, fake_storage: FakeStorage) -> None:
    """Not only those found in other combinations: a job with one combination can still
    name its ancestors from the tree."""
    job = viewer_job(fake_storage)
    page = client.get(
        f"/jobs/{job.id}/tree/", {"combination": MUSCLE, "aligner": "muscle"}
    ).content.decode()
    assert page.count('class="node-popover"') == 2  # Node1 and Node2


def test_the_popover_lists_each_matchs_branch_support(
    client: Client, fake_storage: FakeStorage
) -> None:
    """Owner, 2026-10-30: a column per support measure beside each match, and the node's
    own supports above them. viewer_job's supports.tsv gives every combination's Node2
    UFBoot 88; MUSCLE's match for MAFFT's Node2 is its Node1, which has none."""
    job = viewer_job(fake_storage)
    page = client.get(
        f"/jobs/{job.id}/tree/", {"combination": WAG, "aligner": "mafft"}
    ).content.decode()
    popover = page.split('class="node-popover" data-node="Node2"', 1)[1].split("</div>", 1)[0]
    assert "Branch support: UFBoot 88" in popover
    assert '<th class="num">UFBoot</th>' in popover and "SH-aLRT" not in popover
    rows = popover.split("<tbody>", 1)[1].split("</tr>")
    first = next(r for r in rows if "MODELFINDER" in r.upper() and "Node2" in r)
    assert '<td class="num">88</td>' in first
    muscle = next(r for r in rows if "MUSCLE" in r)
    assert '<td class="num">–</td>' in muscle


def test_a_description_is_shared_by_the_clade(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    # Described without a name: no name appears, the description does, on MUSCLE's Node1.
    name(client, job, "Node2", "", description="Binds DNA.\r\nOrder it first.")
    assert not AncestorName.objects.exists()
    page = viewer(client, job, "Node1")
    assert "Binds DNA.\nOrder it first.</textarea>" in page and "<strong>Node1</strong>" in page
    # Named and described together; clearing the description removes it.
    name(client, job, "Node2", "Anc.AB", description="Binds DNA.")
    assert AncestorName.objects.get().name == "Anc.AB"
    name(client, job, "Node2", "Anc.AB", description="")
    assert AncestorNote.objects.count() == 0 and AncestorName.objects.count() == 1
    # Clearing the name removes it, as "Remove the name" does.
    name(client, job, "Node2", "", description="kept")
    assert not AncestorName.objects.exists() and AncestorNote.objects.get().text == "kept"


def test_descriptions_travel_in_the_zip(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    name(client, job, "Node2", "Anc.AB", description="Binds DNA.")
    data = download(client, job)
    from phylobot import jobs as jobs_mod

    jobs_mod.delete([job])
    upload(client, data)
    imported = Job.objects.get(imported_at__isnull=False)
    assert AncestorNote.objects.get(job=imported).text == "Binds DNA."
