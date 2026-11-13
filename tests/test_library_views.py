"""Ancestral Libraries on the web (phylobot.library_views, step L4): "Add to library" from
the viewers, the top-bar chip, the library page, and that every value is the job's own."""

import json
import re
import uuid

import pytest
from django.contrib.auth import get_user_model
from django.test import Client

from phylobot import (
    ancestors,
    gaps,
    jobs,
    library,
    library_export,
    library_glyphs,
    library_views,
    sites_excel,
)
from phylobot.models import Job, Library, LibraryExport
from phylobot.modes import local_user
from tests.conftest import FakeStorage

KEY = "mafft/iqtree/LG"
TREE = "((a:0.1,b:0.2)Node2:0.3,c:0.4,(d:0.1,e:0.1)Node3:0.5)Node1;"


def finished_job(fake_storage: FakeStorage, name: str = "Ime2 test", **params: object) -> Job:
    """A finished protein job, one combination, tree ((a,b)Node2,c,(d,e)Node3)Node1."""
    job = Job.objects.create(owner=local_user(), name=name, params=dict(params))
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combo = {"key": KEY, "aligner": "mafft", "tree": "iqtree", "model": "LG",
             "label": "MAFFT · IQ-TREE · LG"}  # fmt: skip
    files = {
        "summary.json": json.dumps({"primary": KEY, "combinations": [combo]}),
        f"{KEY}/tree.nwk": TREE,
        f"{KEY}/alignment.fasta": "".join(f">{n}\nMK\n" for n in "abcde"),
        f"{KEY}/ancestors.fasta": ">Node1\nMK\n>Node2\nMK\n>Node3\nMK\n",
        f"{KEY}/ancestors.gaps-fitch.fasta": ">Node1\nM-\n>Node2\nMK\n>Node3\nM-\n",
        # Node2's second site is uncertain: R at 0.35 (AltAll's 0.2 counts it).
        f"{KEY}/ancestors_posteriors.tsv": "node\tsite\tK\tM\tR\n"
        + "".join(
            f"{n}\t1\t0.0\t1.0\t0.0\n{n}\t2\t{k}\t0.0\t{r}\n"
            for n, k, r in [("Node1", 0.9, 0.1), ("Node2", 0.65, 0.35), ("Node3", 0.99, 0.01)]
        ),
    }
    for file, text in files.items():
        fake_storage.objects[prefix + file] = text.encode()
    return job


def uncertain_counts(page: str) -> list[str]:
    """The uncertain-site counts on a library page (each row's cell, then its summary's),
    from the markup only: the page's script names the same attribute."""
    markup = page.rsplit("<script>", 1)[0]
    return [c.split(">", 1)[1].split("<", 1)[0] for c in markup.split('data-count="uncertain"')[1:]]


def add(client: Client, job: Job, *nodes: str, **extra: str):  # type: ignore[no-untyped-def]
    data = {"combination": KEY, "gaps": "original", "node": list(nodes)} | extra
    return client.post(f"/jobs/{job.id}/library/add/", data)


@pytest.mark.django_db
def test_add_from_the_tree_viewer_starts_a_library(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    page = client.get(
        f"/jobs/{job.id}/tree/", {"combination": KEY, "node": "Node2"}
    ).content.decode()
    assert "Add Node2 to a new library" in page
    assert '<input type="hidden" name="node" value="Node2" form="library-add-node">' in page
    assert 'id="library-add-node"' in page and "library-chip" not in page

    back = f"/jobs/{job.id}/tree/?combination={KEY}&node=Node2"
    response = add(client, job, "Node2", next=back)
    assert response.status_code == 302 and response["Location"] == back
    record = Library.objects.get()
    assert (record.name, record.code, record.owner, record.job) == (
        "Library 1",
        "IME2TE",
        local_user(),
        job,
    )
    assert record.spec_version == library.SPEC_VERSION and isinstance(record.spec["seed"], int)
    (target,) = record.targets
    assert (target["combination"], target["node"], target["version"], target["below"]) == (
        KEY,
        "Node2",
        "original",
        2,
    )
    # The partition is the node's identity under any rooting: gaps.node_partitions'.
    expected = {label: parts for parts, label in gaps.node_partitions(TREE).items()}["Node2"]
    assert {frozenset(p) for p in target["partition"]} == expected

    page = client.get(back).content.decode()
    assert ">Add Node2 to</button>" in page
    # The button now offers the job's libraries, the current one first, and a new one.
    menu = page.split('<select name="library" form="library-add-node"', 1)[1].split("</select>", 1)[
        0
    ]
    assert f'<option value="{record.id}" selected>Library 1</option>' in menu
    assert '<option value="new">a new library</option>' in menu


@pytest.mark.django_db
def test_a_partition_names_a_node_where_a_split_names_a_branch(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node1", "Node2")
    partitions = [t["partition"] for t in Library.objects.get().targets]
    # Node1 and Node2 share a branch, so one split would fit both; their partitions differ,
    # and the top node gets one too (a split of it would be empty).
    assert partitions[0] != partitions[1] and all(len(p) == 3 for p in partitions)


@pytest.mark.django_db
def test_adding_a_lineage_and_the_same_node_twice(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node1")
    add(client, job, "Node1", "Node3")  # a lineage overlapping what is already there
    add(client, job, "Node3")  # a second click
    add(client, job, "Node3", gaps="fitch")  # the same node in another reconstruction
    (record,) = Library.objects.all()
    assert [library.target_key(t) for t in record.targets] == [
        (KEY, "Node1", "original"),
        (KEY, "Node3", "original"),
        (KEY, "Node3", "fitch"),
    ]


@pytest.mark.django_db
def test_add_refuses_anything_that_is_not_the_jobs_own(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    someone = get_user_model().objects.create(username="someone")
    theirs = Job.objects.create(owner=someone, name="theirs")
    their_library = Library.objects.create(owner=someone, job=theirs, name="L", code="L")
    cases: list[dict[str, str | list[str]]] = [
        {"combination": "../../../etc"},
        {"combination": "mafft/iqtree/WAG"},  # a real-looking key this job lacks
        {"gaps": "grasp"},  # a version this combination does not have
        {"gaps": "../x"},
        {"node": ["Node9"]},
        {"node": ["a"]},  # an extant sequence has no reconstruction
        {"node": ["<script>"]},
        {"node": []},
        {"library": "not-a-uuid"},
        {"library": str(uuid.uuid4())},
        {"library": str(their_library.id)},
    ]
    for bad in cases:
        data: dict[str, str | list[str]] = {"combination": KEY, "gaps": "original"}
        data |= {"node": ["Node2"]} | bad
        assert client.post(f"/jobs/{job.id}/library/add/", data).status_code == 404, bad
    assert client.post(f"/jobs/{theirs.id}/library/add/", {}).status_code == 404
    assert client.get(f"/jobs/{job.id}/library/add/").status_code == 405
    assert not Library.objects.filter(owner=local_user()).exists()
    assert their_library.targets == []


@pytest.mark.django_db
def test_add_needs_a_finished_job(client: Client, fake_storage: FakeStorage) -> None:
    job = Job.objects.create(owner=local_user(), name="queued")
    assert add(client, job, "Node1").status_code == 404


@pytest.mark.django_db
def test_add_returns_only_to_the_jobs_own_pages(client: Client, fake_storage: FakeStorage) -> None:
    job = finished_job(fake_storage)
    home = f"/jobs/{job.id}/"
    for bad in ["https://evil.example/", "//evil.example/", f"/jobs/{job.id}//evil.example",
                f"/jobs/{job.id}/\\evil", "/jobs/", f"/jobs/{uuid.uuid4()}/tree/", ""]:  # fmt: skip
        assert add(client, job, "Node2", next=bad)["Location"] == home, bad


@pytest.mark.django_db
def test_library_page_rename_remove_new_and_delete(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node1", "Node2")
    record = Library.objects.get()
    url = f"/jobs/{job.id}/library/{record.id}/"
    page = client.get(url).content.decode()
    assert "Ancestors (2)" in page and "MAFFT · IQ-TREE · LG" in page
    assert f"/jobs/{job.id}/ancestors/?aligner=mafft&amp;combination=" in page
    # An "On this page" sidebar, and the two tables fold.
    assert 'data-toc aria-label="On this page"' in page and 'href="#ancestors"' in page
    assert 'aria-controls="ancestor-table-body"' in page
    # Cross-combination disagreement (phylobot.library_disagreement): the Disputed column and,
    # with one combination per clade, the prompt to add more.
    assert ">Disputed " in page and 'data-count="disputed"' in page
    assert "Add the same ancestor from more than one aligner or model" in page

    client.post(url + "rename/", {"name": "  Stability   ladder ", "code": "ime-2!x9"})
    record.refresh_from_db()
    assert (record.name, record.code) == ("Stability ladder", "IME2X9")  # cleaned, not refused
    client.post(url + "rename/", {"name": "", "code": "???"})
    record.refresh_from_db()
    assert (record.name, record.code) == ("Stability ladder", "LIB")
    # The title's pencil posts the name alone; the code and notes stay as they were.
    notes = "Thermostability set.\r\nOrder by Friday."
    client.post(url + "rename/", {"code": "LIB", "notes": notes})
    client.post(url + "rename/", {"name": "Stability ladder"})
    record.refresh_from_db()
    assert (record.name, record.code) == ("Stability ladder", "LIB")
    assert record.notes == "Thermostability set.\nOrder by Friday."
    page = client.get(url).content.decode()
    assert "data-rename" in page and "Code and Description" in page
    assert "Thermostability set.\nOrder by Friday.</textarea>" in page

    response = client.post(url + "remove/", {"combination": KEY, "node": "Node1",
                                             "version": "original", "next": url})  # fmt: skip
    assert response["Location"] == url
    record.refresh_from_db()
    assert [t["node"] for t in record.targets] == ["Node2"]

    # A new library becomes the current one: the viewers' buttons follow it.
    response = client.post(f"/jobs/{job.id}/library/new/")
    second = Library.objects.exclude(id=record.id).get()
    assert response["Location"] == f"/jobs/{job.id}/library/{second.id}/"
    assert second.name == "Library 2"
    tree = client.get(f"/jobs/{job.id}/tree/", {"combination": KEY, "node": "Node3"}).content
    assert f'<option value="{second.id}" selected>Library 2'.encode() in tree

    response = client.post(url + "delete/")
    assert response["Location"] == f"/jobs/{job.id}/"
    record.refresh_from_db()
    assert record.deleted_at is not None  # the record is kept
    assert client.get(url).status_code == 404
    client.post(f"/jobs/{job.id}/library/new/")
    assert Library.objects.filter(name="Library 3").exists()  # names never repeat


@pytest.mark.django_db
def test_others_libraries_are_404(client: Client, fake_storage: FakeStorage) -> None:
    job = finished_job(fake_storage)
    someone = get_user_model().objects.create(username="someone")
    theirs = Library.objects.create(owner=someone, job=job, name="Theirs", code="T")
    url = f"/jobs/{job.id}/library/{theirs.id}/"
    assert client.get(url).status_code == 404
    for action in ("rename/", "remove/", "delete/"):
        assert client.post(url + action, {"name": "mine now"}).status_code == 404
    theirs.refresh_from_db()
    assert theirs.name == "Theirs" and theirs.deleted_at is None
    # Another job's library is not reachable through this job's URL either.
    mine_elsewhere = finished_job(fake_storage, name="other")
    add(client, mine_elsewhere, "Node1")
    elsewhere = Library.objects.get(job=mine_elsewhere)
    assert client.get(f"/jobs/{job.id}/library/{elsewhere.id}/").status_code == 404


@pytest.mark.django_db
def test_job_page_lists_libraries_where_there_can_be_some(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Ancestral Libraries" in page and "the first one starts a library" in page
    add(client, job, "Node2")
    page = client.get(f"/jobs/{job.id}/").content.decode()
    assert "Library 1</a>" in page and "IME2TE · 1 ancestor" in page and "current" in page
    trees_only = finished_job(fake_storage, name="trees", stop_after="trees")
    assert b"Ancestral Libraries" not in client.get(f"/jobs/{trees_only.id}/").content
    queued = Job.objects.create(owner=local_user(), name="queued")
    assert b"Ancestral Libraries" not in client.get(f"/jobs/{queued.id}/").content


@pytest.mark.django_db
def test_the_bar_carries_no_library_chip(client: Client, fake_storage: FakeStorage) -> None:
    """The chip is gone from the bar (2026-11-11): the viewers' "Add to" menus name the
    library, and the job page and Libraries list lead to it."""
    job = finished_job(fake_storage)
    add(client, job, "Node2")
    assert b"library-chip" not in client.get(f"/jobs/{job.id}/").content
    assert b"library-chip" not in client.get("/jobs/").content


@pytest.mark.django_db
def test_library_page_draws_each_ancestor(client: Client, fake_storage: FakeStorage) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node1", "Node2")
    add(client, job, "Node3", gaps="fitch")
    record = Library.objects.get()
    url = f"/jobs/{job.id}/library/{record.id}/"
    page = client.get(url).content.decode()
    # One viewer row per target, each able to open its sequence below.
    assert page.count('"anchors": ["target-') == 1 and "Node3 \\u00b7 fitch" in page
    # Counts per row: only Node2's site at 0.2 (Node1's R is 0.1; Fitch masked Node3's).
    assert uncertain_counts(page) == ["0", "0", "1", "1", "0", "0"]  # each row's cell, summary
    assert 'data-copy="&gt;Node2 MAFFT · IQ-TREE · LG original\nMK\n"' in page
    assert "Node3's sequence: 1 residue, " in page
    # Site maps are fetched when opened, not carried by the page.
    markup = page.rsplit("<script>", 1)[0]
    assert 'class="sm-site' not in markup
    assert markup.count('class="site-map" role="list" data-src=') == 3
    at_low = client.get(url, {"altall": "0.1"}).content.decode()
    assert uncertain_counts(at_low) == ["1", "1", "1", "1", "0", "0"]


@pytest.mark.django_db
def test_library_page_marks_lengths_against_the_extant_sequences(
    client: Client, fake_storage: FakeStorage
) -> None:
    """Both tables flag a length outside the extant range (every extant row is MK): the
    ancestors at 2 residues pass, Fitch's one-residue Node3 is flagged short."""
    job = finished_job(fake_storage)
    add(client, job, "Node1", "Node2")
    add(client, job, "Node3", gaps="fitch")
    record = Library.objects.get()
    page = client.get(f"/jobs/{job.id}/library/{record.id}/").content.decode()
    ancestors_table = page.split('class="jobs library-table"', 1)[1].split("</table>", 1)[0]
    assert ancestors_table.count('class="length-flag') == 1
    assert 'class="length-flag short"' in ancestors_table
    assert "1 residues shorter than the shortest extant sequence." in ancestors_table
    members_table = page.split('class="jobs members-table"', 1)[1].split("</table>", 1)[0]
    flagged = members_table.count('class="length-flag short"')
    assert flagged >= 1 and flagged == members_table.count("1&nbsp;(!)")
    assert "vs extant" not in page and "<svg" not in members_table
    # One histogram of them all, between the ancestors' viewer and their table.
    chart = page.split('class="plot length-chart"', 1)[1].split("</figure>", 1)[0]
    assert "Extant (n = 5): mean 2 · median 2 · sd 0" in chart
    assert chart.count('class="length-mark"') == 2  # Node1 and Node2 share a length
    assert ">Node1, Node2</text>" in chart and ">Node3 · fitch</text>" in chart
    assert "Sequence length" in chart and "outside the extant range" not in chart
    assert "Node3 (MAFFT · IQ-TREE · LG, " in chart and 'data-plot-download="svg"' in chart
    assert page.index("length-chart") < page.index('class="jobs library-table"')


@pytest.mark.django_db
def test_library_page_suggests_the_defaults_to_a_thin_job(
    client: Client, fake_storage: FakeStorage
) -> None:
    thin = finished_job(fake_storage, aligners=["mafft"], modelfinder=True, gap_methods=["fitch"])
    add(client, thin, "Node1")
    page = client.get(f"/jobs/{thin.id}/library/{Library.objects.get().id}/").content.decode()
    assert "This job ran 1 aligner, 1 model and 1 gap method," in re.sub(r"\s+", " ", page)
    assert f'href="/jobs/new/?from={thin.id}">Duplicate it</a>' in page
    rich = finished_job(
        fake_storage, aligners=["mafft", "muscle"], modelfinder=True, fixed_models=True,
        models=["LG"], gap_methods=["fitch", "parsimony"],
    )  # fmt: skip
    add(client, rich, "Node1")
    record = Library.objects.get(job=rich)
    page = client.get(f"/jobs/{rich.id}/library/{record.id}/").content.decode()
    assert "Duplicate it</a>" not in page


@pytest.mark.django_db
def test_site_maps_are_fetched_by_target_key(client: Client, fake_storage: FakeStorage) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node2")
    record = Library.objects.get()
    url = f"/jobs/{job.id}/library/{record.id}/sequence/"
    target = {"combination": KEY, "node": "Node2", "version": "original"}
    fragment = client.get(url, target | {"altall": "0.2"}).content.decode()
    assert 'data-top="35" data-alts="K65 R35"' not in fragment  # K is the ML state
    assert 'data-top="35" data-alts="R35"' in fragment and "sm-tile sm-alt q35" in fragment
    assert "sm-alt" not in client.get(url, target | {"altall": "0.4"}).content.decode()
    # Only the library's own targets, and only the user's library.
    for bad in [{"node": "Node1"}, {"version": "fitch"}, {"combination": "x"}, {"node": ""}]:
        assert client.get(url, target | bad).status_code == 404, bad
    someone = get_user_model().objects.create(username="someone")
    theirs = Library.objects.create(owner=someone, job=job, name="T", code="T", spec=record.spec)
    assert client.get(f"/jobs/{job.id}/library/{theirs.id}/sequence/", target).status_code == 404
    # A reload after a palette change draws the sequences that were open, inline.
    page = client.get(f"/jobs/{job.id}/library/{record.id}/", {"open": "target-1"}).content
    assert b'data-target="target-1" data-drawn="20" open' in page and b'class="sm-site' in page


@pytest.mark.django_db
def test_library_page_controls_take_only_known_values(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node2")
    url = f"/jobs/{job.id}/library/{Library.objects.get().id}/"
    for query in [{"palette": "<script>"}, {"altall": "nan"}, {"altall": "9"}, {"palette": "../x"}]:
        response = client.get(url, query)
        assert response.status_code == 200, query
        assert b"<script>alert" not in response.content
    page = client.get(url, {"palette": "viridis", "altall": "0.9"}).content.decode()
    assert '<option value="viridis" selected>' in page
    assert 'name="altall" value="0.5"' in page  # clamped to the largest meaningful threshold
    # The slider shows the same value and carries no name: the number field submits.
    assert "Uncertainty Threshold" in page
    slider = page.split('class="threshold-slider"', 1)[1].split(">", 1)[0]
    assert 'value="0.5"' in slider and "name=" not in slider


@pytest.mark.django_db
def test_library_page_survives_a_target_the_job_no_longer_has(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node2")
    record = Library.objects.get()
    record.spec, _ = library.add_targets(
        record.spec, [library.target_entry("mafft/iqtree/WAG", "Node2", "original", [], 0)]
    )
    record.save()
    response = client.get(f"/jobs/{job.id}/library/{record.id}/")
    assert response.status_code == 200
    assert "No longer in this job's results." in response.content.decode()
    # Only the target it still has is drawn; the other is listed as missing.
    assert response.content.count(b'"anchors": ["target-1"]') == 1


@pytest.mark.django_db
def test_my_libraries_lists_only_mine_on_live_jobs(
    client: Client, fake_storage: FakeStorage
) -> None:
    assert b"No libraries yet" in client.get("/libraries/").content
    first = finished_job(fake_storage, name="first job")
    second = finished_job(fake_storage, name="second job")
    add(client, first, "Node1")
    add(client, second, "Node1", "Node2")
    someone = get_user_model().objects.create(username="someone")
    Library.objects.create(owner=someone, job=first, name="Not mine", code="X")
    page = client.get("/libraries/").content.decode()
    assert 'href="/libraries/" aria-current="page">Libraries' in page
    assert "first job" in page and "second job" in page and "Not mine" not in page
    assert page.index("second job") < page.index("first job")  # the last changed first
    Job.objects.filter(id=second.id).update(deleted_at=first.created_at)
    page = client.get("/libraries/").content.decode()
    assert "second job" not in page and "first job" in page
    assert b"Libraries</a>" in client.get("/about/").content  # in the top bar everywhere


@pytest.mark.django_db
def test_my_libraries_lists_libraries_with_their_jobs_beneath(
    client: Client, fake_storage: FakeStorage
) -> None:
    thin = finished_job(fake_storage, name="thin job", aligners=["mafft"], gap_methods=["fitch"])
    rich = finished_job(
        fake_storage, name="rich job", aligners=["mafft", "muscle"], tree_builders=["iqtree"]
    )
    add(client, thin, "Node1")
    add(client, rich, "Node2")
    client.post(f"/jobs/{rich.id}/library/new/")  # a second library on the same job
    add(client, rich, "Node3")
    page = client.get("/libraries/").content.decode()
    # One card per library, the last changed first, each with its job as a sub-item.
    assert page.count('class="library-card"') == 3
    assert page.count('class="card-jobs"') == 3
    assert "2 aligners · 1 tree search" in page and "1 aligner · 1 gap method" in page
    cards = page[page.index('class="library-cards"') :]  # past the "Added ..." messages
    assert cards.index("Node3") < cards.index("Node2") < cards.index("Node1")
    # No "this job is thin" notice: the inventory already says how many aligners there are.
    assert "One aligner only" not in page
    # Everything is a draft until a design is saved, and no member count is claimed.
    assert page.count('class="badge draft">Draft<') == 3  # the filter's option is not a badge
    assert "members last exported" not in page

    # Filters: by job, and by a state nothing is in yet.
    only = client.get(f"/libraries/?job={rich.id}").content.decode()
    assert only.count('class="library-card"') == 2 and "thin job" not in only
    assert "No libraries match" in client.get("/libraries/?state=exported").content.decode()
    # An unknown state is ignored rather than refused, as the viewers' parameters are.
    loose = client.get("/libraries/?state=nonsense").content.decode()
    assert loose.count('class="library-card"') == 3


@pytest.mark.django_db
def test_my_libraries_filters_by_the_job_s_tags(client: Client, fake_storage: FakeStorage) -> None:
    """Libraries have no tags of their own: they filter by their job's (owner decision)."""
    job = finished_job(fake_storage, name="tagged job")
    add(client, job, "Node1")
    other = finished_job(fake_storage, name="untagged job")
    add(client, other, "Node1")
    job.add_tags(["lysin"])
    page = client.get("/libraries/?tag=lysin").content.decode()
    assert "tagged job" in page and "untagged job" not in page
    assert client.get("/libraries/?tag=%23lysin").content.decode().count("tagged job") >= 1
    assert "No libraries match" in client.get("/libraries/?tag=missing").content.decode()


@pytest.mark.django_db
def test_my_libraries_lists_deleted_records_beneath(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage, name="a job")
    add(client, job, "Node1")
    record = Library.objects.get()
    client.post(f"/jobs/{job.id}/library/{record.id}/delete/")
    page = client.get("/libraries/").content.decode()
    assert "Deleted Libraries" in page and "No libraries yet" in page
    assert page.count('class="library-card"') == 0


@pytest.mark.django_db
def test_library_page_explains_why_other_jobs_libraries_are_not_destinations(
    client: Client, fake_storage: FakeStorage
) -> None:
    """A target names a combination, node and version of its own job's results, so a
    library on another job cannot hold it. The dialog says so rather than looking short."""
    job = finished_job(fake_storage, name="this job")
    add(client, job, "Node1")
    record = Library.objects.get()
    page = client.get(f"/jobs/{job.id}/library/{record.id}/").content.decode()
    assert "Only this job's libraries are listed" not in page  # no others anywhere yet
    other = finished_job(fake_storage, name="another job")
    add(client, other, "Node1")
    page = client.get(f"/jobs/{job.id}/library/{record.id}/").content.decode()
    assert "Only this job's libraries are listed" in page
    assert "Your other libraries (1)" in page


@pytest.mark.django_db
def test_library_page_designs_the_set_and_saves_it(
    client: Client, fake_storage: FakeStorage
) -> None:
    """The Design section turns targets into members: a live count, a budget, a strategy
    and the coverage readout, with the controls in the address until Save writes them."""
    job = finished_job(fake_storage)
    add(client, job, "Node2", "Node3")
    record = Library.objects.get()
    page = client.get(f"/jobs/{job.id}/library/{record.id}/").content.decode()
    assert "Design the Set" in page
    # Two targets, AltAll on and the default ten samples: 2 x (1 + 1 + 10) candidates.
    assert "24 candidates" in page
    assert f'value="{library.DEFAULT_BUDGET}"' in page
    assert "distinct sequence" in page and 'class="design-coverage' in page

    # The address governs before anything is saved, so a view is shareable.
    fewer = client.get(f"/jobs/{job.id}/library/{record.id}/?samples=0&budget=3").content.decode()
    assert "4 candidates" in fewer  # 2 x (ML + AltAll)
    assert library.design_from(Library.objects.get().spec).axes.samples == library.DEFAULT_SAMPLES
    # "I can order" alone stays with the library the moment it is changed: a reload that
    # sends only the view (the Controls panel's) used to put the saved number back.
    assert library.design_from(Library.objects.get().spec).budget == 3
    again = client.get(f"/jobs/{job.id}/library/{record.id}/?extant=a").content.decode()
    assert 'name="budget" id="budget-number" value="3"' in again
    assert Library.objects.get().spec["targets"] == record.spec["targets"]

    # Save writes the spec; nothing else in it changes.
    before = Library.objects.get().spec
    response = client.post(
        f"/jobs/{job.id}/library/{record.id}/design/",
        {"budget": "8", "strategy": "spanning", "samples": "3", "subs": "1",
         "alt_all": "1", "altall": "0.25",
         "next": f"/jobs/{job.id}/library/{record.id}/"},
    )  # fmt: skip
    assert response.status_code == 302
    saved = library.design_from(Library.objects.get().spec)
    assert saved.budget == 8 and saved.strategy == "spanning"
    assert saved.axes.samples == 3 and saved.axes.max_substitutions == 1
    assert saved.axes.alt_all and saved.axes.threshold == 0.25
    assert Library.objects.get().spec["targets"] == before["targets"]
    assert Library.objects.get().spec["seed"] == before["seed"]


@pytest.mark.django_db
def test_library_design_follows_the_uncertainty_threshold(
    client: Client, fake_storage: FakeStorage
) -> None:
    """The threshold is what makes a site uncertain, so it decides how many sequences
    there are to choose from. The design form carries it, and the page's script keeps that
    copy in step with the slider: without it, Update sent the value the page was loaded
    with and the threshold appeared to snap back (owner, 2026-10-08)."""
    job = finished_job(fake_storage)
    add(client, job, "Node2", "Node3")
    record = Library.objects.get()
    counts = {}
    for threshold in ("0.05", "0.2", "0.5"):
        page = client.get(
            f"/jobs/{job.id}/library/{record.id}/?altall={threshold}&subs=1"
        ).content.decode()
        found = re.search(r"(\d+) candidates?", page)
        assert found is not None
        counts[threshold] = int(found.group(1))
        # The form sends back the threshold it was drawn at, not the stored one.
        assert f'name="altall" value="{threshold}" data-design-threshold' in page
    # A higher threshold leaves fewer sites uncertain, so fewer single-site alternates.
    # Raising it never adds sequences, and the fixture's one borderline site drops out
    # between 0.2 and 0.5, so the count falls there.
    assert counts["0.05"] >= counts["0.2"] > counts["0.5"]


def viewer_data(page: str, frame: str) -> dict:  # type: ignore[type-arg]
    script = rf'<script id="aln-data-{frame}" type="application/json">(.*?)</script>'
    found = re.search(script, page, re.S)
    assert found is not None, f"no viewer {frame}"
    return json.loads(found.group(1))  # type: ignore[no-any-return]


@pytest.mark.django_db
def test_members_are_drawn_as_rows_beneath_their_ancestors(
    client: Client, fake_storage: FakeStorage
) -> None:
    """L5c: the Design section draws the members as alignment-viewer rows, one per
    member, each tied to its row in the table, with its changes from the ML sequence
    marked. Which members are drawn is part of the saved view."""
    job = finished_job(fake_storage)
    add(client, job, "Node2")
    add(client, job, "Node3", gaps="fitch")
    record = Library.objects.get()
    url = f"/jobs/{job.id}/library/{record.id}/"
    page = client.get(url).content.decode()
    data = viewer_data(page, "m0")
    assert data["names"] == ["IME2TE-C1-N2-ML", "IME2TE-C1-N2-ALT", "IME2TE-C2-N3-ML-F"]
    assert all(data["member"]) and not any(data["dim"])
    for anchor in data["anchors"]:
        assert f'<tr id="{anchor}" class="picked">' in page
    # AltAll changes Node2's one uncertain site; a gap version's ML is its own reference.
    assert data["marked"] == [[], [1], []]
    assert data["rows"] == ["MK", "MR", "M-"]
    # The ancestors' viewer has no member rows.
    assert not any(viewer_data(page, "0")["member"])
    assert 'value="selected" form="design-form"' in page

    client.post(
        f"{url}design/",
        {"budget": "5", "strategy": "spanning", "samples": "10", "subs": "0", "alt_all": "1",
         "altall": "0.15", "members": "all", "next": url},
    )  # fmt: skip
    assert library.view_from(Library.objects.get().spec).members == "all"
    assert 'name="members" value="all"' in client.get(url).content.decode()


@pytest.mark.django_db
def test_members_viewer_downloads_as_a_coloured_workbook(
    client: Client, fake_storage: FakeStorage
) -> None:
    """ "Download Excel" beside the members viewer gives its rows as coloured cells, boxed
    where a member leaves the ML sequence, with the members table and what was drawn."""
    from io import BytesIO

    from openpyxl import load_workbook

    job = finished_job(fake_storage)
    add(client, job, "Node2")
    record = Library.objects.get()
    url = f"/jobs/{job.id}/library/{record.id}/"
    page = client.get(f"{url}?palette=viridis").content.decode()
    link = re.search(r'<a href="([^"]+)" data-excel', page)
    assert link and "export=xlsx" in link.group(1) and "palette=viridis" in link.group(1)
    response = client.get(link.group(1).replace("&amp;", "&"))
    assert response["Content-Type"].startswith("application/vnd.openxmlformats")
    assert 'filename="IME2TE.members.xlsx"' in response["Content-Disposition"]
    book = load_workbook(BytesIO(response.content))
    assert book.sheetnames == ["Sequences", "Probabilities", "Members", "About"]
    grid = [list(r) for r in book["Sequences"].iter_rows()]
    assert [c.value for c in grid[1][:3]] == ["Alignment column", 1, 2]
    assert [c.value for c in grid[2][:3]] == ["IME2TE-C1-N2-ML", "M", "K"]
    assert [c.value for c in grid[3][:3]] == ["IME2TE-C1-N2-ALT", "M", "R"]
    assert grid[3][2].border.left.style == "medium"  # AltAll's change, boxed
    assert grid[2][2].border.left.style is None
    assert grid[2][1].fill.start_color.rgb.endswith(
        sites_excel.hex_of(ancestors.rgb(1.0, "viridis"))
    )
    numbers = [list(r) for r in book["Probabilities"].iter_rows(values_only=True)]
    assert numbers[2][1] == 1.0
    members = list(book["Members"].iter_rows(values_only=True))
    assert members[0][0] == "Name" and members[1][0] == "IME2TE-C1-N2-ML"
    about = {r[0]: r[1] for r in book["About"].iter_rows(max_row=6, values_only=True)}
    assert about["Drawn"] == "The selected set"


def test_member_viewers_draw_the_cut_only_when_asked() -> None:
    """The selected set by default; every member on request, the budget's cut dimmed.
    A member is drawn on its own ancestor's entry, never on another target's."""
    target = library_glyphs.Entry(
        "Node2", "target-1", "Node2", "MKAL", [[("M", 1.0)], [("K", 0.6), ("R", 0.4)],
        [("A", 1.0)], [("L", 1.0)]], "mafft", "MAFFT", {"a": "MKAL"},
    )  # fmt: skip
    viewers = library_views.Viewers([], library_glyphs.Layout("mafft"), [], [], {0: target})
    source = library.Source(combination=1, node="Node2")

    def row(n: int, sequence: str, picked: bool) -> dict[str, object]:
        member = library.Member(sequence=sequence, kind="sample", variant=f"S{n}", source=source)
        return {"member": member, "name": f"X-S{n}", "picked": picked, "anchor": f"member-{n}"}

    elsewhere = row(9, "MRAL", True)
    elsewhere["member"] = library.Member(  # a target the page could not draw
        sequence="MRAL", kind="sample", variant="S9", source=library.Source(2, "Node3")
    )
    rows = [row(1, "MKAL", True), row(2, "MRAL", False), elsewhere]
    frames, _ = library_views.member_viewers(rows, viewers, "selected", "coolwarm", 20)
    ((frame, data),) = frames
    assert data["names"] == ["X-S1"] and frame.json_id == "aln-data-m0"
    frames, _ = library_views.member_viewers(rows, viewers, "all", "coolwarm", 20)
    ((frame, data),) = frames
    assert data["names"] == ["X-S1", "X-S2"] and data["dim"] == [False, True]
    assert data["marked"] == [[], [1]]
    assert library_views.member_viewers([], viewers, "all", "coolwarm", 20) == ([], "own")


@pytest.mark.django_db
def test_library_design_refuses_more_candidates_than_it_will_generate(
    client: Client, fake_storage: FakeStorage
) -> None:
    """`count_candidates` runs before anything is generated, so a wide set says so
    instead of hanging."""
    job = finished_job(fake_storage)
    add(client, job, "Node2", "Node3")
    record = Library.objects.get()
    # Samples are clamped per ancestor, so the limit is reached by having two of them.
    page = client.get(f"/jobs/{job.id}/library/{record.id}/?samples=10000").content.decode()
    assert "over the" in page and "this page will generate" in page
    assert "distinct sequence" not in page


@pytest.mark.django_db
def test_library_page_colours_every_sequence_by_the_chosen_scheme(
    client: Client, fake_storage: FakeStorage
) -> None:
    """One Colour choice paints the viewer frames and the site maps alike. The palette
    belongs to the posterior scheme, so it is disabled under the others rather than left
    looking as though it still applied (owner, 2026-10-08)."""
    job = finished_job(fake_storage)
    add(client, job, "Node2")
    record = Library.objects.get()
    page = client.get(f"/jobs/{job.id}/library/{record.id}/").content.decode()
    assert '"scheme": "posterior"' in page  # the viewers' data, and the default
    assert 'name="palette" data-palette>' in page  # not disabled
    assert "Posterior probability" in page

    chem = client.get(f"/jobs/{job.id}/library/{record.id}/?scheme=chemistry").content.decode()
    assert '"scheme": "chemistry"' in chem
    assert 'name="palette" data-palette disabled' in chem
    assert "Acidic D E" in chem  # the scheme's own key replaces the palette's
    # An unknown scheme falls back rather than being refused, as the viewers' parameters do.
    odd = client.get(f"/jobs/{job.id}/library/{record.id}/?scheme=nonsense").content.decode()
    assert '"scheme": "posterior"' in odd


@pytest.mark.django_db
def test_saving_keeps_the_view_as_well_as_the_design(
    client: Client, fake_storage: FakeStorage
) -> None:
    """One Save stores both, so reopening a library brings back how it was left (owner,
    2026-10-08). The address still wins, so a shared link shows what its sender saw."""
    job = finished_job(fake_storage)
    add(client, job, "Node2")
    record = Library.objects.get()
    assert "Not saved yet" in client.get(f"/jobs/{job.id}/library/{record.id}/").content.decode()
    client.post(
        f"/jobs/{job.id}/library/{record.id}/design/",
        {"budget": "48", "strategy": "spanning", "samples": "10", "subs": "0", "alt_all": "1",
         "altall": "0.15", "palette": "viridis", "scheme": "chemistry", "zoom": "fit",
         "aln": "mafft", "next": f"/jobs/{job.id}/library/{record.id}/"},
    )  # fmt: skip
    spec = Library.objects.get().spec
    view = library.view_from(spec)
    assert view.palette == "viridis" and view.scheme == "chemistry"
    assert library.saved_at(spec)  # the moment, for "Saved ... ago"
    assert library.design_from(spec).budget == 48

    page = client.get(f"/jobs/{job.id}/library/{record.id}/").content.decode()
    assert '"scheme": "chemistry"' in page and "Saved" in page
    # The address still wins over the saved view.
    other = client.get(f"/jobs/{job.id}/library/{record.id}/?scheme=none").content.decode()
    assert '"scheme": "none"' in other


@pytest.mark.django_db
def test_duplicating_a_library_copies_its_ancestors_and_design(
    client: Client, fake_storage: FakeStorage
) -> None:
    """Several configurations are several libraries, not a second saved thing inside one
    (owner decision 2026-10-08). The seed is copied, so only what the copy is changed to
    differs."""
    job = finished_job(fake_storage)
    add(client, job, "Node2", "Node3")
    record = Library.objects.get()
    client.post(
        f"/jobs/{job.id}/library/{record.id}/design/",
        {"budget": "12", "strategy": "per_node", "samples": "4", "subs": "0", "alt_all": "1",
         "altall": "0.2", "next": f"/jobs/{job.id}/library/{record.id}/"},
    )  # fmt: skip
    response = client.post(f"/jobs/{job.id}/library/{record.id}/duplicate/")
    assert response.status_code == 302
    copy = Library.objects.exclude(id=record.id).get()
    assert copy.name == f"{record.name} (copy)" and copy.job_id == job.id
    assert copy.spec["seed"] == Library.objects.get(id=record.id).spec["seed"]
    assert [library.target_key(t) for t in copy.targets] == [
        library.target_key(t) for t in Library.objects.get(id=record.id).targets
    ]
    assert library.design_from(copy.spec) == library.design_from(
        Library.objects.get(id=record.id).spec
    )
    # It is a library of its own: changing one leaves the other alone.
    client.post(
        f"/jobs/{job.id}/library/{copy.id}/design/",
        {"budget": "99", "strategy": "supported", "samples": "0", "subs": "0",
         "next": f"/jobs/{job.id}/library/{copy.id}/"},
    )  # fmt: skip
    assert library.design_from(Library.objects.get(id=record.id).spec).budget == 12
    assert library.design_from(Library.objects.get(id=copy.id).spec).budget == 99


@pytest.mark.django_db
def test_exporting_freezes_what_the_page_shows_and_leaves_a_receipt(
    client: Client, fake_storage: FakeStorage
) -> None:
    """An export writes its four documents and a copy of the spec, numbered, and keeps
    them: a set that was ordered can be read back (step L5b)."""
    job = finished_job(fake_storage)
    add(client, job, "Node2", "Node3")
    record = Library.objects.get()
    page = client.get(f"/jobs/{job.id}/library/{record.id}/").content.decode()
    assert "Export (" in page and "Exports" not in page  # the button, but no section yet

    taken = client.post(
        f"/jobs/{job.id}/library/{record.id}/export/",
        {"altall": "0.2", "budget": "3", "strategy": "supported", "samples": "2",
         "subs": "0", "alt_all": "1", "next": f"/jobs/{job.id}/library/{record.id}/"},
    )  # fmt: skip
    assert taken.status_code == 302
    receipt = LibraryExport.objects.get()
    assert receipt.number == 1 and 1 <= receipt.member_count <= 3  # at most the budget
    # The receipt freezes what was exported, not whatever the library last saved.
    assert library.design_from(receipt.spec).budget == 3
    assert library.design_from(Library.objects.get().spec).budget == library.DEFAULT_BUDGET
    written = set(fake_storage.objects) & {
        f"{receipt.storage_prefix}{name}" for name in library_export.FILES
    }
    assert len(written) == len(library_export.FILES)

    # Every file downloads, under a name that says which library and which export.
    for name in library_export.FILES:
        got = client.get(f"/jobs/{job.id}/library/{record.id}/exports/1/{name}")
        assert got.status_code == 200 and got.content
        assert f"-1-{name}" in got["Content-Disposition"]
    assert client.get(f"/jobs/{job.id}/library/{record.id}/exports/1/passwd").status_code == 404
    gone = client.get(f"/jobs/{job.id}/library/{record.id}/exports/9/{library_export.FASTA}")
    assert gone.status_code == 404
    # The FASTA holds what was selected; the manifest lists every member, flagged.
    ordered = client.get(f"/jobs/{job.id}/library/{record.id}/exports/1/{library_export.FASTA}")
    assert ordered.content.decode().count(">") == receipt.member_count
    table = client.get(
        f"/jobs/{job.id}/library/{record.id}/exports/1/{library_export.MANIFEST}"
    ).content.decode()
    assert len(table.splitlines()) - 1 >= receipt.member_count
    flags = [line.split("\t")[3] for line in table.splitlines()[1:]]
    assert flags.count("yes") == receipt.member_count

    # A second export is numbered, never overwriting the first.
    client.post(
        f"/jobs/{job.id}/library/{record.id}/export/",
        {"altall": "0.2", "budget": "2", "strategy": "supported", "samples": "0",
         "subs": "0", "alt_all": "1"},
    )  # fmt: skip
    assert [e.number for e in LibraryExport.objects.order_by("number")] == [1, 2]
    assert "Exports" in client.get(f"/jobs/{job.id}/library/{record.id}/").content.decode()


@pytest.mark.django_db
def test_restoring_an_export_puts_its_design_back(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node2")
    record = Library.objects.get()
    client.post(
        f"/jobs/{job.id}/library/{record.id}/export/",
        {"altall": "0.2", "budget": "7", "strategy": "per_node", "samples": "1",
         "subs": "0", "alt_all": "1"},
    )  # fmt: skip
    client.post(
        f"/jobs/{job.id}/library/{record.id}/design/",
        {"budget": "99", "strategy": "supported", "samples": "5", "subs": "0", "altall": "0.3"},
    )
    assert library.design_from(Library.objects.get().spec).budget == 99
    response = client.post(f"/jobs/{job.id}/library/{record.id}/exports/1/restore/")
    assert response.status_code == 302
    back = library.design_from(Library.objects.get().spec)
    assert back.budget == 7 and back.strategy == "per_node" and back.axes.samples == 1
    # The export itself is untouched, and the library keeps its own seed.
    assert library.design_from(LibraryExport.objects.get().spec).budget == 7
    assert Library.objects.get().spec["seed"] == record.spec["seed"]


@pytest.mark.django_db
def test_deleting_a_job_purges_its_libraries_exports(
    client: Client, fake_storage: FakeStorage
) -> None:
    """A library's exports live under their own prefix, so the job's purge has to name
    them: the ancestors are going, and an export that outlived them would promise a
    provenance nothing could check. The records stay, as deleted jobs' own do."""
    job = finished_job(fake_storage)
    add(client, job, "Node2")
    record = Library.objects.get()
    client.post(
        f"/jobs/{job.id}/library/{record.id}/export/",
        {"altall": "0.2", "budget": "2", "strategy": "supported", "samples": "0",
         "subs": "0", "alt_all": "1"},
    )  # fmt: skip
    assert any(k.startswith(f"libraries/{record.id}/") for k in fake_storage.objects)
    jobs.delete([job])
    assert jobs.purge(Job.objects.get(id=job.id))
    assert not any(k.startswith(f"libraries/{record.id}/") for k in fake_storage.objects)
    assert not any(k.startswith(job.storage_prefix) for k in fake_storage.objects)
    assert LibraryExport.objects.count() == 1 and Library.objects.count() == 1


def two_aligner_job(fake_storage: FakeStorage) -> Job:
    """A finished job with MAFFT and MUSCLE alignments of the same five sequences."""
    job = finished_job(fake_storage, name="two aligners")
    prefix = f"jobs/{job.id}/results/"
    summary = json.loads(fake_storage.objects[prefix + "summary.json"])
    other = "muscle/iqtree/LG"
    summary["combinations"].append(
        {
            "key": other,
            "aligner": "muscle",
            "aligner_label": "MUSCLE",
            "tree": "iqtree",
            "model": "LG",
            "label": "MUSCLE · IQ-TREE · LG",
        }  # fmt: skip
    )
    summary["combinations"][0]["aligner_label"] = "MAFFT"
    files = {
        "summary.json": json.dumps(summary),
        f"{other}/tree.nwk": TREE,
        f"{other}/alignment.fasta": "".join(f">{n}\nM-K\n" for n in "abcde"),
        f"{other}/ancestors.fasta": ">Node1\nM-K\n>Node2\nM-K\n>Node3\nM-K\n",
    }
    for file, text in files.items():
        fake_storage.objects[prefix + file] = text.encode()
    return job


@pytest.mark.django_db
def test_library_page_layouts(client: Client, fake_storage: FakeStorage) -> None:
    job = two_aligner_job(fake_storage)
    add(client, job, "Node2")
    add(client, job, "Node2", combination="muscle/iqtree/LG")
    url = f"/jobs/{job.id}/library/{Library.objects.get().id}/"
    page = client.get(url).content.decode()
    # Two alignments: by default each on its own frame, never stacked as if aligned.
    assert "MAFFT alignment · 2 columns" in page and "MUSCLE alignment · 3 columns" in page
    assert '<option value="own" selected>' in page
    # The "Map ancestors onto" layout menu lives in the Controls panel (the "aln" select,
    # where the agree/disagree tracks appear); separate by default.
    assert 'name="aln"' in page and "Map ancestors onto" in page
    for option in ('value="mafft"', 'value="muscle"', 'value="sequence"'):
        assert option in page
    for name in "abcde":  # every extant sequence is offered, for showing or laying out on
        assert f'name="extant" value="{name}"' in page

    page = client.get(url, {"aln": "mafft"}).content.decode()
    assert "On MAFFT columns · 2" in page and "Aligners agree" in page
    assert "MUSCLE Node2" in page and "MAFFT Node2" in page
    # One control maps both viewers: the Compare members menu is gone, and with the
    # ancestors on MAFFT the members viewer (group "m") is on MAFFT too.
    assert 'name="maln"' not in page and "Compare members" not in page
    assert page.count("On MAFFT columns") >= 2  # ancestors' viewer and members' viewer
    # Extant sequences join the ancestors as rows.
    page = client.get(url, {"aln": "mafft", "extant": ["c", "a"]}).content.decode()
    assert '"names": ["c", "a", "MAFFT Node2"' in page
    page = client.get(url, {"aln": "sequence", "extant": ["c"]}).content.decode()
    assert "On c · 2 residues" in page and "Inserted residues" in page
    for bad in ["../x", "<script>", "nonsense", ""]:
        response = client.get(url, {"aln": bad, "extant": bad})
        assert (
            response.status_code == 200
            and '<option value="own" selected>' in response.content.decode()
        )


@pytest.mark.django_db
def test_library_page_zoom(client: Client, fake_storage: FakeStorage) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node2")
    url = f"/jobs/{job.id}/library/{Library.objects.get().id}/"
    assert '"cell": "fit"' in client.get(url).content.decode()
    assert '"cell": 3.5' in client.get(url, {"zoom": "3.5"}).content.decode()
    for bad in ["<script>", "abc", "-5"]:
        page = client.get(url, {"zoom": bad}).content.decode()
        assert '"cell": "fit"' in page or '"cell": 0.05' in page


@pytest.mark.django_db
def test_add_to_a_chosen_library_or_a_new_one(client: Client, fake_storage: FakeStorage) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node1")
    first = Library.objects.get()
    add(client, job, "Node2", library="new")
    second = Library.objects.exclude(id=first.id).get()
    assert [t["node"] for t in second.targets] == ["Node2"]
    add(client, job, "Node3", library=str(first.id))  # not the current one
    first.refresh_from_db()
    assert [t["node"] for t in first.targets] == ["Node1", "Node3"]


def keys(*nodes: str, version: str = "original") -> list[str]:
    return [json.dumps([KEY, node, version]) for node in nodes]


@pytest.mark.django_db
def test_remove_copy_and_move_the_ticked_ancestors(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node1", "Node2", "Node3")
    source = Library.objects.get()
    url = f"/jobs/{job.id}/library/{source.id}/targets/"
    page = f"/jobs/{job.id}/library/{source.id}/"

    # Move two to a new library; the page returns to where it was.
    response = client.post(
        url, {"action": "move", "to": "new", "target": keys("Node1", "Node2"), "next": page}
    )
    assert response["Location"] == page
    dest = Library.objects.exclude(id=source.id).get()
    source.refresh_from_db()
    assert [t["node"] for t in source.targets] == ["Node3"]
    assert [t["node"] for t in dest.targets] == ["Node1", "Node2"]
    assert dest.targets[0]["partition"]  # copied whole, partition included

    # Copy: one new, one already there, and one the page lists but the library lost.
    client.post(url, {"action": "copy", "to": str(dest.id), "target": keys("Node3", "Node1")})
    dest.refresh_from_db()
    source.refresh_from_db()
    assert [t["node"] for t in dest.targets] == ["Node1", "Node2", "Node3"]
    assert [t["node"] for t in source.targets] == ["Node3"]  # a copy leaves it in place
    response = client.post(
        url, {"action": "copy", "to": str(dest.id), "target": keys("Node3")}, follow=True
    )
    assert "Copied 1 ancestor to Library 2 (1 already there)." in response.content.decode()

    client.post(url, {"action": "remove", "target": keys("Node3")})
    source.refresh_from_db()
    assert source.targets == []


@pytest.mark.django_db
def test_bulk_actions_refuse_what_is_not_the_users(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node1")
    source = Library.objects.get()
    url = f"/jobs/{job.id}/library/{source.id}/targets/"
    elsewhere = finished_job(fake_storage, name="other")
    add(client, elsewhere, "Node1")
    other_job_library = Library.objects.get(job=elsewhere)
    someone = get_user_model().objects.create(username="someone")
    theirs = Library.objects.create(owner=someone, job=job, name="T", code="T")
    target = keys("Node1")
    for bad in [
        {"action": "drop", "target": target},
        {"action": "remove", "target": ["{not json"]},
        {"action": "remove", "target": ['["only", "two"]']},
        {"action": "remove", "target": ["[1, 2, 3]"]},
        {"action": "copy", "target": target},  # no destination
        {"action": "copy", "target": target, "to": "not-a-uuid"},
        {"action": "copy", "target": target, "to": str(other_job_library.id)},  # not this job
        {"action": "move", "target": target, "to": str(theirs.id)},
    ]:
        assert client.post(url, bad).status_code == 404, bad
    assert client.post(f"/jobs/{job.id}/library/{theirs.id}/targets/",
                       {"action": "remove", "target": target}).status_code == 404  # fmt: skip
    source.refresh_from_db()
    other_job_library.refresh_from_db()
    assert len(source.targets) == 1 and len(other_job_library.targets) == 1
    # Copying a library onto itself is refused politely, and changes nothing.
    response = client.post(
        url, {"action": "copy", "to": str(source.id), "target": target}, follow=True
    )
    assert "Choose another library" in response.content.decode()
    source.refresh_from_db()
    assert len(source.targets) == 1


@pytest.mark.django_db
def test_copy_to_several_libraries_at_once(client: Client, fake_storage: FakeStorage) -> None:
    job = finished_job(fake_storage)
    add(client, job, "Node1", "Node2")
    source = Library.objects.get()
    url = f"/jobs/{job.id}/library/{source.id}/targets/"
    # Copy into two new libraries in one action.
    client.post(url, {"action": "copy", "to": ["new", "new"], "target": keys("Node1")})
    # "new" twice is one destination: the action names the libraries it wrote to.
    made = Library.objects.exclude(id=source.id)
    assert made.count() == 1
    client.post(url, {"action": "copy", "to": ["new"], "target": keys("Node2")})
    first, second = Library.objects.exclude(id=source.id).order_by("created_at")
    client.post(url, {"action": "copy", "to": [str(first.id), str(second.id)],
                      "target": keys("Node1", "Node2")})  # fmt: skip
    first.refresh_from_db()
    second.refresh_from_db()
    assert [t["node"] for t in first.targets] == ["Node1", "Node2"]
    assert [t["node"] for t in second.targets] == ["Node2", "Node1"]
    source.refresh_from_db()
    assert len(source.targets) == 2  # a copy leaves them in place
    # Moving to several libraries is not one action.
    assert client.post(url, {"action": "move", "to": [str(first.id), str(second.id)],
                             "target": keys("Node1")}).status_code == 404  # fmt: skip
    assert client.post(url, {"action": "copy", "target": keys("Node1")}).status_code == 404
