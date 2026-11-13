"""The mutational path page (phylobot.path_views): a path picked on the tree viewer opens
on a page of its own, compared across the job's combinations, and can be pinned."""

import json
import re
from io import BytesIO
from typing import Any
from urllib.parse import parse_qs, urlsplit

import pytest
from django.contrib.auth.models import User
from django.test import Client
from openpyxl import load_workbook

from phylobot.models import Job, SavedPath
from tests.conftest import FakeStorage
from tests.test_views import viewer_job

pytestmark = pytest.mark.django_db

MAFFT, WAG, MUSCLE = "mafft/iqtree/modelfinder", "mafft/iqtree/WAG", "muscle/iqtree/modelfinder"
# viewer_job: MAFFT's tree ((a,b)Node2,c<i>)Node1, MUSCLE's the same numbered the other way
# round; every combination's ancestors are Node1 MKK and Node2 MRK, so along the same
# path MAFFT reads MKK -> MRK -> a (MK-) and MUSCLE reads MRK -> MKK -> a.
PATH = {"path_from": "Node1", "path_to": "a", "gaps": "original"}


def page(client: Client, job: Job, htmx: bool = False, **query: Any) -> Any:
    headers = {"HX-Request": "1"} if htmx else {}
    return client.get(f"/jobs/{job.id}/path/", PATH | query, headers=headers)


def frames(html: str) -> list[dict[str, Any]]:
    """The alignment viewers' data, in page order."""
    found = re.findall(r'<script id="aln-data-p\d+" type="application/json">(.*?)</script>', html)
    return [json.loads(text) for text in found]


def test_tree_viewer_hands_the_path_to_its_page(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    url = f"/jobs/{job.id}/tree/"
    html = client.get(url, {"path_from": "Node1"}).content.decode()
    # The ends are picked on the tree (path_pick.js): each branch names its ends, each
    # tip its name, and the form says where a picked path opens.
    assert '<input type="hidden" name="path_from" value="Node1">' in html
    assert '<input type="hidden" name="path_to" value="">' in html
    assert html.count("data-path-pick") == 1 and f'data-path-url="/jobs/{job.id}/path/"' in html
    assert 'data-child="Node2" data-parent="Node1"' in html
    assert 'data-child="c&lt;i&gt;" data-parent="Node1"' in html  # escaped
    assert 'class="tree-leaf" data-leaf="a"' in html
    # A traced path is highlighted on the tree and linked to its page; no table here.
    query = {"path_from": "Node1", "path_to": "a", "trace": "1", "lengths": "true"}
    response = client.get(url, query, headers={"HX-Request": "1"})
    html = response.content.decode()
    assert html.count('stroke="#e08214" stroke-width="4" stroke-linecap="round"') == 2  # elbows
    assert "path-table" not in html and "K2R" not in html
    link = re.search(r'<a class="button with-icon" href="([^"]+)">', html)
    assert link and "Open Node1 → a" in html
    opened = parse_qs(urlsplit(link.group(1).replace("&amp;", "&")).query)
    assert opened["path_from"] == ["Node1"] and opened["path_to"] == ["a"]
    assert opened["combination"] == [MAFFT] and "lengths=true" in opened["back"][0]
    assert "path_from=Node1&path_to=a" in response["HX-Push-Url"]
    # an end that is not below the start is dropped; Clear drops the path
    html = client.get(url, query | {"path_from": "Node2", "path_to": "c<i>"}).content.decode()
    assert "Open Node2" not in html
    response = client.get(url, query | {"path_clear": "1"})
    assert response.status_code == 302 and "path_" not in response["Location"]
    assert client.get(url, query | {"export": "xlsx"}).status_code == 200  # moved to the page


def test_path_page_compares_the_path_across_alignments(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = viewer_job(fake_storage)
    response = page(client, job, extant=["b", "nope"], back="lengths=true&next=//evil.example")
    html = response.content.decode()
    assert response.status_code == 200 and "<h1>Mutational Path</h1>" in html
    assert "<strong>Node1 → Node2 → a</strong>" in html
    # the way back: the tree as it was drawn, the path highlighted, nothing else carried
    back = re.search(r'id="path-back"><a href="([^"]+)">← Return to Tree Viewer', html)
    assert back
    tree = parse_qs(urlsplit(back.group(1).replace("&amp;", "&")).query)
    assert tree["lengths"] == ["true"] and tree["trace"] == ["1"] and "next" not in tree
    assert tree["path_from"] == ["Node1"] and tree["path_to"] == ["a"]
    # Checked against MUSCLE's best (MAFFT's own best is the one traced): there Node1 ->
    # Node2 is R -> K, so K2R is this combination's alone; losing the K into a is shared.
    assert "Checked against 1 other combination" in html
    assert re.search(r'class="change substitution lvl-none" title="[^"]*different: MUSCLE[^"]*'
                     r'\(R→K\)[^"]*">K2R<sub>0/1</sub>', html)  # fmt: skip
    assert re.search(r'class="change loss lvl-all"[^>]*>K3-<sub>1/1</sub>', html)
    assert 'id="branch-0"' in html and 'href="#branch-1"' in html
    # the line: three stops, a band per level present, the tip marked extant
    assert html.count("<circle") == 3 and 'fill="#eb6834"' in html and 'fill="#1f66c1"' in html
    assert "Real branch lengths" in html and 'font-style="italic">a</text>' in html
    # the support table, a branch per row
    assert "Support Across Combinations" in html and ">0/1</td>" in html and ">1/2</td>" in html
    # the traced alignment's viewer: the path's rows, then the chosen extant sequence
    (data,) = frames(html)
    assert data["names"] == ["Node1", "Node2", "a", "b"] and data["scheme"] == "changes"
    assert data["rows"] == ["MKK", "MRK", "MK-", "M-K"]
    assert data["marked"] == [[], [1], [1, 2], []] and data["onPath"] == [True, True, True, False]
    assert data["notes"][1]["1"].startswith("K2R on Node1 → Node2\nthe same in 0 of 1")
    assert data["tracks"] == []
    assert "Select Path</h2>" in html and 'id="help-path-line"' in html
    assert "The path along branches from Node1 to a," in html
    # the old table is still here, for this combination, with the extant row
    assert "path-table" in html and "For comparison" in html and "<em>b</em>" in html
    assert "nope" not in html
    response = page(client, job, extant=["b"], htmx=True)
    assert "path_from=Node1&path_to=a" in response["HX-Push-Url"]
    assert "extant=b" in response["HX-Push-Url"] and "hx-swap-oob" in response.content.decode()
    assert "<h1>" not in response.content.decode()


def test_path_page_shows_other_alignments_beside_it(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = viewer_job(fake_storage)
    # "+" adds a combination once: the address then names it
    response = page(client, job, add=MUSCLE)
    assert response.status_code == 302 and "add=" not in response["Location"]
    assert "also=muscle%2Fiqtree%2Fmodelfinder" in response["Location"]
    response = page(client, job, add_best="1", htmx=True)
    assert "also=muscle%2Fiqtree%2Fmodelfinder" in response["HX-Push-Url"]
    html = page(client, job, also=[MUSCLE, "nope/x/y", MAFFT], mutations="support").content.decode()
    main, other = frames(html)
    assert other["names"] == ["Node2", "Node1", "a"]  # the same ancestors, as MUSCLE names them
    assert other["rows"] == ["MRK", "MKK", "MK-"] and other["sync"] == [0, 1, 2]
    assert other["match"] == [0, 1, 2] and main["match"] is None and main["hoverEvents"]
    # coloured by level: MUSCLE's own R2K is not on the traced path, the lost K is shared
    assert main["levelColours"][3] == "#eb6834" and main["levels"][1] == {"1": 3}
    assert "in every alignment × model combo" in html and "in most combos" in html
    assert other["levels"][1] == {"1": 5} and other["levels"][2] == {"2": 0}
    assert "1 of the traced path's 3 comparable changes is made here too" in html
    assert 'name="also" value="muscle/iqtree/modelfinder"' in html and "Trace here" in html
    # a control named "remove" would hide the form's remove(), and HTMX could not swap it
    assert 'name="drop" value="muscle/iqtree/modelfinder"' in html and 'name="remove"' not in html
    # side by side every ancestor sits beside itself (here the trees agree: no empty rows)
    html = page(client, job, also=[MUSCLE], layout="side", mutations="red").content.decode()
    main, other = frames(html)
    assert 'class="path-frames side"' in html and html.count('class="aln-overview"') == 2
    assert 'name="layout" value="side" checked' in html
    assert main["scheme"] == "posterior" and main["markColour"] == "#d7191c"
    assert len(main["rows"]) == len(other["rows"]) == 3
    response = page(client, job, also=[MUSCLE, WAG], drop=MUSCLE)
    assert "also=mafft%2Fiqtree%2FWAG" in response["Location"]
    assert "muscle" not in response["Location"].split("also=")[1]
    # "Trace here": the path carried onto MUSCLE's tree by clade, MAFFT kept beside it
    response = client.get(
        f"/jobs/{job.id}/path/",
        PATH | {"combination": MUSCLE, "aligner": "muscle", "of": MAFFT, "also": [MAFFT]},
    )
    moved = parse_qs(urlsplit(response["Location"]).query)
    assert moved["combination"] == [MUSCLE] and moved["path_from"] == ["Node2"]
    assert moved["path_to"] == ["a"] and moved["also"] == [MAFFT] and "of" not in moved


def test_path_page_without_a_path(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    html = client.get(f"/jobs/{job.id}/path/").content.decode()
    assert "No path is traced" in html and "path-frames" not in html
    # a descendant that is not below the start is no path
    html = page(client, job, path_from="Node2", path_to="c<i>").content.decode()
    assert "No path is traced" in html and "c<i>" not in html
    assert page(client, job, export="xlsx").status_code == 200
    assert client.get(f"/jobs/{job.id}/path/", {"export": "xlsx"}).status_code == 200
    # someone else's job is not found
    stranger = Job.objects.create(owner=User.objects.create(username="other"), params={})
    assert client.get(f"/jobs/{stranger.id}/path/").status_code == 404
    unfinished = Job.objects.create(owner=job.owner, params={})
    assert "no ancestral sequences" in client.get(f"/jobs/{unfinished.id}/path/").content.decode()


def test_path_page_with_codons(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)  # ancestors Node1 MKK, Node2 MRK; extant a MK-
    html = page(client, job).content.decode()
    assert 'name="silent"' not in html  # protein: no codon options
    # A codon reconstruction: Node1 -> Node2 changes K2R (AAA -> AGA) and, silently,
    # K3K (AAA -> AAG).
    results = f"jobs/{job.id}/results/{MAFFT}/"
    fake_storage.objects[results + "alignment.cds.fasta"] = (
        b">a\nATGAGA---\n>b\nATG---AAG\n>c<i>\nATGAAAAAG\n"
    )
    fake_storage.objects[results + "ancestors.cds.fasta"] = (
        b">Node1\nATGAAAAAA\n>Node2\nATGAGAAAG\n"
    )
    html = page(client, job).content.decode()
    assert '<span class="change substitution lvl-none" title="column 2 · AAA → AGA' in html
    assert ">K3K</span>" in html and "(silent)" in html  # a silent change is not compared
    assert "1 amino-acid change (non-synonymous), 1 silent (synonymous)" in html
    assert 'class="silent"' in html and 'title="AAG"' in html
    html = page(client, job, silent="hide").content.decode()
    assert ">K3K</span>" not in html and "1 silent (synonymous), hidden" in html
    response = page(client, job, cells="codons", htmx=True)
    html = response.content.decode()
    assert 'class="path-table codons"' in html and ">AGA" in html
    assert "silent=show&cells=codons" in response["HX-Push-Url"]
    link = re.search(r'<a href="([^"]+)">Download Excel</a>', html)
    assert link
    book = load_workbook(BytesIO(client.get(link.group(1).replace("&amp;", "&")).content))
    changes = [[c.value for c in row] for row in book["Changes"].iter_rows()]
    assert "Codons" in changes[0]
    assert changes[2][2:4] == ["K3K", "synonymous"] and "AAA → AAG" in changes[2]
    # codeml's free ratios ran: dN/dS beside each branch (undefined without synonymous
    # changes), whichever way codeml's tree lists the two ends.
    assert "branch-rate" not in html
    fake_storage.objects[results + "codeml/FR/branches.tsv"] = (
        b"parent\tchild\tt\tN\tS\tdN_dS\tdN\tdS\tN_dN\tS_dS\n"
        b"Node2\tNode1\t0.1\t300\t100\t0.42\t0.01\t0.024\t3\t2.4\n"
        b"Node2\ta\t0.1\t300\t100\t999\t0.02\t0\t6\t0\n"
    )
    html = page(client, job).content.decode()
    assert ">dN/dS 0.42</span>" in html and "free-ratio estimates" in html
    assert "dN/dS undefined: fewer than one synonymous change expected" in html


def test_path_excel_and_probabilities(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    # probabilities for every site of both ancestors (the fixture's pp file has one)
    fake_storage.objects[f"jobs/{job.id}/results/{MAFFT}/ancestors_pp.tsv"] = (
        "node\tsite\tstate\tpp\n"
        + "".join(f"{n}\t{s}\tM\t0.{s}{s}\n" for n in ("Node1", "Node2") for s in (1, 2, 3))
    ).encode()
    html = page(client, job, extant=["b"]).content.decode()
    assert "K<small>0.22</small>" in html  # Node1, column 2: K at 0.22
    (data,) = frames(html)
    assert data["pp"][:2] == [[11, 22, 33], [11, 22, 33]] and data["pp"][2] == []  # extant
    link = re.search(r'<a href="([^"]+)">Download Excel</a>', html)
    assert link and "export=xlsx" in link.group(1)
    response = client.get(link.group(1).replace("&amp;", "&"))
    assert 'filename="path.Node1-a.xlsx"' in response["Content-Disposition"]
    book = load_workbook(BytesIO(response.content))
    assert book.sheetnames == ["Path", "Probabilities", "Changes", "About"]
    grid = [[c.value for c in row] for row in book["Path"].iter_rows()]
    assert grid[0][0] == "Alignment column" and grid[2][0] == "Node1"
    assert grid[4][0] == "a" and grid[5][0] == "For comparison" and grid[6][0] == "b"
    probs = [[c.value for c in row] for row in book["Probabilities"].iter_rows()]
    assert probs[2][1] == 0.22 and probs[6][1] is None  # extant: no probability
    changes = [[c.value for c in row] for row in book["Changes"].iter_rows()]
    assert changes[1][:4] == ["Node1 → Node2", 2, "K2R", "substitution"]
    # each change with what the compared combinations say
    assert changes[0][-3:] == ["Same in", "Compared with", "Support"]
    assert changes[1][-3:-1] == [0, 1] and "different: MUSCLE" in str(changes[1][-1])
    assert book["Path"]["B4"].border.left.style == "medium"  # Node2 changed at column 2


def test_saved_paths(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    html = page(client, job, also=[MUSCLE], layout="side", back="lengths=true").content.decode()
    view = re.search(r'name="view" value="([^"]+)" form="path-save"', html)
    assert view and "Save path" in html and "None yet" in html
    posted = view.group(1).replace("&amp;", "&")
    assert "back=" not in posted  # where one came from is not part of the path
    url = f"/jobs/{job.id}/path/save/"
    response = client.post(url, {"view": posted, "name": "  "})
    assert response.status_code == 302 and not SavedPath.objects.exists()
    response = client.post(url, {"view": posted, "name": "  To   <b>a</b> "})
    saved = SavedPath.objects.get()
    assert (saved.name, saved.combination, saved.start, saved.end) == (
        "To <b>a</b>", MAFFT, "Node1", "a",
    )  # fmt: skip
    assert saved.view["also"] == [MUSCLE] and saved.view["layout"] == "side"
    html = client.get(response["Location"]).content.decode()
    assert "To &lt;b&gt;a&lt;/b&gt;" in html and "<b>a</b>" not in html and "Save again" in html
    # listed on the tree viewer too, and it opens the view it was saved with
    tree = client.get(f"/jobs/{job.id}/tree/").content.decode()
    link = re.search(r'<li><a href="([^"]+)">To &lt;b&gt;a', tree)
    assert link and 'aria-expanded="true" aria-controls="mutational-path-body"' in tree
    assert "Mutational Paths (1 saved)" in tree
    opened = client.get(link.group(1).replace("&amp;", "&")).content.decode()
    assert 'class="path-frames side"' in opened and len(frames(opened)) == 2
    # saving the same path again renames it; another path is another entry
    client.post(url, {"view": posted, "name": "Renamed"})
    assert SavedPath.objects.get().name == "Renamed"
    # a path that is not one of the job's is refused, as is anyone else's job
    assert (
        client.post(url, {"view": "path_from=Node2&path_to=nope", "name": "x"}).status_code == 404
    )
    stranger = Job.objects.create(owner=User.objects.create(username="other"), params={})
    assert (
        client.post(f"/jobs/{stranger.id}/path/save/", {"view": posted, "name": "x"}).status_code
        == 404
    )
    remove = f"/jobs/{job.id}/path/saved/{saved.id}/delete/"
    assert client.get(remove).status_code == 405
    response = client.post(remove, {"next": "https://evil.example/"})
    assert response["Location"] == f"/jobs/{job.id}/tree/" and not SavedPath.objects.exists()
