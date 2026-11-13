"""The structure page (phylobot.structure_views): attaching structures to a job, and its
ancestors painted on them, on one panel or several. The network is always mocked."""

import io
import json
import re
import zipfile
from typing import Any
from urllib.parse import parse_qs, urlsplit

import pytest
from django.core.cache import cache
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import Client

from phylobot import jobs
from phylobot.homologs import remote as homolog_remote
from phylobot.models import Job, Structure
from phylobot.modes import local_user
from phylobot.structures import colours
from tests.conftest import FakeStorage
from tests.test_structures import cif_text, pdb_text

pytestmark = pytest.mark.django_db

KEY = "mafft/iqtree/LG"
ALIGNMENT = {
    "a_P12931": "MKTAYIAKQRQISFVKSHFSRQ",
    "b": "MKTAYIAKQ--ISFVKSHFSRQ",
    "c": "MRTAYLAKQRQISFVRSHFSRQ",
}
ANCESTORS = {"Node1": "MKTAYIAKQRQISFVKSHFSRQ", "Node2": "MKTAYLAKQRQISFIKSHFSRQ"}
# Chain A is a_P12931's residues 3-18, numbered from 10; chain B is like nothing in the job.
STRUCTURE = pdb_text({"A": "TAYIAKQRQISFVKSH", "B": "GGGGGGGGGG"}, start=10)
MODEL = cif_text(ANCESTORS["Node2"], title="Model of Node2")  # as AlphaFold 3 writes one


@pytest.fixture(autouse=True)
def _fresh_cache() -> None:
    cache.clear()


def fasta(rows: dict[str, str]) -> str:
    return "".join(f">{name}\n{row}\n" for name, row in rows.items())


def structure_job(fake_storage: FakeStorage) -> Job:
    job = Job.objects.create(owner=local_user(), params={})
    assert jobs.mark(job.id, Job.Status.RUNNING) and jobs.mark(job.id, Job.Status.SUCCEEDED)
    prefix = f"jobs/{job.id}/results/"
    combo = {"key": KEY, "aligner": "mafft", "tree": "iqtree", "model": "LG", "rank_bic": 1}
    states = sorted(set("".join(ANCESTORS.values())))
    posteriors = "node\tsite\t" + "\t".join(states) + "\n"
    for node, row in ANCESTORS.items():
        for site, state in enumerate(row, start=1):
            values = ("0.8" if s == state else "0.01" for s in states)
            posteriors += f"{node}\t{site}\t" + "\t".join(values) + "\n"
    files = {
        "alignment.fasta": fasta(ALIGNMENT),
        "tree.nwk": "((a_P12931,b)Node2,c)Node1;",
        "ancestors.fasta": fasta(ANCESTORS),
        "ancestors_posteriors.tsv": posteriors,
    }
    fake_storage.objects[prefix + "summary.json"] = json.dumps(
        {"primary": KEY, "combinations": [combo]}
    ).encode()
    for name, text in files.items():
        fake_storage.objects[f"{prefix}{KEY}/{name}"] = text.encode()
    fake_storage.map_ancestors(job)
    return job


def upload(client: Client, job: Job, *files: tuple[str, bytes], **fields: str) -> Any:
    sent = [SimpleUploadedFile(name, data) for name, data in files or [("my model.pdb", STRUCTURE)]]
    return client.post(f"/jobs/{job.id}/structure/upload/", {"file": sent, **fields})


def data_of(html: str) -> dict[str, Any]:
    found = re.search(r'<script id="structure-data" type="application/json">(.*?)</script>', html)
    assert found, "the page has no structure data"
    loaded: dict[str, Any] = json.loads(found.group(1))
    return loaded


def query_of(url: str) -> dict[str, list[str]]:
    return parse_qs(urlsplit(url).query)


def test_a_job_without_structures_offers_to_attach_some(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = structure_job(fake_storage)
    html = client.get(f"/jobs/{job.id}/structure/").content.decode()
    assert "<h1>Ancestors on Structure</h1>" in html
    assert "No structure is attached to this job yet." in html and "structure-data" not in html
    # the job's sequence names that carry a UniProt accession are offered
    assert 'data-structure-suggest="P12931" title="a_P12931"' in html
    # the section sits in the viewer's own columns, folded like its sections
    assert '<details class="viewer-section" id="structures" open>' in html
    # an ancestor to model elsewhere: its sequence to copy, gaps removed
    assert f'data-structure-copy="{ANCESTORS["Node1"]}"' in html
    # the ancestor viewer and the job page lead here
    viewer = client.get(f"/jobs/{job.id}/ancestors/").content.decode()
    assert "View on Structure</a>" in viewer and f"/jobs/{job.id}/structure/?combination=" in viewer
    assert "Attach a structure</a>" in client.get(f"/jobs/{job.id}/").content.decode()


def test_an_uploaded_structure_is_checked_then_painted_with_an_ancestor(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = structure_job(fake_storage)
    refused = upload(client, job, ("x.pdb", b"<html><script>alert(1)</script></html>"))
    assert refused.status_code == 302 and not Structure.objects.exists()
    response = upload(client, job, label="  Anc   model ", predicted="on")
    structure = Structure.objects.get()
    assert query_of(response["Location"]) == {"ps": [str(structure.id)], "pw": ["a"]}
    assert (structure.label, structure.source, structure.predicted) == ("Anc model", "upload", True)
    assert structure.format == "pdb"
    assert structure.chains == [{"id": "A", "residues": 16}, {"id": "B", "residues": 10}]
    assert fake_storage.objects[f"jobs/{job.id}/structures/{structure.id}.pdb"] == STRUCTURE

    served = client.get(f"/jobs/{job.id}/structure/{structure.id}/file")
    assert served.content == STRUCTURE and served["Content-Type"].startswith("text/plain")
    assert served["X-Content-Type-Options"] == "nosniff"

    html = client.get(f"/jobs/{job.id}/structure/").content.decode()
    data = data_of(html)
    (panel,) = data["panels"]
    assert panel["file"] == f"/jobs/{job.id}/structure/{structure.id}/file"
    assert panel["format"] == "pdb" and panel["matrix"] is None and data["scheme"] == "pp"
    assert panel["keys"][0] == "A:10" and len(panel["keys"]) == 26 and panel["cols"][0] == 2
    assert "column 3 · Node1: T (probability 0.80)" in panel["tips"][0]
    assert "not in the alignment" in panel["tips"][-1] and panel["cols"][-1] is None
    assert panel["colours"][-1] == colours.UNMAPPED != panel["colours"][0]
    assert "read through <strong>a_P12931</strong>" in html and "16 of 16 residues aligned" in html
    assert 'id="structure-stage-0" hx-preserve="true"' in html
    # unticked boxes still say so (a hidden "off" before each), and the "?" explains the schemes
    assert '<input type="hidden" name="surface" value="off">' in html
    assert "What the Colours Show" in html and "Red sites" in html

    # two sequences: what changed between them, on the structure and in the table
    response = client.get(
        f"/jobs/{job.id}/structure/", {"node": "Node1", "to": "Node2", "labels": "on"},
        headers={"HX-Request": "1"},
    )  # fmt: skip
    html = response.content.decode()
    pushed = query_of(response["HX-Push-Url"])
    assert pushed["to"] == ["Node2"] and pushed["scheme"] == ["changes"]
    assert pushed["ps"] == [str(structure.id)] and pushed["pw"] == ["a"] and pushed["pc"] == ["A"]
    assert "<h1>" not in html  # the viewer alone
    data = data_of(html)
    (panel,) = data["panels"]
    assert [panel["keys"][i] for i in panel["marked"]] == ["A:13", "A:22"]
    assert data["labels"] == "on" and panel["labels_text"] == {"3": "I6L", "12": "V15I"}
    assert panel["colours"][3] == colours.SUBSTITUTED and panel["colours"][0] == colours.SAME
    assert 'data-structure-col="5">I6L</button>' in html
    assert "2 sites differ between <strong>Node1</strong> and" in html
    # an extant sequence as the second one, two branches down: the path's steps are offered
    response = client.get(f"/jobs/{job.id}/structure/", {"to": "b", "scheme": "steps"})
    html = response.content.decode()
    assert "R10-" in html and "Q11-" in html and "They happened along 2 branches." in html
    data = data_of(html)
    assert data["scheme"] == "steps" and "I6L on Node1 → Node2" in data["panels"][0]["tips"][3]
    assert "L6I on Node2 → b" in data["panels"][0]["tips"][3]  # and back again
    html = client.get(f"/jobs/{job.id}/structure/", {"to": "<nope>"}).content.decode()
    assert data_of(html)["scheme"] == "pp" and "<nope>" not in html

    dropped = client.post(f"/jobs/{job.id}/structure/{structure.id}/drop/")
    assert dropped.status_code == 302 and not Structure.objects.exists()
    assert not any("structures/" in key for key in fake_storage.objects)


def test_sequences_and_structures_are_mixed_on_panels(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = structure_job(fake_storage)
    # two files at once: a crystal structure, and a model of Node2 as AlphaFold writes one
    upload(client, job, ("crystal.pdb", STRUCTURE), ("fold_node2_model_0.cif", MODEL))
    crystal, model = Structure.objects.order_by("format").reverse()
    assert (crystal.format, model.format) == ("pdb", "cif")
    assert model.predicted and not crystal.predicted  # the mmCIF says it is a model
    assert model.title == "Model of Node2" and model.label == "fold_node2_model_0"
    assert f"jobs/{job.id}/structures/{model.id}.cif" in fake_storage.objects
    url = f"/jobs/{job.id}/structure/"
    both = {"node": "Node1", "to": "Node2"}

    # each end of the branch on the structure closest to it: Node2 on its own model
    response = client.get(url, both | {"preset": "closest"})
    assert response.status_code == 302  # a one-off button: the address gets the view
    view = query_of(response["Location"])
    assert view["pw"] == ["a", "b"] and view["ps"][1] == str(model.id) and "preset" not in view
    html = client.get(response["Location"]).content.decode()
    assert "This structure is Node2 itself." in html and "is this sequence" in html

    # the two ends on different structures: laid on each other, sharing a camera
    asked = both | {"ps": [str(crystal.id), str(model.id)], "pw": ["a", "b"], "pc": ["", ""]}
    html = client.get(url, asked).content.decode()
    first, second = data_of(html)["panels"]
    assert first["matrix"] is None and len(second["matrix"]) == 12 and second["linked"]
    assert (first["format"], second["format"]) == ("pdb", "cif")
    assert "Laid on panel 1 by 16 of 16 shared residues (RMSD 0.0 Å)" in html
    assert 'id="structure-stage-1" hx-preserve="true"' in html
    # the strips scroll together by alignment column unless unlocked
    assert 'name="scroll" value="on" checked' in html and 'data-key="A:13" data-col="5"' in html
    assert data_of(client.get(url, asked | {"scroll": "off"}).content.decode())["scroll"] == "off"
    # the same change in both, each on its own residue, named first to second in both
    assert [first["keys"][i] for i in first["marked"]] == ["A:13", "A:22"]
    assert [second["keys"][i] for i in second["marked"]] == ["A:6", "A:15"]
    assert set(first["labels_text"].values()) == set(second["labels_text"].values())
    assert second["title"] == "Node2 on fold_node2_model_0" and second["cols"][5] == 5
    assert "<td>A ILE 13</td><td>A LEU 6</td>" in html
    # the panels can be freed from each other
    html = client.get(url, asked | {"link": "off"}).content.decode()
    assert data_of(html)["panels"][1]["matrix"] is None and "Laid on panel 1" not in html

    # one ancestor on every structure; a panel added, a panel dropped
    view = query_of(client.get(url, {"preset": "every"})["Location"])
    assert sorted(view["ps"]) == sorted([str(crystal.id), str(model.id)]) and view["pw"] == [
        "a",
        "a",
    ]
    view = query_of(client.get(url, both | {"ps": [str(crystal.id)], "add": "1"})["Location"])
    assert view["ps"] == [str(crystal.id)] * 2 and view["pw"] == ["a", "b"]  # the other end
    view = query_of(client.get(url, asked | {"drop": "0"})["Location"])
    assert view["ps"] == [str(model.id)] and view["pw"] == ["b"]
    # an unknown structure or sequence in the address is dropped
    view = query_of(client.get(url, {"ps": ["nope", str(model.id)], "pw": ["a", "b"], "add": "x"})[
        "Location"
    ])  # fmt: skip
    assert view["ps"][0] == str(model.id) and set(view["pw"]) == {"a"}


def test_a_predictors_zip_gives_its_top_ranked_model(
    client: Client, fake_storage: FakeStorage
) -> None:
    job = structure_job(fake_storage)
    packed = io.BytesIO()
    with zipfile.ZipFile(packed, "w") as archive:
        for rank in range(3):
            archive.writestr(f"fold_node2/fold_node2_model_{rank}.cif", MODEL)
            archive.writestr(f"fold_node2/fold_node2_full_data_{rank}.json", "{}")
        archive.writestr("../../evil.cif", MODEL)  # a name is only ever read, never a path
        archive.writestr("terms_of_use.md", "terms")
    upload(client, job, ("fold_node2.zip", packed.getvalue()))
    (model,) = Structure.objects.all()
    assert (model.label, model.format, model.predicted) == ("fold_node2_model_0", "cif", True)
    assert [k for k in fake_storage.objects if "structures/" in k] == [
        f"jobs/{job.id}/structures/{model.id}.cif"
    ]
    for bad in (b"not a zip", b"PK\x03\x04 broken"):
        upload(client, job, ("x.zip", bad))
    empty = io.BytesIO()
    with zipfile.ZipFile(empty, "w") as archive:
        archive.writestr("notes.txt", "nothing")
    upload(client, job, ("empty.zip", empty.getvalue()))
    assert Structure.objects.count() == 1


def test_attaching_by_identifier_sends_only_the_identifiers(
    client: Client, fake_storage: FakeStorage, monkeypatch: pytest.MonkeyPatch
) -> None:
    job = structure_job(fake_storage)
    asked: list[str] = []

    def answer(url: str, data: object = None) -> bytes:
        asked.append(url)
        if "1BAD" in url:
            raise homolog_remote.SearchError("the search service answered 404: nothing")
        if url.endswith(".pdb"):
            return STRUCTURE
        return json.dumps({"struct": {"title": "A kinase <b>domain</b>"}}).encode()

    monkeypatch.setattr(homolog_remote, "request", answer)
    url = f"/jobs/{job.id}/structure/attach/"
    client.post(url, {"identifier": "2src"})  # no consent: nothing is sent
    client.post(url, {"identifier": "2SRC, ../../etc/passwd", "consent": "on"})
    client.post(url, {"pick": ["rcsb:2SRC", "rcsb:https://evil.example"], "consent": "on"})
    assert asked == [] and not Structure.objects.exists()
    response = client.post(url, {"identifier": "2src", "consent": "on"})
    structure = Structure.objects.get()
    assert asked == [
        "https://files.rcsb.org/download/2SRC.pdb",
        "https://data.rcsb.org/rest/v1/core/entry/2SRC",
    ]
    assert (structure.label, structure.source, structure.accession) == ("2SRC", "rcsb", "2SRC")
    assert not structure.predicted and str(structure.id) in response["Location"]
    # several ticked in the lookup's results: the new ones fetched, the rest reported
    del asked[:]
    response = client.post(
        url, {"pick": ["rcsb:2SRC", "rcsb:1OPL", "rcsb:1BAD", "rcsb:1FMK"], "consent": "on"}
    )
    assert [u for u in asked if u.endswith(".pdb")] == [
        f"https://files.rcsb.org/download/{entry}.pdb" for entry in ("1OPL", "1BAD", "1FMK")
    ]
    assert sorted(s.label for s in Structure.objects.all()) == ["1FMK", "1OPL", "2SRC"]
    assert "ps=" not in response["Location"]  # several: the view stays as it was
    html = client.get(f"/jobs/{job.id}/structure/").content.decode()
    assert "Attached 1OPL, 1FMK." in html and "Already attached: 2SRC." in html
    assert "Not attached: 1BAD (RCSB has no entry 1BAD)." in html
    assert "A kinase &lt;b&gt;domain&lt;/b&gt;" in html and "<b>domain" not in html


def test_finding_structures_by_name(
    client: Client, fake_storage: FakeStorage, monkeypatch: pytest.MonkeyPatch
) -> None:
    job = structure_job(fake_storage)
    asked: list[str] = []
    entry = {
        "primaryAccession": "P12931", "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "organism": {"scientificName": "Homo sapiens"}, "sequence": {"length": 536},
        "uniProtKBCrossReferences": [
            {"database": "PDB", "id": "2SRC", "properties": [
                {"key": "Method", "value": "X-ray"}, {"key": "Resolution", "value": "1.50 A"}]},
            {"database": "AlphaFoldDB", "id": "P12931"}],
    }  # fmt: skip

    def answer(url: str, data: object = None) -> bytes:
        asked.append(url)
        return json.dumps({"results": [entry]}).encode()

    monkeypatch.setattr(homolog_remote, "request", answer)
    url = f"/jobs/{job.id}/structure/find/"
    assert "Tick this box" in client.get(url, {"protein": "SRC"}).content.decode()
    assert asked == []  # nothing leaves without the consent tick
    found = client.get(url, {"protein": 'SRC" OR x', "species": "Homo sapiens", "consent": "on"})
    html = found.content.decode()
    assert len(asked) == 1 and "xref_pdb" in asked[0] and "%22SRC+OR+x%22" in asked[0]
    # every structure is a box to tick: any number are attached together
    assert '<input type="checkbox" name="pick" value="rcsb:2SRC"' in html
    assert '<input type="checkbox" name="pick" value="alphafold:P12931">' in html
    assert html.count("data-structure-attach") == 1
