"""The tree viewer's "Node labels" dropdown (hide all / show all / show named only): it
governs the node name labels (NodeN and user names) drawn beside each node; branch support
labels keep their own control."""

import pytest
from django.test import Client

from phylobot.models import Job
from tests.conftest import FakeStorage
from tests.test_views import viewer_job

pytestmark = pytest.mark.django_db

WAG = "mafft/iqtree/WAG"


def page(client: Client, job: Job, **extra: str) -> str:
    query = {"combination": WAG, "aligner": "mafft", **extra}
    return client.get(f"/jobs/{job.id}/tree/", query).content.decode()


def name(client: Client, job: Job, node: str, text: str) -> None:
    client.post(
        f"/jobs/{job.id}/ancestors/name/",
        {
            "combination": WAG,
            "aligner": "mafft",
            "gaps": "original",
            "node": node,
            "ancestor_name": text,
        },  # fmt: skip
    )


def test_the_dropdown_is_offered_with_three_options(
    client: Client, fake_storage: FakeStorage
) -> None:
    html = page(client, job := viewer_job(fake_storage))
    assert job and 'name="node_labels"' in html
    for value in ("all", "named", "none"):
        assert f'value="{value}"' in html
    assert 'value="all" selected' in html  # the default


def test_show_all_hide_all_and_named_only(client: Client, fake_storage: FakeStorage) -> None:
    job = viewer_job(fake_storage)
    name(client, job, "Node2", "Anc.AB")
    named = 'font-weight="bold">Anc.AB'  # the drawn custom name
    plain = 'paint-order="stroke">Node1'  # Node1's own (unnamed) node label

    shown = page(client, job)  # default is "all"
    assert named in shown and plain in shown

    only_named = page(client, job, node_labels="named")
    assert named in only_named and plain not in only_named  # unnamed node's label hidden

    none = page(client, job, node_labels="none")
    assert named not in none and plain not in none

    # The choice round-trips in the control (so a shared URL reopens the same way).
    assert 'value="named" selected' in only_named


def test_a_crowded_tree_starts_with_named_labels_only(
    client: Client, fake_storage: FakeStorage, monkeypatch: pytest.MonkeyPatch
) -> None:
    """NodeN beside every node of a large tree collides with its neighbours, so unless the
    reader asked, a crowded tree labels only the ancestors they named; asking for all still
    gives all, and keeps it in the URL the viewer's links carry."""
    from phylobot import treeview

    job = viewer_job(fake_storage)
    name(client, job, "Node2", "Anc.AB")
    monkeypatch.setattr(treeview, "CROWDED_NODES", 1)
    plain = 'paint-order="stroke">Node1'
    crowded = page(client, job)
    assert 'value="named" selected' in crowded
    assert 'font-weight="bold">Anc.AB' in crowded and plain not in crowded
    assert "hover one for its name" in crowded
    asked = page(client, job, node_labels="all")
    assert 'value="all" selected' in asked and plain in asked
