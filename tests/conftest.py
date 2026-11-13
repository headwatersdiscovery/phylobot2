import os
import re
from typing import Any

import pytest

# Settings require a secret key; tests that don't touch services run on the host.
os.environ.setdefault("DJANGO_SECRET_KEY", "test-only-not-secret")
os.environ.setdefault("PHYLOBOT_MODE", "local")


RESULT_KEY = re.compile(r"^jobs/([0-9a-f-]{36})/results/(.+)$")


class Objects(dict[str, bytes]):
    """The fake's objects. A file stored under a job's results/ is also recorded in its
    manifest (phylobot.manifest), as the worker does when it uploads, so tests that put
    results in place see them on every page; `indexing = False` turns that off."""

    indexing = True

    def __setitem__(self, key: str, data: bytes) -> None:
        super().__setitem__(key, data)
        if self.indexing and (found := RESULT_KEY.match(key)):
            from django.utils import timezone

            from phylobot.models import Job, ResultFile

            if Job.objects.filter(id=found[1]).exists():
                ResultFile.objects.update_or_create(
                    job_id=found[1], path=found[2], defaults={"size": len(data)}
                )
                Job.objects.filter(id=found[1], results_indexed_at__isnull=True).update(
                    results_indexed_at=timezone.now()
                )

    def __delitem__(self, key: str) -> None:
        super().__delitem__(key)
        if self.indexing and (found := RESULT_KEY.match(key)):
            from phylobot.models import ResultFile

            ResultFile.objects.filter(job_id=found[1], path=found[2]).delete()

    def pop(self, key: str, *default: Any) -> Any:
        if key in self:
            data = self[key]
            del self[key]
            return data
        return default[0] if default else super().pop(key)


class FakeStorage:
    """In-memory stand-in for phylobot.storage in the web views."""

    def __init__(self) -> None:
        self.objects = Objects()
        self.measured: dict[tuple[str, str], dict[tuple[str, str], Any]] = {}

    def put_bytes(self, key: str, data: bytes) -> None:
        self.objects[key] = data

    def put_file(self, key: str, path: Any) -> None:
        self.objects[key] = path.read_bytes()

    def get_file(self, key: str, path: Any) -> None:
        path.write_bytes(self.objects[key])

    def delete_key(self, key: str) -> None:
        self.objects.pop(key, None)

    def list_keys(self, prefix: str) -> list[str]:
        return sorted(k for k in self.objects if k.startswith(prefix))

    def list_sizes(self, prefix: str) -> dict[str, int]:
        return {k: len(self.objects[k]) for k in self.list_keys(prefix)}

    def exists(self, key: str) -> bool:
        return key in self.objects

    def read_bytes(self, key: str) -> bytes:
        return self.objects[key]

    def size(self, key: str) -> int:
        return len(self.objects[key])

    def read_range(self, key: str, start: int, length: int) -> bytes:
        return self.objects[key][start : start + length]

    def stream(self, key: str) -> Any:
        yield self.objects[key]

    def copy_key(self, source: str, dest: str) -> None:
        self.objects[dest] = self.objects[source]

    def copy_prefix(self, source: str, dest: str) -> int:
        keys = self.list_keys(source)
        for key in keys:
            self.objects[dest + key[len(source) :]] = self.objects[key]
        return len(keys)

    def prefix_size(self, prefix: str) -> int:
        return sum(len(self.objects[k]) for k in self.list_keys(prefix))

    def map_ancestors(self, job: Any) -> None:
        """Write the job's ancestor_map.tsv from its stored trees, as the pipeline's
        ANCESTOR_MAP does (each combination's tree.rooted.nwk if it has one)."""
        from phylobot.pipeline.ancestor_map import build

        prefix = f"{job.storage_prefix}results/"
        lines = []
        for key in self.list_keys(prefix):
            relative = key[len(prefix) :]
            combination, _, name = relative.rpartition("/")
            if name == "tree.rooted.nwk" or (
                name == "tree.nwk" and f"{prefix}{combination}/tree.rooted.nwk" not in self.objects
            ):
                kind = "rooted" if name == "tree.rooted.nwk" else "unrooted"
                lines.append(f"{combination}\t{kind}\t{self.objects[key].decode().strip()}")
        self.objects[prefix + "ancestor_map.tsv"] = build("\n".join(lines)).encode()

    def store_properties(
        self,
        results: str,
        combination: str,
        node: str,
        version: str,
        summary: dict[str, Any],
        draws: dict[str, list[float]],
    ) -> None:
        """Add one node's measurements to a combination's packed property files under
        `results` (a job's results/ prefix), as PROPERTIES writes them."""
        import tempfile
        from pathlib import Path

        from phylobot import property_store

        found = self.measured.setdefault((results, combination), {})
        found[(version, node)] = property_store.Measured(version, node, draws, summary)
        keys = sorted({k for m in found.values() for k in m.values})
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            draws_n = max((len(v) for m in found.values() for v in m.values.values()), default=1)
            property_store.write(out, list(found.values()), keys, draws_n)
            for name in (property_store.DRAWS, property_store.SUMMARIES):
                key = f"{results}{combination}/{property_store.FOLDER}/{name}"
                self.objects[key] = (out / name).read_bytes()
        from phylobot import property_views

        property_views._summaries.cache_clear()  # finished jobs' never change; tests' do

    def delete_prefix(self, prefix: str) -> int:
        doomed = [k for k in self.objects if k.startswith(prefix)]
        for key in doomed:
            del self.objects[key]
        return len(doomed)


@pytest.fixture(autouse=True)
def no_open_tree(monkeypatch: pytest.MonkeyPatch, tmp_path_factory: pytest.TempPathFactory) -> None:
    """No test reaches the Open Tree of Life, NCBI or UniProt: the web tier's clients
    fail unless a test gives them a recording (species_web.make_transport,
    make_accession_transport), and cache nowhere shared."""
    from phylobot import species_web

    def refuse() -> Any:
        raise AssertionError("a test tried to reach the Open Tree of Life")

    monkeypatch.setattr(species_web, "make_transport", refuse)
    monkeypatch.setattr(species_web, "make_accession_transport", refuse)
    monkeypatch.setattr(species_web, "CACHE_DIR", tmp_path_factory.mktemp("opentree"))
    species_web.reset_client()


@pytest.fixture(autouse=True)
def no_real_storage(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> None:
    """No fast test reaches the real object store: without this, a test that submits or
    deletes a job through the views wrote into the developer's own Garage bucket (and
    failed wherever Garage was not running). Use the `fake_storage` fixture instead;
    golden and e2e tests drive the real stack on purpose."""
    if any(request.node.get_closest_marker(m) for m in ("golden", "golden_full", "e2e")):
        return
    from phylobot import storage

    def refuse() -> Any:
        raise AssertionError("a test tried to reach real object storage: use fake_storage")

    monkeypatch.setattr(storage, "client", refuse)


@pytest.fixture(autouse=True)
def fixed_machine(monkeypatch: pytest.MonkeyPatch) -> None:
    """Every test's Docker has 8 CPUs and 16 GiB (phylobot.machine.detect asks it at each
    job launch); tests of detection itself replace this."""
    from phylobot import machine

    monkeypatch.setattr(machine, "docker_info", lambda: (8, 16 * 1024**3))


@pytest.fixture
def fake_storage(monkeypatch: pytest.MonkeyPatch) -> FakeStorage:
    from phylobot import property_views

    property_views._summaries.cache_clear()
    fake = FakeStorage()
    for name in (
        "put_bytes",
        "put_file",
        "get_file",
        "delete_key",
        "list_keys",
        "list_sizes",
        "exists",
        "read_bytes",
        "read_range",
        "size",
        "stream",
        "copy_key",
        "copy_prefix",
        "delete_prefix",
        "prefix_size",
    ):
        monkeypatch.setattr(f"phylobot.storage.{name}", getattr(fake, name))
    return fake
