"""A small client for the Open Tree of Life API, v3 (https://api.opentreeoflife.org).

Our own client (urllib), not the `opentree` package (unmaintained since 2022) and not
ete (GPL-3, which must not be imported into our BSD-3 package). Only species names and
NCBI taxids are sent, never sequences; the caller asks the user first. Open Tree
publishes no rate limit, so requests are paced (`MIN_GAP_SECONDS`) and answers can be
cached on disk (`Cache`), the induced subtree keyed by the synthesis release.

Every call goes through a `Transport`, so tests replay recorded answers
(tests/fixtures/opentree/) without the network.

Behaviour found by testing (2026-10-16, synthesis opentree16.1), not documented:
- one pruned id (e.g. a "barren" or "merged" taxon) fails the whole induced-subtree
  request with HTTP 400 and `{"unknown": {"ott<id>": "pruned_ott_id"}}`;
- a non-monophyletic ("broken") taxon is answered by its stand-in, in `broken`
  (`{"ott632176": "ott985906"}`: Crustacea is drawn as Pancrustacea);
- a requested higher taxon is an internal node if other requested taxa lie inside it,
  and a tip otherwise;
- name matching in a narrow context silently drops names, so the context is always
  "All life".
"""

import hashlib
import json
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

BASE_URL = "https://api.opentreeoflife.org/v3"
TIMEOUT = 60  # seconds without data
TOTAL_TIMEOUT = 180  # seconds for a whole answer
MAX_BYTES = 20_000_000  # an induced subtree of a few thousand tips is well under this
CHUNK = 1 << 16
MIN_GAP_SECONDS = 1.0
MAX_NAMES = 5000  # per match_names request (Open Tree accepts 10,000)
CACHE_DAYS = 30
CONTEXT = "All life"


class OpenTreeError(Exception):
    pass


@dataclass(frozen=True)
class Answer:
    status: int  # HTTP status
    body: Any  # parsed JSON


# (endpoint, JSON payload) -> Answer. Raises OpenTreeError when the service cannot be
# reached or does not answer in JSON.
Transport = Callable[[str, dict[str, Any]], Answer]


class Http:
    """The real transport, paced: at most one request per `MIN_GAP_SECONDS`."""

    def __init__(self, base_url: str = BASE_URL, min_gap: float = MIN_GAP_SECONDS):
        self.base_url = base_url
        self.min_gap = min_gap
        self.last = 0.0

    def __call__(self, endpoint: str, payload: dict[str, Any]) -> Answer:
        wait = self.last + self.min_gap - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        req = urllib.request.Request(
            f"{self.base_url}/{endpoint}",
            data=json.dumps(payload).encode(),
            headers={"Content-Type": "application/json", "User-Agent": "PhyloBot"},
        )
        try:
            try:
                with urllib.request.urlopen(req, timeout=TIMEOUT) as response:
                    status, raw = response.status, read_bounded(response)
            except urllib.error.HTTPError as exc:
                status, raw = exc.code, read_bounded(exc)
        except (urllib.error.URLError, TimeoutError) as exc:
            raise OpenTreeError(f"could not reach the Open Tree of Life: {exc}") from exc
        finally:
            self.last = time.monotonic()
        try:
            return Answer(status, json.loads(raw))
        except ValueError:
            raise OpenTreeError(
                f"the Open Tree of Life answered {status} without JSON: {raw[:200]!r}"
            ) from None


def read_bounded(response: Any) -> bytes:
    deadline = time.monotonic() + TOTAL_TIMEOUT
    parts: list[bytes] = []
    size = 0
    while chunk := response.read(CHUNK):
        parts.append(chunk)
        size += len(chunk)
        if size > MAX_BYTES:
            raise OpenTreeError(f"the Open Tree of Life answer is over {MAX_BYTES} bytes")
        if time.monotonic() > deadline:
            raise OpenTreeError(f"the Open Tree of Life took over {TOTAL_TIMEOUT} s to answer")
    return b"".join(parts)


class Cache:
    """Answers on disk, one JSON file per (endpoint, payload), kept `days` days."""

    def __init__(self, directory: Path, days: float = CACHE_DAYS):
        self.directory = directory
        self.seconds = days * 86400

    def path(self, endpoint: str, payload: dict[str, Any]) -> Path:
        key = json.dumps([endpoint, payload], sort_keys=True)
        return self.directory / f"{hashlib.sha256(key.encode()).hexdigest()}.json"

    def get(self, endpoint: str, payload: dict[str, Any]) -> Any | None:
        path = self.path(endpoint, payload)
        try:
            if time.time() - path.stat().st_mtime > self.seconds:
                return None
            return json.loads(path.read_text())
        except OSError, ValueError:
            return None

    def put(self, endpoint: str, payload: dict[str, Any], body: Any) -> None:
        self.directory.mkdir(parents=True, exist_ok=True)
        path = self.path(endpoint, payload)
        partial = path.with_suffix(".part")
        partial.write_text(json.dumps(body))
        partial.replace(path)


@dataclass(frozen=True)
class Taxon:
    ott_id: int
    name: str
    rank: str
    flags: tuple[str, ...]
    in_synthesis: bool  # False: suppressed from the synthetic tree
    ncbi: tuple[int, ...]  # NCBI taxids among its sources


@dataclass(frozen=True)
class Match:
    taxon: Taxon
    synonym: bool  # matched one of the taxon's synonyms (e.g. an older name)
    approximate: bool
    score: float


@dataclass(frozen=True)
class Induced:
    newick: str  # labels in Open Tree's "name_and_id" format
    broken: dict[int, str]  # requested ott id -> the node id standing in for it
    pruned: tuple[int, ...]  # requested ids missing from the synthesis (left out)
    studies: tuple[str, ...]  # supporting studies


def taxon(entry: dict[str, Any]) -> Taxon:
    ncbi = []
    for source in entry.get("tax_sources") or []:
        prefix, _, number = str(source).partition(":")
        if prefix == "ncbi" and number.isdigit():
            ncbi.append(int(number))
    return Taxon(
        ott_id=int(entry["ott_id"]),
        name=str(entry["name"]),
        rank=str(entry.get("rank") or ""),
        flags=tuple(str(f) for f in entry.get("flags") or []),
        in_synthesis=not entry.get("is_suppressed_from_synth", False),
        ncbi=tuple(ncbi),
    )


class Client:
    def __init__(self, transport: Transport | None = None, cache: Cache | None = None):
        self.transport = transport or Http()
        self.cache = cache
        self._synthesis: dict[str, str] | None = None

    def request(self, endpoint: str, payload: dict[str, Any], version: str = "") -> Answer:
        """One call; with a `version` (the release the answer depends on), answered from
        the cache when it can be, and successful answers cached under that version."""
        key = {**payload, "_version": version}
        if version and self.cache and (body := self.cache.get(endpoint, key)) is not None:
            return Answer(200, body)
        answer = self.transport(endpoint, payload)
        if version and self.cache and answer.status == 200:
            self.cache.put(endpoint, key, answer.body)
        return answer

    def synthesis(self) -> dict[str, str]:
        """The synthesis release ("synth_id", "taxonomy_version", "date_created").
        Asked once per client and never cached: it versions the other answers."""
        if self._synthesis is None:
            answer = self.request("tree_of_life/about", {})
            if answer.status != 200 or not isinstance(answer.body, dict):
                raise OpenTreeError(message("tree_of_life/about", answer))
            self._synthesis = {
                key: str(answer.body.get(key) or "")
                for key in ("synth_id", "taxonomy_version", "date_created")
            }
            if not self._synthesis["synth_id"]:
                raise OpenTreeError("the Open Tree of Life did not name its synthesis release")
        return self._synthesis

    def match_names(self, names: list[str], approximate: bool = False) -> dict[str, list[Match]]:
        """{name: matches}; exact matches unless `approximate` (suggestions)."""
        result: dict[str, list[Match]] = {name: [] for name in names}
        unique = sorted(set(names))
        version = self.synthesis()["taxonomy_version"]
        for start in range(0, len(unique), MAX_NAMES):
            payload = {
                "names": unique[start : start + MAX_NAMES],
                "do_approximate_matching": approximate,
                "context_name": CONTEXT,
            }
            answer = self.request("tnrs/match_names", payload, version)
            if answer.status != 200:
                raise OpenTreeError(message("tnrs/match_names", answer))
            for entry in answer.body.get("results") or []:
                name = entry.get("name")
                if name in result:
                    result[name] = [
                        Match(
                            taxon=taxon(m["taxon"]),
                            synonym=bool(m.get("is_synonym")),
                            approximate=bool(m.get("is_approximate_match")),
                            score=float(m.get("score") or 0),
                        )
                        for m in entry.get("matches") or []
                    ]
        return result

    def taxon_from_ncbi(self, taxid: int) -> Taxon | None:
        """The Open Tree taxon with this NCBI taxid; None if there is none."""
        payload = {"source_id": f"ncbi:{int(taxid)}"}
        answer = self.request("taxonomy/taxon_info", payload, self.synthesis()["taxonomy_version"])
        if answer.status == 400 and "No taxon" in str(answer.body.get("message", "")):
            return None
        if answer.status != 200:
            raise OpenTreeError(message("taxonomy/taxon_info", answer))
        return taxon(answer.body)

    def taxon_by_ott(self, ott_id: int) -> Taxon | None:
        """The Open Tree taxon with this OTT id; None if there is none."""
        payload = {"ott_id": int(ott_id)}
        answer = self.request("taxonomy/taxon_info", payload, self.synthesis()["taxonomy_version"])
        if answer.status == 400:
            return None
        if answer.status != 200:
            raise OpenTreeError(message("taxonomy/taxon_info", answer))
        return taxon(answer.body)

    def parent(self, ott_id: int) -> Taxon | None:
        """The taxon's parent in the Open Tree taxonomy (its lineage's first entry)."""
        payload = {"ott_id": int(ott_id), "include_lineage": True}
        answer = self.request("taxonomy/taxon_info", payload, self.synthesis()["taxonomy_version"])
        if answer.status != 200 or not isinstance(answer.body, dict):
            return None
        lineage = answer.body.get("lineage") or []
        return taxon(lineage[0]) if lineage and isinstance(lineage[0], dict) else None

    def induced_subtree(self, ott_ids: list[int], attempts: int = 3) -> Induced:
        """The synthetic tree induced by these taxa. Ids Open Tree has pruned from the
        synthesis fail the whole request; they are dropped and the request retried."""
        wanted = sorted(set(ott_ids))
        pruned: list[int] = []
        version = self.synthesis()["synth_id"]
        for _ in range(attempts):
            if len(wanted) < 2:
                raise OpenTreeError("fewer than two of the species are in the Open Tree synthesis")
            payload = {"ott_ids": wanted, "label_format": "name_and_id"}
            answer = self.request("tree_of_life/induced_subtree", payload, version)
            unknown = answer.body.get("unknown") if isinstance(answer.body, dict) else None
            if answer.status == 400 and isinstance(unknown, dict) and unknown:
                gone = {ott_number(k) for k in unknown}
                if None in gone or not gone <= set(wanted):
                    raise OpenTreeError(message("tree_of_life/induced_subtree", answer))
                pruned.extend(sorted(i for i in gone if i is not None))
                wanted = [i for i in wanted if i not in gone]
                continue
            if answer.status != 200:
                raise OpenTreeError(message("tree_of_life/induced_subtree", answer))
            newick = answer.body.get("newick")
            if not isinstance(newick, str) or not newick.strip().endswith(";"):
                raise OpenTreeError("the Open Tree of Life answered without a tree")
            broken: dict[int, str] = {}
            for key, value in (answer.body.get("broken") or {}).items():
                if (number := ott_number(key)) is not None:
                    broken[number] = str(value)
            return Induced(
                newick=newick.strip(),
                broken=broken,
                pruned=tuple(pruned),
                studies=tuple(str(s) for s in answer.body.get("supporting_studies") or []),
            )
        raise OpenTreeError("the Open Tree of Life kept refusing taxa; giving up")


def ott_number(node_id: str) -> int | None:
    """ "ott123" -> 123; anything else (e.g. "mrcaott1ott2") -> None."""
    digits = str(node_id).removeprefix("ott")
    return int(digits) if digits.isdigit() and str(node_id).startswith("ott") else None


def message(endpoint: str, answer: Answer) -> str:
    detail = answer.body.get("message") if isinstance(answer.body, dict) else None
    return f"the Open Tree of Life ({endpoint}) answered {answer.status}: {str(detail)[:300]}"
