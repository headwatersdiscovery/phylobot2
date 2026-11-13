"""A job's results as its pages read them: its own files, plus, for a job that extends
another (`phylobot.extend`), the parent's files for the alignments it inherits.

Nothing is copied. A child that adds MUSCLE to a MAFFT job ran only MUSCLE; its pages
read `mafft/...` from the parent's results and `muscle/...` from its own. A child that
adds a model or a tree source ran only the combinations its parent lacked, so a folder
like `mafft/` holds both jobs' files: ownership is per file, the job's own first, then
the parent's. A child of a child resolves through its parent, so ownership composes.

The tables every alignment shares are joined here: summary.json and
model_comparison.tsv from the joined combinations, re-ranked by BIC within each
alignment as `pipeline/summary.py` ranks one job's (a model added to an inherited
alignment changes that alignment's ranking), and sources.tsv and versions.tsv as the
rows of both. Reconciled jobs are never extended (their species tree is resolved from
every tree, `extend.refusal`).

Every page that reads a job's results goes through `results`, so an extension looks like
one job that ran every aligner, and its ZIP holds every file.
"""

import json
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any

from phylobot import manifest, storage
from phylobot.models import Job
from phylobot.pipeline.export import FIT_FIELDS
from phylobot.pipeline.summary import COMPARISON_FIELDS, compare

SUMMARY = "summary.json"
COMPARISON = "model_comparison.tsv"
# rows of both
JOINED_TABLES = ("model_comparison.tsv", "sources.tsv", "versions.tsv", "ancestor_map.tsv")
READERS = 16  # concurrent reads in Results.read_many


@dataclass
class Results:
    """A job's result files by path relative to results/, wherever each one lives."""

    job: Job
    keys: set[str] = field(default_factory=set)
    summary: dict[str, Any] = field(default_factory=dict)
    # An inherited file -> the job whose results/ holds it (absent: this job).
    owners: dict[str, Job] = field(default_factory=dict)
    comparison: list[dict[str, str]] = field(default_factory=list)  # the joined ranking
    # The parent's tables, for the joined ones (empty unless the job inherits).
    inherited_tables: dict[str, str] = field(default_factory=dict)
    inherited: list[str] = field(default_factory=list)  # the alignments it inherits
    lost: bool = False  # it inherits from a job whose files are gone
    memo: dict[str, Any] = field(default_factory=dict)  # what readers parsed, per request

    def owner(self, relative: str) -> Job:
        """The job holding a file, or a folder (a combination's key) and its files."""
        if relative in self.owners:
            return self.owners[relative]
        if relative in self.keys:
            return self.job
        folder = relative.rstrip("/") + "/"
        if any(k.startswith(folder) for k in self.keys if k not in self.owners):
            return self.job
        return next((j for k, j in self.owners.items() if k.startswith(folder)), self.job)

    def prefix(self, relative: str) -> str:
        """results/ of the job that holds `relative` (a file or a combination's key)."""
        return f"{self.owner(relative).storage_prefix}results/"

    def key(self, relative: str) -> str:
        return self.prefix(relative) + relative

    def read_bytes(self, relative: str) -> bytes:
        if relative == SUMMARY and self.inherited:
            return (json.dumps(self.summary, indent=2, ensure_ascii=False) + "\n").encode()
        if relative == COMPARISON and self.comparison:
            rows = ["\t".join(COMPARISON_FIELDS)]
            rows += ["\t".join(r.get(f, "") for f in COMPARISON_FIELDS) for r in self.comparison]
            return ("\n".join(rows) + "\n").encode()
        if relative in self.inherited_tables:
            own = ""
            if relative in self.own:
                own = storage.read_bytes(self.key(relative)).decode()
            return join_tables(self.inherited_tables[relative], own).encode()
        return storage.read_bytes(self.key(relative))

    def read(self, relative: str) -> str:
        return self.read_bytes(relative).decode()

    def read_range(self, relative: str, start: int, length: int) -> bytes:
        """Part of one stored file (never a joined table)."""
        return storage.read_range(self.key(relative), start, length)

    def read_many(self, relatives: list[str]) -> dict[str, str]:
        """Several files at once, read concurrently (a page that reads one small file
        per combination waits for 285 round trips otherwise)."""
        with ThreadPoolExecutor(max_workers=READERS) as pool:
            return dict(zip(relatives, pool.map(self.read, relatives), strict=True))

    def sizes(self, under: str = "") -> dict[str, int]:
        """Bytes of each file (relative to results/), optionally only under a folder,
        from the manifest of the job holding it."""
        out = {p: n for p, n in manifest.sizes(self.job, under).items() if p not in self.owners}
        for holder in {j.id: j for j in self.owners.values()}.values():
            found = manifest.sizes(holder, under)
            out |= {p: n for p, n in found.items() if self._held_by(p, holder)}
        return out

    def _held_by(self, relative: str, holder: Job) -> bool:
        return relative in self.owners and self.owners[relative].id == holder.id

    @property
    def own(self) -> set[str]:
        """The files this job holds itself."""
        return {k for k in self.keys if k not in self.owners}


def inheritance(job: Job) -> tuple[Job | None, list[str]]:
    """The job whose alignments' results this one reads, and those alignments' keys."""
    of = job.params.get("inherits") or {}
    folders = [str(k) for k in of.get("alignments") or []]
    if not folders or job.extends_id is None or job.imported_at is not None:
        return None, []
    return job.extends, folders


def listing(job: Job) -> set[str]:
    """The job's own result files, from its manifest (never a storage listing)."""
    return manifest.paths(job)


def results(job: Job, keys: set[str] | None = None) -> Results:
    """The job's results, joined with what it inherits. `keys` is its own listing when
    the caller has one already."""
    own = listing(job) if keys is None else set(keys)
    out = Results(job, set(own))
    if job.status == Job.Status.SUCCEEDED and SUMMARY in own:
        out.summary = json.loads(storage.read_bytes(f"{job.storage_prefix}results/{SUMMARY}"))
    parent, folders = inheritance(job)
    if parent is None or job.status != Job.Status.SUCCEEDED:
        return out
    out.inherited = folders
    if parent.purged_at is not None or parent.deleted_at is not None:
        out.lost = True
        return out
    above = results(parent)
    if not above.summary:
        out.lost = True
        return out
    for key in above.keys:
        if key.split("/", 1)[0] in folders and key not in own:  # its own copy wins
            out.keys.add(key)
            out.owners[key] = above.owner(key)
    out.summary = join_summaries(above.summary, out.summary, set(folders))
    out.comparison = rerank(out.summary)
    for name in JOINED_TABLES:
        if name in above.keys:
            out.inherited_tables[name] = filter_rows(above.read(name), set(folders))
            out.keys.add(name)
    return out


def rerank(summary: dict[str, Any]) -> list[dict[str, str]]:
    """The joined combinations ranked as one job's SUMMARY ranks them (BIC within each
    alignment and data type), written back into the summary's entries; [] if a fit is
    missing (the joined table then stands as it is)."""
    combos = summary.get("combinations", [])
    fits = []
    for c in combos:
        fit = c.get("fit") or {}
        if any(f not in fit for f in FIT_FIELDS):
            return []
        values = {f: str(fit[f]) for f in FIT_FIELDS}
        fits.append({"aligner": c["aligner"], "tree": c["tree"], "model_key": c["model"], **values})
    rows = compare(fits)
    by_key = {f"{r['aligner']}/{r['tree']}/{r['model_key']}": r for r in rows}
    for c in combos:
        c["rank_bic"] = int(by_key[c["key"]]["rank_bic"])
        c["delta_bic"] = float(by_key[c["key"]]["delta_bic"])
    return rows


def cached(job: Job) -> Results:
    """`results(job)`, once per Job instance (one request): for readers that are handed
    only the job, file by file."""
    found = getattr(job, "_joined_results", None)
    if found is None:
        found = results(job)
        job._joined_results = found  # type: ignore[attr-defined]
    return found


def join_summaries(
    parent: dict[str, Any], own: dict[str, Any], folders: set[str]
) -> dict[str, Any]:
    """The parent's entries for the inherited alignments with the job's own, sorted by
    key as one job's summary is. The parent's primary combination stays the primary."""
    mine = {c["key"] for c in own.get("combinations", [])}
    aligned = {a["key"] for a in own.get("alignments", [])}
    alignments = [a for a in parent.get("alignments", []) if a["key"] in folders - aligned]
    combos = [
        dict(c) for c in parent.get("combinations", [])
        if c["aligner"] in folders and c["key"] not in mine
    ]  # fmt: skip
    return {
        **own,
        "alignments": sorted(alignments + own.get("alignments", []), key=lambda a: a["key"]),
        "combinations": sorted(
            combos + [dict(c) for c in own.get("combinations", [])], key=lambda c: c["key"]
        ),
        "primary": parent.get("primary") if combos else own.get("primary", ""),
    }


def filter_rows(text: str, folders: set[str]) -> str:
    """A shared table's header and the rows that belong to the inherited alignments:
    model_comparison.tsv's by its aligner column, ancestor_map.tsv's by its combination's,
    sources.tsv's alignment rows by key (its tree rows and every versions.tsv row hold
    for the whole job)."""
    lines = [line for line in text.splitlines() if line]
    if not lines:
        return ""
    header = lines[0].split("\t")
    if header[0] == "aligner":  # model_comparison.tsv
        return "\n".join([lines[0]] + [r for r in lines[1:] if r.split("\t")[0] in folders])
    if header[0] == "combination":  # ancestor_map.tsv
        rows = [r for r in lines[1:] if r.split("/", 1)[0] in folders]
        return "\n".join([lines[0], *rows])
    if header[:2] == ["key", "kind"]:  # sources.tsv
        rows = [
            r for r in lines[1:] if r.split("\t")[1] != "alignment" or r.split("\t")[0] in folders
        ]
        return "\n".join([lines[0], *rows])
    return "\n".join(lines)


def join_tables(inherited: str, own: str) -> str:
    """Two tables with one header: the inherited rows first, then the job's own, each
    row once (versions.tsv and sources.tsv repeat what both jobs ran)."""
    a = [line for line in inherited.splitlines() if line]
    b = [line for line in own.splitlines() if line]
    if not a:
        return own
    if not b:
        return "\n".join(a) + "\n"
    rows = list(dict.fromkeys(a[1:] + b[1:]))
    return "\n".join([a[0], *rows]) + "\n"


def dependents(job: Job) -> list[Job]:
    """Jobs that read some of their results from this one, directly or through another
    (not deleted): what deleting it would take from them."""
    out = []
    for child in job.extensions.filter(deleted_at__isnull=True).order_by("created_at"):
        if inheritance(child)[0] is not None:
            out += [child, *dependents(child)]
    return out
