"""How a combination's property measurements are stored.

PROPERTIES writes, per combination, under `<combination>/properties/`:

    draws.parquet    every draw: one row per (gap version, node, draw), a float64 column
                     per property; sorted, zstd, row groups of NODES_PER_GROUP nodes
    summaries.json   {"nodes": {version: {node: summary}}}: each node's summary as the
                     panel reads it (n, mean, sd, interval, f(ML), shift, cv, ...)
    index.json       what was asked and measured

Three files per combination, where a file per node, gap version and property made a
285-combination job 488,205 objects.
The per-draw values are kept, not only the summaries, so a histogram can be rebinned and
a figure rebuilt from what was published.

The readers fetch one node's draws with ranged reads: Parquet's footer, then only the
row groups whose statistics can hold that node. pyarrow is imported only where draws are
read or written, never by importing this module.
"""

from __future__ import annotations

import io
import json
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

FOLDER = "properties"
DRAWS = "draws.parquet"
SUMMARIES = "summaries.json"
INDEX = "index.json"
# One-node row groups double the file (tiny groups compress badly: the storage report);
# a group of several nodes is still one small read for the panel.
NODES_PER_GROUP = 8
COMPRESSION = "zstd"


def draws_path(combination: str) -> str:
    """Relative to results/."""
    return f"{combination}/{FOLDER}/{DRAWS}"


def summaries_path(combination: str) -> str:
    return f"{combination}/{FOLDER}/{SUMMARIES}"


@dataclass
class Measured:
    """One node and gap version's measurements, as PROPERTIES makes them."""

    version: str
    node: str
    values: dict[str, list[float]]  # property key -> draws, in draw order
    summary: dict[str, Any]


def write(out: Path, measured: Iterable[Measured], keys: list[str], draws: int) -> None:
    """draws.parquet and summaries.json under `out` (the combination's properties/)."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    out.mkdir(parents=True, exist_ok=True)
    found = sorted(measured, key=lambda m: (m.version, m.node))
    columns: dict[str, list[Any]] = {"version": [], "node": [], "draw": []}
    columns |= {key: [] for key in keys}
    nodes: dict[str, dict[str, Any]] = {}
    for m in found:
        n = max((len(v) for v in m.values.values()), default=0)
        columns["version"] += [m.version] * n
        columns["node"] += [m.node] * n
        columns["draw"] += range(1, n + 1)
        for key in keys:
            values = m.values.get(key, [])
            columns[key] += values + [None] * (n - len(values))  # a property not drawn: null
        nodes.setdefault(m.version, {})[m.node] = m.summary
    schema = pa.schema(
        [("version", pa.string()), ("node", pa.string()), ("draw", pa.int32())]
        + [(key, pa.float64()) for key in keys]
    )
    table = pa.table(columns, schema=schema)
    pq.write_table(
        table,
        out / DRAWS,
        compression=COMPRESSION,
        row_group_size=max(1, draws) * NODES_PER_GROUP,
        write_statistics=True,
    )
    (out / SUMMARIES).write_text(json.dumps({"nodes": nodes}, indent=1, sort_keys=True) + "\n")


def summary(text: str, node: str, version: str) -> dict[str, Any] | None:
    """One node's summary from summaries.json, None if it was not measured there."""
    try:
        found = json.loads(text)["nodes"][version][node]
    except KeyError, TypeError, ValueError:
        return None
    return found if isinstance(found, dict) and found.get("properties") else None


def measured_nodes(text: str, version: str) -> set[str]:
    """The nodes summaries.json has for a gap version."""
    try:
        found = json.loads(text)["nodes"][version]
    except KeyError, TypeError, ValueError:
        return set()
    return set(found) if isinstance(found, dict) else set()


class RangedFile(io.RawIOBase):
    """A read-only, seekable file over ranged reads of a stored object, for pyarrow:
    reading one node touches the footer and that node's row groups, not the file."""

    def __init__(self, size: int, read_range: Callable[[int, int], bytes]) -> None:
        self.size = size
        self.position = 0
        self.read_range = read_range

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.position

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        base = {io.SEEK_SET: 0, io.SEEK_CUR: self.position, io.SEEK_END: self.size}[whence]
        self.position = max(0, base + offset)
        return self.position

    def readinto(self, buffer: Any) -> int:
        n = min(len(buffer), self.size - self.position)
        if n <= 0:
            return 0
        data = self.read_range(self.position, n)
        buffer[: len(data)] = data
        self.position += len(data)
        return len(data)


def read_draws(source: Any, node: str, version: str, keys: list[str]) -> dict[str, list[float]]:
    """One node's draws per property from draws.parquet (a path or a file object).
    Properties the file lacks are left out (their rows are then drawn without bars)."""
    import pyarrow.parquet as pq

    file = pq.ParquetFile(source)
    names = set(file.schema_arrow.names)
    wanted = [k for k in keys if k in names]
    if not wanted:
        return {}
    groups = []
    for i in range(file.metadata.num_row_groups):
        group = file.metadata.row_group(i)
        if _may_hold(group, "version", version) and _may_hold(group, "node", node):
            groups.append(i)
    if not groups:
        return {}
    table = file.read_row_groups(groups, columns=["version", "node", "draw", *wanted])
    rows = table.to_pydict()
    picked = [
        i for i, (v, n) in enumerate(zip(rows["version"], rows["node"], strict=True))
        if v == version and n == node
    ]  # fmt: skip
    picked.sort(key=lambda i: rows["draw"][i])
    out = {k: [float(rows[k][i]) for i in picked if rows[k][i] is not None] for k in wanted}
    return {k: values for k, values in out.items() if values}


def _may_hold(group: Any, column: str, value: str) -> bool:
    """Whether a row group's statistics allow `value` in `column` (True if unknown)."""
    for i in range(group.num_columns):
        chunk = group.column(i)
        if chunk.path_in_schema != column:
            continue
        stats = chunk.statistics
        if stats is None or not stats.has_min_max:
            return True
        return bool(stats.min <= value <= stats.max)
    return True
