"""A job's cells, frame and ancestors, run through the decomposition.

`run` reads a finished job's published files through `read` (a relative path to its text,
None when absent: the staged folder in the pipeline, a dict in the tests) and returns the
files to publish under `uncertainty/`. `phylobot.uncertainty` does the arithmetic; this
module decides what the cells are and what each one says at each reference column.

Pipeline code: it runs in phylobot-tools, so it imports nothing from Django.
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from phylobot import alignment_map, analogs, ancestors, results
from phylobot import uncertainty as u
from phylobot.gaps import MASKING
from phylobot.trees import load_tree

Read = Callable[[str], str | None]
POSTERIORS = "ancestors_posteriors.tsv"

OUT = "uncertainty"
MODES = ("posterior", "map")
MODE_LABELS = {"posterior": "Posterior", "map": "Most probable sequences"}
PARTS = ("all", "presence", "residue")
FACTOR_LABELS = {
    "alignment": "Alignment",
    "tree": "Tree",
    "model": "Model",
    "gaps": "Gap placement",
}
WEIGHTINGS = ("flat", "bic")
# Gap versions a run may use (§2): the masks by default; the rest only when asked for.
OPT_IN = (ancestors.ORIGINAL, "prank", "grasp")
VERSIONS = (*MASKING, *OPT_IN)
POINT_ESTIMATES = ("prank", "grasp")  # their own residues, no posteriors
DRIVER_FLOOR = 0.005  # sites varying less than this between cells get no drivers row
FORMAT = 1


class SettingsError(ValueError):
    pass


@dataclass(frozen=True)
class Settings:
    """What a run was asked for (§2.1): the weighting, multipliers per factor level, the
    gap versions (empty: the masks each combination has) and the reference combination
    (empty: the viewers' default)."""

    weighting: str = "flat"
    multipliers: dict[str, dict[str, float]] = field(default_factory=dict)
    versions: tuple[str, ...] = ()
    reference: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "weighting": self.weighting,
            "multipliers": {
                f: dict(sorted(m.items())) for f, m in sorted(self.multipliers.items())
            },
            "versions": list(self.versions),
            "reference": self.reference,
        }

    @classmethod
    def parse(cls, data: object) -> Settings:
        """Settings from JSON, refusing anything but the shapes and names it may hold:
        these come from a job's params, so they are checked again here."""
        if not isinstance(data, dict) or set(data) - {
            "weighting",
            "multipliers",
            "versions",
            "reference",
        }:
            raise SettingsError("settings must be an object with known keys")
        weighting = data.get("weighting", "flat")
        if weighting not in WEIGHTINGS:
            raise SettingsError(f"unknown weighting {weighting!r}")
        multipliers: dict[str, dict[str, float]] = {}
        raw = data.get("multipliers") or {}
        if not isinstance(raw, dict):
            raise SettingsError("multipliers must be an object")
        for factor, levels in raw.items():
            if factor not in u.FACTORS or not isinstance(levels, dict):
                raise SettingsError(f"unknown factor {factor!r}")
            multipliers[factor] = {}
            for level, value in levels.items():
                if (
                    not isinstance(level, str)
                    or not isinstance(value, int | float)
                    or isinstance(value, bool)
                ):
                    raise SettingsError("a multiplier is a level name and a number")
                if not (0 <= value <= 1000 and math.isfinite(value)):
                    raise SettingsError("multipliers run from 0 to 1000")
                multipliers[factor][level] = float(value)
        versions = data.get("versions") or []
        if not isinstance(versions, list) or any(v not in VERSIONS for v in versions):
            raise SettingsError("unknown gap version")
        reference = data.get("reference") or ""
        if not isinstance(reference, str):
            raise SettingsError("reference must be a combination key")
        return cls(weighting, multipliers, tuple(dict.fromkeys(versions)), reference)


@dataclass(frozen=True)
class Cell:
    combination: dict[str, Any]
    version: str
    weight: float

    @property
    def key(self) -> str:
        return f"{self.combination['key']}@{self.version}"

    @property
    def levels(self) -> tuple[str, str, str, str]:
        c = self.combination
        return (c["aligner"], c["tree"], c["model"], self.version)


def version_files(version: str) -> str:
    return ancestors.VERSIONS[version][1]


def eligible(summary: dict[str, Any], read: Read) -> list[dict[str, Any]]:
    """The combinations with ancestors and posteriors, in the viewers' order, without
    fixed-model rows that repeat ModelFinder's pick (§2: the same fit twice)."""
    order = results.combination_order(summary.get("primary", ""))
    found = [
        c for c in summary.get("combinations", [])
        if read(f"{c['key']}/ancestors.fasta") is not None and has_posteriors(read, c["key"])
    ]  # fmt: skip
    kept: list[dict[str, Any]] = []
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for c in found:
        groups.setdefault((c["aligner"], c.get("data", "amino")), []).append(c)
    for group in groups.values():
        kept += results.merge_modelfinder_duplicates(group, {})[0]
    return sorted(kept, key=order)


def has_posteriors(read: Read, key: str) -> bool:
    return read(f"{key}/ancestors_posteriors.tsv") is not None


def default_reference(combinations: list[dict[str, Any]]) -> str:
    """The viewers' default: the first alignment's best combination by BIC."""
    if not combinations:
        return ""
    first = [c for c in combinations if c["aligner"] == combinations[0]["aligner"]]
    ranked = [c for c in first if isinstance(c.get("rank_bic"), int)]
    return str(min(ranked, key=lambda c: c["rank_bic"])["key"] if ranked else first[0]["key"])


def bic_of(combination: dict[str, Any]) -> float | None:
    try:
        value = float(combination.get("fit", {}).get("bic", ""))
    except TypeError, ValueError:
        return None
    return value if math.isfinite(value) else None


def base_weights(combinations: list[dict[str, Any]], settings: Settings) -> dict[str, float]:
    """{combination: weight} before gap versions share it (§2.1). BIC: Schwarz weights
    within each alignment and data type, each of which then weighs the same; a group
    with a combination lacking a BIC falls back to flat within it. Multipliers apply on
    top, a zero leaving the level out."""
    weights = {c["key"]: 1.0 for c in combinations}
    if settings.weighting == "bic":
        groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
        for c in combinations:
            groups.setdefault((c["aligner"], c.get("data", "amino")), []).append(c)
        for group in groups.values():
            bics = [bic_of(c) for c in group]
            if any(b is None for b in bics):
                share = [1.0] * len(group)
            else:
                best = min(b for b in bics if b is not None)
                share = [math.exp(-((b or 0.0) - best) / 2) for b in bics]
            total = sum(share)
            for c, s in zip(group, share, strict=True):
                weights[c["key"]] = len(combinations) / len(groups) * s / total
    for c in combinations:
        for factor, level in zip(u.FACTORS[:3], (c["aligner"], c["tree"], c["model"]), strict=True):
            weights[c["key"]] *= settings.multipliers.get(factor, {}).get(level, 1.0)
    return weights


def versions_of(read: Read, key: str, settings: Settings) -> list[str]:
    """The gap versions a combination contributes: those asked for that it published,
    or (none asked for) its masks, or with no masks at all its original reconstruction."""
    have = [v for v in VERSIONS if read(f"{key}/{version_files(v)}") is not None]
    if settings.versions:
        chosen = [v for v in settings.versions if v in have]
    else:
        chosen = [v for v in MASKING if v in have] or [ancestors.ORIGINAL]
    gaps = settings.multipliers.get("gaps", {})
    return [v for v in chosen if gaps.get(v, 1.0) > 0]


def cells_of(combinations: list[dict[str, Any]], read: Read, settings: Settings) -> list[Cell]:
    weights = base_weights(combinations, settings)
    cells = []
    for c in combinations:
        if weights[c["key"]] <= 0:
            continue
        versions = versions_of(read, c["key"], settings)
        gaps = settings.multipliers.get("gaps", {})
        for version in versions:
            share = weights[c["key"]] / len(versions) * gaps.get(version, 1.0)
            cells.append(Cell(c, version, share))
    return cells


# --- Frame -------------------------------------------------------------------------


@dataclass(frozen=True)
class Frame:
    """Reference columns mapped into one alignment: `source[r]` is the column of that
    alignment holding most of reference column r's residues (-1: none), `agreement[r]`
    the Jaccard index of the two columns' residues."""

    source: np.ndarray
    agreement: np.ndarray


def frame(reference: dict[str, str], other: dict[str, str]) -> Frame:
    width = len(next(iter(reference.values()), ""))
    if reference == other:
        return Frame(np.arange(width), np.ones(width))
    mapped = alignment_map.map_columns(reference, other)
    source = np.full(width, -1)
    agreement = np.zeros(width)
    for column, target in mapped.target.items():
        source[column] = target
        agreement[column] = mapped.agreement[column]
    return Frame(source, agreement)


# --- What a cell says at each column ------------------------------------------------


def posterior_table(text: str, node: str) -> np.ndarray:
    """[sites, 20] `node`'s posteriors in AMINO order, each row summing to 1 (the TSV
    rounds to five places)."""
    header, *rows = ancestors.posterior_lines(text, node)
    states = header.split("\t")[2:]
    order = [states.index(a) if a in states else None for a in u.AMINO]  # absent: never
    by_site = {}
    for row in rows:
        fields = row.split("\t")
        values = [float(x) for x in fields[2:]]
        by_site[int(fields[1])] = [0.0 if i is None else values[i] for i in order]
    table = np.array([by_site[s] for s in sorted(by_site)], dtype=float).reshape(-1, 20)
    total = table.sum(axis=1, keepdims=True)
    return np.asarray(np.divide(table, total, out=np.full_like(table, 1 / 20), where=total > 0))


class Combination:
    """One combination's files, read once and parsed per node on demand."""

    def __init__(self, info: dict[str, Any], read: Read) -> None:
        self.info = info
        self.key = info["key"]
        self.read = read
        self._fastas: dict[str, dict[str, str]] = {}
        self._ml: dict[str, np.ndarray] | None = None
        index = read(f"{self.key}/{ancestors.index_name(POSTERIORS)}")
        self._index = ancestors.read_index(index) if index else None
        self._whole: str | None = None
        self._bytes = b""  # the table, for slicing by the index

    def posteriors(self, node: str) -> np.ndarray | None:
        if self._whole is None:
            self._whole = self.read(f"{self.key}/{POSTERIORS}") or ""
            self._bytes = self._whole.encode() if self._index else b""
        text: str | None = self._whole
        if self._index is not None:
            text = self._index.text(node, lambda start, n: self._bytes[start : start + n])
        if not text:
            return None
        table = posterior_table(text, node)
        return table if len(table) else None

    def fasta(self, version: str) -> dict[str, str]:
        if version not in self._fastas:
            text = self.read(f"{self.key}/{version_files(version)}") or ""
            self._fastas[version] = ancestors.parse_fasta(text)
        return self._fastas[version]

    def ml_gaps(self, node: str) -> np.ndarray | None:
        if self._ml is None:
            self._ml = {}
            text = self.read(f"{self.key}/gaps-ml.tsv") or ""
            rows: dict[str, dict[int, float]] = {}
            for line in text.splitlines()[1:]:
                name, site, p = line.split("\t")
                rows.setdefault(name, {})[int(site)] = float(p)
            self._ml = {n: np.array([s[k] for k in sorted(s)]) for n, s in rows.items()}
        return self._ml.get(node)


def cell_states(
    combo: Combination, version: str, node: str
) -> tuple[np.ndarray, np.ndarray] | None:
    """(P(gap) [sites], residues [sites, 20]) on the combination's own columns, or None
    when it has no reconstruction of the node in this version."""
    posterior = combo.posteriors(node)
    if posterior is None:
        return None
    width = len(posterior)
    gap = np.zeros(width)
    if version in ("parsimony", "fitch", *POINT_ESTIMATES):
        row = combo.fasta(version).get(node)
        if row is None or len(row) != width:
            return None
        gap = np.array([c == u.GAP for c in row], dtype=float)
        if version in POINT_ESTIMATES:  # their own residues: point masses
            posterior = np.full((width, 20), 1 / 20)
            for i, c in enumerate(row):
                if c in u.AMINO:
                    posterior[i] = 0.0
                    posterior[i, u.AMINO.index(c)] = 1.0
    elif version == "ml":
        p = combo.ml_gaps(node)
        if p is None or len(p) != width:
            return None
        gap = p
    return gap, posterior


# --- The run -----------------------------------------------------------------------


def preorder(newick: str) -> list[str]:
    tree = load_tree(newick)
    return [n.label for n in tree.preorder_internal_node_iter() if n.label]


@dataclass
class NodeResult:
    node: str
    present: int
    of: int
    leaves: int
    modes: dict[str, u.Split]
    pairs: dict[str, np.ndarray]  # mode -> [cells, cells] (NaN where either lacks the node)
    drivers: list[tuple[int, str, str, str, float]]  # column, factor, level, state, p
    agreement: np.ndarray


def tsv(rows: Iterator[list[Any]] | list[list[Any]]) -> str:
    return "".join("\t".join(str(x) for x in row) + "\n" for row in rows)


def number(value: float, places: int = 5) -> str:
    """Fixed places with trailing zeros dropped: most values in a run are 0, written "0"."""
    text = f"{value:.{places}f}".rstrip("0").rstrip(".")
    return "0" if text in ("", "-0") else text


SITE_VALUES = ("total", "within", *u.FACTORS, "other")


def site_values(part: u.Components) -> list[np.ndarray]:
    """A part's per-site arrays in SITE_VALUES order."""
    fixed = {"total": part.total, "within": part.within, "other": part.other}
    return [fixed[k] if k in fixed else part.shares[k] for k in SITE_VALUES]


# sites/<n>.tsv: one row per mode and reference column, the presence and the residue
# part side by side (the whole is their sum, so it is not written).
SITE_HEADER = [
    "mode", "column", "agreement",
    *(f"{part}_{k}" for part in ("presence", "residue") for k in SITE_VALUES),
]  # fmt: skip


@dataclass
class Prepared:
    summary: dict[str, Any]
    combinations: list[dict[str, Any]]
    cells: list[Cell]
    reference: str
    by_clade: bool
    places: dict[str, dict[str, analogs.Place]]
    reference_alignment: dict[str, str]
    frames: dict[str, Frame]  # combination -> its alignment's frame
    order: list[str]


def prepare(read: Read, settings: Settings) -> Prepared:
    summary = json.loads(read("summary.json") or "{}")
    combinations = eligible(summary, read)
    if not combinations:
        raise SettingsError("the job has no ancestors with posteriors")
    keys = [c["key"] for c in combinations]
    reference = settings.reference or default_reference(combinations)
    if reference not in keys:
        raise SettingsError(f"reference {reference!r} is not a combination with ancestors")
    # By clade when every combination has an outgroup-rooted tree (analog_views.trees).
    by_clade = all(read(f"{k}/tree.rooted.nwk") is not None for k in keys)
    tree_file = "tree.rooted.nwk" if by_clade else "tree.nwk"
    places = {k: analogs.places(read(f"{k}/{tree_file}") or "") for k in keys}
    alignments: dict[str, dict[str, str]] = {}
    for k in keys:
        alignments[k] = ancestors.parse_fasta(read(f"{k}/alignment.fasta") or "")
    reference_alignment = alignments[reference]
    frames: dict[str, Frame] = {}
    by_text: dict[str, Frame] = {}
    for k in keys:
        text = json.dumps(alignments[k], sort_keys=True)
        if text not in by_text:
            by_text[text] = frame(reference_alignment, alignments[k])
        frames[k] = by_text[text]
    cells = [c for c in cells_of(combinations, read, settings) if c.weight > 0]
    if not cells:
        raise SettingsError("every combination was weighted out")
    order = [n for n in preorder(read(f"{reference}/{tree_file}") or "") if n in places[reference]]
    return Prepared(
        summary,
        combinations,
        cells,
        reference,
        by_clade,
        places,
        reference_alignment,
        frames,
        order,
    )


def node_result(
    p: Prepared, read_combo: dict[str, Combination], node: str, matches: dict[str, str]
) -> NodeResult:
    width = len(next(iter(p.reference_alignment.values()), ""))
    count = len(p.cells)
    gap = np.zeros((count, width))
    residues = np.full((count, width, 20), 1 / 20)
    weights = np.zeros((count, width))
    for i, cell in enumerate(p.cells):
        key = cell.combination["key"]
        name = matches.get(key)
        if name is None:
            continue
        states = cell_states(read_combo[key], cell.version, name)
        if states is None:
            continue
        g, q = states
        source = p.frames[key].source
        ok = (source >= 0) & (source < len(g))
        gap[i, ok] = g[source[ok]]
        residues[i, ok] = q[source[ok]]
        weights[i, ok] = cell.weight
    levels = [c.levels for c in p.cells]
    point_gap, point = u.most_probable(gap, residues)
    modes = {
        "posterior": u.categorical(gap, residues, weights, levels),
        "map": u.categorical(point_gap, point, weights, levels),
    }
    valid = weights > 0
    pairs = {}
    for mode, (g, q) in {"posterior": (gap, residues), "map": (point_gap, point)}.items():
        present = 1.0 - g
        same = g[:, None] * g[None] + present[:, None] * present[None] * np.einsum(
            "csk,dsk->cds", q, q
        )
        both = valid[:, None] & valid[None]
        found = np.where(both, 1.0 - same, 0.0).sum(-1)
        lacking = ~both.any(-1)
        found[lacking] = np.nan
        pairs[mode] = found
    drivers: list[tuple[int, str, str, str, float]] = []
    split = modes["posterior"].all
    between = split.total - split.within
    varying = np.flatnonzero(between >= DRIVER_FLOOR)
    for f, factor in enumerate(u.FACTORS):
        if len({lv[f] for lv, w in zip(levels, weights.any(1), strict=True) if w}) < 2:
            continue
        for level, (state, prob) in u.level_states(gap, residues, weights, levels, f).items():
            drivers += [
                (int(s) + 1, factor, level, u.STATES[state[s]], float(prob[s])) for s in varying
            ]
    combos_present = {
        c.combination["key"] for c, w in zip(p.cells, weights.any(1), strict=True) if w
    }
    agreement = np.ones(width)
    for key in combos_present:
        agreement = np.minimum(agreement, p.frames[key].agreement)
    return NodeResult(
        node,
        present=len(combos_present),
        of=len({c.combination["key"] for c in p.cells}),
        leaves=len(p.places[p.reference][node].below),
        modes=modes,
        pairs=pairs,
        drivers=sorted(drivers),
        agreement=agreement,
    )


def run(read: Read, settings: Settings) -> dict[str, str]:
    """{path under uncertainty/: text} for every file of §4."""
    p = prepare(read, settings)
    combos = {c["key"]: Combination(c, read) for c in p.combinations}
    matched = analogs.index(p.places, p.reference, p.by_clade)
    node_rows: list[list[Any]] = [
        ["mode", "node", "part", "present", "of", "leaves", "total", "within", *u.FACTORS, "other"]
    ]
    files: dict[str, str] = {}
    index_rows = [["node", "file"]]
    totals: dict[str, dict[str, list[dict[str, float]]]] = {
        m: {part: [] for part in PARTS} for m in MODES
    }
    pair_sums = {m: np.zeros((len(p.cells), len(p.cells))) for m in MODES}
    pair_counts = {m: np.zeros((len(p.cells), len(p.cells))) for m in MODES}
    for number_, node in enumerate(p.order, start=1):
        matches = {p.reference: node} | {a.combination: a.node for a in matched.get(node, [])}
        result = node_result(p, combos, node, matches)
        name = f"{number_}.tsv"
        index_rows.append([node, name])
        site_rows: list[list[Any]] = [SITE_HEADER]
        for mode in MODES:
            for part, comp in result.modes[mode].parts().items():
                sums = comp.summed()
                totals[mode][part].append(sums)
                node_rows.append([
                    mode, node, part, result.present, result.of, result.leaves,
                    *(number(sums[k], 4) for k in ("total", "within", *u.FACTORS, "other")),
                ])  # fmt: skip
            split = result.modes[mode]
            columns = [site_values(c) for c in (split.presence, split.residue)]
            for i in range(len(result.agreement)):
                values = [number(float(v[i]), 4) for part in columns for v in part]
                site_rows.append([mode, i + 1, number(float(result.agreement[i]), 3), *values])
            found = result.pairs[mode]
            seen = ~np.isnan(found)
            pair_sums[mode][seen] += found[seen]
            pair_counts[mode][seen] += 1
        files[f"sites/{name}"] = tsv(site_rows)
        files[f"drivers/{name}"] = tsv(
            [["column", "factor", "level", "state", "p"]]
            + [[c, f, lv, st, number(pr, 3)] for c, f, lv, st, pr in result.drivers]
        )
        # Cells by their position in summary.json's "cells".
        pair_rows: list[list[Any]] = [["mode", "a", "b", "differences"]]
        for mode in MODES:
            found = result.pairs[mode]
            for i in range(len(p.cells)):
                for j in range(len(p.cells)):
                    if not np.isnan(found[i, j]):
                        pair_rows.append([mode, i, j, number(float(found[i, j]), 3)])
        files[f"pairs/{name}"] = tsv(pair_rows)
    for folder in ("sites", "drivers", "pairs"):
        files[f"{folder}/{ancestors.SPLIT_INDEX}"] = tsv(index_rows)
    files["nodes.tsv"] = tsv(node_rows)
    files["summary.json"] = (
        json.dumps(summary_of(p, totals, pair_sums, pair_counts), indent=1) + "\n"
    )
    return files


def summary_of(
    p: Prepared,
    totals: dict[str, dict[str, list[dict[str, float]]]],
    pair_sums: dict[str, np.ndarray],
    pair_counts: dict[str, np.ndarray],
) -> dict[str, Any]:
    keys = ("total", "within", *u.FACTORS, "other")
    weight = sum(c.weight for c in p.cells)
    factors = []
    for f, factor in enumerate(u.FACTORS):
        levels: dict[str, dict[str, Any]] = {}
        for cell in p.cells:
            level = cell.levels[f]
            entry = levels.setdefault(
                level, {"key": level, "label": level_label(cell, factor), "weight": 0.0, "cells": 0}
            )
            entry["weight"] += cell.weight / weight
            entry["cells"] += 1
        factors.append({
            "key": factor, "label": FACTOR_LABELS[factor], "varied": len(levels) > 1,
            "levels": [dict(e, weight=round(e["weight"], 6)) for e in levels.values()],
        })  # fmt: skip
    job: dict[str, Any] = {}
    for mode in MODES:
        job[mode] = {}
        for part in PARTS:
            rows = totals[mode][part]
            job[mode][part] = {
                k: round(float(np.mean([r[k] for r in rows])) if rows else 0.0, 4) for k in keys
            }
    pairs = {}
    for mode in MODES:
        mean = np.divide(
            pair_sums[mode],
            pair_counts[mode],
            out=np.full_like(pair_sums[mode], np.nan),
            where=pair_counts[mode] > 0,
        )
        pairs[mode] = [[None if np.isnan(x) else round(float(x), 3) for x in row] for row in mean]
    return {
        "format": FORMAT,
        "reference": p.reference,
        "by_clade": p.by_clade,
        "columns": len(next(iter(p.reference_alignment.values()), "")),
        "nodes": len(p.order),
        "combinations": len({c.combination["key"] for c in p.cells}),
        "factors": factors,
        "cells": [
            {
                "key": c.key,
                "combination": c.combination["key"],
                "version": c.version,
                "levels": list(c.levels),
                "weight": round(c.weight / weight, 6),
                "label": c.combination.get("label", c.combination["key"]),
            }
            for c in p.cells
        ],  # fmt: skip
        "job": job,
        "pairs": pairs,
        "point_estimates": sorted({c.version for c in p.cells} & set(POINT_ESTIMATES)),
    }


def level_label(cell: Cell, factor: str) -> str:
    c = cell.combination
    match factor:
        case "alignment":
            return str(c.get("aligner_label") or c["aligner"])
        case "tree":
            return str(c.get("tree_label") or c["tree"])
        case "model":  # ModelFinder's pick can differ by alignment: the level is the choice
            return str(c.get("model_label") or c["model"])
    return ancestors.short_version_label(cell.version)
