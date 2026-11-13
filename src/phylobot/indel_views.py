"""The indel-length histograms: thumbnails on the job page, one chart in the alignment
viewer (server SVG; counting is `phylobot.indels`).

Every alignment is measured on its best tree by BIC, so a job's page and its viewer agree.
Charts carry colours as SVG attributes (dark.css remaps them; the viewer's chart downloads).
The axes are constants: thumbnails always run 1..THUMBNAIL_LAST, whichever job they show,
so any two can be compared; a row of them also shares one y scale.
"""

from __future__ import annotations

from collections.abc import Callable, Collection
from dataclasses import dataclass
from typing import Any

from django.core.cache import cache

from phylobot import catalog, indels
from phylobot.alignview import parse_fasta
from phylobot.models import Job
from phylobot.results import data_type

THUMBNAIL_LAST = 20
CHART_LAST = 40
CACHE_SECONDS = 3600


@dataclass(frozen=True)
class Measured:
    """One alignment's lengths and the tree they were read on."""

    aligner: str
    label: str
    lengths: indels.Lengths
    combination: str  # "" without a tree


@dataclass(frozen=True)
class Bar:
    x: float
    y: float
    width: float
    height: float
    title: str
    longer: bool  # the last bin: everything past the axis


@dataclass(frozen=True)
class Chart:
    measured: Measured
    width: int
    height: int
    left: int
    right: int
    top: int
    bottom: int
    bars: list[Bar]
    x_ticks: list[tuple[float, str]]
    y_ticks: list[tuple[float, str]]
    expected: str  # the geometric of the same mean, as SVG points; "" on thumbnails

    @property
    def center(self) -> float:
        return (self.left + self.right) / 2


def best_tree(summary: dict[str, Any], keys: Collection[str], aligner: str) -> tuple[str, str]:
    """The alignment's best combination by BIC (amino-acid models if it has them, whose
    ranks are their own) and its tree, outgroup-rooted if it can be; ("", "") if it has none."""
    combos = [
        c
        for c in summary.get("combinations", [])
        if c["aligner"] == aligner and f"{c['key']}/tree.nwk" in keys
    ]
    amino = [c for c in combos if data_type(c) == "amino"] or combos
    if not amino:
        return "", ""
    best = min(
        amino,
        key=lambda c: (c.get("rank_bic", 0), not catalog.is_modelfinder(c["model"]), c["key"]),
    )
    key = str(best["key"])
    rooted = f"{key}/tree.rooted.nwk"
    return key, rooted if rooted in keys else f"{key}/tree.nwk"


def measure(
    job: Job,
    aligner: str,
    label: str,
    alignment_path: str,
    summary: dict[str, Any],
    keys: Collection[str],
    read: Callable[[str], str],
) -> Measured:
    """Cached per job, alignment and tree: a finished job's results never change, and Fitch
    on every column takes 0.03-0.15 s on the paper's families."""
    combination, tree_path = best_tree(summary, keys, aligner)
    key = f"indels:v1:{job.id}:{aligner}:{combination}"
    found: indels.Lengths | None = cache.get(key)
    if found is None:
        rows = parse_fasta(read(alignment_path))
        found = _lengths(rows, read(tree_path) if tree_path else "")
        cache.set(key, found, CACHE_SECONDS)
    return Measured(aligner, label, found, combination if found.on_tree else "")


def _lengths(rows: dict[str, str], newick: str) -> indels.Lengths:
    if newick:
        try:
            measured = indels.on_tree(newick, rows)
            if measured.on_tree:
                return measured
        except KeyError:  # a tree of other names: the alignment still has its gaps
            pass
    return indels.in_alignment(rows)


def shares(measured: Measured, last: int) -> list[float]:
    total = measured.lengths.events
    return [n / total if total else 0.0 for n in measured.lengths.binned(last)]


def y_scale(tallest: float) -> float:
    """The next 20% above the tallest bar, so a row's bars read against round figures."""
    steps = int(tallest * 5 - 1e-9) + 1
    return max(min(steps, 5), 1) / 5


def chart(measured: Measured, last: int, scale: float, thumbnail: bool) -> Chart:
    width, height = (200, 104) if thumbnail else (760, 280)
    left, right = (30, width - 4) if thumbnail else (52, width - 12)
    top, bottom = (6, height - 18) if thumbnail else (30, height - 44)
    step = (right - left) / (last + 1)
    span = bottom - top
    gap = 1.0 if thumbnail else 2.0
    counts = measured.lengths.binned(last)
    bars = []
    for i, (n, share) in enumerate(zip(counts, shares(measured, last), strict=True)):
        if not n:
            continue
        name = f"Longer than {last}" if i == last else f"Length {i + 1}"
        bar_height = max(share / scale * span, 1.0)
        bars.append(Bar(
            x=round(left + i * step + gap / 2, 2),
            y=round(bottom - bar_height, 2),
            width=round(step - gap, 2),
            height=round(bar_height, 2),
            title=f"{name}: {n} ({100 * share:.0f}%)",
            longer=i == last,
        ))  # fmt: skip
    every = 5 if last <= 20 else 10
    numbers = [1, *range(every, last, every)]  # the last bin's label stands where `last` would
    x_ticks = [(round(left + (k - 0.5) * step, 2), str(k)) for k in numbers]
    x_ticks.append((round(left + (last + 0.5) * step, 2), f">{last}"))
    fractions = (1.0,) if thumbnail else (0.25, 0.5, 0.75, 1.0)
    y_ticks = [(round(bottom - f * span, 2), f"{100 * f * scale:.0f}%") for f in fractions]
    expected = ""
    if not thumbnail and measured.lengths.mean > 1:
        points = [
            f"{left + (i + 0.5) * step:.1f},{bottom - min(p / scale, 1.0) * span:.1f}"
            for i, p in enumerate(indels.geometric(measured.lengths.mean, last))
        ]
        expected = " ".join(points)
    return Chart(measured, width, height, left, right, top, bottom, bars, x_ticks, y_ticks,
                 expected)  # fmt: skip


def thumbnails(measured: list[Measured]) -> list[Chart]:
    """One per alignment on the same axes: 1..THUMBNAIL_LAST always, and one y scale."""
    tallest = max((max(shares(m, THUMBNAIL_LAST)) for m in measured), default=0.0)
    scale = y_scale(tallest)
    return [chart(m, THUMBNAIL_LAST, scale, thumbnail=True) for m in measured]


def large(measured: Measured) -> Chart:
    return chart(measured, CHART_LAST, y_scale(max(shares(measured, CHART_LAST))), thumbnail=False)


def tsv(measured: Measured) -> str:
    """The lengths behind a chart, every length rather than the chart's bins."""
    lengths = measured.lengths
    lines = [
        f"# {'on tree ' + measured.combination if lengths.on_tree else 'no tree: distinct gaps'}",
        f"# interior events {lengths.events}; end gaps left out {lengths.end_gaps}",
        "length\tevents",
    ]
    lines += [f"{k}\t{lengths.counts[k]}" for k in sorted(lengths.counts)]
    return "\n".join(lines) + "\n"
