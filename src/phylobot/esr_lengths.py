"""ESR by branch length: the job page's "By branch length" section.

Accuracy falls on long pendant branches by construction (a distant tip is harder to
reconstruct), so the question here is whether *calibration* changes with them: the gap
between the mean posterior and accuracy, and the expected calibration error (ECE), per
quarter of the held-out sequences by pendant length; a scatter of each sequence's gap
against its pendant; and, with both levels, E1's log score minus E0's per sequence (where
the sequence informed its own reconstruction through the fit, E0 is too kind to it).

Two practical figures sit beside them: wrong calls with a posterior of 0.9 or more, per
sequence, and how many wrong calls AltAll (every residue with a posterior of at least
0.2) would have caught, i.e. where the true residue had at least 0.2.

Everything comes from the job's esr_sites.tsv and esr_tips.tsv at view time, so finished
jobs get it too (their results never change). The sites file is several MB; it is tallied
once per process into per-sequence counts and kept (`tally`).
"""

import functools
import math
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

from phylobot import esr, storage

BINS = 10  # as the pipeline's reliability bins (pipeline.esr.BINS)
CONFIDENT = 0.9  # "confidently wrong": the ML residue had at least this posterior
ALTALL = 0.2  # AltAll's usual threshold (Eick et al. 2017)
MIN_QUARTER_TIPS = 8  # fewer held-out sequences than this: no quarters (two per quarter)
MIN_SCATTER_TIPS = 3
TALLIES_CACHED = 64
SITE_COLUMNS = ("tip", "level", "true", "ml", "pp_ml", "pp_true")
QUARTERS = ("Shortest quarter", "Second quarter", "Third quarter", "Longest quarter")


@dataclass
class Tally:
    """One sequence at one level, or a group of them: counts to sum, figures to read."""

    tips: int = 0
    sites: int = 0
    correct: int = 0
    pp: float = 0.0  # summed posterior of the ML residue
    log: float = 0.0  # summed log posterior of the true residue
    confident_wrong: int = 0
    caught: int = 0  # wrong calls whose true residue had at least ALTALL
    bin_sites: list[int] = field(default_factory=lambda: [0] * BINS)
    bin_correct: list[int] = field(default_factory=lambda: [0] * BINS)
    bin_pp: list[float] = field(default_factory=lambda: [0.0] * BINS)

    def add(self, other: Tally) -> None:
        self.tips += other.tips
        self.sites += other.sites
        self.correct += other.correct
        self.pp += other.pp
        self.log += other.log
        self.confident_wrong += other.confident_wrong
        self.caught += other.caught
        for i in range(BINS):
            self.bin_sites[i] += other.bin_sites[i]
            self.bin_correct[i] += other.bin_correct[i]
            self.bin_pp[i] += other.bin_pp[i]

    @property
    def wrong(self) -> int:
        return self.sites - self.correct

    @property
    def accuracy(self) -> float:
        return self.correct / self.sites if self.sites else math.nan

    @property
    def mean_pp(self) -> float:
        return self.pp / self.sites if self.sites else math.nan

    @property
    def gap(self) -> float:
        return self.mean_pp - self.accuracy

    @property
    def log_score(self) -> float:
        return self.log / self.sites if self.sites else math.nan

    @property
    def ece(self) -> float:
        """Expected calibration error: |accuracy - mean posterior| per bin of the ML
        residue's posterior, weighted by the bin's sites."""
        if not self.sites:
            return math.nan
        return sum(abs(self.bin_correct[i] - self.bin_pp[i]) for i in range(BINS)) / self.sites

    @property
    def confident_wrong_per_tip(self) -> float:
        return self.confident_wrong / self.tips if self.tips else math.nan

    @property
    def caught_share(self) -> float:
        return self.caught / self.wrong if self.wrong else math.nan


Tallies = dict[tuple[str, str], Tally]  # (tip, level) -> its tally


def tally_sites(lines: Iterable[str]) -> Tallies:
    """esr_sites.tsv's rows summed per sequence and level. Reads only the columns it needs;
    a file without them gives nothing."""
    out: Tallies = {}
    rows = iter(lines)
    header = next(rows, "").rstrip("\n").split("\t")
    if not all(c in header for c in SITE_COLUMNS):
        return {}
    tip_i, level_i, true_i, ml_i, pp_i, true_pp_i = (header.index(c) for c in SITE_COLUMNS)
    for line in rows:
        cells = line.rstrip("\n").split("\t")
        if len(cells) < len(header):
            continue
        key = (cells[tip_i], cells[level_i])
        t = out.get(key)
        if t is None:
            t = out[key] = Tally(tips=1)
        p, p_true = float(cells[pp_i]), float(cells[true_pp_i])
        hit = cells[true_i] == cells[ml_i]
        b = min(int(p * BINS), BINS - 1)
        t.sites += 1
        t.correct += hit
        t.pp += p
        t.log += math.log(max(p_true, esr.LOG_FLOOR))
        t.bin_sites[b] += 1
        t.bin_correct[b] += hit
        t.bin_pp[b] += p
        if not hit:
            t.confident_wrong += p >= CONFIDENT
            t.caught += p_true >= ALTALL
    return out


@functools.lru_cache(maxsize=TALLIES_CACHED)
def tally(key: str) -> Tallies:
    """A finished job's esr_sites.tsv, tallied once (its results never change)."""
    return tally_sites(storage.read_bytes(key).decode().splitlines())


def pendants(tips_text: str) -> dict[str, float]:
    """Each held-out sequence's pendant branch length, from esr_tips.tsv."""
    lines = tips_text.splitlines()
    if not lines:
        return {}
    header = lines[0].split("\t")
    if "tip" not in header or "pendant" not in header:
        return {}
    tip_i, pendant_i = header.index("tip"), header.index("pendant")
    out = {}
    for line in lines[1:]:
        cells = line.split("\t")
        try:
            out[cells[tip_i]] = float(cells[pendant_i])
        except IndexError, ValueError:
            continue
    return out


# --- the table ------------------------------------------------------------------------


@dataclass
class Row:
    level: str
    label: str
    low: float  # the group's shortest and longest pendant
    high: float
    median: float
    t: Tally


def quarters(tallies: Tallies, lengths: dict[str, float], levels: list[str]) -> list[Row]:
    """Per level, the held-out sequences in four groups by pendant length (as equal in
    number as they divide), then all of them. Without enough sequences, only "All"."""
    out = []
    for level in levels:
        tips = sorted((lengths[tip], tip) for tip, lvl in tallies if lvl == level
                      and tip in lengths)  # fmt: skip
        if not tips:
            continue
        n = len(tips)
        groups: list[tuple[str, list[tuple[float, str]]]] = []
        if n >= MIN_QUARTER_TIPS:
            groups = [(QUARTERS[q], tips[q * n // 4 : (q + 1) * n // 4]) for q in range(4)]
        groups.append(("All", tips))
        for label, members in groups:
            total = Tally()
            for _, tip in members:
                total.add(tallies[(tip, level)])
            values = [length for length, _ in members]
            out.append(Row(level, label, values[0], values[-1], _median(values), total))
    return out


def _median(values: list[float]) -> float:
    n = len(values)
    return values[n // 2] if n % 2 else (values[n // 2 - 1] + values[n // 2]) / 2


# --- the scatters ---------------------------------------------------------------------

SCATTER_W, SCATTER_H = 420, 300
LEFT, RIGHT, TOP, BOTTOM = 58, 16, 14, 44
COLOURS = {"e0": "#2a78d6", "e1": "#eb6834"}  # as esr_views.COLOURS: dark.css remaps them


@dataclass
class Dot:
    x: float
    y: float
    title: str


@dataclass
class Series:
    level: str
    colour: str
    dots: list[Dot]
    line: str = ""  # the quarters' figure at each quarter's median pendant
    legend_y: float = 0.0


@dataclass
class Scatter:
    key: str
    y_label: str
    width: int = SCATTER_W
    height: int = SCATTER_H
    left: int = LEFT
    top: int = TOP
    right_edge: int = SCATTER_W - RIGHT
    bottom_edge: int = SCATTER_H - BOTTOM
    zero_y: float = 0.0
    x_ticks: list[dict[str, Any]] = field(default_factory=list)
    y_ticks: list[dict[str, Any]] = field(default_factory=list)
    series: list[Series] = field(default_factory=list)

    @property
    def legend_x(self) -> int:
        return self.left + 12

    @property
    def x_middle(self) -> float:
        return (self.left + self.right_edge) / 2

    @property
    def y_middle(self) -> float:
        return (self.top + self.bottom_edge) / 2


def _log_axis(values: list[float]) -> tuple[float, float, list[float]]:
    """A log10 range covering the values (zero lengths clamped to 1/10,000 of the
    longest) and its ticks: decades, with 2 and 5 between them over a short range."""
    hi = max(values)
    positive = [v for v in values if v > 0]
    lo = max(min(positive) if positive else hi, hi / 1e4) if hi > 0 else 1e-3
    hi = hi if hi > 0 else 1.0
    a, b = math.floor(math.log10(lo)), math.ceil(math.log10(hi))
    if b == a:
        b += 1
    steps = (1, 2, 5) if b - a <= 3 else (1,)
    ticks = [m * 10.0**e for e in range(a, b + 1) for m in steps]
    lo_t = max(t for t in ticks if t <= lo)
    hi_t = min(t for t in ticks if t >= hi)
    return math.log10(lo_t), math.log10(hi_t), [t for t in ticks if lo_t <= t <= hi_t]


def _symmetric_axis(values: list[float], step: float) -> tuple[float, list[float]]:
    """A range ±m around 0 covering the values, ticked at 0, ±k and (mostly) ±2k, where
    k is a whole number of steps."""
    biggest = max((abs(v) for v in values if not math.isnan(v)), default=0.0)
    n = max(1, math.ceil(biggest / step - 1e-9))  # steps needed
    per_tick = math.ceil(n / 2)
    k = per_tick * step
    count = math.ceil(n / per_tick)  # 1 or 2 ticks each side
    return count * k, [round(i * k, 10) for i in range(-count, count + 1)]


def _number(v: float) -> str:
    if not v:
        return "0"
    return f"{v:+.3f}".rstrip("0") if abs(v) < 0.01 else f"{v:+.2f}"


def _length(v: float) -> str:
    return f"{v:.3g}"


def scatter(
    key: str,
    y_label: str,
    points: dict[str, list[tuple[float, float, str]]],
    lines: dict[str, list[tuple[float, float]]],
    step: float,
) -> Scatter | None:
    """Per level, one dot per sequence (pendant on a log axis, a signed figure), the
    quarters' figure as a line over them, and 0 drawn across."""
    every = [p for dots in points.values() for p in dots]
    if len(every) < MIN_SCATTER_TIPS:
        return None
    out = Scatter(key, y_label)
    lo, hi, x_ticks = _log_axis([x for x, _, _ in every])
    m, y_ticks = _symmetric_axis([y for _, y, _ in every] + [y for ln in lines.values()
                                                             for _, y in ln], step)  # fmt: skip
    w, h = out.right_edge - out.left, out.bottom_edge - out.top

    def x(v: float) -> float:
        v = min(max(math.log10(v) if v > 0 else lo, lo), hi)
        return out.left + (v - lo) / (hi - lo) * w

    def y(v: float) -> float:
        return out.bottom_edge - (v + m) / (2 * m) * h

    out.zero_y = y(0)
    out.x_ticks = [{"x": x(t), "label": _length(t)} for t in x_ticks]
    out.y_ticks = [{"y": y(t), "label": _number(t)} for t in y_ticks]
    for level in points:
        dots = [Dot(round(x(px), 1), round(y(py), 1), title) for px, py, title in points[level]]
        line = " ".join(f"{'M' if i == 0 else 'L'}{x(px):.1f},{y(py):.1f}"
                        for i, (px, py) in enumerate(lines.get(level, [])))  # fmt: skip
        legend = out.top + 12 + 20 * len(out.series)
        out.series.append(Series(level, COLOURS[level], dots, line, legend))
    return out


# --- the section ----------------------------------------------------------------------


def _figure(v: float, digits: int = 3, signed: bool = False) -> str:
    if math.isnan(v):
        return "–"
    return f"{v:+.{digits}f}" if signed else f"{v:.{digits}f}"


def _percent(v: float) -> str:
    return "–" if math.isnan(v) else f"{round(100 * v)}%"


def context(tallies: Tallies, tips_text: str, levels: list[str]) -> dict[str, Any] | None:
    """The section: per level the headline additions, the quarters table, the gap
    scatter and (with both levels) E1 minus E0. None when there is nothing to show."""
    lengths = pendants(tips_text)
    levels = [lvl for lvl in levels if any(level == lvl for _, level in tallies)]
    if not lengths or not levels:
        return None
    rows = quarters(tallies, lengths, levels)
    table = [{
        "level": r.level.upper(),
        "label": r.label,
        "pendant": f"{_length(r.low)}–{_length(r.high)}" if r.low != r.high else _length(r.low),
        "tips": r.t.tips,
        "sites": r.t.sites,
        "accuracy": _figure(r.t.accuracy),
        "mean_pp": _figure(r.t.mean_pp),
        "gap": _figure(r.t.gap, signed=True),
        "ece": _figure(r.t.ece),
        "confident_wrong": _figure(r.t.confident_wrong_per_tip, 1),
        "caught": _percent(r.t.caught_share),
        "all": r.label == "All",
    } for r in rows]  # fmt: skip
    overall = {r.level: r.t for r in rows if r.label == "All"}
    headline = {level: {
        # The summary's log score floored a printed 0 at 1e-300 before 2026-11-02; this
        # one uses the posterior's printed precision (metrics.LOG_FLOOR).
        "log_score": t.log_score,
        "ece": _figure(t.ece),
        "confident_wrong": _figure(t.confident_wrong_per_tip, 1),
        "caught": _percent(t.caught_share),
    } for level, t in overall.items()}  # fmt: skip
    quartered = [r for r in rows if r.label != "All"]

    gap_points = {level: [(lengths[tip], tallies[(tip, lvl)].gap,
                           f"{tip}: pendant {_length(lengths[tip])}, accuracy "
                           f"{_figure(tallies[(tip, lvl)].accuracy)}, mean posterior "
                           f"{_figure(tallies[(tip, lvl)].mean_pp)}, gap "
                           f"{_figure(tallies[(tip, lvl)].gap, signed=True)}")
                          for tip, lvl in tallies if lvl == level and tip in lengths]
                  for level in levels}  # fmt: skip
    gap_lines = {level: [(r.median, r.t.gap) for r in quartered if r.level == level]
                 for level in levels}  # fmt: skip
    gaps = scatter("gap", "Mean posterior − accuracy", gap_points, gap_lines, 0.05)

    difference = None
    if "e0" in levels and "e1" in levels:
        both = [tip for tip, lvl in tallies if lvl == "e0" and (tip, "e1") in tallies
                and tip in lengths]  # fmt: skip
        diffs = {tip: tallies[(tip, "e1")].log_score - tallies[(tip, "e0")].log_score
                 for tip in both}  # fmt: skip
        points = {"e1": [(lengths[tip], d, f"{tip}: pendant {_length(lengths[tip])}, E1 − E0 "
                          f"log score {_figure(d, 4, signed=True)} a site")
                         for tip, d in diffs.items()]}  # fmt: skip
        e0 = {r.label: r for r in quartered if r.level == "e0"}
        line = [(r.median, r.t.log_score - e0[r.label].t.log_score) for r in quartered
                if r.level == "e1" and r.label in e0]  # fmt: skip
        difference = scatter("difference", "E1 − E0 log score a site", points, {"e1": line},
                             0.005)  # fmt: skip
    return {
        "table": table,
        "quartered": bool(quartered),
        "headline": headline,
        "gaps": gaps,
        "difference": difference,
        "confident": f"{CONFIDENT:g}",
        "altall": f"{ALTALL:g}",
        "min_tips": MIN_QUARTER_TIPS,
    }


def for_job(
    prefix: str, keys: set[str], tips_text: str, levels: list[str]
) -> dict[str, Any] | None:
    """The section for a finished ESR job whose results are under `prefix`."""
    sites_key = "esr/esr_sites.tsv"
    if sites_key not in keys or not tips_text:
        return None
    try:
        tallies = tally(prefix + sites_key)
    except Exception:
        return None  # absent, never broken
    return context(tallies, tips_text, [lvl for lvl in esr.LEVELS if lvl in levels])
