"""The Uncertainty page's tables as downloads (CSV and Excel): every site of one ancestor,
and every ancestor of the run.

Pure over a run's parsed files (`uncertainty_views.rows_of`): the page shows the sites
that vary most, the download has all of them. A table is a header and rows of plain
values (numbers stay numbers, so a spreadsheet sorts and filters them).
"""

from __future__ import annotations

import csv
import io
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from phylobot import uncertainty_run as run_mod

COMPONENTS = (
    ("within", "Posterior"),
    ("alignment", "Alignment"),
    ("tree", "Tree"),
    ("model", "Model"),
    ("gaps", "Gap placement"),
    ("other", "Other"),
)
VALUES = ("total", *(k for k, _ in COMPONENTS))


@dataclass(frozen=True)
class Table:
    title: str  # the sheet's name
    header: list[str]
    rows: list[list[Any]]


def site_table(
    site_rows: list[dict[str, str]],
    drivers: list[dict[str, str]],
    mode: str,
    factors: list[dict[str, Any]],
) -> Table:
    """One row per reference column: agreement, the whole site's split, its two parts'
    totals, each part's split, then what each level of every varied choice reconstructs
    there ("K 0.93"; blank where the column hardly varies between the choices)."""
    varied = [f for f in factors if f.get("varied")]
    levels = [
        (f["key"], lv["key"], f"{f['label']}: {lv['label']}") for f in varied for lv in f["levels"]
    ]
    said: dict[tuple[int, str, str], str] = {}
    for d in drivers:
        said[(int(d["column"]), d["factor"], d["level"])] = f"{d['state']} {float(d['p']):.2f}"
    # Every number is a variance of the column's one-hot state, so one unit throughout:
    # the chance two reconstructions chosen at random differ there (`about`'s "Unit").
    header = ["Column", "Alignment agreement", "Variance"]
    header += [f"Variance: {label}" for _, label in COMPONENTS]
    header += ["Variance: residue or gap", "Variance: which residue"]
    for name in ("residue or gap", "which residue"):
        header += [f"Variance, {name}: {label}" for _, label in COMPONENTS]
    header += [f"Reconstruction, {label}" for _, _, label in levels]
    rows = []
    for r in site_rows:
        if r["mode"] != mode:
            continue
        column = int(r["column"])
        parts = {p: {k: float(r[f"{p}_{k}"]) for k in VALUES} for p in ("presence", "residue")}
        whole = {k: parts["presence"][k] + parts["residue"][k] for k in VALUES}
        row: list[Any] = [column, float(r["agreement"]), round(whole["total"], 5)]
        row += [round(whole[k], 5) for k, _ in COMPONENTS]
        row += [parts["presence"]["total"], parts["residue"]["total"]]
        for p in ("presence", "residue"):
            row += [parts[p][k] for k, _ in COMPONENTS]
        row += [said.get((column, factor, level), "") for factor, level, _ in levels]
        rows.append(row)
    rows.sort(key=lambda row: row[0])
    return Table("Sites", header, rows)


def ancestor_table(node_rows: list[dict[str, str]], mode: str) -> Table:
    """One row per ancestor, in the run's order (the reference tree, root first)."""
    header = ["Ancestor", "Present in", "Of combinations", "Sequences below", "Total"]
    header += ["Residue or gap", "Which residue", *(label for _, label in COMPONENTS)]
    by_node: dict[str, dict[str, dict[str, str]]] = {}
    for r in node_rows:
        if r["mode"] == mode:
            by_node.setdefault(r["node"], {})[r["part"]] = r
    rows = []
    for node, parts in by_node.items():
        whole = parts["all"]
        rows.append([
            node, int(whole["present"]), int(whole["of"]), int(whole["leaves"]),
            float(whole["total"]), float(parts["presence"]["total"]),
            float(parts["residue"]["total"]), *(float(whole[k]) for k, _ in COMPONENTS),
        ])  # fmt: skip
    return Table("Ancestors", header, rows)


def about(
    job_name: str, summary: dict[str, Any], mode: str, settings: dict[str, Any]
) -> list[tuple[str, str]]:
    weights = settings.get("weighting", "flat")
    multipliers = settings.get("multipliers") or {}
    return [
        ("Job", job_name),
        ("Compared", run_mod.MODE_LABELS.get(mode, mode)),
        ("Unit", "Variance of the one-hot state (20 residues and the gap). Per column it is "
                 "the chance that two reconstructions chosen at random differ there (0 to "
                 "0.95); per ancestor, summed over columns, the expected number of positions "
                 "at which they differ. Every part (posterior, each choice) is in the same "
                 "unit, and the parts add up to the total."),
        ("Reference", str(summary.get("reference", ""))),
        ("Ancestors matched", "by clade on outgroup-rooted trees" if summary.get("by_clade")
         else "by the subtrees that meet at them (no outgroup)"),
        ("Weights", weights + (f"; multipliers {multipliers}" if multipliers else "")),
        ("Reconstructions", ", ".join(c.get("key", "") for c in summary.get("cells", []))),
        ("Not varied", ", ".join(
            f["label"] for f in summary.get("factors", []) if not f.get("varied")
        ) or "none"),
        ("Reconstruction columns", "Each choice's reconstructions pooled: most probable "
                                   "state and its "
                           "probability, '-' a gap; blank where the column varies by less than "
                           f"{run_mod.DRIVER_FLOOR} between the choices"),
        ("Method", "Law of total variance on the one-hot state, Shapley shares"),
    ]  # fmt: skip


def as_csv(table: Table) -> str:
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(table.header)
    writer.writerows(table.rows)
    return out.getvalue()


def as_xlsx(table: Table, notes: Sequence[tuple[str, str]]) -> bytes:
    """The table on its own sheet (header bold, frozen), and an "About" sheet."""
    from openpyxl import Workbook  # web tier only: the pipeline image has no openpyxl
    from openpyxl.cell import WriteOnlyCell
    from openpyxl.styles import Alignment, Font, PatternFill

    book = Workbook(write_only=True)
    sheet = book.create_sheet(table.title)
    sheet.freeze_panes = "B2"
    bold = Font(bold=True)
    shade = PatternFill("solid", start_color="F4F6F9", end_color="F4F6F9")
    wrap = Alignment(wrap_text=True, vertical="bottom")
    heading = []
    for name in table.header:
        cell = WriteOnlyCell(sheet, value=name)
        cell.font, cell.fill, cell.alignment = bold, shade, wrap
        heading.append(cell)
    sheet.append(heading)
    for row in table.rows:
        sheet.append(row)
    notes_sheet = book.create_sheet("About")
    notes_sheet.column_dimensions["A"].width = 22
    notes_sheet.column_dimensions["B"].width = 100
    for key, value in notes:
        label = WriteOnlyCell(notes_sheet, value=key)
        label.font = bold
        notes_sheet.append([label, value])
    out = io.BytesIO()
    book.save(out)
    return out.getvalue()
