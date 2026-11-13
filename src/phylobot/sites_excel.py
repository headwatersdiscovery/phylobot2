"""Excel workbooks of the viewers' tables, coloured as on the page: the ancestor
viewer's site table (`workbook`), a traced mutational path (`path_workbook`), a
library's members viewer (`library_workbook`) and a job page's results table
(`results_workbook`).

Sheet "Sites": one row per alignment column; the chosen extant sequences' residues and
residue numbers; then each shown ancestral state as a pair of cells (state, posterior
probability as a number, so it sorts and filters), both filled with the tile's colour.
Sheet "About": what was shown (job, alignment, model, node, gap version, palette), the
node's summary statistics, and a colour key.
"""

from functools import cache
from io import BytesIO
from typing import Any

from openpyxl import Workbook
from openpyxl.cell import Cell, WriteOnlyCell
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

from phylobot import ancestors, library_glyphs, paths

GAP_FILL = "ECEFF3"  # as td.tile.gap on the page
HEADER_FILL = "F4F6F9"


CENTRED = Alignment(horizontal="center")
HEADER_ALIGN = Alignment(wrap_text=True, vertical="bottom", horizontal="center")


# One style object per colour: openpyxl matches styles cell by cell, so fresh objects
# for every cell make a thousand-row table several times slower.
@cache
def fill(hex_colour: str) -> PatternFill:
    return PatternFill("solid", start_color=hex_colour, end_color=hex_colour)


@cache
def font(hex_colour: str, bold: bool = False) -> Font:
    return Font(color=hex_colour, bold=bold)


def hex_of(rgb: tuple[int, int, int]) -> str:
    return "{:02X}{:02X}{:02X}".format(*rgb)


def tile_style(tile: ancestors.Tile, palette: str) -> tuple[PatternFill | None, Font]:
    if tile.pp is not None:
        rgb = ancestors.rgb(tile.pp, palette)
        return fill(hex_of(rgb)), font(ancestors.text_colour(rgb).lstrip("#").upper())
    if tile.kind == "gap":
        return fill(GAP_FILL), font("5D6675")
    return None, font("111111", bold=True)  # PRANK/GRASP residue: no probability, uncoloured


def workbook(context: dict[str, Any], about: list[tuple[str, str]]) -> bytes:
    """`context` is the viewer's (views.ancestor_context); `about` labels what it shows.
    Written in openpyxl's write-only (streaming) mode, several times faster."""
    rows: list[ancestors.SiteRow] = context["rows"]
    palette: str = context["palette"]
    book = Workbook(write_only=True)
    sheet = book.create_sheet("Sites")

    def cell(
        value: Any,
        fill_: PatternFill | None = None,
        font_: Font | None = None,
        alignment: Alignment | None = CENTRED,
        number_format: str = "",
    ) -> Cell:
        c = WriteOnlyCell(sheet, value=value)
        if fill_ is not None:
            c.fill = fill_
        if font_ is not None:
            c.font = font_
        if alignment is not None:
            c.alignment = alignment
        if number_format:
            c.number_format = number_format
        return c

    tiles = max((len(row.tiles) for row in rows), default=1)
    calibrated = bool(context.get("calibration"))  # a T column after the site's
    first_tile = 2 + calibrated + 2 * len(context["extants"])
    # Layout settings come before the first row in write-only mode.
    sheet.freeze_panes = f"{get_column_letter(first_tile)}2"  # header and site columns stay
    sheet.column_dimensions["A"].width = 9
    for column in range(2, first_tile):
        sheet.column_dimensions[get_column_letter(column)].width = 12
    for column in range(first_tile, first_tile + 2 * tiles):
        sheet.column_dimensions[get_column_letter(column)].width = 7

    header = [f"Site in {context['aligner_label']}"] + (["T"] if calibrated else [])
    for name in context["extants"]:
        header += [f"{name} residue", f"{name} site"]
    for k in range(1, tiles + 1):
        header += [f"State {k}", f"PP {k}"]
    heading = fill(HEADER_FILL), font("111111", bold=True), HEADER_ALIGN
    sheet.append([cell(h, *heading) for h in header])
    for row in rows:
        cells = [cell(row.column)]
        if calibrated:
            cells.append(cell(row.temperature, number_format="0.00"))
        for residue, number in row.extants:
            cells += [cell(residue), cell(number)]
        for tile in row.tiles:
            background, style_font = tile_style(tile, palette)
            pp = round(tile.pp, 4) if tile.pp is not None else None
            state = f"{tile.state} ({tile.note})" if tile.note else tile.state  # codon (aa)
            cells += [cell(state, background, style_font),
                      cell(pp, background, style_font, number_format="0.00")]  # fmt: skip
        sheet.append(cells)

    info = book.create_sheet("About")
    info.column_dimensions["A"].width = 26
    info.column_dimensions["B"].width = 70
    bold = font("111111", bold=True)
    for label, value in about:
        info.append([cell(label, font_=bold, alignment=None), value])
    info.append([])
    info.append([cell("Colour key", font_=bold, alignment=None),
                 "posterior probability (as the viewer's palette)"])  # fmt: skip
    for pp in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        rgb = ancestors.rgb(pp, palette)
        text = font(ancestors.text_colour(rgb).lstrip("#").upper())
        info.append(["", cell(pp, fill(hex_of(rgb)), text, alignment=None, number_format="0.0")])
    info.append(["", cell("gap", fill(GAP_FILL), alignment=None)])

    out = BytesIO()
    book.save(out)
    return out.getvalue()


EDGE = Side(style="medium", color="111111")
CHANGED = Border(left=EDGE, right=EDGE, top=EDGE, bottom=EDGE)
DASHED = Side(style="mediumDashed", color="111111")
SILENT = Border(left=DASHED, right=DASHED, top=DASHED, bottom=DASHED)  # silent change


def path_workbook(
    table: paths.PathTable,
    palette: str,
    about: list[tuple[str, str]],
    support: dict[tuple[int, int], tuple[int, int, str]] | None = None,
) -> bytes:
    """Sheets: "Path" (residues, coloured by posterior probability, changed cells boxed),
    "Probabilities" (the same cells' probabilities, as numbers), "Changes" (one row per
    change, branch by branch), and "About". `support` (phylobot.path_compare, from the
    path page) adds to each change how many of the compared combinations make it too:
    {(branch index, alignment column): (the same in, compared with, in words)}."""
    book = Workbook(write_only=True)
    bold = font("111111", bold=True)
    heading = fill(HEADER_FILL)

    def write_grid(title: str, numbers: bool) -> None:
        sheet = book.create_sheet(title)
        sheet.freeze_panes = "B3"
        sheet.column_dimensions["A"].width = 30
        for column in range(2, len(table.columns) + 2):
            sheet.column_dimensions[get_column_letter(column)].width = 6 if numbers else 4

        def cell(value: Any, fill_: PatternFill | None = None, font_: Font | None = None,
                 border: Border | None = None, number_format: str = "") -> Cell:  # fmt: skip
            c = WriteOnlyCell(sheet, value=value)
            c.alignment = CENTRED
            if fill_ is not None:
                c.fill = fill_
            if font_ is not None:
                c.font = font_
            if border is not None:
                c.border = border
            if number_format:
                c.number_format = number_format
            return c

        first = table.rows[0].name
        sheet.append([cell("Alignment column", heading, bold)]
                     + [cell(c, heading, bold) for c in table.columns])  # fmt: skip
        sheet.append([cell(f"Site in {first}", heading, bold)]
                     + [cell(n, heading) for n in table.first_numbers])  # fmt: skip
        compared = False
        for row in table.rows:
            if not row.on_path and not compared:
                sheet.append([cell("For comparison", font_=bold)])
                compared = True
            cells = [cell(row.name, font_=bold)]
            for c in row.cells:
                look: tuple[PatternFill | None, Font | None]
                if c.pp is not None:
                    rgb = ancestors.rgb(c.pp, palette)
                    look = fill(hex_of(rgb)), font(ancestors.text_colour(rgb).lstrip("#").upper())
                elif c.char == "–":
                    look = fill(GAP_FILL), font("5D6675")
                else:
                    look = None, None
                value: Any = (c.pp if c.pp is not None else None) if numbers else c.char
                cells.append(cell(
                    value, fill_=look[0], font_=look[1],
                    border=CHANGED if c.changed else SILENT if c.silent else None,
                    number_format="0.00" if numbers else "",
                ))  # fmt: skip
            sheet.append(cells)

    write_grid("Path", numbers=False)
    write_grid("Probabilities", numbers=True)

    changes = book.create_sheet("Changes")
    changes.append([
        WriteOnlyCell(changes, value=h) for h in (
            "Branch", "Alignment column", "Change", "Kind", "Before", "After",
            "PP before", "PP after", *(["Codons"] if table.codons else []),
            *(["Same in", "Compared with", "Support"] if support else []),
        )
    ])  # fmt: skip
    for index, branch in enumerate(table.branches):
        for change in branch.changes:
            found = (support or {}).get((index, change.column), (None, None, None))
            changes.append([
                f"{branch.parent} → {branch.child}", change.column, change.label, change.kind,
                change.before, change.after, change.pp_before, change.pp_after,
                *([change.codons] if table.codons else []),
                *(found if support else []),
            ])  # fmt: skip

    info = book.create_sheet("About")
    info.column_dimensions["A"].width = 22
    info.column_dimensions["B"].width = 80
    for label, value in about:
        info.append([label, value])

    def counts(b: paths.BranchChanges) -> str:
        silent = f", {b.synonymous} silent" if table.codons else ""
        return f"{b.substitutions} substitutions{silent}, {b.gains} gained, {b.losses} lost"

    info.append(["Branches", "; ".join(f"{b.parent} → {b.child}: {counts(b)}"
                                        for b in table.branches)])  # fmt: skip
    legend = "numbered by alignment column; K300R: K became R; -300R: R gained; K300-: K lost."
    if table.codons:
        legend += " L300L: silent (same amino acid, another codon); dashed cells."
    info.append(["Changes", legend + " Boxed cells differ from the row above."])

    out = BytesIO()
    book.save(out)
    return out.getvalue()


def _faint(rgb: tuple[int, int, int], alpha: float = 0.35) -> tuple[int, int, int]:
    """A colour as the viewer draws a dimmed row: at `alpha` over white."""
    r, g, b = (round(255 - (255 - c) * alpha) for c in rgb)
    return r, g, b


def library_workbook(
    frames: list[library_glyphs.Frame],
    members: list[dict[str, Any]],
    palette: str,
    about: list[tuple[str, str]],
) -> bytes:
    """The library page's members viewer (L5c) as a workbook. Sheets: "Sequences" (each
    frame's rows as on the page, one cell per column or residue, coloured by posterior
    probability, members' changes from the ML sequence boxed, members the budget left out
    faded), "Probabilities" (the same cells' probabilities, as numbers), "Members" (the
    table beneath the viewer) and "About"."""
    book = Workbook(write_only=True)
    bold = font("111111", bold=True)
    heading = fill(HEADER_FILL)
    units = max((f.units for f in frames), default=0)

    def write_grid(title: str, numbers: bool) -> None:
        sheet = book.create_sheet(title)
        sheet.freeze_panes = "B1"
        sheet.column_dimensions["A"].width = 30
        for column in range(2, units + 2):
            sheet.column_dimensions[get_column_letter(column)].width = 6 if numbers else 3

        def cell(value: Any, fill_: PatternFill | None = None, font_: Font | None = None,
                 border: Border | None = None, number_format: str = "",
                 alignment: Alignment | None = CENTRED) -> Cell:  # fmt: skip
            c = WriteOnlyCell(sheet, value=value)
            if alignment is not None:
                c.alignment = alignment
            if fill_ is not None:
                c.fill = fill_
            if font_ is not None:
                c.font = font_
            if border is not None:
                c.border = border
            if number_format:
                c.number_format = number_format
            return c

        for n, frame in enumerate(frames):
            if n:
                sheet.append([])
            label = f"{frame.title} · {frame.members} member{'s' if frame.members != 1 else ''}"
            sheet.append([cell(label, font_=bold, alignment=None)])
            what = "Alignment column" if frame.unit == "column" else "Residue"
            sheet.append([cell(what, heading, bold)]
                         + [cell(i, heading, bold) for i in range(1, frame.units + 1)])  # fmt: skip
            for row in frame.rows:
                marked = set(row.marked)
                name_font = font("8A929E") if row.dim else bold if not row.member else None
                name = f"{row.name} (left out by the budget)" if row.dim else row.name
                cells = [cell(name, font_=name_font, alignment=None)]
                for i, (state, pp) in enumerate(zip(row.states, row.pp, strict=True)):
                    look: tuple[PatternFill | None, Font | None]
                    if pp >= 0:
                        rgb = ancestors.rgb(pp / 100, palette)
                        text = ancestors.text_colour(rgb).lstrip("#").upper()
                        if row.dim:
                            rgb, text = _faint(rgb), "8A929E"
                        look = fill(hex_of(rgb)), font(text)
                    elif state == "-":
                        look = fill(GAP_FILL), font("5D6675")
                    else:
                        look = None, font("8A929E" if row.dim else "111111")
                    value: Any = (pp / 100 if pp >= 0 else None) if numbers else state
                    cells.append(cell(
                        value, fill_=look[0], font_=look[1],
                        border=CHANGED if i in marked else None,
                        number_format="0.00" if numbers else "",
                    ))  # fmt: skip
                sheet.append(cells)

    write_grid("Sequences", numbers=False)
    write_grid("Probabilities", numbers=True)

    table = book.create_sheet("Members")
    table.freeze_panes = "A2"
    for column, width in zip("ABCDEFGHI", (30, 22, 9, 9, 9, 9, 9, 14, 12), strict=True):
        table.column_dimensions[column].width = width
    names = ("Name", "Variant", "Support", "Weight", "Mean pp", "Length", "Changes", "From",
             "In the set")  # fmt: skip
    table.append([WriteOnlyCell(table, value=h) for h in names])
    for row in members:
        member = row["member"]
        mean_pp = row["mean_pp"]
        table.append([
            row["name"], row["variant_label"], member.support, round(member.weight, 3),
            round(mean_pp, 3) if mean_pp is not None else None, member.length,
            row["differences"], member.source.node, "yes" if row["picked"] else "no",
        ])  # fmt: skip

    info = book.create_sheet("About")
    info.column_dimensions["A"].width = 22
    info.column_dimensions["B"].width = 80
    for label, value in about:
        info.append([label, value])
    info.append(["Cells", "each member on its ancestor's columns, coloured by that ancestor's "
                 "posterior probability for the state it carries; boxed cells differ from "
                 "the ML sequence; faded rows were left out by the budget."])  # fmt: skip
    info.append([])
    info.append(["Colour key", "posterior probability (as the page's palette)"])
    for pp in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        rgb = ancestors.rgb(pp, palette)
        c = WriteOnlyCell(info, value=pp)
        c.fill = fill(hex_of(rgb))
        c.font = font(ancestors.text_colour(rgb).lstrip("#").upper())
        c.number_format = "0.0"
        info.append(["", c])

    out = BytesIO()
    book.save(out)
    return out.getvalue()


def results_workbook(
    columns: list[str], rows: list[list[str]], about: list[tuple[str, str]], name: str = "Models"
) -> bytes:
    """A job page's results table (one alignment, one data type): a sheet `name` ("Models",
    a row per combination) with its statistics as numbers, unrounded; sheet "About"."""
    book = Workbook(write_only=True)
    sheet = book.create_sheet(name)
    sheet.freeze_panes = "C2"
    for column, name in enumerate(columns, start=1):
        wide = name in ("Combination", "Model", "Fitted model", "Tree source")
        sheet.column_dimensions[get_column_letter(column)].width = 28 if wide else 12
    bold, heading = font("111111", bold=True), fill(HEADER_FILL)

    def header(value: str) -> Cell:
        c = WriteOnlyCell(sheet, value=value)
        c.font, c.fill, c.alignment = bold, heading, HEADER_ALIGN
        return c

    sheet.append([header(name) for name in columns])
    for row in rows:
        sheet.append([_number(value) for value in row])
    info = book.create_sheet("About")
    info.column_dimensions["A"].width = 22
    info.column_dimensions["B"].width = 80
    for label, value in about:
        info.append([label, value])
    out = BytesIO()
    book.save(out)
    return out.getvalue()


def _number(text: str) -> Any:
    """A statistic as a number, so it sorts and filters; anything else as text."""
    try:
        value = float(text)
    except ValueError:
        return text
    return int(value) if value.is_integer() and "." not in text and "e" not in text else value
