// Sortable tables: a header button sorts its column, ascending first, then descending
// on a second click. Delegated from the document, so it keeps working after HTMX
// replaces the job's status block. Numbers sort as numbers; a cell without one ("–")
// always goes last. A cell's data-sort, when it has one, is sorted instead of its text
// (a duration's seconds, a name without its hints). A table.grouped sorts its <tbody>s
// as units, each by its first row (a job and its stage rows stay together). Without
// this script the table keeps the server's order.
(() => {
  const value = (row, index) => {
    const cell = row.cells[index];
    if (!cell) return "";
    return (cell.dataset.sort !== undefined ? cell.dataset.sort : cell.textContent).trim();
  };
  // What moves when the table sorts: [key row, element to move] pairs.
  const units = (table) => table.classList.contains("grouped")
    ? [...table.tBodies].filter((b) => b.rows.length).map((b) => [b.rows[0], b])
    : [...(table.tBodies[0] ? table.tBodies[0].rows : [])].map((r) => [r, r]);

  document.addEventListener("click", (event) => {
    const button = event.target.closest && event.target.closest("th.sortable button.sort");
    if (!button) return;
    const th = button.closest("th");
    const table = th.closest("table");
    const items = units(table);
    if (!items.length) return;
    const index = [...th.parentElement.children].indexOf(th);
    const ascending = th.getAttribute("aria-sort") !== "ascending";
    const numeric = items.every(([row]) => {
      const text = value(row, index);
      return text === "–" || text === "" || Number.isFinite(parseFloat(text));
    });
    const direction = ascending ? 1 : -1;
    items.sort(([a], [b]) => {
      const x = value(a, index), y = value(b, index);
      if (numeric) {
        const p = parseFloat(x), q = parseFloat(y);
        if (!Number.isFinite(p) || !Number.isFinite(q)) {
          return Number.isFinite(p) ? -1 : Number.isFinite(q) ? 1 : 0;
        }
        return (p - q) * direction;
      }
      return x.localeCompare(y, undefined, { numeric: true }) * direction;
    });
    const parent = items[0][1].parentElement;
    parent.append(...items.map(([, element]) => element));
    table.querySelectorAll("th[aria-sort]").forEach((other) => other.removeAttribute("aria-sort"));
    th.setAttribute("aria-sort", ascending ? "ascending" : "descending");
  });
})();

// A search box over a table (input[data-table-search] naming the table's id): hides the
// rows (or, in a table.grouped, the <tbody>s, by their first row) whose text does not
// contain every word typed. Case-insensitive; an empty box shows everything.
(() => {
  document.addEventListener("input", (event) => {
    const box = event.target.closest && event.target.closest("input[data-table-search]");
    if (!box) return;
    const table = document.getElementById(box.dataset.tableSearch);
    if (!table) return;
    const words = box.value.toLowerCase().split(/\s+/).filter(Boolean);
    const grouped = table.classList.contains("grouped");
    const items = grouped ? [...table.tBodies] : [...(table.tBodies[0] ? table.tBodies[0].rows : [])];
    let shown = 0;
    for (const item of items) {
      const row = grouped ? item.rows[0] : item;
      const text = row ? row.textContent.toLowerCase() : "";
      const match = words.every((word) => text.includes(word));
      item.hidden = !match;
      if (match) shown += 1;
    }
    const count = document.querySelector(`[data-table-search-count="${CSS.escape(box.dataset.tableSearch)}"]`);
    if (count) count.textContent = words.length ? `${shown} of ${items.length}` : "";
  });
})();
