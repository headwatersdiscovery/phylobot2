// The alignment editor's edit bar (phylobot.alignment_editor_views). alignment.js draws
// the alignment and reports the selection ("aln:select": columns and rows as shown,
// 0-based) and clicks on the Name header ("aln:sort"). An edit fills one of the form's
// one-off fields with 1-based positions as shown and submits the form; the server turns
// them into the stored positions of the hidden cut/drop/sort fields and swaps the editor.
//
// Undo: before each edit the current cut/drop/sort are pushed on a stack per alignment
// version (kept in sessionStorage, so it survives a reload; an overwrite starts a new
// version, whose positions differ); Undo puts the last back and submits.
//
// "Save (overwrite)" asks first (delete.js's confirmDialog), then posts the save form.
(() => {
  const FIELDS = ["cut", "drop", "sort"];
  const ONE_OFF = ["cut_columns", "keep_columns", "drop_rows", "keep_rows", "sort_toggle"];
  const MAX_UNDO = 200;
  let selection = { columns: [], rows: [] };

  const editorForm = () => document.querySelector("form[data-editor]");
  const key = (form) => `phylobot.alignment-undo.${form.dataset.editor}.v${form.dataset.version}`;
  function stack(form) {
    try { return JSON.parse(sessionStorage.getItem(key(form)) || "[]"); } catch { return []; }
  }
  function keep(form, entries) {
    try { sessionStorage.setItem(key(form), JSON.stringify(entries.slice(-MAX_UNDO))); } catch { /* private mode */ }
    syncUndo(form, entries.length);
  }
  const state = (form) => Object.fromEntries(FIELDS.map((f) => [f, form.elements[f].value]));
  function remember(form) {
    const entries = stack(form);
    const now = state(form);
    const last = entries[entries.length - 1];
    if (!last || FIELDS.some((f) => last[f] !== now[f])) entries.push(now);
    keep(form, entries);
  }
  function syncUndo(form, count) {
    const button = form.querySelector("[data-undo]");
    if (button) button.disabled = count === 0;
  }

  // "1-3,7" (1-based) from sorted 0-based positions.
  function ranges(positions) {
    const parts = [];
    let start = null, prev = null;
    for (const p of positions) {
      if (start === null) { start = prev = p; continue; }
      if (p === prev + 1) { prev = p; continue; }
      parts.push(start === prev ? `${start + 1}` : `${start + 1}-${prev + 1}`);
      start = prev = p;
    }
    if (start !== null) parts.push(start === prev ? `${start + 1}` : `${start + 1}-${prev + 1}`);
    return parts.join(",");
  }

  function submitEdit(field, positions) {
    const form = editorForm();
    if (!form || !positions.length) return;
    remember(form);
    for (const f of ONE_OFF) form.elements[f].value = "";
    form.elements[field].value = ranges(positions);
    form.requestSubmit();
  }
  function undo() {
    const form = editorForm();
    if (!form) return;
    const entries = stack(form);
    const previous = entries.pop();
    if (!previous) return;
    keep(form, entries);
    for (const f of FIELDS) form.elements[f].value = previous[f];
    for (const f of ONE_OFF) form.elements[f].value = "";
    form.requestSubmit();
  }
  function clearSelection() {
    const wrap = editorForm()?.querySelector(".aln-wrap");
    if (wrap) wrap.dispatchEvent(new CustomEvent("aln:update", { detail: { clear: true } }));
  }

  // The bar: what is selected, and the buttons for it.
  function describe(form) {
    const { columns, rows } = selection;
    const label = form.querySelector("[data-selection]");
    const actions = form.querySelector("[data-selection-actions]");
    const show = (selector, on) => { form.querySelector(selector).hidden = !on; };
    const plural = (k, one, many) => `${k} ${k === 1 ? one : many}`;
    show("[data-cut-columns]", columns.length > 0);
    show("[data-keep-columns]", columns.length > 0);
    show("[data-drop-rows]", rows.length > 0);
    show("[data-keep-rows]", rows.length > 0);
    actions.hidden = !columns.length && !rows.length;
    if (columns.length) {
      label.textContent = `${plural(columns.length, "column", "columns")} selected (${ranges(columns)})`;
      form.querySelector("[data-cut-columns]").textContent = `Delete ${plural(columns.length, "column", "columns")}`;
    } else if (rows.length) {
      label.textContent = `${plural(rows.length, "sequence", "sequences")} selected`;
      form.querySelector("[data-drop-rows]").textContent = `Remove ${plural(rows.length, "sequence", "sequences")}`;
    } else {
      label.textContent = label.dataset.idle || label.textContent;
    }
  }

  document.addEventListener("aln:select", (event) => {
    const form = editorForm();
    if (!form || !form.contains(event.target)) return;
    selection = event.detail;
    describe(form);
  });
  document.addEventListener("aln:sort", (event) => {
    const form = editorForm();
    if (!form || !form.contains(event.target)) return;
    remember(form);
    for (const f of ONE_OFF) form.elements[f].value = "";
    form.elements.sort_toggle.value = "1";
    form.requestSubmit();
  });
  document.addEventListener("click", (event) => {
    const form = editorForm();
    if (!form) return;
    const target = event.target.closest("button, a");
    if (!target || !form.parentElement.contains(target)) return;
    if (target.matches("[data-overwrite]")) {
      const save = target.closest("form");
      confirmDialog({
        title: `Overwrite “${target.dataset.name}”?`,
        body: "Its rows are replaced by the alignment as edited here, and this can't be undone. "
          + "Jobs that used it keep their own copy, as do copies saved from it.",
        action: "Overwrite",
        danger: true,
      }, () => {
        save.elements.mode.value = "overwrite";
        save.submit();
      });
    } else if (target.matches("[data-undo]")) undo();
    else if (target.matches("[data-cut-columns]")) submitEdit("cut_columns", selection.columns);
    else if (target.matches("[data-keep-columns]")) submitEdit("keep_columns", selection.columns);
    else if (target.matches("[data-drop-rows]")) submitEdit("drop_rows", selection.rows);
    else if (target.matches("[data-keep-rows]")) submitEdit("keep_rows", selection.rows);
    else if (target.matches("[data-clear]")) clearSelection();
    else if (target.matches("[data-edit]")) remember(form); // before HTMX sends it
  });
  document.addEventListener("keydown", (event) => {
    const form = editorForm();
    if (!form || event.target.closest("input, textarea, select, [contenteditable]")) return;
    if ((event.metaKey || event.ctrlKey) && !event.shiftKey && event.key.toLowerCase() === "z") {
      event.preventDefault();
      undo();
    } else if (event.key === "Delete" || event.key === "Backspace") {
      if (selection.columns.length) submitEdit("cut_columns", selection.columns);
      else if (selection.rows.length) submitEdit("drop_rows", selection.rows);
      else return;
      event.preventDefault();
    } else if (event.key === "Escape") {
      clearSelection();
    }
  });

  function ready() {
    const form = editorForm();
    if (!form) return;
    selection = { columns: [], rows: [] };
    const label = form.querySelector("[data-selection]");
    label.dataset.idle = label.textContent.trim().replace(/\s+/g, " ");
    syncUndo(form, stack(form).length);
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", ready);
  else ready();
  document.addEventListener("htmx:afterSettle", ready);
})();
