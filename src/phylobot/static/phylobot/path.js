// The mutational path page (phylobot.path_views; _path_viewer.html). The page is one GET
// form that HTMX swaps in place; this script adds what the server cannot:
//
// - The alignment viewers scroll together. Each is an alignment.js viewer on its own
//   columns; a viewer other than the traced one carries `sync` in its data, the traced
//   alignment's column to hold beside each of its own (held to an extant sequence's
//   residues, or to the residues the alignments share). Scrolling one moves the others
//   to the same place in the sequences, and side by side the rows move together too.
// - Pointing at a column boxes it, and the column(s) of the same site in the others.
// - One zoom for every viewer, in place, its value kept in the form and the address.
// - What the reader had folded or open, where they had scrolled the alignments and where
//   the page was all survive the swaps, so changing a setting never moves the page.
// - Fields that belong to another form (naming a path to save it) and the extant search
//   must not reload the viewer when they change.
(() => {
  const MIN_CELL = 0.05, MAX_CELL = 20;  // alignment.js, library_glyphs
  const cellAt = (v) => MIN_CELL * (MAX_CELL / MIN_CELL) ** (Number(v) / 100);
  const sliderAt = (cell) => 100 * Math.log(cell / MIN_CELL) / Math.log(MAX_CELL / MIN_CELL);
  const form = () => document.getElementById("path-viewer");
  const $$ = (selector, root = document) => [...root.querySelectorAll(selector)];

  // --- The viewers, in page order (the traced one first).
  const viewers = () => $$("#path-viewer .path-frame .aln-root").map((root) => {
    const wrap = root.querySelector(".aln-wrap");
    const source = root.querySelector('script[id^="aln-data"]');
    if (!wrap.pathData) wrap.pathData = JSON.parse(source.textContent);
    return { wrap, scroller: wrap.querySelector(".aln-scroller"), data: wrap.pathData };
  });
  const cellOf = (v) => Number(v.wrap.dataset.cell) || 12;
  // A viewer's column (fractional) <-> the traced alignment's column.
  const toMain = (v, column) => {
    const sync = v.data.sync;
    if (!sync || !sync.length) return column;
    const i = Math.max(0, Math.min(sync.length - 1, Math.floor(column)));
    const next = sync[Math.min(sync.length - 1, i + 1)];
    return sync[i] + (next - sync[i]) * (column - Math.floor(column));
  };
  const fromMain = (v, column) => {
    const sync = v.data.sync;
    if (!sync || !sync.length) return column;
    let lo = 0, hi = sync.length - 1;  // the first of its columns at or past `column`
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (sync[mid] < column) lo = mid + 1; else hi = mid;
    }
    return lo;
  };

  // Only the viewer being scrolled leads; the others follow without answering back.
  let leader = null, release = 0;
  const lead = (scroller) => {
    leader = scroller;
    clearTimeout(release);
    release = setTimeout(() => { leader = null; }, 120);
  };
  const place = (all, column) => all.forEach((v) => {
    v.scroller.scrollLeft = fromMain(v, column) * cellOf(v);
  });
  document.addEventListener("scroll", (event) => {
    const scroller = event.target;
    if (!scroller.classList || !scroller.classList.contains("aln-scroller")) return;
    if (!scroller.closest("#path-viewer .path-frame")) return;
    if (leader && leader !== scroller) return;  // moved by us
    lead(scroller);
    const all = viewers();
    const from = all.find((v) => v.scroller === scroller);
    if (!from) return;
    const side = Boolean(scroller.closest(".path-frames.side"));
    const column = toMain(from, scroller.scrollLeft / cellOf(from));
    for (const v of all) {
      if (v === from) continue;
      v.scroller.scrollLeft = fromMain(v, column) * cellOf(v);
      if (side) v.scroller.scrollTop = scroller.scrollTop;
    }
  }, true);

  // --- The same site in every viewer: the column under the pointer is boxed, and in each
  // other viewer the column(s) holding the same residues (`data.match`: per column, the
  // traced alignment's column, -1 where it shares no residue with it).
  document.addEventListener("aln:hover", (event) => {
    const all = viewers();
    const from = all.find((v) => v.wrap === event.target);
    if (!from) return;
    const column = event.detail.column;
    const match = from.data.match;
    const main = column < 0 ? -1 : match ? match[column] : column;
    for (const v of all) {
      let columns = [];
      if (column >= 0 && v === from) columns = [column];
      if (column >= 0 && main >= 0) {
        if (!v.data.match) columns = [main];
        else {
          const same = [];
          v.data.match.forEach((m, c) => { if (m === main) same.push(c); });
          columns = v === from ? [...new Set([column, ...same])] : same;
        }
      }
      v.wrap.dispatchEvent(new CustomEvent("aln:update", { detail: { columns } }));
    }
  });

  // --- Zoom.
  const placeSlider = () => {
    const f = form();
    const slider = f && f.querySelector(".path-zoom .zoom-slider");
    const first = viewers()[0];
    if (slider && first) slider.value = String(Math.round(sliderAt(cellOf(first))));
  };
  let queued = null;
  const zoom = (value) => {
    const idle = queued === null;
    queued = value;
    if (!idle) return;
    requestAnimationFrame(() => {
      const cell = queued;
      queued = null;
      const all = viewers();
      const first = all[0];
      const at = first ? toMain(first, first.scroller.scrollLeft / cellOf(first)) : 0;
      all.forEach((v) => v.wrap.dispatchEvent(new CustomEvent("aln:update", { detail: { cell } })));
      const text = cell === "fit" ? "fit" : cell.toFixed(2);
      const field = form().querySelector('.path-zoom input[name="zoom"]');
      if (field) field.value = text;
      const url = new URL(window.location.href);
      url.searchParams.set("zoom", text);
      history.replaceState(history.state, "", url);
      requestAnimationFrame(() => {  // "fit" and the new widths are drawn: line them up again
        if (first) lead(first.scroller);
        place(all, at);
        placeSlider();
      });
    });
  };
  document.addEventListener("input", (event) => {
    if (event.target.matches("#path-viewer .path-zoom .zoom-slider")) zoom(cellAt(event.target.value));
    if (event.target.matches("#path-viewer .extant-search")) {
      const wanted = event.target.value.trim().toLowerCase();
      $$(".extant-list label", event.target.closest("fieldset")).forEach((label) => {
        label.hidden = Boolean(wanted) && !label.textContent.toLowerCase().includes(wanted);
      });
    }
  });
  document.addEventListener("click", (event) => {
    if (event.target.closest("#path-viewer .path-zoom .zoom-fit")) zoom("fit");
  });

  // --- Changes that must not reload the viewer: the zoom slider and the search (no
  // server value), and the fields of the page-level forms.
  const quiet = (target) => target.matches && target.matches(
    "#path-viewer .zoom-slider, #path-viewer .extant-search, #path-viewer [form]");
  document.addEventListener("change", (event) => {
    if (quiet(event.target)) event.stopPropagation();
  }, true);
  document.addEventListener("keydown", (event) => {
    if (event.key !== "Enter" || !quiet(event.target)) return;
    event.preventDefault();  // Enter in a name field saves; elsewhere it does nothing
    const owner = event.target.getAttribute("form");
    if (owner) document.getElementById(owner).requestSubmit();
  }, true);
  // --- Across a swap (and a reload): which sections are folded, the place in the
  // alignments, and the place on the page, so that changing a setting never moves it.
  // The folds are read from the page just before each request and put back the moment
  // the new form is in (afterSwap, before anything is painted). They are never taken
  // from "toggle" events during a swap: a <details> that arrives open fires one too,
  // which would overwrite what the reader chose. They are kept for the tab as well
  // (sessionStorage), so "Trace here", a saved path or a reload opens the page as it was
  // left. The page is then held where it was: the control just used stays put.
  const STORE = "phylobot.path.sections";
  const sections = () => $$("#path-viewer details[id]");
  const snapshot = () => Object.fromEntries(sections().map((d) => [d.id, d.open]));
  const keepForTab = (states) => {
    const lasting = { ...states };
    delete lasting["path-extant"];  // a menu: open only while it is being used
    try { sessionStorage.setItem(STORE, JSON.stringify(lasting)); } catch { /* private mode */ }
  };
  const restore = (states) => {
    if (!states) return;
    for (const d of sections()) if (d.id in states) d.open = states[d.id];
  };
  let kept = null, folds = null, hold = null;
  document.addEventListener("toggle", (event) => {  // folded by hand: for the next visit too
    if (!folds && event.target.matches && event.target.matches("#path-viewer details[id]")) keepForTab(snapshot());
  }, true);
  // The same control after the swap: by name, and by value where several share a name.
  const same = (field) => {
    const named = $$(`#path-viewer [name="${CSS.escape(field.name)}"]`);
    return (field.byValue ? named.find((el) => el.value === field.value) : named[0]) || null;
  };
  const holdPlace = () => {
    if (!hold) return;
    const now = hold.field ? same(hold.field) : null;
    if (now && now.getClientRects().length) window.scrollBy(0, now.getBoundingClientRect().top - hold.top);
    else window.scrollTo(window.scrollX, hold.y);
  };
  document.addEventListener("htmx:beforeRequest", (event) => {
    if (event.target.id !== "path-viewer") return;
    const first = viewers()[0];
    kept = first ? toMain(first, first.scroller.scrollLeft / cellOf(first)) : null;
    folds = snapshot();
    keepForTab(folds);
    const used = document.activeElement;
    hold = { y: window.scrollY };
    if (used && used.name && used.closest("#path-viewer")) {
      hold.field = { name: used.name, value: used.value, byValue: ["radio", "checkbox", "submit"].includes(used.type) };
      hold.top = used.getBoundingClientRect().top;
    }
  });
  document.addEventListener("htmx:afterSwap", (event) => {
    if (event.target.id !== "path-viewer") return;
    restore(folds);
    holdPlace();
  });
  document.addEventListener("htmx:afterSettle", (event) => {
    if (event.target.id !== "path-viewer") return;
    restore(folds);
    requestAnimationFrame(() => {  // after alignment.js has drawn (and sized) the new viewers
      const all = viewers();
      if (kept !== null && all.length) {
        lead(all[0].scroller);
        place(all, kept);
      }
      holdPlace();
      kept = null;
      hold = null;
      setTimeout(() => { folds = null; }, 0);  // once the arriving sections' own toggles have fired
      placeSlider();
    });
  });
  const start = () => {
    try { restore(JSON.parse(sessionStorage.getItem(STORE) || "null")); } catch { /* nothing kept */ }
    requestAnimationFrame(placeSlider);
  };
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();
})();
