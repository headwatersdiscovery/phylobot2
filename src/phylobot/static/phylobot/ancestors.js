// Ancestor viewer: hovering the probability trace shows the site under the cursor (its
// number in the ancestor, alignment column, and top three states); clicking jumps to
// that site's row in the table (each plot's PNG/SVG downloads are plots.js). Delegated,
// so it keeps working after HTMX swaps.
(() => {
  const parsed = new WeakMap();
  function sites(svg) {
    const source = document.getElementById(svg.dataset.sites);
    if (!parsed.has(source)) parsed.set(source, JSON.parse(source.textContent));
    return parsed.get(source);
  }
  function siteAt(svg, event) {
    const box = svg.getBoundingClientRect();
    const x = ((event.clientX - box.left) / box.width) * svg.viewBox.baseVal.width;
    const list = sites(svg);
    const fraction = (x - Number(svg.dataset.left)) / Number(svg.dataset.width);
    const i = Math.round(fraction * (list.length - 1));
    return i >= 0 && i < list.length ? { i, site: list[i], x } : null;
  }
  function hide(figure) {
    figure.querySelector(".pp-tip").hidden = true;
    figure.querySelector(".guide").setAttribute("visibility", "hidden");
  }
  document.addEventListener("mousemove", (event) => {
    const svg = event.target.closest && event.target.closest("svg.pp-trace");
    document.querySelectorAll("figure.trace").forEach((f) => { if (!f.contains(svg)) hide(f); });
    if (!svg) return;
    const figure = svg.closest("figure");
    const hit = siteAt(svg, event);
    if (!hit) { hide(figure); return; }
    const [number, column, ...states] = hit.site;
    const guide = figure.querySelector(".guide");
    const list = sites(svg);
    const left = Number(svg.dataset.left), width = Number(svg.dataset.width);
    const gx = left + (list.length > 1 ? hit.i / (list.length - 1) : 0.5) * width;
    guide.setAttribute("x1", gx);
    guide.setAttribute("x2", gx);
    guide.setAttribute("visibility", "visible");
    const tip = figure.querySelector(".pp-tip");
    const pairs = [];
    for (let k = 0; k < states.length; k += 2) {
      if (states[k]) pairs.push(`${states[k]} ${states[k + 1].toFixed(2)}`);
    }
    // (the branch-site plot is drawn along the alignment: its sites are columns)
    const where = number === column ? `Alignment column ${column}` : `Site ${number} · alignment column ${column}`;
    tip.textContent = `${where}\n${pairs.join(" · ")}`;
    tip.hidden = false;
    const box = figure.getBoundingClientRect();
    const px = event.clientX - box.left;
    tip.style.left = `${Math.min(px + 12, box.width - tip.offsetWidth - 4)}px`;
    tip.style.top = `${event.clientY - box.top + 12}px`;
  });
  document.addEventListener("click", (event) => {
    const svg = event.target.closest && event.target.closest("svg.pp-trace");
    if (!svg) return;
    const hit = siteAt(svg, event);
    const row = hit && document.getElementById(`site-${hit.site[1]}`);
    if (!row) return;
    row.scrollIntoView({ block: "center", behavior: "smooth" });
    row.classList.remove("flash");
    void row.offsetWidth;  // restart the highlight animation
    row.classList.add("flash");
  });
})();

// AltAll threshold: the slider and the number field stay in step, and the count of
// ambiguous sites is recomputed in the browser as the slider moves (the server sends
// every site's best alternative, sorted, so this is a binary search, not a scan). The
// slider has no name, so the number field is still what the form submits; releasing the
// slider fires the form's change trigger, which refreshes the download links.
(() => {
  const parsed = new WeakMap();
  function probabilities(panel) {
    const source = panel.querySelector("#altall-pps");
    if (!source) return [];
    if (!parsed.has(source)) parsed.set(source, JSON.parse(source.textContent));
    return parsed.get(source);
  }
  // `pps` is sorted high to low: the first index below `threshold` is the count at it.
  function countAtLeast(pps, threshold) {
    let low = 0;
    let high = pps.length;
    while (low < high) {
      const mid = (low + high) >> 1;
      if (pps[mid] >= threshold) low = mid + 1;
      else high = mid;
    }
    return low;
  }
  function show(panel, value) {
    const count = panel.querySelector(".altall-count");
    if (!count) return;
    const n = countAtLeast(probabilities(panel), value);
    const unit = count.dataset.unit || "site";
    count.querySelector("strong").textContent = String(n);
    count.lastChild.textContent = ` ambiguous ${unit}${n === 1 ? "" : "s"} of ${count.dataset.total}`;
  }
  document.addEventListener("input", (event) => {
    const panel = event.target.closest && event.target.closest(".altall");
    if (!panel) return;
    const slider = panel.querySelector(".altall-slider");
    const field = panel.querySelector('input[name="altall"]');
    if (!slider || !field) return;
    const value = Number(event.target.value);
    if (!Number.isFinite(value)) return;
    if (event.target === slider) field.value = slider.value;
    else slider.value = field.value;
    show(panel, value);
  });
})();

// The AltAll threshold drawn on the support plots: a vertical line on the histogram of
// the alternatives' probabilities (sites to its right are the ones AltAll switches) and
// a level on the top-three trace (where the 2nd-best line rises above it). Both are SVG
// attributes, so the PNG/SVG downloads carry the mark exactly as shown. Delegated and
// re-applied after HTMX swaps, like the rest of this file.
(() => {
  function place(value) {
    const mark = document.querySelector(".altall-histogram .altall-mark");
    if (mark) {
      const svg = mark.closest("svg");
      const left = Number(svg.dataset.left);
      const width = Number(svg.dataset.width);
      if (Number.isFinite(left) && Number.isFinite(width)) {
        const x = (left + value * width).toFixed(1);
        mark.setAttribute("x1", x);
        mark.setAttribute("x2", x);
      }
    }
    const level = document.querySelector(".pp-trace .altall-level");
    if (level) {
      const svg = level.closest("svg");
      const top = Number(svg.dataset.top);
      const height = Number(svg.dataset.plotHeight);
      if (Number.isFinite(top) && Number.isFinite(height)) {
        const y = (top + (1 - value) * height).toFixed(1);
        level.setAttribute("y1", y);
        level.setAttribute("y2", y);
      }
    }
  }
  function current() {
    const field = document.querySelector('.altall input[name="altall"]');
    return field ? Number(field.value) : NaN;
  }
  const show = () => {
    const value = current();
    if (Number.isFinite(value)) place(value);
  };
  document.addEventListener("input", (event) => {
    if (!event.target.closest || !event.target.closest(".altall")) return;
    const value = Number(event.target.value);
    if (Number.isFinite(value)) place(value);
  });
  document.addEventListener("DOMContentLoaded", show);
  document.addEventListener("htmx:afterSettle", show);
  show();
})();

// Collapsed sections stay collapsed. The viewer's form — the sections with it — is
// replaced on every menu change, so without this a reader who folded "Sites" away would
// have it unfold again the moment they picked another node. Remembered per section in
// this browser, as the tree viewer's caption and the job page's steps panel are, and
// wrapped in try/catch because private mode makes localStorage throw.
(() => {
  const KEY = "phylobot.ancestorSections";
  function closed() {
    try {
      return new Set(JSON.parse(localStorage.getItem(KEY) || "[]"));
    } catch {
      return new Set();
    }
  }
  function remember(set) {
    try { localStorage.setItem(KEY, JSON.stringify([...set])); } catch { /* private mode */ }
  }
  const apply = () => {
    const shut = closed();
    for (const section of document.querySelectorAll("details.viewer-section[id]")) {
      // A section the page itself starts folded (data-folded: it holds nothing yet) stays
      // folded until the reader opens it, remembered as "+id".
      section.open = section.hasAttribute("data-folded") ? shut.has(`+${section.id}`) : !shut.has(section.id);
    }
  };
  document.addEventListener("toggle", (event) => {
    const section = event.target;
    if (!section.matches || !section.matches("details.viewer-section[id]")) return;
    const shut = closed();
    if (section.hasAttribute("data-folded")) {
      if (section.open) shut.add(`+${section.id}`);
      else shut.delete(`+${section.id}`);
    } else if (section.open) shut.delete(section.id);
    else shut.add(section.id);
    remember(shut);
  }, true);
  document.addEventListener("DOMContentLoaded", apply);
  document.addEventListener("htmx:afterSettle", apply);
  apply();
})();

// Naming and describing the ancestor: Save turns amber once the name or the description
// differs from what was saved, so it is clear the change still needs a click.
(() => {
  document.addEventListener("input", (event) => {
    const box = event.target.closest && event.target.closest(".ancestor-name");
    if (!box) return;
    const changed = [...box.querySelectorAll("input[name=ancestor_name], textarea[name=description]")]
      .some((field) => field.value !== field.defaultValue);
    const save = box.querySelector("[data-name-save]");
    if (save) save.classList.toggle("dirty", changed);
  });
})();

// "Copy sequence": the ancestor's sequence as plain letters (its residues are links, one
// per site, so the text is taken without the markup's whitespace), and a brief "Copied".
document.addEventListener("click", async (event) => {
  const button = event.target.closest && event.target.closest("[data-copy-sequence]");
  if (!button) return;
  const sequence = button.closest("details, section").querySelector("[data-sequence]");
  const label = button.querySelector("span") || button;
  const was = label.textContent;
  try {
    await navigator.clipboard.writeText(sequence.textContent.replace(/\s+/g, ""));
    label.textContent = "Copied";
  } catch {
    getSelection().selectAllChildren(sequence);
    label.textContent = "Selected: press Ctrl+C";
  }
  setTimeout(() => { label.textContent = was; }, 1600);
});
