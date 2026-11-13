// The structure page's 3D views (structure_views, _structure_viewer.html), drawn with
// 3Dmol.js. The server decides every colour (structures.colours) and sends them in
// #structure-data: for each panel one colour per residue, keyed "chain:number" as
// structures.pdb.key writes them. This script only paints, and keeps the camera where
// the reader left it.
//
// A panel's stage (#structure-stage-<n>) is kept across the form's HTMX swaps
// (hx-preserve), so a new ancestor or colour scheme repaints the structure that is
// already loaded. Everything outside the stages is replaced with each swap: its buttons
// are reached by delegation.
//
// Panels of different structures share one frame: the server sends each a rigid
// transform onto the first panel's structure (structures.superpose), applied to the
// atoms as they are loaded, and the viewers are then linked so one camera serves all.
(function () {
  "use strict";

  const OTHER = "#b9bec6"; // chains and residues the page has no colour for
  const BACKBONE = new Set(["N", "C", "O", "OXT"]);
  const texts = new Map(); // file address -> a promise of its text
  let stages = []; // per panel: {index, el, viewer, loaded, alpha, panel}
  let data = null;
  let picked = null; // {panel, key}

  function keyOf(atom) {
    const chain = (atom.chain || "").trim() || "_";
    return chain + ":" + atom.resi + (atom.icode || "").trim();
  }

  function stageElement(index) {
    return document.getElementById("structure-stage-" + index);
  }

  function status(index, text) {
    const stage = stageElement(index);
    const el = stage && stage.querySelector("[data-structure-status]");
    if (el) {
      el.textContent = text;
      el.hidden = !text;
    }
  }

  function readData() {
    const el = document.getElementById("structure-data");
    if (!el) return null;
    try {
      return JSON.parse(el.textContent);
    } catch (e) {
      return null;
    }
  }

  function hex(colour) {
    // "rgb(12, 34, 56)" (a computed style) as "#0c2238"
    const m = /rgba?\((\d+),\s*(\d+),\s*(\d+)/.exec(colour || "");
    if (!m) return colour || "#ffffff";
    return "#" + [m[1], m[2], m[3]].map((v) => Number(v).toString(16).padStart(2, "0")).join("");
  }

  function background(index) {
    if (data.background === "white") return "#ffffff";
    if (data.background === "black") return "#000000";
    return hex(getComputedStyle(stageElement(index)).backgroundColor);
  }

  function elementColour(atom, carbon) {
    if (atom.elem === "C") return carbon;
    const table = ($3Dmol.elementColors && $3Dmol.elementColors.rasmol) || {};
    const found = table[atom.elem];
    return found === undefined ? 0x999999 : found;
  }

  function label(viewer, atom, text) {
    viewer.addLabel(text, {
      position: { x: atom.x, y: atom.y, z: atom.z },
      fontSize: 12,
      fontColor: "white",
      backgroundColor: "black",
      backgroundOpacity: 0.7,
      inFront: true,
    });
  }

  function linked(stage) {
    return data.link === "on" && stage.panel.linked && stages.length > 1;
  }

  function paint(stage) {
    const viewer = stage.viewer;
    const panel = stage.panel;
    if (!viewer || !panel || stage.loaded === "") return;
    const colourOf = new Map();
    panel.keys.forEach((key, i) => colourOf.set(key, panel.colours[i]));
    const marked = new Set(panel.marked.map((i) => panel.keys[i]));
    const known = (atom) => colourOf.has(keyOf(atom));
    const base = (atom) => colourOf.get(keyOf(atom)) || OTHER;
    const sticks = (atom) => elementColour(atom, base(atom));
    const protein = { predicate: known };

    viewer.removeAllSurfaces();
    viewer.removeAllLabels();
    viewer.removeAllShapes();
    viewer.setStyle({}, {});
    // What the page has no residue for (nucleic acids, short peptides): plain grey.
    viewer.setStyle({ hetflag: false, predicate: (a) => !known(a) }, { cartoon: { color: OTHER } });
    if (data.style === "stick") {
      viewer.setStyle(protein, { stick: { colorfunc: sticks, radius: 0.2 } });
    } else if (data.style === "sphere") {
      viewer.setStyle(protein, { sphere: { colorfunc: base } });
    } else {
      const cartoon = { colorfunc: base };
      if (data.style === "trace") cartoon.style = "trace";
      viewer.setStyle(protein, { cartoon: cartoon });
      if (data.sidechains !== "none") {
        const every = data.sidechains === "all";
        viewer.addStyle(
          {
            predicate: (a) =>
              known(a) && !BACKBONE.has(a.atom) && (every || marked.has(keyOf(a))),
          },
          { stick: { colorfunc: sticks, radius: 0.22 } }
        );
      }
    }
    if (data.ligands === "on") {
      viewer.addStyle(
        { hetflag: true, predicate: (a) => !known(a) && a.resn !== "HOH" && a.resn !== "DOD" },
        { stick: { radius: 0.2 }, sphere: { scale: 0.28 } }
      );
    }
    if (data.surface === "on") {
      viewer.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.85, colorfunc: base }, protein, protein);
    }
    // Residues the structure has no place for: a marker between the two it has.
    (panel.between || []).forEach((gap) => {
      const a = stage.alpha.get(gap.after);
      const b = stage.alpha.get(gap.before);
      if (!a || !b) return;
      viewer.addSphere({
        center: { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, z: (a.z + b.z) / 2 },
        radius: 1.1,
        color: data.betweenColour,
      });
    });
    if (data.labels === "on") {
      Object.keys(panel.labels_text || {}).forEach((i) => {
        const atom = stage.alpha.get(panel.keys[Number(i)]);
        if (atom) label(viewer, atom, panel.labels_text[i]);
      });
    }
    viewer.setBackgroundColor(background(stage.index));
    // Linked panels follow the first: it alone spins, and the others are set to its view.
    const follows = linked(stage) && stage.index > 0 && stages[0] && linked(stages[0]);
    viewer.spin(data.spin === "on" && !follows ? "y" : false);
    paintPicked(stage);
    viewer.render();
  }

  // The keys of `stage` that are the picked residue: itself in its own panel, and in the
  // others the residues of the same alignment column (the focus chain's first).
  function pickedKeys(stage) {
    if (!picked) return [];
    if (stage.index === picked.panel) return [picked.key];
    const from = stages[picked.panel];
    if (!from || !from.panel) return [];
    const column = from.panel.cols[from.panel.keys.indexOf(picked.key)];
    if (column === null || column === undefined) return [];
    const keys = [];
    stage.panel.cols.forEach((c, i) => {
      if (c === column) keys.push(stage.panel.keys[i]);
    });
    const focus = keys.filter((key) => key.split(":")[0] === stage.panel.focus);
    return focus.length ? focus : keys;
  }

  function tipOf(stage, key) {
    const i = stage.panel ? stage.panel.keys.indexOf(key) : -1;
    return i < 0 ? "" : stage.panel.tips[i];
  }

  function paintPicked(stage) {
    const keys = new Set(pickedKeys(stage));
    if (!keys.size) return;
    stage.viewer.addStyle(
      { predicate: (a) => keys.has(keyOf(a)) },
      { stick: { radius: 0.3, colorscheme: "default" } }
    );
    keys.forEach((key) => {
      const atom = stage.alpha.get(key);
      // Not over a change's own label ("T222K"), where those are drawn.
      const at = String(stage.panel.keys.indexOf(key));
      const named = data.labels === "on" && (stage.panel.labels_text || {})[at];
      if (atom && !named) label(stage.viewer, atom, atom.resn + " " + key.split(":")[1]);
    });
  }

  // The picked residue outside the 3D views: outlined in each panel's strip, and
  // described under the panels.
  function showPicked() {
    document.querySelectorAll(".sx.picked").forEach((el) => el.classList.remove("picked"));
    const line = document.querySelector("[data-structure-picked]");
    const lines = [];
    stages.forEach((stage) => {
      const strip = document.querySelector('[data-structure-strip="' + stage.index + '"]');
      pickedKeys(stage).forEach((key, n) => {
        const text = tipOf(stage, key);
        if (text && n === 0) {
          lines.push(stages.length > 1 ? "Panel " + (stage.index + 1) + ": " + text : text);
        }
        const cell = strip && strip.querySelector('.sx[data-key="' + CSS.escape(key) + '"]');
        if (!cell) return;
        cell.classList.add("picked");
        const scroller = cell.closest(".sx-scroll");
        const left = cell.offsetLeft - scroller.clientWidth / 2;
        scroller.scrollTo({ left: Math.max(0, left), behavior: "smooth" });
      });
    });
    if (line) line.textContent = lines.join("\n");
  }

  // The scroll lock: the strips are held at the same alignment column. The column at the
  // left edge of the strip being scrolled is found, and every other strip is put where
  // its own residue of that column (or the next it has) is at its left edge.
  function followScroll(source) {
    if (!data || data.scroll !== "on") return;
    if (source.scrollLeft === source.syncedTo) return; // moved by this function: not a lead
    const cells = Array.from(source.querySelectorAll(".sx[data-col]"));
    const first = cells.find((cell) => cell.offsetLeft + cell.offsetWidth > source.scrollLeft);
    if (!first) return;
    const column = Number(first.dataset.col);
    const within = source.scrollLeft - first.offsetLeft;
    document.querySelectorAll("[data-structure-strip] .sx-scroll").forEach((other) => {
      if (other === source) return;
      const twin = Array.from(other.querySelectorAll(".sx[data-col]")).find(
        (cell) => Number(cell.dataset.col) >= column
      );
      if (!twin) return;
      const most = other.scrollWidth - other.clientWidth;
      const left = Math.max(0, Math.min(most, twin.offsetLeft + within));
      other.syncedTo = left;
      other.scrollLeft = left;
      other.syncedTo = other.scrollLeft; // as the browser rounded it
    });
  }

  function fullScreen() {
    return document.documentElement.classList.contains("structure-full");
  }

  // Full screen is the whole page's, with the Structures section laid over it (site.css,
  // html.structure-full): the section is replaced by every HTMX swap, and a full-screen
  // element that is replaced drops out of full screen, so the controls could not come.
  function setFullScreen(on) {
    document.documentElement.classList.toggle("structure-full", on);
    const button = document.querySelector("[data-structure-full]");
    if (button) button.textContent = on ? "Exit full screen" : "Full screen";
    const view = document.getElementById("view");
    if (on && view) view.open = true;
    stages.forEach((stage) => stage.viewer && stage.viewer.resize());
  }

  function pick(panel, key, zoom) {
    const same = picked && picked.panel === panel && picked.key === key;
    picked = same && !zoom ? null : { panel: panel, key: key };
    stages.forEach(paint);
    showPicked();
    if (!zoom || !picked) return;
    // The residue with its neighbours along the chain, not the residue alone. Linked
    // panels are zoomed by the one that leads them; the others each by themselves.
    stages.forEach((stage) => {
      if (linked(stage) && stage.index !== picked.panel && linked(stages[picked.panel])) return;
      const keys = pickedKeys(stage);
      if (!keys.length) return;
      const at = stage.panel.keys.indexOf(keys[0]);
      const near = new Set(stage.panel.keys.slice(Math.max(0, at - 7), at + 8));
      stage.viewer.zoomTo({ predicate: (a) => near.has(keyOf(a)) }, 600);
    });
  }

  function pickColumn(column) {
    for (const stage of stages) {
      if (!stage.panel || !stage.panel.cols) continue;
      const found = [];
      stage.panel.cols.forEach((c, i) => {
        if (c === column) found.push(stage.panel.keys[i]);
      });
      const focus = found.filter((key) => key.split(":")[0] === stage.panel.focus);
      const key = focus[0] || found[0];
      if (key) {
        pick(stage.index, key, true);
        return true;
      }
    }
    return false;
  }

  function showTip(stage, atom) {
    const tip = stage.el.parentElement.querySelector("[data-structure-tip]");
    if (!tip) return;
    tip.textContent =
      tipOf(stage, keyOf(atom)) || atom.resn + " " + atom.resi + " · chain " + atom.chain;
    tip.hidden = false;
    const box = stage.el.getBoundingClientRect();
    tip.style.left = Math.max(4, Math.min(stage.pointer.x - box.left + 14, box.width - 260)) + "px";
    tip.style.top = stage.pointer.y - box.top + 14 + "px";
  }

  function hideTip(stage) {
    const tip = stage.el.parentElement.querySelector("[data-structure-tip]");
    if (tip) tip.hidden = true;
  }

  // Every viewer that shares the frame shows what the others show.
  function link() {
    const together = stages.filter((stage) => stage.viewer && stage.loaded && linked(stage));
    stages.forEach((stage) => {
      if (stage.viewer) stage.viewer.linkedViewers = [];
    });
    if (together.length < 2) return;
    together.forEach((stage) => {
      together.forEach((other) => {
        if (other !== stage) stage.viewer.linkViewer(other.viewer);
      });
    });
    const lead = together[0];
    together.slice(1).forEach((stage) => stage.viewer.setView(lead.viewer.getView(), true));
  }

  function load(stage, text, wanted) {
    const viewer = stage.viewer;
    const panel = stage.panel;
    viewer.linkedViewers = [];
    viewer.removeAllModels();
    viewer.addModel(text, panel.format === "cif" ? "cif" : "pdb");
    const m = panel.matrix;
    if (m) {
      // Onto the first panel's structure: x' = R x + t, before anything is drawn.
      viewer.selectedAtoms({}).forEach((atom) => {
        const x = atom.x;
        const y = atom.y;
        const z = atom.z;
        atom.x = m[0] * x + m[1] * y + m[2] * z + m[9];
        atom.y = m[3] * x + m[4] * y + m[5] * z + m[10];
        atom.z = m[6] * x + m[7] * y + m[8] * z + m[11];
      });
    }
    stage.alpha = new Map();
    viewer.selectedAtoms({ atom: "CA" }).forEach((atom) => {
      const key = keyOf(atom);
      if (!stage.alpha.has(key)) stage.alpha.set(key, atom);
    });
    viewer.setHoverDuration(60);
    viewer.setHoverable({}, true, (atom) => showTip(stage, atom), () => hideTip(stage));
    viewer.setClickable({}, true, (atom) => pick(stage.index, keyOf(atom), false));
    stage.loaded = wanted;
    if (picked && picked.panel === stage.index) picked = null;
    paint(stage);
    // A panel that follows the first takes its view; any other frames its own structure.
    const lead = stages[0];
    const follows = stage.index > 0 && linked(stage) && lead && lead.loaded && linked(lead);
    if (!follows) viewer.zoomTo();
    viewer.render();
    link();
  }

  function fetchText(file) {
    if (!texts.has(file)) {
      const asked = fetch(file, { credentials: "same-origin" }).then((response) => {
        if (!response.ok) throw new Error(String(response.status));
        return response.text();
      });
      asked.catch(() => texts.delete(file));
      texts.set(file, asked);
    }
    return texts.get(file);
  }

  function release(stage) {
    // A dropped panel's WebGL context is given back at once: a browser allows few.
    try {
      const canvas = stage.el.querySelector("canvas");
      const gl = canvas && (canvas.getContext("webgl2") || canvas.getContext("webgl"));
      const lose = gl && gl.getExtension("WEBGL_lose_context");
      if (lose) lose.loseContext();
    } catch (e) {
      /* nothing to release */
    }
  }

  function refresh() {
    data = readData();
    if (!data || !data.panels) {
      stages.forEach(release);
      stages = [];
      return;
    }
    if (typeof window.$3Dmol === "undefined") {
      data.panels.forEach((_, i) => status(i, "The 3D viewer could not start in this browser."));
      return;
    }
    const before = stages;
    stages = data.panels.map((panel, index) => {
      const stageEl = stageElement(index);
      const el = stageEl && stageEl.querySelector("[data-structure-canvas]");
      let stage = before.find((s) => s.el === el);
      if (!el) return { index: index, el: null, viewer: null, loaded: "", panel: panel };
      if (!stage) {
        stage = { el: el, viewer: null, loaded: "", alpha: new Map(), pointer: { x: 0, y: 0 } };
        try {
          stage.viewer = $3Dmol.createViewer(el, { backgroundColor: "white", antialias: true });
        } catch (e) {
          stage.viewer = null;
        }
        if (stage.viewer) {
          const viewer = stage.viewer;
          // A resized viewer keeps its own view: the panels are put in step again.
          new ResizeObserver(() => {
            viewer.resize();
            link();
          }).observe(el);
          el.addEventListener("mousemove", (e) => (stage.pointer = { x: e.clientX, y: e.clientY }));
          el.addEventListener("mouseleave", () => hideTip(stage));
        }
      }
      stage.index = index;
      stage.panel = panel;
      return stage;
    });
    before.filter((old) => !stages.includes(old) && old.el).forEach(release);
    if (picked && !stages[picked.panel]) picked = null;

    stages.forEach((stage) => {
      if (!stage.el) return;
      if (!stage.viewer) {
        status(stage.index, "The 3D viewer needs WebGL, which this browser does not offer.");
        return;
      }
      // The same file under another transform is another thing to draw.
      const wanted = stage.panel.file + "|" + JSON.stringify(stage.panel.matrix);
      if (stage.loaded === wanted) {
        paint(stage);
        return;
      }
      stage.wanted = wanted;
      status(stage.index, "Loading the structure…");
      fetchText(stage.panel.file)
        .then((text) => {
          if (stage.wanted !== wanted || !stages.includes(stage)) return; // chosen again meanwhile
          load(stage, text, wanted);
          status(stage.index, "");
        })
        .catch(() => status(stage.index, "The structure's file could not be loaded."));
    });
    link();
    showPicked();
    if (fullScreen()) setFullScreen(true); // the button is new with every swap
  }

  // One picture of every panel, side by side, each under its title.
  function saveImage() {
    const shown = stages.filter((stage) => stage.viewer && stage.loaded);
    if (!shown.length) return;
    const loaded = shown.map(
      (stage) =>
        new Promise((resolve) => {
          const image = new Image();
          image.onload = () => resolve({ image: image, title: stage.panel.title || "" });
          image.onerror = () => resolve(null);
          image.src = stage.viewer.pngURI();
        })
    );
    Promise.all(loaded).then((found) => {
      const pictures = found.filter(Boolean);
      if (!pictures.length) return;
      const band = 34;
      const height = Math.max(...pictures.map((p) => p.image.height));
      const sheet = document.createElement("canvas");
      sheet.width = pictures.reduce((sum, p) => sum + p.image.width, 0);
      sheet.height = height + band;
      const pen = sheet.getContext("2d");
      const dark = data.background === "black";
      pen.fillStyle = dark ? "#000000" : "#ffffff";
      pen.fillRect(0, 0, sheet.width, sheet.height);
      pen.fillStyle = dark ? "#ffffff" : "#111111";
      pen.font = "16px system-ui, sans-serif";
      let x = 0;
      pictures.forEach((p) => {
        pen.fillText(p.title, x + 10, 22, p.image.width - 20);
        pen.drawImage(p.image, x, band);
        x += p.image.width;
      });
      const anchor = document.createElement("a");
      anchor.href = sheet.toDataURL("image/png");
      anchor.download = "structure.png";
      anchor.click();
    });
  }

  function copy(button) {
    const text = button.dataset.structureCopy;
    const said = button.textContent;
    const done = () => {
      button.textContent = "copied";
      setTimeout(() => (button.textContent = said), 1200);
    };
    if (navigator.clipboard) navigator.clipboard.writeText(text).then(done, () => {});
  }

  // The lookup's results: how many structures are ticked, and whether any is.
  function countTicked(form) {
    const ticked = form.querySelectorAll('input[name="pick"]:checked').length;
    const button = form.querySelector("[data-structure-attach]");
    const note = form.querySelector("[data-structure-ticked]");
    if (button) {
      button.disabled = !ticked;
      button.textContent = ticked
        ? "Attach " + ticked + " structure" + (ticked === 1 ? "" : "s")
        : "Attach the ticked structures";
    }
    if (note) note.hidden = ticked > 0;
  }

  document.addEventListener("change", (event) => {
    const form = event.target.closest && event.target.closest("[data-structure-found]");
    if (form) countTicked(form);
  });

  document.addEventListener("click", (event) => {
    const target = event.target instanceof Element ? event.target : null;
    if (!target) return;
    const cell = target.closest(".sx[data-key]");
    const go = target.closest("[data-structure-col]");
    const copied = target.closest("[data-structure-copy]");
    const suggested = target.closest("[data-structure-suggest]");
    if (cell) {
      const strip = cell.closest("[data-structure-strip]");
      pick(Number(strip.dataset.structureStrip), cell.dataset.key, true);
    } else if (go) {
      const panels = document.querySelector("[data-structure-panels]");
      if (pickColumn(Number(go.dataset.structureCol)) && panels) {
        panels.scrollIntoView({ block: "center", behavior: "smooth" });
      }
    } else if (copied) {
      copy(copied);
    } else if (suggested) {
      const box = suggested.closest("form").querySelector('input[name="identifier"]');
      const have = box.value.split(/[\s,;]+/).filter(Boolean);
      if (!have.includes(suggested.dataset.structureSuggest)) {
        box.value = have.concat(suggested.dataset.structureSuggest).join(", ");
      }
    } else if (target.closest("[data-structure-reset]")) {
      picked = null;
      stages.forEach(paint);
      showPicked();
      // The first panel frames its structure and the linked ones follow; the rest each
      // frame their own.
      stages.forEach((stage) => {
        const follows = stage.index > 0 && linked(stage) && stages[0] && linked(stages[0]);
        if (stage.viewer && stage.loaded && !follows) stage.viewer.zoomTo({}, 400);
      });
    } else if (target.closest("[data-structure-image]")) {
      saveImage();
    } else if (target.closest("[data-structure-full]")) {
      const root = document.documentElement;
      if (fullScreen()) {
        if (document.fullscreenElement) document.exitFullscreen();
        setFullScreen(false);
      } else if (root.requestFullscreen) {
        root.requestFullscreen().then(() => setFullScreen(true), () => {});
      }
    }
  });
  document.addEventListener("fullscreenchange", () => {
    if (!document.fullscreenElement) setFullScreen(false); // left with Escape
    stages.forEach((stage) => stage.viewer && stage.viewer.resize());
  });
  // Scroll events do not bubble: they are caught on their way down.
  document.addEventListener(
    "scroll",
    (event) => {
      const el = event.target;
      if (el instanceof Element && el.matches("[data-structure-strip] .sx-scroll")) followScroll(el);
    },
    true
  );
  document.addEventListener("DOMContentLoaded", refresh);
  document.addEventListener("htmx:afterSettle", refresh);
})();
