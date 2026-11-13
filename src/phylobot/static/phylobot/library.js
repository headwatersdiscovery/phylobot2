// The library page (phylobot.library_views): its alignment viewers (alignment.js), the
// site maps beneath each row, the Controls panel, bulk actions and notices.
//
// Threshold, palette and zoom act in place and keep their values in the address, so the
// page never jumps while they are used. Changing the alignment or the extant sequences
// needs the server; those reloads keep the open sequences and the scroll position.
(() => {
  const $ = (selector, root = document) => root.querySelector(selector);
  const $$ = (selector, root = document) => [...root.querySelectorAll(selector)];
  const page = $(".library-page");
  if (!page) return;
  const setParam = (key, value) => {
    const url = new URL(location.href);
    url.searchParams.set(key, value);
    history.replaceState(null, "", url);
  };

  // Notices appear as a toast along the bottom: an action's result is read where the
  // user is working, and the page never scrolls to show it.
  let timer = null;
  const toast = (text, kind = "error") => {
    const box = $("#library-toast");
    box.className = `toast ${kind}`;
    box.textContent = text;
    box.hidden = false;
    clearTimeout(timer);
    timer = setTimeout(() => { box.hidden = true; }, 6000);
  };
  // Django's messages (after an action) are shown the same way.
  $$(".library-messages .message").forEach((message) => {
    toast(message.textContent.trim(), message.classList.contains("error") ? "error" : "success");
    message.remove();
  });

  // A reload (an action, or another alignment) comes back to the same place: the same
  // sequences open, the same scroll position.
  const scrollKey = `scroll:${location.pathname}`;
  const keepPlace = (form) => {
    if (form && form.method.toLowerCase() === "get") {
      $$("details.site-map-details[open]").forEach((details) => {
        const keep = document.createElement("input");
        Object.assign(keep, { type: "hidden", name: "open", value: details.dataset.target });
        form.appendChild(keep);
      });
    }
    try { sessionStorage.setItem(scrollKey, String(window.scrollY)); } catch {}
  };
  try {
    const saved = sessionStorage.getItem(scrollKey);
    if (saved !== null) {
      sessionStorage.removeItem(scrollKey);
      requestAnimationFrame(() => window.scrollTo(0, Number(saved)));
    }
  } catch {}

  const viewers = () => $$(".library-frame .aln-wrap");
  const updateViewers = (patch) => viewers().forEach((wrap) => {
    wrap.dispatchEvent(new CustomEvent("aln:update", { detail: patch }));
  });

  // --- Colours. The same interpolation as the server (ancestors.rgb / colour), so the
  // tiles the page draws match the ones it was sent.
  const stops = JSON.parse($("#palette-stops").textContent);
  const tint = (pct, palette) => {
    const scale = stops[palette] || Object.values(stops)[0];
    const p = Math.min(Math.max(pct / 100, 0), 1);
    let i = 1;
    while (i < scale.length - 1 && p > scale[i][0]) i++;
    const [low, a] = scale[i - 1], [high, b] = scale[i];
    const t = high > low ? (p - low) / (high - low) : 0;
    const rgb = a.map((x, k) => Math.round(x + (b[k] - x) * t));
    const luminance = (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]) / 255;
    return { fill: `rgb(${rgb.join(", ")})`, text: luminance > 0.55 ? "#111111" : "#ffffff" };
  };
  // The site maps' tiles name a class per whole percentage; a new palette rewrites them,
  // which recolours every tile on the page at once without touching the markup.
  const paintTiles = (palette) => {
    const rules = [];
    for (let pct = 0; pct <= 100; pct++) {
      const { fill, text } = tint(pct, palette);
      rules.push(`.site-map .q${pct}{background-color: ${fill}; color: ${text}}`);
    }
    $("#tile-colours").textContent = rules.join("");
    $$(".legend-tile").forEach((tile) => {
      const { fill, text } = tint(Number(tile.dataset.pp) * 100, palette);
      tile.style.backgroundColor = fill;
      tile.style.color = text;
    });
  };

  // --- Threshold. Compared as a whole percentage, as the server does
  // (library_glyphs.percent), so the counts here are the server's. A site's alternatives
  // travel in data-alts ("R35 Q25"); each sparkline's ticks in data-ticks ("x:35 ...").
  let current = Math.max(Math.round(Number($('input[name="altall"]').value) * 100), 1);
  const drawMap = (details) => {
    $$(".sm-site[data-top]", details).forEach((site) => {
      site.classList.toggle("uncertain", Number(site.dataset.top) >= current);
      if (!site.dataset.alts) return;
      $$(".sm-alt", site).forEach((tile) => tile.remove());
      for (const code of site.dataset.alts.split(" ")) {
        const pct = Number(code.slice(1));
        if (pct < current) continue;
        const tile = document.createElement("span");
        tile.className = `sm-tile sm-alt q${pct}`;
        tile.textContent = code[0];
        site.appendChild(tile);
      }
    });
    details.dataset.drawn = String(current);
  };
  const applyThreshold = (t) => {
    current = Math.max(Math.round(t * 100), 1);
    $$('tbody[id^="target-"]').forEach((body) => {
      const svg = $("svg[data-ticks]", body);
      if (!svg) return;
      const ticks = svg.dataset.ticks
        ? svg.dataset.ticks.split(" ").map((c) => c.split(":").map(Number)) : [];
      const shown = ticks.filter(([, top]) => top >= current);
      const base = Number(svg.dataset.baseline);
      $("path.ticks", svg).setAttribute("d", shown.map(([x]) => `M${x} ${base}V${base - 6}`).join(""));
      $$('[data-count="uncertain"]', body).forEach((cell) => { cell.textContent = shown.length; });
    });
    $$("details.site-map-details[open]").forEach((details) => {
      if ($(".site-map", details).dataset.loaded) drawMap(details);
    });
    $$("[data-threshold-text]").forEach((e) => { e.textContent = (current / 100).toFixed(2); });
    updateViewers({ threshold: current });
    setParam("altall", (current / 100).toFixed(2));
  };

  // --- "Updating library…". Shown while the server is working out a design, whether the
  // page is refreshing its Design section through HTMX or submitting and reloading.
  //
  // A design can come back in a blink or take seconds, so the overlay waits before
  // appearing and then stays long enough to be read: flashing it for 150 ms would be
  // worse than not showing it at all.
  const BUSY_WAIT = 250;    // ms before an update is worth announcing
  const BUSY_LEAST = 500;   // ms it stays once it has appeared
  const BUSY_LIMIT = 60000; // a missed event must not strand it over the page
  const busyBox = $("#library-busy");
  let busyTimer = null, busyGuard = null, busyShown = 0;
  const hideBusy = () => {
    clearTimeout(busyTimer);
    clearTimeout(busyGuard);
    busyTimer = busyGuard = null;
    if (!busyBox || busyBox.hidden) return;
    const left = BUSY_LEAST - (Date.now() - busyShown);
    setTimeout(() => { if (!busyTimer) busyBox.hidden = true; }, Math.max(0, left));
  };
  const showBusy = (label) => {
    if (!busyBox) return;
    if (label) $(".busy-label", busyBox).textContent = label;
    if (busyTimer || !busyBox.hidden) return;
    busyTimer = setTimeout(() => {
      busyTimer = null;
      busyShown = Date.now();
      busyBox.hidden = false;
    }, BUSY_WAIT);
    clearTimeout(busyGuard);
    busyGuard = setTimeout(hideBusy, BUSY_LIMIT);
  };
  // An `outerHTML` swap detaches the form that asked for it, so htmx retargets the end of
  // the request: `htmx:afterRequest` arrives at `body`, and the swap events at the new
  // `#design`. Matching the *form* there left the overlay up for ever. The design refresh
  // is the only htmx on this page — the site maps use `fetch` — so any of these means it
  // is done, and each hides it whichever arrives first.
  document.body.addEventListener("htmx:beforeRequest", (event) => {
    if (event.target && event.target.id === "design-form") showBusy("Updating library…");
  });
  [
    "htmx:afterRequest", "htmx:afterSettle", "htmx:afterSwap",
    "htmx:responseError", "htmx:sendError", "htmx:timeout", "htmx:swapError",
  ].forEach((name) => document.body.addEventListener(name, hideBusy));
  // Submissions that leave the page: the overlay stands until the new one arrives, and a
  // reader who comes back with the Back button finds it gone (bfcache included).
  document.addEventListener("submit", (event) => {
    const form = event.target;
    if (!(form instanceof HTMLElement) || form.hasAttribute("hx-get")) return;
    if (!form.closest(".library-page, #library-panel, #design, .library-settings, main")) return;
    const from = event.submitter;
    showBusy((from && from.dataset.busy) || form.dataset.busy || "Updating library…");
  });
  window.addEventListener("pageshow", () => {
    clearTimeout(busyTimer);
    clearTimeout(busyGuard);
    busyTimer = busyGuard = null;
    if (busyBox) busyBox.hidden = true;
  });

  // --- The Controls panel.
  const controls = $("#library-controls");
  const frame = (fn) => { // at most one update per animation frame while a slider moves
    let pending = false;
    return (...args) => {
      if (pending) return;
      pending = true;
      requestAnimationFrame(() => { pending = false; fn(...args); });
    };
  };
  const thresholdSlider = $(".threshold-slider", controls);
  const thresholdField = $('input[name="altall"]', controls);
  const applySoon = frame((value) => applyThreshold(value));
  const clamp = (value) => Math.min(
    Math.max(Number(value) || 0, Number(thresholdSlider.min)), Number(thresholdSlider.max),
  );
  thresholdSlider.addEventListener("input", () => {
    thresholdField.value = Number(thresholdSlider.value).toFixed(2);
    applySoon(clamp(thresholdSlider.value));
    syncDesign(clamp(thresholdSlider.value));
  });
  thresholdField.addEventListener("input", () => {
    thresholdSlider.value = thresholdField.value;
    applySoon(clamp(thresholdField.value));
    syncDesign(clamp(thresholdField.value));
  });

  // The Design section is computed by the server from the same threshold, so when the
  // slider settles the section is refreshed in place (HTMX swaps it; the page keeps its
  // scroll). Its forms carry a copy of the threshold, which must be kept in step or an
  // Update would send the value the page was loaded with.
  const designForms = () => $$('#design-form, #design-save');
  // Palette, scheme and zoom act in place, so the Save form's copy of them — and the
  // design form's, which redraws the members' viewer — would otherwise still hold what
  // the page was loaded with.
  const rememberView = (what, value) => {
    $$(`[data-save-${what}]`).forEach((field) => { field.value = value; });
  };
  let designTimer = null;
  const syncDesign = (value) => {
    designForms().forEach((form) => {
      const field = $('input[name="altall"]', form);
      if (field) field.value = Number(value).toFixed(2);
    });
    const form = $("#design-form");
    if (!form || !window.htmx) return;
    clearTimeout(designTimer);
    designTimer = setTimeout(() => {
      const shown = $("[data-design-threshold]");
      if (shown && Number(shown.value) === Number(thresholdField.value)) {
        window.htmx.trigger(form, "submit");
      }
    }, 400);  // after the slider stops, not on every pixel of its travel
  };

  const BUDGET_SLIDER_BASE = 192;  // the slider's default reach (library_views.BUDGET_SLIDER_BASE)
  const budgetFor = (root) => {
    const slider = $(".budget-slider", root);
    const field = $("#budget-number", root);
    if (!slider || !field) return;
    const maxLabel = $(".budget-slider-max", root);
    // Default max 192; once the entered value passes it, the slider spans twice that value
    // so the thumb never pins at the end. The end label follows.
    const fitMax = (value) => {
      const max = value > BUDGET_SLIDER_BASE ? value * 2 : BUDGET_SLIDER_BASE;
      if (Number(slider.max) !== max) {
        slider.max = String(max);
        if (maxLabel) maxLabel.textContent = max;
      }
    };
    fitMax(Number(field.value) || 0);
    slider.addEventListener("input", () => { field.value = slider.value; });
    slider.addEventListener("change", () => {
      field.value = slider.value;
      field.dispatchEvent(new Event("change", { bubbles: true }));
    });
    field.addEventListener("input", () => {
      fitMax(Number(field.value) || 0);
      slider.value = field.value;
    });
  };
  budgetFor(document);
  document.body.addEventListener("htmx:afterSwap", (event) => {
    if (event.target && event.target.id === "design") budgetFor(event.target);
  });

  // --- Colour scheme. One choice paints every sequence on the page: the viewer frames
  // (alignment.js draws them) and the site map under each ancestor (q-classes for the
  // posterior scheme, a colour per residue for the others).
  const schemeField = $("[data-scheme]", controls);
  // Each scheme's key, from the same table the server draws its legend from, so the panel
  // can swap it without a round trip (phylobot.alignment_views.LEGENDS).
  const SCHEME_KEYS = {
    chemistry: [["Acidic D E", "#eaa57b"], ["Basic K R H", "#7fb1d4"],
                ["Aliphatic A V L I M", "#f5ee9a"], ["Aromatic F W Y", "#e3b3cc"],
                ["Polar S T N Q", "#7fcbb4"], ["Special G P C", "#b0d9f2"]],
    hydropathy: [["Hydrophilic (\u22124.5)", "#6f9fcf"], ["0", "#ffffff"],
                 ["Hydrophobic (+4.5)", "#ec9e4a"]],
    structure: [["Helix former", "#eaa57b"], ["Strand former", "#7fb1d4"],
                ["Turn former", "#7fcbb4"]],
  };
  const showSchemeLegend = () => {
    const keys = SCHEME_KEYS[scheme] || [];
    let list = $("[data-scheme-legend]", controls);
    if (!keys.length) { if (list) list.remove(); return; }
    if (!list) {
      list = document.createElement("ul");
      list.className = "panel-field scheme-legend";
      list.setAttribute("data-scheme-legend", "");
      $("[data-posterior-legend]", controls).after(list);
    }
    list.innerHTML = keys.map(([label, colour]) =>
      `<li><span class="swatch" style="background: ${colour}"></span>${label}</li>`).join("");
  };
  let scheme = schemeField ? schemeField.value : "posterior";
  const residueColour = () => (window.alignmentResidueColours || {})[scheme];
  const paintSiteMaps = () => {
    const colour = residueColour();
    $$(".site-map .sm-tile").forEach((tile) => {
      if (!colour) {  // the posterior scheme: the q-classes carry the colour
        tile.style.backgroundColor = "";
        tile.style.color = "";
        return;
      }
      const fill = colour(tile.textContent.trim());
      tile.style.backgroundColor = fill || "transparent";
      tile.style.color = fill ? "#111" : "";
    });
  };
  if (schemeField) {
    schemeField.addEventListener("change", () => {
      scheme = schemeField.value;
      const posterior = scheme === "posterior";
      updateViewers({ scheme });
      paintSiteMaps();
      setParam("scheme", scheme);
      rememberView("scheme", scheme);
      const paletteField = $("[data-palette]", controls);
      paletteField.disabled = !posterior;
      paletteField.closest(".panel-field").classList.toggle("inactive", !posterior);
      const legend = $("[data-posterior-legend]", controls);
      if (legend) legend.classList.toggle("inactive", !posterior);
      showSchemeLegend();
    });
  }
  paintSiteMaps();  // sequences already open when the page loads under a residue scheme

  const palette = $("[data-palette]", controls);
  palette.addEventListener("change", () => {
    paintTiles(palette.value);
    updateViewers({ palette: stops[palette.value] });
    setParam("palette", palette.value);
    rememberView("palette", palette.value);
  });

  // Zoom: from the whole alignment in view ("Fit") to one letter per site. The scale is
  // exponential, so each step of the slider is the same proportional change.
  // A slider sits on the ancestors' viewer and another on the members' (the Design
  // section, which HTMX replaces, so they are found when used). They set one zoom for
  // every viewer and move together; its value rides in the forms, so a reload that
  // needs the server, or a refreshed Design section, keeps it.
  const zoomField = $('input[name="zoom"]', controls);
  const MIN_CELL = 0.05, MAX_CELL = 20;
  const cellAt = (v) => MIN_CELL * (MAX_CELL / MIN_CELL) ** (Number(v) / 100);
  const sliderAt = (cell) => 100 * Math.log(cell / MIN_CELL) / Math.log(MAX_CELL / MIN_CELL);
  // Where each slider sits: at the zoom, or under "fit" at the width its own frame fitted.
  const placeSliders = () => $$(".zoom-slider").forEach((slider) => {
    const wrap = $(".aln-wrap", slider.closest(".library-frame"));
    const cell = zoomField.value && zoomField.value !== "fit"
      ? Number(zoomField.value) : Number(wrap && wrap.dataset.cell) || 1;
    slider.value = String(Math.round(sliderAt(cell)));
  });
  const applyZoom = frame((cell) => {
    updateViewers({ cell });
    const value = cell === "fit" ? "fit" : cell.toFixed(2);
    zoomField.value = value;
    setParam("zoom", value);
    rememberView("zoom", value);
    requestAnimationFrame(placeSliders);  // the other slider follows; "fit" is now drawn
  });
  document.addEventListener("input", (event) => {
    if (event.target.matches(".zoom-slider")) applyZoom(cellAt(event.target.value));
  });
  document.addEventListener("click", (event) => {
    if (event.target.closest(".zoom-fit")) applyZoom("fit");
  });
  placeSliders();
  document.body.addEventListener("htmx:afterSettle", (event) => {
    if (event.target && event.target.id === "design") requestAnimationFrame(placeSliders);
  });

  // "Download Excel" gives what the page shows. Palette, scheme, zoom and the design
  // change the address in place, so its link is built from the address when clicked.
  document.addEventListener("click", (event) => {
    const link = event.target.closest("[data-excel]");
    if (!link) return;
    const query = new URLSearchParams(window.location.search);
    query.set("export", "xlsx");
    link.href = `${window.location.pathname}?${query}`;
  });

  // The alignment and the extant sequences need the server: submit, keeping the place.
  // The design's own fields (wells, variants, strategy) ride along, or the reload would
  // show the saved design again in place of the one on the page.
  const carryDesign = () => {
    const design = $("#design-form");
    if (!design) return;
    const own = new Set(new FormData(controls).keys());
    for (const [name, value] of new FormData(design)) {
      if (own.has(name)) continue;
      const field = document.createElement("input");
      Object.assign(field, { type: "hidden", name, value });
      controls.append(field);
    }
  };
  $$("[data-reload]", controls).forEach((control) => control.addEventListener("change", () => {
    keepPlace(controls);
    carryDesign();
    controls.requestSubmit();
  }));

  const search = $(".extant-search", controls);
  if (search) {
    search.addEventListener("input", () => {
      const wanted = search.value.trim().toLowerCase();
      $$(".extant-list label", controls).forEach((label) => {
        label.hidden = Boolean(wanted) && !label.textContent.toLowerCase().includes(wanted);
      });
    });
  }

  // --- Site maps: fetched when a row is first opened (the page would otherwise carry
  // every residue of every ancestor), and redrawn in place as the threshold moves.
  const loadMap = async (details) => {
    const map = $(".site-map", details);
    if (map.dataset.loaded) {
      if (Number(details.dataset.drawn) !== current) drawMap(details);
      return;
    }
    map.dataset.loaded = "1";
    const asked = current; // the slider may move while the map is on its way
    try {
      const response = await fetch(`${map.dataset.src}&altall=${(asked / 100).toFixed(2)}`);
      if (!response.ok) throw new Error(String(response.status));
      map.innerHTML = await response.text();
      paintSiteMaps();  // a map fetched under a residue scheme arrives in q-classes
      details.dataset.drawn = String(asked);
      if (asked !== current) drawMap(details);
    } catch {
      delete map.dataset.loaded;
      map.innerHTML = '<p class="error">Could not load this sequence; close and open it to try again.</p>';
    }
  };
  $$("details.site-map-details").forEach((details) => {
    details.addEventListener("toggle", () => { if (details.open) loadMap(details); });
    if (details.open) loadMap(details);
  });
  // Hover read-outs: the fixed part comes from the server, the alternatives from data-alts.
  $$(".site-map").forEach((map) => map.addEventListener("mouseover", (event) => {
    const site = event.target.closest(".sm-site[data-alts]");
    if (!site || site.dataset.titled) return;
    const alts = site.dataset.alts.split(" ")
      .map((c) => `${c[0]} ${(Number(c.slice(1)) / 100).toFixed(2)}`);
    site.title = `${site.title}, ${alts.join(", ")}`;
    site.dataset.titled = "1";
  }));
  const openTarget = (anchor) => {
    const details = $(`details[data-target="${anchor}"]`);
    if (!details) return;
    details.open = true;
    details.scrollIntoView({ block: "nearest" });
  };
  $$("[data-expand]").forEach((button) => button.addEventListener("click", () => {
    $$("details.site-map-details").forEach((d) => { d.open = button.dataset.expand === "open"; });
  }));
  // Clicking an ancestor's row in a viewer opens its sequence below; clicking a member's
  // finds its row in the members table. Delegated, since the members' viewer arrives
  // with every refresh of the Design section.
  const findMember = (anchor) => {
    const row = document.getElementById(anchor);
    if (!row) return;
    $$(".members-table tr.found").forEach((r) => r.classList.remove("found"));
    row.classList.add("found");
    row.scrollIntoView({ block: "center", behavior: "smooth" });
  };
  document.addEventListener("aln:click", (event) => {
    const frameEl = event.target.closest(".library-frame");
    if (!frameEl) return;
    const data = JSON.parse($('script[id^="aln-data"]', frameEl).textContent);
    const anchor = (data.anchors || [])[event.detail.row];
    if (!anchor) return;
    if (anchor.startsWith("member-")) findMember(anchor);
    else openTarget(anchor);
  });
  // And the other way: a member's name in the table shows its row in the viewer.
  document.addEventListener("click", (event) => {
    const button = event.target.closest("[data-member]");
    if (!button) return;
    const anchor = button.dataset.member;
    $$(".library-frame-m").forEach((frameEl) => {
      const data = JSON.parse($('script[id^="aln-data"]', frameEl).textContent);
      const wrap = $(".aln-wrap", frameEl);
      const found = (data.anchors || []).includes(anchor);
      wrap.dispatchEvent(new CustomEvent("aln:update", { detail: { highlight: found ? anchor : "" } }));
      if (found) frameEl.scrollIntoView({ block: "nearest", behavior: "smooth" });
    });
  });
  // Which members the viewer draws: part of the design form (form=), refreshed in place.
  document.addEventListener("change", (event) => {
    const choice = event.target.closest("[data-members-shown]");
    if (!choice) return;
    rememberView("members", choice.value);
    setParam("members", choice.value);
    const form = $("#design-form");
    if (form && window.htmx) window.htmx.trigger(form, "submit");
    else if (form) form.requestSubmit();
  });
  $$("[data-copy]").forEach((button) => button.addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText(button.dataset.copy);
      const was = button.textContent;
      button.textContent = "Copied";
      setTimeout(() => { button.textContent = was; }, 1500);
    } catch { toast("Could not copy; select the sequence instead."); }
  }));

  // --- The Controls panel slides in from the right, as the New Job page's does.
  const panel = $("#library-panel"), tab = $("#library-panel-tab");
  const showPanel = (open) => {
    panel.classList.toggle("open", open);
    panel.inert = !open;
    tab.hidden = open;
    tab.setAttribute("aria-expanded", String(open));
    document.body.classList.toggle("library-panel-open", open);
    try { localStorage.setItem("phylobot_library_panel2", open ? "open" : "closed"); } catch {}
  };
  $(".job-panel-close", panel).addEventListener("click", () => showPanel(false));
  tab.addEventListener("click", () => showPanel(true));
  // Key bumped to _panel2 so the controls open by default: the panel holds the extant-sequence
  // chooser and the layout menu, and a stale "closed" from before was hiding them.
  let wanted = "open";
  try { wanted = localStorage.getItem("phylobot_library_panel2") || "open"; } catch {}
  showPanel(wanted === "open" && window.innerWidth > 900);

  // --- Bulk actions on the ticked ancestors.
  const targetsForm = $("#targets-form");
  const ticked = () => $$('input[name="target"]:checked');
  $$("[data-select-all]").forEach((all) => all.addEventListener("change", () => {
    $$('input[name="target"]').forEach((box) => { box.checked = all.checked; });
  }));

  const dialog = $("#library-move");
  const dialogList = $(".library-choices", dialog);
  const dialogSearch = $(".library-choice-search", dialog);
  let action = "copy";
  const openDialog = (which) => {
    action = which;
    const copying = which === "copy";
    $(".library-move-title", dialog).textContent =
      `${copying ? "Copy" : "Move"} ${ticked().length} ancestor${ticked().length === 1 ? "" : "s"}`;
    $(".library-move-hint", dialog).textContent = copying
      ? "Choose one or more libraries to copy them into."
      : "Choose the library to move them to; they leave this one.";
    $(".library-move-go", dialog).textContent = copying ? "Copy" : "Move";
    $$("input", dialogList).forEach((input) => {
      input.type = copying ? "checkbox" : "radio";
      input.checked = false;
    });
    dialogSearch.value = "";
    $$("label", dialogList).forEach((label) => { label.hidden = false; });
    dialog.showModal();
    dialogSearch.focus();
  };
  dialogSearch.addEventListener("input", () => {
    const wantedText = dialogSearch.value.trim().toLowerCase();
    $$("label", dialogList).forEach((label) => {
      label.hidden = Boolean(wantedText) && !label.textContent.toLowerCase().includes(wantedText);
    });
  });
  $(".library-move-go", dialog).addEventListener("click", (event) => {
    event.preventDefault();
    const chosen = $$("input:checked", dialogList).map((input) => input.value);
    if (!chosen.length) {
      $(".library-move-error", dialog).textContent = "Choose a library first.";
      return;
    }
    $$('input[name="to"]', targetsForm).forEach((input) => input.remove());
    for (const value of chosen) {
      const field = document.createElement("input");
      Object.assign(field, { type: "hidden", name: "to", value });
      targetsForm.appendChild(field);
    }
    $("[data-action-field]", targetsForm).value = action;
    dialog.close();
    keepPlace(null);
    targetsForm.requestSubmit();
  });

  $$("[data-bulk]").forEach((button) => button.addEventListener("click", (event) => {
    event.preventDefault();
    const count = ticked().length;
    if (!count) return toast("Tick the ancestors to act on first.");
    $(".library-move-error", dialog).textContent = "";
    if (button.dataset.bulk === "remove") {
      confirmDialog({
        title: `Remove ${count} ancestor${count === 1 ? "" : "s"}?`,
        body: "They leave this library. The job and its results are not touched.",
        action: "Remove",
        danger: true,
      }, () => {
        $("[data-action-field]", targetsForm).value = "remove";
        keepPlace(null);
        targetsForm.requestSubmit();
      });
      return;
    }
    openDialog(button.dataset.bulk);
  }));

  $("#delete-library").addEventListener("click", (event) => {
    confirmDialog({
      title: `Delete “${event.currentTarget.dataset.name}”?`,
      body: "It leaves this job's libraries. The job and its results are not touched.",
      action: "Delete",
      danger: true,
    }, () => $("#delete-library-form").submit());
  });
})();
