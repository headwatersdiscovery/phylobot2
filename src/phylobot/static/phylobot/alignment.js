// Alignment viewer (phylobot.alignment_views): the grid and tracks drawn in one canvas
// from the page's JSON. Only the visible window is drawn (pilot jobs reach 500 x 5000
// cells), so the canvas is the size of the view and a transparent scroller on top of it
// provides native scrolling; names (left) and tracks (top) stay in place. Re-initialised
// after each HTMX swap.
//
// Other pages embed it too (the library page, phylobot.library_views): each viewer is an
// `.aln-root` holding its own `script#aln-data-<n>`, and an "aln:update" event on its
// `.aln-wrap` (detail: palette, threshold, cell) redraws it in place. The "posterior"
// scheme colours ancestor rows by `data.pp` (whole percentages) through `data.palette`,
// outlines sites whose strongest alternative (`data.top`) reaches `data.threshold`, and
// fades the cells in `data.faded`. `data.cell` may be "fit" (the whole width in view).
// Library members (`data.member`) are rows of their own: `data.marked` outlines their
// changes from the ML sequence, `data.dim` dims a whole row, and an "aln:update" with
// `highlight` (a row's anchor) scrolls to that row and rings it.
//
// The mutational path page (phylobot.path_views) draws a path's ancestors, oldest first:
// `data.marked` is then each row's sites that differ from the row above. Its "changes"
// scheme leaves the path's rows (`data.onPath`) grey and colours only those sites: by
// posterior probability, or through `data.levels` (per row, column -> index into
// `data.levelColours`). Under "posterior" the marks are rings in `data.markColour`,
// `data.markWidth` wide. `data.blank` rows are kept empty (an ancestor another
// combination has and this one lacks), `data.notes` adds a line to a cell's tooltip, and
// `data.namesWidth` fixes the names' gutter so stacked viewers line up. With
// `data.hoverEvents` the viewer reports the column under the pointer ("aln:hover", detail:
// column, -1 on leaving), and an "aln:update" with `columns` boxes those columns from top
// to bottom: the page uses the two to show the same site in every viewer.
//
// The alignment editor (phylobot.alignment_editor_views) sets `data.editable`: dragging
// across the ruler or the grid selects columns (Shift extends, Ctrl/Cmd adds), clicking or
// dragging over names selects sequences, and the names' header ("Name", with the order in
// `data.sort`) asks for a sort. The viewer reports its selection in an "aln:select" event
// (detail: columns, rows: positions shown, 0-based) and an "aln:sort" event; an
// "aln:update" with `clear` empties the selection. The editor applies edits itself.
// `data.origin` gives each column's number before the editor's cuts (for the tooltip).
(() => {
  const GAP = "-";
  const AMINO = "ARNDCQEGHILKMFPSTWYV";
  const RULER = 18, TRACK = 44, TRACK_GAP = 6;
  const MIN_CELL = 0.05, MAX_CELL = 20; // the zoom slider's ends (library_glyphs too)
  const MONO = "ui-monospace, SFMono-Regular, Menlo, monospace";
  const SERIES = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"];

  // Colours shared with the legends (alignment_views.LEGENDS).
  const CLUSTAL = { blue: "#80a0f0", red: "#f01505", green: "#15c015", pink: "#f08080",
    magenta: "#c048c0", orange: "#f09048", cyan: "#15a4a4", yellow: "#c0c000" };
  const CHEMISTRY = {};
  [["DE", "#eaa57b"], ["KRH", "#7fb1d4"], ["AVLIM", "#f5ee9a"], ["FWY", "#e3b3cc"],
   ["STNQ", "#7fcbb4"], ["GPC", "#b0d9f2"]].forEach(([set, colour]) => {
    for (const r of set) CHEMISTRY[r] = colour;
  });
  const KD = { A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5, Q: -3.5, E: -3.5, G: -0.4, H: -3.2,
    I: 4.5, L: 3.8, K: -3.9, M: 1.9, F: 2.8, P: -1.6, S: -0.8, T: -0.7, W: -0.9, Y: -1.3, V: 4.2 };
  const CF = { A: [1.42, 0.83, 0.66], R: [0.98, 0.93, 0.95], N: [0.67, 0.89, 1.56],
    D: [1.01, 0.54, 1.46], C: [0.70, 1.19, 1.19], Q: [1.11, 1.10, 0.98], E: [1.51, 0.37, 0.74],
    G: [0.57, 0.75, 1.56], H: [1.00, 0.87, 0.95], I: [1.08, 1.60, 0.47], L: [1.21, 1.30, 0.59],
    K: [1.14, 0.74, 1.01], M: [1.45, 1.05, 0.60], F: [1.13, 1.38, 0.60], P: [0.57, 0.55, 1.52],
    S: [0.77, 0.75, 1.43], T: [0.83, 1.19, 0.96], W: [1.08, 1.37, 0.96], Y: [0.69, 1.47, 1.14],
    V: [1.06, 1.70, 0.50] };
  const STRUCTURE = [[234, 165, 123], [127, 177, 212], [127, 203, 180]]; // helix, strand, turn
  const GAP_DIFF = "#c8ccd4";

  function mix(rgb, t) { // blend towards white: t = 0 white, 1 the colour
    return `rgb(${rgb.map((c) => Math.round(255 + (c - 255) * t)).join(",")})`;
  }
  function hydropathy(r) {
    const v = KD[r];
    if (v === undefined) return null;
    const t = Math.min(1, Math.abs(v) / 4.5);
    return mix(v < 0 ? [111, 159, 207] : [236, 158, 74], t);
  }
  function structure(r) {
    const p = CF[r];
    if (!p) return null;
    const k = p.indexOf(Math.max(...p));
    if (p[k] < 1.05) return null;
    return mix(STRUCTURE[k], Math.min(1, 0.35 + (p[k] - 1.05) / 0.6));
  }

  // Clustal X (Jalview's rules): a residue is coloured only if its column passes one of
  // its thresholds. "KR:60" = K and R together above 60%; "K,R,Q:80" = any one above 80%.
  const RULES = {
    A: ["blue", ["WLVIMAFCHP:60"]], I: ["blue", ["WLVIMAFCHP:60"]], L: ["blue", ["WLVIMAFCHP:60"]],
    M: ["blue", ["WLVIMAFCHP:60"]], F: ["blue", ["WLVIMAFCHP:60"]], W: ["blue", ["WLVIMAFCHP:60"]],
    V: ["blue", ["WLVIMAFCHP:60"]],
    K: ["red", ["KR:60", "K,R,Q:80"]], R: ["red", ["KR:60", "K,R,Q:80"]],
    N: ["green", ["N:50", "N,Y:85"]],
    Q: ["green", ["KR:60", "QE:50", "Q,E,K,R:85"]],
    E: ["magenta", ["KR:60", "QE:50", "E,Q,D:85"]],
    D: ["magenta", ["KR:60", "K,R,Q:85", "ED:50"]],
    G: ["orange", ["G:0"]], P: ["yellow", ["P:0"]],
    H: ["cyan", ["WLVIMAFCHP:60", "W,Y,A,C,P,Q,F,H,I,L,M,V:85"]],
    Y: ["cyan", ["WLVIMAFCHP:60", "W,Y,A,C,P,Q,F,H,I,L,M,V:85"]],
    S: ["green", ["WLVIMAFCHP:60", "TS:50", "S,T:85"]],
    T: ["green", ["WLVIMAFCHP:60", "TS:50", "S,T:85"]],
  };
  function passes(rule, share) {
    const [set, pct] = rule.split(":");
    const limit = Number(pct) / 100;
    if (set.includes(",")) return set.split(",").some((r) => (share[r] || 0) > limit);
    let sum = 0;
    for (const r of set) sum += share[r] || 0;
    return limit === 0 ? sum > 0 : sum > limit;
  }

  // Where the user had scrolled, kept across HTMX swaps (a new colour scheme, track,
  // zoom or order re-renders the viewer): in columns and the top sequence's name, not
  // pixels, so it survives zoom and row-order changes.
  let current = null, saved = null;
  const SELECTED = "rgba(47, 111, 222, 0.24)", SELECTED_RULER = "rgba(47, 111, 222, 0.6)";
  document.addEventListener("htmx:beforeRequest", (event) => {
    if (current && event.target.id === "alignment-viewer") saved = current();
  });

  function init(root) {
    const source = root.querySelector('script[id^="aln-data"]');
    const wrap = root.querySelector(".aln-wrap");
    if (!source || !wrap || wrap.dataset.ready) return;
    wrap.dataset.ready = "1";
    const data = JSON.parse(source.textContent);
    const canvas = wrap.querySelector("canvas");
    const scroller = wrap.querySelector(".aln-scroller");
    const sizer = wrap.querySelector(".aln-sizer");
    const tip = wrap.querySelector(".aln-tip");
    const ctx = canvas.getContext("2d");
    const style = getComputedStyle(wrap);
    const colours = {
      bg: style.getPropertyValue("--aln-bg").trim() || "#fff",
      fg: style.getPropertyValue("--aln-fg").trim() || "#1d2330",
      muted: style.getPropertyValue("--aln-muted").trim() || "#8a929e",
      line: style.getPropertyValue("--aln-line").trim() || "#d8dde5",
      band: style.getPropertyValue("--aln-band").trim() || "#f3f5f8",
      selected: style.getPropertyValue("--aln-selected").trim() || SELECTED,
      quiet: style.getPropertyValue("--aln-quiet").trim() || "#e4e7ec",
    };
    const editable = Boolean(data.editable);
    let selCols = new Set(), selRows = new Set();
    let anchorCol = -1, anchorRow = -1, drag = null;
    const unit = data.codon ? 3 : 1; // characters per site
    let cell = Number(data.cell) || 12; // "fit" is resolved once the width is known
    let cw = cell * (data.codon ? 2.2 : 1); // a codon tile is wider, not 3x
    // Rows are as tall as the page says, or (the alignment viewer, whose zoom is
    // continuous) follow the cell width: the old menu's steps, 3 -> 4, 6 -> 9, 11 -> 15
    // and 15 -> 19 px, joined by straight lines.
    const rowFor = (c) => Math.round(c <= 3 ? 4 : c <= 6 ? 4 + (c - 3) * 5 / 3
      : c <= 11 ? 9 + (c - 6) * 6 / 5 : Math.min(c + 4, 24));
    let rh = data.row;
    const rows = data.rows;
    const n = rows.length;
    const width = rows.reduce((m, r) => Math.max(m, r.length / unit), 0);
    const tracksH = RULER + data.tracks.length * (TRACK + TRACK_GAP);
    ctx.font = `12px ${MONO}`;
    const widestName = Math.max(80, ...data.names.map((s) => ctx.measureText(s).width)) + 18;
    // The tracks' names share the gutter with the sequence names (beside each track's
    // scale, 40px): where the names are short, make room for them rather than cut
    // "Conservation (Jensen-Shannon)" down to "Conserva…".
    ctx.font = "11px system-ui, sans-serif";
    const widestTrack = Math.min(215, Math.max(0, ...data.tracks.map((t) => ctx.measureText(t.label).width)) + 46);
    ctx.font = `12px ${MONO}`;
    const namesW = Number(data.namesWidth) || Math.min(280, Math.max(widestName, widestTrack));
    const fitCell = () => Math.max(0.05,
      (wrap.clientWidth - namesW - 2) / Math.max(width, 1) / (data.codon ? 2.2 : 1));
    let fitting = false; // "fit" follows the viewer's width as the window changes
    const setCell = (value) => {
      fitting = value === "fit";
      cell = fitting ? fitCell() : Math.max(0.05, Number(value) || 12);
      cw = cell * (data.codon ? 2.2 : 1);
      if (data.rowFromCell) rh = rowFor(cell);
      wrap.dataset.cell = String(cell); // so a zoom slider can show where "fit" landed
    };
    setCell(data.cell);

    // The posterior scheme: palette colours per whole percentage, and text that reads on them.
    const faded = (data.faded || []).map((columns) => new Set(columns));
    const marked = (data.marked || []).map((columns) => new Set(columns));
    const dim = data.dim || [];
    const member = data.member || [];
    const onPath = data.onPath || [], blank = data.blank || [];
    const levels = data.levels || [], notes = data.notes || [];
    const levelColours = data.levelColours || [];
    const markColour = data.markColour || colours.fg, markWidth = Number(data.markWidth) || 1.5;
    const CHANGED = "#6b7280"; // a changed site with no probability to colour it by
    const textOn = (hex) => { // black or white, whichever reads on a #rrggbb fill
      const v = parseInt(hex.slice(1), 16);
      return (0.299 * (v >> 16) + 0.587 * ((v >> 8) & 255) + 0.114 * (v & 255)) / 255 > 0.55 ? "#111111" : "#ffffff";
    };
    let highlighted = -1;
    let boxed = [], hovered = -1; // columns boxed by the page; the column last reported
    const tints = new Map();
    function tint(pct) {
      if (tints.has(pct)) return tints.get(pct);
      const stops = data.palette || [[0, [215, 48, 39]], [1, [69, 117, 180]]];
      const p = Math.min(Math.max(pct / 100, 0), 1);
      let i = 1;
      while (i < stops.length - 1 && p > stops[i][0]) i++;
      const [lo, a] = stops[i - 1], [hi, b] = stops[i];
      const t = hi > lo ? (p - lo) / (hi - lo) : 0;
      const rgb = a.map((x, k) => Math.round(x + (b[k] - x) * t));
      const light = (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]) / 255 > 0.55;
      const out = { fill: `rgb(${rgb.join(",")})`, text: light ? "#111111" : "#ffffff" };
      tints.set(pct, out);
      return out;
    }
    const ppAt = (r, c) => {
      const row = data.pp && data.pp[r];
      const v = row ? row[c] : null;
      return v === null || v === undefined || v < 0 ? null : v;
    };

    const site = (r, c) => rows[r].substr(c * unit, unit);
    const aminoAt = (r, c) => (data.codon ? (data.amino[r][c] || GAP) : rows[r][c] || GAP);
    // Per-column consensus (extant rows only), of the drawn states and of amino acids.
    const consensus = new Map();
    const extantRows = data.ancestor.flatMap((ancestor, r) => (ancestor ? [] : [r]));
    function columnInfo(c) {
      if (consensus.has(c)) return consensus.get(c);
      // Residues only (Clustal X, identity), and with gaps as a state (the consensus that
      // differences compare with: a column of mostly gaps has a gap as its consensus).
      const counts = {}, amino = {};
      let gaps = 0;
      for (const r of extantRows) {
        const s = site(r, c), a = aminoAt(r, c);
        if (a === GAP) { gaps++; continue; }
        counts[s] = (counts[s] || 0) + 1;
        amino[a] = (amino[a] || 0) + 1;
      }
      let best = null, top = 0;
      for (const [s, k] of Object.entries(counts)) if (k > top) { best = s; top = k; }
      let bestAmino = null, topAmino = 0;
      for (const [a, k] of Object.entries(amino)) if (k > topAmino) { bestAmino = a; topAmino = k; }
      const share = {};
      for (const [a, k] of Object.entries(amino)) share[a] = k / data.extant;
      const gapped = gaps > top;
      const info = { state: gapped ? GAP.repeat(unit) : best, amino: gapped ? GAP : bestAmino,
        residueAmino: bestAmino, aminoShare: topAmino / data.extant, share };
      consensus.set(c, info);
      return info;
    }
    const referenceState = (c) => (data.reference >= 0 ? site(data.reference, c) : columnInfo(c).state);
    const referenceAmino = (c) => (data.reference >= 0 ? aminoAt(data.reference, c) : columnInfo(c).amino);
    const isGap = (s) => s.includes(GAP);

    // A cell's fill (null: none) and whether it is drawn as a dot (same as reference).
    function paint(r, c) {
      const s = site(r, c), a = aminoAt(r, c);
      const changed = marked[r] ? marked[r].has(c) : false;
      if (data.scheme === "changes") { // a path's rows: grey, but for what changed
        if (!onPath[r]) return { fill: null };
        if (!changed) return { fill: isGap(s) ? null : colours.quiet, quiet: true };
        const level = levelColours.length && levels[r] ? levels[r][c] : undefined;
        if (level !== undefined) return { fill: levelColours[level], text: textOn(levelColours[level]) };
        const pct = levelColours.length ? null : ppAt(r, c);
        return pct === null ? { fill: CHANGED, text: "#ffffff" } : tint(pct);
      }
      if (isGap(s)) {
        if (changed && data.markColour) return { fill: GAP_DIFF }; // a loss, so its ring shows
        if (data.scheme === "differences" || data.scheme === "synonymous") {
          const ref = referenceState(c);
          return { fill: ref && !isGap(ref) && r !== data.reference ? GAP_DIFF : null };
        }
        return { fill: null };
      }
      switch (data.scheme) {
        case "clustalx": {
          const rule = RULES[a];
          if (a === "C") {
            const share = columnInfo(c).share;
            if ((share.C || 0) >= 1) return { fill: CLUSTAL.pink };
            return { fill: passes("WLVIMAFCHP:60", share) ? CLUSTAL.blue : null };
          }
          if (!rule) return { fill: null };
          const share = columnInfo(c).share;
          return { fill: rule[1].some((x) => passes(x, share)) ? CLUSTAL[rule[0]] : null };
        }
        case "posterior": {
          const pct = ppAt(r, c);
          if (pct === null) return { fill: null };
          const { fill, text } = tint(pct);
          return { fill, text, faded: faded[r] ? faded[r].has(c) : false };
        }
        case "chemistry": return { fill: CHEMISTRY[a] || null };
        case "hydropathy": return { fill: hydropathy(a) };
        case "structure": return { fill: structure(a) };
        case "identity": {
          const info = columnInfo(c);
          if (a !== info.residueAmino) return { fill: null };
          const p = info.aminoShare;
          return { fill: p > 0.8 ? "#6464ff" : p > 0.6 ? "#9999ff" : p > 0.4 ? "#ccccff" : null };
        }
        case "differences": {
          if (r === data.reference) return { fill: CHEMISTRY[a] || null };
          const ref = data.codon ? referenceAmino(c) : referenceState(c);
          if (a === ref) return { fill: null, dot: true };
          return { fill: ref === GAP || ref === null ? GAP_DIFF : CHEMISTRY[a] || GAP_DIFF };
        }
        case "synonymous": {
          if (r === data.reference) return { fill: null };
          const ref = referenceState(c);
          if (s === ref) return { fill: null, dot: true };
          if (!ref || isGap(ref)) return { fill: GAP_DIFF };
          return { fill: a === referenceAmino(c) ? "#a8d8f2" : "#f0c56b" };
        }
        default: return { fill: null };
      }
    }

    // The overview strip (.aln-overview, where the page has one): the whole alignment
    // drawn into the strip's width, one sample per pixel, with a box over the part in view.
    // Residues take the scheme's colour (a neutral grey where the scheme leaves them
    // plain) and gaps the background, so the gapped regions show at any size. It appears
    // only when the view does not already hold every column and row.
    const overview = root.querySelector(".aln-overview");
    const OVERVIEW_H = 44;
    let overviewStale = true;
    const overviewSpan = () => Math.max(1, overview.clientWidth - namesW);
    function drawOverview() {
      overviewStale = false;
      const oc = overview.querySelector("canvas");
      const octx = oc.getContext("2d");
      const dpr = window.devicePixelRatio || 1;
      const w = overview.clientWidth;
      oc.width = Math.round(w * dpr);
      oc.height = Math.round(OVERVIEW_H * dpr);
      oc.style.width = `${w}px`;
      oc.style.height = `${OVERVIEW_H}px`;
      octx.fillStyle = colours.bg;
      octx.fillRect(0, 0, oc.width, oc.height);
      const x0 = Math.round(namesW * dpr), pw = oc.width - x0, ph = oc.height;
      const across = Math.max(1, Math.min(pw, width)), down = Math.max(1, Math.min(ph, n));
      for (let j = 0; j < down; j++) {
        const r = Math.min(n - 1, Math.floor((j + 0.5) * n / down));
        const top = Math.floor(j * ph / down), bottom = Math.floor((j + 1) * ph / down);
        for (let i = 0; i < across; i++) {
          const c = Math.min(width - 1, Math.floor((i + 0.5) * width / across));
          if (isGap(site(r, c))) continue;
          octx.fillStyle = paint(r, c).fill || colours.line;
          const left = Math.floor(i * pw / across);
          octx.fillRect(x0 + left, top, Math.floor((i + 1) * pw / across) - left, bottom - top);
        }
      }
      octx.setTransform(dpr, 0, 0, dpr, 0, 0);
      octx.fillStyle = colours.muted;
      octx.font = "11px system-ui, sans-serif";
      octx.textBaseline = "middle";
      octx.fillText("Overview", 6, OVERVIEW_H / 2);
      octx.fillStyle = colours.line;
      octx.fillRect(namesW - 1, 0, 1, OVERVIEW_H);
      octx.setTransform(1, 0, 0, 1, 0, 0);
    }
    function placeOverview() {
      if (!overview) return;
      const totalW = width * cw, totalH = n * rh;
      const viewW = scroller.clientWidth - namesW, viewH = scroller.clientHeight - tracksH;
      const needed = totalW > viewW + 1 || totalH > viewH + 1;
      if (overview.hidden === needed) {
        overview.hidden = !needed;
        overviewStale = true;
      }
      if (!needed) return;
      if (overviewStale) drawOverview();
      const span = overviewSpan(), box = overview.querySelector(".aln-overview-window").style;
      box.left = `${namesW + scroller.scrollLeft / totalW * span}px`;
      box.width = `${Math.min(1, viewW / totalW) * span}px`;
      box.top = `${scroller.scrollTop / totalH * OVERVIEW_H}px`;
      box.height = `${Math.min(1, viewH / totalH) * OVERVIEW_H}px`;
    }
    if (overview) {
      // Dragging the box keeps hold of it where it was grabbed; a click elsewhere centres
      // the view there.
      let grab = null;
      const at = (event) => {
        const box = overview.getBoundingClientRect();
        return [(event.clientX - box.left - namesW) / overviewSpan(), (event.clientY - box.top) / OVERVIEW_H];
      };
      const centre = () => [
        (scroller.scrollLeft + (scroller.clientWidth - namesW) / 2) / (width * cw),
        (scroller.scrollTop + (scroller.clientHeight - tracksH) / 2) / (n * rh)];
      const go = (event) => {
        const [fx, fy] = at(event);
        scroller.scrollLeft = (fx - grab[0]) * width * cw - (scroller.clientWidth - namesW) / 2;
        scroller.scrollTop = (fy - grab[1]) * n * rh - (scroller.clientHeight - tracksH) / 2;
      };
      overview.addEventListener("pointerdown", (event) => {
        if (event.button !== 0) return;
        event.preventDefault();
        const [fx, fy] = at(event), [cx, cy] = centre();
        const halfW = (scroller.clientWidth - namesW) / (width * cw) / 2;
        const halfH = (scroller.clientHeight - tracksH) / (n * rh) / 2;
        const inside = Math.abs(fx - cx) <= halfW && Math.abs(fy - cy) <= halfH;
        grab = inside ? [fx - cx, fy - cy] : [0, 0];
        overview.setPointerCapture(event.pointerId);
        go(event);
      });
      overview.addEventListener("pointermove", (event) => { if (grab) go(event); });
      const release = () => { grab = null; };
      overview.addEventListener("pointerup", release);
      overview.addEventListener("pointercancel", release);
    }

    function sizeCanvas() {
      const dpr = window.devicePixelRatio || 1;
      // Room for a horizontal scrollbar where the system draws one, and the border.
      const chrome = 18 + 2;
      let height;
      if (data.height === "fit") { // the window, below the sticky top bar
        const bar = document.querySelector("header");
        const top = bar && getComputedStyle(bar).position === "sticky" ? bar.offsetHeight : 0;
        height = Math.max(200, window.innerHeight - top - 16);
      } else {
        height = tracksH + Number(data.height) * rh + chrome;
      }
      height = Math.min(height, tracksH + n * rh + chrome);
      wrap.style.height = `${height}px`;
      sizer.style.width = `${namesW + width * cw}px`;
      sizer.style.height = `${tracksH + n * rh}px`;
      const w = scroller.clientWidth, h = scroller.clientHeight;
      canvas.width = Math.round(w * dpr);
      canvas.height = Math.round(h * dpr);
      canvas.style.width = `${w}px`;
      canvas.style.height = `${h}px`;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }

    function range(track) {
      let lo = track.low, hi = track.high;
      if (lo === null || hi === null) {
        let min = Infinity, max = -Infinity;
        for (const s of track.series) for (const v of s.values) {
          if (v !== null) { min = Math.min(min, v); max = Math.max(max, v); }
        }
        if (min === Infinity) { min = 0; max = 1; }
        if (lo === null) lo = min;
        if (hi === null) hi = max;
      }
      if (hi === lo) hi = lo + 1;
      return [lo, hi];
    }
    const ranges = data.tracks.map(range);
    const fmt = (v) => (Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2));

    function draw() {
      const w = scroller.clientWidth, h = scroller.clientHeight;
      const sx = scroller.scrollLeft, sy = scroller.scrollTop;
      ctx.fillStyle = colours.bg;
      ctx.fillRect(0, 0, w, h);
      const c0 = Math.max(0, Math.floor(sx / cw));
      const c1 = Math.min(width, Math.ceil((sx + w - namesW) / cw));
      const r0 = Math.max(0, Math.floor(sy / rh));
      const r1 = Math.min(n, Math.ceil((sy + h - tracksH) / rh));
      const x = (c) => namesW + c * cw - sx;
      const y = (r) => tracksH + r * rh - sy;
      const letters = cell >= 9;
      const outline = data.scheme === "posterior" && data.top && cw >= 4;

      // Grid.
      ctx.save();
      ctx.beginPath();
      ctx.rect(namesW, tracksH, w - namesW, h - tracksH);
      ctx.clip();
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.font = `${Math.min(data.codon ? 11 : 13, rh - 3)}px ${MONO}`;
      for (let r = r0; r < r1; r++) {
        if (data.ancestor[r]) {
          ctx.fillStyle = colours.band;
          ctx.fillRect(namesW, y(r), w - namesW, rh);
        }
        if (blank[r]) continue;
        if (dim[r]) ctx.globalAlpha = 0.35;
        for (let c = c0; c < c1; c++) {
          const p = paint(r, c);
          if (p.fill) {
            ctx.fillStyle = p.fill;
            if (p.faded) ctx.globalAlpha = dim[r] ? 0.15 : 0.35;
            ctx.fillRect(x(c), y(r), cw, rh);
            ctx.globalAlpha = dim[r] ? 0.35 : 1;
          }
          if (marked[r] && marked[r].has(c)) { // a member's or a path's change: at any zoom
            if (data.scheme === "changes") { /* the colour is the mark */ } else if (cw >= 4) {
              ctx.strokeStyle = markColour;
              ctx.lineWidth = markWidth;
              ctx.strokeRect(x(c) + markWidth / 2, y(r) + markWidth / 2, cw - markWidth, rh - markWidth);
            } else { // too narrow to ring: a band along the top, leaving the colour beneath
              ctx.fillStyle = markColour;
              ctx.fillRect(x(c), y(r) + 1, Math.max(1, cw), 4);
            }
          } else if (outline && data.top[r] && data.top[r][c] >= data.threshold) {
            ctx.strokeStyle = colours.fg;
            ctx.lineWidth = 1;
            ctx.strokeRect(x(c) + 0.5, y(r) + 0.5, cw - 1, rh - 1);
          }
          if (!letters) continue;
          const s = site(r, c);
          if (p.dot) {
            ctx.fillStyle = colours.muted;
            ctx.fillText("·", x(c) + cw / 2, y(r) + rh / 2);
          } else {
            ctx.fillStyle = p.quiet ? colours.muted : p.text && !p.faded ? p.text
              : p.fill && !p.faded ? (p.fill === "#6464ff" ? "#fff" : "#111") : isGap(s) ? colours.muted : colours.fg;
            ctx.fillText(data.codon && cw < 22 ? aminoAt(r, c) : s, x(c) + cw / 2, y(r) + rh / 2 + 0.5);
          }
        }
        ctx.globalAlpha = 1;
      }
      if (editable && (selCols.size || selRows.size)) {
        ctx.fillStyle = colours.selected;
        for (let c = c0; c < c1; c++) if (selCols.has(c)) ctx.fillRect(x(c), tracksH, cw, h - tracksH);
        for (let r = r0; r < r1; r++) if (selRows.has(r)) ctx.fillRect(namesW, y(r), w - namesW, rh);
      }
      if (boxed.length) { // the same site in every viewer: a tinted, outlined column
        const bottom = Math.min(h, y(n));
        for (const c of boxed) {
          if (c < c0 || c >= c1) continue;
          const bw = Math.max(cw, 2);
          ctx.fillStyle = "rgba(224, 130, 20, 0.16)";
          ctx.fillRect(x(c), tracksH, bw, bottom - tracksH);
          ctx.strokeStyle = "#e08214";
          ctx.lineWidth = 2;
          ctx.strokeRect(x(c) + 1, tracksH + 1, bw - 2, bottom - tracksH - 2);
        }
      }
      if (highlighted >= r0 && highlighted < r1) {
        ctx.strokeStyle = "#d55e00";
        ctx.lineWidth = 2;
        ctx.strokeRect(namesW + 1, y(highlighted) + 1, w - namesW - 2, rh - 2);
      }
      ctx.restore();

      // Tracks and ruler (fixed at the top, scrolled sideways with the grid).
      ctx.save();
      ctx.fillStyle = colours.bg;
      ctx.fillRect(0, 0, w, tracksH);
      ctx.beginPath();
      ctx.rect(namesW, 0, w - namesW, tracksH);
      ctx.clip();
      ctx.font = `10px ${MONO}`;
      ctx.textAlign = "left";
      ctx.textBaseline = "top";
      // Below a pixel per column (zoomed right out), labels at least ~48 px apart.
      const nice = [100, 200, 500, 1000, 2000, 5000, 10000];
      const step = cw >= 9 ? 10 : cw >= 5 ? 20 : cw >= 1 ? 50 : nice.find((k) => k * cw >= 48) || 20000;
      if (editable && selCols.size) {
        ctx.fillStyle = SELECTED_RULER;
        for (let c = c0; c < c1; c++) if (selCols.has(c)) ctx.fillRect(x(c), 0, cw, RULER);
      }
      ctx.fillStyle = colours.muted;
      for (let c = Math.floor(c0 / step) * step; c < c1; c += step) {
        if (c === 0) continue;
        ctx.fillRect(x(c - 1) + cw / 2, RULER - 5, 1, 5);
        ctx.fillText(String(c), x(c - 1) + cw / 2 + 2, 2);
      }
      data.tracks.forEach((track, t) => {
        const top = RULER + t * (TRACK + TRACK_GAP) + TRACK_GAP / 2;
        const [lo, hi] = ranges[t];
        const yv = (v) => top + TRACK - ((v - lo) / (hi - lo)) * TRACK;
        ctx.fillStyle = colours.band;
        ctx.fillRect(namesW, top, w - namesW, TRACK);
        if (lo < 0 && hi > 0) {
          ctx.fillStyle = colours.line;
          ctx.fillRect(namesW, yv(0), w - namesW, 1);
        }
        if (track.series.length === 1) {
          ctx.fillStyle = SERIES[0];
          const base = yv(Math.max(lo, Math.min(hi, 0)));
          const values = track.series[0].values;
          for (let c = c0; c < c1; c++) {
            const v = values[c];
            if (v === null || v === undefined) continue;
            const yy = yv(Math.max(lo, Math.min(hi, v)));
            ctx.fillRect(x(c) + (cw > 3 ? 0.5 : 0), Math.min(yy, base), Math.max(1, cw - (cw > 3 ? 1 : 0)), Math.max(1, Math.abs(base - yy)));
          }
        } else {
          track.series.forEach((series, k) => {
            ctx.strokeStyle = SERIES[k % SERIES.length];
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            let open = false;
            for (let c = Math.max(0, c0 - 1); c < Math.min(width, c1 + 1); c++) {
              const v = series.values[c];
              if (v === null || v === undefined) { open = false; continue; }
              const px = x(c) + cw / 2, py = yv(Math.max(lo, Math.min(hi, v)));
              if (open) ctx.lineTo(px, py); else ctx.moveTo(px, py);
              open = true;
            }
            ctx.stroke();
          });
        }
        if (track.text && cell >= 9) {
          ctx.fillStyle = colours.fg;
          ctx.textAlign = "center";
          ctx.font = `11px ${MONO}`;
          for (let c = c0; c < c1; c++) ctx.fillText(track.text[c] || "", x(c) + cw / 2, top + 2);
          ctx.textAlign = "left";
          ctx.font = `10px ${MONO}`;
        }
      });
      ctx.restore();

      // Names (fixed at the left) and the track labels above them.
      ctx.save();
      ctx.fillStyle = colours.bg;
      ctx.fillRect(0, 0, namesW, h);
      ctx.beginPath();
      ctx.rect(0, tracksH, namesW, h - tracksH);
      ctx.clip();
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.font = `${Math.min(12, Math.max(8, rh - 3))}px system-ui, sans-serif`;
      if (editable && selRows.size) {
        ctx.fillStyle = colours.selected;
        for (let r = r0; r < r1; r++) if (selRows.has(r)) ctx.fillRect(0, y(r), namesW, rh);
      }
      if (rh >= 8) {
        for (let r = r0; r < r1; r++) {
          // Members are designed sequences, not reconstructions: upright, and dimmed with
          // their row when the budget left them out.
          const italic = data.ancestor[r] && !member[r];
          ctx.fillStyle = italic || dim[r] || blank[r] ? colours.muted : colours.fg;
          ctx.font = `${italic ? "italic " : ""}${r === highlighted ? "bold " : ""}${Math.min(12, rh - 3)}px system-ui, sans-serif`;
          ctx.fillText(fit(data.names[r], namesW - 14), 6, y(r) + rh / 2);
        }
      }
      ctx.restore();
      ctx.save();
      ctx.font = "11px system-ui, sans-serif";
      ctx.textBaseline = "top";
      ctx.fillStyle = colours.muted;
      if (editable) { // the names' header: click to sort
        ctx.fillStyle = colours.fg;
        ctx.fillText(`Name ${data.sort === "name" ? "▲" : data.sort === "-name" ? "▼" : "⇅"}`, 6, 3);
      } else ctx.fillText("Column", 6, 3);
      data.tracks.forEach((track, t) => {
        const top = RULER + t * (TRACK + TRACK_GAP) + TRACK_GAP / 2;
        const [lo, hi] = ranges[t];
        ctx.fillStyle = colours.fg;
        ctx.textAlign = "left";
        ctx.fillText(fit(track.label, namesW - 40), 6, top + 2);
        ctx.fillStyle = colours.muted;
        ctx.textAlign = "right";
        ctx.fillText(fmt(hi), namesW - 6, top);
        ctx.fillText(fmt(lo), namesW - 6, top + TRACK - 11);
        if (track.series.length > 1) {
          ctx.textAlign = "left";
          track.series.slice(0, 3).forEach((s, k) => {
            ctx.fillStyle = SERIES[k % SERIES.length];
            ctx.fillText(fit(`— ${s.label}`, namesW - 44), 10, top + 14 + k * 10);
          });
        }
      });
      ctx.fillStyle = colours.line;
      ctx.fillRect(namesW - 1, 0, 1, h);
      ctx.fillRect(0, tracksH - 1, w, 1);
      ctx.restore();
      placeOverview();
    }
    function fit(text, room) {
      if (ctx.measureText(text).width <= room) return text;
      let t = text;
      while (t.length > 1 && ctx.measureText(`${t}…`).width > room) t = t.slice(0, -1);
      return `${t}…`;
    }

    // Where the pointer is: a cell, a track value, or the ruler.
    function hit(event) {
      const box = scroller.getBoundingClientRect();
      const px = event.clientX - box.left, py = event.clientY - box.top;
      if (px > scroller.clientWidth || py > scroller.clientHeight) return null;
      if (px < namesW) { // the editor's names and their header
        if (!editable || px < 0 || py < 0) return null;
        if (py < RULER) return { header: true, px, py };
        if (py < tracksH) return null;
        const r = Math.floor((py - tracksH + scroller.scrollTop) / rh);
        return r >= 0 && r < n ? { r, name: true, px, py } : null;
      }
      const c = Math.floor((px - namesW + scroller.scrollLeft) / cw);
      if (c < 0 || c >= width) return null;
      if (py < tracksH) {
        const t = Math.floor((py - RULER) / (TRACK + TRACK_GAP));
        return { c, track: py >= RULER && t < data.tracks.length ? t : -1, px, py };
      }
      const r = Math.floor((py - tracksH + scroller.scrollTop) / rh);
      return r >= 0 && r < n ? { c, r, px, py } : null;
    }
    function residueNumber(r, c) {
      if (isGap(site(r, c))) return null;
      let k = 0;
      for (let i = 0; i <= c; i++) if (!isGap(site(r, i))) k++;
      return k;
    }
    function describe(h) {
      if (h.header) return data.sort === "name" ? "Sorted A to Z: click for Z to A" : "Click to sort by name";
      if (h.name) {
        return `${data.names[h.r]}\nClick to select; Shift-click for a range, ${
          navigator.platform.startsWith("Mac") ? "⌘" : "Ctrl"}-click to add`;
      }
      const origin = data.origin && data.origin.length ? data.origin[h.c] : 0;
      const lines = [`Column ${h.c + 1}${origin && origin !== h.c + 1 ? ` (column ${origin} before the cuts)` : ""}`];
      if (h.r !== undefined && blank[h.r]) return "This combination's tree has no such ancestor";
      if (h.r !== undefined) {
        const s = site(h.r, h.c), number = residueNumber(h.r, h.c);
        const kind = member[h.r] ? (dim[h.r] ? " (member, left out by the budget)" : " (member)")
          : data.ancestor[h.r] ? " (ancestor)" : "";
        lines.unshift(data.names[h.r] + kind);
        let what = isGap(s) ? "gap" : `${s}${number ? ` (residue ${number})` : ""}`;
        if (data.codon && !isGap(s)) what += ` → ${aminoAt(h.r, h.c)}`;
        lines.push(what);
        const pct = ppAt(h.r, h.c);
        if (pct !== null) {
          const alts = data.alts && data.alts[h.r] && data.alts[h.r][h.c];
          const more = alts ? alts.split(" ").map((x) => `${x[0]} ${(Number(x.slice(1)) / 100).toFixed(2)}`) : [];
          lines.push([`${s} ${(pct / 100).toFixed(2)}`, ...more].join(", "));
          if (faded[h.r] && faded[h.r].has(h.c)) lines.push("aligners disagree on this column");
          if (marked[h.r] && marked[h.r].has(h.c) && member[h.r]) lines.push("differs from the ML sequence");
        }
        const note = notes[h.r] && notes[h.r][h.c];
        if (note) lines.push(note);
      } else if (h.track >= 0) {
        const track = data.tracks[h.track];
        lines.unshift(track.label);
        track.series.forEach((s) => {
          const v = s.values[h.c];
          const shown = v === null || v === undefined ? "—" : fmt(v);
          lines.push(track.series.length > 1 ? `${s.label}: ${shown}` : shown);
        });
        if (track.text && track.text[h.c]) lines.push(`Consensus ${track.text[h.c]}`);
      }
      return lines.join("\n");
    }
    const report = (column) => {
      if (!data.hoverEvents || column === hovered) return;
      hovered = column;
      wrap.dispatchEvent(new CustomEvent("aln:hover", { bubbles: true, detail: { column } }));
    };
    scroller.addEventListener("mousemove", (event) => {
      const h = hit(event);
      report(h && h.c !== undefined ? h.c : -1);
      if (!h) { tip.hidden = true; scroller.style.cursor = ""; return; }
      tip.textContent = describe(h);
      tip.hidden = false;
      const left = Math.min(h.px + 14, scroller.clientWidth - tip.offsetWidth - 4);
      tip.style.left = `${Math.max(0, left)}px`;
      tip.style.top = `${h.py + 16}px`;
      scroller.style.cursor = h.header || h.name || data.ancestorViewer ? "pointer"
        : editable ? "cell" : "";
    });
    scroller.addEventListener("mouseleave", () => { tip.hidden = true; report(-1); });
    scroller.addEventListener("click", (event) => {
      const h = hit(event);
      if (h && data.rowEvents && h.r !== undefined) {
        wrap.dispatchEvent(new CustomEvent("aln:click", { bubbles: true,
          detail: { row: h.r, column: h.c, name: data.names[h.r], ancestor: data.ancestor[h.r] } }));
        return;
      }
      if (!h || !data.ancestorViewer) return;
      let url = data.ancestorViewer;
      if (h.r !== undefined && data.ancestor[h.r]) url += `&node=${encodeURIComponent(data.names[h.r])}`;
      window.location.href = `${url}#site-${h.c + 1}`;
    });

    // The editor's selection: columns by dragging across the ruler or grid, sequences by
    // clicking or dragging over their names. Selecting one kind clears the other, so a
    // Delete key has one meaning.
    function emitSelection() {
      const sorted = (set) => [...set].sort((a, b) => a - b);
      wrap.dispatchEvent(new CustomEvent("aln:select", { bubbles: true,
        detail: { columns: sorted(selCols), rows: sorted(selRows) } }));
    }
    const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
    function dragTo(clientX, clientY) {
      const box = scroller.getBoundingClientRect();
      let at;
      if (drag.kind === "columns") {
        at = clamp(Math.floor((clientX - box.left - namesW + scroller.scrollLeft) / cw), 0, width - 1);
      } else {
        at = clamp(Math.floor((clientY - box.top - tracksH + scroller.scrollTop) / rh), 0, n - 1);
      }
      const chosen = new Set(drag.base);
      for (let i = Math.min(drag.anchor, at); i <= Math.max(drag.anchor, at); i++) chosen.add(i);
      if (drag.kind === "columns") selCols = chosen; else selRows = chosen;
      redraw();
      emitSelection();
    }
    function autoScroll() { // while dragging past an edge of the view
      if (!drag) return;
      const box = scroller.getBoundingClientRect();
      const { clientX, clientY } = drag.pointer;
      let moved = false;
      if (drag.kind === "columns") {
        const left = box.left + namesW, right = box.left + scroller.clientWidth;
        const step = clientX > right - 12 ? 1 : clientX < left + 12 ? -1 : 0;
        if (step) { scroller.scrollLeft += step * Math.max(cw, 8); moved = true; }
      } else {
        const top = box.top + tracksH, bottom = box.top + scroller.clientHeight;
        const step = clientY > bottom - 8 ? 1 : clientY < top + 8 ? -1 : 0;
        if (step) { scroller.scrollTop += step * rh; moved = true; }
      }
      if (moved) dragTo(clientX, clientY);
      drag.timer = requestAnimationFrame(autoScroll);
    }
    if (editable) {
      scroller.addEventListener("mousedown", (event) => {
        if (event.button !== 0) return;
        const h = hit(event);
        if (!h) return;
        if (h.header) {
          wrap.dispatchEvent(new CustomEvent("aln:sort", { bubbles: true }));
          return;
        }
        if (h.c === undefined && !h.name) return; // a track's label area
        event.preventDefault();
        tip.hidden = true;
        const add = event.metaKey || event.ctrlKey;
        const kind = h.name ? "rows" : "columns";
        const at = h.name ? h.r : h.c;
        const existing = kind === "rows" ? selRows : selCols;
        if (kind === "rows") selCols = new Set(); else selRows = new Set();
        if (add && existing.has(at) && !event.shiftKey) { // Ctrl/Cmd-click a chosen one: unchoose
          existing.delete(at);
          redraw();
          emitSelection();
          return;
        }
        let anchor = kind === "rows" ? anchorRow : anchorCol;
        if (!event.shiftKey || anchor < 0) anchor = at;
        if (kind === "rows") anchorRow = anchor; else anchorCol = anchor;
        drag = { kind, anchor, base: add ? new Set(existing) : new Set(),
                 pointer: { clientX: event.clientX, clientY: event.clientY } };
        dragTo(event.clientX, event.clientY);
        drag.timer = requestAnimationFrame(autoScroll);
      });
      function onMove(event) {
        if (!wrap.isConnected) { // replaced by an HTMX swap
          window.removeEventListener("mousemove", onMove);
          window.removeEventListener("mouseup", onUp);
          return;
        }
        if (!drag) return;
        drag.pointer = { clientX: event.clientX, clientY: event.clientY };
        dragTo(event.clientX, event.clientY);
      }
      function onUp() {
        if (!drag) return;
        cancelAnimationFrame(drag.timer);
        drag = null;
      }
      window.addEventListener("mousemove", onMove);
      window.addEventListener("mouseup", onUp);
    }

    let pending = false;
    const redraw = () => {
      if (pending) return;
      pending = true;
      requestAnimationFrame(() => { pending = false; draw(); });
    };
    scroller.addEventListener("scroll", redraw, { passive: true });
    new ResizeObserver(() => {
      if (fitting) setCell("fit");
      overviewStale = true;
      sizeCanvas();
      draw();
    }).observe(wrap);
    const onResize = () => {
      if (!wrap.isConnected) { window.removeEventListener("resize", onResize); return; }
      if (fitting) setCell("fit");
      sizeCanvas();
      draw();
    };
    window.addEventListener("resize", onResize); // "fit" follows the window's height
    // Pages that embed the viewer change its palette, threshold or zoom in place. A zoom
    // keeps the column at the middle of the view where it was.
    wrap.addEventListener("aln:update", (event) => {
      const patch = event.detail || {};
      if (patch.scheme) data.scheme = patch.scheme;  // the column cache is scheme-free
      if (patch.palette) { data.palette = patch.palette; tints.clear(); }
      if (patch.scheme || patch.palette) overviewStale = true;
      if (patch.threshold !== undefined) data.threshold = patch.threshold;
      if (patch.columns !== undefined) boxed = patch.columns;
      if (patch.clear) {
        selCols = new Set();
        selRows = new Set();
        emitSelection();
      }
      if (patch.highlight !== undefined) {
        highlighted = (data.anchors || []).indexOf(patch.highlight);
        if (highlighted >= 0) { // into view, a few rows from the top
          const view = scroller.clientHeight - tracksH;
          const top = highlighted * rh, now = scroller.scrollTop;
          if (top < now || top + rh > now + view) scroller.scrollTop = Math.max(0, top - 2 * rh);
        }
      }
      if (patch.cell !== undefined) {
        const view = scroller.clientWidth - namesW;
        const middle = (scroller.scrollLeft + view / 2) / cw;
        const rows = scroller.clientHeight - tracksH;
        const centre = (scroller.scrollTop + rows / 2) / rh;
        setCell(patch.cell);
        sizeCanvas();
        scroller.scrollLeft = Math.max(0, middle * cw - view / 2);
        const after = scroller.clientHeight - tracksH;
        if (data.rowFromCell) scroller.scrollTop = Math.max(0, centre * rh - after / 2);
      }
      draw();
    });
    sizeCanvas();
    if (saved) { // the view this one replaced: same first column and top sequence
      scroller.scrollLeft = saved.column * cw;
      const top = data.names.indexOf(saved.name);
      scroller.scrollTop = top >= 0 ? (top + saved.offset) * rh : 0;
      // The editor's selection survives a change of colours, tracks or height, not an edit.
      if (editable && saved.state === data.stateKey) {
        selCols = new Set(saved.columns.filter((c) => c < width));
        selRows = new Set(saved.rows.filter((r) => r < n));
        if (selCols.size || selRows.size) setTimeout(emitSelection);
      }
      saved = null;
    }
    draw();
    current = () => {
      const row = Math.min(n - 1, Math.floor(scroller.scrollTop / rh));
      return { column: scroller.scrollLeft / cw, name: data.names[row],
               offset: scroller.scrollTop / rh - row,
               state: data.stateKey, columns: [...selCols], rows: [...selRows] };
    };

    // The viewer page's zoom slider (the library page drives its own, in library.js):
    // in place, without asking the server, on the library's exponential scale from
    // MIN_CELL to MAX_CELL pixels per column, so each step is the same proportional
    // change. Its value rides in the form's hidden field and the address, so the next
    // swap or a reload keeps it.
    const zoom = root.querySelector(".aln-zoom");
    if (zoom) {
      const field = zoom.querySelector('input[name="zoom"]');
      const slider = zoom.querySelector(".zoom-slider");
      const cellAt = (v) => MIN_CELL * (MAX_CELL / MIN_CELL) ** (Number(v) / 100);
      const sliderAt = (c) => 100 * Math.log(c / MIN_CELL) / Math.log(MAX_CELL / MIN_CELL);
      const place = () => { slider.value = String(Math.round(sliderAt(cell))); };
      let queued = null;
      const apply = (value) => {
        const first = queued === null;
        queued = value;
        if (!first) return;
        requestAnimationFrame(() => {
          const v = queued;
          queued = null;
          wrap.dispatchEvent(new CustomEvent("aln:update", { detail: { cell: v } }));
          field.value = v === "fit" ? "fit" : v.toFixed(2);
          const url = new URL(window.location.href);
          url.searchParams.set("zoom", field.value);
          history.replaceState(history.state, "", url);
          if (v === "fit") place();
        });
      };
      // A range input fires "change", which would send the whole form to the server.
      slider.addEventListener("change", (event) => event.stopPropagation());
      slider.addEventListener("input", () => apply(cellAt(slider.value)));
      zoom.querySelector(".zoom-fit").addEventListener("click", () => apply("fit"));
      place();
    }

    const png = root.querySelector(".aln-png");
    if (png) {
      png.addEventListener("click", (event) => {
        event.preventDefault();
        const link = document.createElement("a");
        link.download = `${png.dataset.stem}.png`;
        link.href = canvas.toDataURL("image/png");
        link.click();
      });
    }
  }

  const start = () => {
    const roots = document.querySelectorAll(".aln-root");
    if (roots.length) roots.forEach(init);
    else init(document);
  };
  // The schemes that colour a residue on its own, for pages that draw sequences outside a
  // viewer canvas (the library's site maps). Column-dependent schemes — Clustal X,
  // identity, differences — are deliberately absent: they need a column to judge.
  window.alignmentResidueColours = {
    chemistry: (r) => CHEMISTRY[r] || null,
    hydropathy,
    structure,
    none: () => null,
  };

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();
  document.addEventListener("htmx:afterSettle", start);
})();
