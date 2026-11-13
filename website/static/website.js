// phylobot.com only (website/). The app's own scripts do the backdrop and the copy buttons;
// this adds what a landing page with no server behind it needs.
(() => {
  const root = document.documentElement;
  const still = matchMedia("(prefers-reduced-motion: reduce)").matches;
  const store = {
    get(key) { try { return localStorage.getItem(key); } catch { return null; } },
    set(key, value) { try { localStorage.setItem(key, value); } catch { /* private window */ } },
  };

  // ---- skins: the menu sets the attribute (every skin's stylesheet is loaded), and the
  // choice is kept in this browser; the script in <head> applies it before the first paint.
  const picker = document.querySelector("[data-skin-picker]");
  if (picker) {
    picker.value = root.dataset.skin;
    picker.addEventListener("change", () => {
      root.dataset.skin = picker.value;
      root.dataset.backdrop = picker.selectedOptions[0].dataset.backdrop;   // backdrop.js watches
      store.set("phylobot_skin", picker.value);
    });
  }

  // ---- the section links: shown once the title has scrolled away, marking where you are.
  const nav = document.querySelector("[data-site-nav]");
  const hero = document.querySelector(".title-page");
  if (nav && hero && "IntersectionObserver" in window) {
    new IntersectionObserver(([entry]) => { nav.hidden = entry.isIntersecting; }, { rootMargin: "-80px 0px 0px 0px" })
      .observe(hero);
    const links = new Map();
    for (const link of nav.querySelectorAll("a[href^='#']")) {
      const target = document.getElementById(link.getAttribute("href").slice(1));
      if (target && target !== document.getElementById("top")) links.set(target, link);
    }
    const marker = new IntersectionObserver((entries) => {
      for (const entry of entries) {
        if (!entry.isIntersecting) continue;
        for (const link of links.values()) link.removeAttribute("aria-current");
        links.get(entry.target).setAttribute("aria-current", "true");
      }
    }, { rootMargin: "-45% 0px -50% 0px" });
    for (const target of links.keys()) marker.observe(target);
  }

  // ---- the pipeline, played -------------------------------------------------------
  //
  // An invented family of a hundred-odd sequences on its tree, taken through the five steps
  // the list beside it names: unaligned, aligned, the tree drawn, an ancestor at each node,
  // and the changes along one path. Real families are this size or larger, so the picture is
  // too: every sequence is a row, most of them a hair's height, and six in the middle are
  // enlarged as under a lens (the rows beside them swell towards it) so the residues can be
  // read. Only those six and their five ancestors are worked out by hand; the rest are
  // generated around them, the same on every visit. The ancestors are consistent with the tips (each branch's
  // changes are listed in PATH for the one path shown); the support values are invented.
  function pipeline(demo) {
    const stage = demo.querySelector("[data-demo-stage]");
    const buttons = [...demo.querySelectorAll("[data-step]")];
    const TIPS = [
      { name: "seq 49", seq: "MKVISEGDWKTA" },
      { name: "seq 50", seq: "MKVI--GDWRTS" },
      { name: "seq 51", seq: "MRALSDGNWRTA" },
      { name: "seq 52", seq: "MRVLSDGNWRSA" },
      { name: "seq 53", seq: "MRVLTDGDFRTA" },
      { name: "seq 54", seq: "MRVLSDGDFRTG" },
    ];
    const NODES = {
      n1: { label: "Node 52", x: 62, seq: "MKVLSEGDWRTA", of: ["n2", "n3"], support: [1, .62, .97, .71, .93, .58, 1, .83, .88, .95, .97, .9] },
      n2: { label: "Node 53", x: 160, seq: "MKVISEGDWRTA", of: [0, 1], support: [1, .95, .99, .93, .9, .88, 1, .99, .99, .84, .98, .8] },
      n3: { label: "Node 55", x: 112, seq: "MRVLSDGDWRTA", of: ["n4", "n5"], support: [1, .96, .97, .98, .95, .94, 1, .74, .79, .99, .97, .96] },
      n4: { label: "Node 56", x: 176, seq: "MRVLSDGNWRTA", of: [2, 3], support: [1, .99, .86, .99, .99, .98, 1, .94, .99, .99, .87, .97] },
      n5: { label: "Node 59", x: 176, seq: "MRVLSDGDFRTA", of: [4, 5], support: [1, .99, .99, .99, .9, .98, 1, .97, .95, .99, .98, .86] },
    };
    // Node 52 to seq 53: the branches on the way, and what changed on each.
    const PATH = [["n1", "n3", "K2R E6D"], ["n3", "n5", "W9F"], ["n5", 4, "S5T"]];
    const CHANGED = [1, 5, 8, 4];                       // those changes' columns, from 0
    const KIND = { M: "h", V: "h", L: "h", I: "h", A: "h", W: "h", F: "h", S: "p", T: "p", N: "p", G: "p", Q: "p",
                   K: "b", R: "b", H: "b", D: "a", E: "a" };
    const LETTERS = Object.keys(KIND);
    const PITCH = 32, TOP = 4, TIP_X = 250, TILE = 30, FIRST = 58;
    // The rest of the family: how many rows above and below the six, and how tall a row is
    // by its distance from them (the lens), down to a hair.
    const ABOVE = 48, BELOW = 52, LENS = [13, 10, 7.5, 5.5, 4, 3, 2.2, 1.7], THIN = 1.3;
    const DWELL = { 1: 2600, 2: 2800, 3: 3200, 4: 5400, 5: 5400 };
    const SVG = "http://www.w3.org/2000/svg";
    const el = (tag, cls, text) => {
      const node = document.createElement(tag);
      if (cls) node.className = cls;
      if (text !== undefined) node.textContent = text;
      return node;
    };
    const svg = (tag, attrs) => {
      const node = document.createElementNS(SVG, tag);
      for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
      return node;
    };

    // The same "random" family on every visit.
    let seed = 20161107;
    const rand = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
    const relative = (base) => {
      const seq = [...base].map((letter) => (rand() < 0.13 ? LETTERS[Math.floor(rand() * LETTERS.length)] : letter));
      if (rand() < 0.18) { const at = 3 + Math.floor(rand() * 6); seq[at] = "-"; if (rand() < 0.5) seq[at + 1] = "-"; }
      return seq.join("");
    };
    const context = (count, base, top, towardsLens) => {
      const rows = [];
      let y = top;
      for (let i = 0; i < count; i++) {
        const h = LENS[towardsLens ? count - 1 - i : i] ?? THIN;
        rows.push({ y, h, seq: relative(base) });
        y += h;
      }
      return { rows, bottom: y };
    };
    const above = context(ABOVE, "MKVISEGDWRTA", TOP, true);
    const lensTop = above.bottom + 5;
    const below = context(BELOW, "MRVLSDGDWRTA", lensTop + TIPS.length * PITCH + 1, false);
    const HEIGHT = Math.ceil(below.bottom + 18 + 28 + 4);

    const canvas = el("div", "demo-canvas");
    canvas.style.height = HEIGHT + "px";
    const tree = svg("svg", { viewBox: `0 0 690 ${HEIGHT}`, width: 690, height: HEIGHT });
    canvas.append(tree);

    // Where everything sits: tips down the rows, each node midway between its children.
    const tipY = (i) => lensTop + i * PITCH + 14;
    const place = (id) => {
      if (typeof id === "number") return { x: TIP_X, y: tipY(id) };
      const node = NODES[id];
      if (node.y === undefined) node.y = (place(node.of[0]).y + place(node.of[1]).y) / 2;
      return node;
    };
    const depth = (id, d = 0) => { NODES[id].depth = d; NODES[id].of.forEach((c) => typeof c === "string" && depth(c, d + 1)); };
    depth("n1", 1);
    const onPath = new Set(PATH.map(([from, to]) => from + ">" + to));
    // The rest of the tree: each run of rows splits in two, again and again, down to the tips.
    const grow = (rows, lo, hi, level) => {
      if (lo === hi) return { x: TIP_X, y: rows[lo].y + rows[lo].h / 2 };
      const split = lo + 1 + Math.floor(((rand() + rand()) / 2) * (hi - lo));
      const kids = [grow(rows, lo, split - 1, level + 1), grow(rows, split, hi, level + 1)];
      const node = { x: Math.min(234, 26 + level * 13 + rand() * 7), y: (kids[0].y + kids[1].y) / 2 };
      for (const kid of kids) {
        const branch = svg("path", { class: "branch ctx", pathLength: 1, d: `M${node.x},${node.y} V${kid.y} H${kid.x}` });
        branch.style.setProperty("--delay", Math.min(level, 8) * 0.12 + "s");
        tree.append(branch);
      }
      return node;
    };
    const root = { x: 8, y: place("n1").y };
    const clades = [grow(above.rows, 0, ABOVE - 1, 0), grow(below.rows, 0, BELOW - 1, 0)];
    tree.append(svg("path", { class: "branch ctx", pathLength: 1, d: `M0,${root.y} H${root.x}` }));
    for (const clade of clades) {
      tree.append(svg("path", { class: "branch ctx", pathLength: 1, d: `M${root.x},${root.y} V${clade.y} H${clade.x}` }));
    }
    tree.append(svg("path", { class: "branch", pathLength: 1, d: `M${root.x},${root.y} H${NODES.n1.x}` }));
    for (const [id, node] of Object.entries(NODES)) {
      for (const child of node.of) {
        const to = place(child);
        const branch = svg("path", { class: "branch" + (onPath.has(id + ">" + child) ? " on-path" : ""), pathLength: 1,
                                     d: `M${node.x},${node.y} V${to.y} H${to.x}` });
        branch.style.setProperty("--delay", (node.depth + 1) * 0.22 + "s");
        tree.append(branch);
      }
    }
    TIPS.forEach((_, i) => tree.append(svg("circle", { class: "tip-dot", cx: TIP_X, cy: tipY(i), r: 2.5 })));
    const halos = {}, dots = {}, folds = {};
    for (const [id, node] of Object.entries(NODES)) {
      // A little ribbon beside the node: a helix's coils running into a strand's arrow.
      folds[id] = svg("path", { class: "fold-glyph", d: `M${node.x - 30},${node.y - 6} c6,-11 -3,-11 3.6,0 c6,-11 -3,-11 3.6,0 c6,-11 -3,-11 3.6,0 l8,0 m-3.2,-3.2 l3.2,3.2 l-3.2,3.2` });
      tree.append(folds[id]);
      halos[id] = svg("circle", { class: "halo", cx: node.x, cy: node.y, r: 10 });
      dots[id] = svg("circle", { class: "node", cx: node.x, cy: node.y, r: 5 });
      dots[id].addEventListener("click", () => { stop(); if (step < 4) show(4); select(id, true); });
      tree.append(halos[id], dots[id]);
    }

    // The sequences: a tile per column, which knows where it sits unaligned and aligned.
    const tile = (letter, column, packed) => {
      const t = el("span", "res " + (letter === "-" ? "gap" : KIND[letter]), letter === "-" ? "–" : letter);
      t.style.setProperty("--aligned", FIRST + column * TILE + "px");
      t.style.setProperty("--packed", FIRST + packed * TILE + "px");
      return t;
    };
    const tipRows = TIPS.map((tip, i) => {
      const row = el("div", "demo-row tip" + (i === PATH.at(-1)[1] ? " on-path" : ""));
      row.style.top = lensTop + i * PITCH + "px";
      row.append(el("span", "name", tip.name));
      let seen = 0;
      const tiles = [...tip.seq].map((letter, column) => tile(letter, column, letter === "-" ? column : seen++));
      row.append(...tiles);
      canvas.append(row);
      return tiles;
    });
    // The other rows: the same tiles, a hair tall; the ones nearest the lens show letters.
    for (const { y, h, seq } of [...above.rows, ...below.rows]) {
      const row = el("div", "demo-row tip mini");
      row.style.top = y + "px";
      row.style.setProperty("--h", Math.max(0.9, h - (h > 3 ? 1.5 : 0.35)) + "px");
      if (h >= 9) row.style.setProperty("--fs", h * 0.62 + "px");
      let seen = 0;
      row.append(...[...seq].map((letter, column) => {
        const t = tile(letter, column, letter === "-" ? column : seen++);
        if (h < 9) t.textContent = "";
        return t;
      }));
      canvas.append(row);
    }
    const ancestor = el("div", "demo-row ancestor");
    ancestor.style.top = below.bottom + 18 + "px";
    const ancestorName = el("span", "name");
    const ancestorTiles = [];
    ancestor.append(ancestorName);
    for (let column = 0; column < TIPS[0].seq.length; column++) {
      const t = el("span", "res");
      t.style.transform = `translateX(${FIRST + column * TILE}px)`;
      ancestorTiles.push(t);
      ancestor.append(t);
    }
    canvas.append(ancestor);
    for (const [from, to, changes] of PATH) {
      const a = place(from), b = place(to);
      const chip = el("span", "demo-chip", changes);
      chip.style.left = (a.x + b.x) / 2 + "px";
      chip.style.top = b.y + "px";
      canvas.append(chip);
    }
    stage.append(canvas);
    // A narrow window gets the whole picture, smaller, rather than half of it and a scrollbar.
    const fit = () => {
      const room = stage.clientWidth - parseFloat(getComputedStyle(stage).paddingLeft) * 2;
      const scale = Math.min(1, room / 690);
      canvas.style.transform = scale < 1 ? `scale(${scale})` : "";
      canvas.style.marginBottom = scale < 1 ? -HEIGHT * (1 - scale) + "px" : "";
      canvas.style.marginLeft = scale < 1 ? "0" : "";
    };
    fit();
    window.addEventListener("resize", fit);

    // ---- one ancestor in the row beneath: its residues settle one after another ----
    let selected = "", settling = 0;
    function select(id, animate) {
      const node = NODES[id], run = ++settling;
      selected = id;
      for (const key of Object.keys(NODES)) {
        dots[key].classList.toggle("selected", key === id);
        halos[key].classList.toggle("selected", key === id);
        folds[key].classList.toggle("selected", key === id);
      }
      ancestorName.textContent = node.label;
      ancestorTiles.forEach((t, column) => {
        const put = (letter) => {
          t.textContent = letter;
          t.classList.remove("h", "p", "b", "a");
          t.classList.add(KIND[letter]);
        };
        const settle = () => {
          put(node.seq[column]);
          t.style.setProperty("--pp", node.support[column] * 100 + "%");
        };
        if (!animate || still) return settle();
        const until = performance.now() + 120 + column * 45;
        (function flicker() {
          if (run !== settling) return;                       // another node was chosen meanwhile
          if (performance.now() >= until) return settle();
          put(LETTERS[Math.floor(Math.random() * LETTERS.length)]);
          setTimeout(flicker, 60);
        })();
      });
    }

    // ---- the steps ----
    let step = 0, playing = false, timer = 0, tours = [], seen = false, held = false;
    function show(next) {
      const entering = next !== step;
      step = next;
      canvas.dataset.step = step;
      demo.style.setProperty("--dwell", DWELL[step] + "ms");
      buttons.forEach((b) => (Number(b.dataset.step) === step ? b.setAttribute("aria-current", "step") : b.removeAttribute("aria-current")));
      tours.forEach(clearTimeout);
      tours = [];
      const path = step === 5;
      ancestorTiles.forEach((t, column) => t.classList.toggle("changed", path && CHANGED.includes(column)));
      tipRows[PATH.at(-1)[1]].forEach((t, column) => t.classList.toggle("changed", path && CHANGED.includes(column)));
      if (path) select("n1", entering && selected !== "n1");
      else if (step === 4 && entering) {
        select("n1", true);
        // While it plays by itself, look at an ancestor further down the tree every so often.
        if (playing) tours = [setTimeout(() => select("n3", true), 1800), setTimeout(() => select("n5", true), 3600)];
      }
    }
    function schedule() {
      clearTimeout(timer);
      if (playing) timer = setTimeout(() => { show(step % 5 + 1); schedule(); }, DWELL[step]);
    }
    const toggle = el("button", "link");
    toggle.type = "button";
    function play() { playing = true; demo.classList.add("playing"); toggle.textContent = "Pause"; schedule(); }
    function stop() { playing = false; demo.classList.remove("playing"); toggle.textContent = "Play"; clearTimeout(timer); tours.forEach(clearTimeout); }
    toggle.addEventListener("click", () => {
      held = playing;                                   // paused by hand: scrolling back must not restart it
      if (playing) stop();
      else { show(step % 5 + 1); play(); }
    });
    demo.querySelector(".demo-note").append(toggle);
    buttons.forEach((b) => b.addEventListener("click", () => { held = true; stop(); show(Number(b.dataset.step)); }));

    // At rest the whole picture is there (the ancestors). It plays from the start the first
    // time it comes into view, and only while it is in view.
    show(4);
    stop();
    select("n1", false);
    if (still || !("IntersectionObserver" in window)) return;
    new IntersectionObserver(([entry]) => {
      if (!entry.isIntersecting) { if (playing) { stop(); held = false; } return; }
      if (held || playing) return;
      if (!seen) { seen = true; show(1); }
      play();
    }, { threshold: 0.45 }).observe(stage);
  }
  document.querySelectorAll("[data-demo]").forEach(pipeline);

  // ---- screenshots open at full size ----
  let zoom = null;
  document.addEventListener("click", (event) => {
    const shot = event.target.closest(".shot img");
    if (!shot || !window.HTMLDialogElement) return;
    if (!zoom) {
      zoom = document.createElement("dialog");
      zoom.className = "zoom";
      zoom.append(document.createElement("img"));
      zoom.addEventListener("click", () => zoom.close());
      document.body.append(zoom);
    }
    zoom.firstChild.src = shot.currentSrc || shot.src;
    zoom.firstChild.alt = shot.alt;
    zoom.showModal();
  });

  // ---- the tutorial video: the player replaces the panel only when the reader asks for it,
  // so nothing is fetched from YouTube before that. Without this script the panel is a link.
  document.addEventListener("click", (event) => {
    const poster = event.target.closest(".video-poster");
    const video = poster && poster.closest("[data-video]");
    if (!video || !/^[A-Za-z0-9_-]{11}$/.test(video.dataset.video)) return;
    event.preventDefault();
    const frame = document.createElement("iframe");
    frame.src = "https://www.youtube-nocookie.com/embed/" + video.dataset.video + "?autoplay=1&rel=0";
    frame.title = video.dataset.title || "Tutorial video";
    frame.allow = "accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share";
    frame.allowFullscreen = true;
    frame.referrerPolicy = "strict-origin-when-cross-origin";
    poster.replaceWith(frame);
  });

  // ---- folded lists: a link to something inside one opens it, and one button opens every
  // citation, for the reader collecting them.
  const reveal = () => {
    const target = location.hash.length > 1 && document.getElementById(decodeURIComponent(location.hash.slice(1)));
    if (!target) return;
    let opened = false;
    for (let node = target; node; node = node.parentElement) {
      if (node.tagName === "DETAILS" && !node.open) { node.open = true; opened = true; }
    }
    if (opened) target.scrollIntoView();
  };
  window.addEventListener("hashchange", reveal);
  reveal();
  for (const button of document.querySelectorAll("[data-unfold]")) {
    const folds = () => [...document.querySelectorAll(button.dataset.unfold)];
    const label = () => {
      button.textContent = folds().every((fold) => fold.open) ? "Hide the citations" : "Show every citation";
    };
    button.hidden = false;
    button.addEventListener("click", () => {
      const open = !folds().every((fold) => fold.open);
      folds().forEach((fold) => { fold.open = open; });
      label();
    });
    document.addEventListener("toggle", label, true);
    label();
  }

  // ---- the PhyloBot '95 skin's furniture: a hit counter that counts this browser's own
  // visits, and a helper that stays dismissed.
  const visits = (Number(store.get("phylobot_visits")) || 0) + 1;
  store.set("phylobot_visits", String(visits));
  const counter = document.querySelector("[data-visits]");
  if (counter) {
    counter.querySelector(".digits").textContent = String(visits).padStart(6, "0");
    counter.hidden = false;
  }
  const helper = document.querySelector("[data-retro-helper]");
  if (helper && store.get("phylobot_helper") === "dismissed") helper.hidden = true;
  document.addEventListener("click", (event) => {
    if (!event.target.closest("[data-retro-dismiss]") || !helper) return;
    helper.hidden = true;
    store.set("phylobot_helper", "dismissed");
  });
})();
