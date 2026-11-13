// Pick a mutational path on the tree (_tree_viewer.html, _mutational_path.html): "Trace a
// path on the tree" starts picking; the first click names the ancestor to start from,
// hovering then previews the route to whatever is under the pointer, and the second click
// (a descendant, or an ancestor of the first: they are swapped) opens the path on its own
// page (phylobot.path_views). Escape or the button again cancels. The route comes from
// the branches' data-child / data-parent, so the server still checks what it is sent
// (views.path_context).
(() => {
  let picking = null;  // {figure, from, status} while picking

  const figureOf = () => document.querySelector("#tree-viewer figure.tree-figure");
  // Something on the tree that can end a path: an ancestral node (or a folded clade's
  // apex) by data-node, or an extant sequence by data-leaf.
  const pickable = (target) => {
    const el = target.closest && target.closest("a.tree-node[data-node], .tree-clade[data-node], .tree-leaf[data-leaf]");
    if (!el || !el.closest("figure.tree-figure")) return null;
    return { el, name: el.dataset.node || el.dataset.leaf, ancestor: !el.dataset.leaf };
  };
  const parents = (figure) => {
    const up = new Map();
    figure.querySelectorAll("g.tree-branch[data-child]").forEach((g) => {
      if (g.dataset.child) up.set(g.dataset.child, g.dataset.parent);
    });
    return up;
  };
  // The names from `top` down to `bottom`, oldest first, or null if `bottom` is not below.
  const route = (up, top, bottom) => {
    const names = [bottom];
    let at = bottom;
    while (up.has(at)) {
      at = up.get(at);
      names.unshift(at);
      if (at === top) return names;
    }
    return null;
  };
  // A route either way: [from, to] in tree order, or null when neither is below the other.
  const between = (up, a, b) => {
    const down = route(up, a, b);
    if (down) return { from: a, to: b, names: down };
    const back = route(up, b, a);
    return back ? { from: b, to: a, names: back } : null;
  };

  const say = (text) => { if (picking) picking.status.querySelector("span").textContent = text; };
  const clearPreview = (figure) => {
    figure.querySelectorAll(".path-preview-route").forEach((el) => el.remove());
    figure.querySelectorAll(".path-target").forEach((el) => el.classList.remove("path-target"));
  };
  const stop = () => {
    if (!picking) return;
    const { figure, status } = picking;
    clearPreview(figure);
    figure.classList.remove("path-picking");
    figure.querySelectorAll(".path-start").forEach((el) => el.classList.remove("path-start"));
    status.remove();
    document.querySelectorAll("[data-path-pick]").forEach((b) => b.setAttribute("aria-pressed", "false"));
    picking = null;
  };
  const start = () => {
    const figure = figureOf();
    if (!figure) return;
    const status = document.createElement("p");
    status.className = "path-pick-status";
    status.setAttribute("role", "status");
    status.innerHTML = '<span></span> <button type="button" class="link" data-path-pick-cancel>Cancel</button>';
    figure.querySelector(".plot-tools").after(status);
    figure.classList.add("path-picking");
    document.querySelectorAll("[data-path-pick]").forEach((b) => b.setAttribute("aria-pressed", "true"));
    picking = { figure, from: null, status, up: parents(figure) };
    say("Click the ancestor the path starts from.");
    figure.scrollIntoView({ block: "nearest", behavior: "smooth" });
  };

  // Open the path on its own page (phylobot.path_views): this tree's combination and gap
  // version, the two ends, and how the tree is drawn now ("back"), so the page can lead
  // back to it.
  const trace = (from, to) => {
    const form = document.getElementById("tree-viewer");
    if (!form || !form.dataset.pathUrl) return;
    const params = new URLSearchParams();
    for (const name of ["aligner", "combination", "gaps"]) {
      const field = form.elements[name];
      if (field && field.value) params.set(name, field.value);
    }
    params.set("path_from", from);
    params.set("path_to", to);
    const back = new URLSearchParams(window.location.search);
    for (const name of ["aligner", "combination", "gaps", "path_from", "path_to", "trace", "node"]) back.delete(name);
    if (back.toString()) params.set("back", back.toString());
    say(`Tracing ${from} → ${to}…`);
    window.location.href = `${form.dataset.pathUrl}?${params}`;
  };

  document.addEventListener("click", (event) => {
    if (event.target.closest("[data-path-pick]")) {
      event.preventDefault();
      if (picking) stop(); else start();
      return;
    }
    if (!picking) return;
    if (event.target.closest("[data-path-pick-cancel]")) { stop(); return; }
    const hit = pickable(event.target);
    if (!hit) return;
    // A pick, not the node's popover or its link to the ancestor viewer.
    event.preventDefault();
    event.stopPropagation();
    const { figure, up } = picking;
    if (!picking.from) {
      if (!hit.ancestor) { say("A path starts at an ancestor: click a node, not a sequence."); return; }
      picking.from = hit.name;
      hit.el.classList.add("path-start");
      say(`From ${hit.name}: now click a descendant (an ancestor or a sequence below it).`);
      return;
    }
    if (hit.name === picking.from) {  // the start again: choose another
      hit.el.classList.remove("path-start");
      picking.from = null;
      clearPreview(figure);
      say("Click the ancestor the path starts from.");
      return;
    }
    const way = between(up, picking.from, hit.name);
    if (!way) { say(`${hit.name} is not below ${picking.from}: a path runs down the tree. Pick a descendant, or click ${picking.from} again to start elsewhere.`); return; }
    const { from, to } = way;
    stop();
    trace(from, to);
  }, true);  // capture: before tree.js opens a popover

  document.addEventListener("mouseover", (event) => {
    if (!picking || !picking.from) return;
    const { figure, up } = picking;
    clearPreview(figure);
    const hit = pickable(event.target);
    if (!hit || hit.name === picking.from) return;
    const way = between(up, picking.from, hit.name);
    if (!way) return;
    hit.el.classList.add("path-target");
    // Each branch with its elbow from the parent, drawn over the tree as one route.
    const steps = new Set(way.names.slice(1));
    const points = [];
    let layer = null;
    figure.querySelectorAll("g.tree-branch[data-child]").forEach((g) => {
      if (!steps.has(g.dataset.child)) return;
      const edge = g.querySelector("line.edge");
      const [x1, x2, y] = ["x1", "x2", "y1"].map((a) => edge.getAttribute(a));
      points.push(`M ${x1} ${g.dataset.parentY} V ${y} H ${x2}`);
      layer = g.parentNode;
    });
    if (!layer) return;
    const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
    path.setAttribute("class", "path-preview-route");
    path.setAttribute("d", points.join(" "));
    layer.insertBefore(path, layer.querySelector(".tree-leaf, .tree-clade, a.tree-node"));  // under dots and names
  });

  document.addEventListener("keydown", (event) => { if (event.key === "Escape") stop(); });
  document.addEventListener("htmx:beforeSwap", stop);
})();
