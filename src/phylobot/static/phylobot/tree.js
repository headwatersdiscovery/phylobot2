// The tree viewer's node popovers (_tree_viewer.html, `analog_views.popovers`): clicking
// an ancestral node opens its popover (owner, 2026-10-30: hovering made it hard to use),
// with the way to the ancestor viewer, the ancestor's name to set or change, and the same
// ancestor in the job's other combinations. It stays until another node is clicked, a
// click lands outside it, its × or Escape. Without this script the node's link opens
// the ancestor viewer, as before. Listeners sit on the document, so they keep working
// after HTMX swaps the viewer.
(() => {
  let shown = null;

  // The popover for a clicked thing (ancestral node or collapsed-clade triangle), matched
  // by data-node within the figure; `.clade-popover` ones are the triangles'.
  const popoverFor = (el, clade) => {
    const figure = el.closest("figure");
    if (!figure) return null;
    const selector = clade ? ".clade-popover" : ".node-popover:not(.clade-popover)";
    return [...figure.querySelectorAll(selector)]
      .find((p) => p.dataset.node === el.dataset.node) || null;
  };

  const hide = () => {
    if (shown) shown.hidden = true;
    shown = null;
  };

  const show = (anchor, popover) => {
    if (shown && shown !== popover) shown.hidden = true;
    const figure = popover.closest("figure");
    const dot = anchor.querySelector(".dot") || anchor;
    const box = figure.getBoundingClientRect();
    const at = dot.getBoundingClientRect();
    popover.hidden = false;
    // Beside the node, on the right; flipped left if it would leave the figure.
    let left = at.right - box.left + figure.scrollLeft + 10;
    if (left + popover.offsetWidth > figure.scrollLeft + figure.clientWidth) {
      left = Math.max(0, at.left - box.left + figure.scrollLeft - popover.offsetWidth - 10);
    }
    popover.style.left = `${left}px`;
    popover.style.top = `${at.top - box.top + figure.scrollTop - 12}px`;
    shown = popover;
    const first = popover.querySelector("a.button, button, input");
    if (first) first.focus({ preventScroll: true });
  };

  const back = () => window.location.pathname + window.location.search;

  // Name or rename through the page's one POST form (tree.html): the popover sits inside
  // the viewer's GET form, where a form of its own cannot.
  const postName = (popover, remove) => {
    const form = document.getElementById("tree-ancestor-name");
    const input = popover.querySelector("[data-pop-name]");
    if (!form || !input) return;
    if (!remove && !input.reportValidity()) return;
    const set = (name, value) => { form.elements[name].value = value; };
    set("node", popover.dataset.node);
    set("combination", popover.dataset.combination);
    set("aligner", popover.dataset.aligner);
    set("ancestor_name", input.value.trim());
    set("remove", remove ? "1" : "");
    set("return_to", back());
    form.submit();
  };

  // Collapse / expand / restyle a clade, through the page's #tree-clade-collapse form.
  // `style` (the restyle values) is read from the triangle's popover; a bare collapse
  // sends none, so the server applies the defaults.
  const postCollapse = (popover, action) => {
    const form = document.getElementById("tree-clade-collapse");
    if (!form) return;
    const set = (name, value) => { form.elements[name].value = value; };
    set("node", popover.dataset.node);
    set("combination", popover.dataset.combination);
    set("aligner", popover.dataset.aligner);
    set("action", action);
    if (action === "update") {
      const label = popover.querySelector("[data-clade-label]");
      const fill = popover.querySelector("[data-clade-fill-input]");
      const border = popover.querySelector("[data-clade-border-input]");
      const width = popover.querySelector("[data-clade-width]");
      set("label", label ? label.value.trim() : "");
      set("fill", popover.dataset.fill || (fill ? fill.value : ""));
      set("border_colour", popover.dataset.border || (border ? border.value : ""));
      set("border_width", width ? width.value : "");
    } else {
      ["label", "fill", "border_colour", "border_width"].forEach((n) => set(n, ""));
    }
    set("return_to", back());
    form.submit();
  };

  const nodeOf = (target) => target.closest && target.closest("a.tree-node[data-node]");
  const cladeOf = (target) => target.closest && target.closest(".tree-clade[data-node]");

  document.addEventListener("click", (event) => {
    const link = nodeOf(event.target);
    if (link) {
      const popover = popoverFor(link, false);
      if (!popover) return;  // no popover: the link opens the ancestor viewer
      event.preventDefault();
      if (shown === popover) hide();
      else show(link, popover);
      return;
    }
    const triangle = cladeOf(event.target);
    if (triangle) {
      const popover = popoverFor(triangle, true);
      if (!popover) return;
      event.preventDefault();
      if (shown === popover) hide();
      else show(triangle, popover);
      return;
    }
    const popover = event.target.closest && event.target.closest(".node-popover");
    if (!popover) { hide(); return; }
    if (event.target.closest("[data-pop-close]")) { hide(); return; }
    if (event.target.closest("[data-pop-save]")) { postName(popover, false); return; }
    if (event.target.closest("[data-pop-remove]")) { postName(popover, true); return; }
    if (event.target.closest("[data-pop-collapse]")) { postCollapse(popover, "collapse"); return; }
    if (event.target.closest("[data-pop-expand]")) { postCollapse(popover, "expand"); return; }
    if (event.target.closest("[data-pop-restyle]")) { postCollapse(popover, "update"); return; }
    const fill = event.target.closest("[data-clade-fill]");
    if (fill) {
      popover.dataset.fill = fill.dataset.cladeFill;
      const input = popover.querySelector("[data-clade-fill-input]");
      if (input && fill.dataset.cladeFill !== "none") input.value = fill.dataset.cladeFill;
      return;
    }
    const border = event.target.closest("[data-clade-border]");
    if (border) {
      popover.dataset.border = border.dataset.cladeBorder;
      const input = popover.querySelector("[data-clade-border-input]");
      if (input) input.value = border.dataset.cladeBorder;
    }
  });
  // A custom colour from the picker is authoritative over the last swatch.
  document.addEventListener("input", (event) => {
    const t = event.target;
    if (!t.matches) return;
    if (t.matches("[data-clade-fill-input]")) t.closest(".node-popover").dataset.fill = t.value;
    else if (t.matches("[data-clade-border-input]")) t.closest(".node-popover").dataset.border = t.value;
  });
  // Every control in a popover sits inside the viewer's GET form: its change must not reach
  // that form (HTMX would redraw the tree and lose the popover).
  document.addEventListener("change", (event) => {
    if (event.target.closest && event.target.closest(".node-popover")) event.stopPropagation();
  }, true);
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape") { hide(); return; }
    if (event.key !== "Enter" || !event.target.matches) return;
    if (event.target.matches("[data-pop-name]")) {
      event.preventDefault();
      postName(event.target.closest(".node-popover"), false);
    } else if (event.target.matches("[data-clade-label]")) {
      event.preventDefault();
      postCollapse(event.target.closest(".node-popover"), "update");
    }
  });
  document.addEventListener("htmx:beforeSwap", hide);
})();

// The tree viewer's Caption: open or closed as the reader left it (this browser), since
// the viewer is redrawn on every change.
(() => {
  const KEY = "phylobot-tree-caption";
  const apply = () => {
    const caption = document.querySelector("details.tree-caption");
    if (!caption) return;
    try { caption.open = localStorage.getItem(KEY) === "open"; } catch (_) { /* private mode */ }
  };
  document.addEventListener("toggle", (event) => {
    if (!event.target.matches || !event.target.matches("details.tree-caption")) return;
    try { localStorage.setItem(KEY, event.target.open ? "open" : "closed"); } catch (_) { /* ignore */ }
  }, true);
  document.addEventListener("DOMContentLoaded", apply);
  document.addEventListener("htmx:afterSettle", apply);
})();
