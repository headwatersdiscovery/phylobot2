// The tree viewer's overview: the whole tree as a thumbnail pinned at the figure's upper
// right while the page scrolls (as the alignment viewer's overview strip, turned upright),
// with a box over the part of the tree in the window. Drag the box, or click anywhere on
// the thumbnail, to move. It is drawn on a canvas from the tree's own SVG lines, in the
// colours they have on the page (so a traced path and coloured branches show, and the dark
// skins' remapped ink is used), and appears only when the tree is taller than the window.
// Nothing is added to the figure itself, so its PNG and SVG downloads are unchanged.
(() => {
  const WIDTH = 150, MAX_HEIGHT = 420;  // CSS px; the thumbnail keeps the tree's proportions
  let state = null;  // { svg, box, canvas, frame, scale }

  const barHeight = () => {
    const bar = document.querySelector("header");
    return bar && getComputedStyle(bar).position === "sticky" ? bar.offsetHeight : 0;
  };

  function draw() {
    const { svg, canvas } = state;
    const page = svg.getBoundingClientRect();
    if (!page.width || !page.height) return;
    const limit = Math.min(MAX_HEIGHT, (innerHeight - barHeight()) * 0.6);
    const scale = Math.min(WIDTH / page.width, limit / page.height);
    const w = Math.round(page.width * scale), h = Math.round(page.height * scale);
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    canvas.style.width = `${w}px`;
    canvas.style.height = `${h}px`;
    state.scale = scale;
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr * scale, 0, 0, dpr * scale, 0, 0);
    ctx.clearRect(0, 0, page.width, page.height);
    const hair = 1 / scale;  // one thumbnail pixel, in page units
    // Folded clades: their triangles, apex at the left.
    svg.querySelectorAll(".tree-clade polygon").forEach((shape) => {
      const r = shape.getBoundingClientRect();
      const style = getComputedStyle(shape);
      ctx.beginPath();
      ctx.moveTo(r.left - page.left, r.top - page.top + r.height / 2);
      ctx.lineTo(r.right - page.left, r.top - page.top);
      ctx.lineTo(r.right - page.left, r.bottom - page.top);
      ctx.closePath();
      ctx.fillStyle = style.fill === "none" ? "transparent" : style.fill;
      ctx.fill();
      ctx.strokeStyle = style.stroke;
      ctx.lineWidth = hair;
      ctx.stroke();
    });
    // Branches and connectors are level or upright lines; the slanted "cut" marks on long
    // branches and the wide transparent lines that catch the pointer are left out.
    svg.querySelectorAll("line, polyline").forEach((line) => {
      const stroke = line.getAttribute("stroke");
      if (!stroke || stroke === "transparent" || line.closest("[font-size]")) return;
      const r = line.getBoundingClientRect();
      if (line.tagName === "line" && r.width > 1.5 && r.height > 1.5) return;
      const path = stroke === "#e08214" || stroke === "#8e44ad";  // a traced path, a foreground
      ctx.fillStyle = getComputedStyle(line).stroke;
      const thick = hair * (path ? 2.5 : 1);
      if (line.tagName === "polyline") {  // a path's elbow: its upright, then its level part
        ctx.fillRect(r.left - page.left - thick / 2, r.top - page.top, thick, r.height);
        ctx.fillRect(r.left - page.left, r.bottom - page.top - thick / 2, r.width, thick);
      } else if (r.width >= r.height) {
        ctx.fillRect(r.left - page.left, r.top - page.top + r.height / 2 - thick / 2, Math.max(r.width, hair), thick);
      } else {
        ctx.fillRect(r.left - page.left + r.width / 2 - thick / 2, r.top - page.top, thick, Math.max(r.height, hair));
      }
    });
  }

  function place() {
    if (!state || !state.svg.isConnected) return;
    const { svg, box, frame } = state;
    const page = svg.getBoundingClientRect();
    const top = barHeight();
    const needed = page.height > (innerHeight - top) * 1.15;
    if (box.hidden === needed) {
      box.hidden = !needed;
      if (needed) draw();
    }
    if (!needed) return;
    box.parentElement.style.top = `${top + 8}px`;
    const clamp = (v) => Math.max(0, Math.min(1, v));
    const y0 = clamp((top - page.top) / page.height), y1 = clamp((innerHeight - page.top) / page.height);
    const x0 = clamp(-page.left / page.width), x1 = clamp((innerWidth - page.left) / page.width);
    const w = page.width * state.scale, h = page.height * state.scale;
    frame.style.left = `${x0 * w}px`;
    frame.style.width = `${(x1 - x0) * w}px`;
    frame.style.top = `${y0 * h}px`;
    frame.style.height = `${Math.max(0, y1 - y0) * h}px`;
  }

  function init() {
    const box = document.querySelector(".tree-overview");
    const svg = document.querySelector(".tree-figure > svg");
    if (!box || !svg) { state = null; return; }
    if (state && state.box === box) return;
    state = { svg, box, canvas: box.querySelector("canvas"), frame: box.querySelector(".tree-overview-window"), scale: 1 };
    box.hidden = true;
    // Dragging the box keeps hold of it where it was grabbed; a click elsewhere centres
    // the window there.
    let grab = null;
    const go = (event) => {
      const thumb = box.getBoundingClientRect();
      const page = svg.getBoundingClientRect();
      const fy = (event.clientY - thumb.top) / thumb.height - grab;
      const view = innerHeight - barHeight();
      window.scrollTo({ top: window.scrollY + page.top + fy * page.height - barHeight() - view / 2 });
    };
    box.addEventListener("pointerdown", (event) => {
      if (event.button !== 0) return;
      event.preventDefault();
      const thumb = box.getBoundingClientRect();
      const frame = state.frame.getBoundingClientRect();
      const inside = event.clientY >= frame.top && event.clientY <= frame.bottom;
      grab = inside ? (event.clientY - (frame.top + frame.height / 2)) / thumb.height : 0;
      box.setPointerCapture(event.pointerId);
      go(event);
    });
    box.addEventListener("pointermove", (event) => { if (grab !== null) go(event); });
    const release = () => { grab = null; };
    box.addEventListener("pointerup", release);
    box.addEventListener("pointercancel", release);
    place();
  }

  let queued = false;
  const later = () => {
    if (queued) return;
    queued = true;
    requestAnimationFrame(() => { queued = false; place(); });
  };
  addEventListener("scroll", later, { passive: true });
  addEventListener("resize", () => { if (state && !state.box.hidden) draw(); later(); });
  document.addEventListener("DOMContentLoaded", init);
  document.addEventListener("htmx:afterSettle", init);
  // Folding a section above the tree moves it without a scroll.
  document.addEventListener("click", () => setTimeout(later, 60));
  if (document.readyState !== "loading") init();
})();
