// Help popovers (the "?" of _help.html): placed beside the button that opened them rather
// than in the middle of the window, and opened by resting the pointer on the "?" as well as
// by clicking it. A hover closes again when the pointer leaves; a click pins it until the
// reader clicks elsewhere or presses Escape (the popover API's own light dismiss). Delegated
// from the document, so it covers what HTMX swaps in. Without this script the "?" still
// opens its popover on a click, centred.
(() => {
  const OPEN_AFTER = 180;  // ms the pointer rests on a "?" before it opens
  const CLOSE_AFTER = 220;  // ms of grace to travel from the "?" into the popover
  const hovering = matchMedia("(hover: hover) and (pointer: fine)");
  const anchors = new WeakMap();  // popover -> the "?" that opened it (several may share one)
  const hovered = new WeakSet();  // popovers opened by hovering, not pinned by a click
  let opening = 0, closing = 0;

  const popoverOf = (dot) => document.getElementById(dot.getAttribute("popovertarget") || "");
  const isOpen = (pop) => pop.matches(":popover-open");

  const place = (pop) => {
    const dot = anchors.get(pop);
    if (!dot || !dot.isConnected) return;
    const r = dot.getBoundingClientRect(), pad = 8, gap = 6;
    Object.assign(pop.style, { position: "fixed", inset: "auto", margin: "0" });
    const w = pop.offsetWidth, h = pop.offsetHeight;
    const left = Math.min(Math.max(pad, r.left - 12), Math.max(pad, innerWidth - w - pad));
    const below = r.bottom + gap, above = r.top - gap - h;
    const top = below + h <= innerHeight - pad || above < pad ? below : above;
    pop.style.left = `${left}px`;
    pop.style.top = `${Math.max(pad, top)}px`;
  };

  // The toggle event comes after the popover is shown, so its size can be measured.
  document.addEventListener("toggle", (event) => {
    const pop = event.target;
    if (!(pop instanceof HTMLElement) || !pop.classList.contains("help-popover")) return;
    if (event.newState === "open") place(pop);
    else hovered.delete(pop);
  }, true);
  const replace = () => document.querySelectorAll(".help-popover:popover-open").forEach(place);
  addEventListener("scroll", replace, { capture: true, passive: true });
  addEventListener("resize", replace);

  document.addEventListener("click", (event) => {
    const dot = event.target.closest && event.target.closest(".help-dot[popovertarget]");
    const pop = dot && popoverOf(dot);
    if (!pop) return;
    clearTimeout(opening);
    if (hovered.has(pop) && isOpen(pop)) {
      // Already showing from the hover: the click keeps it (the default would close it).
      event.preventDefault();
      hovered.delete(pop);
      return;
    }
    anchors.set(pop, dot);
  }, true);

  document.addEventListener("mouseover", (event) => {
    if (!hovering.matches || !event.target.closest) return;
    const dot = event.target.closest(".help-dot[popovertarget]");
    if (dot) {
      const pop = popoverOf(dot);
      if (!pop || !pop.showPopover) return;
      clearTimeout(closing);
      clearTimeout(opening);
      opening = setTimeout(() => {
        if (isOpen(pop)) {
          if (!hovered.has(pop)) return;  // pinned by a click: leave it where it is
          anchors.set(pop, dot);
          place(pop);
          return;
        }
        anchors.set(pop, dot);
        hovered.add(pop);
        pop.showPopover();
      }, OPEN_AFTER);
    } else if (event.target.closest(".help-popover")) {
      clearTimeout(closing);
    }
  });

  document.addEventListener("mouseout", (event) => {
    if (!event.target.closest) return;
    const dot = event.target.closest(".help-dot[popovertarget]");
    const pop = dot ? popoverOf(dot) : event.target.closest(".help-popover");
    if (!pop) return;
    const to = event.relatedTarget;
    if (to && to.closest && (to.closest(".help-popover") === pop || (dot && to.closest(".help-dot") === dot))) return;
    clearTimeout(opening);
    if (!hovered.has(pop)) return;
    clearTimeout(closing);
    closing = setTimeout(() => { if (hovered.has(pop) && isOpen(pop)) pop.hidePopover(); }, CLOSE_AFTER);
  });
})();
