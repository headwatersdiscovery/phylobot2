// "On this page" sidebars: mark the link of the section being read (aria-current), so the
// contents show where you are. Without this script the sidebar still works as links.
//
// Any <nav data-toc> whose links point at ids on the page is handled, so the About page,
// the ancestor viewer and the job page share one implementation. The viewer is swapped in
// place by HTMX, and the job page's list while the job runs, so the links are looked up
// afresh on every update rather than bound once. A section folded away (fold.js) has no
// layout, so it is never "current".
(() => {
  const update = () => {
    for (const nav of document.querySelectorAll("nav[data-toc]")) {
      const links = new Map();
      for (const a of nav.querySelectorAll("a[href^='#']")) {
        const target = document.getElementById(a.getAttribute("href").slice(1));
        if (target && target.getClientRects().length) links.set(target, a);
      }
      if (!links.size) continue;
      const targets = [...links.keys()];
      // The last heading or section that has passed just below the sticky top bar.
      const line = 90;
      let current = targets[0];
      for (const t of targets) if (t.getBoundingClientRect().top <= line) current = t;
      if (window.innerHeight + window.scrollY >= document.body.scrollHeight - 2) {
        current = targets.at(-1);
      }
      for (const a of nav.querySelectorAll("a[aria-current]")) a.removeAttribute("aria-current");
      links.get(current).setAttribute("aria-current", "true");
    }
  };
  let queued = false;
  const later = () => {
    if (queued) return;
    queued = true;
    requestAnimationFrame(() => { queued = false; update(); });
  };
  window.addEventListener("scroll", later, { passive: true });
  document.addEventListener("DOMContentLoaded", update);
  document.addEventListener("htmx:afterSettle", update);
  document.addEventListener("click", () => setTimeout(later));  // a section folded or opened
})();
