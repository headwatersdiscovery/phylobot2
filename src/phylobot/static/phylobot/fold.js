// Sections that fold away (the job page's Results and each alignment in it): a heading's
// button.fold-toggle shows or hides the element its aria-controls names. Delegated from
// the document, so it keeps working after HTMX replaces the job's status block; what the
// reader folded is remembered for the page and put back after each refresh. Following a
// link to something folded away unfolds it. Without this script everything stays open.
(() => {
  const folded = new Map();  // body id -> hidden?

  const set = (button, open) => {
    const body = document.getElementById(button.getAttribute("aria-controls"));
    if (!body) return;
    button.setAttribute("aria-expanded", String(open));
    body.hidden = !open;
    folded.set(body.id, !open);
  };

  const toggleOf = (body) => document.querySelector(`button.fold-toggle[aria-controls="${CSS.escape(body.id)}"]`);

  // Unfold everything that hides `target`, the target too if it is a folding section.
  const reveal = (target) => {
    for (let node = target; node; node = node.parentElement) {
      if (node.classList && node.classList.contains("fold-body") && node.hidden) {
        const button = toggleOf(node);
        if (button) set(button, true);
      }
    }
    const own = target.querySelector(":scope > h2 > button.fold-toggle, :scope > h3 > button.fold-toggle");
    if (own && own.getAttribute("aria-expanded") === "false") set(own, true);
  };

  document.addEventListener("click", (event) => {
    const button = event.target.closest && event.target.closest("button.fold-toggle");
    if (button) {
      set(button, button.getAttribute("aria-expanded") !== "true");
      return;
    }
    const link = event.target.closest && event.target.closest("a[href^='#']");
    if (link && link.getAttribute("href").length > 1) {
      const target = document.getElementById(link.getAttribute("href").slice(1));
      if (target) reveal(target);  // before the browser scrolls to it
    }
  });

  const restore = () => {
    for (const [id, hidden] of folded) {
      const body = document.getElementById(id);
      const button = body && toggleOf(body);
      if (button) {
        button.setAttribute("aria-expanded", String(!hidden));
        body.hidden = hidden;
      }
    }
  };
  document.addEventListener("htmx:afterSwap", restore);
  const fromHash = () => {
    const target = location.hash.length > 1 && document.getElementById(decodeURIComponent(location.hash.slice(1)));
    if (target) { reveal(target); target.scrollIntoView(); }
  };
  window.addEventListener("hashchange", fromHash);
})();
