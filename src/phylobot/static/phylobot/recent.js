// "Jump Back In" on the home page: the result pages this browser last had open.
//
// Every viewer is a GET form whose address holds its whole state, so remembering a page is
// remembering its address. That is kept here, in localStorage, per signed-in user: the
// server writes nothing for a page view. A page worth returning to names its kind on
// <body data-page-kind> (phylobot.skins.PAGES); one entry is kept per page, and its
// address follows the viewer's own controls (HTMX pushes each change into the history).
// The home page draws the list, after asking the server which records still exist
// (home_views.alive), so a deleted job's pages drop off it.
(() => {
  const KEPT = 12, SHOWN = 8;
  const UUID = /[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}/g;
  const body = document.body;
  const key = "phylobot_recent:" + (body.dataset.user || "");

  function read() {
    try {
      const list = JSON.parse(localStorage.getItem(key) || "[]");
      // Only what this script wrote: an address on this site, and text.
      return Array.isArray(list) ? list.filter((e) => e && typeof e.url === "string"
        && /^\/(?!\/)/.test(e.url) && typeof e.path === "string" && typeof e.title === "string") : [];
    } catch { return []; }
  }
  function write(list) {
    try { localStorage.setItem(key, JSON.stringify(list.slice(0, KEPT))); } catch { /* private window */ }
  }

  // ---- remembering ---------------------------------------------------------------

  // "Tree · Ime2 · PhyloBot" -> "Ime2": the page's own name, without the words the list's
  // kind label already says.
  function title(kind) {
    const dull = new Set(["PhyloBot", "PhyloBot 2", kind, "Find Homologs"]);
    const parts = document.title.split(" · ").map((part) => part.trim()).filter((part) => part && !dull.has(part));
    return parts.join(" · ").slice(0, 120);
  }

  function remember() {
    const kind = body.dataset.pageKind;
    if (!kind || !body.dataset.user) return;
    const path = location.pathname;
    const entry = { path, url: path + location.search, kind, title: title(kind), at: Date.now() };
    write([entry, ...read().filter((e) => e.path !== path)]);
  }
  remember();
  document.addEventListener("htmx:pushedIntoHistory", remember);
  document.addEventListener("htmx:replacedInHistory", remember);

  // ---- the home page's list -----------------------------------------------------

  const panel = document.querySelector("[data-recent-pages]");
  if (!panel) return;
  const list = panel.querySelector("ol");

  function ago(at) {
    const minutes = Math.floor((Date.now() - at) / 60000);
    if (minutes < 1) return "just now";
    if (minutes < 60) return minutes + " min ago";
    if (minutes < 24 * 60) return Math.floor(minutes / 60) + " h ago";
    const days = Math.floor(minutes / (24 * 60));
    if (days === 1) return "yesterday";
    if (days < 14) return days + " days ago";
    return new Date(at).toISOString().slice(0, 10);
  }

  function draw(entries) {
    list.replaceChildren(...entries.slice(0, SHOWN).map((entry) => {
      // Built from text nodes only: titles are the user's own words.
      const item = document.createElement("li");
      const kind = document.createElement("span");
      kind.className = "kind";
      kind.textContent = entry.kind || "Page";
      const link = document.createElement("a");
      link.href = entry.url;
      link.textContent = entry.title || entry.kind || entry.path;
      const when = document.createElement("span");
      when.className = "when";
      when.textContent = Number.isFinite(entry.at) ? ago(entry.at) : "";
      const forget = document.createElement("button");
      forget.type = "button";
      forget.className = "link forget";   // .link: no skin boxes it as a button
      forget.title = "Remove from this list";
      forget.setAttribute("aria-label", "Remove " + link.textContent + " from this list");
      forget.textContent = "×";
      forget.addEventListener("click", () => {
        const kept = read().filter((e) => e.path !== entry.path);
        write(kept);
        draw(kept);
      });
      item.append(kind, link, when, forget);
      return item;
    }));
    panel.hidden = !entries.length;
  }

  panel.querySelector("[data-recent-clear]").addEventListener("click", () => { write([]); draw([]); });

  // Ask which of the records in these addresses still exist; until the answer comes (or if
  // it never does) the list shows what is remembered.
  const entries = read();
  draw(entries);
  const ids = [...new Set(entries.flatMap((entry) => entry.path.match(UUID) || []))];
  if (!ids.length) return;
  const query = ids.map((id) => "id=" + id).join("&");
  fetch(panel.dataset.aliveUrl + "?" + query, { headers: { Accept: "application/json" } })
    .then((response) => (response.ok ? response.json() : null))
    .then((answer) => {
      if (!answer || !Array.isArray(answer.live)) return;
      const live = new Set(answer.live);
      const kept = read().filter((entry) => (entry.path.match(UUID) || []).every((id) => live.has(id)));
      write(kept);
      draw(kept);
    })
    .catch(() => { /* offline: keep what is shown */ });
})();
