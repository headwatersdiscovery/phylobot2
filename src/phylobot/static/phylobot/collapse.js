// Collapsible page sections (the Settings page): each `section.collapsible[id]` gets its
// heading wrapped in a button that opens and closes it. Without this script every
// section stays open. Closed by default; open when it holds an error, when the URL names
// it (#id), or when this browser last left it open (localStorage, best effort).
(() => {
  const KEY = "phylobot_open_sections:" + location.pathname;
  const sections = [...document.querySelectorAll("section.collapsible[id]")];
  if (!sections.length) return;

  let saved = null;
  try {
    saved = JSON.parse(localStorage.getItem(KEY) || "null");
  } catch (e) {
    saved = null;
  }

  const setOpen = (section, open) => {
    section.classList.toggle("collapsed", !open);
    section.querySelector(":scope > h2 > .section-toggle").setAttribute("aria-expanded", String(open));
  };
  const remember = () => {
    const open = sections.filter((s) => !s.classList.contains("collapsed")).map((s) => s.id);
    try {
      localStorage.setItem(KEY, JSON.stringify(open));
    } catch (e) {
      // private window or blocked storage: the page still works, it just forgets
    }
  };
  const target = () => decodeURIComponent(location.hash.slice(1));

  for (const section of sections) {
    const heading = section.querySelector(":scope > h2");
    if (!heading) continue;
    const button = document.createElement("button");
    button.type = "button";
    button.className = "section-toggle";
    button.setAttribute("aria-controls", section.id);
    while (heading.firstChild) button.appendChild(heading.firstChild);
    heading.appendChild(button);
    button.addEventListener("click", () => {
      setOpen(section, section.classList.contains("collapsed"));
      remember();
    });
    const open =
      section.querySelector(".errors, .errorlist") !== null ||
      section.id === target() ||
      (Array.isArray(saved) && saved.includes(section.id));
    setOpen(section, open);
  }

  window.addEventListener("hashchange", () => {
    const section = sections.find((s) => s.id === target());
    if (section) {
      setOpen(section, true);
      remember();
    }
  });
})();
