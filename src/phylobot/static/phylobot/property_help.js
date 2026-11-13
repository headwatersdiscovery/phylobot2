// The (?) beside a property's title: show how it is computed, with its citations.
// Delegated from the document, so it keeps working after an HTMX swap replaces the
// ancestor viewer, and shared by the viewer and the comparison page. Without this
// script the panels are readable and the citations are on the About page.
(() => {
  document.addEventListener("click", (event) => {
    const button = event.target.closest && event.target.closest(".property-help");
    if (!button) return;
    const note = document.getElementById(button.getAttribute("aria-controls"));
    if (!note) return;
    const open = button.getAttribute("aria-expanded") === "true";
    button.setAttribute("aria-expanded", open ? "false" : "true");
    note.hidden = open;
  });
})();
