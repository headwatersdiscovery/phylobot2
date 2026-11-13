// The top bar's menus (Sequences, Validation, Account): each is a <details>, so it opens
// and closes on its own summary without this script. This closes an open one when the
// reader clicks elsewhere or presses Escape, and keeps one open at a time where the
// browser does not yet do that itself from the shared `name`.
(() => {
  const menus = () => document.querySelectorAll("header nav details.nav-menu");
  const close = (except) => menus().forEach((menu) => { if (menu !== except) menu.open = false; });
  document.addEventListener("click", (event) => {
    const inside = event.target.closest && event.target.closest("header nav details.nav-menu");
    close(inside);
  });
  document.addEventListener("keydown", (event) => {
    if (event.key !== "Escape") return;
    const open = document.querySelector("header nav details.nav-menu[open]");
    if (!open) return;
    open.open = false;
    open.querySelector("summary").focus();
  });
  document.addEventListener("toggle", (event) => {
    const menu = event.target;
    if (menu.matches && menu.matches("header nav details.nav-menu") && menu.open) close(menu);
  }, true);
})();
