// PhyloBot '95 skin: the home page's helpful robot. "Don't show this again" hides it and
// sets a cookie, so the server stops drawing it (home.html checks the cookie).
(() => {
  document.addEventListener("click", (event) => {
    const button = event.target.closest("[data-retro-dismiss]");
    if (!button) return;
    document.cookie = "phylobot_helper=dismissed; path=/; max-age=31536000; SameSite=Lax";
    button.closest("[data-retro-helper]").hidden = true;
  });
})();
