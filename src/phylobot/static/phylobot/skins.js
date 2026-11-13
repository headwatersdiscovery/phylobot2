// The skin menu in the top bar, and the Display panel's menus on the Settings page
// (phylobot.skins). Every skin's stylesheet is already loaded (scoped by data-skin), so
// switching only sets the attribute; the cookies let the server render the next page the
// same way. backdrop.js watches the attributes.
(() => {
  const YEAR = "; path=/; max-age=31536000; SameSite=Lax";
  const root = document.documentElement;
  // What the reader picked for the skins that show a backdrop. Kept while they visit one
  // that shows none, so coming back restores it rather than losing the choice.
  let preference = root.dataset.backdropPreference || "none";

  function remember(name, value) {
    document.cookie = name + "=" + encodeURIComponent(value) + YEAR;
  }

  document.addEventListener("change", (event) => {
    const target = event.target;
    if (target.matches("[data-skin-picker]")) {
      root.dataset.skin = target.value;
      remember("phylobot_skin", target.value);
      // Only Clean, Dark and Spectrum show a backdrop; each option says which it is.
      const shows = target.selectedOptions[0].dataset.backdrop !== "none";
      root.dataset.backdrop = shows ? preference : "none";
      // The Display panel, if that is the page we are on: say why the box is empty, and
      // name the skin we have just moved to rather than the one the server rendered.
      document.querySelectorAll("[data-no-backdrop]").forEach((n) => (n.hidden = shows));
      document.querySelectorAll("[data-has-backdrop]").forEach((n) => (n.hidden = !shows));
      const name = document.querySelector("[data-skin-name]");
      if (name) name.textContent = target.selectedOptions[0].textContent;
    } else if (target.matches("[data-look]")) {
      // The backdrop's style, speed and density: each menu names the attribute it sets
      // and the cookie that remembers it, so adding a dial needs no change here.
      root.dataset[target.dataset.attr] = target.value;
      remember(target.dataset.cookie, target.value);
      if (target.dataset.cookie === "phylobot_backdrop") preference = target.value;
    }
  });
})();
