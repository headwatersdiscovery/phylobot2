// The Simulations and ESR pages (form[data-show-when-form]): show each option's settings
// only when it is chosen.
// An element with data-show-when="field=value [field=value ...]" is shown when any of its
// conditions holds; a checkbox's value is "on" when ticked. Re-applied after HTMX swaps
// (the tree and ancestor menus).
(() => {
  const form = document.querySelector("form[data-show-when-form]");
  if (!form) return;

  function current(name) {
    const fields = form.querySelectorAll(`[name="${name}"]`);
    for (const field of fields) {
      if (field.type === "radio") { if (field.checked) return field.value; }
      else if (field.type === "checkbox") return field.checked ? "on" : "";
      else return field.value;
    }
    return "";
  }

  function apply() {
    for (const el of form.querySelectorAll("[data-show-when]")) {
      const show = el.dataset.showWhen.split(/\s+/).some((rule) => {
        const [name, value] = rule.split("=");
        return current(name) === value;
      });
      el.hidden = !show;
    }
    const shape = current("shape");
    for (const help of form.querySelectorAll("#shape-help [data-shape]")) {
      help.hidden = help.dataset.shape !== shape;
    }
  }

  form.addEventListener("change", apply);
  document.body.addEventListener("htmx:afterSwap", apply);
  apply();
})();
