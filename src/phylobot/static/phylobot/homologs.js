// Find homologs, review panel (_homolog_filters.html). The panel is re-rendered by HTMX
// on every change, so all listeners live on the document.
(function () {
  // Sliders paired with number boxes. Dragging updates the box without a request; the
  // request goes when the slider is released ("change") or the typing pauses.
  document.addEventListener("input", (event) => {
    const el = event.target;
    if (el.matches("input[type=range][data-for]")) {
      const box = document.getElementById(el.dataset.for);
      if (box) box.value = el.value;
    } else if (el.matches(".homolog-filters input[type=number]")) {
      const slider = document.querySelector(`input[type=range][data-for="${el.id}"]`);
      if (slider && el.value !== "") slider.value = el.value;
    }
  });

  // A tick that disagrees with the filters becomes an include (ticked, filters remove
  // it) or exclude (unticked, filters keep it). Capture phase: this runs before HTMX
  // sends the form on the same "change".
  document.addEventListener("change", (event) => {
    const box = event.target;
    if (!box.matches("input.pick")) return;
    const holder = box.form.querySelector(".choices");
    holder.querySelectorAll("input[name=include], input[name=exclude]").forEach((h) => {
      if (h.value === box.dataset.id) h.remove();
    });
    if (box.checked !== (box.dataset.default === "1")) {
      const h = document.createElement("input");
      h.type = "hidden";
      h.name = box.checked ? "include" : "exclude";
      h.value = box.dataset.id;
      holder.appendChild(h);
    }
  }, true);

  // The hit table's controls change the form's hidden state, then send it once:
  // a heading's button sets ?sort= (again on the same column reverses it); "Select: all /
  // none" sets the base the ticks apply to and drops the ticks; "Reset to the filters"
  // drops both (its href is the same view without them, for opening in a new tab).
  function send(form) {
    form.querySelector(".choices").dispatchEvent(new Event("change", { bubbles: true }));
  }
  function setHidden(form, name, value) {
    const holder = form.querySelector(".choices");
    holder.querySelectorAll(`input[name=${name}]`).forEach((h) => h.remove());
    if (value === null) return;
    const h = document.createElement("input");
    h.type = "hidden";
    h.name = name;
    h.value = value;
    holder.appendChild(h);
  }
  document.addEventListener("click", (event) => {
    const sort = event.target.closest("button.sort[data-sort]");
    const base = event.target.closest("[data-select-base]");
    const reset = event.target.closest("a[data-reset]");
    if (!sort && !base && !reset) return;
    const form = (sort || base || reset).closest("form");
    if (sort) {
      const key = sort.dataset.sort;
      const current = form.querySelector(".choices input[name=sort]:not([disabled])");
      const now = current ? current.value : "rank";
      setHidden(form, "sort", now === key ? `-${key}` : key);
    } else {
      event.preventDefault();
      setHidden(form, "include", null);
      setHidden(form, "exclude", null);
      setHidden(form, "base", base ? base.dataset.selectBase : null);
    }
    send(form);
  });

  // Remember which hit tables are open, so the next render keeps them open.
  document.addEventListener("toggle", (event) => {
    const details = event.target;
    if (!(details instanceof HTMLDetailsElement) || !details.dataset.table) return;
    const holder = details.closest("form").querySelector(".choices");
    holder.querySelectorAll(`input[name=show][value="${details.dataset.table}"]`).forEach((h) => h.remove());
    if (details.open) {
      const h = document.createElement("input");
      h.type = "hidden";
      h.name = "show";
      h.value = details.dataset.table;
      holder.appendChild(h);
    }
  }, true);
})();
