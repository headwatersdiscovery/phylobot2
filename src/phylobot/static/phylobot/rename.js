// Renaming in place (a job's or a search's title: _job_title.html, _homolog_title.html):
// the pencil shows the form; Cancel or Escape hides it again, with the old name back.
// Saving is an HTMX post that swaps the whole title back, done editing.
(function () {
  function editing(title, on) {
    title.querySelector(".title-view").hidden = on;
    const form = title.querySelector(".title-edit");
    form.hidden = !on;
    const input = form.querySelector("input[name=name]");
    if (on) {
      input.focus();
      input.select();
    } else {
      input.value = input.defaultValue;
    }
    changed(title);
  }
  // A Save beside the name (button.confirm: the library's and alignment's titles) turns
  // amber once the name differs from the saved one, so it is clear it still needs a click.
  function changed(title) {
    const confirm = title.querySelector(".title-edit .confirm");
    const input = title.querySelector(".title-edit input[name=name]");
    if (confirm && input) confirm.classList.toggle("dirty", input.value !== input.defaultValue);
  }
  document.addEventListener("input", (event) => {
    const title = event.target.closest && event.target.closest(".editable-title");
    if (title) changed(title);
  });
  document.addEventListener("click", (event) => {
    const title = event.target.closest(".editable-title");
    if (!title) return;
    if (event.target.closest("[data-rename]")) editing(title, true);
    if (event.target.closest("[data-rename-cancel]")) editing(title, false);
  });
  document.addEventListener("keydown", (event) => {
    const title = event.target.closest && event.target.closest(".editable-title");
    if (title && event.key === "Escape") editing(title, false);
  });
})();
