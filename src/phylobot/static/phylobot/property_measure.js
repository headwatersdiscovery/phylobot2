// The ancestor viewer's "measure this property" dialog (owner decision: the draws, the
// scope and the button were sitting in the panel whether or not anyone wanted them).
//
// Delegated from the document, so it survives the HTMX swaps that replace the viewer.
// Four jobs:
//   - swap the no-script checklist for the buttons it can drive, so the panel degrades
//     to a working form when this script is absent;
//   - open the dialog for the property that was clicked, carrying its key;
//   - stop the dialog's change events reaching the viewer's form, whose hx-trigger is
//     "change": a change there would re-request the viewer and close the dialog
//     mid-decision (the same trap as property_redraw.js);
//   - rescale the estimate as the draws and the scope change. The viewer cannot be asked
//     for it without closing the dialog, so the arithmetic is repeated here:
//     `properties_run.cost_of` is cost-for-one-ancestor-at-1000-draws x ancestors x
//     draws/1000, and `properties_run.cost_text` formats it. Keep the three in step.
(() => {
  function formatCost(seconds) {
    if (!(seconds > 0)) return "—";
    if (seconds < 1) return "under a second";
    const whole = Math.round(seconds);
    const hours = Math.floor(whole / 3600);
    const minutes = Math.floor((whole % 3600) / 60);
    const rest = whole % 60;
    if (hours) return minutes ? `${hours} h ${minutes} min` : `${hours} h`;
    if (minutes) return rest ? `${minutes} min ${rest} s` : `${minutes} min`;
    return `${rest} s`;
  }

  const dialog = () => document.getElementById("measure-dialog");

  function sync(box) {
    const unit = Number(box.dataset.unitCost || 0);
    const draws = Number(box.querySelector("[data-measure-draws]").value || 0);
    const scope = box.querySelector("[data-measure-scope]");
    const targets = Number(scope.selectedOptions[0]?.dataset.targets || 1);
    box.querySelector("[data-measure-targets]").textContent = String(targets);
    box.querySelector("[data-measure-cost]").textContent =
      formatCost(unit * targets * (draws / 1000));
    const estimate = box.querySelector("[data-measure-estimate]");
    if (estimate) estimate.hidden = !(unit > 0);  // nothing to price only if the row has no cost
  }

  // With the script present the buttons drive the dialog, so the plain checklist and its
  // controls (which only exist for readers without it) are taken out of the way.
  function dress() {
    for (const panel of document.querySelectorAll("[data-property-offer]")) {
      panel.querySelectorAll("[data-offer-fallback]").forEach((el) => { el.hidden = true; });
      panel.querySelectorAll(".offer-open").forEach((el) => { el.hidden = false; });
    }
  }
  document.addEventListener("DOMContentLoaded", dress);
  document.addEventListener("htmx:afterSwap", dress);
  dress();

  document.addEventListener("click", (event) => {
    const opener = event.target.closest && event.target.closest(".offer-open");
    if (opener) {
      const box = dialog();
      if (!box) return;
      const redraw = opener.dataset.offerRedraw === "on";
      box.dataset.unitCost = opener.dataset.offerUnitCost || "0";
      box.querySelector("[data-measure-key]").value = opener.dataset.offerKey || "";
      box.querySelector("[data-measure-label]").textContent =
        opener.dataset.offerLabel || "this property";
      // One dialog, two modes: a property already measured is drawn again with a new
      // seed (which `redraw` tells the server), and there is no cost to estimate for
      // it, since its own figures are already on the page.
      box.querySelector("[data-measure-redraw]").value = redraw ? "on" : "";
      box.querySelector("[data-measure-verb]").textContent = redraw ? "Recompute" : "Measure";
      box.querySelector("[data-measure-note]").hidden = redraw;
      box.querySelector("[data-redraw-note]").hidden = !redraw;
      sync(box);
      box.showModal();
      return;
    }
    const closer = event.target.closest && event.target.closest("[data-close-dialog]");
    if (closer && closer.closest("dialog")) closer.closest("dialog").close();
  });

  for (const kind of ["change", "input"]) {
    document.addEventListener(kind, (event) => {
      const box = event.target.closest && event.target.closest("dialog.measure-dialog");
      if (!box) return;
      event.stopPropagation();
      sync(box);
    }, true);
  }
})();
