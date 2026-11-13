// Keep the reader's place on the properties-across-combinations page (owner,
// 2026-10-21): toggling a combination off and on is how the page is read, and a jump
// back to the top after every click made that comparison impossible to see.
//
// The page is one GET form, so every change is a real navigation — which is what keeps
// each view in the address and shareable. So the scroll position is carried across the
// navigation in sessionStorage instead.
//
// The controls submit through `this.form.submit()`, which does **not** fire a submit
// event, so there is nothing to listen for there. The position is saved on `change` and
// on `click` in the *capture* phase, which runs before the element's own onchange
// handler navigates away.
(() => {
  const KEY = "phylobot:compare-scroll";
  const form = () => document.getElementById("compare");

  function save() {
    if (!form()) return;
    try {
      sessionStorage.setItem(KEY, String(window.scrollY));
    } catch {
      /* private window, or site data blocked: the page still works, it just jumps. */
    }
  }

  for (const kind of ["change", "click"]) {
    document.addEventListener(kind, (event) => {
      const el = event.target;
      if (el && el.closest && el.closest("#compare")) save();
    }, true);
  }

  function restore() {
    if (!form()) return;
    let y = null;
    try {
      y = sessionStorage.getItem(KEY);
      sessionStorage.removeItem(KEY);  // one navigation only, so a later visit starts at the top
    } catch {
      return;
    }
    if (y === null) return;
    // The browser restores its own position on back/forward; this is for the forward
    // navigation a control causes, where it would otherwise land at the top.
    window.scrollTo(0, Number(y) || 0);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", restore);
  } else {
    restore();
  }
})();

// How many combinations per alignment to compare (owner, 2026-10-30): more than the
// default reads and draws many more distributions, so it asks first. Backing out puts
// the menu back. A new number shows every combination it brings in, so the old ticks
// (which named the previous set) are cleared before submitting.
(() => {
  function submitTop(select) {
    for (const box of select.form.querySelectorAll('input[name="show"]')) box.checked = false;
    select.form.submit();
  }

  function setUp() {
    const select = document.querySelector('#compare select[name="top"]');
    const dialog = document.getElementById("top-warning");
    if (!select || !dialog) return;
    let previous = select.value;
    const limit = Number(select.dataset.topDefault);

    select.addEventListener("change", () => {
      const n = select.value === "all" ? Number(select.dataset.perAlignment) : Number(select.value);
      if (n <= limit || typeof dialog.showModal !== "function") {
        submitTop(select);
        return;
      }
      const groups = Number(select.dataset.groups) || 1;
      dialog.querySelector("[data-top-count]").textContent = String(n * groups);
      dialog.showModal();
    });
    dialog.querySelector("[data-top-confirm]").addEventListener("click", () => {
      dialog.close();
      previous = select.value;
      submitTop(select);
    });
    const revert = () => {
      select.value = previous;
    };
    dialog.querySelector("[data-top-cancel]").addEventListener("click", () => {
      revert();
      dialog.close();
    });
    dialog.addEventListener("cancel", revert);  // Escape
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", setUp);
  } else {
    setUp();
  }
})();
