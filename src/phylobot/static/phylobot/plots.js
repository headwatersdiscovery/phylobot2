// Figure downloads: a button with data-plot-download="png" or "svg" inside a <figure
// data-filename="..."> saves the figure's own SVG, drawn in the browser from the page's
// markup (plot styles are SVG attributes, so the file looks as shown, without the
// page's CSS and so without the skin). Elements classed "hover" or "guide" are left
// out. Delegated, so it keeps working after HTMX swaps.
(() => {
  function save(blob, filename) {
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    document.body.append(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(link.href), 1000);
  }
  function standalone(svg) {
    const copy = svg.cloneNode(true);
    copy.querySelectorAll(".hover, .guide").forEach((element) => element.remove());
    copy.removeAttribute("class");
    for (const name of [...copy.getAttributeNames()]) {
      if (name.startsWith("data-")) copy.removeAttribute(name);
    }
    return new XMLSerializer().serializeToString(copy);
  }
  document.addEventListener("click", (event) => {
    const button = event.target.closest && event.target.closest("[data-plot-download]");
    if (!button) return;
    const figure = button.closest("figure");
    const svg = figure.querySelector(":scope > svg");
    const text = standalone(svg);
    const name = figure.dataset.filename;
    if (button.dataset.plotDownload === "svg") {
      save(new Blob([text], { type: "image/svg+xml" }), `${name}.svg`);
      return;
    }
    const scale = 3;  // print quality
    const box = svg.viewBox.baseVal;
    const image = new Image();
    image.onload = () => {
      const canvas = document.createElement("canvas");
      canvas.width = box.width * scale;
      canvas.height = box.height * scale;
      canvas.getContext("2d").drawImage(image, 0, 0, canvas.width, canvas.height);
      canvas.toBlob((blob) => save(blob, `${name}.png`), "image/png");
    };
    image.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(text)}`;
  });
})();

// Figure tooltips: in a <figure data-tips>, hovering an element with data-tip shows its
// text in the figure's .pp-tip at the pointer (instant, unlike an SVG <title>). Only
// figures that opt in, since the ancestor viewer's trace drives its own .pp-tip.
(() => {
  function tipOf(event) {
    const figure = event.target.closest && event.target.closest("figure[data-tips]");
    return figure && [figure, figure.querySelector(":scope > .pp-tip")];
  }
  document.addEventListener("pointermove", (event) => {
    const found = tipOf(event);
    if (!found || !found[1]) return;
    const [figure, tip] = found;
    const source = event.target.closest("[data-tip]");
    if (!source) {
      tip.hidden = true;
      return;
    }
    tip.textContent = source.dataset.tip;
    tip.hidden = false;
    const box = figure.getBoundingClientRect();
    const px = event.clientX - box.left;
    tip.style.left = `${Math.max(0, Math.min(px + 12, box.width - tip.offsetWidth - 4))}px`;
    tip.style.top = `${event.clientY - box.top + 14}px`;
  });
  document.addEventListener("pointerout", (event) => {
    const found = tipOf(event);
    if (found && found[1] && !found[0].contains(event.relatedTarget)) found[1].hidden = true;
  });
})();

// Flag tooltips (.length-flag > .flag-tip, shown by CSS on hover or focus): fixed, so a
// table's scroll box cannot clip them, and placed here under the flag, kept on screen.
(() => {
  function place(event) {
    const flag = event.target.closest && event.target.closest(".length-flag");
    const tip = flag && flag.querySelector(".flag-tip");
    if (!tip) return;
    const box = flag.getBoundingClientRect();
    const width = tip.offsetWidth || 352;
    tip.style.left = `${Math.max(8, Math.min(box.left, innerWidth - width - 8))}px`;
    tip.style.top = `${box.bottom + 4}px`;
  }
  document.addEventListener("mouseover", place);
  document.addEventListener("focusin", place);
})();
