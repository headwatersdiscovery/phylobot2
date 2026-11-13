// My usage: hover and keyboard read-outs for the server-drawn charts (phylobot.usage).
// Line charts: a column per point; the crosshair and a dot snap to it. Bars and
// segments: the mark itself. The tooltip only repeats what the tables beneath show.
(() => {
  const tip = document.getElementById("viz-tip");
  if (!tip) return;

  function show(mark, svgX, svgY) {
    tip.replaceChildren();
    const value = document.createElement("strong");
    value.textContent = mark.dataset.tipValue;  // text, never HTML: names are user data
    const label = document.createElement("span");
    label.textContent = mark.dataset.tipLabel;
    tip.append(value, label);
    tip.hidden = false;
    const svg = mark.ownerSVGElement;
    const box = svg.getBoundingClientRect();
    const scale = box.width / svg.viewBox.baseVal.width;
    const x = box.left + window.scrollX + svgX * scale;
    const y = box.top + window.scrollY + svgY * scale;
    const left = Math.min(x - tip.offsetWidth / 2, window.scrollX + document.documentElement.clientWidth - tip.offsetWidth - 8);
    tip.style.left = `${Math.max(left, window.scrollX + 8)}px`;
    tip.style.top = `${y - tip.offsetHeight - 12}px`;
  }

  function hide() { tip.hidden = true; }

  document.querySelectorAll(".viz[data-line]").forEach((chart) => {
    const cross = chart.querySelector(".viz-cross");
    const dot = chart.querySelector(".viz-cross-dot");
    function enter(event) {
      const hit = event.currentTarget;
      const x = hit.dataset.x, y = hit.dataset.y;
      cross.setAttribute("x1", x); cross.setAttribute("x2", x);
      cross.setAttribute("visibility", "visible");
      dot.setAttribute("cx", x); dot.setAttribute("cy", y);
      dot.setAttribute("visibility", "visible");
      show(hit, Number(x), Number(y));
    }
    function leave() {
      cross.setAttribute("visibility", "hidden");
      dot.setAttribute("visibility", "hidden");
      hide();
    }
    chart.querySelectorAll(".viz-hit").forEach((hit) => {
      hit.addEventListener("pointerenter", enter);
      hit.addEventListener("focus", enter);
      hit.addEventListener("pointerleave", leave);
      hit.addEventListener("blur", leave);
    });
  });

  document.querySelectorAll(".viz-fill[data-tip-value]").forEach((mark) => {
    function enter() {
      const box = mark.getBBox();
      show(mark, box.x + box.width / 2, box.y);
    }
    mark.addEventListener("pointerenter", enter);
    mark.addEventListener("focus", enter);
    mark.addEventListener("pointerleave", hide);
    mark.addEventListener("blur", hide);
  });

  // The filters apply as soon as they change; "Custom" shows its dates first.
  const form = document.getElementById("usage-filters");
  const custom = form.querySelector(".usage-custom");
  form.addEventListener("change", (event) => {
    const el = event.target;
    if (el.name === "range" && el.value === "custom") { custom.hidden = false; return; }
    if (el.type === "date" && el.value && Number(el.value.slice(0, 4)) < 1900) return;
    if (el.name === "range") custom.querySelectorAll("input").forEach((d) => { d.disabled = true; });
    form.requestSubmit();
  });
  window.addEventListener("pageshow", () => {
    form.querySelectorAll("input").forEach((el) => { el.disabled = false; });
  });
})();
