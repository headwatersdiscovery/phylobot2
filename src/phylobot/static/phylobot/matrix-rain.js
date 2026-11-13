// Matrix skin: residues falling behind the top bar, drawn on a canvas only while
// data-skin="matrix" (switching skins starts or stops it). The letters are the twenty
// amino acids. Off when the reader prefers reduced motion, and paused in hidden tabs.
(() => {
  const LETTERS = "ACDEFGHIKLMNPQRSTVWY";
  const STEP = 14;           // column width and row height, px
  const FRAME_MS = 70;
  let canvas = null, timer = null, drops = [];

  function resize() {
    const header = canvas.parentElement;
    const ratio = window.devicePixelRatio || 1;
    canvas.width = header.clientWidth * ratio;
    canvas.height = header.clientHeight * ratio;
    canvas.getContext("2d").setTransform(ratio, 0, 0, ratio, 0, 0);
    const columns = Math.ceil(header.clientWidth / STEP);
    const rows = Math.ceil(header.clientHeight / STEP);
    drops = Array.from({ length: columns }, (_, i) => drops[i] ?? -Math.floor(Math.random() * rows * 3));
  }

  function frame() {
    if (document.hidden) return;
    const context = canvas.getContext("2d");
    const width = canvas.width, height = canvas.height;
    context.fillStyle = "rgba(2, 4, 2, 0.22)";   // fade what was drawn before: the trails
    context.fillRect(0, 0, width, height);
    context.font = `${STEP - 2}px ui-monospace, Menlo, monospace`;
    const rows = Math.ceil(canvas.parentElement.clientHeight / STEP);
    drops.forEach((y, i) => {
      if (y >= 0) {
        context.fillStyle = Math.random() < 0.08 ? "#d8ffe2" : "#3dff73";
        context.fillText(LETTERS[Math.floor(Math.random() * LETTERS.length)], i * STEP + 2, (y + 1) * STEP - 3);
      }
      drops[i] = y > rows + Math.random() * 12 ? -Math.floor(Math.random() * rows * 2) : y + 1;
    });
  }

  function start() {
    // The top bar; on phylobot.com, which has none, the element that asks for the rain.
    const header = document.querySelector("body > header") || document.querySelector("[data-matrix-rain]");
    if (canvas || !header || matchMedia("(prefers-reduced-motion: reduce)").matches) return;
    canvas = document.createElement("canvas");
    canvas.className = "matrix-rain";
    canvas.setAttribute("aria-hidden", "true");
    header.prepend(canvas);
    resize();
    window.addEventListener("resize", resize);
    timer = setInterval(frame, FRAME_MS);
  }

  function stop() {
    if (!canvas) return;
    clearInterval(timer);
    window.removeEventListener("resize", resize);
    canvas.remove();
    canvas = null;
    drops = [];
  }

  const update = () => (document.documentElement.dataset.skin === "matrix" ? start() : stop());
  new MutationObserver(update).observe(document.documentElement, { attributes: true, attributeFilter: ["data-skin"] });
  update();
})();
