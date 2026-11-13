// Backdrop for the landing and About pages: real protein structures drifting slowly down
// behind the content. The twenty structures come from the RCSB PDB and are committed under
// static/phylobot/backdrop/ (see its README); the images are opaque JPEGs on white, so the
// white is keyed out here and what is left is tinted to suit the skin.
//
// The reader turns it on or off, and sets its speed and density, under Settings -> Display
// (phylobot.skins sets data-backdrop and the two dials on <html>). Only the Clean, Dark and
// Spectrum skins show one at all; on the others data-backdrop is "none" whatever the reader
// picked. Colours come from the skin as tokens on .backdrop (--backdrop-colors, -count,
// -alpha, -hue-by-x), with the muted fallback below. With prefers-reduced-motion the
// structures are scattered once and never move: the texture without the motion.
(() => {
  const FRAME_MS = 33;                  // ~30 fps; the motion is slow, so that is plenty
  const FALL = [26, 46];                // seconds to cross the viewport
  const FALLBACK = ["#8799a8", "#9a8fa3", "#8fa69a", "#a99a83"];

  let masks = null, pending = null;      // the structure images, loaded once and shared
  const still = matchMedia("(prefers-reduced-motion: reduce)").matches;

  const rand = (lo, hi) => lo + Math.random() * (hi - lo);
  const pick = (list) => list[Math.floor(Math.random() * list.length)];
  const ratio = () => window.devicePixelRatio || 1;

  function read(container) {
    const root = document.documentElement.dataset;
    const style = getComputedStyle(container);
    const number = (name, fallback) => {
      const value = parseFloat(style.getPropertyValue(name));
      return Number.isFinite(value) ? value : fallback;
    };
    const colors = style.getPropertyValue("--backdrop-colors")
      .split(",").map((c) => c.trim()).filter(Boolean);
    return {
      style: document.documentElement.dataset.backdrop || "none",
      colors: colors.length ? colors : FALLBACK,
      count: Math.round(number("--backdrop-count", 15)),
      alpha: number("--backdrop-alpha", 0.13),
      scale: number("--backdrop-scale", 1),
      hueByX: style.getPropertyValue("--backdrop-hue-by-x").trim() === "1",
      // The reader's own dials (Settings -> Display), as multipliers of all of the above.
      speed: Number(root.backdropSpeed) || 1,
      density: Number(root.backdropDensity) || 1,
    };
  }

  // ---- real structures from the PDB ----------------------------------------------

  // The images are JPEGs on white. Everything keeps its distance from white as alpha, so
  // the structure survives with its shading and the page shows through around it.
  function keyWhite(image) {
    const mask = document.createElement("canvas");
    mask.width = image.naturalWidth;
    mask.height = image.naturalHeight;
    const ctx = mask.getContext("2d", { willReadFrequently: true });
    ctx.drawImage(image, 0, 0);
    const pixels = ctx.getImageData(0, 0, mask.width, mask.height);
    const data = pixels.data;
    for (let i = 0; i < data.length; i += 4) {
      const ink = 255 - Math.min(data[i], data[i + 1], data[i + 2]);
      data[i] = data[i + 1] = data[i + 2] = 0;
      data[i + 3] = ink < 24 ? 0 : ink;        // JPEG ringing leaves a faint halo on white
    }
    ctx.putImageData(pixels, 0, 0);
    return mask;
  }

  function load(url) {
    if (masks || pending) return pending;
    // phylobot.com carries the list inline, each image as a data: address
    // (window.phylobotStructures, written by website/build.py): a page opened straight
    // from disk may neither fetch a file nor read the pixels of an image beside it.
    const inline = window.phylobotStructures;
    pending = (inline ? Promise.resolve(inline) : fetch(url).then((response) => response.json()))
      .then((list) => Promise.all(list.map((entry) => new Promise((resolve) => {
        const image = new Image();
        image.addEventListener("load", () => resolve(keyWhite(image)));
        image.addEventListener("error", () => resolve(null));
        image.src = entry.src || url.replace(/[^/]*$/, entry.file);
      }))))
      .then((loaded) => { masks = loaded.filter(Boolean); })
      .catch(() => { masks = []; });
    return pending;
  }

  function tint(mask, color, size) {
    const image = document.createElement("canvas");
    image.width = image.height = Math.ceil(size);
    const ctx = image.getContext("2d");
    ctx.drawImage(mask, 0, 0, image.width, image.height);
    ctx.globalCompositeOperation = "source-in";       // keep the shape, replace the colour
    ctx.fillStyle = color;
    ctx.fillRect(0, 0, image.width, image.height);
    return image;
  }

  // ---- one backdrop ----------------------------------------------------------

  // Each .backdrop element runs its own engine, so the Settings page can preview the
  // dials in a box beside them while the landing page fills the window.
  function createBackdrop(container) {
    let canvas = null, context = null, config = null;
    let shapes = [], timer = 0, last = 0;

    // ---- the falling shapes ---------------------------------------------------------

    // Where a structure of this size may fall, as a fraction of the width (so a resize
    // keeps the layout): anywhere, or to either side of the element named by data-clear
    // (the landing page's logo), clear of it by the structure's tumbling size and drift.
    function across(size, drift) {
      const target = container.dataset.clear && document.querySelector(container.dataset.clear);
      const width = container.clientWidth;
      const box = target && target.getBoundingClientRect();
      if (!box || !box.width || !width) return Math.random();   // absent or hidden by the skin
      const own = container.getBoundingClientRect();
      const margin = size * 0.71 + drift;          // half the diagonal, so rotation stays clear
      const left = Math.max(0, (box.left - own.left - margin) / width);
      const right = Math.max(0, 1 - (box.right - own.left + margin) / width);
      if (left + right < 0.1) return Math.random();  // a phone: no room beside it
      const u = Math.random() * (left + right);
      return u < left ? u : 1 - right + (u - left);
    }

    function spawn(y) {
      const size = rand(84, 160) * config.scale;
      const drift = rand(10, 50);
      const x = across(size, drift);
      const color = config.hueByX
        ? config.colors[Math.min(config.colors.length - 1, Math.floor(x * config.colors.length))]
        : pick(config.colors);
      const image = tint(pick(masks), color, size * ratio());
      const shape = {
        image, x, y,
        width: image.width / ratio(), height: image.height / ratio(),
        alpha: config.alpha * rand(0.55, 1),
        speed: container.clientHeight / rand(FALL[0], FALL[1]),   // px per second
        drift, driftSpeed: rand(0.05, 0.16), phase: rand(0, Math.PI * 2),
        angle: rand(0, Math.PI * 2), spin: rand(-0.06, 0.06),
      };
      if (y === null) shape.y = -shape.height / 2;
      return shape;
    }

    function paint(dt) {
      const width = container.clientWidth, height = container.clientHeight;
      context.clearRect(0, 0, width, height);
      shapes.forEach((shape, i) => {
        shape.y += shape.speed * dt;
        shape.phase += shape.driftSpeed * dt;
        shape.angle += shape.spin * dt;
        if (shape.y - shape.height > height) {
          shapes[i] = spawn(null);
          return;
        }
        context.globalAlpha = shape.alpha;
        context.save();
        context.translate(shape.x * width + Math.sin(shape.phase) * shape.drift, shape.y);
        context.rotate(shape.angle);
        context.drawImage(shape.image, -shape.width / 2, -shape.height / 2,
                          shape.width, shape.height);
        context.restore();
      });
      context.globalAlpha = 1;
    }

    function tick(now) {
      timer = requestAnimationFrame(tick);
      const elapsed = now - last;
      if (elapsed >= 0 && elapsed < FRAME_MS) return;  // a clock that went backwards must not stall it
      last = now;
      // A hidden tab comes back with a huge gap; the speed dial then scales every motion
      // at once -- falling, drift and tumble.
      paint(Math.min(Math.max(elapsed, 0) / 1000, 0.1) * config.speed);
    }

    function resize() {
      canvas.width = container.clientWidth * ratio();
      canvas.height = container.clientHeight * ratio();
      context.setTransform(ratio(), 0, 0, ratio(), 0, 0);
      if (shapes.length) paint(0);     // resizing clears the canvas
    }

    function populate() {
      // The skin's count suits a desktop window; a phone gets proportionally fewer, or the
      // same structures would crowd a third of the width.
      const room = Math.min(1, container.clientWidth / 1300);
      const count = Math.max(2, Math.round(config.count * room * config.density));
      shapes = Array.from({ length: count },
                             () => spawn(rand(-container.clientHeight * 0.3, container.clientHeight)));
      paint(0);                           // the first frame now, not whenever the loop first runs
      if (still) return;
      last = performance.now();
      timer = requestAnimationFrame(tick);
    }

    function start() {
      config = read(container);
      if (config.style === "none") return;
      canvas = document.createElement("canvas");
      container.append(canvas);
      context = canvas.getContext("2d");
      resize();
      window.addEventListener("resize", resize);
      // The images have to arrive first; if they cannot, the page simply has no backdrop.
      load(container.dataset.structures).then(() => {
        if (canvas && masks.length) populate();
      });
    }

    function stop() {
      if (!canvas) return;
      cancelAnimationFrame(timer);
      window.removeEventListener("resize", resize);
      canvas.remove();
      canvas = context = config = null;
      shapes = [];
      timer = 0;
    }

    return { start, stop };
  }

  const instances = [...document.querySelectorAll(".backdrop")].map(createBackdrop);
  if (!instances.length) return;      // every page but the landing, About and Settings pages
  instances.forEach((backdrop) => backdrop.start());
  // The menus change the style, the palette or a dial without reloading (skins.js sets the
  // attributes); the loaded structure images are kept across a restart.
  new MutationObserver(() => instances.forEach((backdrop) => { backdrop.stop(); backdrop.start(); }))
    .observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["data-skin", "data-backdrop", "data-backdrop-speed", "data-backdrop-density"],
    });
})();
