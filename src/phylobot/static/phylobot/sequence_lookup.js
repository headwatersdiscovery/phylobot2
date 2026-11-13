// "Find a sequence by name" on the Find homologs page (homologs.html):
//   1. As-you-type suggestions for the Species and Protein/gene fields, filled into their
//      <datalist>s from /homologs/suggest/ (UniProt). The field keeps working if a request
//      fails — suggestions are a convenience, not a gate.
//   2. "Use this sequence" drops a looked-up sequence into the homolog search box. The
//      results are swapped in by HTMX, so that listener sits on the document.
(() => {
  const form = document.querySelector(".sequence-lookup-form");
  const suggestUrl = form && form.dataset.suggestUrl;

  if (suggestUrl) {
    const DEBOUNCE = 220; // ms after the last keystroke before asking
    const MIN = 2; // characters before asking

    // Replace a datalist's options with freshly fetched ones.
    const fill = (list, items) => {
      list.textContent = "";
      for (const item of items) {
        const option = document.createElement("option");
        option.value = item.value;
        if (item.label && item.label !== item.value) option.label = item.label;
        list.appendChild(option);
      }
    };

    for (const field of form.querySelectorAll("input[data-suggest]")) {
      const list = document.getElementById(field.dataset.suggestList);
      if (!list) continue;
      let timer = null;
      let token = 0; // ignore answers that arrive out of order
      field.addEventListener("input", () => {
        const value = field.value.trim();
        if (timer) clearTimeout(timer);
        if (value.length < MIN) return;
        timer = setTimeout(() => {
          const mine = ++token;
          const params = new URLSearchParams({ kind: field.dataset.suggest, q: value });
          const speciesField = field.dataset.suggestSpecies
            ? document.getElementById(field.dataset.suggestSpecies)
            : null;
          if (speciesField && speciesField.value.trim()) params.set("species", speciesField.value.trim());
          fetch(`${suggestUrl}?${params.toString()}`, { headers: { Accept: "application/json" } })
            .then((response) => (response.ok ? response.json() : { items: [] }))
            .then((data) => {
              if (mine === token) fill(list, (data && data.items) || []);
            })
            .catch(() => {
              /* leave the datalist as it is: suggestions are optional */
            });
        }, DEBOUNCE);
      });
    }
  }

  // Keep the search consent naming the service the user picked (EBI / NCBI). The radio
  // values are "service:database", so the service is the part before the colon.
  const serviceLabel = document.querySelector("[data-search-service]");
  if (serviceLabel) {
    const update = () => {
      const chosen = document.querySelector('input[name="database"]:checked');
      if (chosen) serviceLabel.textContent = chosen.value.split(":")[0].toUpperCase();
    };
    for (const radio of document.querySelectorAll('input[name="database"]')) {
      radio.addEventListener("change", update);
    }
    update();
  }

  document.addEventListener("click", (event) => {
    const button = event.target.closest(".use-sequence");
    if (!button) return;
    const query = document.getElementById("id_query");
    const name = document.getElementById("id_name");
    if (query) {
      query.value = button.dataset.sequence || "";
      query.dispatchEvent(new Event("input", { bubbles: true }));
    }
    if (name && !name.value) name.value = button.dataset.name || "";
    const section = query ? query.closest("section") : null;
    if (section) section.scrollIntoView({ behavior: "smooth", block: "start" });
    if (query) query.focus({ preventScroll: true });
  });
})();
