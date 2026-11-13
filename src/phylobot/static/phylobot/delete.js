// Shared confirmation dialog (_delete_dialog.html) for deleting jobs or searches and for
// running a job again. onConfirm runs only if the user presses the action button.
const DELETE_WORDING = {
  job: {
    one: "Its input and result files will be permanently removed. A record that it ran is kept.",
    many: "Their input and result files will be permanently removed. A record that they ran is kept.",
    plural: "jobs",
  },
  search: {
    one: "Its query and results will be permanently removed. A record that it ran is kept.",
    many: "Their queries and results will be permanently removed. A record that they ran is kept.",
    plural: "searches",
  },
  alignment: {
    one: "Its file will be permanently removed. Jobs that used it keep their own copy.",
    many: "Their files will be permanently removed. Jobs that used them keep their own copies.",
    plural: "alignments",
  },
};

function confirmDialog({ title, body, note, action, danger }, onConfirm) {
  const dialog = document.getElementById("confirm-delete");
  dialog.querySelector("#confirm-title").textContent = title;
  dialog.querySelector("#confirm-body").textContent = body;
  const noteLine = dialog.querySelector("#confirm-running");
  noteLine.hidden = !note;
  noteLine.textContent = note || "";
  const button = dialog.querySelector("#confirm-action");
  button.textContent = action;
  button.className = danger ? "danger" : "primary";
  dialog.returnValue = "";
  dialog.onclose = () => { if (dialog.returnValue === "confirm") onConfirm(); };
  dialog.showModal();
}

// `items` is a list of {name, status, dependents}; `kind` is "job", "search" or
// "alignment". `dependents` names the jobs that read results from a job (extensions,
// phylobot.joined): deleting it takes those results from them, which the note says.
function confirmDelete(items, onConfirm, kind = "job") {
  const words = DELETE_WORDING[kind];
  const one = items.length === 1;
  const unfinished = items.filter((item) => ["queued", "running"].includes(item.status)).length;
  const deleting = new Set(items.map((item) => item.name));
  const readers = [...new Set(items.flatMap((item) => item.dependents || []))]
    .filter((name) => !deleting.has(name));
  const notes = [];
  if (unfinished) {
    notes.push(unfinished === 1 ? `1 unfinished ${kind} will be stopped.`
      : `${unfinished} unfinished ${words.plural} will be stopped.`);
  }
  if (readers.length) {
    notes.push(`${readers.length === 1 ? "This job depends" : "These jobs depend"} on `
      + `${one ? "it" : "them"}, reading results from ${one ? "it" : "them"} instead of `
      + `keeping copies, and will lose those results: ${readers.map((n) => `“${n}”`).join(", ")}.`);
  }
  confirmDialog({
    title: one ? `Delete “${items[0].name}”?` : `Delete ${items.length} ${words.plural}?`,
    body: one ? words.one : words.many,
    note: notes.join(" "),
    action: "Delete",
    danger: true,
  }, onConfirm);
}

