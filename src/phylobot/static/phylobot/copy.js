// Copy buttons beside command snippets: copy the snippet's text, say so briefly.
document.addEventListener("click", async (event) => {
  const button = event.target.closest("button.copy");
  if (!button) return;
  const code = button.closest(".command").querySelector("code");
  try {
    await navigator.clipboard.writeText(code.textContent);
    button.textContent = "Copied";
  } catch {
    button.textContent = "Select and copy";
  }
  setTimeout(() => { button.textContent = "Copy"; }, 1500);
});
