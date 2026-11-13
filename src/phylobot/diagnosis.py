"""Why a job failed, in plain words, from what the job keeps (pure).

The job page's "Why it failed" box (`views.failure_context`) reads the job's one-line
error (set by the worker: an input error, the time limit, an interruption, or a pipeline
failure), its progress (the trace's tasks, with their exit status for newer jobs) and
its stored logs: the console log (pipeline.log, with Nextflow's "Error executing
process" report), Nextflow's own log, and each task's stderr. `diagnose` names the
step that failed, guesses the cause from what those logs commonly say (memory, disk,
the container runtime, a tool's own ERROR lines, ...), and says what to do; for
anything it doesn't recognise it still names the step and quotes its last lines.
Computed per request, so it covers every failed job, old ones included.
"""

import re
from dataclasses import dataclass, field, replace
from typing import Any

from phylobot.executors import GENERIC_FAILURE, INTERRUPTED, task_name
from phylobot.jobs import STOPPED
from phylobot.stages import LABELS

# The same two messages as stored before jobs could be resumed (older jobs keep them).
LEGACY_STOPS = (
    "The job was interrupted. Please resubmit.",
    "The job stopped unexpectedly. Please resubmit.",
)
DONE = {"COMPLETED", "CACHED"}
MAX_QUOTED = 12  # lines of a tool's output quoted
# Lines every task's stderr has on an ARM Mac (amd64 images), never the cause.
NOISE = re.compile(r"requested image's platform|^\s*$")


@dataclass(frozen=True)
class Diagnosis:
    kind: str  # input, interrupted, time_limit, memory, disk, runtime, image, tool, unknown
    title: str
    explanation: str
    advice: str
    step: str = ""  # the failed (or running) step, as the steps graph labels it
    quoted: list[str] = field(default_factory=list)  # the tool's own last lines
    inputs_ok: bool = False  # the input checks had passed
    rerun: bool = False  # running it again as it is may well work
    generic: bool = False  # the worker's own message says nothing more than this does


@dataclass
class Logs:
    """What a failed job left: the console log, Nextflow's log, task stderr by task."""

    pipeline: str = ""
    nextflow: str = ""
    tasks: dict[str, str] = field(default_factory=dict)  # "MODELFINDER (mafft)" -> stderr
    outputs: dict[str, str] = field(default_factory=dict)  # ... -> stdout (IQ-TREE's ERRORs)


def step_label(name: str) -> str:
    """ "MODELFINDER (mafft)" -> "ModelFinder (mafft)"."""
    name = task_name(name)
    process, _, tag = name.partition(" ")
    label = LABELS.get(process, process.replace("_", " ").title())
    return f"{label} {tag}".strip()


def failed_task(progress: dict[str, Any], logs: Logs) -> tuple[str, str]:
    """(the failed task's name, its exit status if known)."""
    for task in progress.get("tasks") or []:
        if task.get("status") == "FAILED":
            return task_name(task["name"]), str(task.get("exit") or "")
    if found := re.search(r"Error executing process > '([^']+)'", logs.pipeline + logs.nextflow):
        name = task_name(found.group(1))
        status = re.search(r"Command exit status:\s*\n\s*(-?\d+)", logs.pipeline)
        status = status or re.search(r"terminated with an error exit status \((-?\d+)\)",
                                     logs.nextflow)  # fmt: skip
        return name, status.group(1) if status else ""
    return "", ""


def running_task(progress: dict[str, Any]) -> str:
    """The step that was running (or about to) when the job stopped."""
    tasks = progress.get("tasks") or []
    for task in tasks:
        if task.get("status") not in DONE and task.get("status") != "FAILED":
            return task_name(task["name"])
    return task_name(progress.get("last") or "")


def inputs_passed(progress: dict[str, Any]) -> bool:
    return any(
        task_name(t.get("name", "")).startswith("VALIDATE") and t.get("status") in DONE
        for t in progress.get("tasks") or []
    )


def command_error(pipeline: str) -> list[str]:
    """Nextflow's "Command error:" block from the console log."""
    found = re.search(r"Command error:\n(.*?)(?:\n\s*\nWork dir:|\n\nTip:|\Z)", pipeline, re.S)
    return found.group(1).splitlines() if found else []


def quoted(lines: list[str]) -> list[str]:
    kept = [ln.rstrip() for ln in lines if not NOISE.search(ln)]
    return kept[-MAX_QUOTED:]


def tool_errors(text: str) -> list[str]:
    """A tool's own error lines (IQ-TREE, RAxML-NG and our steps all write them)."""
    return [
        ln.strip() for ln in text.splitlines()
        if re.match(r"\s*(ERROR|Error|Alignment error|Input error|Species error)\b", ln)
    ]  # fmt: skip


# (kind, pattern, title, explanation, advice, rerun): checked in order on the failed
# task's stderr, the command error and Nextflow's report.
PATTERNS = [
    (
        "memory",
        r"exit status \(137\)|^137$|\bKilled\b|bad_alloc|MemoryError|Cannot allocate memory|"
        r"[Oo]ut of memory|OutOfMemory",
        "The step ran out of memory",
        "The computer (or container) it ran on had less memory than the step needed, and the "
        "system stopped it.",
        "Try fewer or shorter sequences, or fewer models (memory grows with both); on your own "
        "machine, give Docker more memory.",
        False,
    ),
    (
        "disk",
        r"No space left on device|Disk quota exceeded",
        "The disk filled up",
        "There was no room left to write the step's files.",
        "Free some disk space (on your own machine, `make clean` removes regenerable scratch), "
        "then run the job again.",
        True,
    ),
    (
        "cpus",
        r"Requested CPUs are not available|for cpuset cpus",
        "The machine's CPU setting is wrong",
        "PHYLOBOT_DOCKER_CPUSET in .env names CPUs that Docker cannot give this step, so its "
        "container could not start. That is a problem with the machine, not with your job.",
        "Run `make docker-cpus` on the machine, put the value it suggests in .env, restart "
        "the worker, then run the job again.",
        True,
    ),
    (
        "runtime",
        r"IllegalThreadStateException|process hasn't exited|Cannot connect to the Docker "
        r"daemon|Error response from daemon|docker: error|context deadline exceeded",
        "The container runtime stopped responding",
        "Docker, which runs each step, failed or stopped answering while this step ran. That "
        "is a problem with the machine, not with your job.",
        "Run the job again (Duplicate and run again). On your own machine, check that Docker "
        "is running and responsive first.",
        True,
    ),
    (
        "esr_model",
        r"ValueError: fixed model .* gives lnL|ValueError: cannot fix model component",
        "ESR could not fix the model's parameters",
        "ESR reruns IQ-TREE with every parameter of the model fixed at its fitted value, and "
        "first checks that this gives the fit's log-likelihood again. Here it did not (or the "
        "model has a component it cannot fix), so the model was misread and nothing was "
        "reconstructed with it.",
        "Please report it, with the model named above. Another combination of the job, or a "
        "simpler model on an uploaded tree, will usually run.",
        False,
    ),
    (
        "image",
        r"Unable to find image|pull access denied|manifest unknown|failed to pull|"
        r"repository does not exist",
        "A tool's container could not be fetched",
        "The step's pinned container image could not be downloaded.",
        "Check the network connection and run the job again.",
        True,
    ),
]


def diagnose(error: str, progress: dict[str, Any], logs: Logs) -> Diagnosis | None:
    """Why the job failed (None if there is nothing to say: no error recorded)."""
    if not error:
        return None
    found = cause(error, progress, logs)
    return replace(found, generic=error == GENERIC_FAILURE)


def cause(error: str, progress: dict[str, Any], logs: Logs) -> Diagnosis:
    ok = inputs_passed(progress)
    if error.startswith("Input error:"):
        return Diagnosis(
            "input", "A problem with your input", error.removeprefix("Input error:").strip(),
            "Fix it and run the job again (Duplicate and run again keeps your settings).",
            step="Check input",
        )  # fmt: skip
    if error in (INTERRUPTED, STOPPED, *LEGACY_STOPS):
        running = running_task(progress)
        where = f" during “{step_label(running)}”" if running else ""
        return Diagnosis(
            "interrupted", "PhyloBot stopped while your job was running",
            f"The worker that runs jobs was stopped or restarted{where} (for example by a "
            "server update or restart), so the job could not finish. Nothing points to a "
            "problem with the job itself." + (" Your input had passed its checks." if ok else ""),
            "Resume it with the button below; nothing about the job itself needs changing.",
            step=step_label(running) if running else "", inputs_ok=ok, rerun=True,
        )  # fmt: skip
    if error.startswith("The job exceeded the time limit"):
        running = running_task(progress)
        return Diagnosis(
            "time_limit", "The job ran out of time", error,
            "Ask for less work (fewer sequences, aligners, tree sources or models, or fewer "
            "branch-support replicates); the work estimate on the New Job page shows what "
            "each choice costs.",
            step=step_label(running) if running else "", inputs_ok=ok,
        )  # fmt: skip
    name, status = failed_task(progress, logs)
    stderr = logs.tasks.get(name, "")
    if not quoted(stderr.splitlines()):  # nothing but noise: the tool wrote to stdout
        stderr = logs.outputs.get(name, "") or stderr
    report = "\n".join(command_error(logs.pipeline))
    caused = re.search(r"Caused by:\n(.*?)\n\n", logs.pipeline, re.S)
    text = "\n".join([stderr, report, caused.group(1) if caused else "",
                      f"exit status ({status})" if status else ""])  # fmt: skip
    step = step_label(name) if name else ""
    for kind, pattern, title, explanation, advice, rerun in PATTERNS:
        if re.search(pattern, text, re.M):
            lines = quoted(stderr.splitlines() or report.splitlines())
            return Diagnosis(kind, title, explanation, advice, step, lines, ok, rerun)
    if errors := tool_errors(stderr + "\n" + report):
        return Diagnosis(
            "tool", f"{step.split(' (')[0] or 'A step'} reported an error",
            "The program running this step stopped with the message below.",
            "If the message is about your data (sequence names, lengths, characters, too "
            "few sequences), fix it and run again; otherwise please report it, with the log.",
            step, quoted(errors), ok,
        )  # fmt: skip
    if error == GENERIC_FAILURE or name:
        lines = quoted(stderr.splitlines() or report.splitlines())
        status_text = f" (exit status {status})" if status else ""
        return Diagnosis(
            "unknown", "The pipeline failed" + (f" at “{step}”" if step else ""),
            f"The step stopped with an error{status_text} that PhyloBot does not recognise."
            + (" Its last lines are below." if lines else ""),
            "Running it again may work if the cause was temporary; if it fails the same way, "
            "please report it with the log (Download ZIP has everything).",
            step, lines, ok,
        )  # fmt: skip
    return Diagnosis("unknown", "The job failed", error, "", inputs_ok=ok)
