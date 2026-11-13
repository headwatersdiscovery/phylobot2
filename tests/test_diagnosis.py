"""Why a job failed, from its logs (phylobot.diagnosis)."""

from typing import Any

from phylobot.diagnosis import Logs, diagnose
from phylobot.executors import GENERIC_FAILURE, INTERRUPTED
from phylobot.jobs import STOPPED

NOISE = (
    "WARNING: The requested image's platform (linux/amd64) does not match the detected host "
    "platform (linux/arm64/v8) and no specific platform was requested\n"
)


def tasks(*rows: tuple[str, str, str]) -> dict[str, Any]:
    return {"tasks": [{"name": n, "status": s, "exit": e} for n, s, e in rows]}


VALIDATED = ("PHYLOBOT:VALIDATE (1)", "COMPLETED", "0")


def console(process: str, caused: str, status: str = "", error: str = "") -> str:
    """The console log's report of a failed process, as Nextflow writes it."""
    text = f"ERROR ~ Error executing process > '{process}'\n\nCaused by:\n  {caused}\n\n"
    if status:
        text += f"Command executed:\n\n  iqtree3 ...\n\nCommand exit status:\n  {status}\n\n"
    if error:
        text += f"Command output:\n  (empty)\n\nCommand error:\n{error}\n\nWork dir:\n  /w\n"
    return text


def test_input_errors_and_interruptions() -> None:
    found = diagnose("Input error: need at least 4 sequences, found 2", {}, Logs())
    assert (
        found
        and found.kind == "input"
        and found.explanation == "need at least 4 sequences, found 2"
    )
    progress = tasks(VALIDATED, ("PHYLOBOT:MODELFINDER (provided1)", "COMPLETED", "0"),
                     ("PHYLOBOT:TREE_IQTREE (provided1/modelfinder)", "RUNNING", ""))  # fmt: skip
    for error in (INTERRUPTED, STOPPED):
        found = diagnose(error, progress, Logs())
        assert found and found.kind == "interrupted" and found.rerun and found.inputs_ok
        assert found.step == "IQ-TREE search (provided1/modelfinder)"
        assert "Your input had passed its checks" in found.explanation and not found.generic
    late = diagnose("The job exceeded the time limit (24 h).", progress, Logs())
    assert late and late.kind == "time_limit" and not late.rerun


def test_out_of_memory_from_the_exit_status() -> None:
    progress = tasks(VALIDATED, ("PHYLOBOT:TREE_RAXML (mafft/LG+G4)", "FAILED", "137"))
    found = diagnose(GENERIC_FAILURE, progress, Logs(tasks={"TREE_RAXML (mafft/LG+G4)": NOISE}))
    assert found and found.kind == "memory" and found.step == "RAxML-NG search (mafft/LG+G4)"
    assert found.generic and found.quoted == []


def test_container_runtime_and_disk_from_the_logs() -> None:
    # A real failure (2026-10-14): Docker stopped answering during ModelFinder.
    log = console("PHYLOBOT:MODELFINDER (mafft)", "process hasn't exited")
    found = diagnose(GENERIC_FAILURE, tasks(VALIDATED), Logs(pipeline=log))
    assert found and found.kind == "runtime" and found.step == "ModelFinder (mafft)" and found.rerun
    # Docker's own words when PHYLOBOT_DOCKER_CPUSET names CPUs the VM lacks (2026-10-25).
    cpus = console("PHYLOBOT:ALIGN_MAFFT (mafft)", "Process terminated", "125",
                   "docker: Error response from daemon: Requested CPUs are not available - "
                   "requested 0-12, available: 0-7")  # fmt: skip
    found = diagnose(GENERIC_FAILURE, tasks(VALIDATED), Logs(pipeline=cpus))
    assert found and found.kind == "cpus" and found.rerun
    full = console("PHYLOBOT:EXPORT (mafft/iqtree/LG+G4)", "Process terminated", "1",
                   "OSError: [Errno 28] No space left on device")  # fmt: skip
    found = diagnose(GENERIC_FAILURE, {}, Logs(pipeline=full))
    assert found and found.kind == "disk" and found.step == "Export results (mafft/iqtree/LG+G4)"


def test_a_tool_s_own_error_lines_are_quoted() -> None:
    # IQ-TREE writes its ERROR to stdout; stderr holds only noise.
    name = "MODELFINDER (mafft)"
    stdout = "IQ-TREE version 3.1.3\nReading alignment file...\nERROR: Sequence s0004 is too long\n"
    log = console(f"PHYLOBOT:{name}", "Process terminated with an error exit status (2)", "2")
    found = diagnose(GENERIC_FAILURE, tasks(VALIDATED, (f"PHYLOBOT:{name}", "FAILED", "2")),
                     Logs(pipeline=log, tasks={name: NOISE}, outputs={name: stdout}))  # fmt: skip
    assert found and found.kind == "tool" and found.title == "ModelFinder reported an error"
    assert found.quoted == ["ERROR: Sequence s0004 is too long"]


def test_anything_else_still_names_the_step_and_its_last_lines() -> None:
    name = "GAPS_GRASP (mafft/iqtree/LG+G4)"
    stderr = NOISE + "Exception in thread main\n  at bnkit.Main(Main.java:12)\n"
    found = diagnose(GENERIC_FAILURE, tasks(VALIDATED, (f"PHYLOBOT:{name}", "FAILED", "1")),
                     Logs(tasks={name: stderr}))  # fmt: skip
    assert found and found.kind == "unknown" and "exit status 1" in found.explanation
    assert found.step.endswith("(mafft/iqtree/LG+G4)")
    assert found.quoted == ["Exception in thread main", "  at bnkit.Main(Main.java:12)"]
    assert diagnose("", {}, Logs()) is None


def test_an_esr_model_that_would_not_fix() -> None:
    """phylobot-esr fit refuses a fixed model that does not give the fit's lnL back
    (phylobot.esr.Combination); the traceback's last line, as the task's stderr has it."""
    error = (
        "Traceback (most recent call last):\n  ...\nValueError: fixed model "
        "LG+F{0.07,...}+R4{0.25,0.15} gives lnL -4414.9, the fit -4414.5\n"
    )
    progress = tasks(("ESR:ESR_FIT (fit)", "FAILED", "1"))
    caused = "Process `ESR:ESR_FIT (fit)` terminated with an error exit status (1)"
    logs = Logs(pipeline=console("ESR:ESR_FIT (fit)", caused, "1", error),
                tasks={"ESR_FIT (fit)": error})  # fmt: skip
    found = diagnose(GENERIC_FAILURE, progress, logs)
    assert found and found.kind == "esr_model"
