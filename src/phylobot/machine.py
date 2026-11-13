"""The machine the worker runs pipelines on, and each job's share of it.

Nextflow's local executor counts what tasks *reserve* (1 CPU and the step memory each), never
what they use, and on its own it assumes the whole machine is its own. Two jobs at once
would each plan for all of it. So the worker tells every job, at launch, its share:

    CPUs   = usable CPUs // slots
    memory = (memory - reserved) / slots - HEAD_MEMORY

where `slots` is how many jobs can run at once. Each total comes from, in order: what the
worker detects (Docker's own figures), else .env (PHYLOBOT_CPUS, PHYLOBOT_MEMORY_GB),
else the Settings page; the memory to leave for everything else (the web stack, the
desktop) and what each step reserves likewise from the built-in default, .env or the
page. Steps run at once = the share's memory // the step reservation, at most its
CPUs. The share is written into
the attempt's own Nextflow config (resources.config), never params.json: a resume runs
on whatever machine resumes it, and params.json must not change (every cache key would).

Detection runs at every launch and in the worker's maintenance (WorkerMachine, which
the Settings page shows), so a BIOS or Docker settings change is picked up without
anyone re-running anything.
"""

import logging
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from django.conf import settings
from django.utils import timezone

from phylobot import site_settings
from phylobot.models import WorkerMachine

logger = logging.getLogger(__name__)

GIB = 1024**3
MIB = 1024**2
# What a pipeline task reserves unless the Settings page or PHYLOBOT_TASK_MEMORY_GB says
# otherwise: nextflow.config's process.memory, for runs without the worker (pinned by a
# test). The worker's resources.config replaces it with the setting. No job's share may be
# smaller than one task's: a task asking for more than its executor has fails.
DEFAULT_TASK_MEMORY = 2 * GIB
# Each job's Nextflow (a JVM in the runner container) holds this outside its tasks:
# 0.4-0.5 GB measured on a 43-task job, with room for larger DAGs.
HEAD_MEMORY = 1 * GIB
DOCKER_TIMEOUT = 20  # seconds; Docker Desktop sometimes stalls


@dataclass(frozen=True)
class Detected:
    """What the worker found. `cpus` are Docker's minus isolated ones (or those
    PHYLOBOT_DOCKER_CPUSET names); `note` says what went wrong if Docker did not answer."""

    cpus: int
    isolated: str  # isolcpus, as the kernel lists them ("" for none)
    memory_bytes: int
    note: str = ""


@dataclass(frozen=True)
class Totals:
    """The machine as PhyloBot treats it, with where each figure comes from."""

    cpus: int
    cpus_source: str
    memory_bytes: int
    memory_source: str
    reserved_bytes: int
    reserved_source: str
    task_bytes: int = DEFAULT_TASK_MEMORY  # what each step reserves
    task_source: str = "built-in default"


@dataclass(frozen=True)
class Share:
    """One job's share: what Nextflow may reserve at once, and why it is that."""

    cpus: int
    memory_bytes: int
    slots: int
    warning: str = ""
    task_bytes: int = DEFAULT_TASK_MEMORY  # what each step reserves

    @property
    def tasks_at_once(self) -> int:
        return max(1, min(self.cpus, self.memory_bytes // self.task_bytes))


def cpu_list(text: str) -> set[int]:
    """ "0-3,6" -> {0, 1, 2, 3, 6}; the kernel's and Docker's format."""
    out: set[int] = set()
    for part in filter(None, text.strip().split(",")):
        low, _, high = part.partition("-")
        out.update(range(int(low), int(high or low) + 1))
    return out


def usable_cpus(online: int, isolated: str, cpuset: str) -> int:
    """CPUs task containers can really use: those PHYLOBOT_DOCKER_CPUSET names, else
    every online one but the isolated (a task stuck on an isolated CPU shares it)."""
    if cpuset:
        return len(cpu_list(cpuset))
    return max(1, online - len(cpu_list(isolated)))


def docker_info() -> tuple[int, int]:
    """(CPUs, memory in bytes) of the machine Docker runs containers on: the host on
    Linux, the VM on a Mac."""
    out = subprocess.run(
        ["docker", "info", "--format", "{{.NCPU}} {{.MemTotal}}"],
        capture_output=True, text=True, timeout=DOCKER_TIMEOUT, check=True,
    )  # fmt: skip
    cpus, memory = out.stdout.split()
    return int(cpus), int(memory)


def meminfo_bytes(text: str) -> int:
    match = re.search(r"^MemTotal:\s+(\d+) kB", text, re.MULTILINE)
    return int(match.group(1)) * 1024 if match else 0


def detect() -> Detected:
    """The worker shares Docker's kernel (its socket's), so /sys is the same one tasks
    see; Docker's own figures are preferred, and the worker's used if it does not answer."""
    isolated = ""
    if (path := Path("/sys/devices/system/cpu/isolated")).is_file():
        isolated = path.read_text().strip()
    note = ""
    try:
        online, memory = docker_info()
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        note = f"Docker did not answer ({type(error).__name__}); the worker's own figures"
        online = os.cpu_count() or 1
        memory = (
            meminfo_bytes(Path("/proc/meminfo").read_text())
            if Path("/proc/meminfo").is_file()
            else 0
        )
    cpus = usable_cpus(online, isolated, settings.PHYLOBOT_DOCKER_CPUSET)
    return Detected(cpus=cpus, isolated=isolated, memory_bytes=memory, note=note)


def record(found: Detected) -> None:
    """Keep what this worker found, for the Settings page (one row per worker host)."""
    WorkerMachine.objects.update_or_create(
        host=os.uname().nodename,
        defaults={
            "cpus": found.cpus,
            "isolated": found.isolated,
            "memory_bytes": found.memory_bytes,
            "note": found.note,
            "detected_at": timezone.now(),
        },
    )


def latest() -> WorkerMachine | None:
    """The machine a worker last reported (the web tier cannot ask Docker itself)."""
    return WorkerMachine.objects.order_by("-detected_at").first()


def totals(found: Detected | None) -> Totals:
    """Detected, then .env, then the Settings page: each later one wins when set, but
    never above what was detected (the page refuses more; .env, or a machine that has
    since shrunk, is capped here)."""
    cpus, cpus_source = (found.cpus, "detected") if found else (0, "not detected yet")
    memory, memory_source = (found.memory_bytes, "detected") if found else (0, "not detected yet")
    if settings.PHYLOBOT_CPUS:
        cpus, cpus_source = settings.PHYLOBOT_CPUS, "from PHYLOBOT_CPUS"
    if settings.PHYLOBOT_MEMORY_GB:
        memory, memory_source = int(settings.PHYLOBOT_MEMORY_GB * GIB), "from PHYLOBOT_MEMORY_GB"
    if (chosen := site_settings.override("machine_cpus")) is not None:
        cpus, cpus_source = int(chosen), "set on the Settings page"
    if (chosen := site_settings.override("machine_memory_gb")) is not None:
        memory, memory_source = int(float(chosen) * GIB), "set on the Settings page"
    if found is not None and cpus > found.cpus:
        cpus, cpus_source = found.cpus, f"{cpus_source}, capped at the {found.cpus} detected"
    if found is not None and memory > found.memory_bytes:
        memory = found.memory_bytes
        memory_source = f"{memory_source}, capped at the {gb(found.memory_bytes)} detected"
    row = site_settings.stored()
    reserved_gb = site_settings.reserved_memory_gb()
    reserved_source = site_settings.source("reserved_memory_gb", row)
    task_bytes = int(site_settings.task_memory_gb() * GIB)
    task_source = site_settings.source("task_memory_gb", row)
    return Totals(
        cpus, cpus_source, memory, memory_source, int(reserved_gb * GIB), reserved_source,
        task_bytes, task_source,
    )  # fmt: skip


def detected_limits() -> tuple[int, float] | None:
    """(CPUs, memory in GB rounded down to 0.1) the Settings page may not exceed: what a
    worker last detected. None before any worker has looked: nothing to check against."""
    found = latest()
    if found is None:
        return None
    return found.cpus, int(found.memory_bytes / GIB * 10) / 10


def check_overrides(
    cpus: int | None,
    memory_gb: float | None,
    reserved_gb: float | None,
    task_gb: float | None = None,
) -> dict[str, str]:
    """The Settings page's errors, by field: no more CPUs or memory than detected, and a
    reserve and step reservation that leave room for at least one step."""
    errors: dict[str, str] = {}
    limits = detected_limits()
    if limits is not None:
        most_cpus, most_gb = limits
        if cpus is not None and cpus > most_cpus:
            errors["machine_cpus"] = (
                f"This machine has {most_cpus} CPUs for pipeline steps (as detected): at most "
                f"{most_cpus}."
            )
        if memory_gb is not None and memory_gb > most_gb:
            errors["machine_memory_gb"] = (
                f"This machine has {most_gb:g} GB (as detected): at most {most_gb:g}."
            )
    if reserved_gb is not None or task_gb is not None:
        if memory_gb is None:
            memory_gb = settings.PHYLOBOT_MEMORY_GB or (limits[1] if limits else None)
        # A blank field on the page means unset: the .env or built-in value applies.
        step = task_gb if task_gb is not None else settings.PHYLOBOT_TASK_MEMORY_GB
        reserve = reserved_gb if reserved_gb is not None else settings.PHYLOBOT_RESERVED_MEMORY_GB
        room = step + HEAD_MEMORY / GIB
        if memory_gb is not None and reserve > memory_gb - room:
            # Name the field the user set; if both, the reserve (the larger lever).
            field = "reserved_memory_gb" if reserved_gb is not None else "task_memory_gb"
            if field == "reserved_memory_gb":
                errors[field] = (
                    f"That leaves too little of the {memory_gb:g} GB for one job with one step "
                    f"({room:g} GB): at most {max(memory_gb - room, 0):g}."
                )
            else:
                most = max(memory_gb - reserve - HEAD_MEMORY / GIB, 0)
                errors[field] = (
                    f"One step this large does not fit in the {memory_gb:g} GB after the "
                    f"{reserve:g} GB reserve and {HEAD_MEMORY / GIB:g} GB for Nextflow: at "
                    f"most {most:g}."
                )
    return errors


def slots() -> int:
    """How many jobs can run at once here: the worker's limit, and in local mode (one
    user) also that user's own, so a lone job is never held to half a machine."""
    count = site_settings.worker_concurrency()
    if settings.PHYLOBOT_MODE == "local":
        from phylobot.models import limits_for
        from phylobot.modes import local_user

        count = min(count, limits_for(local_user()).max_concurrent_jobs)
    return max(1, count)


def share(machine: Totals, jobs_at_once: int) -> Share:
    """One job's share. Never below one task's reservation (`machine.task_bytes`), or no
    task could start: the warning says the machine is overcommitted instead."""
    jobs_at_once = max(1, jobs_at_once)
    cpus = max(1, machine.cpus // jobs_at_once)
    free = machine.memory_bytes - machine.reserved_bytes
    memory = free // jobs_at_once - HEAD_MEMORY
    task = machine.task_bytes
    warning = ""
    if memory < task:
        warning = (
            f"Too little memory for {jobs_at_once} job{'s' * (jobs_at_once > 1)} at once: "
            f"{gb(free)} after the reserve is less than "
            f"{gb((task + HEAD_MEMORY) * jobs_at_once)}. Each job still gets one "
            "task's worth, so together they may use more than the machine has."
        )
        memory = task
    return Share(
        cpus=cpus, memory_bytes=memory, slots=jobs_at_once, warning=warning, task_bytes=task
    )


def for_launch() -> Share:
    """Detect now, record it, and work out the share for a job starting now."""
    found = detect()
    try:
        record(found)
    except Exception:
        logger.warning("could not record the detected machine", exc_info=True)
    return share(totals(found), slots())


def nextflow_config(job_share: Share) -> str:
    """The attempt's resources.config (passed with -c): this job's share for Nextflow's
    local executor, what each step reserves (the setting, in place of nextflow.config's
    default), and the retry's memory capped to the share. Numbers only."""
    megabytes = int(job_share.memory_bytes // MIB)
    task = int(job_share.task_bytes // MIB)
    return (
        "// Written by the worker for this attempt (phylobot.machine): this job's share of\n"
        f"// the machine, with {job_share.slots} job(s) able to run at once.\n"
        "executor {\n"
        f"    cpus   = {int(job_share.cpus)}\n"
        f"    memory = '{megabytes} MB'\n"
        "}\n"
        "process {\n"
        "    // nextflow.config doubles a task's memory when it is retried after running\n"
        "    // out; never past the job's share, which Nextflow would refuse to schedule.\n"
        f"    memory = {{ [{task}.MB * task.attempt, {megabytes}.MB].min() }}\n"
        "}\n"
    )


def gb(n: int) -> str:
    return f"{n / GIB:.1f} GB"
