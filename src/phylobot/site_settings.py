"""Site-wide settings that the Settings page can change while PhyloBot runs.

Each value comes from, in order:
1. the Settings page (the SiteSettings row), if a value is set there;
2. otherwise the environment variable or built-in default, as Django's settings hold it.

So .env still configures a fresh install, and the page overrides it without restarting
containers. Read values through the functions here, never from django.conf.settings
directly, or a change on the page is silently ignored. Values are read from the database
on every call (one small query): workers and the search runner pick up changes at
their next step.

Not here, on purpose: anything that must be fixed before the app starts or is secret
(PHYLOBOT_MODE, DEBUG, ALLOWED_HOSTS, the secret key, database and storage credentials).
"""

import os
from collections.abc import Callable
from dataclasses import dataclass
from datetime import timedelta
from typing import Any

from django.conf import settings

from phylobot import calibration
from phylobot.models import SiteSettings, UserCalibration

MEGABYTE = 1_000_000


@dataclass(frozen=True)
class Setting:
    field: str  # SiteSettings field
    label: str
    help: str
    group: str
    fallback: Callable[[], Any]  # the environment's or built-in value, in the page's unit
    env: str = ""  # the environment variable that can set it, if any
    unit: str = ""
    low: int = 0  # bounds for numbers on the page
    high: int = 0
    zero: str = ""  # what 0 means, shown in place of "0 <unit>"
    switch: bool = False  # on/off, not a number
    real: bool = False  # a decimal number, not a whole one
    detected: bool = False  # without the environment's value, what the worker detects


def _detected(field: str, env: str, unit: float = 1) -> Callable[[], Any]:
    """The environment's value if set, else what a worker last detected (None before any
    has), in the page's unit."""

    def fallback() -> Any:
        if chosen := getattr(settings, env):
            return chosen
        from phylobot.machine import latest  # it imports this module

        found = latest()
        if found is None:
            return None
        detected = getattr(found, field)
        return detected if unit == 1 else round(detected / unit, 1)

    return fallback


def _limit(name: str) -> Callable[[], Any]:
    return lambda: settings.PHYLOBOT_DEFAULT_LIMITS[name]


SEARCH, USERS, JOBS = "Homolog Search", "Limits for Every User", "Jobs and Worker"
SIMULATIONS, BILLING, IMPORT = "Simulations", "Billing", "Importing Jobs"
CALIBRATION = "Calibrated Posterior Probabilities"
# Built-in defaults, pending the owner's confirmation.
SIMULATION_MAX_REPLICATES = 100
SIMULATION_MAX_MEGARESIDUES = 50
# Built-in limits on an imported ZIP (phylobot.job_import), in the page's units.
IMPORT_DEFAULTS = {
    "import_max_megabytes": 2_000,
    "import_max_unpacked_megabytes": 3_000,
    "import_max_files": 50_000,
    "import_max_record_megabytes": 64,
    "import_max_libraries": 200,
    "import_max_exports": 500,
}


@dataclass(frozen=True)
class ImportLimits:
    """What phylobot.job_import accepts, in bytes and counts."""

    upload_bytes: int
    unpacked_bytes: int
    files: int
    record_bytes: int
    libraries: int
    exports: int


def _import(field: str) -> Callable[[], Any]:
    return lambda: IMPORT_DEFAULTS[field]


def _beta(key: str) -> Callable[[], Any]:
    return lambda: calibration.BUILT_IN.betas[key]


# Coefficients are bounded on the page; T itself is bounded again in phylobot.calibration.
CALIBRATION_BOUND = 5


def _calibration_settings() -> list[Setting]:
    intercept = Setting(
        "calibration_intercept",
        "Intercept (ln T at the reference site)",
        "The ancestor viewer's calibrated probabilities raise each state's probability to "
        "1/T and renormalise, with ln T = intercept + the sum of coefficient x covariate "
        "below; the ML state never changes. The reference site has every covariate at 0: "
        "branches of 0.1 substitutions per site, the mean rate, a conserved column, no "
        "gaps. Positive values flatten (less confident), negative sharpen. The built-in "
        "intercept is ln(1.24) ≈ 0.215, a global temperature fitted on a benchmark; the "
        "covariate terms are 0. These are the "
        "site's defaults: each user can set their own over them on My calibration.",
        CALIBRATION,
        lambda: calibration.BUILT_IN.intercept,
        low=-CALIBRATION_BOUND,
        high=CALIBRATION_BOUND,
        real=True,
    )
    return [intercept] + [
        Setting(
            f"calibration_{c.key}",
            f"Coefficient: {c.label.lower()}",
            f"{c.meaning} 0 at {c.reference}.",
            CALIBRATION,
            _beta(c.key),
            low=-CALIBRATION_BOUND,
            high=CALIBRATION_BOUND,
            real=True,
        )
        for c in calibration.COVARIATES
    ]


SETTINGS = [
    Setting(
        "search_email",
        "Contact email for searches",
        "EBI and NCBI require one with every search, and write to it if a search causes "
        "trouble: yours when PhyloBot runs on your machine, the operator's when hosted. "
        "Searches are refused while it is empty.",
        SEARCH,
        lambda: settings.PHYLOBOT_SEARCH_EMAIL,
        env="PHYLOBOT_SEARCH_EMAIL",
    ),
    Setting(
        "searches_per_day",
        "Searches per user per day",
        "New searches a user may start in any 24 hours.",
        SEARCH,
        lambda: settings.PHYLOBOT_SEARCHES_PER_DAY,
        low=0,
        high=10_000,
    ),
    Setting(
        "searches_in_flight_per_user",
        "Searches waiting at a service, per user",
        "A user's further searches wait in the queue. EBI asks that at most 30 wait at "
        "once from all of PhyloBot.",
        SEARCH,
        lambda: settings.PHYLOBOT_SEARCHES_IN_FLIGHT_PER_USER,
        low=1,
        high=30,
    ),
    Setting(
        "search_time_limit_minutes",
        "Abandon a search after",
        "Counted from when it was started, including time waiting in the queue.",
        SEARCH,
        lambda: round(settings.PHYLOBOT_SEARCH_TIME_LIMIT.total_seconds() / 60),
        unit="minutes",
        low=5,
        high=24 * 60,
    ),
    Setting(
        "max_sequences",
        "Sequences per job",
        "Default for every user; a user's own limits (Django admin, User limits) win.",
        USERS,
        _limit("max_sequences"),
        low=4,
        high=100_000,
    ),
    Setting(
        "max_length",
        "Residues per sequence",
        "Default for every user; a user's own limits win.",
        USERS,
        _limit("max_length"),
        low=10,
        high=100_000,
    ),
    Setting(
        "max_upload_megabytes",
        "Upload size per file",
        "Default for every user; a user's own limits win. At most 50 MB, the server's "
        "ceiling on any request.",
        USERS,
        lambda: round(settings.PHYLOBOT_DEFAULT_LIMITS["max_upload_bytes"] / MEGABYTE),
        unit="MB",
        low=1,
        high=50,
    ),
    Setting(
        "max_concurrent_jobs",
        "Jobs running at once, per user",
        "Default for every user; a user's own limits win. Further jobs wait in the queue.",
        USERS,
        _limit("max_concurrent_jobs"),
        low=1,
        high=100,
    ),
    Setting(
        "job_time_limit_hours",
        "Stop a job after",
        "A running job is stopped and failed after this long. Applies to jobs that start "
        "after the change.",
        JOBS,
        lambda: round(settings.PHYLOBOT_JOB_TIME_LIMIT.total_seconds() / 3600),
        unit="hours",
        low=1,
        high=24 * 30,
    ),
    Setting(
        "worker_concurrency",
        "Jobs the worker runs at once",
        "Across all users, per worker (one runs here). Each job is given an equal share "
        "of the CPUs and memory below, so more at once means less each. Takes effect "
        "within seconds; running jobs are not stopped if you lower it.",
        JOBS,
        lambda: settings.PHYLOBOT_WORKER_CONCURRENCY,
        env="PHYLOBOT_WORKER_CONCURRENCY",
        low=1,
        high=64,
    ),
    Setting(
        "machine_cpus",
        "CPUs for pipeline steps",
        "What the worker detects from Docker, less any isolated CPUs. Lower it to keep "
        "CPUs free for other work. Read when each job starts.",
        JOBS,
        _detected("cpus", "PHYLOBOT_CPUS"),
        env="PHYLOBOT_CPUS",
        low=1,
        high=1024,
        detected=True,
    ),
    Setting(
        "machine_memory_gb",
        "Memory of this machine",
        "What the worker detects from Docker (on a Mac, Docker's VM). Lower it if other "
        "programs need more than the reserve below. Read when each job starts.",
        JOBS,
        _detected("memory_bytes", "PHYLOBOT_MEMORY_GB", 1024**3),
        env="PHYLOBOT_MEMORY_GB",
        unit="GB",
        low=1,
        high=65_536,
        real=True,
        detected=True,
    ),
    Setting(
        "reserved_memory_gb",
        "Memory to leave for everything else",
        "Never given to pipeline steps: for the web portal, Postgres, Garage, and whatever "
        "else runs on this machine (a desktop takes 1-2 GB). The rest is shared between "
        "the jobs that can run at once. Read when each job starts.",
        JOBS,
        lambda: settings.PHYLOBOT_RESERVED_MEMORY_GB,
        env="PHYLOBOT_RESERVED_MEMORY_GB",
        unit="GB",
        low=0,
        high=65_536,
        real=True,
    ),
    Setting(
        "task_memory_gb",
        "Memory each pipeline step reserves",
        "Nextflow runs as many steps at once as a job's share holds at this much each, "
        "whatever they really use (tens of MB for most families; the trace's peak_rss "
        "shows it). Lower it to run more steps at once; a step that runs out is retried "
        "once with twice as much, within its job's share. Measuring properties (2 GB, or "
        "4 GB with predictors) and uncertainty runs (2 GB) keep their own. Read when each "
        "job starts.",
        JOBS,
        lambda: settings.PHYLOBOT_TASK_MEMORY_GB,
        env="PHYLOBOT_TASK_MEMORY_GB",
        unit="GB",
        low=1,
        high=64,
        real=True,
    ),
    Setting(
        "keep_failed_work_days",
        "Keep failed jobs' work folders for",
        "A failed job's folder on the worker holds every step it finished, which is what "
        "lets Resume skip them; for big jobs it can take many gigabytes. After this many "
        "days since the job failed the worker removes it, and a later resume starts from "
        "the beginning. 0 = keep until the job is deleted.",
        JOBS,
        lambda: settings.PHYLOBOT_KEEP_FAILED_WORK_DAYS,
        unit="days",
        low=0,
        high=3650,
        zero="until deleted",
    ),
    Setting(
        "archive_failed_work",
        "Archive failed jobs' work folders to object storage",
        "Copies each failed job's folder to the S3 bucket (Garage when local), so it can "
        "still be resumed if the worker's disk is lost, as when a cloud machine is "
        "replaced. Pointless on a laptop, where Garage shares that disk. The archive goes "
        "with the folder (the period above, or deleting the job) and once the job succeeds.",
        JOBS,
        lambda: settings.PHYLOBOT_ARCHIVE_FAILED_WORK,
        env="PHYLOBOT_ARCHIVE_FAILED_WORK",
        switch=True,
    ),
    Setting(
        "failed_task_rate_percent",
        "Bill failed steps at",
        "Percent of the normal rate for the CPU time of steps that failed or were stopped "
        "(0 = not billed, 100 = billed in full). Completed steps are always billed in "
        "full, and a resumed job's reused steps only once. Each attempt keeps the rate in "
        "force when it ended.",
        BILLING,
        lambda: settings.PHYLOBOT_FAILED_TASK_RATE_PERCENT,
        unit="%",
        low=0,
        high=100,
    ),
    Setting(
        "simulation_max_replicates",
        "Replicates per simulation",
        "Each replicate is one simulated alignment with its true tree and ancestors.",
        SIMULATIONS,
        lambda: SIMULATION_MAX_REPLICATES,
        low=1,
        high=1000,
    ),
    Setting(
        "simulation_max_megaresidues",
        "Residues written per simulation",
        "Replicates x (2 x sequences - 1) x sites, counting the true ancestors: what a "
        "simulation stores. AliSim itself is fast; storage is the cost.",
        SIMULATIONS,
        lambda: SIMULATION_MAX_MEGARESIDUES,
        unit="million",
        low=1,
        high=10_000,
    ),
    Setting(
        "import_max_megabytes",
        "Largest ZIP to import",
        "Import job on My Jobs (local mode only, for now). The ZIP as uploaded; the web "
        "server holds it in a temporary file while it is checked.",
        IMPORT,
        _import("import_max_megabytes"),
        unit="MB",
        low=1,
        high=4_000,  # see the next setting
    ),
    Setting(
        "import_max_unpacked_megabytes",
        "Its files, uncompressed",
        "Their total size, counted as they are read, so a ZIP that claims less (a "
        '"zip bomb") is stopped there. This is what the job takes in storage. At most '
        "4000 MB: the import runs inside the web request, which the server stops after "
        "30 seconds, and 1000 MB took 7 s on a laptop (results compress about 7x, so a "
        "small ZIP can hold much more).",
        IMPORT,
        _import("import_max_unpacked_megabytes"),
        unit="MB",
        low=1,
        high=4_000,
    ),
    Setting(
        "import_max_files",
        "Files in a ZIP",
        "Inputs, results and library exports together.",
        IMPORT,
        _import("import_max_files"),
        low=10,
        high=10_000_000,
    ),
    Setting(
        "import_max_record_megabytes",
        "Size of the job's record",
        "phylobot-job.json in the ZIP: the job's settings, progress, libraries and the "
        "list of its files, read into memory whole.",
        IMPORT,
        _import("import_max_record_megabytes"),
        unit="MB",
        low=1,
        high=1_000,
    ),
    Setting(
        "import_max_libraries",
        "Libraries per imported job",
        "",
        IMPORT,
        _import("import_max_libraries"),
        low=0,
        high=100_000,
    ),
    Setting(
        "import_max_exports",
        "Exports per library",
        "Each export is a receipt plus its files.",
        IMPORT,
        _import("import_max_exports"),
        low=0,
        high=100_000,
    ),
    *_calibration_settings(),
]
BY_FIELD = {s.field: s for s in SETTINGS}


def stored() -> SiteSettings | None:
    return SiteSettings.objects.filter(pk=1).first()


def value(field: str, row: SiteSettings | None = None) -> Any:
    """The value in effect, in the page's unit."""
    setting = BY_FIELD[field]
    chosen = getattr(row if row is not None else stored(), field, None)
    if chosen is None or chosen == "":
        return setting.fallback()
    return chosen


def source(field: str, row: SiteSettings | None) -> str:
    """Where the value in effect comes from, for the page."""
    chosen = getattr(row, field, None)
    if chosen is not None and chosen != "":
        return "set here"
    setting = BY_FIELD[field]
    if setting.env and os.environ.get(setting.env, "").strip():
        return f"from {setting.env}"
    if setting.detected:
        return "detected" if setting.fallback() is not None else "not detected yet"
    return "built-in default"


# Typed accessors for the code that uses them. Without a value set on the page, each
# returns Django's setting exactly (no unit conversion), so tests and .env behave as before.


def override(field: str, row: SiteSettings | None = None) -> Any:
    """The value set on the page, or None."""
    chosen = getattr(row if row is not None else stored(), field, None)
    return None if chosen == "" else chosen


def search_email() -> str:
    chosen = override("search_email")
    return str(chosen) if chosen is not None else str(settings.PHYLOBOT_SEARCH_EMAIL)


def searches_per_day() -> int:
    chosen = override("searches_per_day")
    return int(chosen) if chosen is not None else int(settings.PHYLOBOT_SEARCHES_PER_DAY)


def searches_in_flight_per_user() -> int:
    chosen = override("searches_in_flight_per_user")
    if chosen is not None:
        return int(chosen)
    return int(settings.PHYLOBOT_SEARCHES_IN_FLIGHT_PER_USER)


def search_time_limit() -> timedelta:
    chosen = override("search_time_limit_minutes")
    if chosen is not None:
        return timedelta(minutes=int(chosen))
    limit: timedelta = settings.PHYLOBOT_SEARCH_TIME_LIMIT
    return limit


def job_time_limit() -> timedelta:
    chosen = override("job_time_limit_hours")
    if chosen is not None:
        return timedelta(hours=int(chosen))
    limit: timedelta = settings.PHYLOBOT_JOB_TIME_LIMIT
    return limit


def worker_concurrency() -> int:
    chosen = override("worker_concurrency")
    return int(chosen) if chosen is not None else int(settings.PHYLOBOT_WORKER_CONCURRENCY)


def reserved_memory_gb() -> float:
    chosen = override("reserved_memory_gb")
    return float(chosen) if chosen is not None else float(settings.PHYLOBOT_RESERVED_MEMORY_GB)


def task_memory_gb() -> float:
    chosen = override("task_memory_gb")
    return float(chosen) if chosen is not None else float(settings.PHYLOBOT_TASK_MEMORY_GB)


def failed_task_rate_percent() -> int:
    chosen = override("failed_task_rate_percent")
    if chosen is not None:
        return int(chosen)
    return int(settings.PHYLOBOT_FAILED_TASK_RATE_PERCENT)


def keep_failed_work() -> timedelta | None:
    """How long a failed job's local folder is kept; None = until the job is deleted."""
    chosen = override("keep_failed_work_days")
    days = int(chosen) if chosen is not None else int(settings.PHYLOBOT_KEEP_FAILED_WORK_DAYS)
    return timedelta(days=days) if days else None


def archive_failed_work() -> bool:
    chosen = override("archive_failed_work")
    return bool(chosen) if chosen is not None else bool(settings.PHYLOBOT_ARCHIVE_FAILED_WORK)


def simulation_max_replicates() -> int:
    chosen = override("simulation_max_replicates")
    return int(chosen) if chosen is not None else SIMULATION_MAX_REPLICATES


def simulation_max_residues() -> int:
    chosen = override("simulation_max_megaresidues")
    return 1_000_000 * (int(chosen) if chosen is not None else SIMULATION_MAX_MEGARESIDUES)


def default_limits() -> dict[str, int]:
    """Every user's limits unless their own (UserLimits) say otherwise."""
    row = stored()
    limits = {k: int(v) for k, v in settings.PHYLOBOT_DEFAULT_LIMITS.items()}
    for field in ("max_sequences", "max_length", "max_concurrent_jobs"):
        if (chosen := override(field, row)) is not None:
            limits[field] = int(chosen)
    if (megabytes := override("max_upload_megabytes", row)) is not None:
        limits["max_upload_bytes"] = int(megabytes) * MEGABYTE
    return limits


# Calibration coefficients come in three layers, each field on its own: the user's own
# (UserCalibration, My calibration), else the site's defaults (the Settings page), else
# the built-in ones (calibration.BUILT_IN). Every field in the order of the formula:
CALIBRATION_TERMS = ("intercept", *calibration.KEYS)


def built_in_term(term: str) -> float:
    built_in = calibration.BUILT_IN
    return built_in.intercept if term == "intercept" else built_in.betas[term]


def site_calibration(row: SiteSettings | None = None) -> dict[str, float | None]:
    """The site's defaults as set on the Settings page (None = the built-in value)."""
    row = row if row is not None else stored()
    return {t: override(f"calibration_{t}", row) for t in CALIBRATION_TERMS}


def own_calibration(user: Any) -> UserCalibration | None:
    if user is None or not getattr(user, "is_authenticated", False):
        return None
    return UserCalibration.objects.filter(user_id=user.pk).first()


def calibration_coefficients(user: Any = None) -> calibration.Coefficients:
    """The coefficients in force for `user` (None: the site's). The source says which
    layer they come from: anything but the built-in ones is no longer what was fitted."""
    row = stored()
    site = site_calibration(row)
    own = own_calibration(user)
    mine = {t: getattr(own, t) if own is not None else None for t in CALIBRATION_TERMS}
    values = {}
    for t in CALIBRATION_TERMS:
        chosen = mine[t] if mine[t] is not None else site[t]
        values[t] = float(chosen) if chosen is not None else built_in_term(t)
    if any(v is not None for v in mine.values()):
        assert own is not None
        source = f"your own (My calibration, {own.updated_at:%Y-%m-%d})"
    elif any(v is not None for v in site.values()):
        when = f", {row.updated_at:%Y-%m-%d}" if row is not None and row.updated_at else ""
        source = f"this site's defaults (set on the Settings page{when})"
    else:
        return calibration.BUILT_IN
    return calibration.Coefficients(
        intercept=values["intercept"],
        betas={k: values[k] for k in calibration.KEYS},
        source=source,
    )


def import_limits() -> ImportLimits:
    row = stored()

    def chosen(field: str) -> int:
        set_here = override(field, row)
        return int(set_here) if set_here is not None else IMPORT_DEFAULTS[field]

    return ImportLimits(
        upload_bytes=chosen("import_max_megabytes") * MEGABYTE,
        unpacked_bytes=chosen("import_max_unpacked_megabytes") * MEGABYTE,
        files=chosen("import_max_files"),
        record_bytes=chosen("import_max_record_megabytes") * MEGABYTE,
        libraries=chosen("import_max_libraries"),
        exports=chosen("import_max_exports"),
    )
