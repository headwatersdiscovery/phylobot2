"""Django settings. Every value comes from the environment; see .env.example."""

import os
import re
from datetime import timedelta
from pathlib import Path

from django.core.exceptions import ImproperlyConfigured

BASE_DIR = Path(__file__).resolve().parent


def env_bool(name: str, default: bool = False) -> bool:
    return os.environ.get(name, str(default)).strip().lower() in {"1", "true", "yes", "on"}


def env_list(name: str, default: str = "") -> list[str]:
    return [item.strip() for item in os.environ.get(name, default).split(",") if item.strip()]


SECRET_KEY = os.environ["DJANGO_SECRET_KEY"]
DEBUG = env_bool("DJANGO_DEBUG")
ALLOWED_HOSTS = env_list("DJANGO_ALLOWED_HOSTS", "localhost,127.0.0.1")

# Deployment mode, required. "local": no accounts, every request is the built-in local
# user; only reachable from this machine. "cloud": real accounts, login required.
PHYLOBOT_MODE = os.environ.get("PHYLOBOT_MODE", "")
LOCAL_HOSTS = {"localhost", "127.0.0.1", "::1", "[::1]"}
DEV_SECRET_MARKERS = ("local-dev-only", "not-secret", "change-me")
# Garage access key IDs start "GK" (AWS ones "AKIA"/"ASIA"): such a key in cloud mode is a
# local .env carried over, and the hosted site's storage is S3.
GARAGE_KEY_PREFIX = "GK"

if PHYLOBOT_MODE not in {"local", "cloud"}:
    raise ImproperlyConfigured("PHYLOBOT_MODE must be set to 'local' or 'cloud'")
if PHYLOBOT_MODE == "local" and not set(ALLOWED_HOSTS) <= LOCAL_HOSTS:
    raise ImproperlyConfigured(
        "PHYLOBOT_MODE=local serves every request as the local user, so it must not be "
        f"reachable from other machines; DJANGO_ALLOWED_HOSTS may only list {sorted(LOCAL_HOSTS)}"
    )
if PHYLOBOT_MODE == "cloud":
    if DEBUG:
        raise ImproperlyConfigured("PHYLOBOT_MODE=cloud refuses to run with DJANGO_DEBUG on")
    if len(SECRET_KEY) < 50 or any(marker in SECRET_KEY for marker in DEV_SECRET_MARKERS):
        raise ImproperlyConfigured("PHYLOBOT_MODE=cloud needs a real DJANGO_SECRET_KEY")
    if os.environ.get("AWS_ACCESS_KEY_ID", "").startswith(GARAGE_KEY_PREFIX):
        raise ImproperlyConfigured(
            "PHYLOBOT_MODE=cloud refuses a local Garage key; "
            "on AWS leave AWS_ACCESS_KEY_ID unset (IAM role), or set your own key"
        )

INSTALLED_APPS = [
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    "phylobot",
]


def middleware_for(mode: str) -> list[str]:
    """Local mode signs every request in as the local user; cloud mode requires login."""
    identity = {
        "local": "phylobot.modes.LocalUserMiddleware",
        "cloud": "django.contrib.auth.middleware.LoginRequiredMiddleware",
    }[mode]
    return [
        "phylobot.timing.TimingMiddleware",
        "django.middleware.security.SecurityMiddleware",
        "whitenoise.middleware.WhiteNoiseMiddleware",
        "django.contrib.sessions.middleware.SessionMiddleware",
        "django.middleware.common.CommonMiddleware",
        "django.middleware.csrf.CsrfViewMiddleware",
        "django.contrib.auth.middleware.AuthenticationMiddleware",
        identity,
        "django.contrib.messages.middleware.MessageMiddleware",
        "django.middleware.clickjacking.XFrameOptionsMiddleware",
    ]


MIDDLEWARE = middleware_for(PHYLOBOT_MODE)
LOGIN_URL = "login"
LOGIN_REDIRECT_URL = "home"
LOGOUT_REDIRECT_URL = "login"

ROOT_URLCONF = "phylobot.urls"
WSGI_APPLICATION = "phylobot.wsgi.application"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [],
        "OPTIONS": {
            # Local mode reads templates fresh on every request, so edits show at once
            # (gunicorn only reloads on Python changes). Cloud mode caches them.
            "loaders": (
                ["django.template.loaders.app_directories.Loader"]
                if PHYLOBOT_MODE == "local"
                else [
                    (
                        "django.template.loaders.cached.Loader",
                        ["django.template.loaders.app_directories.Loader"],
                    )
                ]
            ),
            "context_processors": [
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
                "phylobot.skins.context",
                "phylobot.library_views.context",
            ],
        },
    },
]

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.postgresql",
        "NAME": os.environ.get("POSTGRES_DB", "phylobot"),
        "USER": os.environ.get("POSTGRES_USER", "phylobot"),
        "PASSWORD": os.environ.get("POSTGRES_PASSWORD", ""),
        "HOST": os.environ.get("POSTGRES_HOST", "localhost"),
        "PORT": os.environ.get("POSTGRES_PORT", "5432"),
        "CONN_MAX_AGE": 60,
        "CONN_HEALTH_CHECKS": True,
    }
}

# Object storage (Garage locally, S3 in production). Credentials come from the standard
# AWS_* environment variables that boto3 reads.
S3_ENDPOINT_URL = os.environ.get("S3_ENDPOINT_URL") or None
S3_BUCKET = os.environ.get("S3_BUCKET", "phylobot")

# Site-wide default limits; per-user overrides live in UserLimits (Django admin).
PHYLOBOT_DEFAULT_LIMITS = {
    "max_sequences": 500,
    "max_length": 5000,
    "max_upload_bytes": 10_000_000,
    # One job already runs its steps side by side; more at once mostly compete (owner
    # decision 2026-10-26, was 3).
    "max_concurrent_jobs": 1,
}
# Worker (manage.py run_worker). PROJECT_ROOT: this checkout's absolute path, mounted at
# the same path in the worker, the Nextflow runner and task containers (set by make).
PHYLOBOT_PROJECT_ROOT = os.environ.get("PROJECT_ROOT", "")
# Keep in step with the runner service image in docker-compose.yml (a test checks).
PHYLOBOT_NEXTFLOW_IMAGE = "nextflow/nextflow:26.04.6"
PHYLOBOT_WORKER_CONCURRENCY = int(os.environ.get("PHYLOBOT_WORKER_CONCURRENCY", "2"))
# CPUs every pipeline task container may use (docker --cpuset-cpus), e.g. "0-6"; empty
# means all. OrbStack isolates its VM's last CPU, and emulated IQ-TREE gets stuck there.
# `make docker-cpus` suggests a value for this machine; see .env.example.
PHYLOBOT_DOCKER_CPUSET = os.environ.get("PHYLOBOT_DOCKER_CPUSET", "").strip()
if PHYLOBOT_DOCKER_CPUSET and not re.fullmatch(r"\d+(-\d+)?(,\d+(-\d+)?)*", PHYLOBOT_DOCKER_CPUSET):
    raise ImproperlyConfigured(
        f"PHYLOBOT_DOCKER_CPUSET must look like 0-6 or 0,2-5, not {PHYLOBOT_DOCKER_CPUSET!r}"
    )
# The machine as PhyloBot treats it (phylobot.machine): unset (0), what the worker detects
# from Docker. Set either to use less, or when detection is wrong; the Settings page wins
# over both. The reserve is memory left for everything but pipeline steps.
PHYLOBOT_CPUS = int(os.environ.get("PHYLOBOT_CPUS") or 0)
PHYLOBOT_MEMORY_GB = float(os.environ.get("PHYLOBOT_MEMORY_GB") or 0)
PHYLOBOT_RESERVED_MEMORY_GB = float(os.environ.get("PHYLOBOT_RESERVED_MEMORY_GB") or 2)
# What each pipeline step reserves (phylobot.machine): it sets how many steps run at once.
# The default, 2, is nextflow.config's own (pinned by a test).
PHYLOBOT_TASK_MEMORY_GB = float(os.environ.get("PHYLOBOT_TASK_MEMORY_GB") or 2)
PHYLOBOT_JOB_TIME_LIMIT = timedelta(hours=24)  # owner decision: a compute guardrail
# Failed steps are billed at this percent of the normal rate (0 = not billed); the admin
# sets it on the Settings page. 100 keeps what cpu_seconds always counted.
PHYLOBOT_FAILED_TASK_RATE_PERCENT = 100
# Days the worker keeps a failed job's local folder (its resume cache); 0 = until the job
# is deleted. Set on the Settings page.
PHYLOBOT_KEEP_FAILED_WORK_DAYS = 0
# Archive failed jobs' work folders to object storage (phylobot.work_archive), so a resume
# survives losing the worker's disk. Off by default: locally Garage shares that disk.
PHYLOBOT_ARCHIVE_FAILED_WORK = os.environ.get("PHYLOBOT_ARCHIVE_FAILED_WORK", "") == "1"
PHYLOBOT_HEARTBEAT_SECONDS = 10
PHYLOBOT_STALE_AFTER = timedelta(minutes=10)  # well above the heartbeat interval

# Homolog searches (phylobot.searches, run_searches). EBI and NCBI require a contact
# email with every search: the person running this PhyloBot (see .env.example).
PHYLOBOT_SEARCH_EMAIL = os.environ.get("PHYLOBOT_SEARCH_EMAIL", "").strip()
# Placeholders pending owner decision: searches per user per day, a user's
# searches waiting at a service at once, and how long a search may wait in all.
PHYLOBOT_SEARCHES_PER_DAY = 25
PHYLOBOT_SEARCHES_IN_FLIGHT_PER_USER = 3
PHYLOBOT_SEARCH_TIME_LIMIT = timedelta(hours=2)

# Hard ceiling on any request body, above every per-user upload limit.
DATA_UPLOAD_MAX_MEMORY_SIZE = 50_000_000
FILE_UPLOAD_MAX_MEMORY_SIZE = 5_000_000

AUTH_PASSWORD_VALIDATORS = [
    {"NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator"},
    {"NAME": "django.contrib.auth.password_validation.MinimumLengthValidator"},
    {"NAME": "django.contrib.auth.password_validation.CommonPasswordValidator"},
    {"NAME": "django.contrib.auth.password_validation.NumericPasswordValidator"},
]

LANGUAGE_CODE = "en-us"
TIME_ZONE = "UTC"
USE_I18N = True
USE_TZ = True

STATIC_URL = "static/"
WHITENOISE_USE_FINDERS = True  # serve app static files without a collectstatic step
STATIC_ROOT = None  # set when deployment adds a collectstatic step

DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"
