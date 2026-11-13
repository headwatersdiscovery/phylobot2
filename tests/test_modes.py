import os
import subprocess
import sys
from unittest import mock

import pytest
from django.contrib.auth import get_user_model
from django.test import Client, override_settings

from phylobot.modes import LOCAL_USERNAME, local_user
from phylobot.settings import middleware_for

CLOUD = {"PHYLOBOT_MODE": "cloud", "MIDDLEWARE": middleware_for("cloud")}
REAL_KEY = "k" * 60
# Variables the settings check; the test container's own values are dropped.
SETTINGS_ENV = ("PHYLOBOT_", "DJANGO_", "AWS_")


def settings_error(**env: str) -> str:
    """Import settings in a fresh interpreter with `env`; return the error ('' if none)."""
    base = {k: v for k, v in os.environ.items() if not k.startswith(SETTINGS_ENV)}
    proc = subprocess.run(
        [sys.executable, "-c", "import phylobot.settings"],
        env={**base, "DJANGO_SECRET_KEY": "local-dev-only", **env},
        capture_output=True,
        text=True,
    )
    return proc.stderr.strip().splitlines()[-1] if proc.returncode else ""


def test_mode_is_required() -> None:
    assert "PHYLOBOT_MODE must be set" in settings_error()
    assert "PHYLOBOT_MODE must be set" in settings_error(PHYLOBOT_MODE="production")


def test_local_mode_refuses_public_hosts() -> None:
    assert settings_error(PHYLOBOT_MODE="local") == ""
    assert "must not be reachable" in settings_error(
        PHYLOBOT_MODE="local", DJANGO_ALLOWED_HOSTS="localhost,phylobot.example.org"
    )
    assert "must not be reachable" in settings_error(
        PHYLOBOT_MODE="local", DJANGO_ALLOWED_HOSTS="*"
    )


# The shape of a Garage key (GK + 24 hex), generated per install by scripts/env-secrets.sh.
GARAGE_KEY = "GK" + "0" * 24


def test_cloud_mode_refuses_dev_settings() -> None:
    assert "real DJANGO_SECRET_KEY" in settings_error(PHYLOBOT_MODE="cloud")
    assert "DJANGO_DEBUG" in settings_error(
        PHYLOBOT_MODE="cloud", DJANGO_SECRET_KEY=REAL_KEY, DJANGO_DEBUG="true"
    )
    assert "local Garage key" in settings_error(
        PHYLOBOT_MODE="cloud", DJANGO_SECRET_KEY=REAL_KEY, AWS_ACCESS_KEY_ID=GARAGE_KEY
    )
    assert settings_error(PHYLOBOT_MODE="cloud", DJANGO_SECRET_KEY=REAL_KEY) == ""
    assert settings_error(PHYLOBOT_MODE="local", AWS_ACCESS_KEY_ID=GARAGE_KEY) == ""


@pytest.mark.django_db
def test_local_mode_acts_as_local_user(client: Client) -> None:
    assert client.get("/").status_code == 200
    client.get("/")
    users = get_user_model().objects.filter(username=LOCAL_USERNAME)
    assert users.count() == 1
    assert not users.get().has_usable_password()
    assert client.get("/admin/").status_code == 200  # local user can use the admin


@pytest.mark.django_db
def test_local_user_is_stable() -> None:
    assert local_user().pk == local_user().pk


@pytest.mark.django_db
@override_settings(**CLOUD)
def test_cloud_mode_requires_login(client: Client) -> None:
    response = client.get("/")
    assert response.status_code == 302
    assert response["Location"].startswith("/accounts/login/")
    assert client.get("/accounts/login/").status_code == 200
    assert not get_user_model().objects.filter(username=LOCAL_USERNAME).exists()

    get_user_model().objects.create_user("alice", password="pw-for-tests-only")
    assert client.login(username="alice", password="pw-for-tests-only")
    assert client.get("/").status_code == 200


@pytest.mark.django_db
@override_settings(**CLOUD)
def test_healthz_needs_no_login(client: Client) -> None:
    with mock.patch("phylobot.storage.bucket_reachable", return_value=True):
        response = client.get("/healthz")
    assert response.status_code == 200
    assert response.json() == {"postgres": "ok", "storage": "ok"}
    with mock.patch("phylobot.storage.bucket_reachable", return_value=False):
        assert client.get("/healthz").status_code == 503


def test_template_caching_by_mode() -> None:
    """Local mode picks up template edits live; cloud mode caches templates."""
    code = "import phylobot.settings as s; print(s.TEMPLATES[0]['OPTIONS']['loaders'])"
    base = {k: v for k, v in os.environ.items() if not k.startswith(SETTINGS_ENV)}
    for mode, key, cached in (("local", "local-dev-only", False), ("cloud", REAL_KEY, True)):
        env = {**base, "PHYLOBOT_MODE": mode, "DJANGO_SECRET_KEY": key}
        out = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True)
        assert ("cached.Loader" in out.stdout) is cached, out.stderr
