"""Site settings: the page's values override the environment and defaults, blank falls
back, the code reads them live, and only the right people see the page."""

from datetime import timedelta
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.test import Client, override_settings

from phylobot import calibration, searches, site_settings
from phylobot.models import SiteSettings, UserCalibration, UserLimits, limits_for
from phylobot.modes import local_user
from phylobot.settings import middleware_for

pytestmark = pytest.mark.django_db


def save(**fields: Any) -> None:
    row = SiteSettings.objects.filter(pk=1).first() or SiteSettings(pk=1)
    for name, value in fields.items():
        setattr(row, name, value)
    row.save()


def test_without_the_page_everything_comes_from_django_settings(settings: Any) -> None:
    settings.PHYLOBOT_SEARCH_EMAIL = "env@example.org"
    settings.PHYLOBOT_JOB_TIME_LIMIT = timedelta(minutes=30)  # kept exactly, not rounded
    assert site_settings.search_email() == "env@example.org"
    assert site_settings.job_time_limit() == timedelta(minutes=30)
    assert site_settings.default_limits() == settings.PHYLOBOT_DEFAULT_LIMITS


def test_the_page_overrides_and_blank_falls_back(settings: Any) -> None:
    settings.PHYLOBOT_SEARCH_EMAIL = "env@example.org"
    save(search_email="page@example.org", job_time_limit_hours=3, worker_concurrency=5)
    assert site_settings.search_email() == "page@example.org"
    assert site_settings.job_time_limit() == timedelta(hours=3)
    assert site_settings.worker_concurrency() == 5
    save(search_email="", job_time_limit_hours=None)
    assert site_settings.search_email() == "env@example.org"
    assert site_settings.job_time_limit() == settings.PHYLOBOT_JOB_TIME_LIMIT


def test_default_limits_then_each_users_own(settings: Any) -> None:
    save(max_sequences=40, max_upload_megabytes=2)
    user = User.objects.create_user("someone")
    limits = limits_for(user)
    assert (limits.max_sequences, limits.max_upload_bytes) == (40, 2_000_000)
    assert limits.max_length == settings.PHYLOBOT_DEFAULT_LIMITS["max_length"]
    UserLimits.objects.create(user=user, max_sequences=7)
    assert limits_for(user).max_sequences == 7  # the user's own limit wins


def test_search_email_set_on_the_page_enables_searching(settings: Any) -> None:
    settings.PHYLOBOT_SEARCH_EMAIL = ""
    with pytest.raises(searches.LimitError, match="Settings page"):
        searches.check_can_search(local_user())
    save(search_email="page@example.org")
    searches.check_can_search(local_user())  # no restart needed


def test_sources_shown_on_the_page(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PHYLOBOT_SEARCH_EMAIL", "env@example.org")
    assert site_settings.source("search_email", None) == "from PHYLOBOT_SEARCH_EMAIL"
    assert site_settings.source("searches_per_day", None) == "built-in default"
    save(searches_per_day=3)
    assert site_settings.source("searches_per_day", site_settings.stored()) == "set here"


# --- the page ----------------------------------------------------------------------


def test_local_mode_settings_page(client: Client) -> None:
    home = client.get("/").content.decode()
    assert ">Settings</a>" in home and "Admin Settings" not in home
    page = client.get("/settings/").content.decode()
    assert "Contact email for searches" in page and "built-in default" in page
    # Local mode has no accounts: Billing is hidden and the per-user framing is dropped.
    assert ">Billing</h2>" not in page
    assert "Limits Across Jobs" in page and "Limits for Every User" not in page
    assert "Searches per day" in page and "Searches per user per day" not in page
    assert "Searches waiting at a service" in page
    assert "Searches waiting at a service, per user" not in page
    # Save lives in the sticky "On this page" sidebar and still submits the settings form.
    assert 'form="settings-form"' in page and "settings-save" in page

    response = client.post(
        "/settings/", {"search_email": "me@example.org", "searches_per_day": "5"}
    )
    assert response.status_code == 302
    row = SiteSettings.objects.get()
    assert (row.search_email, row.searches_per_day, row.worker_concurrency) == (
        "me@example.org",
        5,
        None,
    )
    assert row.updated_by == local_user()
    page = client.get("/settings/").content.decode()
    assert "set here" in page and "me@example.org" in page

    client.post("/settings/", {"search_email": "", "searches_per_day": ""})
    row = SiteSettings.objects.get()
    assert (row.search_email, row.searches_per_day) == ("", None)


@pytest.mark.parametrize(
    "fields",
    [
        {"searches_in_flight_per_user": "31"},  # EBI: at most 30 at once
        {"max_upload_megabytes": "51"},  # the server's request ceiling
        {"worker_concurrency": "0"},
        {"search_email": "not an email"},
    ],
)
def test_settings_page_refuses_out_of_bounds(client: Client, fields: dict[str, str]) -> None:
    response = client.post("/settings/", fields)
    assert response.status_code == 200 and "errors" in response.content.decode()
    assert not SiteSettings.objects.exists()


@override_settings(PHYLOBOT_MODE="cloud", MIDDLEWARE=middleware_for("cloud"))
def test_cloud_mode_settings_are_for_staff_only(client: Client) -> None:
    assert client.get("/settings/").status_code == 302  # login first
    client.force_login(User.objects.create_user("someone"))
    assert client.get("/settings/").status_code == 404
    assert client.post("/settings/", {"searches_per_day": "1"}).status_code == 404
    assert "settings" not in client.get("/jobs/").content.decode().split("</nav>")[0].lower()
    client.force_login(User.objects.create_user("admin", is_staff=True))
    assert client.get("/settings/").status_code == 200
    assert "Admin Settings" in client.get("/jobs/").content.decode()


def test_failed_work_retention_page(client: Client) -> None:
    import re

    from phylobot.models import Job

    page = re.sub(r"\s+", " ", client.get("/settings/").content.decode())
    assert "Keep failed jobs&#x27; work folders for" in page
    assert "<strong>until deleted</strong>" in page and "Kept now: none" in page
    owner = local_user()
    Job.objects.create(owner=owner, status=Job.Status.FAILED, work_bytes=3 * 1024**3)
    Job.objects.create(owner=owner, status=Job.Status.FAILED, work_bytes=1024**3)
    Job.objects.create(owner=owner, status=Job.Status.SUCCEEDED, work_bytes=10**12)  # not kept
    assert site_settings.keep_failed_work() is None
    client.post("/settings/", {"keep_failed_work_days": "14"})
    assert site_settings.keep_failed_work() == timedelta(days=14)
    page = re.sub(r"\s+", " ", client.get("/settings/").content.decode())
    assert "<strong>14</strong> days" in page and "Kept now: 4.0 GB in 2 failed jobs" in page
    client.post("/settings/", {"keep_failed_work_days": "0"})
    assert site_settings.keep_failed_work() is None  # 0 = until the job is deleted
    # A zero is shown, not "(empty)". (failed_task_rate_percent is a Billing setting, which
    # local mode hides, so this uses searches_per_day instead.)
    save(searches_per_day=0)
    page = re.sub(r"\s+", " ", client.get("/settings/").content.decode())
    assert "<strong>0</strong>" in page


def test_archive_switch(client: Client) -> None:
    import re

    assert site_settings.archive_failed_work() is False
    page = re.sub(r"\s+", " ", client.get("/settings/").content.decode())
    assert "<strong>off</strong>" in page and "Archived now: none" in page
    client.post("/settings/", {"archive_failed_work": "1"})
    assert SiteSettings.objects.get().archive_failed_work is True
    assert site_settings.archive_failed_work() is True
    page = client.get("/settings/").content.decode()
    assert '<option value="1" selected>On</option>' in page
    client.post("/settings/", {"archive_failed_work": "0"})
    assert site_settings.archive_failed_work() is False  # off on the page beats .env
    client.post("/settings/", {"archive_failed_work": ""})
    assert SiteSettings.objects.get().archive_failed_work is None


# --- calibration coefficients ------------------------------------------------------


def test_calibration_coefficients_built_in_until_set_on_the_page() -> None:
    assert site_settings.calibration_coefficients() == calibration.BUILT_IN
    save(calibration_intercept=0.4, calibration_branch_length=-0.25)
    chosen = site_settings.calibration_coefficients()
    assert chosen.intercept == 0.4 and chosen.betas["branch_length"] == -0.25
    assert chosen.betas["entropy"] == calibration.BUILT_IN.betas["entropy"]
    assert chosen.source.startswith("this site's defaults") and not chosen.fitted
    save(calibration_intercept=None, calibration_branch_length=None)
    assert site_settings.calibration_coefficients() == calibration.BUILT_IN


def test_settings_page_takes_decimal_coefficients(client: Client) -> None:
    page = client.get("/settings/").content.decode()
    assert "Calibrated Posterior Probabilities" in page and "Coefficient: site rate" in page
    chosen = {"calibration_intercept": "0.405", "calibration_entropy": "-1.5"}
    response = client.post("/settings/", chosen)
    assert response.status_code == 302
    row = SiteSettings.objects.get()
    assert (row.calibration_intercept, row.calibration_entropy) == (0.405, -1.5)
    assert row.calibration_site_rate is None
    response = client.post("/settings/", {"calibration_intercept": "6"})  # out of bounds
    assert response.status_code == 200 and "errors" in response.content.decode()


def test_a_users_own_coefficients_beat_the_sites_field_by_field() -> None:
    save(calibration_intercept=0.4, calibration_entropy=1.0)
    user = User.objects.create_user("someone")
    assert site_settings.calibration_coefficients(user).intercept == 0.4  # nothing of theirs
    UserCalibration.objects.create(user=user, intercept=-0.2)
    mine = site_settings.calibration_coefficients(user)
    assert mine.intercept == -0.2 and mine.betas["entropy"] == 1.0  # the rest the site's
    assert mine.source.startswith("your own")
    assert site_settings.calibration_coefficients().intercept == 0.4  # nobody else's
    other = User.objects.create_user("other")
    assert site_settings.calibration_coefficients(other).intercept == 0.4


def test_settings_page_resets_only_the_calibration_defaults(client: Client) -> None:
    save(calibration_intercept=0.4, calibration_entropy=1.0, searches_per_day=3)
    page = client.get("/settings/").content.decode()
    assert "Reset to built-in defaults" in page and "No effect (T = 1)" in page
    assert client.post("/settings/calibration-defaults/", {"to": "none"}).status_code == 302
    row = SiteSettings.objects.get()
    assert {getattr(row, f"calibration_{t}") for t in site_settings.CALIBRATION_TERMS} == {0.0}
    assert row.searches_per_day == 3  # untouched
    coefficients = site_settings.calibration_coefficients()
    assert coefficients.identity and coefficients.source.startswith("this site's defaults")
    client.post("/settings/calibration-defaults/", {"to": "default"})
    assert site_settings.calibration_coefficients() == calibration.BUILT_IN
    assert client.post("/settings/calibration-defaults/", {"to": "x"}).status_code == 404


def test_my_calibration_save_reset_and_no_effect(client: Client) -> None:
    save(calibration_intercept=0.4)
    page = client.get("/calibration/").content.decode()
    assert "My Calibration" in page and "Reset to site defaults" in page
    assert client.post("/calibration/", {"action": "save", "intercept": "0.7"}).status_code == 302
    own = UserCalibration.objects.get(user=local_user())
    assert own.intercept == 0.7 and own.entropy is None
    assert site_settings.calibration_coefficients(local_user()).intercept == 0.7
    client.post("/calibration/", {"action": "none"})
    assert site_settings.calibration_coefficients(local_user()).identity
    client.post("/calibration/", {"action": "default"})
    assert site_settings.calibration_coefficients(local_user()).intercept == 0.4  # the site's
    response = client.post("/calibration/", {"action": "save", "intercept": "9"})
    assert response.status_code == 200 and "errors" in response.content.decode()


@override_settings(PHYLOBOT_MODE="cloud", MIDDLEWARE=middleware_for("cloud"))
def test_cloud_users_set_their_own_but_not_the_sites(client: Client) -> None:
    client.force_login(User.objects.create_user("someone"))
    assert client.get("/calibration/").status_code == 200
    client.post("/calibration/", {"action": "save", "entropy": "0.5"})
    assert UserCalibration.objects.get().entropy == 0.5
    assert client.post("/settings/calibration-defaults/", {"to": "none"}).status_code == 404
    assert not SiteSettings.objects.exists()
