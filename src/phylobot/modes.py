"""Deployment-mode behaviour. See PHYLOBOT_MODE in settings."""

from collections.abc import Callable

from django.contrib.auth import get_user_model
from django.contrib.auth.models import User
from django.http import HttpRequest, HttpResponse

LOCAL_USERNAME = "local"


def local_user() -> User:
    """The single user every request acts as in local mode (created on first use)."""
    user: User
    user, created = get_user_model().objects.get_or_create(
        username=LOCAL_USERNAME,
        defaults={"is_staff": True, "is_superuser": True},  # the admin works locally
    )
    if created:
        user.set_unusable_password()
        user.save(update_fields=["password"])
    return user


class LocalUserMiddleware:
    """Local mode only: treat every request as the local user. Never used in cloud mode;
    settings refuse local mode unless the site is reachable only from this machine."""

    def __init__(self, get_response: Callable[[HttpRequest], HttpResponse]) -> None:
        self.get_response = get_response

    def __call__(self, request: HttpRequest) -> HttpResponse:
        request.user = local_user()
        return self.get_response(request)
