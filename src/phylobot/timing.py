"""Request timing: how long each request took and how much of it was the database.

Every response gets a `Server-Timing` header in local mode (the browser's developer tools
show it under Network -> Timing); requests slower than SLOW_SECONDS are logged in both
modes, with their database share, so a slow page can be traced from `docker compose logs
web` without reproducing it.
"""

import logging
import time
from collections.abc import Callable
from typing import Any

from django.conf import settings
from django.db import connection
from django.http import HttpRequest, HttpResponse

logger = logging.getLogger(__name__)
SLOW_SECONDS = 0.3


class TimingMiddleware:
    def __init__(self, get_response: Callable[[HttpRequest], HttpResponse]) -> None:
        self.get_response = get_response

    def __call__(self, request: HttpRequest) -> HttpResponse:
        queries = [0, 0.0]  # count, seconds

        def count(execute: Callable[..., Any], *args: Any) -> Any:
            start = time.perf_counter()
            try:
                return execute(*args)
            finally:
                queries[0] += 1
                queries[1] += time.perf_counter() - start

        start = time.perf_counter()
        with connection.execute_wrapper(count):
            response = self.get_response(request)
        total = time.perf_counter() - start
        if settings.PHYLOBOT_MODE == "local":
            db = f'db;dur={queries[1] * 1000:.1f};desc="{queries[0]} queries"'
            response["Server-Timing"] = f"total;dur={total * 1000:.1f}, {db}"
        if total > SLOW_SECONDS:
            logger.warning(
                "slow request: %s %s %.0f ms (database %.0f ms, %d queries)",
                request.method, request.get_full_path()[:200], total * 1000,
                queries[1] * 1000, queries[0],
            )  # fmt: skip
        return response
