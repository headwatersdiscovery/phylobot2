"""Follow homolog searches at EBI and NCBI: python manage.py run_searches

One loop moves every due search one step on (phylobot.searches): submit, check, fetch.
Searches mostly wait on the service, so no thread waits on each; one runner follows
any number. Several runners may run at once (claims are locked). Needs no Docker and no
pipeline: only the database, object storage, and outbound HTTPS.
"""

import logging
import signal
import time
from types import FrameType
from typing import Any

from django.core.management.base import BaseCommand
from django.db import close_old_connections

from phylobot import searches

logger = logging.getLogger(__name__)

IDLE_SECONDS = 2
ERROR_SECONDS = 15


class Command(BaseCommand):
    help = "Follow homolog searches at the remote services."

    def handle(self, *args: Any, **options: Any) -> None:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
        stopping = False

        def stop(signum: int, frame: FrameType | None) -> None:
            nonlocal stopping
            logger.info("stopping after the current step")
            stopping = True

        signal.signal(signal.SIGTERM, stop)
        signal.signal(signal.SIGINT, stop)
        logger.info("following homolog searches")
        while not stopping:
            close_old_connections()
            try:
                due = searches.claim_due()
            except Exception:
                # The database is unreachable or not yet migrated (e.g. started before
                # the web container applied a migration): wait and try again, never exit.
                logger.exception("could not claim searches; retrying in %d s", ERROR_SECONDS)
                time.sleep(ERROR_SECONDS)
                continue
            for search in due:
                if stopping:
                    break  # the rest are picked up again when their lease runs out
                try:
                    searches.advance(search)
                except Exception:
                    # A bug, not service trouble: log it and let the lease expire, so
                    # one bad search cannot stop the others.
                    logger.exception("search %s: unexpected error", search.id)
            if not due:
                time.sleep(IDLE_SECONDS)
