"""Run queued jobs: python manage.py run_worker

Claims jobs from the Job table (phylobot.jobs) and runs each in a thread with the
LocalDockerExecutor, up to site_settings.worker_concurrency() at once (read every loop,
so the Settings page changes it live). On SIGTERM/SIGINT it stops claiming, stops
running pipelines, and marks their jobs failed.
"""

import logging
import os
import shutil
import signal
import socket
import subprocess
import threading
import time
import uuid
from pathlib import Path
from types import FrameType
from typing import Any

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import close_old_connections

from phylobot import jobs, machine, manifest, site_settings, work_archive
from phylobot.executors import (
    LocalDockerExecutor,
    expire_work_folders,
    measure_work_folders,
    record_unrecorded,
    run_job,
)
from phylobot.models import Job

logger = logging.getLogger("phylobot.worker")


def is_uuid(text: str) -> bool:
    try:
        uuid.UUID(text)
    except ValueError:
        return False
    return True


POLL_SECONDS = 5


class Command(BaseCommand):
    help = "Run queued PhyloBot jobs."

    def handle(self, *args: Any, **options: Any) -> None:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
        root = Path(settings.PHYLOBOT_PROJECT_ROOT)
        if not (root / "pipeline" / "main.nf").is_file():
            raise CommandError(
                f"PROJECT_ROOT ({root or 'unset'}) must be this checkout, mounted at the "
                "same path as on the host; start the worker with `make up`."
            )
        name = f"{socket.gethostname()}-{os.getpid()}"
        executor = LocalDockerExecutor()
        stop = threading.Event()
        running: dict[str, threading.Thread] = {}
        containers: dict[str, str] = {}  # job id -> its current attempt's container

        def shutdown(signum: int, frame: FrameType | None) -> None:
            logger.info("stopping: no new jobs; stopping %d running", len(running))
            stop.set()
            for job_id in running:
                subprocess.run(["docker", "stop", containers[job_id]], capture_output=True)

        signal.signal(signal.SIGTERM, shutdown)
        signal.signal(signal.SIGINT, shutdown)
        logger.info("worker %s: up to %d jobs", name, site_settings.worker_concurrency())
        self.detect()

        last_sweep = 0.0
        archiver: threading.Thread | None = None  # archives failed jobs' folders
        while not stop.is_set() or running:
            for job_id, thread in list(running.items()):
                if not thread.is_alive():
                    del running[job_id]
                    containers.pop(job_id, None)
            if stop.is_set():
                time.sleep(1)
                continue
            close_old_connections()
            if time.monotonic() - last_sweep > 60:
                if swept := jobs.sweep_stale(settings.PHYLOBOT_STALE_AFTER, set(running)):
                    logger.info("marked %d stale jobs failed", swept)
                if recorded := record_unrecorded(root):
                    logger.info("recorded the compute of %d stopped attempts", recorded)
                if purged := jobs.purge_deleted():
                    logger.info("removed files of %d deleted jobs", purged)
                if measured := jobs.measure_unmeasured():
                    logger.info("measured the storage of %d finished jobs", measured)
                if listed := manifest.index_unindexed():
                    logger.info("indexed the result files of %d finished jobs", listed)
                self.remove_local_folders(root, set(running))
                if expired := expire_work_folders(root, set(running)):
                    logger.info("removed the work folders of %d expired failed jobs", expired)
                measure_work_folders(root, set(running))
                work_archive.tidy()
                self.detect(log=False)
                if archiver is None or not archiver.is_alive():
                    archiver = threading.Thread(
                        target=self.archive, args=(root, set(running)), daemon=True
                    )
                    archiver.start()
                last_sweep = time.monotonic()
            job: Job | None = None
            if len(running) < site_settings.worker_concurrency():
                job = jobs.claim_next(name)
            if job is None:
                time.sleep(POLL_SECONDS)
                continue
            logger.info("job %s: claimed", job.id)
            thread = threading.Thread(target=self.run_one, args=(job, executor, stop), daemon=True)
            running[str(job.id)] = thread
            containers[str(job.id)] = executor.container(job)
            thread.start()

    @staticmethod
    def detect(log: bool = True) -> None:
        """Record the machine for the Settings page (each launch detects it again)."""
        try:
            found = machine.detect()
            machine.record(found)
        except Exception:
            logger.warning("could not detect the machine", exc_info=True)
            return
        if log:
            logger.info("machine: %d CPUs, %s%s", found.cpus, machine.gb(found.memory_bytes),
                        f" ({found.note})" if found.note else "")  # fmt: skip

    @staticmethod
    def remove_local_folders(root: Path, running: set[str]) -> None:
        """Remove kept local folders (failed jobs) once their job is deleted or gone."""
        folders = {p.name: p for p in (root / "jobs").glob("*") if p.is_dir() and is_uuid(p.name)}
        candidates = set(folders) - running
        if not candidates:
            return
        keep = {
            str(job_id)
            for job_id in Job.objects.filter(
                id__in=candidates, deleted_at__isnull=True
            ).values_list("id", flat=True)
        }
        for name in candidates - keep:
            shutil.rmtree(folders[name], ignore_errors=True)
            logger.info("removed local folder of deleted job %s", name)

    @staticmethod
    def archive(root: Path, running: set[str]) -> None:
        try:
            work_archive.archive_pending(root, running)
        except Exception:
            logger.exception("archiving failed jobs' work folders")
        finally:
            close_old_connections()

    @staticmethod
    def run_one(job: Job, executor: LocalDockerExecutor, stop: threading.Event) -> None:
        try:
            run_job(job, executor, stop)
        finally:
            close_old_connections()
