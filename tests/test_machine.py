"""The machine and each job's share of it (phylobot.machine), the Settings page's "This
machine", and the admin's Work Records (phylobot.work_records)."""

import re
import subprocess
from pathlib import Path
from typing import Any

import pytest
from django.contrib.auth.models import User
from django.test import Client, override_settings

from phylobot import machine, site_settings, work_records
from phylobot.models import Job, JobAttempt, SiteSettings, WorkerMachine
from phylobot.modes import local_user
from phylobot.settings import middleware_for

pytestmark = pytest.mark.django_db
ROOT = Path(__file__).resolve().parents[1]
GIB = 1024**3


def save(**fields: Any) -> None:
    row = SiteSettings.objects.filter(pk=1).first() or SiteSettings(pk=1)
    for name, value in fields.items():
        setattr(row, name, value)
    row.save()


def test_usable_cpus_leave_out_isolated_ones_or_take_the_cpuset() -> None:
    assert machine.cpu_list("0-3,6") == {0, 1, 2, 3, 6}
    assert machine.usable_cpus(12, "", "") == 12
    assert machine.usable_cpus(8, "7", "") == 7  # OrbStack isolates its last CPU
    assert machine.usable_cpus(8, "7", "0-5") == 6  # the cpuset is what tasks get


def test_detect_asks_docker_and_falls_back_to_the_worker(monkeypatch: pytest.MonkeyPatch) -> None:
    found = machine.detect()  # conftest's Docker: 8 CPUs, 16 GiB
    assert (found.cpus, found.memory_bytes, found.note) == (8, 16 * GIB, "")

    def stalled() -> tuple[int, int]:
        raise subprocess.TimeoutExpired("docker", machine.DOCKER_TIMEOUT)

    monkeypatch.setattr(machine, "docker_info", stalled)
    found = machine.detect()
    assert found.cpus >= 1 and "Docker did not answer" in found.note


def test_meminfo() -> None:
    assert machine.meminfo_bytes("MemTotal:        9443584 kB\nMemFree: 1 kB\n") == 9443584 * 1024


def test_totals_detected_then_env_then_the_page(settings: Any) -> None:
    found = machine.Detected(cpus=12, isolated="", memory_bytes=9 * GIB)
    totals = machine.totals(found)
    assert (totals.cpus, totals.cpus_source) == (12, "detected")
    assert totals.memory_bytes == 9 * GIB and totals.reserved_bytes == 2 * GIB
    settings.PHYLOBOT_CPUS, settings.PHYLOBOT_MEMORY_GB = 10, 8.0
    totals = machine.totals(found)
    assert (totals.cpus, totals.memory_bytes) == (10, 8 * GIB)
    assert totals.memory_source == "from PHYLOBOT_MEMORY_GB"
    save(machine_cpus=6, machine_memory_gb=4.5, reserved_memory_gb=0.5)
    totals = machine.totals(found)
    assert (totals.cpus, totals.memory_bytes, totals.reserved_bytes) == (
        6,
        int(4.5 * GIB),
        GIB // 2,
    )
    assert totals.cpus_source == "set on the Settings page"


def test_a_job_gets_its_share_of_what_the_reserve_leaves() -> None:
    totals = machine.Totals(12, "", 15 * GIB, "", 2 * GIB, "")
    alone = machine.share(totals, 1)
    assert (alone.cpus, alone.memory_bytes) == (12, 12 * GIB)  # 15 - 2 - 1 for Nextflow
    assert alone.tasks_at_once == 6  # 2 GB reserved per step
    two = machine.share(totals, 2)
    assert (two.cpus, two.memory_bytes, two.warning) == (6, int(5.5 * GIB), "")
    # Never below one step's reservation, or no step could start; said so.
    crowded = machine.share(totals, 5)
    assert crowded.memory_bytes == machine.DEFAULT_TASK_MEMORY
    assert "Too little memory" in crowded.warning


def test_a_smaller_step_reservation_runs_more_steps_at_once() -> None:
    # The laptop that prompted the setting: 12 CPUs, 13.9 GiB, held to 5 steps at 2 GB.
    def laptop(task: int = machine.DEFAULT_TASK_MEMORY) -> machine.Totals:
        return machine.Totals(12, "", int(13.93 * GIB), "", 2 * GIB, "", task)

    assert machine.share(laptop(), 1).tasks_at_once == 5
    one_gb = machine.share(laptop(GIB), 1)
    assert one_gb.tasks_at_once == 10 and one_gb.task_bytes == GIB
    # Never more steps than CPUs, however small the reservation.
    assert machine.share(laptop(GIB // 2), 1).tasks_at_once == 12
    # Crowded: each job still gets one step's worth of the setting.
    crowded = machine.share(laptop(GIB), 12)
    assert crowded.memory_bytes == GIB and "Too little memory" in crowded.warning


def test_the_step_reservation_comes_from_settings_then_env(settings: Any) -> None:
    found = machine.Detected(12, "", 15 * GIB)
    assert (machine.totals(found).task_bytes, machine.totals(found).task_source) == (
        machine.DEFAULT_TASK_MEMORY,
        "built-in default",
    )
    settings.PHYLOBOT_TASK_MEMORY_GB = 1.5
    assert machine.totals(found).task_bytes == int(1.5 * GIB)
    save(task_memory_gb=1.0)
    totals = machine.totals(found)
    assert (totals.task_bytes, totals.task_source) == (GIB, "set here")
    text = machine.nextflow_config(machine.share(totals, 1))
    assert "[1024.MB * task.attempt, " in text  # what each step reserves, and its retry


def test_jobs_at_once_in_local_mode_is_also_the_users_limit(settings: Any) -> None:
    settings.PHYLOBOT_MODE = "local"
    save(worker_concurrency=4)
    assert machine.slots() == 1  # the default per-user limit
    save(max_concurrent_jobs=3)
    assert machine.slots() == 3
    settings.PHYLOBOT_MODE = "cloud"
    assert machine.slots() == 4


def test_the_default_is_one_job_at_a_time_per_user() -> None:
    assert site_settings.default_limits()["max_concurrent_jobs"] == 1


def test_resources_config_is_numbers_only() -> None:
    text = machine.nextflow_config(machine.Share(cpus=6, memory_bytes=6 * GIB, slots=2))
    assert "cpus   = 6" in text and "memory = '6144 MB'" in text
    assert "[2048.MB * task.attempt, 6144.MB].min()" in text
    text = machine.nextflow_config(
        machine.Share(cpus=6, memory_bytes=6 * GIB, slots=2, task_bytes=GIB)
    )
    assert "[1024.MB * task.attempt, 6144.MB].min()" in text
    code = [line for line in text.splitlines() if not line.strip().startswith("//")]
    assert not re.search(r"[\"$`\\\\]", "\n".join(code))  # nothing a shell could read


def test_task_memory_matches_nextflow_config() -> None:
    config = (ROOT / "pipeline" / "nextflow.config").read_text()
    match = re.search(r"memory\s*=\s*\{ (\d+)\.GB \* task\.attempt \}", config)
    assert match and int(match.group(1)) * GIB == machine.DEFAULT_TASK_MEMORY
    # The built-in default of the setting is the same (no PHYLOBOT_TASK_MEMORY_GB in tests).
    assert site_settings.task_memory_gb() * GIB == machine.DEFAULT_TASK_MEMORY
    assert "task.exitStatus == 137 && task.attempt == 1 ? 'retry'" in config


def test_record_keeps_one_row_per_host() -> None:
    machine.record(machine.Detected(12, "", 9 * GIB))
    machine.record(machine.Detected(12, "", 15 * GIB))
    row = WorkerMachine.objects.get()
    assert machine.latest() == row and row.memory_bytes == 15 * GIB


def test_settings_page_shows_this_machine(client: Client) -> None:
    page = client.get("/settings/").content.decode()
    assert "No worker has reported this machine yet" in page
    assert "not detected yet" in page
    machine.record(machine.Detected(12, "", 9 * GIB))
    page = client.get("/settings/").content.decode()
    assert "This machine" in page and "Each job gets" in page
    assert "12 CPUs, 6.0" in page  # 9 - 2 reserved - 1 for Nextflow
    assert site_settings.source("machine_cpus", None) == "detected"
    assert site_settings.value("machine_cpus") == 12  # a count, not 12.0
    assert site_settings.value("machine_memory_gb") == 9.0
    # Overrides go through the same form as every other setting.
    response = client.post("/settings/", {"reserved_memory_gb": "1.5", "machine_cpus": "10"})
    assert response.status_code == 302
    assert site_settings.reserved_memory_gb() == 1.5
    page = client.get("/settings/").content.decode()
    assert "10 CPUs, 6.5" in page
    assert "up to 3 steps at once" in page and "Reserved by each step" in page
    response = client.post("/settings/", {"reserved_memory_gb": "1.5", "task_memory_gb": "1"})
    assert response.status_code == 302 and site_settings.task_memory_gb() == 1.0
    page = client.get("/settings/").content.decode()
    assert "up to 6 steps at once" in page  # 9 - 1.5 - 1 = 6.5 GB, at 1 GB a step
    assert client.post("/settings/", {"task_memory_gb": "0.5"}).status_code == 200  # below 1
    assert "/settings/work-records/" in page


# Work Records


def finished(name: str, **fields: Any) -> Job:
    defaults: dict[str, Any] = {
        "owner": local_user(), "name": name, "status": Job.Status.SUCCEEDED,
        "params": {"aligners": ["mafft"], "models": ["LG+G4"], "stop_after": "trees"},
        "work_units": 2.0, "cpu_seconds": 300.0,
    }  # fmt: skip
    from django.utils import timezone

    job = Job.objects.create(**(defaults | fields))
    Job.objects.filter(pk=job.pk).update(finished_at=timezone.now())
    job.refresh_from_db()
    return job


def task(name: str, status: str, realtime: str, peak_rss: str = "") -> dict[str, str]:
    return {"name": f"PHYLOBOT:{name}", "status": status, "realtime": realtime, "cpus": "1",
            "peak_rss": peak_rss}  # fmt: skip


TASKS = [
    task("VALIDATE (1)", "COMPLETED", "2s"),
    task("ALIGN_MAFFT (mafft)", "CACHED", "1m", "30.5 MB"),
    task("TREE_IQTREE (mafft/LG+G4)", "FAILED", "9m"),
    task("TREE_IQTREE (mafft/LG+G4)", "COMPLETED", "4m", "1.5 GB"),
]


def test_stages_count_each_step_once_beside_the_estimate() -> None:
    job = finished("sized", input_size={"sequences": 20, "length": 400.0},
                   progress={"tasks": TASKS})  # fmt: skip
    stages = {s.key: s for s in work_records.stages(job)}
    assert stages["setup"].cpu_seconds == 2 and stages["align"].cpu_seconds == 60  # reused counts
    trees = stages["trees"]
    assert (trees.tasks, trees.cpu_seconds) == (1, 240)  # the failed try is the job's only
    assert trees.peak_rss_bytes == int(1.5 * GIB) and trees.estimated_units
    assert stages["align"].peak_rss_bytes == int(30.5 * 1024**2)
    # Without a recorded size: CPU time, no estimate per stage.
    old = finished("old", progress={"tasks": TASKS})
    assert all(s.estimated_units is None for s in work_records.stages(old))


def test_work_records_page_and_csv(client: Client) -> None:
    job = finished("sized", input_size={"sequences": 20, "length": 400.0},
                   progress={"tasks": TASKS})  # fmt: skip
    JobAttempt.objects.create(job=job, number=1, cpus=6, memory_bytes=6 * GIB)
    Job.objects.create(owner=local_user(), name="still running", status=Job.Status.RUNNING,
                       cpu_seconds=5.0)  # fmt: skip
    page = client.get("/settings/work-records/").content.decode()
    assert "sized" in page and "still running" not in page
    assert "20 seq × 400" in page and "given 6 CPUs" in page
    assert "By Stage" in page and "Trees" in page and "1.5\xa0GB" in page  # filesizeformat
    # Sortable and searchable by job (sort.js): each job's <tbody> moves as one.
    assert 'class="results work-records grouped" id="work-records-table"' in page
    assert 'data-table-search="work-records-table"' in page
    assert '<button type="button" class="sort">CPU time used</button>' in page
    assert f'data-sort="{job.cpu_seconds}"' in page
    csv = client.get("/settings/work-records/?format=csv")
    assert csv["Content-Type"] == "text/csv"
    lines = csv.content.decode().splitlines()
    assert lines[0].split(",") == work_records.CSV_COLUMNS
    assert any(",trees,1," in line for line in lines[1:])
    assert client.get("/settings/work-records/?limit=nonsense").status_code == 200


@override_settings(PHYLOBOT_MODE="cloud", MIDDLEWARE=middleware_for("cloud"))
def test_work_records_are_the_admins_only(client: Client) -> None:
    finished("someone's")
    assert client.get("/settings/work-records/").status_code == 302  # login first
    client.force_login(User.objects.create_user("plain"))
    assert client.get("/settings/work-records/").status_code == 404
    client.force_login(User.objects.create_user("admin", is_staff=True))
    assert "someone" in client.get("/settings/work-records/").content.decode()


def test_trace_memory() -> None:
    assert work_records.memory("830.5 MB") == int(830.5 * 1024**2)
    assert work_records.memory("12 KB") == 12 * 1024
    assert work_records.memory("-") is None


def test_command_line_jobs_record_their_estimate_and_size(fake_storage: Any) -> None:
    from django.core.management import call_command

    call_command("submit_job", str(ROOT / "src" / "phylobot" / "examples" / "ime2.fasta"))
    job = Job.objects.get()
    assert job.input_size["sequences"] == 7 and job.work_units


def test_overrides_may_not_exceed_what_was_detected(client: Client) -> None:
    # Before any worker has looked there is nothing to check against.
    assert machine.check_overrides(100, 500.0, None) == {}
    machine.record(machine.Detected(12, "", int(9.006 * GIB)))
    errors = machine.check_overrides(100, 9.5, None)
    assert "at most 12" in errors["machine_cpus"] and "at most 9" in errors["machine_memory_gb"]
    assert machine.check_overrides(12, 9.0, 6.0) == {}  # all of it is allowed
    # The reserve must leave one job room for one step (2 GB) and Nextflow (1 GB).
    assert "at most 6" in machine.check_overrides(None, None, 6.5)["reserved_memory_gb"]
    assert "at most 1" in machine.check_overrides(None, 4.0, 1.5)["reserved_memory_gb"]
    # A smaller step leaves more for the reserve; a step too large to fit is refused.
    assert machine.check_overrides(None, 4.0, 1.5, 1.0) == {}
    errors = machine.check_overrides(None, 9.0, None, 7.5)
    assert "at most 6" in errors["task_memory_gb"]  # 9 - 2 reserve - 1 Nextflow
    response = client.post("/settings/", {"machine_cpus": "100"})
    page = response.content.decode()
    assert response.status_code == 200 and "at most 12" in page
    assert not SiteSettings.objects.exists()
    assert 'max="12"' in client.get("/settings/").content.decode()  # the browser's limit too


def test_env_values_above_the_machine_are_capped(settings: Any) -> None:
    settings.PHYLOBOT_CPUS, settings.PHYLOBOT_MEMORY_GB = 64, 512.0
    totals = machine.totals(machine.Detected(12, "", 9 * GIB))
    assert (totals.cpus, totals.memory_bytes) == (12, 9 * GIB)
    assert totals.cpus_source == "from PHYLOBOT_CPUS, capped at the 12 detected"


def test_settings_sections_collapse(client: Client) -> None:
    page = client.get("/settings/").content.decode()
    assert page.count('class="step collapsible"') >= 6  # every group
    assert 'id="display"' in page and "phylobot/collapse.js" in page
    assert "Open all" not in page and "Close all" not in page
