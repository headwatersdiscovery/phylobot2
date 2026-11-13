from typing import Any

from django.db import migrations, models
import django.db.models.deletion


def first_attempts(apps: Any, schema_editor: Any) -> None:
    """Every job that has started gets its attempt 1, from the Job row's own figures.
    Older jobs never split failed from completed steps, so all of their CPU time is
    recorded as completed and billable, as cpu_seconds always counted it."""
    Job = apps.get_model("phylobot", "Job")
    JobAttempt = apps.get_model("phylobot", "JobAttempt")
    rows = []
    for job in Job.objects.filter(started_at__isnull=False).iterator():
        terminal = job.status in ("succeeded", "failed", "cancelled")
        rows.append(
            JobAttempt(
                job_id=job.id,
                number=1,
                worker=job.worker,
                started_at=job.started_at,
                finished_at=job.finished_at,
                status=job.status if terminal else "",
                error=job.error,
                completed_cpu_seconds=job.cpu_seconds,
                failed_cpu_seconds=0 if job.cpu_seconds is not None else None,
                billable_cpu_seconds=job.cpu_seconds,
            )
        )
    JobAttempt.objects.bulk_create(rows, batch_size=500)
    Job.objects.filter(cpu_seconds__isnull=False).update(
        billable_cpu_seconds=models.F("cpu_seconds")
    )


class Migration(migrations.Migration):
    dependencies = [
        ("phylobot", "0015_job_storage_bytes"),
    ]

    operations = [
        migrations.AddField(
            model_name="job",
            name="attempt",
            field=models.PositiveIntegerField(default=1),
        ),
        migrations.AddField(
            model_name="job",
            name="billable_cpu_seconds",
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="failed_task_rate_percent",
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
        migrations.CreateModel(
            name="JobAttempt",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("number", models.PositiveIntegerField()),
                ("worker", models.CharField(blank=True, max_length=100)),
                ("started_at", models.DateTimeField(blank=True, null=True)),
                ("finished_at", models.DateTimeField(blank=True, null=True)),
                ("status", models.CharField(blank=True, choices=[("queued", "Queued"), ("running", "Running"), ("succeeded", "Succeeded"), ("failed", "Failed"), ("cancelled", "Cancelled")], max_length=16)),
                ("error", models.TextField(blank=True)),
                ("completed_cpu_seconds", models.FloatField(blank=True, null=True)),
                ("failed_cpu_seconds", models.FloatField(blank=True, null=True)),
                ("tasks_run", models.PositiveIntegerField(blank=True, null=True)),
                ("tasks_reused", models.PositiveIntegerField(blank=True, null=True)),
                ("failed_rate_percent", models.PositiveIntegerField(blank=True, null=True)),
                ("billable_cpu_seconds", models.FloatField(blank=True, null=True)),
                ("job", models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name="attempts", to="phylobot.job")),
            ],
            options={"ordering": ["job", "number"]},
        ),
        migrations.AddConstraint(
            model_name="jobattempt",
            constraint=models.UniqueConstraint(fields=("job", "number"), name="unique_job_attempt"),
        ),
        migrations.RunPython(first_attempts, migrations.RunPython.noop),
    ]
