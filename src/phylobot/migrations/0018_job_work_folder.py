from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("phylobot", "0017_job_resume_cache"),
    ]

    operations = [
        migrations.AddField(
            model_name="job",
            name="work_bytes",
            field=models.BigIntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="job",
            name="work_expired_at",
            field=models.DateTimeField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="keep_failed_work_days",
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
    ]
