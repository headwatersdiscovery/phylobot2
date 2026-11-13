from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("phylobot", "0018_job_work_folder"),
    ]

    operations = [
        migrations.AddField(
            model_name="job",
            name="work_archive_attempt",
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="job",
            name="work_archive_bytes",
            field=models.BigIntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="archive_failed_work",
            field=models.BooleanField(blank=True, null=True),
        ),
    ]
