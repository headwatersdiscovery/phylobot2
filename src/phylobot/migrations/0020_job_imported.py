from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("phylobot", "0019_job_work_archive"),
    ]

    operations = [
        migrations.AddField(
            model_name="job",
            name="imported_at",
            field=models.DateTimeField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="job",
            name="imported_from",
            field=models.JSONField(blank=True, default=dict),
        ),
    ]
