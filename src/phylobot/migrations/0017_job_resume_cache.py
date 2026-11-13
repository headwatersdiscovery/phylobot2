from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("phylobot", "0016_job_attempts"),
    ]

    operations = [
        migrations.AddField(
            model_name="job",
            name="resume_cache",
            field=models.BooleanField(default=False),
        ),
    ]
