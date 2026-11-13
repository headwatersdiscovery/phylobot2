import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("phylobot", "0032_alignment_tracks"),
    ]

    operations = [
        migrations.CreateModel(
            name="SavedPath",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True, primary_key=True, serialize=False, verbose_name="ID"
                    ),
                ),
                ("name", models.CharField(max_length=80)),
                ("combination", models.CharField(max_length=255)),
                ("start", models.CharField(max_length=255)),
                ("end", models.CharField(max_length=255)),
                ("view", models.JSONField(default=dict)),
                ("created_at", models.DateTimeField(auto_now_add=True)),
                ("updated_at", models.DateTimeField(auto_now=True)),
                (
                    "job",
                    models.ForeignKey(
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="saved_paths",
                        to="phylobot.job",
                    ),
                ),
            ],
            options={
                "ordering": ["name", "id"],
                "constraints": [
                    models.UniqueConstraint(
                        fields=("job", "combination", "start", "end"), name="unique_saved_path"
                    )
                ],
            },
        ),
    ]
