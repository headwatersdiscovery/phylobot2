import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("phylobot", "0023_site_settings_calibration"),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name="UserCalibration",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("intercept", models.FloatField(blank=True, null=True)),
                ("branch_length", models.FloatField(blank=True, null=True)),
                ("tip_distance", models.FloatField(blank=True, null=True)),
                ("site_rate", models.FloatField(blank=True, null=True)),
                ("entropy", models.FloatField(blank=True, null=True)),
                ("gap_fraction", models.FloatField(blank=True, null=True)),
                ("updated_at", models.DateTimeField(auto_now=True)),
                (
                    "user",
                    models.OneToOneField(
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="calibration",
                        to=settings.AUTH_USER_MODEL,
                    ),
                ),
            ],
        ),
    ]
