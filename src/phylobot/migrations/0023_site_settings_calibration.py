from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("phylobot", "0022_job_extends"),
    ]

    operations = [
        migrations.AddField(
            model_name="sitesettings",
            name="calibration_intercept",
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="calibration_branch_length",
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="calibration_tip_distance",
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="calibration_site_rate",
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="calibration_entropy",
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="calibration_gap_fraction",
            field=models.FloatField(blank=True, null=True),
        ),
    ]
