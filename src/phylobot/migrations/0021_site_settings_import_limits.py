from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("phylobot", "0020_job_imported"),
    ]

    operations = [
        migrations.AddField(
            model_name="sitesettings",
            name="import_max_megabytes",
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="import_max_unpacked_megabytes",
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="import_max_files",
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="import_max_record_megabytes",
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="import_max_libraries",
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="sitesettings",
            name="import_max_exports",
            field=models.PositiveIntegerField(blank=True, null=True),
        ),
    ]
