import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('phylobot', '0028_task_memory'),
    ]

    operations = [
        migrations.CreateModel(
            name='CladeCollapse',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('clade', models.CharField(max_length=16)),
                ('parts', models.CharField(max_length=16)),
                ('label', models.CharField(blank=True, default='', max_length=60)),
                ('fill', models.CharField(default='#cfd8e6', max_length=16)),
                ('border_colour', models.CharField(default='#2457a6', max_length=16)),
                ('border_width', models.FloatField(default=1.6)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('job', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='clade_collapses', to='phylobot.job')),
            ],
            options={
                'constraints': [models.UniqueConstraint(fields=('job', 'clade'), name='unique_clade_collapse')],
            },
        ),
    ]
