# Generated by Django 5.2.4 on 2025-07-04 13:46

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('nucleotide', '0002_alter_nucleotide_seq'),
    ]

    operations = [
        migrations.AddField(
            model_name='nucleotide',
            name='description',
            field=models.TextField(default='', max_length=250, verbose_name='Name'),
        ),
        migrations.AddField(
            model_name='nucleotide',
            name='entrez_id',
            field=models.CharField(default='', max_length=100, verbose_name='Entrez ID'),
        ),
        migrations.AddField(
            model_name='nucleotide',
            name='name',
            field=models.CharField(default='', max_length=100, verbose_name='Name'),
        ),
        migrations.AlterField(
            model_name='nucleotide',
            name='seq',
            field=models.TextField(default='', verbose_name='Sequence'),
        ),
    ]
