from django.db import models

#from eweb.core.fields import SequenceField


class Nucleotide(models.Model):

    entrez_id = models.CharField("Entrez ID", max_length=100, default='')
    name = models.CharField("Name", max_length=100, default='')
    description = models.TextField("Name", max_length=250, default='')
    seq = models.TextField("Sequence", default='')
