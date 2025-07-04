from django.db import models

#from eweb.core.fields import SequenceField


class Nucleotide(models.Model):

    seq = models.TextField("Sequence", blank=True, default='')
