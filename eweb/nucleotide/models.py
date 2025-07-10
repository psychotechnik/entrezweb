from pathlib import Path

from Bio import SeqIO
from django.db import models


def nucleotide_filepath(instance, filename):
    return Path("nucleotide-fasta") / str(instance.entrez_id) / filename


class Nucleotide(models.Model):

    entrez_id = models.CharField(
        "Entrez ID",
        db_index=True,
        max_length=100,
        unique=True,
    )
    title = models.CharField("Name", max_length=250)
    extra = models.CharField("Extra", max_length=250)
    seq_length = models.IntegerField("Sequence Length", default=0)

    fasta_file = models.FileField(
        verbose_name="Nuncleotide Seq as fasta file",
        upload_to=nucleotide_filepath,
        blank=True,
        null=True,
    )
    class Meta:
        ordering = ('entrez_id',)
        verbose_name = "Nucleotide"

    def __str__(self):
        return f"Nucleotide {self.name} [{self.entrez_id}]"

    def __repr__(self):
        return f'<Nucleotide name="{self.name}" entrez_id="{self.entrez_id}">'

    @property
    def seq(self):
        if not Path(self.fasta_file.file.name).exists():
            raise ValueError("sequence fasta files does not exist on dist.")

        with open(self.fasta_file.file.name, "r") as handle:
            return str(SeqIO.read(handle, "fasta").seq)

#import gzip
#from Bio import SeqIO
#with gzip.open("ls_orchid.gbk.gz", "rt") as handle:
#    print(sum(len(r) for r in SeqIO.parse(handle, "gb")))
