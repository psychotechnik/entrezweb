from Bio import Entrez

from django.core.management.base import BaseCommand

from eweb.nucleotide.models import Nucleotide


class Command(BaseCommand):
    help = "fetch record from Entrez"

    #def add_arguments(self, parser):
    #    parser.add_argument('sample', nargs='+')

    def handle(self, *args, **options):

        id = "30271926"

        handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="xml")
        recs = Entrez.read(handle)
        for r in recs:
            #for k, v in r.items():
            #    print(k, v)
            #print(r["TSeq_sequence"])
            n, created = Nucleotide.objects.get_or_create(
                seq=r["TSeq_sequence"],
            )
            if created: 
                print(f"created nucleotide id: {id}")
