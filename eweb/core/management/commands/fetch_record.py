from Bio import Entrez, SeqIO

from django.core.management.base import BaseCommand

from eweb.nucleotide.models import Nucleotide


class Command(BaseCommand):
    help = "fetch record from Entrez"

    #def add_arguments(self, parser):
    #    parser.add_argument('sample', nargs='+')

    def handle(self, *args, **options):

        id = "30271926"

        with Entrez.efetch(
            db="nucleotide", rettype="fasta", retmode="fasta", id=id
        ) as handle:
            seq_record = SeqIO.read(handle, "fasta")
        print("%s with %i features" % (seq_record.id, len(seq_record.features)))

        #import ipdb;ipdb.set_trace()

        n, created = Nucleotide.objects.get_or_create(
            entrez_id=seq_record.id,
            name=seq_record.name,
            description=seq_record.description,
            seq=str(seq_record.seq)
        )
        if created: 
            print(f"created nucleotide {n.entrez_id}")

        """
        stream = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="xml")
        recs = Entrez.read(stream)
        for r in recs:
            #for k, v in r.items():
            #    print(k, v)
            #print(r["TSeq_sequence"])
            n, created = Nucleotide.objects.get_or_create(
                seq=r["TSeq_sequence"],
            )
            if created: 
                print(f"created nucleotide id: {id}")
        stream.close()
        """
