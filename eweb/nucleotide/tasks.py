from Bio import Entrez, SeqIO
from celery import shared_task

from eweb.nucleotide.models import Nucleotide


@shared_task
def download_nucleotide_task(seq_id: str):
    assert seq_id, "Sequence ID was not provided."
    print(f"download task: {seq_id=}")

    #id = "224589800"

    with Entrez.efetch(
        db="nucleotide", rettype="fasta", retmode="fasta", id=seq_id
    ) as handle:
        seq_record = SeqIO.read(handle, "fasta")

    print("%s with %i features" % (seq_record.id, len(seq_record.features)))

    n, created = Nucleotide.objects.get_or_create(
        entrez_id=seq_id,
        name=seq_record.name,
        description=seq_record.description,
        seq=str(seq_record.seq)
    )
    if created: 
        print(f"created nucleotide {n.entrez_id}")
    else:
        print(f"nucleotide {n.entrez_id} has already been downloaded")
