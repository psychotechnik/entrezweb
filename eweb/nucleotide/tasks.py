from Bio import Entrez, SeqIO
from celery import shared_task
from django.core.files.base import ContentFile

from eweb.nucleotide.models import Nucleotide


@shared_task
def download_nucleotide_task(seq_id: str):
    assert seq_id, "Sequence ID was not provided."
    print(f"download task: {seq_id=}")

    #id = "224589800"

    #http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=30271926&rettype=fasta

    with Entrez.efetch(
        db="nucleotide",
        rettype="fasta",
        retmode="text",
        id=seq_id,
    ) as handle:
        seq_record = SeqIO.read(handle, "fasta")

    try:
        n = Nucleotide.objects.get(
            entrez_id=seq_id,
        )
    except Nucleotide.DoesNotExist:
        n = Nucleotide.objects.create(
            entrez_id=seq_id,
            name=seq_record.name,
            description=seq_record.description,
            seq=str(seq_record.seq),
        )
        n.fasta_file.save(
            f"{seq_id}.fasta",
            ContentFile(seq_record.format("fasta")),
        )
        print(f"created nucleotide {n.entrez_id}")
    else:
        print(f"nucleotide {n.entrez_id} has already been downloaded")
