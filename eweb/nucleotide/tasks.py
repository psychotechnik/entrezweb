from Bio import Entrez, SeqIO
from celery import shared_task
from django.core.files.base import ContentFile

from eweb.nucleotide.models import Nucleotide


@shared_task
def download_nucleotide_task(seq_id: str):
    assert seq_id, "Sequence ID was not provided."
    print(f"download task: {seq_id=}")

    #id = "224589800"

    with Entrez.esummary(db="nucleotide", id=seq_id) as handle:
        record = Entrez.read(handle)[0]

    title = record["Title"]
    extra = record["Extra"]
    length = record["Length"].numerator

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
            title=title,
            extra=extra,
            seq_length=length,
        )
        n.fasta_file.save(
            f"{seq_id}.fasta",
            ContentFile(seq_record.format("fasta")),
        )
        print(f"created nucleotide {n.entrez_id}")
    else:
        print(f"nucleotide {n.entrez_id} has already been downloaded")
