from Bio import Entrez
from django.conf import settings

#from django.shortcuts import get_object_or_404
from django_celery_results.models import TaskResult
from ninja import Router

from eweb.nucleotide.tasks import download_nucleotide_task

from .models import Nucleotide

Entrez.email = settings.ADMINS[0][1]

router = Router()


@router.get('/')
def list_(request, seq_ids: str | None = None):

    results = []
    if seq_ids:
        seq_ids: list[str] = seq_ids.split(",")
        try:
            with Entrez.esummary(db="nucleotide", id=seq_ids) as handle:
                for record in Entrez.read(handle):
                    try:
                        n = Nucleotide.objects.get(entrez_id=record["Id"])
                        print(f"getting from db {n.entrez_id}")
                        results.append({
                            "entrez_id": n.entrez_id,
                            "title": n.title,
                            "extra": n.extra,
                            "seq_length": n.seq_length,
                            "downloaded": True,
                        })
                    except Nucleotide.DoesNotExist:
                        results.append({
                            "entrez_id": record["Id"],
                            "title": record["Title"],
                            "extra": record["Extra"],
                            "seq_length": record["Length"].numerator,
                            "downloaded": False,
                        })
        except RuntimeError:
            pass
            # FIXME
            # handle invalid ids
            #RuntimeError: Invalid uid 1230000000000000000000000000 at position= 32

        return results

    return [{
        "entrez_id": n.entrez_id,
        "title": n.title,
        "extra": n.extra,
        "seq_length": n.seq_length,
        "downloaded": True,
    } for n in Nucleotide.objects.all()]

@router.get('/{seq_id}/')
def get_(request, seq_id: int):
    try:
        nucleotide = Nucleotide.objects.get(entrez_id=seq_id)
    except Nucleotide.DoesNotExist:
        return {}

    return {
        "entrez_id": nucleotide.entrez_id,
        "title": nucleotide.title,
        "extra": nucleotide.extra,
        "seq_length": nucleotide.seq_length,
    }

@router.get('/download/{seq_id}/')
def download(request, seq_id: int):
    result = download_nucleotide_task.delay(seq_id)
    return {"task_id": result.id}


@router.get('/download-progress/{task_id}/')
def download_progress(request, task_id: str):
    try:
        task_result = TaskResult.objects.get(task_id=task_id)
    except TaskResult.DoesNotExist:
        return {"status": "in-progres"}

    return {"status": task_result.status }

