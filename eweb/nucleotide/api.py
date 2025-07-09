from ninja import Router

from eweb.nucleotide.tasks import download_nucleotide_task
from django.shortcuts import get_object_or_404
from django_celery_results.models import TaskResult

from .models import Nucleotide

router = Router()

@router.get('/')
def list_nucleotides(request):
    return [
        {
            "entrez_id": n.entrez_id,
            "name": n.name,
            "description": n.description
        }
        for n in Nucleotide.objects.all()
    ]

@router.get('/{seq_id}/')
def nucleotide_details(request, seq_id: int):
    nucleotide = Nucleotide.objects.get(entrez_id=seq_id)
    return {
        "entrez_id": nucleotide.entrez_id,
        "name": nucleotide.name,
        "description": nucleotide.description,
        "sequence": nucleotide.seq,
    }

@router.get('/download/{seq_id}/')
def nucleotide_download(request, seq_id: int):
    result = download_nucleotide_task.delay(seq_id)
    return {"task_id": result.id}


@router.get('/download-progress/{task_id}/')
def nucleotide_download_progress(request, task_id: str):
    task_result = TaskResult.objects.get(task_id=task_id)
    return {"status": task_result.status }

