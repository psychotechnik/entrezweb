import re

from Bio import Entrez
from django.conf import settings

#from django.shortcuts import get_object_or_404
from django_celery_results.models import TaskResult
from ninja import Router

from eweb.nucleotide.tasks import download_nucleotide_task

from .data import (
    build_seq_row,
    chars_per_part,
    num_of_columns,
    seq_parts,
)
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


@router.get('/seq-table/{seq_id}/')
def get_seq_table(request, seq_id: int, seq_search_query: str):
    try:
        nucleotide = Nucleotide.objects.get(entrez_id=seq_id)
    except Nucleotide.DoesNotExist:
        return {"result": "nucleotide not found"}

    #nucleotide.entrez_id,
    #nucleotide.title,
    #nucleotide.extra,
    #nucleotide.seq_length,

    rows = []
    parts = seq_parts(nucleotide.seq)
    row_count = len(parts) // num_of_columns
    marker_left, marker_right = 1, chars_per_part * num_of_columns

    """
    query_matches = [match.start() for match in matches]
    {'seq_search_query': 'ATATTAGGTT'}

    TGCATGCCTA|GTGCACCTAC
    100 TGCATGCCTA
    110 GTGCACCTAC
    16872 GTGCACCTAC
    """

    query_matches = []
    if seq_search_query:
        matches = re.finditer(seq_search_query.replace(' ', ''), nucleotide.seq)
        for match in matches:
            query_matches.append((match.start(), match.group()))
        print(query_matches)

    for row_index in range(0, row_count):
        if 0: #query_matches:
            highlight_positions = []
            #for match_index, match_seq in query_matches:
            #    print(f"{match_index=} {match_seq=}")
            #    for position in range(len(match_seq)):
            #        highlight_positions.append(match_index+position)
            #seq_row = build_seq_row(parts, marker_left, marker_right, highlight_positions)
            #row_as_str = render_block_to_string(
            #    'includes/seq-row.html',
            #    'block1', 
            #    {"nucleotide": nucleotide, "seq_row": seq_row}
            #)
            #rows_as_strs.append(row_as_str)

        elif (row_index < 1): # or (row_index > (row_count - 6)):

            seq_row = build_seq_row(parts, marker_left, marker_right, markup_style="terminal")
            print(seq_row.model_dump)
            #import ipdb;ipdb.set_trace()

            #row_as_str = render_block_to_string(
            #    'includes/seq-row.html',
            #    'block1', 
            #    {"nucleotide": nucleotide, "seq_row": seq_row}
            #)
            rows.append(seq_row)

        marker_left+=chars_per_part*num_of_columns
        marker_right+=chars_per_part*num_of_columns 

    return rows

    #if request.htmx:
    #    block_as_string = render_block_to_string(
    #        'includes/seq-table.html',
    #        'block1',
    #        {"nucleotide": nucleotide, "rows_as_strs": rows_as_strs}
    #    )
