import re

from Bio import Entrez
from django.conf import settings

#from django.shortcuts import get_object_or_404
from django_celery_results.models import TaskResult
from ninja import Router

from eweb.nucleotide.tasks import download_nucleotide_task

from eweb import (
    seq_parts,
    chars_per_part,
    num_of_columns,
    top_rows_num,
    buttom_rows_num,
    span_red,
    span_close,
)

from .data import build_seq_row
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
def get_seq_table(request, seq_id: int, seq_search_query: str | None = None):
    try:
        nucleotide = Nucleotide.objects.get(entrez_id=seq_id)
    except Nucleotide.DoesNotExist:
        return {"result": "nucleotide not found"}

    rows = []
    parts = seq_parts(nucleotide.seq)
    row_count = len(parts) // num_of_columns
    marker_left, marker_right = 1, chars_per_part * num_of_columns
    """
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

    match_position_indexes = [m[0] for m in query_matches]
    print(f"{match_position_indexes=}")
    for row_index in range(0, row_count):
        if query_matches:
            highlight_positions = []
            for match_index, match_seq in query_matches:
                #print(f"{match_index=} {match_seq=}")
                for position in range(len(match_seq)):
                    highlight_positions.append(match_index+position)

            if list(filter(lambda x: (marker_left -1) <= x <= (marker_right - 1), match_position_indexes)):

                if highlight_positions:
                    print(f"{match_position_indexes=}")
                    print(f"api: {highlight_positions=} {marker_left=} {marker_right=} ")

                seq_row = build_seq_row(
                    parts,
                    marker_left,
                    marker_right,
                    markup_style="terminal",
                    highlight_positions=highlight_positions,
                )
                rows.append(seq_row)
        elif (row_index <= 25) or (row_index > (row_count - 25)):
            seq_row = build_seq_row(
                parts,
                marker_left,
                marker_right,
                markup_style="terminal"
            )
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
