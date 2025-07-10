import re

from django.http import (
    HttpRequest,
    HttpResponse,
    HttpResponseRedirect,
)
from django.shortcuts import get_object_or_404, render
from django.urls import reverse
from django.views.decorators.http import require_GET, require_http_methods, require_POST
from django_celery_results.models import TaskResult
from django_htmx.middleware import HtmxDetails
from render_block import render_block_to_string

from eweb.nucleotide.models import Nucleotide
from eweb.nucleotide.tasks import download_nucleotide_task

from .data import (
    build_seq_row,
    chars_per_part,
    num_of_columns,
    seq_parts,
)
from .forms import DownloadNucleotideForm, SearchForm


# Typing pattern recommended by django-stubs:
# https://github.com/typeddjango/django-stubs#how-can-i-create-a-httprequest-thats-guaranteed-to-have-an-authenticated-user
class HtmxHttpRequest(HttpRequest):
    htmx: HtmxDetails



@require_GET
def index(request):
    first_nuc = Nucleotide.objects.first() 
    context = {
        "nucleotides": Nucleotide.objects.order_by("entrez_id") ,
    }
    if first_nuc:
        context["seq_table_url"] = reverse("seq-table", args=[first_nuc.entrez_id])

    return render(request, 'index.html', context)

@require_http_methods(["GET", "POST"])
def seq_table(request: HtmxHttpRequest, seq_id: str) -> HttpResponse:
    nucleotide = get_object_or_404(Nucleotide, entrez_id=seq_id)
    seq_search_query = None
    query_matches = []

    if request.method == "POST":
        form = SearchForm(request.POST)
        if form.is_valid():
            print(f"{form.cleaned_data}")
            #return HttpResponseRedirect("/")
            seq_search_query = form.cleaned_data.get("seq_search_query")
        else:
            print("form invalid")

    rows_as_strs = []

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
    if seq_search_query:
        matches = re.finditer(seq_search_query.replace(' ', ''), nucleotide.seq)
        for match in matches:
            query_matches.append((match.start(), match.group()))
        print(query_matches)

    for row_index in range(0, row_count):
        if query_matches:
            highlight_positions = []
            for match_index, match_seq in query_matches:
                print(f"{match_index=} {match_seq=}")
                for position in range(len(match_seq)):
                    highlight_positions.append(match_index+position)
            seq_row = build_seq_row(parts, marker_left, marker_right, highlight_positions)
            row_as_str = render_block_to_string(
                'includes/seq-row.html',
                'block1', 
                {"nucleotide": nucleotide, "seq_row": seq_row}
            )
            rows_as_strs.append(row_as_str)

        elif (row_index <= 2): # or (row_index > (row_count - 6)):
            #if row_index > num_of_columns and row_index < (row_count - 6):
            #    print("sep")
            #    row_as_str = '<tr class="border-b dark:border-gray-700">sep</tr>'
            #else:
            seq_row = build_seq_row(parts, marker_left, marker_right)
            row_as_str = render_block_to_string(
                'includes/seq-row.html',
                'block1', 
                {"nucleotide": nucleotide, "seq_row": seq_row}
            )
            rows_as_strs.append(row_as_str)

        marker_left+=chars_per_part*num_of_columns
        marker_right+=chars_per_part*num_of_columns 

    if request.htmx:
        block_as_string = render_block_to_string(
            'includes/seq-table.html',
            'block1',
            {
                "nucleotide": nucleotide,
                "rows_as_strs": rows_as_strs,
                "seq_table_url": reverse("seq-table", args=[nucleotide.entrez_id])
                #f"/nucleotides/seq-table/{nucleotide.entrez_id}/" 
            }
        )
        return HttpResponse(block_as_string)
    #else:
    #    return render(request, "index.html", ...)


    #return render(request, 'index.html')

@require_POST
def download_nucleotide(request: HtmxHttpRequest) -> HttpResponse:
    seq_id = None
    if request.method == "POST" and request.htmx:
        form = DownloadNucleotideForm(request.POST)
        if form.is_valid():
            print(f"{form.cleaned_data}")
            seq_id = form.cleaned_data.get("seq_id")
            result = download_nucleotide_task.delay(seq_id)
            #res.get(timeout=1)#
            block_as_string = render_block_to_string(
                'includes/seq-download-progress.html',
                'block1',
                {"seq_id": seq_id, "task_id": result.id}
            )
            return HttpResponse(block_as_string)
        else:
            print("form invalid")

    return HttpResponseRedirect("/")


@require_GET
def download_nucleotide_progress(request: HtmxHttpRequest, task_id: str, seq_id: str) -> HttpResponse:
    task_result = get_object_or_404(TaskResult, task_id=task_id)
    print(f"{task_id=} {seq_id=} {task_result.status=}")

    if task_result.status == "SUCCESS":
        nucleotide = get_object_or_404(Nucleotide, entrez_id=seq_id)
        rows_as_strs = []
        marker_left, marker_right = 1, chars_per_part * num_of_columns
        parts = seq_parts(nucleotide.seq)
        row_count = len(parts) // num_of_columns

        for row_index in range(0, row_count):
            if (row_index <= 2): # or (row_index > (row_count - 6)):
                seq_row = build_seq_row(parts, marker_left, marker_right)
                row_as_str = render_block_to_string(
                    'includes/seq-row.html',
                    'block1', 
                    {"nucleotide": nucleotide, "seq_row": seq_row}
                )
                rows_as_strs.append(row_as_str)

            marker_left+=chars_per_part*num_of_columns
            marker_right+=chars_per_part*num_of_columns 

        block_as_string = render_block_to_string(
            'includes/seq-table.html',
            'block1',
            {
                "nucleotide": nucleotide,
                "rows_as_strs": rows_as_strs,
                "close_modal": True,
                "seq_table_url": reverse("seq-table", args=[nucleotide.entrez_id]),
                #"seq_table_url": f"/nucleotides/seq-table/{nucleotide.entrez_id}/" 
            }
        )
        return HttpResponse(block_as_string)

    block_as_string = render_block_to_string(
        'includes/seq-download-progress.html',
        'block1',
        {"seq_id": seq_id, "task_id": task_id}
    )
    return HttpResponse(block_as_string)



