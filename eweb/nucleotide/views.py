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

from eweb import (
    buttom_rows_num,
    chars_per_part,
    num_of_columns,
    seq_parts,
    span_close,
    span_red,
    top_rows_num,
    Timer,

)
from eweb.nucleotide.data import build_seq_row
from eweb.nucleotide.models import Nucleotide
from eweb.nucleotide.tasks import download_nucleotide_task

from .forms import DownloadNucleotideForm, SearchForm


# Typing pattern recommended by django-stubs:
# https://github.com/typeddjango/django-stubs#how-can-i-create-a-httprequest-thats-guaranteed-to-have-an-authenticated-user
class HtmxHttpRequest(HttpRequest):
    htmx: HtmxDetails



@require_GET
def index(request):
    #first_nuc = Nucleotide.objects.get(entrez_id="3008860383")
    first_nuc = Nucleotide.objects.order_by('seq_length').first() 
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
            seq_search_query = form.cleaned_data.get("seq_search_query")
        else:
            print("form invalid")
            #import ipdb;ipdb.set_trace()

            block_as_string = render_block_to_string(
                'includes/seq-table.html',
                'block1',
                {
                    "form_errors": form.errors,
                    #"nucleotide": nucleotide,
                    #"rows_as_strs": rows_as_strs,
                    #"seq_table_url": reverse("seq-table", args=[nucleotide.entrez_id])
                    #f"/nucleotides/seq-table/{nucleotide.entrez_id}/" 
                }
            )
            return HttpResponse(block_as_string)

    rows_as_strs = []
    #parts = seq_parts(nucleotide.seq)
    #row_count = len(parts) // num_of_columns
    if not seq_search_query:
        row_count = top_rows_num
        marker_left, marker_right = 1, chars_per_part * num_of_columns
        sub_seq = nucleotide.seq[:(top_rows_num*num_of_columns*chars_per_part)-1] 
        for row_index in range(0, row_count):
            if (row_index <= top_rows_num) or (row_index > (row_count - buttom_rows_num)):
                idx_start = marker_left-1
                idx_end = num_of_columns*chars_per_part+marker_left
                seq_row = build_seq_row(
                    seq_parts(sub_seq[idx_start:idx_end]),
                    marker_left,
                    marker_right
                )
                row_as_str = render_block_to_string(
                    'includes/seq-row.html',
                    'block1', 
                    {"nucleotide": nucleotide, "seq_row": seq_row}
                )
                rows_as_strs.append(row_as_str)
                idx_start = idx_end+1
                idx_end = num_of_columns*chars_per_part + idx_start

            marker_left+=chars_per_part*num_of_columns
            marker_right+=chars_per_part*num_of_columns 

        #if request.htmx:
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

    # SEARCH QUERY
    matches = re.finditer(seq_search_query.replace(' ', ''), nucleotide.seq)
    for match in matches:
        query_matches.append((match.start(), match.group()))

    print(f"{query_matches=}")

    highlight_positions = []
    for match_row_index, match_seq in query_matches:
        for position in range(len(match_seq)):
            highlight_positions.append(match_row_index+position)

    if highlight_positions:
        print(f"view: {highlight_positions=}")

    match_indexes = [m[0] for m in query_matches]
    print(f"{match_indexes=}")
    if len(match_indexes) > 25:
        block_as_string = render_block_to_string(
            'includes/seq-table.html',
            'block1',
            {"too_many_results": len(match_indexes) }
        )
        return HttpResponse(block_as_string)

    #query_matches=[
    #    (23, 'GTTGCTC'), 
    #    (3276, 'GTTGCTC')
    #]
    #highlight_positions=[23, 24, 25, 26, 27, 28, 29, 3276, 3277, 3278, 3279, 3280, 3281, 3282]

    # CATTTAAAGCAGTGTGTAAAGAGAC
    row_count = 0
    ml, mr = 1, chars_per_part * num_of_columns
    for row_index in range(0, len(nucleotide.seq) // (num_of_columns * chars_per_part)):
        if list(filter(lambda x: ml-1 <= x < mr, match_indexes)):
            row_count += 1
        ml+=chars_per_part*num_of_columns
        mr+=chars_per_part*num_of_columns 
    print(f"{row_count=}")
    t = Timer()
    t.start()

    marker_left, marker_right = 1, chars_per_part*num_of_columns
    all_rows_count = len(nucleotide.seq) // (chars_per_part*num_of_columns)
    for row_index in range(0, all_rows_count):
        #if list(filter(lambda x: marker_left-1 <= x < marker_right, match_indexes)):
        for match_index in match_indexes:
            if marker_left-1 <= match_index <= marker_right-1:
                #print(f"{match_index=} {marker_left=} {marker_right=}")
                idx_start = marker_left-1
                idx_end = marker_right
                seq_row = build_seq_row(
                    seq_parts(nucleotide.seq[idx_start:idx_end]),
                    marker_left,
                    marker_right,
                    highlight_positions=highlight_positions,
                )
                row_as_str = render_block_to_string(
                    'includes/seq-row.html',
                    'block1', 
                    {"nucleotide": nucleotide, "seq_row": seq_row}
                )
                rows_as_strs.append(row_as_str)
                idx_start = idx_end+1
                idx_end = num_of_columns*chars_per_part+idx_start

        marker_left+=chars_per_part*num_of_columns
        marker_right+=chars_per_part*num_of_columns 
    t.stop()

    #if request.htmx:
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
    #return render(request, 'index.html')

@require_POST
def download_nucleotide(request: HtmxHttpRequest) -> HttpResponse:
    seq_id = None
    if request.method == "POST" and request.htmx:
        form = DownloadNucleotideForm(request.POST)
        if form.is_valid():
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
                seq_row = build_seq_row( parts, marker_left, marker_right)
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



