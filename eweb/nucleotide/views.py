import re
import textwrap
from dataclasses import dataclass

from django import forms
from django.http import (
    HttpRequest,
    HttpResponse,
    #HttpResponseRedirect,
)
from django.shortcuts import get_object_or_404, render
from django.views.decorators.http import require_GET, require_http_methods
from render_block import render_block_to_string
from django_htmx.middleware import HtmxDetails

from eweb.nucleotide.models import Nucleotide


# Typing pattern recommended by django-stubs:
# https://github.com/typeddjango/django-stubs#how-can-i-create-a-httprequest-thats-guaranteed-to-have-an-authenticated-user
class HtmxHttpRequest(HttpRequest):
    htmx: HtmxDetails


@dataclass
class SeqPart:
    index_start: int
    index_end: int
    seq: str
    seq_markup: str


@dataclass
class SeqRow:
    marker_left: int
    marker_right: int
    seq_parts: list[SeqPart]


@require_GET
def index(request):
    return render(request, 'index.html')



class SearchForm(forms.Form):
    seq_search_query = forms.CharField(label="Query", max_length=100)


@require_http_methods(["GET", "POST"])
def seq_table(request: HtmxHttpRequest, uid: str) -> HttpResponse:
    nucleotide = get_object_or_404(Nucleotide, entrez_id=uid)
    chars_per_part = 10
    num_of_columns = 5

    seq_search_query = None
    query_matches = None
    span_red = '<span class="text-red-600">'
    span_close = '</span>'

    if request.method == "POST":
        form = SearchForm(request.POST)
        if form.is_valid():
            print(f"{form.cleaned_data}")
            #return HttpResponseRedirect("/")
            seq_search_query = form.cleaned_data.get("seq_search_query")
        else:
            print("form invalid")

    if seq_search_query:
        matches = re.finditer(seq_search_query, nucleotide.seq)
        #for match in matches:
        #    print(match.start(), match.group())
        query_matches = [match.start() for match in matches]

    seq_parts: list[str] = textwrap.TextWrapper(width=chars_per_part).\
        wrap(text=nucleotide.seq)
    rows_as_strs = []
    marker_left, marker_right = 1, chars_per_part * num_of_columns
    row_count = len(seq_parts) // num_of_columns
    for row_index in range(0, row_count):
        if  (row_index < 5) or (row_index > (row_count - 6)):
            seq_index = marker_left - 1
            print(f"{row_index=} {marker_left=} {marker_right=} {row_count=}")

            row1 = seq_parts[seq_index // chars_per_part]
            row1_str = f'{span_red}{row1}{span_close}'
            index_start = marker_left
            index_end = index_start + chars_per_part - 1
            part1 = SeqPart(index_start=index_start, index_end=index_end, seq=row1, seq_markup=row1_str)

            row2 = seq_parts[seq_index // chars_per_part + 1]
            row2_str = f'{span_red}{row2}{span_close}'
            index_start = index_end + 1
            index_end = index_start + chars_per_part - 1
            part2 = SeqPart(index_start=index_start, index_end=index_end, seq=row2, seq_markup=row2_str)

            row3 = seq_parts[seq_index // chars_per_part + 2]
            row3_str = f'{span_red}{row3}{span_close}'
            index_start = index_end + 1
            index_end = index_start + chars_per_part - 1
            part3 = SeqPart(index_start=index_start, index_end=index_end, seq=row3, seq_markup=row3_str)

            row4 = seq_parts[seq_index // chars_per_part + 3]
            row4_str = f'{span_red}{row4}{span_close}'
            index_start = index_end + 1
            index_end = index_start + chars_per_part - 1
            part4 = SeqPart(index_start=index_start, index_end=index_end, seq=row4, seq_markup=row4_str)

            row5 = seq_parts[seq_index // chars_per_part + 4]
            row5_str = f'{span_red}{row5}{span_close}'
            index_start = index_end + 1
            index_end = index_start + chars_per_part - 1
            part5 = SeqPart(index_start=index_start, index_end=index_end, seq=row5, seq_markup=row5_str)

            seq_row = SeqRow(marker_left, marker_right, seq_parts=[part1, part2, part3, part4, part5])

            #if row_index > num_of_columns and row_index < (row_count - 6):
            #    print("sep")
            #    row_as_str = '<tr class="border-b dark:border-gray-700">sep</tr>'
            #else:
            row_as_str = render_block_to_string(
                    'includes/seq-row.html',
                    'block1', 
                    {"nucleotide": nucleotide, "seq_row": seq_row}
                )
            rows_as_strs.append(row_as_str)
        marker_left+=50
        marker_right+=50

    if request.htmx:
        block_as_string = render_block_to_string(
            'includes/seq-table.html',
            'block1',
            {"nucleotide": nucleotide, "rows_as_strs": rows_as_strs}
        )
        return HttpResponse(block_as_string)
    #else:
    #    return render(request, "index.html", ...)


    #return render(request, 'index.html')
