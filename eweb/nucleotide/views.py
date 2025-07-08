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
    row_seq_parts: list[SeqPart]


@require_GET
def index(request):
    context = {
        "nucleotides": Nucleotide.objects.order_by("entrez_id"),
    }
    return render(request, 'index.html', context)



class SearchForm(forms.Form):
    seq_search_query = forms.CharField(label="Query", max_length=100)


chars_per_part = 10
num_of_columns = 5
span_red = '<span class="text-red-600">'
span_close = '</span>'


def build_seq_row(
    seq_parts: list[str],
    marker_left,
    marker_right,
    highlight_positions=None,
):
    if not highlight_positions:
        highlight_positions = []
    print(f"{highlight_positions=}")
    seq_markup_values = []
    if marker_left > 1:
        parts_index_start = marker_left // chars_per_part
    else:
        parts_index_start = marker_left - 1
    parts_index_end = marker_right // chars_per_part
    print(f"{marker_left=} {marker_right=} {parts_index_start=} {parts_index_end=}")

    row_str = "".join(seq_parts[parts_index_start:parts_index_end])
    print(f"{row_str=}")
    assert len(row_str) == 50, f"row seq str len: {len(row_str)}"

    for i, val in enumerate(row_str):
        if i in highlight_positions:
            seq_markup_values.append(f'{span_red}{val}{span_close}')
        else:
            seq_markup_values.append(val)

    print(f"{seq_markup_values=}")
    assert len(seq_markup_values) == 50, f"seq markup len: {len(seq_markup_values)}"

    row_seq_parts = []
    markup_start_index = 0
    markup_end_index = chars_per_part 
    index_start = marker_left
    index_end = index_start + chars_per_part - 1
    for c in range(0, num_of_columns):
        print(f"{markup_start_index=} {markup_end_index=}")
        seq_index = marker_left - 1
        col = seq_parts[seq_index // chars_per_part + c]
        seq_part = SeqPart(
            index_start=index_start,
            index_end=index_end,
            seq=col,
            seq_markup=seq_markup_values[markup_start_index:markup_end_index],
        )
        assert len(seq_part.seq_markup) == 10, f"seq markup len: {len(seq_part.seq_markup)}"
        print(seq_part)
        print()
        row_seq_parts.append(seq_part)

        index_start = index_end + 1
        index_end = index_start + chars_per_part - 1

        markup_start_index = markup_end_index
        markup_end_index = markup_start_index + chars_per_part 
    print()
    print()
    return SeqRow(marker_left, marker_right, row_seq_parts)


@require_http_methods(["GET", "POST"])
def seq_table(request: HtmxHttpRequest, uid: str) -> HttpResponse:
    nucleotide = get_object_or_404(Nucleotide, entrez_id=uid)
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

    seq_parts: list[str] = textwrap.TextWrapper(width=chars_per_part).\
        wrap(text=nucleotide.seq)
    rows_as_strs = []
    marker_left, marker_right = 1, chars_per_part * num_of_columns
    row_count = len(seq_parts) // num_of_columns

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
            seq_row = build_seq_row(seq_parts, marker_left, marker_right, highlight_positions)
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
            seq_row = build_seq_row(seq_parts, marker_left, marker_right)
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
            {"nucleotide": nucleotide, "rows_as_strs": rows_as_strs}
        )
        return HttpResponse(block_as_string)
    #else:
    #    return render(request, "index.html", ...)


    #return render(request, 'index.html')
