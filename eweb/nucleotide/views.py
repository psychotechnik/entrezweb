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
    return render(request, 'index.html')



class SearchForm(forms.Form):
    seq_search_query = forms.CharField(label="Query", max_length=100)


@require_http_methods(["GET", "POST"])
def seq_table(request: HtmxHttpRequest, uid: str) -> HttpResponse:
    nucleotide = get_object_or_404(Nucleotide, entrez_id=uid)
    chars_per_part = 10
    num_of_columns = 5

    seq_search_query = None
    query_matches = []
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

    seq_parts: list[str] = textwrap.TextWrapper(width=chars_per_part).\
        wrap(text=nucleotide.seq)
    rows_as_strs = []
    marker_left, marker_right = 1, chars_per_part * num_of_columns
    row_count = len(seq_parts) // num_of_columns

    def build_seq_row(highlight_positions=None):
        if not highlight_positions:
            highlight_positions = []
        print(f"{highlight_positions=}")
        #['ATATTAGGTT', 'TTTACCTACC', 'CAGGAAAAGC', 'CAACCAACCT', 'CGATCTCTTG']
        seq_markup_values = []
        row_str = "".join(seq_parts[marker_left-1: marker_left-1 // chars_per_part + 4])
        for i, val in enumerate(row_str):
            if i in highlight_positions:
                seq_markup_values.append(f'{span_red}{val}{span_close}')
            else:
                seq_markup_values.append(val)

        for seq_markup in seq_markup_values:
            print(seq_markup)
        row_seq_parts = []
        markup_start_index=0
        markup_end_index=num_of_columns-1
        for c in range(0, num_of_columns):
            print(f"{markup_start_index=} {markup_end_index}")
            index_start = marker_left
            index_end = index_start + chars_per_part - 1

            seq_index = marker_left - 1
            col = seq_parts[seq_index // chars_per_part + c]
            seq_part = SeqPart(
                index_start=index_start,
                index_end=index_end,
                seq=col,
                seq_markup=seq_markup_values[markup_start_index:markup_end_index],
            )
            print(seq_part)
            markup_start_index=markup_end_index
            markup_end_index=markup_start_index+chars_per_part

        return SeqRow(marker_left, marker_right, row_seq_parts)

        """
        seq_index = marker_left - 1
        print(f"{row_index=} {marker_left=} {marker_right=} {row_count=}")

        col1 = seq_parts[seq_index // chars_per_part]
        col1_str = ""
        for char_val in col1:
            if col1.index(char_val) in highlight_positions:
                col1_str+=f'{span_red}{char_val}{span_close}'
            else:
                col1_str+=char_val

        index_start = marker_left
        index_end = index_start + chars_per_part - 1
        part1 = SeqPart(index_start=index_start, index_end=index_end, seq=col1, seq_markup=col1_str)
        """

        #return SeqRow(marker_left, marker_right, seq_parts=[
        #    part1
        #])


        """
        seq_index = marker_left - 1
        print(f"{row_index=} {marker_left=} {marker_right=} {row_count=}")

        col1 = seq_parts[seq_index // chars_per_part]
        col1_str = ""
        for char_val in col1:
            if col1.index(char_val) in highlight_positions:
                col1_str+=f'{span_red}{char_val}{span_close}'
            else:
                col1_str+=char_val

        index_start = marker_left
        index_end = index_start + chars_per_part - 1
        part1 = SeqPart(index_start=index_start, index_end=index_end, seq=col1, seq_markup=col1_str)

        col2 = seq_parts[seq_index // chars_per_part + 1]
        col2_str = ""
        for char_val in col2:
            if col2.index(char_val) in highlight_positions:
                col2_str+=f'{span_red}{char_val}{span_close}'
            else:
                col2_str+=char_val
        index_start = index_end + 1
        index_end = index_start + chars_per_part - 1
        part2 = SeqPart(index_start=index_start, index_end=index_end, seq=col2, seq_markup=col2_str)

        col3 = seq_parts[seq_index // chars_per_part + 2]
        col3_str = ""
        for char_val in col3:
            if col3.index(char_val) in highlight_positions:
                col3_str+=f'{span_red}{char_val}{span_close}'
            else:
                col3_str+=char_val
        index_start = index_end + 1
        index_end = index_start + chars_per_part - 1
        part3 = SeqPart(index_start=index_start, index_end=index_end, seq=col3, seq_markup=col3_str)

        col4 = seq_parts[seq_index // chars_per_part + 3]
        col4_str = ""
        for char_val in col4:
            if col4.index(char_val) in highlight_positions:
                col4_str+=f'{span_red}{char_val}{span_close}'
            else:
                col4_str+=char_val
        index_start = index_end + 1
        index_end = index_start + chars_per_part - 1
        part4 = SeqPart(index_start=index_start, index_end=index_end, seq=col4, seq_markup=col4_str)

        col5 = seq_parts[seq_index // chars_per_part + 4]
        col5_str = ""
        for char_val in col5:
            if col5.index(char_val) in highlight_positions:
                col5_str+=f'{span_red}{char_val}{span_close}'
            else:
                col5_str+=char_val
        index_start = index_end + 1
        index_end = index_start + chars_per_part - 1
        part5 = SeqPart(index_start=index_start, index_end=index_end, seq=col5, seq_markup=col5_str)
        """

        """
        col2 = seq_parts[seq_index // chars_per_part + 1]
        #col2_str = f'{span_red}{col2}{span_close}'
        col2_str = col2
        index_start = index_end + 1
        index_end = index_start + chars_per_part - 1
        part2 = SeqPart(index_start=index_start, index_end=index_end, seq=col2, seq_markup=col2_str)

        col3 = seq_parts[seq_index // chars_per_part + 2]
        #col3_str = f'{span_red}{col3}{span_close}'
        col3_str = col3
        index_start = index_end + 1
        index_end = index_start + chars_per_part - 1
        part3 = SeqPart(index_start=index_start, index_end=index_end, seq=col3, seq_markup=col3_str)

        col4 = seq_parts[seq_index // chars_per_part + 3]
        #col4_str = f'{span_red}{col4}{span_close}'
        col4_str = col4
        index_start = index_end + 1
        index_end = index_start + chars_per_part - 1
        part4 = SeqPart(index_start=index_start, index_end=index_end, seq=col4, seq_markup=col4_str)

        col5 = seq_parts[seq_index // chars_per_part + 4]
        #col5_str = f'{span_red}{col5}{span_close}'
        col5_str = col5
        index_start = index_end + 1
        index_end = index_start + chars_per_part - 1
        part5 = SeqPart(index_start=index_start, index_end=index_end, seq=col5, seq_markup=col5_str)
        """

        #return SeqRow(marker_left, marker_right, seq_parts=[part1, part2, part3, part4, part5])

    """
    query_matches = [match.start() for match in matches]
    {'seq_search_query': 'ATATTAGGTT'}
    [0]
    row_index=0 marker_left=1 marker_right=50 row_count=595
    row_index=1 marker_left=51 marker_right=100 row_count=595
    row_index=2 marker_left=101 marker_right=150 row_count=595
    row_index=3 marker_left=151 marker_right=200 row_count=595
    row_index=4 marker_left=201 marker_right=250 row_count=595

    TGCATGCCTA|GTGCACCTAC
    100 TGCATGCCTA
    110 GTGCACCTAC
    16872 GTGCACCTAC

    """
    if seq_search_query:
        matches = re.finditer(seq_search_query, nucleotide.seq)
        for match in matches:
            query_matches.append((match.start(), match.group()))
        print(query_matches)
    for row_index in range(0, row_count):
        if query_matches:
            highlight_positions = []
            if list(
                filter(lambda m: marker_left-1 <= m < marker_right-1, [m[0] for m in query_matches])
            ):
                for match_index, match_seq in query_matches:
                    print(f"{match_index=} {match_seq=} {marker_left-1=} {marker_right-1=}")
                    #highlight_positions.append(match_index)
                    for position in range(len(match_seq)):
                        highlight_positions.append(match_index+position)

                print(f"{highlight_positions=}")

            if highlight_positions:
                seq_row = build_seq_row(highlight_positions)
                row_as_str = render_block_to_string(
                    'includes/seq-row.html',
                    'block1', 
                    {"nucleotide": nucleotide, "seq_row": seq_row}
                )
                rows_as_strs.append(row_as_str)

        elif (row_index < 5) or (row_index > (row_count - 6)):
            #if row_index > num_of_columns and row_index < (row_count - 6):
            #    print("sep")
            #    row_as_str = '<tr class="border-b dark:border-gray-700">sep</tr>'
            #else:
            seq_row = build_seq_row()
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
