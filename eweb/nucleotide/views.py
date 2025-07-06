import textwrap
from dataclasses import dataclass

from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.http import require_GET
from render_block import render_block_to_string

from eweb.nucleotide.models import Nucleotide


@dataclass
class SeqRow:
    marker_left: int
    marker_right: int
    row1: str
    row2: str
    row3: str
    row4: str
    row5: str


@require_GET
def index(request):
    return render(request, 'index.html')


@require_GET
def seq_table(request, uid: str):

    nucleotide = get_object_or_404(Nucleotide, entrez_id=uid)
    seq_parts: list[str] = textwrap.TextWrapper(width=10).\
        wrap(text=nucleotide.seq)
    rows_as_strs = []
    marker_left, marker_right = 1, 50
    row_count = len(seq_parts) // 5
    for row_index in range(0, row_count):
        if  (row_index < 11): # or (row_index > (row_count - 10)):
            seq_index = marker_left - 1
            print(f"{row_index=} {marker_left=} {marker_right=} {row_count=}")
            seq_row = SeqRow(
                marker_left,
                marker_right, 
                seq_parts[seq_index // 10],
                seq_parts[seq_index // 10+1],
                seq_parts[seq_index // 10+2],
                seq_parts[seq_index // 10+3],
                seq_parts[seq_index // 10+4],
            )
            #if row_index > 10 and row_index < (row_count - 10):
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

    block_as_string = render_block_to_string(
        'includes/seq-table.html',
        'block1',
        {"nucleotide": nucleotide, "rows_as_strs": rows_as_strs}
    )
    #print(block_as_string)
    #import ipdb;ipdb.set_trace()
    return HttpResponse(block_as_string)


    #return render(request, 'index.html')
