import textwrap
from dataclasses import dataclass

from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.http import require_GET
from render_block import render_block_to_string

from eweb.nucleotide.models import Nucleotide


@dataclass
class SeqRow:
    start: int
    end: int
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
    start, end = 0, 5
    for row_index in range(0, len(seq_parts) // 5):
        #print(row, seq_parts[start:end])
        seq_row = SeqRow(
            start,
            end, 
            seq_parts[start],
            seq_parts[start+1],
            seq_parts[start+2],
            seq_parts[start+4],
            seq_parts[start+3],
        )
        start+=5
        end+=5

        row_as_str = render_block_to_string(
                'includes/seq-row.html',
                'block1', 
                {"nucleotide": nucleotide, "seq_row": seq_row}
            )
        rows_as_strs.append(row_as_str)

    block_as_string = render_block_to_string(
        'includes/seq-table.html',
        'block1',
        {"nucleotide": nucleotide, "rows_as_strs": rows_as_strs}
    )
    print(block_as_string)
    #import ipdb;ipdb.set_trace()
    return HttpResponse(block_as_string)


    #return render(request, 'index.html')
