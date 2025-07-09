import textwrap

import pydantic
from typing_extensions import Self


chars_per_part = 10
num_of_columns = 5
span_red = '<span class="text-red-600">'
span_close = '</span>'

#seq_parts: list[str] = textwrap.TextWrapper(width=chars_per_part).\
#    wrap(text=nucleotide.seq)

seq_parts = lambda seq: textwrap.TextWrapper(width=chars_per_part).\
    wrap(text=seq)



class SeqPart(pydantic.BaseModel):
    index_start: int
    index_end: int
    seq: str
    seq_markup: list[str]

    @pydantic.model_validator(mode='after')
    def check_index_values(self) -> Self:
        if self.index_start >= self.index_end:
            raise ValueError('end index should be greater than start index')
        return self

    @pydantic.model_validator(mode='after')
    def check_seq_length(self) -> Self:
        if len(self.seq) != ((self.index_end - self.index_start) + 1):
            raise ValueError(
                f'sequence length [{len(self.seq)}] != {(self.index_end - self.index_start) + 1}'
            )
        return self


class SeqRow(pydantic.BaseModel):
    marker_left: int
    marker_right: int
    row_seq_parts: list[SeqPart]

    @pydantic.model_validator(mode='after')
    def check_markers(self) -> Self:
        if self.marker_left != self.row_seq_parts[0].index_start:
            raise ValueError(
                f'left marker {self.marker_left} does not equal to the first row start index [{self.row_seq_parts[0].index_start}]'
            )
        if self.marker_right != self.row_seq_parts[-1].index_end:
            raise ValueError(
                f'right marker [{self.marker_left}] does not equal to the last row end index [{self.row_seq_parts[-1].index_end}]'
            )
        return self


def build_seq_row(
    parts: list[str],
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

    row_str = "".join(parts[parts_index_start:parts_index_end])
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
        col = parts[seq_index // chars_per_part + c]
        print(f"{markup_start_index=} {markup_end_index}")
        seq_markup = seq_markup_values[markup_start_index:markup_end_index]
        print(f"{seq_markup=}")

        seq_part = SeqPart(
            index_start=index_start,
            index_end=index_end,
            seq=col,
            seq_markup=seq_markup,
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
    seq_row = SeqRow(
        marker_left=marker_left,
        marker_right=marker_right,
        row_seq_parts=row_seq_parts
    )
    print(f"{seq_row.marker_left=} {seq_row.marker_right=}")
    for seq_row_part in seq_row.row_seq_parts:
        print(seq_row_part)

    return seq_row


