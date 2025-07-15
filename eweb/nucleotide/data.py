import pydantic
from typing_extensions import Self

from eweb import (
    seq_parts,
    chars_per_part,
    num_of_columns,
    top_rows_num,
    buttom_rows_num,
    span_red,
    span_close,
    Timer,
)



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

    #@pydantic.model_validator(mode='after')
    #def check_seq_length(self) -> Self:
    #    if len(self.seq) != ((self.index_end - self.index_start) + 1):
    #        raise ValueError(
    #            f'sequence length [{len(self.seq)}] != {(self.index_end - self.index_start) + 1}'
    #        )
    #    return self


class SeqRow(pydantic.BaseModel):
    marker_left: int
    marker_right: int
    row_seq_parts: list[SeqPart]

    """
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
    """


def build_seq_row(
    parts: list[str],
    marker_left: int,
    marker_right: int,
    markup_style: str = "html",
    highlight_positions=None,
):
    #t = Timer()
    #t.start()

    if not highlight_positions:
        highlight_positions = []
    seq_markup_values = []
    row_str = "".join(parts) #[parts_index_start:parts_index_end])
    print(f"{row_str=}")
    print(parts)
    #assert len(row_str) == 50, f"row seq str len: {len(row_str)}"
    print(f"{highlight_positions=}")
    for i, val in enumerate(row_str):
        if marker_left-1+i in highlight_positions:
            seq_markup_values.append(
                f'{span_red}{val}{span_close}' if markup_style == "html" else f'[red]{val}'
            )
        else:
            seq_markup_values.append(val)
    #assert len(seq_markup_values) == 50, f"seq markup len: {len(seq_markup_values)}"
    row_seq_parts = []
    m_start_idx = 0
    m_end_idx = chars_per_part-1
    idx_start = marker_left-1
    idx_end = idx_start+chars_per_part-1
    for c in range(0, num_of_columns):
        seq_markup = seq_markup_values[m_start_idx:m_end_idx]
        #assert len(seq_part.seq_markup) == 10, f"seq markup len: {len(seq_part.seq_markup)}"
        row_seq_parts.append(
            SeqPart(
                index_start=idx_start,
                index_end=idx_end,
                seq=parts[c],
                #seq=parts[seq_index // chars_per_part + c] ,
                seq_markup=seq_markup,
            )
        )
        idx_start = idx_end+1
        idx_end = idx_start+chars_per_part-1
        m_start_idx = m_end_idx
        m_end_idx = m_start_idx + chars_per_part 
    seq_row = SeqRow(
        marker_left=marker_left,
        marker_right=marker_right,
        row_seq_parts=row_seq_parts
    )
    #t.stop()
    return seq_row


