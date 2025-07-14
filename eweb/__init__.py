import textwrap

chars_per_part = 10
num_of_columns = 5
top_rows_num = 25
buttom_rows_num = 25
span_red = '<span class="text-red-600">'
span_close = '</span>'

seq_parts = lambda seq: textwrap.TextWrapper(
    width=chars_per_part
).wrap(text=seq)

