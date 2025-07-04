from django.urls import path

from eweb.nucleotide.views import (
    seq_table,
)

urlpatterns = [
    path("<str:uid>/", seq_table, name="seq-table"),
]
