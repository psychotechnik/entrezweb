from django.urls import path

from eweb.nucleotide.views import (
    seq_table,
    download_nucleotide,
    download_nucleotide_progress,
)

urlpatterns = [
    path("download-nucleotide/", download_nucleotide, name="download-nucleotide"),
    path("download-nucleotide/<str:task_id>/<str:seq_id>/", download_nucleotide_progress, name="download-nucleotide-progress"),
    path("seq-table/<str:seq_id>/", seq_table, name="seq-table"),
]
