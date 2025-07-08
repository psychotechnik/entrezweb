from django import forms


class DownloadNucleotideForm(forms.Form):
    seq_id = forms.CharField(label="UID", max_length=100)


class SearchForm(forms.Form):
    seq_search_query = forms.CharField(label="Query", max_length=100)



