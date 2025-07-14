from django import forms
from django.core.exceptions import ValidationError


class DownloadNucleotideForm(forms.Form):
    seq_id = forms.CharField(label="UID", max_length=100,)


class SearchForm(forms.Form):
    seq_search_query = forms.CharField(label="Query", max_length=100, required=False)

    def clean_seq_search_query(self):
        data = self.cleaned_data["seq_search_query"]
        if data.strip() and not data.strip().isalpha():
            raise ValidationError("search query contains non-alphabetical characters")

        return data



