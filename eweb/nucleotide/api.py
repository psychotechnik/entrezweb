from ninja import Router

from .models import Nucleotide

router = Router()

@router.get('/')
def list_nucleotides(request):
    return [
        {
            "entrez_id": n.entrez_id,
            "name": n.name,
            "description": n.description
        }
        for n in Nucleotide.objects.all()
    ]

@router.get('/{nucleotide_id}')
def nucleotide_details(request, nucleotide_id: int):
    nucleotide = Nucleotide.objects.get(id=nucleotide_id)
    return {
        "entrez_id": nucleotide.entrez_id,
        "name": nucleotide.name,
        "description": nucleotide.description,
        "sequence": nucleotide.seq,
    }
