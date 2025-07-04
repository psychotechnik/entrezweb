from ninja import NinjaAPI

from eweb.nucleotide.api import router as nucleotide_router

api = NinjaAPI()

api.add_router("/nucleotides/", nucleotide_router)
