from Bio import Entrez, SeqIO

from django.core.management.base import BaseCommand

from eweb.nucleotide.models import Nucleotide
from eweb.nucleotide.tasks import download_nucleotide_task

"""

<Id>30271926</Id>
<Item Name="Caption" Type="String">NC_004718</Item>
<Item Name="Title" Type="String">SARS coronavirus Tor2, complete genome</Item>
<Item Name="Extra" Type="String">gi|30271926|ref|NC_004718.3|[30271926]</Item>
<Item Name="Gi" Type="Integer">30271926</Item>
<Item Name="CreateDate" Type="String">2003/04/14</Item>
<Item Name="UpdateDate" Type="String">2020/11/20</Item>
<Item Name="Flags" Type="Integer">512</Item>
<Item Name="TaxId" Type="Integer">227984</Item>
<Item Name="Length" Type="Integer">29751</Item>
<Item Name="Status" Type="String">live</Item>
<Item Name="ReplacedBy" Type="String"/>
<Item Name="Comment" Type="String">
<![CDATA[ ]]>
</Item>
<Item Name="AccessionVersion" Type="String">NC_004718.3</Item>
"""

"""
<TSeq>
<TSeq_seqtype value="nucleotide"/>
<TSeq_accver>NC_004718.3</TSeq_accver>
<TSeq_taxid>227984</TSeq_taxid>
<TSeq_orgname>SARS coronavirus Tor2</TSeq_orgname>
<TSeq_defline>SARS coronavirus Tor2, complete genome</TSeq_defline>
<TSeq_length>29751</TSeq_length>
<TSeq_sequence></TSeq_sequence>
</TSeq>
"""


class Command(BaseCommand):
    help = "fetch record from Entrez"

    #def add_arguments(self, parser):
    #    parser.add_argument('sample', nargs='+')

    def handle(self, *args, **options):

        id = "30271926"
        #id = "224589800"
        res = download_nucleotide_task.delay(id)
        res.get(timeout=1)#

