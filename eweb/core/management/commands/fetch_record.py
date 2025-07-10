from Bio import Entrez, SeqIO

from django.core.management.base import BaseCommand
from django.core.files.base import ContentFile
from django.conf import settings

from eweb.nucleotide.models import Nucleotide

"""
    all_seq_ids = [
        "3008860385",
        "3008860383",
        "3008860381",
        "3008860379",
        "3008860377",
        "3008860375",
        "3008860373",
        "3008860371",
        "3008860369",
        "3008860367",
        "3008860365",
        "3008860363",
        "3008860361",
        "3008860359",
        "3008713457",
        "3008713446",
        "3008712886",
        "3008712874",
        "3008711791",
    ]

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

        #handle = Entrez.esearch(db=db_type, term=search_term, retmax=max_results)
        #rec_list = Entrez.read(handle)
        #handle.close()

        #stream = Entrez.efetch(db=dbToSearch, id=ids, rettype=returnType)
        #records = list(SeqIO.parse(stream, returnType)) # Seq object, can treat like string - See chapter 3 - https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec17
        #stream.close()

        Entrez.email = settings.ADMINS[0][1]

        #  title:  NC_004718.3 SARS coronavirus Tor2, complete genome
        #  extra:  gi|30271926|ref|NC_004718.3|[30271926]
        #  # NC_004718.3
        #
        return_type = "gb" # gb fasta
        return_mode = "text" # text fasta

        if 0: 
            search_term = "coronavirus"
            max_results = 10
            with Entrez.esearch(db="nucleotide", term=search_term, retmax=max_results) as handle:
                rec_list = Entrez.read(handle)
                for rec in rec_list:
                    print(rec)
            ids = rec_list['IdList']
            print(ids)
            with Entrez.efetch(
                db="nucleotide",
                rettype=return_type,
                retmode=return_mode,
                id=ids[0],
            ) as handle:
                seq_record = SeqIO.read(handle, "fasta")

        if 0:
            id = "224589800"
            with Entrez.esummary(db="nucleotide", id=id) as handle:
                record = Entrez.read(handle)
                import ipdb;ipdb.set_trace()

            #record[0]['Length'].numerator -> 29751
            #IntegerElement(29751, attributes={})

        if 1: 
            #https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch?db=nucleotide&rettype=fasta&retmode=text&id=30271926
            #id = "3008860387"
            seq_id = "30271926"
            #id = "224589800"
            with Entrez.efetch(
                db="nucleotide", rettype="fasta", retmode="text", id=seq_id
            ) as handle:
                seq_record = SeqIO.read(handle, "fasta")

            try:
                n = Nucleotide.objects.get(entrez_id=seq_id)
            except Nucleotide.DoesNotExist:
                n = Nucleotide.objects.create(
                    entrez_id=seq_id,
                    name=seq_record.name,
                    description=seq_record.description,
                )
                n.fasta_file.save(
                    f"{seq_id}.fasta",
                    ContentFile(str(seq_record.format("fasta"))),
                )
                import ipdb;ipdb.set_trace()


        """
        n, created = Nucleotide.objects.get_or_create(
            entrez_id=id,
            name=seq_record.name,
            description=seq_record.description,
            seq=str(seq_record.seq)
        )
        if created: 
            print(f"created nucleotide {n.entrez_id}")
        """

        """
        stream = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="xml")
        recs = Entrez.read(stream)
        for r in recs:
            #for k, v in r.items():
            #    print(k, v)
            #print(r["TSeq_sequence"])
            n, created = Nucleotide.objects.get_or_create(
                seq=r["TSeq_sequence"],
            )
            if created: 
                print(f"created nucleotide id: {id}")
        stream.close()
        """
