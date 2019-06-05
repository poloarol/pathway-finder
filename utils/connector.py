"""utils/connector.py Connects to Genbank DB."""

from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from dataclasses import dataclass
from typing import List

import urllib

Entrez.email = 'adjon081@uottawa.ca'


@dataclass
class BioConnect:
    """Provides methods to download and run blast on the NCBI servers."""

    expect: int
    hit: int
    db = 'nucleotide'
    fileFormat = 'XML'
    flavour = 'tblastn'
    service = 'psi'
    blastdb = 'nr'

    def load(self, accession):
        """Download Genbank record from NCBI-Genbank."""
        try:
            handle = Entrez.efetch(db=self.db, id=accession, rettype="gbwithparts", retmode='txt')
            return handle
        except urllib.error.HTTPError as error:
            print(error.read())

    def bioBlast(self, seq) -> List:
        """Run blast on the NCBI servers."""
        handle = NCBIWWW.qblast(self.flavour, self.blastdb, seq, self.fileFormat, self.hit, self.expect, self.service)  # noqa
        record = NCBIXML.parse(handle)
        numList = list()

        for rec in record:
            for align in rec.alignments:
                numList.append(align.accession)

        return numList
