"""utils/connector.py Connects to Genbank DB."""

from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from dataclasses import dataclass
from typing import List

import urllib
import ssl

Entrez.email = 'adjon081@uottawa.ca'
ssl._create_default_https_context = ssl._create_unverified_context


@dataclass
class BioConnect:
    """
    Provides methods to download and run blast on the NCBI servers.

    ...

    Attributes
    ----------
    expect : int
        Expected value of blast output
    hit : int
        Hit size of blast output i.e. How many organism are expected
    db : str
        NCBI DataBase to be used. Set to Nucleotide
    fileFormat : str
        Expected format of blast output. Set to XML
    service : str
        Type of blast to be performed. Set to PSI Blast
    blastdb: str
        The blast db to be used. Set to Non-Redundant

    Methods
    -------
    load(accession=str)
        Provides access to the NCBI DB via the E-Utils params to
        obtain files
    bioBlast(seq: str)
        Performs a online psi-blast using a provided sequence

    """

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
            handle = Entrez.efetch(db=self.db, id=accession, rettype="gbwithparts", retmode='txt')  # noqa
            return handle
        except urllib.error.HTTPError as error:
            print(error.read())

    def bioBlast(self, seq) -> List:
        """Run blast on the NCBI servers."""
        handle = NCBIWWW.qblast(self.flavour, self.blastdb, seq, format_type=self.fileFormat, hitlist_size=self.hit, expect=self.expect, service=self.service)  # noqa
        record = NCBIXML.parse(handle)
        numList = list()

        for rec in record:
            for align in rec.alignments:
                numList.append(align.accession)

        return numList
