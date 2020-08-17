"""Read and Write GB files."""

from .structures import Genome
from .structures import Organism
from .structures import Gene
from .structures import GENE

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

from typing import List
from typing import Dict
from typing import Tuple


import os
import errno
import textwrap


class ReadGB:
    """
    Read GB file and extract information

    ...

    Attributes
    ----------
    flag : bool
        Used to detect error in parsing Gb file
    _record : File Object
        Result of the SeqIO.parse() from biopython
    GENOME : Genome()
        Organisms genome to store genes as key-value pairs

    Methods
    -------
    provideOrg()
        Gives information about the organism
    readfile()
        Parses the GB file, extract useful features and builds genome

    """

    def __init__(self, genbank):
        """Initialize the reader and store gb."""

        self.flag = False

        try:
            self._record = SeqIO.read(genbank, 'genbank')
        except Exception:
            self.flag = True
        self.GENOME = Genome()
        try:
            self._org = Organism(self._record.annotations['organism'], self._record.annotations['accessions'][0])  # noqa
        except KeyError:
            raise KeyError('Organim and Accession keys are nor present.')

    def provideOrg(self):
        """Provide indentifier about organism created."""
        return self._org.info()

    def readfile(self):
        """Read the gb file and return a new Genome."""
        if not self.flag:
            try:
                prot_id: str = 'N/A'
                locus: str = 'N/A'
                product: str = 'N/A'
                gene: str = 'N/A'
                translation: str = 'N/A'
                description: str = 'N/A'
                for i, recs in enumerate(self._record.features):
                    if recs.type == 'CDS':
                        if 'translation' in recs.qualifiers:
                            translation = recs.qualifiers['translation'][0]
                        if 'gene' in recs.qualifiers:
                            gene = recs.qualifiers['gene'][0]
                        if 'locus_tag' in recs.qualifiers:
                            locus = recs.qualifiers['locus_tag'][0]
                        if 'product' in recs.qualifiers:
                            product = recs.qualifiers['product'][0]
                        if 'protein_id' in recs.qualifiers:
                            prot_id = recs.qualifiers['protein_id'][0]
                        if 'description' in recs.qualifiers:
                            description = recs.qualifiers['description'][0]
                        strand: int = int(recs.strand)
                        location: Tuple = (recs.location.start.position, recs.location.end.position)  # noqa
                        g: GENE = Gene(gene, locus, product, prot_id, translation, description, location, strand)  # noqa
                        self.GENOME.addGene(g)
            except KeyError:
                raise KeyError('CDS key not found.')
            return self.GENOME


class Writer:
    """
    Provides a methods to write GB files

    Attributes
    ..........

    _genes (list) : List of all of all inputed genes
    _records (list) : Stores all gb records

    Methods
    .......

    produce_records : Takes a gene info (tuple) and produces
                      a GB record
    parse : invokes the produce_records method and produces a list of
            Gb records
    write : writes a GB record to the file system
    writeGB : provided with a list of GB records, invokes write

    """

    def __init__(self, genes: List):
        """Initialize variables"""
        self._genes: List = genes

    def parse(self):
        """Use the provided list of subsectioned genes, it provides a ist or SeqIO records."""  # noqa
        sequences: str = ""
        for genes in self._genes:
            record = list()
            for gene in genes:
                translation = '\n'.join(textwrap.wrap(gene.trans, 60))
                line: str = "{0} - {1}\n{2}".format(gene.prot_id, gene.protein, translation)
                print(line)
            #     value = self.produce_records(gene)
            #     record.append(value)
            # self._records.append(record)
