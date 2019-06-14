"""Read and Write GB files."""

from .structures import Genome
# from .structures import Organism
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
            # self._org = Organism(self._record.annotations['organism'], self._record.annotations['accessions'][0])  # noqa
            pass
        except KeyError:
            raise KeyError('Organim and Accession keys are nor present.')

    def provideOrg(self) -> Tuple:
        """Provide indentifier about organism created."""
        return self._org.info()

    def readfile(self) -> Dict:
        """Read the gb file and return a new Genome."""
        if not self.flag:
            try:
                prot_id: str = 'n/a'
                locus: str = 'n/a'
                product: str = 'n/a'
                gene: str = 'n/a'
                translation: str = 'n/a'
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
                        strand: int = int(recs.strand)
                        location: Tuple = (recs.location.start.position, recs.location.end.position)  # noqa
                        g: GENE = Gene(gene, locus, product, prot_id, translation, location, strand)  # noqa
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
        self._records: List = list()

    def parse(self) -> List:
        """Use the provided list of subsectioned genes, it provides a ist or SeqIO records."""
        record = list()
        for genes in self._genes:
            for gene in genes:
                record.append(self.produce_records(gene))
            self._records.append(record)

    def produce_records(self, gene):
        """Produce a record of a gene in Gb format."""
        seq = Seq(gene[4], IUPAC.protein)
        record = SeqRecord(seq, id=gene.prot_id, name=gene.gene, description='N/A', annotations={'product': gene.product, 'locus_tag': gene.locus})  # noqa
        loc = SeqFeature(FeatureLocation(gene.loc[0], gene.loc[1]), type='CDS', strand=gene.strand)  # noqa
        record.features.append(loc)
        return record

    def writeGB(self) -> None:
        """Go through records list and write records by calling write method."""
        count = 0
        for record in self._records:
            name = 'expl' + str(count) + '.gb'
            _dir_ = os.path.join(os.path.expanduser('~'), 'Desktop/output/', name)  # noqa
            if not os.path.exists(os.path.dirname(_dir_)):
                try:
                    os.makedirs(os.path.dirname(_dir_))
                except OSError as exc:  # Guard against race condition
                    if exc.errno != errno.EEXIST:
                        raise
            self.write(_dir_, record)
            count = count + 1

    def write(self, dir, record):
        """Write GB files to file system."""
        with open(dir, 'w+') as handle:
            SeqIO.write(record, handle, 'genbank')
