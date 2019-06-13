"""utils/writer.py Used to write GB files in output directory."""

# from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

from typing import List

import os
import errno


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
                    print(_dir_)
                    os.makedirs(os.path.dirname(_dir_))
                    self.write(_dir_, record)
                    count = count + 1
                except OSError as exc:  # Guard against race condition
                    if exc.errno != errno.EEXIST:
                        raise('Unable to create file directory')

    def write(self, dir, record):
        """Write GB files to file system."""
        with open(dir, 'w+') as handle:
            SeqIO.write(record, handle, 'genbank')
