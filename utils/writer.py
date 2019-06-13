"""utils/writer.py Used to write GB files in output directory."""

# from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

from typing import List

import os


class Writer:
    """Provides a methods to write GB files."""

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

    def write(self) -> None:
        """Write GB files to file system."""
        for records in self._records:
            count = 0
            name = 'expl' + str(count) + '.gb'
            _dir_ = os.path.join(os.path.expanduser('~'),'/pathway-finder/output', name)
            with open(_dir_, 'w+') as handle:
                for record in records:
                    SeqIO.write(record, handle, 'genbank')
            count = count + 1
