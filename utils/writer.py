"""utils/writer.py Used to write GB files in output directory."""

# from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

from typing import List


class Writer:
    """Provides a methods to write GB files."""

    def __init__(self, genes: List):
        self._genes: List = genes

    def parse(self) -> List:
        """
        Using the provided list of subsectioned genes, it provides a
        list or SeqIO records.
        """
        records = list()
        for genes in self._genes:
            for gene in genes:
                records.append(self.produce_records(gene))

        return records

    def produce_records(self, genes):
        """Produces a record of a gene in Gb format."""
        seq = Seq(gene[4], IUPAC.protein)
        record = SeqRecord(seq, id=gene[3], name=gene[0], description='N/A', annotations=[gene[1], gene[2]])  # noqa
        loc = SeqFeature(FeatureLocation(gene[5][0], gene[5][1]), type='CDS', strand=gene[6])  # noqa
        record.features[loc]
        return record
