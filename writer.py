"""utils/writer.py Used to write GB files in output directory."""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

from typings import List


class Writer:
    """Provides a methods to write GB files."""

    def __init__(self, genes: List):
        self._genes: List = genes

    def parse():
        """Create a Sequence and Record object."""
        pass
