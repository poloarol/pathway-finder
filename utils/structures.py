""" Definition of the various Data Structures used in the program. """
from dataclasses import dataclass
from dataclasses import field
from typing import Tuple
from typing import NewType
from typing import Set
from typing import Dict
from typing import List
from itertools import cycle
from itertools import islice
# from collections import namedtuple


@dataclass
class Gene:
    """Allow to capture the identifiers of a gene."""

    gene: str
    locus: str
    product: str
    prot_id: str
    location: Tuple

    def keys(self) -> Set:
        # f = namedtuple('Info', 'Gene, Locus, Product, Protein')  # noqa
        return set(self.gene, self.locus, self.product, self.prot_id)

    def values(self) -> Tuple:
        f = (self.gene, self.locus, self.product, self.prot_id, self.location)
        return f


GENE = NewType('GENE', Gene)


@dataclass
class Genome:
    """Allow to simulate bacterial genome, using a dictionary."""

    GENOME: Dict[Set, GENE] = field(default_factory=dict)

    def addGene(self, gene):
        self.GENOME.update(gene.info(), gene.values())

    def findGene(self, ident: str) -> List:
        """ Allows to find a list of all genes having a certain identifier. """

        geneList = list()
        return geneList

    def build(self, ident: str, bp: int) -> List:
        """builds the genomic pathway based on identifier."""

        geneList: List = self.findGene(ident)
        if geneList:
            right: List = self.rbuild(geneList[0], len(geneList), bp)
            left: List = self.lbuild(geneList[0], len(geneList), bp)
            return [a + b for a, b in zip(right, left)]
        return geneList

    def rbuild(self, value: set, size: int, bp: int) -> List:
        """ Generates a list of genes on the right of the core gene."""
        right: List = list()
        keys: List = self.GENOME.keys()
        indices = keys.index(value) if size == 1 else [i for i, x in enumerate(keys) if x == value]  # noqa
        # To speed it up, NumPy can be used (https://stackoverflow.com/questions/6294179/how-to-find-all-occurrences-of-an-element-in-a-list)  # noqa
        # as explained
        kcycle = cycle(keys)  # itertool.cycle, to cycle over the list until desired length or queried gene is met again.  # noqa

        if isinstance(indices, list):
            for i in len(indices):
                start = islice(kcycle, indices, None)
                right.append(self.paths(start, bp))
        else:
            start = islice(kcycle, indices, None)
            right = self.paths(start, bp)
        return right

    def lbuild(self, value: set, size: int, bp: int) -> List:
        left: List = list()
        keys: List = self.GENOME.keys()
        keys.reverse()
        indices = keys.index(value) if size == 1 else [i for i, x in enumerate(keys) if x == value]  # noqa
        kcycle = cycle(keys)

        if isinstance(indices, list):
            for i in len(indices):
                start = islice(kcycle, indices, None)
                left.append(self.paths(start, bp))
        else:
            start = islice(kcycle, indices, None)
            left = self.paths(start, bp)
        return left

    def paths(self, start, bp) -> List:
        path = list()
        length = 0

        while length < bp:
            gene: GENE = next(start)
            info: Tuple = self.GENOME[gene]
            size: int = info.values()[4][1] - info.values()[4][1]  # calculate the size of the gene  # noqa
            length = length + size
            path.append(info)
