"""Definition of the various Data Structures used in the program. """
from dataclasses import dataclass
from dataclasses import field
from typing import Tuple
from typing import NewType
from typing import Dict
from typing import List
from itertools import cycle
from itertools import islice
from Levenshtein import ratio
# from collections import namedtuple


@dataclass
class Organism:
    """Store organism information."""

    sciName: str
    accession: str

    def info(self) -> Tuple:
        """Provide infomation pertaining to indentifying organism."""
        return (self.sciName, self.accession)


@dataclass
class Gene:
    """Allow to capture the identifiers of a gene."""

    gene: str
    locus: str
    product: str
    prot_id: str
    trans: str
    location: Tuple
    strand: int

    def keys(self) -> Tuple:
        """Obtain all keys."""
        # f = namedtuple('Info', 'Gene, Locus, Product, Protein')  # noqa
        return (self.gene, self.locus, self.product, self.prot_id)

    def values(self) -> Tuple:
        """Obtain all genomic information."""
        f = (self.gene, self.locus, self.product, self.prot_id, self.trans, self.location, self.strand)  # noqa
        return f


GENE = NewType('GENE', Gene)


@dataclass
class Genome:
    """Allow to simulate bacterial genome, using a dictionary."""

    GENOME: Dict[Tuple, GENE] = field(default_factory=dict)
    core: str = ''

    def addGene(self, gene):
        """Add  a new to the genome."""
        self.GENOME[gene.keys()] = gene.values()

    def findGene(self, ident: str) -> List:
        """Find a gene by its identifier."""
        geneList = list()
        for keys in self.GENOME:
            if ident in keys:
                geneList.append(keys)
        self.setCore(geneList[0])
        return geneList

    def findCoreGeneBySimilarity(self, seq: str, similarity: float) -> List:
        """Determine the core gene in blast output by percentage similarity of seq compared."""
        geneList: List = list()
        for keys in self.GENOME:
            print(self.GENOME[keys][4])
            # if val >= similarity:
            #     geneList.append(self.GENOME[keys])
            return geneList


    def setCore(self, ident: Tuple) -> None:
        self.core = ident

    def getCore(self) -> Tuple:
        return self.GENOME[self.core][4]

    def build(self, ident: str, bp: int) -> List:
        """Build a genomic pathway based on identifier."""
        geneList: List = self.findGene(ident)
        if geneList:
            right: List = self.rbuild(geneList[0], len(geneList), bp)
            left: List = self.lbuild(geneList[0], len(geneList), bp)
            return [a + b for a, b in zip(right, left)]
        return geneList

    def rbuild(self, value: set, size: int, bp: int) -> List:
        """Build the right genomic pathway."""
        right: List = list()
        keys: List = list(self.GENOME)
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
        """Build the left genomic pathway."""
        left: List = list()
        keys: List = list(self.GENOME)
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
        """Navigate via the genome in a cycle."""
        path = list()
        length = 0

        while length < bp:
            gene: GENE = next(start)
            info: Tuple = self.GENOME[gene]
            size: int = int(info[5][1]) - int(info[5][0])  # calculate the size of the gene  # noqa
            length = length + size
            path.append(info)

        return path
