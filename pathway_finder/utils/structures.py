"""Definition of the various Data Structures used in the program. """
from dataclasses import dataclass
from dataclasses import field
from typing import Tuple, NamedTuple
from typing import NewType
from typing import Dict
from typing import List
from itertools import cycle
from itertools import islice
from collections import namedtuple

import Levenshtein


@dataclass
class Organism:
    """
    Store organism information

    ...

    Attributes
    ----------
    sciName: str
        Scientific name of the organism
    accession: str
        Accession number of the organims

    Methods
    -------
    info()
        Provides information about the organism in question
    """

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
    desc: str
    location: Tuple
    strand: int

    def keys(self) -> Tuple:
        """Obtain all keys."""
        # f = namedtuple('Info', 'Gene, Locus, Product, Protein')  # noqa
        return (self.gene, self.locus, self.product, self.prot_id)

    def values(self) -> NamedTuple:
        """Obtain all genomic information."""
        Gene: NamedTuple = namedtuple('Gene', 'gene locus product prot_id trans desc loc strand')  # noqa
        f = Gene(gene=self.gene, locus=self.locus, product=self.product, prot_id=self.prot_id, trans=self.trans, desc=self.desc,loc=self.location, strand=self.strand)  # noqa
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
        genes = list()
        for keys in self.GENOME:
            if ident in keys:
                genes.append(keys)
        return genes

    def findCoreGeneBySimilarity(self, seq: str, similarity: float):
        """Determine the core gene in blast output by percentage similarity of seq compared."""  # noqa
        genes: List = list()
        for keys in self.GENOME:
            val = Levenshtein.ratio(self.GENOME[keys][4], seq)
            if val >= similarity:
                genes.append(keys)
        return genes

    def setCore(self, ident: Tuple) -> None:
        """Public method to set Genome core gene."""
        self.__setCore__(ident)

    def __setCore__(self, ident: Tuple) -> None:
        """Private method to set Genome core gene."""
        self.core = ident

    def getCore(self) -> str:
        """Provides the core genome's core gene"""
        return self.GENOME[self.core][4]

    def build(self, ident: str, bp: int) -> List:
        """Build a genomic pathway based on identifier."""
        genes: List = self.findGene(ident)
        self.setCore(genes[0])
        if genes:
            right: List = self.rbuild(genes[0], len(genes), bp)
            left: List = self.lbuild(genes[0], len(genes), bp)
            right.extend(left)
            return right
        return genes

    def rbuild(self, value: set, size: int, bp: int) -> List:
        """Build the right genomic pathway."""
        right: List = list()
        keys: List = list(self.GENOME)
        indices = keys.index(value) if size == 1 else [i for i, x in enumerate(keys) if x == value]  # noqa
        # TODO: To speed it up, NumPy can be used
        # (https://stackoverflow.com/questions/6294179/how-to-find-all-occurrences-of-an-element-in-a-list)
        # as explained
        kcycle = cycle(keys)  # itertool.cycle, to cycle over the list until desired length or queried gene is met again.  # noqa

        if isinstance(indices, list):
            for i in indices:
                start = islice(kcycle, indices, None)
                right.extend(self.paths(start, bp))
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
            for i in indices:
                start = islice(kcycle, indices, None)
                next(start)
                left.extend(self.paths(start, bp))
        else:
            start = islice(kcycle, indices, None)
            next(start)
            left = self.paths(start, bp)
        return left

    def buildsimilarity(self, value: set, bp: int):
        """After setting core gene by similarity, use this to build pathway."""
        right: List = self.rbuild(value, 1, bp)
        left: List = self.lbuild(value, 1, bp)
        right.extend(left)
        return right

    def paths(self, start, bp) -> List:
        """Navigate via the genome in a cycle."""
        path = list()
        length = 0

        while length <= bp:
            gene: GENE = next(start)
            info: Tuple = self.GENOME[gene]
            size: int = int(info.loc[1]) - int(info.loc[0])  # calculate the size of the gene  # noqa
            length = length + size
            path.append(info)

        return path
