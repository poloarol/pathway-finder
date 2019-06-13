"""Entry point of program."""

from utils import connector
from utils import reader
from utils import writer

from typing import List
from itertools import chain
from dataclasses import dataclass
from typing import List

import time
import sys


class ReadFile():
    """
    Read and Process GB file

    ...

    Attributes
    ----------
    _gb : str
        accession number of GB file, to download file

    readGB : Object
        Provides methods to read  GB files

    Methods
    -------
    getGenome():
        provides the complete genome of organism
    """

    def __init__(self, gbfile):
        """Initialize the reader."""
        self._gb = gbfile
        self.readGB = reader.ReadGB(self._gb)

    def getGenome(self):
        """Obtain the genome built from passed gb file."""
        self.GENOME = self.readGB.readfile()


@dataclass
class Finder:

    accession: str
    expect: int = 10
    hit: int = 100
    coreGene: str
    bp: int = 2500
    similarity: float = 0.75
    ReadGB = reader.ReadGB
    Writer = writer.Writer
    bconnect = connector.BioConnect(expect, hit)

    def finder(self) -> List:
        """ Calls other classes to generate a directory of similar pathways as the queried one."""
        gb = self.bconnect.load(self.accession)
        rb = self.ReadGB(gb)
        genome = rb.readfile()

        if not genome:
            sys.exit()

        pathway = genome.build(self.coreGene, self.bp)
        pathways = self.flatten(pathway)
        coregene = genome.getCore()
        output = self.bconnect.bioBlast(coregene)
        bpathways = self.repProcedure(output, self.bp, coregene, self.similarity)  # noqa
        bpathways.append(pathways)

        return bpathways

    def produce(self, pathways: List) -> None:
        """Create a gb files of similar pathways compared to queried gene."""
        writer = self.Writer(pathways)
        writer.parse()
        writer.writeGB()

    def repProcedure(self, items: List, bp: int, coreGene: str, similarity: float) -> List:  # noqa
        bpathway: List = list()
        counter: int = 0
        for item in items:
            bconnect = connector.BioConnect(self.expect, self.hit)
            gbfile = bconnect.load(item)
            rb = self.ReadGB(gbfile)
            genome = rb.readfile()

            if not genome:
                sys.exit()

            genes: List = genome.findCoreGeneBySimilarity(coreGene, similarity)
            if genes:
                for gene in genes:
                    genome.setCore(gene)
                    path = self.flatten(genome.buildsimilarity(gene, bp))
                    bpathway.append(path)
            counter = counter + 1
            if(counter % 3 == 0):  # Used because of NCBI's policy on requests without API key. with API key, change to 10  # noqa
                time.sleep(3)
        return bpathway

    def flatten(self, path: List) -> List:
        # TODO: This is just hack, I'll need to work on
        # the data structure output to remove this extra step
        """ Flattens the list i.e. removes nested list in output"""
        return [*chain.from_iterable(x if isinstance(x[0], tuple) else [x] for x in path)]  # noqa


# finder = Finder(accession='CP013839.1', hit=10, expect=100, coreGene="MGAS23530_0009", bp=5000, similarity=0.75)
# finder = Finder(accession='CP013839.1', coreGene="MGAS23530_0009")

# find small genome and build unittest based on that
# Mycoplasma genitalium
# Vibrio cholerae
# Escherichia coli
# Bacilus subtilis
