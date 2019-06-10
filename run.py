"""Entry point of program."""

from reader import ReadGB
from utils import connector

from typing import List
from itertools import chain

import time

BioConnect = connector.BioConnect
bconnect = None

class ReadFile():
    """Read and Process GB file."""

    def __init__(self, gbfile):
        """Initialize the reader."""
        self._gb = gbfile
        self.readGB = ReadGB(self._gb)
        # self.GENOME: Dict = dict()

    def getGenome(self):
        """Obtain the genome built from passed gb file."""
        self.GENOME = self.readGB.readfile()


def main(gbfile, coreGene, bp, similarity):
    bconnect = BioConnect(10, 100)
    gbfile = bconnect.load(gbfile)
    rb = ReadGB(gbfile)
    genome = rb.readfile()
    pathways = genome.build(coreGene, bp)
    coregene = genome.getCore()
    output = bconnect.bioBlast(coregene)
    bpathways = repProcedure(output, bp, coregene, similarity)
    pathways.append(bpathways)
    pathways = list(chain.from_iterable(pathways))
    for pathway in pathways:
        print("==============================================================================")
        print(pathway)
        print("==============================================================================")


def repProcedure(items: List, bp: int, coreGene: str, similarity: int) -> List:
    bpathway: List = list()
    counter: int = 0
    for item in items:
        bconnect = BioConnect(10, 100)
        gbfile = bconnect.load(item)
        rb = ReadGB(gbfile)
        genome = rb.readfile()
        genes: List = genome.findCoreGeneBySimilarity(coreGene, similarity)
        if genes:
            for gene in genes:
                genome.setCore(gene)
                bpathway.append(genome.buildsimilarity(gene, bp))
        counter = counter + 1
        if(counter % 3 == 0):  # Used because of NCBI's policy on requests without API key. with API key, change to 10
            time.sleep(3)
    return bpathway

if __name__ == '__main__':
    main('CP013839.1', "MGAS23530_0009", 5000, 0.75)

# find small genome and build unittest based on that
# Mycoplasma genitalium
# Vibrio cholerae
# Escherichia coli
# Bacilus subtilis
