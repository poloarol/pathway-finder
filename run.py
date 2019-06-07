"""Entry point of program."""

from reader import ReadGB
from utils import connector

from typing import List

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
    pathway = genome.build(coreGene, bp)
    coregene = genome.getCore()
    output = bconnect.bioBlast(coregene)
    pathways = repProcedure(output, bp, coregene, similarity)
    pathway.append(pathways)
    print(pathway)


def repProcedure(items: List, bp: int, coreGene: str, similarity: int) -> List:
    bpathway: List = list()
    for item in items:
        bconnect = BioConnect(10, 100)
        gbfile = bconnect.load(item)
        rb = ReadGB(gbfile)
        genome = rb.readfile()
        coreGeneList = genome.findCoreGeneBySimilarity(coreGene, similarity)
        for i in coreGeneList:
            genome.setCore(i)
            pathway = genome.build(coreGene, bp)
            bpathway.append(pathway)
    return bpathway

if __name__ == '__main__':
    main('CP013839.1', "MGAS23530_0009", 5000, 0.75)

# find small genome and build unittest based on that
# Mycoplasma genitalium
# Vibrio cholerae
# Escherichia coli
# Bacilus subtilis
