"""Entry point of program."""

from utils import connector
from utils import reader

from typing import List

import time
import sys

BioConnect = connector.BioConnect
bconnect = None
ReadGB = reader.ReadGB


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
        self.readGB = ReadGB(self._gb)

    def getGenome(self):
        """Obtain the genome built from passed gb file."""
        self.GENOME = self.readGB.readfile()


def main(gbfile, coreGene, bp, similarity):
    bconnect = BioConnect(10, 100)
    gbfile = bconnect.load(gbfile)
    rb = ReadGB(gbfile)
    genome = rb.readfile()
    if not genome:
        sys.exit()
    pathways = genome.build(coreGene, bp)
    coregene = genome.getCore()
    output = bconnect.bioBlast(coregene)
    bpathways = repProcedure(output, bp, coregene, similarity)
    pathways.append(bpathways)
    # count = 0
    # for pathway in pathways:
    #     if count < 5:
    #         print("===============================================")
    #         result = [x if isinstance(x[0], tuple) else [x] for x in pathway]
    #         print(result)
    #         print("===============================================")
    #         count = count + 1


def repProcedure(items: List, bp: int, coreGene: str, similarity: int) -> List:
    bpathway: List = list()
    counter: int = 0
    for item in items:
        bconnect = BioConnect(10, 100)
        gbfile = bconnect.load(item)
        rb = ReadGB(gbfile)
        genome = rb.readfile()
        if not genome:
            sys.exit()
        genes: List = genome.findCoreGeneBySimilarity(coreGene, similarity)
        if genes:
            for gene in genes:
                genome.setCore(gene)
                bpathway.append(genome.buildsimilarity(gene, bp))
        counter = counter + 1
        if(counter % 3 == 0):  # Used because of NCBI's policy on requests without API key. with API key, change to 10  # noqa
            time.sleep(3)
    return bpathway


if __name__ == '__main__':
    main('CP013839.1', "MGAS23530_0009", 5000, 0.75)

# find small genome and build unittest based on that
# Mycoplasma genitalium
# Vibrio cholerae
# Escherichia coli
# Bacilus subtilis
