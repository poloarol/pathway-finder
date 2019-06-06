"""Entry point of program."""

from reader import ReadGB
from utils.connector import BioConnect

from typing import Dict


class ReadFile():
    """Read and Process GB file."""

    def __init__(self, gbfile):
        self._gb = gbfile
        self.readGB = ReadGB(self._gb)
        # self.GENOME: Dict = dict()

    def getGenome(self):
        self.GENOME = self.readGB.readfile()


def main(gbfile, coreGene, bp):
    bconnect = BioConnect(10, 100)
    gbfile = bconnect.load(gbfile)
    rb = ReadGB(gbfile)
    genome = rb.readfile()
    pathway = genome.build(coreGene, bp)
    coregene = genome.getCore()
    output = bconnect.bioBlast(coregene)
    print(output)
    # print(output)


if __name__ == '__main__':
    main('AP018392.1', "SONE68_0006", 15000)

# find small genome and build unittest based on that
# Mycoplasma genitalium
# Vibrio cholerae
# Escherichia coli
# Bacilus subtilis