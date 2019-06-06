"""Entry point of program."""

from reader import ReadGB
from utils import connector

BioConnect = connector.BioConnect


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
    main('CP013839.1', "MGAS23530_0009", 5000)

# find small genome and build unittest based on that
# Mycoplasma genitalium
# Vibrio cholerae
# Escherichia coli
# Bacilus subtilis
