"""Entry point of program."""

from reader import ReadGB
from utils.connector import BioConnect

from typing import Dict


class ReadFile():
    """Read and Process GB file."""

    def __init__(self, gbfile):
        self._gb = gbfile
        self.readGB = ReadGB(self._gb)
        self.GENOME: Dict = dict()

    def getGenome(self):
        self.GENOME = self.readGB.readfile()


def main(gbfile):
    bconnect = BioConnect(10, 100)
    gbfile = bconnect.load(gbfile)
    rb = ReadGB(gbfile)
    rb.readfile()
    # output = BioConnect.bioBlast(genome.getCore())


if __name__ == '__main__':
    main('AL590842.1')
