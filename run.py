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
    gb = BioConnect.load(gbfile)
    rb = ReadFile(gbfile)
    genome = rb.getGenome()
    output = BioConnect.bioBlast(genome.getCore())
