"""Entry point of program."""

from utils import connector
from utils import processor

from typing import List
from itertools import chain
from dataclasses import dataclass
from typing import List

import time
import sys


BioConnect = connector.BioConnect
Writer = processor.Writer
ReadGB = processor.ReadGB


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


@dataclass
class Finder:
    """
    Provides an interface to call all methods necessary to subset a Gb file
    and search for identical looking pathways as the looked upon one.

    Attributes
    ..........

    accession (str) : Accession number (identifer) for GB file
    coreGene (str) : Gene of interest to build pathway around
    expect (int) = 10 : Expected value of blast output
    hit (int) = 100 : size of blast output i.e. number of organisms
    bp (int) = 2500 : number of base pairs to look for around core gene for
    pathway building
    similarity (float) = 0.75 : allows to find identical genes to core gene in
    other genomes

    Methods
    .......

    finder : Finds identical pathways using the E-Utils and Blast API
    (appropriate functions in code)
    flatten (list) : Flattens a list by removing nested list within it
    repProcedure (list, bp, core gene, similarity) : Parses the blast output
    and finds idetical genes within it
    produce (list) : Produces a GB files of identical files from results of
    finder

    """

    accession: str
    coreGene: str
    expect: int = 10
    hit: int = 100
    bp: int = 2500
    similarity: float = 0.75
    bconnect = BioConnect(expect, hit)

    def finder(self) -> List:
        """ Calls other classes to generate a directory of similar pathways as the queried one."""
        gb = self.bconnect.load(self.accession)
        rb = ReadGB(gb)
        genome = rb.readfile()

        if not genome:
            sys.exit()

        pathway: List = genome.build(self.coreGene, self.bp)
        pathway = self.flatten(pathway)
        coregene: str = genome.getCore()
        output: List = self.bconnect.bioBlast(coregene)
        bpathways: List = self.repProcedure(output, self.bp, coregene, self.similarity)  # noqa
        bpathways.append(pathway)

        return bpathways

    def produce(self, pathways: List) -> None:
        """Create a gb files of similar pathways compared to queried gene."""
        writer = Writer(pathways)
        writer.parse()
        writer.writeGB()

    def repProcedure(self, items: List, bp: int, coreGene: str, similarity: float) -> List:  # noqa
        bpathway: List = list()
        counter: int = 0
        for item in items:
            bconnect = BioConnect(self.expect, self.hit)
            gbfile = bconnect.load(item)
            rb = ReadGB(gbfile)
            genome = rb.readfile()

            if not genome:
                sys.exit()

            genes: List = genome.findCoreGeneBySimilarity(coreGene, similarity)
            if genes:
                for gene in genes:
                    genome.setCore(gene)
                    path: List = self.flatten(genome.buildsimilarity(gene, bp))
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
