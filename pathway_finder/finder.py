"""Entry point of program."""

from utils import connector
from utils import processor

from typing import List
from itertools import chain
from dataclasses import dataclass

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

    email: str
    seq: str = None
    accession: str = None
    coreGene: str = None
    bp: int = 2500
    similarity: float = 0.5
    bconnect = None

    def __post_init__(self):
        self.bconnect = BioConnect(self.email)

    def finder(self):
        """ Calls other classes to generate a directory of similar pathways as the queried one."""  # noqa
        gb = None
        pathway = None
        coregene = None
        output = None
        pathways = None

        if self.accession:
            gb = self.bconnect.load(self.accession)
            rb = ReadGB(gb)
            genome = rb.readfile()
            organism = rb.provideOrg()
            pathway: List = genome.build(self.coreGene, self.bp)
            genome.setCore(self.coreGene)
            output: List = self.bconnect.bioBlast(self.coreGene)
            coregene = genome.getCore()
            pathways: List = self.repProcedure(output, self.bp, coregene, self.similarity)  # noqa
            pathway.append(pathways)
        else:
            output: List = self.bconnect.bioBlast(self.seq)
            pathways: List = self.repProcedure(output, self.bp, self.seq, self.similarity)  # noqa

        return pathways

    def repProcedure(self, items: List, bp: int, coreGene: str, similarity: float):  # noqa
        pathway: List = list()
        counter: int = 0
        for item in items:
            try:
                gbfile = self.bconnect.load(item)
                rb = ReadGB(gbfile)
                genome = rb.readfile()
                organism = rb.provideOrg()

                if not genome:
                    sys.exit()
                
                genes: List = genome.findCoreGeneBySimilarity(coreGene, similarity)
                if genes:
                    for gene in genes:
                        path: List = genome.buildsimilarity(gene, bp)
                        if path:
                            pathway.append(path)
                counter = counter + 1
                if(counter % 3 == 0):  # Used because of NCBI's policy on requests without API key. with API key, change to 10  # noqa
                    time.sleep(3)
            except:
                print("failed to download {0}".format(item))
        return pathway

    def flatten(self, path: List):
        # TODO: This is just hack, I'll need to work on
        # the data structure output to remove this extra step
        """ Flattens the list i.e. removes nested list in output"""
        return [*chain.from_iterable(x if isinstance(x[0], tuple) else [x] for x in path)]  # noqa

# finder = Finder('adjon081@uottawa.ca', accession="KK037233.1", coreGene="EWM62968.1", bp=500, similarity=0.8)  # noqa
# paths = finder.finder()
# for path in paths:
#     if path:
#         print('-----------------------------------------\n')
#         print(path)


# finder.write_fasta(paths)

# finder = Finder(accession='CP013839.1', coreGene="MGAS23530_0009")
# paths = finder.finder()
# paths.produce()
# using protein accession number finder = Finder(coreGene="AAD07482.1")
# using a sequence finder = Finder(seq="AUGTTTYRRSTVVVVALLISSTUCCYTADQ")

# find small genome and build unittest based on that
# Mycoplasma genitalium
# Vibrio cholerae
# Escherichia coli
# Bacilus subtilis
