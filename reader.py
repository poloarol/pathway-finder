"""Read GB file."""

from utils.structures import Gene
from utils.structures import Genome
from utils.structures import Organism
from utils.structures import GENE

from Bio import SeqIO

from typing import Dict, Tuple


class ReadGB:
    """Read GB file and extract information."""

    def __init__(self, genbank):
        """Initialize the reader and store gb."""
        self._record = SeqIO.read(genbank, 'genbank')
        self.GENOME = Genome()
        try:
            # self._org = Organism(self._record.source, self._record.accessions)
            pass
        except KeyError:
            # Throw exception if file is empty
            pass

    def provideOrg(self) -> Tuple:
        """Provide indentifier about organism created."""
        return self._org.info()

    def readfile(self) -> Dict:
        """Read the gb file and return a new Genome."""
        try:
            prot_id: str = 'n/a'
            locus: str = 'n/a'
            product: str = 'n/a'
            gene: str = 'n/a'
            translation: str = 'n/a'
            for i, recs in enumerate(self._record.features):
                if recs.type == 'CDS':
                    if 'translation' in recs.qualifiers:
                        translation = recs.qualifiers['translation'][0]
                    if 'gene' in recs.qualifiers:
                        gene = recs.qualifiers['gene'][0]
                    if 'locus_tag' in recs.qualifiers:
                        locus = recs.qualifiers['locus_tag'][0]
                    if 'product' in recs.qualifiers:
                        product = recs.qualifiers['product'][0]
                    if 'protein_id' in recs.qualifiers:
                        prot_id = recs.qualifiers['protein_id'][0]
                    strand: int = int(recs.strand)
                    location: Tuple = (recs.location.start.position, recs.location.end.position)
                    g: GENE = Gene(gene, locus, product, prot_id, translation, location, strand)
                    self.GENOME.addGene(g)
        except KeyError:
            # if no CDS is found in the file or all features are empty.
            pass
        return self.GENOME
