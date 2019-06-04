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
            self._org: Tuple = tuple(self._record.annotation, self._record.id)
        except:
            # Throw exception if file is empty
            pass
    
    def provideOrg(self) -> Tuple:
        """Provide indentifier about organism created."""
        return self._org
    
    def readfile(self) -> Dict:
        """Read the gb file and return a new Genome."""
        try:
            prot_id: str = 'N/A'
            locus: str = 'N/A'
            product: str = 'N/A'
            gene: str = 'N/A'
            for feature in self._record.features:
                for "CDS" in feature.type:
                    for value in feature.qualifiers.keys():
                        if 'protein_id' in feature.qualifiers.keys():
                            # feature.qualifiers['protein-id'][0]
                            pass
                        if 'product' in feature.qualifiers.keys():
                            pass
                        if 'gene' in feature.qualifiers.keys():
                            pass
                        if 'locus_tag' in feature.qualifiers.keys():
                            pass
                        strand: int = int(feature.strand, base=10)
                        location: Tuple = (feature.location.start.position, feature.location.end.position)
                        gene: GENE = Gene(gene, locus, product, prot_id, location, strand)
                        self.GENOME.addGene(gene)
        return self.GENOME