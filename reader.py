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
        self.GENOME: Dict = Genome()
        try:
            self._org = Organism(self._record.annotations['organism'], self._record.annotations['accessions'])  # noqa
        except KeyError:
            # Throw exception if file is empty
            pass

    def provideOrg(self) -> Tuple:
        """Provide indentifier about organism created."""
        return self._org.info()

    def readfile(self) -> Dict:
        """Read the gb file and return a new Genome."""
        try:
            prot_id: str = 'N/A'
            locus: str = 'N/A'
            product: str = 'N/A'
            gene: str = 'N/A'
            translation: str = 'N/A'
            for feature in self._record.features:
                if "CDS" in feature.type:
                    for value in feature.qualifiers.keys():
                        if 'protein_id' in feature.qualifiers.keys():
                            prot_id = feature.qualifiers['protein-id'][0]
                        if 'product' in feature.qualifiers.keys():
                            product = feature.qualifiers['product'][0]
                        if 'gene' in feature.qualifiers.keys():
                            gene = feature.qualifiers['gene'][0]
                        if 'locus_tag' in feature.qualifiers.keys():
                            locus = feature.qualifiers['locus_tag'][0]
                        if 'translation' in feature.qualifiers.keys():
                            translation = feature.qualifiers['translation'][0]
                        strand: int = int(feature.strand, base=10)
                        location: Tuple = (feature.location.start.position, feature.location.end.position)  # noqa
                        gene: GENE = Gene(gene, locus, product, prot_id, translation, location, strand)  # noqa
                        self.GENOME.addGene(gene)
        except KeyError:
            # if no CDS is found in the file or all features are empty.
            pass
            print(1)
        return self.GENOME
