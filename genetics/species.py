import genetics.genesequence
import genetics.phenotype


class Species:

    __slots__ = ['phenotypes', 'genome_len']

    def __init__(self, phenotypes: genetics.phenotype.Phenotype, genome_len: int):
        self.phenotypes = phenotypes
        self.genome_len = genome_len

    def phenotype_of(self, gene_seq: genetics.genesequence.SimpleGeneSequence):
        return self.phenotypes.from_genome(gene_seq)

    def __eq__(self, other):
        if not isinstance(other, Species):
            return NotImplementedError

        return self.phenotypes == other.phenotypes and self.genome_len == other.genome_len