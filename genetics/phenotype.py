import genetics.genesequence

from abc import ABC


class Phenotype(ABC):
    def __init__(self):
        pass

    def from_genome(self, sequence: genetics.genesequence.SimpleGeneSequence):
        raise NotImplementedError
