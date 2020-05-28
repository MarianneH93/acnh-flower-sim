from __future__ import annotations
from enum import Enum

import genetics.gene


class SimpleAllele(Enum):
    """
    Represents a single allele, either recessive or dominant. Can be combined
    into genes.
    TODO: Should there be different truth values for the different alleles?
    """
    RECESSIVE = 0
    DOMINANT = 1

    def __add__(self, other: SimpleAllele):
        if not isinstance(other, SimpleAllele):
            raise TypeError

        if not self == other:
            return genetics.gene.SimpleGene.HETEROZYGOUS
        if self == SimpleAllele.RECESSIVE:
            return genetics.gene.SimpleGene.RECESSIVE
        return genetics.gene.SimpleGene.DOMINANT