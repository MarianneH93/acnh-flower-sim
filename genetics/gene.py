from __future__ import annotations
from enum import Enum
from typing import Collection

import genetics.allele
import genetics.genesequence
import genetics.cross


class SimpleGene(Enum):
    """
    Represents the two alleles that make up a gene.

    It is simpler to model this as an enum since there are only three
    possibilities.
    """
    RECESSIVE = 0
    HETEROZYGOUS = 1
    DOMINANT = 2

    def __add__(self, other):
        if not isinstance(other, SimpleGene):
            raise TypeError

        genes = (self, other)

        return genetics.genesequence.SimpleGeneSequence(genes)

    @property
    def alleles(self) -> Collection[genetics.allele.SimpleAllele]:
        result = set()
        if not self == SimpleGene.RECESSIVE:
            result.add(genetics.allele.SimpleAllele.DOMINANT)
        if not self == SimpleGene.DOMINANT:
            result.add(genetics.allele.SimpleAllele.RECESSIVE)

        return result

    def cross(self, other: SimpleGene):

        results = genetics.cross.GeneCrossResults(
            [genetics.genesequence.SimpleGeneSequence(self),
             genetics.genesequence.SimpleGeneSequence(other)])

        for allele1 in self.alleles:
            for allele2 in other.alleles:
                results.add_genesequence(genetics.genesequence.SimpleGeneSequence(allele1 + allele2))

        return results