from __future__ import annotations

import genetics.species
import genetics.genesequence


class Specimen:

    __slots__ = ['species', 'genome', 'phenotype']

    def __init__(self, species: genetics.species.Species,
                 genome: genetics.genesequence.SimpleGeneSequence):
        self.species = species
        self.genome = genome
        self.phenotype = species.phenotype_of(genome)

    def __str__(self):
        return "{species} {genome} {phenotype}".format(species=self.species,
                                                       genome=self.genome,
                                                       phenotype=self.phenotype)

    def cross(self, other: Specimen):
        if not self.species == other.species:
            raise TypeError

        return self.genome.cross(other.genome)