from enum import Enum

import acnh.species
from genetics.specimen import Specimen
from genetics.gene import SimpleGene
from genetics.genesequence import SimpleGeneSequence
from genetics.familytree import FamilyTree


class Species(Enum):
    ROSE = 0
    TULIP = 1
    PANSY = 2
    COSMO = 3
    LILY = 4
    HYACINTH = 5
    WINDFLOWER = 6
    MUM = 7


program_setting = Species.ROSE

rose_seed_genomes = [(SimpleGene.RECESSIVE, SimpleGene.RECESSIVE, SimpleGene.HETEROZYGOUS, SimpleGene.RECESSIVE),
                (SimpleGene.RECESSIVE, SimpleGene.HETEROZYGOUS, SimpleGene.RECESSIVE, SimpleGene.RECESSIVE),
                (SimpleGene.DOMINANT, SimpleGene.RECESSIVE, SimpleGene.RECESSIVE, SimpleGene.HETEROZYGOUS)]
rose_seed_genomes = [SimpleGeneSequence(seq) for seq in rose_seed_genomes]

tulip_seed_genomes = [SimpleGeneSequence.from_ordinal(1, 3),
                      SimpleGeneSequence.from_ordinal(6, 3),
                      SimpleGeneSequence.from_ordinal(19, 3)]

cosmo_seed_genomes = [SimpleGeneSequence.from_ordinal(1, 3),
                      SimpleGeneSequence.from_ordinal(7, 3),
                      SimpleGeneSequence.from_ordinal(18, 3)]

mum_seed_genomes = [
    SimpleGeneSequence.from_ordinal(1, 3),
    SimpleGeneSequence.from_ordinal(6, 3),
    SimpleGeneSequence.from_ordinal(18, 3)
]

if program_setting == Species.ROSE:
    seeds = rose_seed_genomes
    species = acnh.species.rose
elif program_setting == Species.TULIP:
    seeds = tulip_seed_genomes
    species = acnh.species.tulip
elif program_setting == Species.COSMO:
    seeds = cosmo_seed_genomes
    species = acnh.species.cosmo
elif program_setting == Species.MUM:
    seeds = mum_seed_genomes
    species = acnh.species.mum

familytree = FamilyTree(seeds, species)

for i in range(3 ** species.genome_len):
    parent1 = SimpleGeneSequence.from_ordinal(i, species.genome_len)
    for j in range(i, 3 ** species.genome_len):
        parent2 = SimpleGeneSequence.from_ordinal(j, species.genome_len)

        cross = parent1.cross(parent2)
        familytree.add(cross, prune=False, verbose=False)

familytree.shortest_path()
familytree.highest_probability()

for i in range(3 ** species.genome_len):
    seq = SimpleGeneSequence.from_ordinal(i, species.genome_len)
    print("Target: {} {}".format(Specimen(species, seq).phenotype, i))
    print("HIGHEST PROBABILITY")
    familytree.print_highest_prob(seq)
    print("SHORTEST PATH")
    familytree.print_shortest_path(seq)
    print()
    print()