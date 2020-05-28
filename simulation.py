from collections import Counter

import acnh.species
from genetics.specimen import Specimen
from genetics.gene import SimpleGene
from genetics.genesequence import SimpleGeneSequence
from genetics.familytree import FamilyTree

seed_genomes = [(SimpleGene.RECESSIVE, SimpleGene.RECESSIVE, SimpleGene.HETEROZYGOUS, SimpleGene.RECESSIVE),
                (SimpleGene.RECESSIVE, SimpleGene.HETEROZYGOUS, SimpleGene.RECESSIVE, SimpleGene.RECESSIVE),
                (SimpleGene.DOMINANT, SimpleGene.RECESSIVE, SimpleGene.RECESSIVE, SimpleGene.HETEROZYGOUS)]
seed_genomes = [SimpleGeneSequence(seq) for seq in seed_genomes]

familytree = FamilyTree(seed_genomes)

last_gen = set(seed_genomes)
next_gen = set()
gen = 1
while True:
    print("Generation {}".format(gen))
    for speca in last_gen:
        for specb in last_gen:
            cross = speca.cross(specb)

            new_phenos = Counter()
            for genseq in cross:
                if genseq not in last_gen:
                    next_gen.add(genseq)
                specimen = Specimen(acnh.species.rose, genseq)
                new_phenos[specimen.phenotype] += 1

            familytree.add(cross)

    if last_gen.union(next_gen) == last_gen:
        print("No new genotypes, stopping")
        break
    last_gen = last_gen.union(next_gen)
    print("{} new genotypes, {} total".format(len(next_gen), len(last_gen)))
    next_gen.clear()
    gen += 1