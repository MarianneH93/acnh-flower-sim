import argparse

import acnh.species
from genetics.genesequence import SimpleGeneSequence


species_map = {'rose': acnh.species.rose,
               'tulip': acnh.species.tulip,
               'pansy': acnh.species.pansy,
               'cosmo': acnh.species.cosmo,
               #'lily': acnh.species.lily,
               #'hyacinth': acnh.species.hyacinth,
               #'windflower': acnh.species.windflower,
               'mum': acnh.species.mum}


def main():
    parser = argparse.ArgumentParser(description='Calculate the possible results of a crossbreed.')
    parser.add_argument('species', metavar='S', type=str)
    parser.add_argument('parent1', metavar='P1', type=int)
    parser.add_argument('parent2', metavar='P2', type=int)
    args = parser.parse_args()

    try:
        species = species_map[args.species]
    except KeyError:
        print('Invalid species: {species}'.format(args.species))
        return

    parent1 = SimpleGeneSequence.from_ordinal(args.parent1, species.genome_len)
    parent2 = SimpleGeneSequence.from_ordinal(args.parent2, species.genome_len)

    cross = parent1.cross(parent2)
    for child in cross:
        phenotype = species.phenotype_of(child)
        print('{pheno} {ord}: {prob}'.format(pheno=phenotype,
                                             ord=child.ordinal,
                                             prob=cross[child]))


if __name__ == '__main__':
    main()