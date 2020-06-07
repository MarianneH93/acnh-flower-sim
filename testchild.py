import argparse
from collections import Counter, Sequence

import acnh.species
from genetics.genesequence import SimpleGeneSequence
from genetics.cross import GeneCrossResults


class Test:
    def __init__(self, test_cross: SimpleGeneSequence, species, cross1, cross2):
        self.test_cross: SimpleGeneSequence = test_cross
        self.species = species
        self.cross1 = cross1
        self.cross2 = cross2


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
    parser.add_argument('subjects', metavar='subjects', nargs='+', action='append', type=int)
    args = parser.parse_args()

    try:
        species = species_map[args.species]
    except KeyError:
        print('Invalid species: {species}'.format(args.species))
        return

    subjects = [SimpleGeneSequence.from_ordinal(subject, species.genome_len) for subject in args.subjects[0]]

    for ordinal in range(3 ** species.genome_len):
        test_cross = SimpleGeneSequence.from_ordinal(ordinal, species.genome_len)

        cross_results = [subject.cross(test_cross) for subject in subjects]
        phenos: Sequence[Counter] = []
        for cross in cross_results:
            new_counter = Counter()
            for child in cross:
                new_counter[species.phenotype_of(child)] += cross.probability(child)
            phenos.append(new_counter)

        unique_phenos = []
        for pheno_index1 in range(len(phenos)):
            unique = set(phenos[pheno_index1].keys())
            for pheno_index2 in range(len(phenos)):
                if pheno_index1 == pheno_index2:
                    continue
                unique = unique.difference(set(phenos[pheno_index2].keys()))
            prob = 0.0
            for unique_pheno in unique:
                prob += phenos[pheno_index1][unique_pheno]
            unique_phenos.append((unique, prob))

        has_unique = False
        for item in unique_phenos:
            if item[0]:
                has_unique = True
                break

        if has_unique:
            print("{colour} {ordinal}".format(colour=species.phenotype_of(test_cross),
                                              ordinal=test_cross.ordinal))
            for index in range(len(unique_phenos)):
                if unique_phenos[index][0]:
                    print("Will only produce {phenos} with {colour} {ordinal}, probability={prob}".format(
                        phenos=unique_phenos[index][0],
                        colour=species.phenotype_of(subjects[index]),
                        ordinal=subjects[index].ordinal,
                        prob=unique_phenos[index][1]
                    ))

if __name__ == '__main__':
    main()