from __future__ import annotations
from typing import Set
from collections import Counter, Sequence

import genetics.genesequence


class GeneCrossResults:
    """
    Acts similar to a Counter, keeping track of the number of times that a
    particular sequence appears in cross results.
    """
    def __init__(self, parents: Sequence[genetics.genesequence.SimpleGeneSequence]):
        if len(parents) < 2:
            parents = parents + parents

        self.parents: Sequence[genetics.genesequence.SimpleGeneSequence] = list(parents)
        self.contents: Counter[genetics.genesequence.SimpleGeneSequence] = Counter()

    def __len__(self):
        return len(self.contents)

    def __contains__(self, item):
        return item in self.contents

    def __iter__(self):
        for item in self.contents:
            yield item

    def __getitem__(self, key: genetics.genesequence.SimpleGeneSequence):
        return self.contents[key]

    def __add__(self, other) -> GeneCrossResults:
        """
        Add multiple GeneCrossResults together, giving the odds for all
        possible permutations.
        :param other:
        :return:
        """

        if not isinstance(other, GeneCrossResults):
            raise TypeError

        parents = []

        for i in range(len(self.parents)):
            parents.append(self.parents[i] + other.parents[i])

        results = GeneCrossResults(self.parents)

        for sequence1 in self.contents:
            for sequence2 in other.contents:
                results.add_genesequence(sequence1+sequence2,
                                         num=self.contents[sequence1] * other.contents[sequence2])

        return results

    def add_genesequence(self, sequence: genetics.genesequence.SimpleGeneSequence, num: int = 1):
        self.contents[sequence] += num

    @property
    def total(self):
        return sum(self.contents.values())

    def probability(self, sequence: genetics.genesequence.SimpleGeneSequence):
        return self[sequence] / self.total

    def __str__(self):
        result = ""

        for sequence, count in self.contents.most_common():
            result += "{}: {}%\n".format(sequence, self.probability(sequence)*100)

        return result