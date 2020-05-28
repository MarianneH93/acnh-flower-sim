from __future__ import annotations
from collections import Counter
from typing import Sequence, Union

import genetics.gene


class SimpleGeneSequence:
    def __init__(self, gene_sequence: Union[Sequence[genetics.gene.SimpleGene], genetics.gene.SimpleGene]):
        self.sequence: Sequence[genetics.gene.SimpleGene]
        if isinstance(gene_sequence, Sequence):
            self.sequence = tuple(gene_sequence)
        elif isinstance(gene_sequence, genetics.gene.SimpleGene):
            self.sequence = (gene_sequence, )

    @classmethod
    def from_ordinal(cls, ordinal: int, seq_len: int):
        """
        Alternative method of creating gene sequences. Uses an int from
        0 (completely recessive) to MAX (all dominant).
        :param ordinal: The gene sequence to generate.
        :param seq_len: Number of simple genes.
        :return: A new gene sequence.
        """
        divisor = ordinal
        genes = []
        while len(genes) < seq_len:
            quotient = divisor // 3
            remainder = divisor % 3

            gene = genetics.gene.SimpleGene(remainder)
            genes.insert(0, gene)
            divisor = quotient

        if divisor > 0:
            raise ValueError

        return cls(genes)

    def __add__(self, other):
        if not isinstance(other, SimpleGeneSequence):
            raise TypeError

        return SimpleGeneSequence(self.sequence + other.sequence)

    def __eq__(self, other: SimpleGeneSequence):
        if not isinstance(other, SimpleGeneSequence):
            raise TypeError

        return self.sequence == other.sequence

    def __iter__(self):
        for gene in self.sequence:
            yield gene

    def __repr__(self):
        genes = ", ".join([str(gene) for gene in self.sequence])
        return "SimpleGeneSequence({})".format(genes)

    def __hash__(self):
        return hash(self.sequence)

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        return self.sequence[key]

    def __add__(self, other):
        return SimpleGeneSequence(self.sequence + other.sequence)

    def cross(self, other: SimpleGeneSequence):
        if not len(self) == len(other):
            raise NotImplementedError

        working_copy = None

        for i in range(len(self)):
            cross = self.sequence[i].cross(other.sequence[i])
            if working_copy is None:
                working_copy = cross
            else:
                working_copy = working_copy + cross

        return working_copy

    @property
    def ordinal(self):
        result = 0
        for gene in self.sequence:
            result *= 3
            result += gene.value

        return result