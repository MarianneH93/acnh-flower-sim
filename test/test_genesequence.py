import unittest

from genetics.gene import SimpleGene
from genetics.genesequence import SimpleGeneSequence


class TestGeneSequence(unittest.TestCase):
    def test_cross(self):
        seqa = SimpleGeneSequence((SimpleGene.HETEROZYGOUS, SimpleGene.HETEROZYGOUS))
        cross = seqa.cross(seqa)
        raise NotImplementedError

    def test_from_ordinal(self):
        seq = SimpleGeneSequence.from_ordinal(0, 4)
        self.assertEqual(len(seq), 4)
        for gene in seq:
            self.assertEqual(gene, SimpleGene.RECESSIVE)

        seq = SimpleGeneSequence.from_ordinal(1, 2)
        self.assertEqual(len(seq), 2)
        self.assertEqual(seq[0], SimpleGene.RECESSIVE)
        self.assertEqual(seq[1], SimpleGene.HETEROZYGOUS)

if __name__ == '__main__':
    unittest.main()