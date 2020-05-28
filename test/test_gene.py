import unittest

from genetics.allele import SimpleAllele
from genetics.gene import SimpleGene
from genetics.genesequence import SimpleGeneSequence


class TestGene(unittest.TestCase):
    def test_hash(self):
        self.assertIsInstance(hash(SimpleGene.RECESSIVE), int)

    def test_add(self):
        with self.assertRaises(TypeError):
            SimpleGene.RECESSIVE + SimpleAllele.RECESSIVE

        self.assertIsInstance(SimpleGene.RECESSIVE + SimpleGene.DOMINANT, SimpleGeneSequence)

    def test_alleles(self):
        self.assertIn(SimpleAllele.RECESSIVE, SimpleGene.RECESSIVE.alleles)
        self.assertNotIn(SimpleAllele.DOMINANT, SimpleGene.RECESSIVE.alleles)
        self.assertIn(SimpleAllele.RECESSIVE, SimpleGene.HETEROZYGOUS.alleles)
        self.assertIn(SimpleAllele.DOMINANT, SimpleGene.HETEROZYGOUS.alleles)
        self.assertIn(SimpleAllele.DOMINANT, SimpleGene.DOMINANT.alleles)
        self.assertNotIn(SimpleAllele.RECESSIVE, SimpleGene.DOMINANT.alleles)

    def test_cross(self):
        result = SimpleGene.HETEROZYGOUS.cross(SimpleGene.HETEROZYGOUS)

        self.assertEqual(result[SimpleGeneSequence((SimpleGene.HETEROZYGOUS,))], 2)
        self.assertEqual(result[SimpleGeneSequence((SimpleGene.DOMINANT,))], 1)
        self.assertEqual(result[SimpleGeneSequence((SimpleGene.RECESSIVE,))], 1)

        print(result)

if __name__ == '__main__':
    unittest.main()