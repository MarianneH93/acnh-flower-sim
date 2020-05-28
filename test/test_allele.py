import unittest

from genetics.gene import SimpleGene
from genetics.allele import SimpleAllele


class TestAllele(unittest.TestCase):

    def test_eq(self):
        self.assertEqual(SimpleAllele.RECESSIVE, SimpleAllele.RECESSIVE)
        self.assertEqual(SimpleAllele.DOMINANT, SimpleAllele.DOMINANT)
        self.assertNotEqual(SimpleAllele.RECESSIVE, SimpleAllele.DOMINANT)

    def test_hash(self):
        self.assertIsInstance(hash(SimpleAllele.RECESSIVE), int)
        self.assertIsInstance(hash(SimpleAllele.DOMINANT), int)

    def test_add(self):
        with self.assertRaises(TypeError):
            SimpleAllele.RECESSIVE + 1

        self.assertIsInstance(SimpleAllele.RECESSIVE + SimpleAllele.DOMINANT, SimpleGene)
        self.assertEqual(SimpleAllele.RECESSIVE + SimpleAllele.RECESSIVE, SimpleGene.RECESSIVE)
        self.assertEqual(SimpleAllele.RECESSIVE + SimpleAllele.DOMINANT, SimpleGene.HETEROZYGOUS)
        self.assertEqual(SimpleAllele.DOMINANT + SimpleAllele.DOMINANT, SimpleGene.DOMINANT)

if __name__ == '__main__':
    unittest.main()