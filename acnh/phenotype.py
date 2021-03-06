from enum import Enum, unique
from typing import List
import genetics.genesequence
from genetics.phenotype import Phenotype


@unique
class Colour(Enum):
    WHITE = 0
    PURPLE = 1
    YELLOW = 2
    RED = 3
    PINK = 4
    ORANGE = 5
    BLACK = 6
    BLUE = 7
    GREEN = 8


""" Maps for each species from an ordinal (list index) to a phenotype. """


class SimpleMapPhenotype:

    def __init__(self, pheno_map: List[Phenotype]):
        self.map = pheno_map

    def from_genome(self, sequence: genetics.genesequence.SimpleGeneSequence):
        return self.map[sequence.ordinal]


rose_map = [
    Colour.WHITE, Colour.WHITE, Colour.WHITE,
    Colour.WHITE, Colour.WHITE, Colour.WHITE,
    Colour.PURPLE, Colour.PURPLE, Colour.PURPLE,
    Colour.YELLOW, Colour.YELLOW, Colour.YELLOW,
    Colour.WHITE, Colour.WHITE, Colour.WHITE,
    Colour.PURPLE, Colour.PURPLE, Colour.PURPLE,
    Colour.YELLOW, Colour.YELLOW, Colour.YELLOW,
    Colour.YELLOW, Colour.YELLOW, Colour.YELLOW,
    Colour.WHITE, Colour.WHITE, Colour.WHITE,
    Colour.RED, Colour.PINK, Colour.WHITE,
    Colour.RED, Colour.PINK, Colour.WHITE,
    Colour.RED, Colour.PINK, Colour.PURPLE,
    Colour.ORANGE, Colour.YELLOW, Colour.YELLOW,
    Colour.RED, Colour.PINK, Colour.WHITE,
    Colour.RED, Colour.PINK, Colour.PURPLE,
    Colour.ORANGE, Colour.YELLOW, Colour.YELLOW,
    Colour.ORANGE, Colour.YELLOW, Colour.YELLOW,
    Colour.RED, Colour.PINK, Colour.WHITE,
    Colour.BLACK, Colour.RED, Colour.PINK,
    Colour.BLACK, Colour.RED, Colour.PINK,
    Colour.BLACK, Colour.RED, Colour.PINK,
    Colour.ORANGE, Colour.ORANGE, Colour.YELLOW,
    Colour.RED, Colour.RED, Colour.WHITE,
    Colour.BLACK, Colour.RED, Colour.PURPLE,
    Colour.ORANGE, Colour.ORANGE, Colour.YELLOW,
    Colour.ORANGE, Colour.ORANGE, Colour.YELLOW,
    Colour.BLUE, Colour.RED, Colour.WHITE
]

rose = SimpleMapPhenotype(rose_map)

tulip_map = [
    Colour.WHITE, Colour.WHITE, Colour.WHITE,
    Colour.YELLOW, Colour.YELLOW, Colour.WHITE,
    Colour.YELLOW, Colour.YELLOW, Colour.YELLOW,
    Colour.RED, Colour.PINK, Colour.WHITE,
    Colour.ORANGE, Colour.YELLOW, Colour.YELLOW,
    Colour.ORANGE, Colour.YELLOW, Colour.YELLOW,
    Colour.BLACK, Colour.RED, Colour.RED,
    Colour.BLACK, Colour.RED, Colour.RED,
    Colour.PURPLE, Colour.PURPLE, Colour.PURPLE
]

tulip = SimpleMapPhenotype(tulip_map)

pansy_map = [
    Colour.WHITE, Colour.WHITE, Colour.BLUE,
    Colour.YELLOW, Colour.YELLOW, Colour.BLUE,
    Colour.YELLOW, Colour.YELLOW, Colour.YELLOW,
    Colour.RED, Colour.RED, Colour.BLUE,
    Colour.ORANGE, Colour.ORANGE, Colour.ORANGE,
    Colour.YELLOW, Colour.YELLOW, Colour.YELLOW,
    Colour.RED, Colour.RED, Colour.PURPLE,
    Colour.RED, Colour.RED, Colour.PURPLE,
    Colour.ORANGE, Colour.ORANGE, Colour.PURPLE
]

pansy = SimpleMapPhenotype(pansy_map)

cosmo_map = [
    Colour.WHITE, Colour.WHITE, Colour.WHITE,
    Colour.YELLOW, Colour.YELLOW, Colour.WHITE,
    Colour.YELLOW, Colour.YELLOW, Colour.YELLOW,
    Colour.PINK, Colour.PINK, Colour.PINK,
    Colour.ORANGE, Colour.ORANGE, Colour.PINK,
    Colour.ORANGE, Colour.ORANGE, Colour.ORANGE,
    Colour.RED, Colour.RED, Colour.RED,
    Colour.ORANGE, Colour.ORANGE, Colour.RED,
    Colour.BLACK, Colour.BLACK, Colour.RED
]

cosmo = SimpleMapPhenotype(cosmo_map)

lily_map = [
    Colour.WHITE, Colour.WHITE, Colour.WHITE,
    Colour.YELLOW, Colour.WHITE, Colour.WHITE,
    Colour.YELLOW, Colour.YELLOW, Colour.WHITE,
    Colour.RED, Colour.PINK, Colour.WHITE,
    Colour.ORANGE, Colour.YELLOW, Colour.YELLOW,
    Colour.ORANGE, Colour.YELLOW, Colour.YELLOW,
    Colour.BLACK, Colour.RED, Colour.PINK,
    Colour.BLACK, Colour.RED, Colour.PINK,
    Colour.ORANGE, Colour.ORANGE, Colour.WHITE
]

lily = SimpleMapPhenotype(lily_map)

hyacinth_map = [
    Colour.WHITE, Colour.WHITE, Colour.BLUE,
    Colour.YELLOW, Colour.YELLOW, Colour.WHITE,
    Colour.YELLOW, Colour.YELLOW, Colour.YELLOW,
    Colour.RED, Colour.PINK, Colour.WHITE,
    Colour.ORANGE, Colour.YELLOW, Colour.YELLOW,
    Colour.ORANGE, Colour.YELLOW, Colour.YELLOW,
    Colour.RED, Colour.RED, Colour.RED,
    Colour.BLUE, Colour.BLUE, Colour.RED,
    Colour.PURPLE, Colour.PURPLE, Colour.PURPLE
]

hyacinth = SimpleMapPhenotype(hyacinth_map)

windflower_map = [
    Colour.WHITE, Colour.WHITE, Colour.BLUE,
    Colour.ORANGE, Colour.ORANGE, Colour.BLUE,
    Colour.ORANGE, Colour.ORANGE, Colour.ORANGE,
    Colour.RED, Colour.RED, Colour.BLUE,
    Colour.PINK, Colour.PINK, Colour.PINK,
    Colour.ORANGE, Colour.ORANGE, Colour.ORANGE,
    Colour.RED, Colour.RED, Colour.PURPLE,
    Colour.RED, Colour.RED, Colour.PURPLE,
    Colour.PINK, Colour.PINK, Colour.PURPLE
]

windflower = SimpleMapPhenotype(windflower_map)

mum_map = [
    Colour.WHITE, Colour.WHITE, Colour.PURPLE,
    Colour.YELLOW, Colour.YELLOW, Colour.WHITE,
    Colour.YELLOW, Colour.YELLOW, Colour.YELLOW,
    Colour.PINK, Colour.PINK, Colour.PINK,
    Colour.YELLOW, Colour.RED, Colour.PINK,
    Colour.PURPLE, Colour.PURPLE, Colour.PURPLE,
    Colour.RED, Colour.RED, Colour.RED,
    Colour.PURPLE, Colour.PURPLE, Colour.RED,
    Colour.GREEN, Colour.GREEN, Colour.RED
]

mum = SimpleMapPhenotype(mum_map)