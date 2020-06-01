from __future__ import annotations
from collections import Counter
from typing import Dict, Collection, FrozenSet, Set

from genetics.genesequence import SimpleGeneSequence
import genetics.cross
import genetics.species
from genetics.specimen import Specimen


class ParentNode:
    __slots__ = ['parents', 'children', 'specimen', 'probability']

    def __init__(self, parents: Collection[genetics.genesequence.SimpleGeneSequence],
                 existing: Set[ChildNode]=[]):
        if len(parents) > 2:
            raise NotImplementedError

        for parent in parents:
            if not isinstance(parent, genetics.genesequence.SimpleGeneSequence):
                raise TypeError

        self.parents: FrozenSet[ChildNode]
        temp_set = set()
        for parent in parents:
            new_child = ChildNode(parent)
            if existing:
                if new_child in existing:
                    for child_node in existing:
                        if child_node == new_child:
                            temp_set.add(child_node)
                            break
                else:
                    temp_set.add(new_child)
                    existing.add(new_child)
            else:
                temp_set.add(new_child)
        self.parents = frozenset(temp_set)
        self.children: Dict[ChildNode, float] = dict()

    def __hash__(self):
        return hash(self.parents)

    def __eq__(self, other):
        if not isinstance(other, ParentNode):
            raise NotImplementedError

        return self.parents == other.parents

    @classmethod
    def zeroth_gen(cls, children):
        result = cls([])
        for child in children:
            child_node = ChildNode(child)
            result.children[child_node] = 1.0
            child_node.shortest_path_parents.add(result)
            child_node.highest_prob_parents.add(result)

        return result

    @property
    def generation(self):
        """
        Same as the generation of the child nodes.
        :return:
        """
        if len(self.parents):
            return max(parent.generation for parent in self.parents) + 1
        return 0

    def ancestors_in_shortest_path(self):
        """
        Includes the parents in the current node as well.
        :return:
        """
        ancestors = set()
        queue = list(self.parents)
        while queue:
            child_node = queue.pop()
            ancestors.update((child_node.specimen,))
            ancestors.update(child_node.ancestors_in_shortest_path)

        return ancestors

class ChildNode:
    def __init__(self, specimen: genetics.genesequence.SimpleGeneSequence):
        if not isinstance(specimen, genetics.genesequence.SimpleGeneSequence):
            raise TypeError
        self.specimen = specimen
        self.shortest_path_parents: Set[ParentNode] = set()
        self.highest_prob_parents: Set[ParentNode] = set()

    def __hash__(self):
        return hash(self.specimen)

    def __eq__(self, other: ChildNode):
        if not isinstance(other, ChildNode):
            raise NotImplementedError

        return self.specimen == other.specimen

    @property
    def generation(self):
        if len(self.shortest_path_parents):
            return list(self.shortest_path_parents)[0].generation
        else:
            return 10000
    
    @property
    def ancestors_in_shortest_path(self) -> Collection[ParentNode]:
        ancestors: Set[ChildNode] = set()
        queue: Collection[ParentNode] = list(self.shortest_path_parents)

        while queue:
            parent_node = queue.pop()
            for child_node in parent_node.parents:
                if child_node.specimen not in ancestors:
                    ancestors.add(child_node.specimen)
                    queue += child_node.shortest_path_parents

        return ancestors


class FamilyTree:
    def __init__(self, first_gen: Collection[genetics.genesequence.SimpleGeneSequence],
                 species: genetics.species.Species):
        """
        Contains the first generation and a set of all genomes that exist
        somewhere in the family tree. The current implementation is used to
        solve for the shortest path from the first gen to all genomes so should
        not include duplicates (could cause loops).
        :param first_gen:
        """
        self.parent_nodes: Set[ParentNode] = set()
        self.child_nodes: Set[ChildNode] = set()
        first_node = ParentNode.zeroth_gen(first_gen)
        self.parent_nodes.add(first_node)
        for child_node in first_node.children:
            self.child_nodes.add(child_node)

        self.species = species

    def __getitem__(self, item: genetics.genesequence.SimpleGeneSequence):
        if not isinstance(item, genetics.genesequence.SimpleGeneSequence):
            raise NotImplementedError

        return self.child_nodes[ChildNode(item)]

    def add(self, results: genetics.cross.GeneCrossResults,
            prune: bool = False,
            verbose: bool = True):

        for result in results.parents:
            if not len(result) == self.species.genome_len:
                raise NotImplementedError
        if verbose:
            parent_specimens = [Specimen(self.species, parent) for parent in results.parents]
            if len(parent_specimens) == 1:
                parent_specimens.append(parent_specimens[0])
            print("Adding new ParentNode: {} + {}".format(str(parent_specimens[0]), str(parent_specimens[1])))

        new_parent_node = ParentNode(results.parents, existing=self.child_nodes)
        if new_parent_node not in self.parent_nodes:
            self.parent_nodes.add(new_parent_node)
        else:
            for node in self.parent_nodes:
                if node == new_parent_node:
                    new_parent_node = node
                    break
            parent_specimens = [Specimen(self.species, parent) for parent in results.parents]
            if len(parent_specimens) == 1:
                parent_specimens.append(parent_specimens[0])
            if verbose:
                print("Found ParentNode: {} + {}".format(str(parent_specimens[0]), str(parent_specimens[1])))

        phenotype_counter = Counter()
        for child in results:
            phenotype_counter[genetics.specimen.Specimen(self.species, child).phenotype] += 1

        if verbose:
            print("Phenotypes: {}".format(phenotype_counter))

        for item in results:
            if not len(item) == self.species.genome_len:
                raise NotImplementedError

            specimen = genetics.specimen.Specimen(self.species, item)

            new_node = ChildNode(item)
            if new_node not in self.child_nodes:
                if verbose:
                    print("Adding new ChildNode: {}".format(specimen))
                self.child_nodes.add(new_node)
            else:
                # Retrieve existing node from the set
                for node in self.child_nodes:
                    if node == new_node:
                        if verbose:
                            print("Found existing child node: {}".format(Specimen(self.species, node.specimen)))
                        new_node = node
                        break

            # Adjust highest probability
            probability = results.probability(new_node.specimen)
            new_parent_node.children[new_node] = probability

    @staticmethod
    def whole_collection_has_shortest_path(child_nodes: ChildNode):
        for node in child_nodes:
            if not node.shortest_path_parents:
                return False

        return True

    def count_phenotypes(self, parent_node: ParentNode):
        counter = Counter()
        for child in parent_node.children:
            specimen = Specimen(self.species, child.specimen)
            counter[specimen.phenotype] += 1
        return counter

    def update_shortest_path(self, child: ChildNode, parent: ParentNode):
        min = -1
        if child.shortest_path_parents:
            for parent_node in child.shortest_path_parents:
                if min == -1 or parent_node.generation < min:
                    min = parent_node.generation

            for parent_node in child.shortest_path_parents:
                if min < parent_node.generation:
                    child.shortest_path_parents.discard(parent_node)

        if min == -1:
            child.shortest_path_parents.add(parent)
        elif parent.generation < min:
            child.shortest_path_parents.clear()
            child.shortest_path_parents.add(parent)
        else:
            return False

        return True

    def update_highest_prob_path(self, child: ChildNode, parent: ParentNode):
        prob = 0.0
        if child.highest_prob_parents:
            for parent_node in child.highest_prob_parents:
                if parent_node.children[child] > prob:
                    prob = parent_node.children[child]
            for parent_node in child.highest_prob_parents:
                if parent_node.children[child] < prob:
                    child.highest_prob_parents.discard(parent_node)

        if child in parent.parents:
            return

        if parent.children[child] > prob:
            child.highest_prob_parents.clear()
            child.highest_prob_parents.add(parent)
        elif parent.children[child] == prob:
            child.highest_prob_parents.add(parent)

    def shortest_path(self, prune=True):

        loop_counter = 1
        while True:
            stop_flag = True
            print("Loop {}".format(loop_counter))
            loop_counter += 1
            for parent_node in self.parent_nodes:
                if self.whole_collection_has_shortest_path(parent_node.parents):
                    phenotypes = self.count_phenotypes(parent_node)

                    for child in parent_node.children:
                        specimen = Specimen(self.species, child.specimen)
                        if not prune or phenotypes[specimen.phenotype] < 2:
                            if self.update_shortest_path(child, parent_node):
                                stop_flag = False
            if stop_flag:
                return

    def highest_probability(self):
        for parent_node in self.parent_nodes:
            for child_node in parent_node.children:
                self.update_highest_prob_path(child_node, parent_node)

    def print_shortest_path(self, seq: genetics.genesequence.SimpleGeneSequence):
        childnode = ChildNode(seq)
        for node in self.child_nodes:
            if childnode == node:
                childnode = node

        parent_node: ParentNode
        for parent_node in childnode.shortest_path_parents:
            parent: ChildNode
            print("PARENTS")
            for parent in parent_node.parents:
                specimen = Specimen(self.species, parent.specimen)
                print("{} {}".format(specimen.phenotype, parent.specimen.ordinal))
            print("CHILDREN")
            for childnode in parent_node.children:
                specimen = Specimen(self.species, childnode.specimen)
                print("{} {} {}".format(parent_node.children[childnode],
                                        specimen.phenotype, childnode.specimen.ordinal))


    def print_highest_prob(self, seq: genetics.genesequence.SimpleGeneSequence):
        if not len(seq) == self.species.genome_len:
            raise NotImplementedError

        node = ChildNode(seq)
        if node in self.child_nodes:
            # Get matching node
            for child_node in self.child_nodes:
                if node == child_node:
                    node = child_node
                    break

            parent_nodes = node.highest_prob_parents
            for parent in parent_nodes:
                print("PARENTS")
                for specimen in parent.parents:
                    print(genetics.specimen.Specimen(self.species, specimen.specimen).phenotype,
                          specimen.specimen.ordinal)
                print("CHILDREN")
                for item in sorted(parent.children.items(), key=lambda item: item[1], reverse=True):
                    print(parent.children[item[0]], genetics.specimen.Specimen(self.species,
                                                                               item[0].specimen).phenotype,
                          item[0].specimen.ordinal)
        print()

    @property
    def genotypes(self) -> Set[SimpleGeneSequence]:
        result = set()
        for child_node in self.child_nodes:
            result.add(child_node.specimen)

        return result