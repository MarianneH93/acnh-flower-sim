from __future__ import annotations
from typing import Dict, Collection, FrozenSet, Set

import genetics.genesequence
import genetics.cross


class ParentNode:
    __slots__ = ['parents', 'children', 'specimen', 'probability']

    def __init__(self, parents: Collection[genetics.genesequence.SimpleGeneSequence]):
        if len(parents) < 1 or len(parents) > 2:
            raise NotImplementedError

        for parent in parents:
            if not isinstance(parent, genetics.genesequence.SimpleGeneSequence):
                raise TypeError

        self.parents: FrozenSet[ChildNode]
        temp_set = set()
        for parent in parents:
            temp_set.add(ChildNode(parent))

        self.parents = frozenset(temp_set)

        self.children: Dict[ChildNode, float] = dict()

    def __hash__(self):
        return hash(self.parents)

    def __eq__(self, other):
        if not isinstance(other, ParentNode):
            raise NotImplementedError

        return self.parents == other.parents

    @property
    def generation(self):
        """
        Same as the generation of the child nodes.
        :return:
        """
        if len(self.parents):
            return max(parent.generation for parent in self.parents) + 1

        return 1


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
            return max([node.generation for node in self.shortest_path_parents])

        return 1


class FamilyTree:
    def __init__(self, first_gen: Collection[genetics.genesequence.SimpleGeneSequence]):
        """
        Contains the first generation and a set of all genomes that exist
        somewhere in the family tree. The current implementation is used to
        solve for the shortest path from the first gen to all genomes so should
        not include duplicates (could cause loops).
        :param first_gen:
        """
        self.parent_nodes: Set[ParentNode] = set()
        self.child_nodes: Set[ChildNode] = set()
        for seq in first_gen:
            child_node = ChildNode(seq)
            self.child_nodes.add(child_node)

    def __getitem__(self, item: genetics.genesequence.SimpleGeneSequence):
        if not isinstance(item, genetics.genesequence.SimpleGeneSequence):
            raise NotImplementedError

        return self.child_nodes[ChildNode(item)]

    def add(self, results: genetics.cross.GeneCrossResults):

        new_parent_node = ParentNode(results.parents)
        if new_parent_node not in self.parent_nodes:
            self.parent_nodes.add(new_parent_node)
        else:
            for node in self.parent_nodes:
                if node == new_parent_node:
                    new_parent_node = node

        for item in results:
            new_node = ChildNode(item)
            if new_node not in self.child_nodes:
                self.child_nodes.add(new_node)
            else:
                # Retrieve existing node from the set
                for node in self.child_nodes:
                    if node == new_node:
                        new_node = node

            probability = results.probability(new_node.specimen)
            new_parent_node.children[new_node] = probability
            if new_node.highest_prob_parents:
                if list(new_node.highest_prob_parents)[0].children[new_node] < probability:
                    new_node.highest_prob_parents.clear()
                    new_node.highest_prob_parents.add(new_parent_node)
                elif list(new_node.highest_prob_parents)[0].children[new_node] == probability:
                    new_node.highest_prob_parents.add(new_parent_node)
            else:
                new_node.highest_prob_parents.add(new_parent_node)

            generation = new_parent_node.generation
            if new_node.shortest_path_parents:
                if list(new_node.shortest_path_parents)[0].generation > generation:
                    new_node.shortest_path_parents.clear()
                    new_node.shortest_path_parents.add(new_parent_node)
                elif list(new_node.shortest_path_parents)[0].generation == probability:
                    new_node.shortest_path_parents.add(new_parent_node)
            else:
                new_node.shortest_path_parents.add(new_parent_node)