#!/usr/bin/env python
# encoding: utf-8


"""
Multiple ways of creating a graph representation of a
L{Metabolism<pyMetabolism.metabolism.Metabolism>} instance.
Generation of a stoichiometric matrix is also possible.

@author: Nikolaus Sonnenschein
@author: Moritz Beber
@contact: niko(at)foo.org
@contact: moritz(at)foo.org
@copyright: Jacobs University Bremen. All rights reserved.
@since: 2009-09-11
@attention: Graph instances and L{Metabolism<pyMetabolism.metabolism.Metabolism>}
            instances are not interlinked:
            Adding or removing reactions from one does not affect the other.
@todo: complete docs
"""


#import random
import logging
import networkx


from pyMetabolism import OptionsManager, new_property
from pyMetabolism.metabolism import Compound, Reaction, Compartment, Metabolism, CompartCompound


class BipartiteMetabolicNetwork(networkx.DiGraph):
    """
    A bipartite graph representation of a L{Metabolism<pyMetabolism.metabolism.Metabolism>}
    instance.

    @todo: many things
    """

    _counter = 0

    def __init__(self, reactions=None, name='',
             split=True, edge_info=False,* args, ** kwargs):
        """
        @param reactions: A list of L{Reaction<pyMetabolism.metabolism.Reaction>}s from
                          which the graph can be populated.
        @type reactions: C{list}

        @param split: Whether to have separate nodes for forward and backward
            direction of reversible reactions.

        @param edge_info: Whether edges should carry 

        @param name: A simple string to identify the graph if desired.
        @type name: C{str}
        """
        super(BipartiteMetabolicNetwork, self).__init__(name=name, *args, ** kwargs)
        self._options = OptionsManager()
        if name:
            self._name = str(name)
        else:
            self.__class__._counter += 1
            self._name = "BipartiteMetabolicNetwork-%d" % self.__class__._counter
        self._logger = logging.getLogger("%s.%s.%s"\
             % (self._options.main_logger_name, self.__class__.__name__, self._name))
        # maintain lists of compounds and reactions for easy iteration over them
        self._reactions = set()
        self._compounds = set()
        self._split = bool(split)
        self._edge_info = bool(edge_info)
        if reactions:
            for rxn in reactions:
                self.add_reaction(rxn)

    @new_property
    def name():
        pass

    @new_property
    def network():
        return {"fset": None, "doc": "get method"}

    @new_property
    def reactions():
        return {"fset": None, "doc": "get method"}

    @new_property
    def compounds():
        return {"fset": None, "doc": "get method"}

    @new_property
    def split():
        return {"fset": None, "doc": "get method"}

    def add_compound(self, compound):
        assert isinstance(compound, (Compound, CompartCompound))
        self._compounds.add(compound)
        self.add_node(compound)

    def remove_compound(self, compound):
        assert isinstance(compound, (Compound, CompartCompound))
        self._compounds.remove(compound)
        self.remove_node(compound)

    def add_reaction(self, rxn):
        """
        Add edges to the bipartite representation of
        L{Metabolism<pyMetabolism.metabolism.Metabolism>} as necessary
        to represent a reaction.

        @param rxn: Instance to be added.
        @type rxn: L{Reaction<pyMetabolism.metabolism.Reaction>}
        """
        assert isinstance(rxn, Reaction)
        self._reactions.add(rxn)
        # in case there are no compounds involved we add the reaction node
        # individually
        self.add_node(rxn)
        if self._split and rxn.reversible:
            subs = list(rxn.substrates)
            prods = list(rxn.products)
            coeffs = list(rxn.stoichiometry)
            back = Reaction(rxn.identifier + self._options.rev_reaction_suffix,
                prods.reverse(), subs.reverse(), coeffs.reverse())
            self._reactions.add(back)
            # in case there are no compounds involved we add the reaction node
            # individually
            self.add_node(back)
            for src in prods:
                self.add_edge(src, back)
            for tar in subs:
                self.add_edge(back, tar)
        for src in rxn.substrates:
            self._compounds.add(src)
            self.add_edge(src, rxn)
            if rxn.reversible and not self._split:
                self.add_edge(rxn, src)
        for tar in rxn.products:
            self._compounds.add(tar)
            self.add_edge(rxn, tar)
            if rxn.reversible and not self._split:
                self.add_edge(tar, rxn)

    def remove_reaction(self, rxn):
        assert isinstance(rxn, Reaction)
        self._reactions.remove(rxn)
        self.remove_node(rxn)

    def update_reactions(self):
        """
        Update the substrate, product, stoichiometry information of reactions
        according to changes that were done in the graph.
        """
        raise NotImplementedError

#    def __getattr__(self, name):
#        return type(self._network).__getattribute__(self._network, name)

class MetaboliteCentricNetwork(networkx.MultiDiGraph):
    """

    """

    _counter = 0

    def __init__(self, reactions=None, name='', labelled=True, * args, ** kwargs):
        self.__class__._counter += 1
        if name:
            self.name = str(name)
        else:
            self.name = "MetaboliteCentricNetwork-%d" % self.__class__._counter
        super(MetaboliteCentricNetwork, self).__init__(*args, ** kwargs)
        self.logger = logging.getLogger(\
                                        "pyMetabolism.MetaboliteCentricNetwork.%s" % self.name)
        self.handler = NullHandler
        self.logger.addHandler(self.handler)
        self.labelled = labelled
        if reactions:
            self._populate_graph_on_init(reactions)

    def _populate_graph_on_init(self, reactions):
        for rxn in reactions:
            self.add_reaction(rxn)

    def add_reaction(self, rxn):
        """docstring for add_reaction"""
        for sub in rxn.substrates:
            for prod in rxn.substrates:
                if sub == prod:
                    continue
                if self.labelled:
                    if rxn.reversible and not\
                        (self.has_edge(sub, prod, key=rxn.identifier) or\
                         self.has_edge(prod, sub, key=rxn.identifier)):
                        self.add_edge(sub, prod, key=rxn.identifier)
                        self.add_edge(prod, sub, key=rxn.identifier)
                    elif not self.has_edge(sub, prod, key=rxn.identifier):
                        self.add_edge(sub, prod, key=rxn.identifier)
                else:
                    if rxn.reversible:
                        self.add_edge(sub, prod)
                        self.add_edge(prod, sub)
                    else:
                        self.add_edge(sub, prod)


class ReactionCentricNetwork(networkx.MultiDiGraph):
    """

    """
    def __init__(self, reactions=None, name='', labelled=True, * args, ** kwargs):
        self.__class__._counter += 1
        if name:
            self.name = str(name)
        else:
            self.name = "ReactionCentricNetwork-%d" % self.__class__._counter
        super(MetaboliteCentricNetwork, self).__init__(*args, ** kwargs)
        self.logger = logging.getLogger(\
                                        "pyMetabolism.ReactionCentricNetwork.%s" % self.name)
        self.labelled = labelled
        if reactions:
            self._populate_graph_on_init(reactions)

    def _populate_graph_on_init(self, reactions):
        for rxn1 in reactions:
            for rxn2 in reactions:
                if not rxn1 == rxn2:
                    self.add_reaction_pair(rxn1, rxn2)

    def add_reaction_pair(self, rxn1, rxn2):
        """
        docstring for add_reaction
        """
        def connect_pair(self, rxn1, rxn2, compartment):
            if self.labelled:
                if not self.has_edge(rxn1, rxn2, key=compartment):
                    self.add_edge(rxn1, rxn2, key=compartment)
            else:
                self.add_edge(rxn1, rxn2)
        # find potential targets
        targets = set(rxn1.substrates + rxn1.products).intersection(
                                                                    set(rxn2.substrates + rxn2.products))
        # compartments!
        for compound in targets:
            if rxn1.compartments_dict[compound] != rxn2.compartments[compound]:
                continue
            rxn1_sub = compound in rxn1.substrates
            rxn2_sub = compound in rxn2.substrates
            if not rxn1_sub and rxn2_sub:
                connect_pair(rxn1, rxn2, rxn1.compartments_dict[compound])
            if not rxn2_sub and rxn1_sub:
                connect_pair(rxn1, rxn2, rxn1.compartments_dict[compound])
