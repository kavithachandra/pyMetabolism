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
from pyMetabolism.metabolism_logging import NullHandler
#from pyMetabolism.metabolism import Compound, Reaction, Compartment, Metabolism


class BipartiteMetabolicNetwork(networkx.DiGraph):
    """
    A bipartite graph representation of a L{Metabolism<pyMetabolism.metabolism.Metabolism>}
    instance.

    @todo: many things
    """

    _counter = 0

    def __init__(self, reactions=None, name='', weighted=True, node_info=False,
                 * args, ** kwargs):
        """
        @param reactions: A list of L{Reaction<pyMetabolism.metabolism.Reaction>}s from
                          which the graph can be populated.
        @type reactions: C{list}

        @param name: A simple string to identify the graph if desired.
        @type name: C{str}

        @param weighted: By default the stoichiometric factors are used as edge
            weights, this can be turned off.
        @type weighted: C{bool}

        @param node_info: If true, information will be added to the graph's
            nodes.
            For compound nodes the following information will be tried in
            order:

                1. L{Compound.charge<pyMetabolism.metabolism.Compound.charge>}
                2. L{Compound.mass<pyMetabolism.metabolism.Compound.mass>}

            For reaction nodes
            L{Reaction.rate_constant<pyMetabolism.metabolism.Reaction.rate_constant>} will
            be used.
        @type node_info: C{bool}
        """
        self.__class__._counter += 1
        if name:
            self.name = str(name)
        else:
            self.name = "BipartiteMetabolicNetwork-%d" % self.__class__._counter
        super(BipartiteMetabolicNetwork, self).__init__(*args, ** kwargs)
        self.logger = logging.getLogger(\
                                        "pyMetabolism.BipartiteMetabolicNetwork.%s" % self.name)
        self.handler = NullHandler
        self.logger.addHandler(self.handler)
        self.weighted = weighted
        if reactions:
            self._populate_graph_on_init(reactions, node_info)

    def _populate_graph_on_init(self, reactions, node_info):
        for rxn in reactions:
            self.add_reaction(rxn, node_info)

    def add_reaction(self, rxn, node_info=False):
        """
        Add edges to the bipartite representation of
        L{Metabolism<pyMetabolism.metabolism.Metabolism>} as necessary
        to represent a reaction.

        @param rxn: Instance to be added.
        @type rxn: L{Reaction<pyMetabolism.metabolism.Reaction>}

        @param node_info: If true, information will be added to the graph's
            nodes.
            For compound nodes the following information will be tried in
            order and are availabe as charge and mass attribute respectively.

                1. L{Compound.charge<pyMetabolism.metabolism.Compound.charge>}
                2. L{Compound.mass<pyMetabolism.metabolism.Compound.mass>}

            For reaction nodes
            L{Reaction.rate_constant<pyMetabolism.metabolism.Reaction.rate_constant>}
            will be used.
        @type node_info: C{bool}
        @note: C{node_info} is purely aimed at graph output to a format that supports
            node attributes.
        """
        if node_info and rxn.rate_constant:
            self.add_node(rxn, rate_constant=rxn.rate_constant)
        for compound in rxn.substrates:
            if node_info and compound.charge:
                self.add_node(compound, charge=compound.charge)
            if node_info and compound.mass:
                self.add_node(compound, mass=compound.mass)
            if self.weighted:
                self.add_edge(compound, rxn, \
                              compartment=rxn.compartments_dict[compound], \
                              weight=rxn.stoichiometry_dict[compound])
            else:
                self.add_edge(compound, rxn,
                              compartment=rxn.compartments_dict[compound])
        for compound in rxn.products:
            if node_info and compound.charge:
                self.add_node(compound, charge=compound.charge)
            if node_info and compound.mass:
                self.add_node(compound, mass=compound.mass)
            if self.weighted:
                self.add_edge(rxn, compound, \
                              compartment=rxn.compartments_dict[compound], \
                              weight=rxn.stoichiometry_dict[compound])
            else:
                self.add_edge(rxn, compound,
                              compartment=rxn.compartments_dict[compound])
        if rxn.reversible:
            for compound in rxn.substrates:
                if self.weighted:
                    self.add_edge(rxn, compound, \
                                  compartment=rxn.compartments_dict[compound], \
                                  weight=abs(rxn.stoichiometry_dict[compound]))
                else:
                    self.add_edge(rxn, compound,
                                  compartment=rxn.compartments_dict[compound])
            for compound in rxn.products:
                if self.weighted:
                    self.add_edge(compound, rxn, \
                                  compartment=rxn.compartments_dict[compound], \
                                  weight=-rxn.stoichiometry_dict[compound])
                else:
                    self.add_edge(compound, rxn,
                                  compartment=rxn.compartments_dict[compound])


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
        self.handler = NullHandler
        self.logger.addHandler(self.handler)
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
