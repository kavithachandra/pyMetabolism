"""
Multiple ways of creating a graph representation of a L{Metabolism} instance.
Generation of a stoichiometric matrix is also possible.

@author: Nikolaus Sonnenschein
@author: Moritz Beber
@contact: niko(at)foo.org
@contact: moritz(at)foo.org
@copyright: Jacobs University Bremen. All rights reserved.
@since: 2009-09-11
@attention: Graph instances and L{Metabolism} instances are not interlinked:
            Adding or removing reactions from one does not affect the other.
@todo: set up logging, complete docs
"""


import random
import networkx
from numpy import hstack, vstack, zeros
from pyMetabolism.metabolism import Compound, Reaction


class BipartiteMetabolicNetwork(networkx.DiGraph):
    """
A bipartite graph representation of a L{Metabolism} instance.

@todo: many things
    """
    
    def __init__(self, reactions=None, name='', weighted=True, node_info=False,
        *args, **kwargs):
        """
@param reactions: A list of L{Reaction}s from which the graph can be populated.
@type reactions: C{list}

@param name: A simple string to identify the graph if desired.
@type name: C{str}

@param weighted: By default the stoichiometric factors are used as edge weights,
                this can be turned off.
@type weighted: C{bool}

@param node_info: If true, information will be added to the graph's nodes.
                 For compound nodes the following information will be tried in
                 order:
                 
                    1. L{Compound.charge}
                    2. L{Compound.mass}
                 
                 For reaction nodes L{Reaction.rate} will be used.
@type node_info: C{bool}
        """
        super(BipartiteMetabolicNetwork, self).__init__(name=str(name), *args,
            **kwargs)
        if reactions:
            self._populate_graph_on_init(reactions, weighted, node_info)
    
    def _populate_graph_on_init(self, reactions, weighted, node_info):
        for rxn in reactions:
            self.add_reaction(rxn)
    
    def add_reaction(self, rxn, weighted=True, node_info=False):
        """
Add edges to the bipartite representation of L{Metabolism} as necessary to
represent a reaction.

@param rxn: Instance to be added.
@type rxn: L{Reaction}

@param weighted: By default the stoichiometric factors are used as edge weights,
                this can be turned off.
@type weighted: C{bool}

@param node_info: If true, information will be added to the graph's nodes.
                 For compound nodes the following information will be tried in
                 order:
                 
                    1. L{Compound.charge}
                    2. L{Compound.mass}
                 
                 For reaction nodes L{Reaction.rate} will be used.
@type node_info: C{bool}
@note: C{node_info} is purely aimed at graph output to a format that supports
       node attributes.
        """
        if node_info and rxn.rate:
            self.add_node(rxn, rate=rxn.rate)
        for compound in rxn.substrates:
            if weighted:
                self.add_edge(compound, rxn, rxn.stoichiometry_dict[compound])
            else:
                self.add_edge(compound, rxn)
            if node_info and compound.charge:
                self.add_node(compound, charge=compound.charge)
            if node_info and compound.mass:
                self.add_node(compound, mass=compound.mass)
        for compound in rxn.products:
            if weighted:
                self.add_edge(rxn, compound, rxn.stoichiometry_dict[compound])
            else:
                self.add_edge(rxn, compound)
            if node_info and compound.charge:
                self.add_node(compound, charge=compound.charge)
            if node_info and compound.mass:
                self.add_node(compound, mass=compound.mass)
        if rxn.reversible:
            for compound in rxn.substrates:
                if weighted:
                    self.add_edge(rxn, compound,
                        rxn.stoichiometry_dict[compound])
                else:
                    self.add_edge(rxn, compound)
            for compound in rxn.products:
                if weighted:
                    self.add_edge(compound, rxn,
                        rxn.stoichiometry_dict[compound])
                else:
                    self.add_edge(compound, rxn)
    
    def _double_edges(self, graph):
        """
Finds all bi-directional edges in the C{graph}.

@param graph: Directed graph instance.
@type graph: L{BipartiteMetabolicNetwork}

@return: A list of all bi-directional edges, contains each edge in both
         directions.
@rtype: C{list}
        """
        double = list()
        for edge in graph.edges_iter():
            if graph.has_edge(edge[1], edge[0]):
                double.append(edge)
        return double
    
    def _check_conditions(self, graph, first, second):
        """
Checks if all necessary conditions to flip an edge C{first} with another
C{second} are fulfilled.

@param graph: Directed graph instance.
@type graph: L{BipartiteMetabolicNetwork}

@param first: Directed edge.
@type first: Two C{tuple}

@param second: Directed edge.
@type second: Two C{tuple}

@rtype: C{bool}

        """
        # check if we would flip to existing edges
        if graph.has_edge(first[0], second[1]):
            return False
        if graph.has_edge(second[0], first[1]):
            return False
        # check if we would introduce double edges
        if graph.has_edge(second[1], first[0]):
            return False
        if graph.has_edge(first[1], second[0]):
            return False
        return True
    
    def _switch_double_edge(self, graph, first):
        """
Switching of a double edge is treated separately from a unidirectional edge.

@param graph: Directed graph instance.
@type graph: L{BipartiteMetabolicNetwork}

@param first: Directed edge, known to be bidirectional.
@type first: Two C{tuple}

@rtype: C{bool}
        """
        # find another bi-directional edge
        double_edges = _double_edges(graph)
        # not enough bi-directional edges in the graph
        if len(double_edges) < 2:
            return False
        second = random.choice(double_edges)
        if first[0] in second:
            return False
        if first[1] in second:
            return False
        if _check_conditions(graph, first, second):
            # add new edges
            graph.add_edge(first[0], second[1])
            graph.add_edge(second[1], first[0])
            graph.add_edge(first[1], second[0])
            graph.add_edge(second[0], first[1])
            # remove old edges
            graph.remove_edge(first[0], first[1])
            graph.remove_edge(first[1], first[0])
            graph.remove_edge(second[0], second[1])
            graph.remove_edge(second[1], second[0])
            return True
        else:
            return False
    
    def _select_edges(self, graph):
        """
Select two edges whose sources lie in the same node set.

@param graph: Directed graph instance.
@type graph: L{BipartiteMetabolicNetwork}

@return: Returns a pair of edges.
@rtype: C{tuple}
        """
        # source is a metabolite
        edges = graph.edges()
        first = random.choice(edges)
        second = random.choice(edges)
        if isinstance(first[0], Compound):
            while not isinstance(second[0], Compound):
                second = random.choice(edges)
        else:
            while not isinstance(second[0], Reaction):
                second = random.choice(edges)
        return (first, second)
    
    def _flip_bipartite_edge(self, graph):
        """
Switching of two selected edges according to a set of rules.

@param graph: Directed graph instance.
@type graph: L{BipartiteMetabolicNetwork}

@return: Whether edge flipping was successful.
@rtype: C{bool}
        """
        (first, second) = _select_edges(graph, metbs, rxns)
        # same edges were picked by random algorithm
        if first[0] in second:
            return False
        if first[1] in second:
            return False
        # first edge is bi-directional
        if graph.has_edge(first[1], first[0]):
            return _switch_double_edge(graph, first)
        if graph.has_edge(second[1], second[0]):
            return False
        # check other necessary conditions
        if _check_conditions(graph, first, second):
            # print "First edge = %s -> %s" % (first[0], first[1])
            # print "Second edge = %s -> %s" % (second[0], second[1])
            # print "Adding %s -> %s" % (first[0], second[1])
            # print "Adding %s -> %s" % (second[0], first[1])
            graph.add_edge(first[0], second[1])
            graph.add_edge(second[0], first[1])
            graph.remove_edge(first[0], first[1])
            graph.remove_edge(second[0], second[1])
            return True
        else:
            return False
    
    def randomise(self, flip=100):
        """
Produces a randomised version of the bipartite graph. This is achieved by
switching edges in the graph a number of
times equal to the number of edges present times C{flip}. This function is
intended as a basis for calculating three-node subgraph statistics in metabolic
networks. As such only degrees of nodes, bi-directional link properties, etc.
(see:
R Milo, S Shen-Orr, S Itzkovitz, N Kashtan, D Chklovskii & U Alon,
Network Motifs: Simple Building Blocks of Complex Networks
Science, 298:824-827 (2002)
) are conserved. For larger subgraph statistics also the smaller subgraph
statistics would have to be conserved. The novelty of this function is to
preserve the bipartite nature of the graph.

@param flip: How often each edge should be switched ideally.
@type flip: C{int}

@return: A randomised version of the bipartite graph.
@rtype: L{BipartiteMetabolicNetwork}
        """
        rnd_graph = self.copy()
        num_rnd = self.size() * flip
        initial_double = len(_double_edges(rnd_graph))
        print "Number of bi-directional edges in the graph: %d"\
            % initial_double
        success = 0
        for i in xrange(num_rnd):
            # keep attempting to switch edges until right conditions are met
            if _flip_bipartite_edge(rnd_graph):
                success += 1
        end_double = len(_double_edges(rnd_graph))
        print "Number of bi-directional edges in the graph after randomisation: %d"\
            % end_double
        if (initial_double - end_double) != 0:
            raise nx.NetworkXError
        print "Flip success rate: %f" % (float(success) / float(num_rnd))
        return rnd_graph


class MetaboliteCentricNetwork(networkx.MultiDiGraph):
    """"""
    def __init__(self, reactions=None, name='', weighted=True, *args, **kwargs):
        super(MetaboliteCentricNetwork, self).__init__(name=str(name), *args, **kwargs)
        if reactions:
            self.__populate_graph_on_init(reactions)
    
    def __populate_graph_on_init(self, reactions):
        for reac in reactions:
            self.add_reaction(reac)
    
    def add_reaction(self, reaction):
        """docstring for add_reaction"""
        edges = [(substr, prod, reaction) for prod in reaction.products for substr in reaction.substrates]
        if reaction.reversible:
            self.add_edges_from(edges + [(j, i, d) for i, j, d in edges])
        else:
            self.add_edges_from(edges)

