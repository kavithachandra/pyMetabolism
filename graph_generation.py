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


import random
import logging
from pyMetabolism import OptionsManager
from pyMetabolism.metabolism import Compound, Reaction
from pyMetabolism.graph_view import BipartiteMetabolicNetwork


def prune_network(graph):
    """
    Removes stub reactions (in- or out-degree of 1 and 0 respectively) and
    assures that all other reactions consume and produce something.

    @param graph: L{BipartiteMetabolicNetwork
        <pyMetabolism.graph_view.BipartiteMetabolicNetwork>}

    @rtype: L{BipartiteMetabolicNetwork
        <pyMetabolism.graph_view.BipartiteMetabolicNetwork>}


    """
    for rxn in graph.reactions:
        in_deg = graph.in_degree(rxn)
        out_deg = graph.out_degree(rxn)
        if in_deg == 0:
            if out_deg <= 1:
                graph.remove_reaction(rxn)
            else:
                targets = graph.successors(rxn)
                flips = random.randint(1, len(targets) - 1)
                while (flips > 0):
                    target = random.choice(targets)
                    factor = graph[rxn][target]["factor"]
                    graph.remove_edge(rxn, target)
                    graph.add_edge(target, rxn, factor=factor)
                    targets.remove(target)
                    flips -= 1
        elif out_deg == 0:
            if in_deg <= 1:
                graph.remove_reaction(rxn)
            else:
                targets = graph.predecessors(rxn)
                flips = random.randint(1, len(targets) - 1)
                while (flips > 0):
                    target = random.choice(targets)
                    factor = graph[target][rxn]["factor"]
                    graph.remove_edge(target, rxn)
                    graph.add_edge(rxn, target, factor=factor)
                    targets.remove(target)
                    flips -= 1


def random_metabolic_network(compounds, reactions, reversible, p, seed=None):
    """
    Creates a bipartite graph that models metabolism according to the principles
    of an Erdos-Renyi-like random graph.

    @param compounds: Integer specifying the number of compounds.

    @param reactions: Integer, specifying the number of reactions.

    @param p: Probability for an edge to be set (on average this is the density
        of the graph).

    @rtype: L{BipartiteMetabolicNetwork
        <pyMetabolism.graph_view.BipartiteMetabolicNetwork>}

    @raise ValueError: If either compounds or reactions are not integers.
    """
    if seed:
        random.seed(seed)
    options = OptionsManager()
    logger = logging.getLogger("%s.random_metabolic_graph"\
        % (options.main_logger_name))
    graph = BipartiteMetabolicNetwork(name="random model")
    for i in xrange(compounds):
        graph.add_compound(Compound("%s%d" % (options.compound_prefix, i)))
    for i in xrange(reactions):
        graph.add_reaction(Reaction("%s%d" % (options.reaction_prefix, i), (),\
            (), ()))
    for src in graph.compounds:
        for tar in graph.reactions:
            if random.random() < p:
                logger.debug("Adding edge %s-%s.", src.identifier,\
                    tar.identifier)
                graph.add_edge(src, tar, factor=1)
            # a conditional case here (elif not if) because we cannot determine
            # substrates and products from bidirectional edges
            elif random.random() < p:
                logger.debug("Adding edge %s-%s.", tar.identifier,\
                    src.identifier)
                graph.add_edge(tar, src, factor=1)
    prune_network(graph)
    reactions = list(graph.reactions)
    reversibles = list()
    for i in xrange(reversible):
        rxn = random.choice(reactions)
        reversibles.append(rxn)
        reactions.remove(rxn)
    for rxn in reversibles:
        rev = Reaction(rxn.identifier + options.rev_reaction_suffix, (), (), ())
        graph.add_reaction(rev)
        for cmpd in graph.predecessors(rxn):
            graph.add_edge(rev, cmpd, factor=1)
        for cmpd in graph.successors(rxn):
            graph.add_edge(cmpd, rev, factor=1)
    return graph

def scale_free_metabolic_network(compounds, reactions):
    """
    Uses a Barabasi-Alberts-like preferential attachment algorithm.
    """
    pass