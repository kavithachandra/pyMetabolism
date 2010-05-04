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
    reacs = list(graph.reactions)
    for rxn in reacs:
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
    logger = logging.getLogger("%s.random_metabolic_network"\
        % (options.main_logger_name))
    graph = BipartiteMetabolicNetwork(name="random model")
    for i in xrange(compounds):
        graph.add_compound(Compound("%s%d" % (options.compound_prefix, i)))
    # choose a number of reactions as reversible
    reversibles = random.sample(xrange(reactions), reversible)
    reversibles.sort()
    for i in xrange(reactions):
        if len(reversibles) > 0 and i == reversibles[0]:
            del reversibles[0]
            graph.add_reaction(Reaction("%s%d" % (options.reaction_prefix, i),\
                (), (), (), reversible=True))
        else:
            graph.add_reaction(Reaction("%s%d" % (options.reaction_prefix, i),\
                (), (), (), reversible=False))
    for src in graph.compounds:
        for tar in graph.reactions:
            if random.random() < p:
                logger.debug("Adding edge %s-%s.", src.identifier,\
                    tar.identifier)
                graph.add_edge(src, tar, factor=0)
            # a conditional case here (elif not if) because we cannot determine
            # substrates and products from bidirectional edges
            elif random.random() < p:
                logger.debug("Adding edge %s-%s.", tar.identifier,\
                    src.identifier)
                graph.add_edge(tar, src, factor=0)
    prune_network(graph)
    return graph

def scale_free_metabolic_network(compounds, reactions, reversible, m, n):
    """
    Uses a Barabasi-Alberts-like preferential attachment algorithm. Adopted from
    the networkx implementation.

    @param m: How many compounds a new reaction node links to

    @param n: How many reactions a new compound node links to
    """
    options = OptionsManager()
    logger = logging.getLogger("%s.scale_free_metabolic_network"\
        % (options.main_logger_name))
    graph = BipartiteMetabolicNetwork(name="scale-free model")
    # target nodes for reactions
    rxn_targets = []
    for i in xrange(m):
        comp = Compound("%s%d" % (options.compound_prefix, i))
        graph.add_compound(comp)
        rxn_targets.append(comp)
    # target nodes for compounds
    comp_targets = []
    # choose a number of reactions as reversible
    reversibles = random.sample(xrange(reactions), reversible)
    reversibles.sort()
    for i in xrange(n):
        if len(reversibles) > 0 and i == reversibles[0]:
            del reversibles[0]
            rxn = Reaction("%s%d" % (options.reaction_prefix, i),\
                (), (), (), reversible=True)
        else:
            rxn = Reaction("%s%d" % (options.reaction_prefix, i),\
                (), (), (), reversible=False)
        graph.add_reaction(rxn)
        comp_targets.append(rxn)
    # biased lists for preferential attachment
    repeated_cmpds = []
    repeated_rxns = []
    # current vertices being added
    current_rxn = n
    current_comp = m
    # to prevent reactions from consuming and producing the same compound
    new_targets = list()
    rm_targets = list()
    while (current_comp < compounds or current_rxn < reactions):
        if current_comp < compounds:
            source = Compound("%s%d" % (options.compound_prefix, current_comp))
            graph.add_compound(source)
            for rxn in comp_targets:
                if random.random() < 0.5:
                    graph.add_edge(source, rxn, factor=0)
                else:
                    graph.add_edge(rxn, source, factor=0)
            repeated_rxns.extend(comp_targets)
            repeated_cmpds.extend([source] * n)
            comp_targets = set()
            while len(comp_targets) < n:
                rxn = random.choice(repeated_rxns)
                comp_targets.add(rxn)
            current_comp += 1
        if current_rxn < reactions:
            if len(reversibles) > 0 and current_rxn == reversibles[0]:
                del reversibles[0]
                source = Reaction("%s%d" % (options.reaction_prefix,\
                    current_rxn), (), (), (), reversible=True)
            else:
                source = Reaction("%s%d" % (options.reaction_prefix,\
                    current_rxn), (), (), (), reversible=False)
            graph.add_reaction(source)
            new_targets = list()
            rm_targets = list()
            for comp in rxn_targets:
                # same compound may not be in substrates and products
                if (comp in graph.predecessors(source)) or (comp in graph.successors(source)):
                    cmpd = random.choice(repeated_cmpds)
                    while cmpd in graph.predecessors(source) or cmpd in graph.successors(source) or cmpd in new_targets or cmpd in rxn_targets:
                        cmpd = random.choice(repeated_cmpds)
                    new_targets.append(cmpd)
                    rm_targets.append(comp)
                    continue
                if random.random() < 0.5:
                    graph.add_edge(source, comp, factor=0)
                else:
                    graph.add_edge(comp, source, factor=0)
            for comp in rm_targets:
                rxn_targets.remove(comp)
            for comp in new_targets:
                rxn_targets.add(comp)
                if random.random() < 0.5:
                    graph.add_edge(source, comp, factor=0)
                else:
                    graph.add_edge(comp, source, factor=0)
            repeated_cmpds.extend(rxn_targets)
            repeated_rxns.extend([source] * m)
            rxn_targets = set()
            while len(rxn_targets) < m:
                comp = random.choice(repeated_cmpds)
                rxn_targets.add(comp)
            current_rxn += 1
    prune_network(graph)
    return graph

if __name__ == "__main__":
    from pyMetabolism.graph_algorithms import normed_in_out_degrees, plot_bipartite_network_log_degree_distribution, plot_bipartite_network_degree_distribution
    model = scale_free_metabolic_network(200, 300, 100, 3, 2)
#    model = random_metabolic_network(200, 300, 100, 0.2)
#    (metb_in, metb_out) = normed_in_out_degrees(model, model.compounds)
#    (rxn_in, rxn_out) = normed_in_out_degrees(model, model.reactions)
#    plot_bipartite_network_log_degree_distribution((metb_in, metb_out, rxn_in,\
#        rxn_out), "scale_free_degree_distribution.png", "directed scale-free bipartite model")
#    plot_bipartite_network_degree_distribution((metb_in, metb_out, rxn_in,\
#        rxn_out), "random_degree_distribution.png", "directed random bipartite model")