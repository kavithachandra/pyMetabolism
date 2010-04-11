# encoding: utf-8


"""
Multiple ways of creating a graph representation of a
L{Metabolism<pyMetabolism.metabolism.Metabolism>} instance.
Generation of a stoichiometric matrix is also possible.

@author: Nikolaus Sonnenschein
@author: Moritz Beber
@copyright: Jacobs University Bremen. All rights reserved.
@since: 2010-04-04
@attention: Graph instances and L{Metabolism<pyMetabolism.metabolism.Metabolism>}
            instances are not interlinked:
            Adding or removing reactions from one does not affect the other.
@todo: complete docs, many things
"""


#import networkx
import openopt
import logging
import random
import numpy
import matplotlib.pyplot


from pyMetabolism import OptionsManager
from pyMetabolism.graph_generation import random_metabolic_network
from pyMetabolism.stoichiometry import StoichiometricMatrix
from pyMetabolism.stoichiometry_algorithms import verify_consistency
from pyMetabolism.miscellaneous import lcm, fxrange
from pyMetabolism.exceptions import PyMetabolismError

def normed_in_out_degrees(graph, nodes):
    """
    @param graph: Metabolic Network
    """
    def normed_bincount(sequence):
        """
        """
        bins = numpy.bincount(sequence)
        total = float(sum(bins))
        normed = list()
        for val in bins:
            normed.append(float(val) / total)
        return normed
    
    in_deg_sq = list()
    for (node, deg) in graph.in_degree_iter(nodes):
        in_deg_sq.append(deg)
    in_deg_sq = normed_bincount(in_deg_sq)
    out_deg_sq = list()
    for (node, deg) in graph.out_degree_iter(nodes):
        out_deg_sq.append(deg)
    out_deg_sq = normed_bincount(out_deg_sq)
    return (in_deg_sq, out_deg_sq)

def plot_bipartite_network_degree_distribution((metb_in_deg, metb_out_deg,\
    rxn_in_deg, rxn_out_deg), filename, title):
    """
    """
    line1 = matplotlib.pyplot.plot(metb_in_deg, marker='.', color='red',\
        linestyle='-')
    line2 = matplotlib.pyplot.plot(metb_out_deg, marker='x', color='red',\
        linestyle='-.', linewidth=2)
    line3 = matplotlib.pyplot.plot(rxn_in_deg, marker='.', color='blue',\
        linestyle='-')
    line4 = matplotlib.pyplot.plot(rxn_out_deg, marker='x', color='blue',\
        linestyle='-.', linewidth=2)
    matplotlib.pyplot.figlegend((line1, line2, line3, line4), ("metbs in-degree",
        "metbs out-degree", "rxns in-degree", "rxns out-degree"), 1)
    matplotlib.pyplot.title(title)
    matplotlib.pyplot.ylabel("P(k)")
    matplotlib.pyplot.xlabel("k")
    matplotlib.pyplot.savefig(filename)
    matplotlib.pyplot.show()

def plot_bipartite_network_log_degree_distribution((metb_in_deg, metb_out_deg,\
    rxn_in_deg, rxn_out_deg), filename, title):
    """
    """
    line1 = matplotlib.pyplot.loglog(metb_in_deg, marker='.', color='red',\
        linestyle='-')
    line2 = matplotlib.pyplot.loglog(metb_out_deg, marker='x', color='red',\
        linestyle='-.', linewidth=2)
    line3 = matplotlib.pyplot.loglog(rxn_in_deg, marker='.', color='blue',\
        linestyle='-')
    line4 = matplotlib.pyplot.loglog(rxn_out_deg, marker='x', color='blue',\
        linestyle='-.', linewidth=2)
    matplotlib.pyplot.figlegend((line1, line2, line3, line4), ("metbs in-degree",
        "metbs out-degree", "rxns in-degree", "rxns out-degree"), 1)
    matplotlib.pyplot.title(title)
    matplotlib.pyplot.ylabel("log[P(k)]")
    matplotlib.pyplot.xlabel("log[k]")
    matplotlib.pyplot.savefig(filename)
    matplotlib.pyplot.show()

def randomly_make_consistent_stoichiometry(graph, factors):
    """
    @param graph: L{BipartiteMetabolicNetwork
        <pyMetabolism.graph_view.BipartiteMetabolicNetwork>}

    @param factors: An iterable of desired stoichiometric factors. May contain
        a specific factor multiple times to bias towards selection of that factor.

    @rtype: L{BipartiteMetabolicNetwork
        <pyMetabolism.graph_view.BipartiteMetabolicNetwork>} with consistent
        stoichiometric factors.
    """
    options = OptionsManager()
    def populate_network(graph, factors):
        for e in graph.edges_iter():
            graph.add_edge(e[0], e[1], factor=random.choice(factors))
    count = 0
    while (True):
        count += 1
        options.logger.info("Attempt number %d", count)
        populate_network(graph, factors)
        matrix = StoichiometricMatrix()
        matrix.make_new_from_network(graph)
        (val, vector) = verify_consistency(matrix)
        if val:
            break

def make_consistent_stoichiometry(graph, factors):
    """
    @param graph: L{BipartiteMetabolicNetwork
        <pyMetabolism.graph_view.BipartiteMetabolicNetwork>}

    @param factors: An iterable of desired stoichiometric factors.

    @rtype: L{BipartiteMetabolicNetwork
        <pyMetabolism.graph_view.BipartiteMetabolicNetwork>} with consistent
        stoichiometric factors.
    """
    options = OptionsManager()
#    step = 1. / float(reduce(lcm, factors) * len(graph.compounds))
    step = reduce(lcm, factors) * len(graph.compounds)
    # generate mass vector for compounds
    mass_vector = dict()
#    masses = list(fxrange(step, 1., step))
    masses = range(1, step)
    for comp in graph.compounds:
        mass = random.choice(masses)
        mass_vector[comp] = mass
        masses.remove(mass)
    options.logger.debug(str(mass_vector))
    upper = max(factors)

    def balance_reaction(reaction):
        """
        balance individual reaction
        """
        in_deg = graph.in_degree(reaction)
        out_deg = graph.out_degree(reaction)
        objective = numpy.ones(in_deg + out_deg)
        A_eq = numpy.empty(in_deg + out_deg)
        # substrate mass_vector
        msg = "%s:\n" % reaction.identifier
        for (i, comp) in enumerate(graph.predecessors(reaction)):
            A_eq[i] = -mass_vector[comp]
            msg += " -%s" % str(mass_vector[comp])
        # product mass_vector
        for (i, comp) in enumerate(graph.successors(reaction)):
            A_eq[i + in_deg] = mass_vector[comp]
            msg += " +%s" % str(mass_vector[comp])
        msg += " = 0"
        options.logger.debug(msg)
        b_eq = numpy.zeros(1)
        lb = numpy.ones(in_deg + out_deg)
        ub = numpy.empty(in_deg + out_deg)
        for i in xrange(in_deg + out_deg):
            ub[i] = upper
        problem = openopt.LP(f=objective, Aeq=A_eq, beq=b_eq, lb=lb, ub=ub,\
            intVars=range(in_deg + out_deg))
        result = problem.solve(options.solver, iprint=-5)
        if not result.isFeasible:
            raise PyMetabolismError("Reaction %s cannot be balanced with the"\
                " given mass vector." % reaction)
        msg = "%s:\n" % reaction.identifier
        for (i, comp) in enumerate(graph.predecessors(reaction)):
            graph[comp][reaction]["factor"] =  result.xf[i]
            msg += " %s +" % result.xf[i]
        msg += " ->"
        for (i, comp) in enumerate(graph.successors(reaction)):
            graph[reaction][comp]["factor"] =  result.xf[in_deg + i]
            msg += " %s +" % result.xf[in_deg + i]
        options.logger.debug(msg)
        options.logger.info(result.xf)

    total = float(len(graph.reactions))
    for (i, rxn) in enumerate(graph.reactions):
        balance_reaction(rxn)
        options.logger.info("%.2f %% complete.", float(i) / total * 100.)


if __name__ == "__main__":
    options = OptionsManager()
    options.logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    options.logger.addHandler(handler)
    model = random_metabolic_network(200, 300, 100, 0.2)
#    (metb_in, metb_out) = normed_in_out_degrees(model, model.compounds)
#    (rxn_in, rxn_out) = normed_in_out_degrees(model, model.reactions)
#    plot_bipartite_network_degree_distribution((metb_in, metb_out, rxn_in, rxn_out),
#        "degree_distribution.png", "directed random bipartite model")
    make_consistent_stoichiometry(model, range(1,7))
    matrix = StoichiometricMatrix()
    matrix.make_new_from_network(model)
    options.logger.info(matrix)
    options.logger.info(str(verify_consistency(matrix)))
