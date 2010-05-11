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
@todo: complete docs, many things, add requirements to algorithms
"""


import openopt
import logging
import random
import numpy
import matplotlib.pyplot
from networkx import cycle_basis


from pyMetabolism import OptionsManager
from pyMetabolism.metabolism import Reaction
from pyMetabolism.stoichiometry import StoichiometricMatrix
from pyMetabolism.stoichiometry_algorithms import verify_consistency
from pyMetabolism.metabolism_exceptions import PyMetabolismError
from pyMetabolism.graph_view import BipartiteMetabolicNetwork


def discover_partners(graph, node, discovered, complete):
    """
    @param graph: L{BipartiteMetabolicNetwork}

    @param node: vertex whose neighbours to discover

    @param discovered: C{set} of already discovered reactions

    @param complete: C{set} of already completed reactions
    """
    for suc in graph.successors(node):
        if suc not in complete:
            discovered.add(suc)
    for pre in graph.predecessors(node):
        if pre not in complete:
            discovered.add(pre)

def balance_reaction_by_mass(graph, factors, reaction, mass_vector, upper,\
    logger):
    """
    All compounds have fixed masses and only stoichiometric factors can be varied.

    balance individual reaction
    """
    options = OptionsManager()
    in_deg = graph.in_degree(reaction)
    out_deg = graph.out_degree(reaction)
    # test different objective functions
    # only zeros leads to fast solutions (no objective function) that are
    # close to the upper boundary
    # all equal leads to slower solutions that are mostly one with a few
    # variations to balance the reactions
    # an objective function with entries of 1 / factor leads to very long
    # solution times but much more varying entries
    # best approach so far: no objective function, but use a starting point
    # by randomly picking coefficients from 'factors'
    objective = numpy.ones(in_deg + out_deg)
#        objective = numpy.ones(in_deg + out_deg)
#        for (i, val) in enumerate(objective):
#            objective[i] = val / float(random.choice(factors))
    start = numpy.empty(in_deg + out_deg)
    for i in xrange(len(start)):
        start[i] = float(random.choice(factors))
    A_eq = numpy.empty(in_deg + out_deg)
    # substrate mass_vector
    msg = "%s:\n" % reaction.identifier
    for (i, comp) in enumerate(graph.predecessors(reaction)):
        A_eq[i] = -mass_vector[comp]
        msg += " -%sS%d" % (str(mass_vector[comp]), i)
    # product mass_vector
    for (i, comp) in enumerate(graph.successors(reaction)):
        A_eq[i + in_deg] = mass_vector[comp]
        msg += " +%sS%d" % (str(mass_vector[comp]), i + in_deg)
    msg += " = 0"
    logger.debug(msg)
    b_eq = numpy.zeros(1)
    lb = numpy.ones(in_deg + out_deg)
    ub = numpy.empty(in_deg + out_deg)
    for i in xrange(in_deg + out_deg):
        ub[i] = upper
    problem = openopt.MILP(f=objective, x0=start, Aeq=A_eq, beq=b_eq, lb=lb,\
        ub=ub, intVars=range(in_deg + out_deg))
    result = problem.solve(options.solver, iprint=-10)
    if not result.isFeasible:
        raise PyMetabolismError("Reaction '%s' cannot be balanced with the"\
            " given mass vector." % reaction.identifier)
    for (i, comp) in enumerate(graph.predecessors(reaction)):
        graph[comp][reaction]["factor"] =  result.xf[i]
    for (i, comp) in enumerate(graph.successors(reaction)):
        graph[reaction][comp]["factor"] =  result.xf[in_deg + i]
    logger.debug(result.xf)

def balance_reaction_by_factors(graph, factors, reaction, mass_vector, upper, lower, logger):
    """
    Only part of compound masses and/or factors are fixed. We set the other
    factors by drawing from the distribution and then choosing appropriate masses.
    """
    options = OptionsManager()
    in_deg = graph.in_degree(reaction)
    out_deg = graph.out_degree(reaction)
    # objective function
    # try with factors according to degree
    objective = numpy.ones(in_deg + out_deg)
    # initial conditions
#    start = numpy.ones(in_deg + out_deg)
    # equality constraints
    A_eq = numpy.empty(in_deg + out_deg)
    # substrate coefficients
    msg = "%s:\n" % reaction.identifier
    for (i, comp) in enumerate(graph.predecessors(reaction)):
        if graph[comp][reaction]["factor"] == 0:
            graph[comp][reaction]["factor"] = random.choice(factors)
        A_eq[i] = -graph[comp][reaction]["factor"]
        msg += " -%sM%d" % (str(graph[comp][reaction]["factor"]), i)
    # product coefficients
    for (i, comp) in enumerate(graph.successors(reaction)):
        if graph[reaction][comp]["factor"] == 0:
            graph[reaction][comp]["factor"] = random.choice(factors)
        A_eq[i + in_deg] = graph[reaction][comp]["factor"]
        msg += " +%sM%d" % (str(graph[reaction][comp]["factor"]), i + in_deg)
    msg += " = 0"
    logger.info(msg)
    b_eq = numpy.zeros(1)
    # starting values
    start = numpy.empty(in_deg + out_deg)
    # lower boundaries
    lb = numpy.zeros(in_deg + out_deg)
    for (i, comp) in enumerate(graph.predecessors(reaction)):
        objective[i] = graph.in_degree(comp) + graph.out_degree(comp)
        if mass_vector[comp] != 0.:
            lb[i] = mass_vector[comp]
            start[i] = mass_vector[comp]
        else:
            lb[i] = lower
            start[i] = 1. / float(graph.in_degree(comp) + graph.out_degree(comp))
    for (i, comp) in enumerate(graph.successors(reaction)):
        objective[i + in_deg] = graph.in_degree(comp) + graph.out_degree(comp)
        if mass_vector[comp] != 0.:
            lb[i + in_deg] = mass_vector[comp]
            start[i + in_deg] = mass_vector[comp]
        else:
            lb[i + in_deg] = lower
            start[i + in_deg] = 1. / float(graph.in_degree(comp) + graph.out_degree(comp))
    # upper boundaries
    ub = numpy.empty(in_deg + out_deg)
    for (i, comp) in enumerate(graph.predecessors(reaction)):
        if mass_vector[comp] != 0.:
            ub[i] = mass_vector[comp]
        else:
            ub[i] = upper
    for (i, comp) in enumerate(graph.successors(reaction)):
        if mass_vector[comp] != 0.:
            ub[i + in_deg] = mass_vector[comp]
        else:
            ub[i + in_deg] = upper
    # problem
    problem = openopt.LP(f=objective, Aeq=A_eq, beq=b_eq, lb=lb, ub=ub)
    result = problem.solve(options.solver, iprint=-10)
    if not result.isFeasible:
        raise PyMetabolismError("Reaction '%s' cannot be balanced with the"\
            " given stoichiometric coefficients and pre-existing masses."\
            % reaction.identifier)
    for (i, comp) in enumerate(graph.predecessors(reaction)):
        mass_vector[comp] =  result.xf[i]
    for (i, comp) in enumerate(graph.successors(reaction)):
        mass_vector[comp] =  result.xf[in_deg + i]
    logger.info(result.xf)

def balance_multiple_reactions(graph, factors, rxns, mass_vector, lower, upper, logger):
    options = OptionsManager()
    sub = BipartiteMetabolicNetwork()
    for rxn in rxns:
        msg = "%s:\n" % rxn.identifier
        sub.add_reaction(rxn)
        for cmpd in graph.predecessors(rxn):
            coefficient = random.choice(factors)
            msg += " -%d%s" % (coefficient, cmpd)
            sub.add_compound(cmpd)
            sub.add_edge(cmpd, rxn, factor=coefficient)
            graph[cmpd][rxn]["factor"] = coefficient
        for cmpd in graph.successors(rxn):
            coefficient = random.choice(factors)
            msg += " +%d%s" % (coefficient, cmpd)
            sub.add_compound(cmpd)
            sub.add_edge(rxn, cmpd, factor=coefficient)
            graph[rxn][cmpd]["factor"] = coefficient
        msg += " = 0\n"
        logger.debug(msg)
    matrix = StoichiometricMatrix()
    matrix.make_new_from_network(sub)
    # objective function
    objective = numpy.ones(matrix.num_compounds)
    # equality constraints
    A_eq = numpy.array(numpy.transpose(matrix.matrix), dtype=float, copy=False)
    b_eq = numpy.zeros(matrix.num_reactions)
    # lower boundary
    lb = numpy.empty(matrix.num_compounds)
    for (i, comp) in enumerate(matrix.compounds):
        if mass_vector[comp] != 0.:
            lb[i] = mass_vector[comp]
        else:
            lb[i] = lower
    # upper boundary
    ub = numpy.empty(matrix.num_compounds)
    for (i, comp) in enumerate(matrix.compounds):
        if mass_vector[comp] != 0.:
            ub[i] = mass_vector[comp]
        else:
            ub[i] = upper
    # starting guess
    start = numpy.empty(matrix.num_compounds)
    logger.debug(A_eq)
    for (i, comp) in enumerate(matrix.compounds):
        start[i] = 1. / float(sub.in_degree(comp) + sub.out_degree(comp))
    problem = openopt.LP(f=objective, Aeq=A_eq, beq=b_eq, lb=lb, ub=ub, x0=start)
    result = problem.solve(options.solver, iprint=-10)
    if not result.isFeasible:
        raise PyMetabolismError("Unable to balance reaction system.")
    for (i, comp) in enumerate(matrix.compounds):
        if mass_vector[comp] == 0.:
            mass_vector[comp] = result.xf[i]
    logger.debug(result.xf)



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

def plot_bipartite_network_log_degree_distribution((metb_in_deg,\
    metb_out_deg, rxn_in_deg, rxn_out_deg), filename, title):
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

def allot_consistent_stoichiometry(graph, factors):
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
    logger = logging.getLogger("%s.allot_consistent_stoichiometry"\
        % options.main_logger_name)
    def populate_network(graph, factors):
        for e in graph.edges_iter():
            graph[e[0]][e[1]]["factor"] = random.choice(factors)
    count = 0
    while (True):
        count += 1
        logger.info("Attempt number %d", count)
        populate_network(graph, factors)
        matrix = StoichiometricMatrix()
        matrix.make_new_from_network(graph)
        (yes, result) = verify_consistency(matrix)
        if yes:
            logger.info("Success!")
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
    logger = logging.getLogger("%s.make_consistent_stoichiometry"\
        % options.main_logger_name)
    # generate mass vector for compounds
    mass_vector = dict()
    degrees = list()
    for comp in graph.compounds:
        degrees.append(graph.in_degree(comp) + graph.out_degree(comp))
    degrees.sort()
    for comp in graph.compounds:
#        mass_vector[comp] = round(abs(random.gauss(len(graph.compounds),\
#            len(graph.compounds) / 2)))
        mass_vector[comp] = degrees[-degrees.index(graph.in_degree(comp) + graph.out_degree(comp))]
    logger.debug(str(mass_vector))
    total = float(len(graph.reactions))
    maxi = max(factors)
    for (i, rxn) in enumerate(graph.reactions):
        balance_reaction_by_mass(graph, factors, rxn, mass_vector, maxi, logger)
        logger.info("%.2f %% complete.", float(i + 1) / total * 100.)

def grow_consistent_stoichiometry(graph, factors):
    """

    """
    options = OptionsManager()
    logger = logging.getLogger("%s.grow_consistent_stoichiometry"\
        % options.main_logger_name)
    cmpds = list(graph.compounds)
    # assign a starting point compound and discovered adjacent reactions
    source = random.choice(cmpds)
    discovered = set()
    complete = set()
    discover_partners(graph, source, discovered, complete)
    # almost empty mass vector that will be filled in a bfs type traversal of the
    # network
    mass_vector = dict(zip(cmpds, numpy.zeros(len(cmpds))))
    mass_vector[source] = 1.
    total = float(len(graph.reactions))
    prog = total
    maxi = max(factors)
    mini = 1. / float(len(graph.compounds))
    while len(discovered) > 0:
        balance_mass = True
        rxn = discovered.pop()
        for comp in graph.predecessors(rxn):
            discover_partners(graph, comp, discovered, complete)
            if mass_vector[comp] == 0.:
                balance_mass = False
        for comp in graph.successors(rxn):
            discover_partners(graph, comp, discovered, complete)
            if mass_vector[comp] == 0.:
                balance_mass = False
        if balance_mass:
            try:
                balance_reaction_by_mass(graph, factors, rxn, mass_vector, maxi)
            except PyMetabolismError:
                balance_reaction_by_mass(graph, factors, rxn, mass_vector, numpy.inf)
        else:
            try:
                balance_reaction_by_factors(graph, factors, rxn, mass_vector, numpy.inf, mini)
            except PyMetabolismError:
                # relax coefficients
                raise NotImplementedError
        complete.add(rxn)
        prog -= 1
        logger.info("%.2f %% complete.", (total - prog) / total * 100.)

def propagate_consistent_stoichiometry(graph, factors):
    """
    algorithm that chooses an initial compound in a bipartite metabolic network,
    then from the accompanying initial reactions chooses balancing masses that
    will spread along the graph in a wave-like fashion.
    """
    options = OptionsManager()
    logger = logging.getLogger("%s.propagate_consistent_stoichiometry"\
        % options.main_logger_name)
    cmpds = list(graph.compounds)
    # assign a starting point compound and discovered adjacent reactions
    source = random.choice(cmpds)
    discovered_rxns = set()
    complete_rxns = set()
    discovered_cmpds = set()
    complete_cmpds = set()
    discovered_cmpds.add(source)
    discover_partners(graph, source, discovered_rxns, complete_rxns)
    mass_vector = dict(zip(cmpds, numpy.zeros(len(cmpds))))
    total = float(len(graph.reactions))
    prog = total
    maxi = max(factors)
    mini = 1. / float(len(cmpds))
    while len(discovered_rxns) > 0:
#        try:
        balance_multiple_reactions(graph, factors, discovered_rxns, mass_vector, mini, numpy.inf, logger)
#        except PyMetabolismError:
#            break
        # move compound(s) to completed
        while len(discovered_cmpds) > 0:
            comp = discovered_cmpds.pop()
            complete_cmpds.add(comp)
        # move reaction(s) to completed
        while len(discovered_rxns) > 0:
            rxn = discovered_rxns.pop()
            discover_partners(graph, rxn, discovered_cmpds, complete_cmpds)
            complete_rxns.add(rxn)
        # add reactions from the newly found compounds
        for comp in discovered_cmpds:
            discover_partners(graph, comp, discovered_rxns, complete_rxns)
        logger.info(mass_vector.values())

def cycles_consistent_stoichiometry(graph, factors):
    """
    Attempt to create a consistent stoichiometry by stabilising cycles within the
    network.
    """
    options = OptionsManager()
    logger = logging.getLogger("%s.cycles_consistent_stoichiometry"\
        % options.main_logger_name)
    cmpds = list(graph.compounds)
    un_graph = graph.to_undirected()
    # attempt to get a consistent
    cycles = cycle_basis(un_graph)
    cycles.sort(cmp=lambda x, y: cmp(len(x), len(y)), reverse=True)

def random_fba(graph, num_inputs, num_outputs):
    """
    @param graph: BipartiteMetabolicNetwork

    @param num_inputs: Sequence of integers for a number of compound sources on
        the same set of sink compounds

    @param num_outputs: Size of the sink vector
    """
    options = OptionsManager()
    logger = logging.getLogger("%s.random_fba" % options.main_logger_name)
    assert num_outputs > 0
    if not num_inputs:
        raise PyMetabolismError("No sources specified for random FBA!")
    for input in num_inputs:
        assert input > 0
    inputs = list()
    potential_inputs = list()
    outputs = list()
    results = list()
    cmpds = list(graph.compounds)
    for comp in graph.compounds:
        if len(outputs) < num_outputs and graph.out_degree(comp) == 0:
            outputs.append(comp)
            cmpds.remove(comp)
        elif graph.in_degree(comp) == 0:
            potential_inputs.append(comp)
            cmpds.remove(comp)

    # if necessary fill up 'outputs' with random compounds
    # they should not be source compounds, thus we remove references
    count = num_outputs - len(outputs)
    while (count > 0):
        comp = random.choice(cmpds)
        outputs.append(comp)
        cmpds.remove(comp)
        count -= 1
    logger.info("%d biomass compounds.", len(outputs))

    # perform FBA for varying numbers of source compounds
    matrix = StoichiometricMatrix()
    matrix.make_new_from_network(graph)
    
    # add biomass reaction
    biomass = Reaction("Biomass", tuple(outputs), (), [-1.]*num_outputs)
    
    # equality constraints
    b_eq = numpy.zeros(matrix.num_compounds)

    # start with the largest number of inputs, then reduce
    num_inputs.sort(reverse=True)
    maxi = num_inputs[0]
    for comp in potential_inputs:
            if len(inputs) < maxi:
                inputs.append(comp)
    # similarly to 'outputs', 'inputs' is filled up
    count = maxi - len(inputs)
    while (count > 0):
        comp = random.choice(cmpds)
        inputs.append(comp)
        cmpds.remove(comp)
        count -= 1
    
    for num in num_inputs:
        system = numpy.array(matrix, dtype=float, copy=True)
        system.add_reaction(biomass)
        logger.info("%d uptake compounds.", num)
        input_rxns = list()
        # adjust source numbers
        while len(inputs) > num:
            comp = random.choice(inputs)
            inputs.remove(comp)

        transp = list()
        for comp in inputs:
            rxn = Reaction("Transport_%s" % comp.identifier,\
                (), (comp), (1.), reversible=True)
            system.add_reaction(rxn)
            transp.append(rxn)

        # objective function
        objective = numpy.zeros(system.num_reactions)
        objective[system.reaction_map[biomass]] = 1.

        # equality constraints
        A_eq = system

        # prepare flux balance analysis linear programming problem
        
        # boundaries
        lb = numpy.zeros(system.num_reactions)
        for rxn in system.reactions:
            if rxn.reversible:
                lb[system.reaction_map[rxn]] = -numpy.inf
        ub = numpy.empty(system.num_reactions)
        ub.fill(numpy.inf)
        # constrain transport reactions
        for rxn in transp:
            ub[system.reaction_map[rxn]] = 1.
            lb[system.reaction_map[rxn]] = -1.
                    
        problem = openopt.LP(f=objective, Aeq=A_eq, beq=b_eq, lb=lb, ub=ub)
        result = problem.maximize(options.solver, iprint=-10)
        logger.debug(result.xf)
        for rxn in graph.reactions:
            if rxn == biomass:
                logger.info("%s : %G", rxn.identifier, result.xf[system.reaction_map[rxn]])
            elif rxn in transp:
                logger.info("%s : %G", rxn.identifier, result.xf[system.reaction_map[rxn]])
            else:
                logger.info("%s : %G", rxn.identifier, result.xf[system.reaction_map[rxn]])
        results.append((result.isFeasible, result.xf))
    return results
