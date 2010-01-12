"""
Multiple ways of generating a null model for the statistical analysis of a
L{Metabolism} instance.

@todo: set up logging, complete docs, complete functionality
"""


import os
import sys

import collections as coll
import numpy as np
import networkx as nx
import decimal as dec
import random as rnd
import FuncDesigner as fd
import openopt as oo


def _double_edges(graph):
    """
'_double_edges' returns a list with all bidirectional edges present in the graph.

Parameters:
'graph' = ###networkx.DiGraph###

Returns:
###list### of ###edge###

    """
    double = list()
    for edge in graph.edges_iter():
        if graph.has_edge(edge[1], edge[0]):
            double.append(edge)
    return double

def _check_conditions(graph, first, second):
    """
'_check_conditions' checks if all necessary conditions to flip two edges are
fulfilled.

Parameters:
'graph' = ###networkx.DiGraph###
'first' = ###edge###
'second' = ###edge###

Returns:
###bool###

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

def _switch_double_edge(graph, first):
    """
'_switch_double_edge' attempts to switch a double edge with 'first'.

Parameters:
'graph' = ###networkx.DiGraph###
'first' = ###edge###

Returns:
###bool###

    """
    # find another bi-directional edge
    double_edges = _double_edges(graph)
    # not enough bi-directional edges in the graph
    if len(double_edges) < 2:
        return False
    second = rnd.choice(double_edges)
    if first == second:
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

def _select_edges(graph, metbs, rxns):
    """
'flip_bipartite_edge' performs the actual switching of edges according to a set
of rules.

Parameters:
'graph' = ###networkx.DiGraph###
'metbs' = ###iterable### of ###str###
'rxns' = ###iterable### of ###str###

Returns:
###bool###

    """
    # source is a metabolite
    edges = graph.edges()
    first = rnd.choice(edges)
    second = rnd.choice(edges)
    if first[0] in metbs:
        while not second[0] in metbs:
            second = rnd.choice(edges)
    else:
        while not second[0] in rxns:
            second = rnd.choice(edges)
    return (first, second)

def _flip_bipartite_edge(graph, metbs, rxns):
    """
'flip_bipartite_edge' performs the actual switching of edges according to a set
of rules.

Parameters:
'graph' = ###networkx.DiGraph###
'metbs' = ###iterable### of ###str###
'rxns' = ###iterable### of ###str###

Returns:
###bool###

    """
    (first, second) = _select_edges(graph, metbs, rxns)
    # same edges were picked by random algorithm
    if first == second:
        return False
    # first edge is bi-directional
    if graph.has_edge(first[1], first[0]):
        return _switch_double_edge(graph, first)
    # check other necessary conditions
    if _check_conditions(graph, first, second):
        graph.add_edge(first[0], second[1])
        graph.add_edge(second[0], first[1])
        graph.remove_edge(first[0], first[1])
        graph.remove_edge(second[0], second[1])
        return True
    else:
        return False

def randomise_bipartite(graph, metbs, rxns, flip=100):
    """
'randomise_bipartite' produces a randomised version of the supplied bipartite
graph 'graph'. This is achieved by switching edges in the graph a number of
times equal to the number of edges present times 'flip'. This function is
intended as a basis for calculating three-node subgraph statistics in metabolic
networks. As such only degrees of nodes, bi-directional link properties, etc.
(see:
R Milo, S Shen-Orr, S Itzkovitz, N Kashtan, D Chklovskii & U Alon,
Network Motifs: Simple Building Blocks of Complex Networks
Science, 298:824-827 (2002)
) are conserved. For larger subgraph statistics also the smaller subgraph
statistics would have to be conserved. The novelty of this function is to
preserve the bipartite nature of the graph. The randomised version is returned.

Parameters:
'graph' = ###networkx.DiGraph###
'metbs' = ###iterable### of ###str###
'rxns' = ###iterable### of ###str###
'flip' = ###int###

Returns:
###networkx.DiGraph###

    """
    rnd_graph = graph.copy()
    num_rnd = graph.size() * flip
    print "Number of bi-directional edges in the graph: %d"\
        % len(_double_edges(graph))
    while num_rnd > 0:
        # keep attempting to switch edges until right conditions are met
        if _flip_bipartite_edge(rnd_graph, metbs, rxns):
            num_rnd -= 1
    print "Number of bi-directional edges in the graph after randomisation: %d"\
        % len(_double_edges(graph))
    return rnd_graph


def _bfs(graph, metbs, rxns, mols, disc, status, stoich):
    while True:
        elem = disc.popleft()
        if status[elem]:
            continue
        for node in graph.successors_iter(elem):
            disc.append(node)
            # if node belongs to a reaction system set molarities and weights
        status[elem] = np.bool_(True)


def random_metabolism_bfs(graph, metbs, rxns, mols):
    """
'random_metabolism_bfs' should build a consistent stoichiometric matrix from a
directed bipartite graph object. It requires a list of all nodes considered
metabolites, of all nodes considered reactions, and a list of molarities that
will be randomly chosen to construct proper reactions. Note that molarities can
be biased by inserting certain numbers multiple times. Basically, this function
performs a breadth-first search on the graph assigning molarities from 'mols' to
edge weights whenever possible and otherwise reconciles metabolite masses.
Parameters:
'graph' = ###networkx.DiGraph###
'metbs' = ###list### of ###str###
'rxns' = ###list### of ###str###
'mols' = ###list### of ###int###

    """
    # a map that stores exact masses of metabolites, thus decimal type
    mass = coll.defaultdict(dec.Decimal) # zero initialised decimal values
    # stoichiometric matrix
    stoich = np.zeros((len(metbs),len(rxns), np.int, 'F'))
    # breadth-first search containers
    disc = coll.deque() # deque used as a queue for discovered nodes
    # status of nodes: 'False' = unknown, 'True' = complete
    status = coll.defaultdict(np.bool_)
    # initialise search
    elem = rnd.choice(graph.nodes())
    mass[elem] = dec.Decimal(1)
    disc.append(elem)
    try:
        _bfs(graph, metbs, rxns, mols, disc, status, stoich)
    except IndexError: # discovered queue is empty, i.e. search complete
        pass
    # do some consistency checks for example all masses > 0 and
    # each reaction must satisfy mass substrates = mass products


def scale_free_with_constraints(pars, rxns):
    """
This function should produce a scale-free bipartite graph implementing
metabolites and reactions in a fashion that produces a valid stoichiometric
matrix. This is an adaptation of the algorithm for directed scale-free graphs
(Bollobas et al, ACM-SCIAM, 2003).

    """
    if (sum(pars) != 1.):
        raise nx.NetworkXException
    # inserting a reaction with chosen metabolite pointing to it with a certain
    # weight and possibly recording substrates/products
    
    # inserting a metabolite with a chosen reaction pointing to it with
    # appropriate weight (needs exact definition)
    # possibly reducing weight of existing products
    pass


def random_metabolism_lp(graph, metbs, rxns, mols, seed = None):
    """
'random_metabolism_lp' should build a consistent stoichiometric matrix from a
directed bipartite graph object. It requires a list of all nodes considered
metabolites, of all nodes considered reactions, and a list of molarities that
will be randomly chosen to construct proper reactions. Note that molarities can
be biased by inserting certain numbers multiple times. This function solves the
following linear programming problem: Minimize the sum over all stoichiometric
coefficients, subject to transpose T(Stoichiometric Matrix S) dot Mass
Conservation Vector m = zero 0, where stoichiometric coefficients are part of a
population of natural numbers.
Parameters:
'graph' = ###networkx.DiGraph###
'metbs' = ###list### of ###str###
'rxns' = ###list### of ###str###
'mols' = ###list### of ###int###
'seed' = ###hashable###

    """
    if seed:
        rnd.seed(seed)
    # a map that stores exact masses of metabolites, thus decimal type
    mass = dict()
    # map with metabolite matrix-indeces
    metb_idx = dict()
    for (i, metb) in enumerate(metbs):
        metb_idx[metb] = i
        mass[metb] = rnd.random()
    # map with reaction rates
    rates = dict()
    # map with reaction matrix-indeces
    rxn_idx = dict()
    for (i, rxn) in enumerate(rxns):
        rxn_idx[rxn] = i
        # rates[rxn] = rnd.random() # not required atm
    m = len(metbs)
    n = len(rxns)
    # problem: minimise sum over stoichiometric coefficients
    f = np.ones([n, m], dtype=float)
    # subject to T(S).m = 0
    a_eq = np.empty([n, m], dtype=float)
    for rxn in rxns:
        for metb in metbs:
            a_eq[rxn_idx[rxn], metb_idx[metb]] = mass[metb]
    b_eq = np.zeros(n, dtype=float)
    # where
    lb = np.zeros([n, m], dtype=float)
    for rxn in rxns:
        for metb in graph.predecessors_iter(rxn):
            lb[rxn_idx[rxn], metb_idx[metb]] = -100.
        for metb in graph.successors_iter(rxn):
            lb[rxn_idx[rxn], metb_idx[metb]] = 1.
    ub = np.zeros([n, m], dtype=float)
    for rxn in rxns:
        for metb in graph.predecessors_iter(rxn):
            lb[rxn_idx[rxn], metb_idx[metb]] = -1.
        for metb in graph.successors_iter(rxn):
            lb[rxn_idx[rxn], metb_idx[metb]] = 100.
    # solve
    p = oo.LP(f=f, A=None, Aeq=a_eq, b=None, beq=b_eq, lb=lb, ub=ub)
    #p.debug = 1
    result = p.solve('cvxopt_lp')
    print result.ff
    print result.xf


def main(argv):
    """
First attempt at implementing random metabolism.

    """
    metbs = ['a', 'b', 'c', 'd']
    rxns = ['R', 'Rev', 'R2', 'R3', 'R4']
    bipartite = nx.DiGraph()
    bipartite.add_nodes_from(metbs)
    bipartite.add_nodes_from(rxns)
    bipartite.add_edge('a', 'R')
    bipartite.add_edge('b', 'R')
    bipartite.add_edge('R', 'c')
    bipartite.add_edge('R', 'd')
    bipartite.add_edge('c', 'Rev')
    bipartite.add_edge('d', 'Rev')
    bipartite.add_edge('Rev', 'a')
    bipartite.add_edge('Rev', 'b')
    bipartite.add_edge('R2', 'b')
    bipartite.add_edge('a', 'R2')
    bipartite.add_edge('R3', 'c')
    bipartite.add_edge('a', 'R3')
    bipartite.add_edge('R4', 'd')
    bipartite.add_edge('a', 'R4')
    print bipartite.edges()
    print randomise_bipartite(bipartite, metbs, rxns).edges()
    


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Usage:\n%s <number of metabolites : int> <number of reactions>"\
            % sys.argv[0]
        sys.exit(2)
    else:
        sys.exit(main(sys.argv[1:]))
