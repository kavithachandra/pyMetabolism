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
