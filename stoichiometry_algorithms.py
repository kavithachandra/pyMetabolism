#! /usr/bin/env python
# encoding: utf-8


import numpy
import openopt
from pyMetabolism.stoichiometry import StoichiometricMatrix


def verify_consistency(stoich_matrix):
    """
    Gevorgyan, et al(2008)
    """
    # objective function
    objective = numpy.ones(stoich_matrix.num_metabolites)
    # inequal constraints
    # equal constraints
    A_eq = stoich_matrix.matrix.copy()
    A_eq = A_eq.transpose()
    b_eq = numpy.zeros(stoich_matrix.num_reactions)
    # boundaries
    lb = numpy.ones(stoich_matrix.num_metabolites)
    ub = numpy.empty(stoich_matrix.num_metabolites)
    for (index, val) in numpy.ndenumerate(ub):
        ub[index] = numpy.inf
    problem = openopt.LP(f=objective, Aeq=A_eq, beq=b_eq, lb=lb, ub=ub)
    result = problem.solve("glpk", iprint=-2)
    return (result.isFeasible, result.xf)

def detect_unconserved_metabolites(stoich_matrix, masses):
    """
    Gevorgyan, et al(2008)
    """
    # objective function
    objective = numpy.ones(stoich_matrix.num_metabolites)
    # inequal constraints
    # equal constraints
    A_eq = stoich_matrix.matrix.copy()
    A_eq = A_eq.transpose()
    b_eq = numpy.zeros(stoich_matrix.num_reactions)
    # boundaries
    lb = numpy.zeros(stoich_matrix.num_metabolites)
    ub = numpy.ones(stoich_matrix.num_metabolites)
    # binary variables
    bool_vals = range(stoich_matrix.num_metabolites)
    problem = openopt.MILP(f=objective, Aeq=A_eq, beq=b_eq, lb=lb, ub=ub,\
        intVars=bool_vals)
    result = problem.solve("glpk", iprint=-2, goal="max")
    return (result.isFeasible, result.xf)



if __name__ == '__main__':
    import pyMetabolism.metabolism as metb
    cmp = metb.Compartment("Cytosol", True)
    rxn1 = metb.Reaction("MAP", (metb.CompartCompound(metb.Compound("atp"), cmp),),\
        (metb.CompartCompound(metb.Compound("adp"), cmp), metb.CompartCompound(metb.Compound("p"), cmp)),\
        (1, 1, 1))
    rxn2 = metb.Reaction("FPK", (metb.CompartCompound(metb.Compound("fructose"), cmp), metb.CompartCompound(metb.Compound("atp"), cmp)),\
        (metb.CompartCompound(metb.Compound("g6p"), cmp), metb.CompartCompound(metb.Compound("adp"), cmp)),\
        (1, 1, 1, 1))
    system = metb.Metabolism([rxn1, rxn2])
    print system
    matrix = StoichiometricMatrix()
    matrix.add_stoichiometry_from(system)
    print matrix
    (yes, results) = verify_consistency(matrix)
    print yes
    print results
    (yes, results) = detect_unconserved_metabolites(matrix, results)
    print "unconserved metabolites:"
    for (name, i) in matrix.compound_map.iteritems():
        if results[i] == 0.:
            print name
    print results
    