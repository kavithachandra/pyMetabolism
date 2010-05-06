#! /usr/bin/env python
# encoding: utf-8


from pyMetabolism import OptionsManager
import numpy
import openopt


def verify_consistency(stoich_matrix):
    """
    Gevorgyan, et al(2008)
    """
    options = OptionsManager()
    # objective function
    objective = numpy.ones(stoich_matrix.num_compounds)
    # inequal constraints
    # equal constraints
    try:
        A_eq = numpy.array(stoich_matrix.matrix, dtype=float, copy=False)
    except ValueError:
        options.logger.critical("Could not set stoichiometric matrix!")
    A_eq = A_eq.transpose()
    b_eq = numpy.zeros(stoich_matrix.num_reactions)
    # boundaries
    lb = numpy.ones(stoich_matrix.num_compounds)
    ub = numpy.empty(stoich_matrix.num_compounds)
    ub.fill(numpy.inf)
    problem = openopt.LP(f=objective, Aeq=A_eq, beq=b_eq, lb=lb, ub=ub)
    result = problem.solve(options.solver, iprint=-2)
    options.logger.debug(str(result.xf))
    return (result.isFeasible, result.xf)

def detect_unconserved_metabolites(stoich_matrix, masses):
    """
    Gevorgyan, et al(2008)
    """
    raise NotImplementedError
    options = OptionsManager()
    # objective function
    objective = numpy.ones(stoich_matrix.num_compounds)
    # inequal constraints
    # equal constraints
    A_eq = stoich_matrix.matrix.copy()
    A_eq = A_eq.transpose()
    b_eq = numpy.zeros(stoich_matrix.num_reactions)
    # boundaries
    lb = numpy.zeros(stoich_matrix.num_compounds)
    ub = numpy.ones(stoich_matrix.num_compounds)
#    ub = numpy.array(masses)
    # binary variables
    bool_vals = range(stoich_matrix.num_compounds)
    problem = openopt.MILP(f=objective, Aeq=A_eq, beq=b_eq, lb=lb, ub=ub,\
        intVars=bool_vals)
    result = problem.solve(options.solver, iprint=-2, goal="max")
    options.logger.debug(result.xf)
    return (result.isFeasible, result.xf)
