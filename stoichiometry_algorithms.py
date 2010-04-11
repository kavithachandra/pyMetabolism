#! /usr/bin/env python
# encoding: utf-8


from pyMetabolism import OptionsManager
import numpy
import openopt
from pyMetabolism.stoichiometry import StoichiometricMatrix
import logging


def form_matrix(stoich_matrix):
    """
    transform any stoichiometric matrix into a numpy array
    """

    try:
        matrix = numpy.array(stoich_matrix.matrix, dtype=float)
    except ValueError:
        options.logger.debug("Could not set stoichiometric matrix!")
    return matrix


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
    for (index, val) in numpy.ndenumerate(ub):
        ub[index] = numpy.inf
    problem = openopt.LP(f=objective, Aeq=A_eq, beq=b_eq, lb=lb, ub=ub)
    result = problem.solve(options.solver, iprint=-2)
    options.logger.debug(str(result.xf))
    return (result.isFeasible, result.xf)

def detect_unconserved_metabolites(stoich_matrix, masses):
    """
    Gevorgyan, et al(2008)
    """
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



if __name__ == '__main__':
    options = OptionsManager()
    options.logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    options.logger.addHandler(handler)
    import pyMetabolism.metabolism as metb
    cmp = metb.Compartment("Cytosol", True)
#    metb.CompartCompound(metb.Compound("fructose"), cmp)
    rxn1 = metb.Reaction("MAP", (metb.CompartCompound(metb.Compound("atp"), cmp),),\
        (metb.CompartCompound(metb.Compound("adp"), cmp), metb.CompartCompound(metb.Compound("p"), cmp)),\
        (2, 1, 1), reversible=True)
    rxn2 = metb.Reaction("FPK", (metb.CompartCompound(metb.Compound("adp"), cmp), metb.CompartCompound(metb.Compound("p"), cmp)),\
        (metb.CompartCompound(metb.Compound("atp"), cmp),),\
        (1, 2, 1))
    rxn3 = metb.Reaction("EX", (metb.CompartCompound(metb.Compound("atp"), metb.Compartment("extern", False)),),\
        (metb.CompartCompound(metb.Compound("atp"), cmp),metb.CompartCompound(metb.Compound("fructose"), cmp)),\
        (1, 1, 2))
    system = metb.Metabolism([rxn1, rxn2])
    options.logger.info(system)
    matrix = StoichiometricMatrix()
    matrix.make_new_from_system(system)
    options.logger.info(matrix)
    for x, val in matrix._compound_map.iteritems():
        print x, ":", val
    print
    for x, val in matrix._reaction_map.iteritems():
        print x, ":", val
    print
    (yes, results) = verify_consistency(matrix)
    options.logger.info(yes)
    (yes, results) = detect_unconserved_metabolites(matrix, results)
    options.logger.info("unconserved metabolites:")
    for (index, name) in enumerate(matrix.compounds):
        if results[index] == 0.:
            options.logger.info(name)
    from pyMetabolism.io.sbml import SBMLParser
    import re
    smallModel = './test_data/Ec_core_flux1.xml'
    bigModel = './test_data/iAF1260.xml'
    # parser = SBMLParser(bigModel)
    parser = SBMLParser(smallModel, rprefix='R_', rsuffix='', cprefix='M_', csuffix=re.compile('_.$'))
    metbol = parser.get_metabolic_system(parser)
    print metbol[0].identifier
    print metbol[0]
#    print "M_actp_c(Cytosol) is in reacs: ", metb.CompartCompound(metb.Compound('M_actp_c'), metb.Compartment("Cytosol"))
    s = StoichiometricMatrix()
    s.make_new_from_system(metbol)
    print s
    (yes, results) = verify_consistency(s)
    print yes
    (yes, results) = detect_unconserved_metabolites(s, results)
    print yes
    for (index, name) in enumerate(s.compounds):
        if results[index] == 0.:
            print "unconserved:", name
        elif results[index] == 1.:
            print "conserved:", name
