#!/usr/bin/env python
# encoding: utf-8
"""
openopt.py

Created by Nikolaus Sonnenschein on 2010-05-05.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

import numpy
from openopt import LP
from pyMetabolism import OptionsManager
from pyMetabolism.stoichiometry import StoichiometricMatrix

class fba(object):
    """docstring for fba"""
    def __init__(self, metabolism):
        super(fba, self).__init__()
        self.metabolism = metabolism
        

if __name__ == '__main__':
    from pyMetabolism.io.sbml import SBMLParser
    from pyMetabolism.metabolism import Metabolism, Reaction, Compartment

    import re
    smallModel = '../test_data/Ec_core_flux1.xml'
    bigModel = '../test_data/iAF1260.xml'
    parser = SBMLParser(smallModel, rprefix='R_', rsuffix='', cprefix='M_', csuffix=re.compile('_.$'))
    metbol = parser.get_metabolic_system()
    s = StoichiometricMatrix()
    s.make_new_from_system(metbol)
    bMets =[m for m in s.compounds if m.compartment == Compartment('Extra_organism', ())] 
    # for m in bMets:
    #     print m
    bIndices = list()
    for m in bMets:
        bIndices.append(s.compound_map[m])
    print bIndices
    obj = numpy.zeros(s.num_reactions)
    obj[s.reaction_map[Reaction('Biomass_Ecoli_core_N__w_GAM_', (), (), ())]] = 1.
    print len(s.matrix)
    # for m in bMets:
    #     s.add_reaction(Reaction('R_'+m.identifier+'_Transp', (m,), (), ()))
    # A_eq = numpy.delete(numpy.array(s.matrix, dtype=float), bIndices, 0)
    transpColumn = numpy.zeros((len(s.matrix), len(bIndices)))
    for i, elem in enumerate(bIndices):
        transpColumn[elem,i] = 1
    A_eq = numpy.append(numpy.array(s.matrix, dtype=float), transpColumn, 1)
    print A_eq
    print numpy.shape(A_eq)
    b_eq = numpy.zeros(len(A_eq))
    print numpy.shape(A_eq)[1]
    print len(b_eq)
    print len(numpy.random.randint(-10, -5, numpy.shape(A_eq)[1]))
    print len(numpy.random.randint(5, 10, numpy.shape(A_eq)[1]))
    lp = LP(f=obj, Aeq=A_eq, beq=b_eq, lb=numpy.random.randint(-6, -0, numpy.shape(A_eq)[1]), ub=numpy.random.randint(0, 6, numpy.shape(A_eq)[1]))
    lp.exportToMPS
    # print help(lp.manage)
    # lp.solve('cvxopt_lp')
    # print help(lp.solve)
    r = lp.solve('glpk', goal='max')
    # print dir(r)
    # print r.ff
    # print r.evals
    # print r.isFeasible
    # print r.msg
    # print r.xf
    # print r.duals
    # print r.rf
    # print r.solverInfo
    # problem = openopt.LP(f=objective, Aeq=A_eq, beq=b_eq, lb=lb, ub=ub)