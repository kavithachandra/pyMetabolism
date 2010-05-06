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
    from pyMetabolism.metabolism import Metabolism, Reaction

    import re
    smallModel = '../test_data/Ec_core_flux1.xml'
    bigModel = '../test_data/iAF1260.xml'
    parser = SBMLParser(smallModel, rprefix='R_', rsuffix='', cprefix='M_', csuffix=re.compile('_.$'))
    metbol = parser.get_metabolic_system()
    s = StoichiometricMatrix()
    s.make_new_from_system(metbol)
    # for elem in s.reactions:
    #     print elem.identifier
    print s.reaction_map[Reaction('Biomass_Ecoli_core_N__w_GAM_', (), (), ())]
    # print Reaction('Biomass_Ecoli_core_N__w_GAM_')
    
    obj = numpy.zeros(s.num_reactions)
    obj[s.reaction_map[Reaction('Biomass_Ecoli_core_N__w_GAM_', (), (), ())]] = 1.
    lp = 