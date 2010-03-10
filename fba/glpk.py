#!/usr/bin/env python
# encoding: utf-8
"""
glpk.py

Created by Nikolaus Sonnenschein on 2010-02-02.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

from ifba.glpki.glpki import *
import ifba.GlpkWrap.metabolism as ifba
from ifba.GlpkWrap.util import WriteCplex
from pyMetabolism.io.sbml import SBMLParser
from pyMetabolism.metabolism import Metabolism
from pyMetabolism.stoichiometry import StoichiometricMatrix

class Metabolism2glpk(object):
    """docstring for metabolism2glpk"""
    def __init__(self, metabolism):
        super(Metabolism2glpk, self).__init__()
        self.metabolism = metabolism
        self.stoich = StoichiometricMatrix()
        self.stoich.add_stoichiometry_from(self.metabolism)
    
    def __stoich_to_sparse_triplet(self):
        """docstring for fname"""
        non_zero = self.stoich.matrix.nonzero()
        num_non_zero = len(non_zero[0])
        non_zero_elements = self.stoich.matrix[non_zero]
        ia = intArray(num_non_zero)
        ja = intArray(num_non_zero)
        ar = doubleArray(num_non_zero)
        print non_zero
        for i in range(0, num_non_zero):
            ia[i+1] = non_zero[0][i] + 1
            ja[i+1] = non_zero[1][i] + 1
            ar[i+1] = non_zero_elements[i]
        print num_non_zero
        print ia[1]
        return (num_non_zero, ia, ja, ar)
    
    def convert_to_ifba_metabolism(self):
        """docstring for convert_to_ifba_metabolism"""
        lp = glp_create_prob()
        (num_non_zero, ia, ja, ar) = self.__stoich_to_sparse_triplet()
        glp_add_rows(lp, self.stoich.get_num_rows())
        glp_add_cols(lp, self.stoich.get_num_cols())
        glp_load_matrix(lp, num_non_zero, ia, ja, ar)
        for i in range(self.stoich.get_num_rows()):
            rev_compound_map = dict((v,k) for k, v in self.stoich.compound_map.iteritems())
            glp_set_row_name(lp, i+1, rev_compound_map[i].get_id())
        for i in range(self.stoich.get_num_cols()):
            rev_reaction_map = dict((v,k) for k, v in self.stoich.reaction_map.iteritems())
            glp_set_col_name(lp, i+1, rev_reaction_map[i].get_id())
        fba_object = ifba.Metabolism(lp)
        col_bounds = dict()
        for r in self.metabolism:
            if r.reversible:
                print r
                col_bounds[r.get_id()] = (-1000, 1000)
            else:
                col_bounds[r.get_id()] = (0, 1000)
        fba_object.modifyColumnBounds(col_bounds)
        fba_object.modifyRowBounds(dict([(r, (0., 0.)) for r in fba_object.getRowIDs()]))
        return fba_object
        
    def convert_to_ifba_glpk(self):
        """docstring for convert_to_ifba_metabolism"""
        pass
    

if __name__ == '__main__':
    smallModel = '../test_data/Ec_core_flux1.xml'
    bigModel = '../test_data/iAF1260.xml'
    parser = SBMLParser(smallModel)
    metbol = parser.get_metabolic_system()
    print metbol[0]
    s = StoichiometricMatrix()
    s.add_stoichiometry_from(metbol)
    print s.matrix
    converter = Metabolism2glpk(metbol)
    lp = converter.convert_to_ifba_metabolism()
    print lp
    print lp.getColumnIDs()
    print lp.getObjectiveFunction()
    lp.simplex()
    print lp
    extra_mets = [m.get_id() for m in metbol.get_compounds() if m.compartment == 'Extra_organism']
    lp.freeMetabolites(extra_mets)
    lp.setOptFlag('Max')
    lp.setReactionObjective('R_Biomass_Ecoli_core_N__w_GAM_')
    print lp.getObjectiveFunction()
    print lp.fba()
    print lp
    WriteCplex(lp, 'debug.lp')
    # print lp
    # print lp.getColumnBounds()
    # print lp.getRowBounds()

