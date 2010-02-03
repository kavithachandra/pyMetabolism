#!/usr/bin/env python
# encoding: utf-8
"""
glpk.py

Created by Nikolaus Sonnenschein on 2010-02-02.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

from ifba.glpki.glpki import *
from ifba.GlpkWrap.metabolism import Metabolism, StoichiometricMatrix
from pyMetabolism.io.sbml import SBMLParser

class Metabolism2glpk(object):
    """docstring for metabolism2glpk"""
    def __init__(self, metabolism):
        super(Metabolism2glpk, self).__init__()
        self.metabolism = metabolism
        
    def convert_to_ifba_metabolism(self):
        """docstring for convert_to_ifba_metabolism"""
        lp = glp_create_prob()
        
    def convert_to_ifba_glpk(self):
        """docstring for convert_to_ifba_metabolism"""
        pass
    

if __name__ == '__main__':
    smallModel = '../test_data/Ec_core_flux1.xml'
    bigModel = '../test_data/iAF1260.xml'
    parser = SBMLParser(smallModel)
    metbol = parser.get_metabolic_system()
    converter = Metabolism2glpk(metbol)
    converter.convert_to_ifba_metabolism()
    

