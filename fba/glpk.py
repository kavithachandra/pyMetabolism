#!/usr/bin/env python
# encoding: utf-8
"""
glpk.py

Created by Nikolaus Sonnenschein on 2010-02-02.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

from ifba.glpki.glpki import *
from ifba.GlpkWrap.metabolism import Metabolism
from pyMetabolism.io.sbml import SBMLParser


if __name__ == '__main__':
    smallModel = '../test_data/Ec_core_flux1.xml'
    bigModel = '../test_data/iAF1260.xml'
    parser = SBMLParser(smallModel)
    print 'Compounds:\n'
    for elem in parser.get_compounds():
        print elem
    print
    print 'Reactions:\n'
    for elem in parser.get_reactions():
        print elem
    print parser.get_metabolic_system()

