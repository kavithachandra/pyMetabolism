# !/usr/bin/env python
# encoding: utf-8
"""
sbml.py

Created by Nikolaus Sonnenschein on 2010-01-11.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
#from pyMetabolism.metabolism import Metabolism, Reaction, Compound
import libsbml
try:
    import libsbml
except ImportError, msg:
    print "Please install libsbml to use SBML functionality!"
    print msg

def sbml2metabolism(path):
    """Loads a SBML model and returns a Metabolism object."""
    tmp = libsbml.readSBML('../test_data/BIOMD0000000001.xml')
    model = tmp.getModel()
    print model.getId()
    print model.getListOfSpecies()[0].getId()

def sbml2metabolismSegFault(path):
    """Loads a SBML model and returns a Metabolism object."""
    model = libsbml.readSBML('../test_data/BIOMD0000000001.xml').getModel()
    print model.getId()
    print model.getListOfSpecies()[0].getId()

if __name__ == '__main__':
    print "Executing sbml2metabolism"
    sbml2metabolism('../test_data/iAF1260.xml')
    print "Executing sbml2metabolismSegFault"
    sbml2metabolismSegFault('../test_data/iAF1260.xml')
