# !/usr/bin/env python
# encoding: utf-8
"""
sbml.py

Created by Nikolaus Sonnenschein on 2010-01-11.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from pyMetabolism.metabolism import Metabolism, Reaction, Compound
try:
    import libsbml
except ImportError, msg:
    print "Please install libsbml to use SBML functionality!"
    print msg

def sbml2metabolism(path):
    """Loads a SBML model and returns a Metabolism object."""
    model = libsbml.readSBML('../test_data/iAF1260.xml').getModel()
    reactionsList = model.getListOfReactions()
    return reactionsList



if __name__ == '__main__':
    tmp = sbml2metabolism('../test_data/iAF1260.xml')
    print tmp[0]
    
