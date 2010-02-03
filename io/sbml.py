# !/usr/bin/env python -V
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


class SBMLParser(object):
    """A class implementing methods for parsing a SBML model"""
    def __init__(self, path):
        super(SBMLParser, self).__init__()
        self.sbml_document = libsbml.readSBML(path)
        if self.sbml_document.getNumErrors() > 0:
            raise Exception, "libsbml function getNumErrors return value greater 0!"
        self.model = self.sbml_document.getModel()
        
    def __parse_sbml_reactant(self, sbml_compound):
        """Able to parse entries from getListOfReactants or getListOfProducts
        
        @todo: Check for meta information and parse if available
        """
        return Compound(sbml_compound.getSpecies())

    def __parse_sbml_species(self, sbml_compound):
        """Able to parse entries from getListOfSpecies
        
        @todo: Check for meta information and parse if available
        """
        return Compound(sbml_compound.getId())
    
    def __parse_sbml_reaction(self, sbml_reaction):
        """Able to parse entries from getListOfReactions"""
        identifier = sbml_reaction.getId()
        list_of_reactants = sbml_reaction.getListOfReactants()
        list_of_products = sbml_reaction.getListOfProducts()
        substrates = [self.__parse_sbml_reactant(elem) for elem in list_of_reactants]
        products = [self.__parse_sbml_reactant(elem) for elem in list_of_products]
        stoichiometry = tuple([-elem.getStoichiometry() for elem in list_of_reactants] + [elem.getStoichiometry() for elem in list_of_products])
        return Reaction(identifier, substrates, products, stoichiometry, reversible=sbml_reaction.getReversible())
    
    def get_reactions(self):
        """Returns a list of reactions parsed from the SBML model"""
        list_of_reactions = self.model.getListOfReactions()
        return tuple([self.__parse_sbml_reaction(r) for r in list_of_reactions])
        
    def get_compounds(self):
        """Returns a list of all compounds parsed from the SBML model"""
        return [Compound(elem.getId()) for elem in self.model.getListOfSpecies()]
        
    def get_metabolic_system(self):
        """Returns an instance of Metabolism including all reactions from the parsed SBML model"""
        return Metabolism(self.get_reactions(), name=self.model.getName())


if __name__ == '__main__':
    smallModel = '../test_data/Ec_core_flux1.xml'
    bigModel = '../test_data/iAF1260.xml'
    parser = SBMLParser(bigModel)
    print 'Compounds:\n'
    for elem in parser.get_compounds():
        print elem
    print len(parser.get_compounds())
    print 'Reactions:\n'
    for elem in parser.get_reactions():
        print elem
    print len(parser.get_reactions())
    print parser.get_metabolic_system()
    
