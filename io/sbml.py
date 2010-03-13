# !/usr/bin/env python -V
# encoding: utf-8
"""
sbml.py

Created by Nikolaus Sonnenschein on 2010-01-11.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import re
from pyMetabolism.metabolism import Metabolism, Reaction, Compound, Compartment,\
    CompartCompound


try:
    import libsbml
except ImportError, msg:
    print "Please install libsbml to use SBML functionality!"
    print msg


class SBMLParser(object):
    """A class implementing methods for parsing a SBML model
    
    @todo: implement convenience stuff
    """
    def __init__(self, path):
        super(SBMLParser, self).__init__()
        self.sbml_document = libsbml.readSBML(path)
        if self.sbml_document.getNumErrors() > 0:
            raise Exception, "libsbml function getNumErrors return value greater 0!"
        self.model = self.sbml_document.getModel()
        
    def _parse_sbml_reactant(self, sbml_species_reference):
        """Able to parse entries from getListOfReactants or getListOfProducts
        
        @todo: Check for meta information and parse if available
        """
        sbml_species = self.model.getSpecies(sbml_species_reference.getSpecies())
        return self._parse_sbml_species(sbml_species)

    def _parse_sbml_species(self, sbml_compound):
        """Able to parse entries from getListOfSpecies
        
        @todo: Check for meta information and parse if available
        @todo: Fix suffix and prefix handling
        """
        return CompartCompound(sbml_compound.getId().split('_')[1], \
            Compartment(sbml_compound.getCompartment(), sbml_compound.getConstant()))
    
    def _parse_sbml_reaction(self, sbml_reaction):
        """Able to parse entries from getListOfReactions"""
        identifier = sbml_reaction.getId()
        list_of_reactants = sbml_reaction.getListOfReactants()
        list_of_products = sbml_reaction.getListOfProducts()
        compartments = list()
        substrates = list()
        for elem in list_of_reactants:
            species_tmp = self._parse_sbml_reactant(elem)
            substrates.append(species_tmp)
        products = list()
        for elem in list_of_products:
            species_tmp = self._parse_sbml_reactant(elem)
            products.append(species_tmp)
        stoichiometry = tuple([-elem.getStoichiometry() for elem in list_of_reactants] + [elem.getStoichiometry() for elem in list_of_products])
        # print identifier
        # print stoichiometry, [elem.identifier for elem in compartments]
        # print len(stoichiometry), len(compartments)
        return Reaction(identifier, substrates, products, stoichiometry, reversible=sbml_reaction.getReversible())
    
    def get_reactions(self):
        """Returns a list of reactions parsed from the SBML model"""
        list_of_reactions = self.model.getListOfReactions()
        return tuple([self._parse_sbml_reaction(r) for r in list_of_reactions])
        
    def get_compounds(self):
        """Returns a list of all compounds parsed from the SBML model"""
        return [self._parse_sbml_species(elem) for elem in self.model.getListOfSpecies()]
        
    def get_metabolic_system(self):
        """Returns an instance of Metabolism including all reactions from the parsed SBML model"""
        return Metabolism(self.get_reactions(), name=self.model.getName())


if __name__ == '__main__':
    smallModel = '../test_data/Ec_core_flux1.xml'
    bigModel = '../test_data/iAF1260.xml'
    parser = SBMLParser(smallModel)
    print 'Compounds:\n'
    for elem in parser.get_compounds():
        print elem
    print len(parser.get_compounds())
    print 'Reactions:\n'
    for elem in parser.get_reactions():
        print elem
#    Compound('for')
    system = parser.get_metabolic_system()
    tmp = system.get_compounds()
    print tmp[0]
    for elem in dir(tmp[0]):
        print elem,  getattr(tmp[0], elem)
    for key in Compartment._memory:
        print "compartment", key
        x = CompartCompound('for', Compartment(key, Compartment._memory[key].constant))
        print x.compartment, id(x)
