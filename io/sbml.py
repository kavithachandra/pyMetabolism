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
    
    @ivar path: path to SBML model
    @type path: C{str}

    @ivar rprefix: prefix e.g. 'R_' that should be stripped from reaction identifiers
    @type rprefix: C{str} or compiled regex object e.g. re.compile('R.')

    @ivar rsuffix: suffix e.g. '_Ec_Coli_Core' that should be stripped from reaction identifiers
    @type rsuffix: C{str} or compiled regex object e.g. re.compile('_.*$')

    @ivar cprefix: prefix e.g. 'M_' that should be stripped from compound identifiers
    @type cprefix: C{str} or compiled regex object e.g. re.compile('M.')

    @ivar csuffix: suffix e.g. '_c' that should be stripped from compound identifiers
    @type csuffix: C{str} or compiled regex object e.g. re.compile('_.$')

    @todo: implement convenience stuff
    @todo: Fix suffix and prefix handling
    @todo: Use logging
    @todo: Check for boundary conditions
    """
    def __init__(self, path, rprefix='', rsuffix='', cprefix='', csuffix=''):
        super(SBMLParser, self).__init__()
        self.sbml_document = libsbml.readSBML(path)
        if self.sbml_document.getNumErrors() > 0:
            raise Exception, "libsbml function getNumErrors return value greater 0!"
        self.model = self.sbml_document.getModel()
        self.rprefix = rprefix
        self.rsuffix = rsuffix
        self.cprefix = cprefix
        self.csuffix = csuffix
        
    def _parse_sbml_reactant(self, sbml_species_reference):
        """Able to parse entries from getListOfReactants or getListOfProducts
        
        @todo: Check for meta information and parse if available
        """
        sbml_species = self.model.getSpecies(sbml_species_reference.getSpecies())
        return self._parse_sbml_species(sbml_species)

    def _parse_sbml_species(self, sbml_compound):
        """Able to parse entries from getListOfSpecies
        
        @todo: Check for meta information and parse if available
        """
        comp_id = re.sub(self.csuffix, '', re.sub(self.cprefix, '', sbml_compound.getId()))
        return CompartCompound(Compound(comp_id), Compartment(sbml_compound.getCompartment(), sbml_compound.getConstant()))
    
    def _parse_sbml_reaction(self, sbml_reaction):
        """Able to parse entries from getListOfReactions"""
        identifier = re.sub(self.rsuffix, '', re.sub(self.rprefix, '', sbml_reaction.getId()))
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
        return Reaction(identifier, substrates, products, stoichiometry, reversible=sbml_reaction.getReversible())
    
    def get_reactions(self):
        """Returns a list of reactions parsed from the SBML model"""
        list_of_reactions = self.model.getListOfReactions()
        return tuple([self._parse_sbml_reaction(r) for r in list_of_reactions])
        
    def get_compounds(self):
        """Returns a list of all compounds parsed from the SBML model"""
        return [self._parse_sbml_species(elem) for elem in self.model.getListOfSpecies()]
        
    def get_metabolic_system(self, name=None):
        """Returns an instance of Metabolism including all reactions from the parsed SBML model"""
        if name:
            return Metabolism(self.get_reactions(), name=name)
        else:
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
    system = parser.get_metabolic_system()
    tmp = list(system.compounds)
    print tmp[0]
    for elem in dir(tmp[0]):
        print elem,  getattr(tmp[0], elem)
    for key in Compartment._memory:
        print "compartment", key
        x = CompartCompound(Compound('for'), Compartment(key, Compartment._memory[key].constant))
        print x.compartment, id(x)
    # print 'Compounds:\n'
    # for elem in parser.get_compounds():
    #     print elem
    # print len(parser.get_compounds())
    # print 'Reactions:\n'
    # for elem in parser.get_reactions():
    #     print elem
    system = parser.get_metabolic_system()
    for r in system:
        print r.identifier, '->', r
    # tmp = system.get_compounds()
    # print tmp[0]
    # for elem in dir(tmp[0]):
    #     print elem,  getattr(tmp[0], elem)
    
    parser = SBMLParser(smallModel, rprefix='R_', rsuffix='', cprefix='M_', csuffix=re.compile('_.$'))
    system = parser.get_metabolic_system(name='ModelWithStrippedSuffices')
    for r in system:
        print r.identifier, '->', r
