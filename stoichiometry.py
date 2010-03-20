#!/usr/bin/env python
# encoding: utf-8


"""
A class that aims to be able to extract all information
one can from the stoichiometric matrix of a list of reactions.

@author: Nikolaus Sonnenschein
@author: Moritz Beber
@contact: niko(at)foo.org
@contact: moritz(at)foo.org
@copyright: Jacobs University Bremen. All rights reserved.
@since: 2009-09-11
@todo: plenty matrix operations

"""


import logging
from numpy import hstack
from numpy import shape
from numpy import vstack
from numpy import zeros
from pyMetabolism import OptionsManager
from pyMetabolism import new_property


class StoichiometricMatrix(object):
    """A class representing a stoichiometric matrix.

    Columns represent reactions
    Rows represent compounds
    Coefficients ...
    """

    _counter = 0

    def __init__(self, *args, **kwargs):
        self.__class__._counter += 1
        super(StoichiometricMatrix, self).__init__(*args, **kwargs)
        self._options = OptionsManager()
        self._logger = logging.getLogger("%s.%s.%d"\
            % (self._options.main_logger_name, self.__class__,\
            self.__class__._counter))
        self._matrix = None
        self._compound_map = dict()
        self._reaction_map = dict()
        self._signum = None

    @new_property
    def options():
        return {"fset": None, "doc": "get method"}

    @new_property
    def logger():
        pass

    @new_property
    def matrix():
        return {"fset": None, "doc": "get method"}

    @new_property
    def compound_map():
        return {"fset": None, "doc": "get method"}

    @new_property
    def reaction_map():
        return {"fset": None, "doc": "get method"}

    @new_property
    def signum():
        return {"fset": None}

    @new_property
    def num_metabolites():
        """
        @return: Returns the number of rows
        @rtype: C{int}
        """
        return {"fset": None, "fget": lambda self: self._matrix.shape[0], "doc": "get method"}

    @new_property
    def num_reactions():
        """
        @return: Returns the number of rows
        @rtype: C{int}
        """
        return {"fset": None, "fget": lambda self: self._matrix.shape[1], "doc": "get method"}

    def __str__(self):
        """docstring for __str__"""
        return self._matrix.__str__()

    def add_stoichiometry_from(self, metabolism):
        """
        """
        j = len(self._reaction_map)
        i = len(self._compound_map)
        for reaction in metabolism.reactions:
            if reaction not in self._reaction_map:
                if self._matrix == None:
                    i = self._init_matrix(reaction)
                    j = 1
                    continue
                self._reaction_map[reaction] = j
                j += 1
                self._matrix = hstack((self._matrix, zeros((self._matrix.shape[0],
                                     1))))
                for compound in reaction:
                    if compound not in self._compound_map:
                        self._compound_map[compound] = i
                        i += 1
                        self._matrix = vstack((self._matrix, zeros((1,
                                             self._matrix.shape[1]))))
                    self._matrix[self._compound_map[compound]]\
                        [self._reaction_map[reaction]] = \
                        reaction.stoich_coeff(compound)

    def _init_matrix(self, reaction):
        num = len(reaction)
        self._matrix = zeros((num, 1))
        for i, compound in enumerate(reaction):
            self._compound_map[compound] = i
            self._matrix[i][0] = reaction.stoich_coeff(compound)
        self._reaction_map[reaction] = 0
        return num

if __name__ == '__main__':
    import pyMetabolism.metabolism as metb
    cmp = metb.Compartment("Cytosol", True)
    rxn1 = metb.Reaction("MAP", (metb.CompartCompound(metb.Compound("atp"), cmp),),\
        (metb.CompartCompound(metb.Compound("adp"), cmp), metb.CompartCompound(metb.Compound("p"), cmp)),\
        (1, 1, 1))
    rxn2 = metb.Reaction("FPK", (metb.CompartCompound(metb.Compound("fructose"), cmp), metb.CompartCompound(metb.Compound("atp"), cmp)),\
        (metb.CompartCompound(metb.Compound("g6p"), cmp), metb.CompartCompound(metb.Compound("adp"), cmp)),\
        (1, 1, 1, 1))
    system = metb.Metabolism([rxn1, rxn2])
    print system
    matrix = StoichiometricMatrix()
    matrix.add_stoichiometry_from(system)
    print matrix
    from pyMetabolism.io.sbml import SBMLParser
    from pyMetabolism.metabolism import Compound, CompartCompound, Compartment
    from numpy import *
    import re
    # set_printoptions(threshold='nan')
    smallModel = './test_data/Ec_core_flux1.xml'
    bigModel = './test_data/iAF1260.xml'
    # parser = SBMLParser(bigModel)
    parser = SBMLParser(smallModel, rprefix='R_', rsuffix='', cprefix='M_', csuffix=re.compile('_.$'))
    metbol = parser.get_metabolic_system(parser)
    print metbol[0].identifier
    print metbol[0]
    print "M_actp_c(Cytosol) is in reacs: ", CompartCompound(Compound('M_actp_c'), Compartment("Cytosol"))
    s = StoichiometricMatrix()
    s.add_stoichiometry_from(metbol)
    # print s.matrix[0:7]
    # print shape(s)
    # print [(c, s.compound_map[c]) for c in metbol[0].get_compounds()]
