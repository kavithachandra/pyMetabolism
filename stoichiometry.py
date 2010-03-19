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
from pyMetabolism.metabolism_logging import NullHandler


class StoichiometricMatrix(object):
    """A class representing a stoichiometric matrix.

    Columns represent reactions
    Rows represent compounds
    Coefficients ...
    """

    _counter = 0

    def __init__(self, * args, ** kwargs):
        self.__class__._counter += 1
        super(StoichiometricMatrix, self).__init__(*args, ** kwargs)
        self.logger = logging.getLogger("pyMetabolism.StoichiometricMatrix.%d"\
                                        % self.__class__._counter)
        self.handler = NullHandler
        self.logger.addHandler(self.handler)
        self.matrix = None
        self.compound_map = dict()
        self.reaction_map = dict()

    def __str__(self):
        """docstring for __str__"""
        return self.matrix.__str__()

    def get_num_rows(self):
        """
        @return: Returns the number of rows
        @rtype: C{int}
        """
        return shape(self.matrix)[0]

    def get_num_cols(self):
        """
        @return: Returns the number of rows
        @rtype: C{int}
        """
        return shape(self.matrix)[1]

    def add_stoichiometry_from(self, metabolism):
        """"""
        j = len(self.reaction_map)
        i = len(self.compound_map)
        for reaction in metabolism.reactions:
            if reaction not in self.reaction_map:
                if self.matrix == None:
                    i = self._init_matrix(reaction)
                    j = 1
                    continue
                self.reaction_map[reaction] = j
                j += 1
                self.matrix = hstack((self.matrix, zeros((self.matrix.shape[0],
                                     1))))
                for compound in reaction:
                    if compound not in self.compound_map:
                        self.compound_map[compound] = i
                        i += 1
                        self.matrix = vstack((self.matrix, zeros((1,
                                             self.matrix.shape[1]))))
                    self.matrix[self.compound_map[compound]]\
                        [self.reaction_map[reaction]] = \
                        reaction.get_stoich_coeff(compound)

    def _init_matrix(self, reaction):
        num = len(reaction)
        self.matrix = zeros((num, 1))
        for i, compound in enumerate(reaction):
            self.compound_map[compound] = i
            self.matrix[i][0] = reaction.get_stoich_coeff(compound)
        self.reaction_map[reaction] = 0
        return num

if __name__ == '__main__':
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
