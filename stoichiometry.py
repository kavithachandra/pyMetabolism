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
import copy
from numpy import hstack, vstack, zeros
from pyMetabolism import OptionsManager, new_property


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
        self._matrix = zeros(shape=(0,0))
        self._compound_map = dict()
        self._reaction_map = dict()
#        self._compound_index = dict()
#        self._reaction_index = dict()
#        self._permutation = None

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
    def compounds():
        return {"fset": None, "fget": lambda self: self._compound_map.keys(),\
            "doc": "get method"}

    @new_property
    def reactions():
        return {"fset": None, "fget": lambda self: self._reaction_map.keys(),\
            "doc": "get method"}

    @new_property
    def num_compounds():
        """
        @return: Returns the number of rows
        @rtype: C{int}
        """
        return {"fset": None, "fget": lambda self: self._matrix.shape[0],\
            "doc": "get method"}

    @new_property
    def num_reactions():
        """
        @return: Returns the number of rows
        @rtype: C{int}
        """
        return {"fset": None, "fget": lambda self: self._matrix.shape[1],\
            "doc": "get method"}

    def reaction_vector(self, rxn):
        """
        """
        return self._matrix[:, self._reaction_map[rxn]]

    def compound_vector(self, comp):
        """
        """
        return self._matrix[self._compound_map[comp], :]

    def __str__(self):
        """docstring for __str__"""
        return self._matrix.__str__()

    def __copy__(self):
        raise NotImplementedError

    def make_new_from_system(self, metabolism):
        """
        """
        self._matrix = zeros((len(metabolism.compounds),\
            len(metabolism.reactions)))
        for (i, comp) in enumerate(metabolism.compounds):
            self._compound_map[comp] = i
        for (i, rxn) in enumerate(metabolism.reactions):
            self._reaction_map[rxn] = i
        for rxn in metabolism.reactions:
            for comp in rxn:
                self._matrix[self._compound_map[comp], self._reaction_map[rxn]]\
                    = rxn.stoich_coeff(comp)

    def make_new_from_network(self, graph):
        """
        """
        self._matrix = zeros(shape=(len(graph.compounds), len(graph.reactions)))
        for (i, comp) in enumerate(graph.compounds):
            self._compound_map[comp] = i
        for (i, rxn) in enumerate(graph.reactions):
            self._reaction_map[rxn] = i
        for rxn in graph.reactions:
            self._logger.debug("Reaction: %s", rxn.identifier)
            self._logger.debug("Substrates:")
            for comp in graph.predecessors(rxn):
                self._logger.debug("%s", comp.identifier)
                self._matrix[self._compound_map[comp], self._reaction_map[rxn]]\
                    = -graph[comp][rxn]["factor"]
            self._logger.debug("Products:")
            for comp in graph.successors(rxn):
                self._logger.debug("%s", comp.identifier)
                self._matrix[self._compound_map[comp], self._reaction_map[rxn]]\
                    = graph[rxn][comp]["factor"]

    def add_reaction(self, reaction):
        """
        @todo: doc
        """
        j = self.num_reactions
        i = self.num_compounds
        self._reaction_map[reaction] = j
        self._matrix = hstack((self._matrix, zeros((self.num_compounds,
                             1))))
        for compound in reaction:
            if compound not in self._compound_map:
                self._compound_map[compound] = i
                i += 1
                self._matrix = vstack((self._matrix, zeros((1,
                                     self.num_reactions))))
            self._matrix[self._compound_map[compound],\
                self._reaction_map[reaction]] = \
                reaction.stoich_coeff(compound)

    def remove_reaction(self, reaction):
        """
        @todo: doc
        """
        raise NotImplementedError
        # have to adjust crappy indices
        participants = list()
        column = self.reaction_vector(reaction)
        for (i, val) in enumerate(column):
            if val != 0.:
                participants.append(column[i])
