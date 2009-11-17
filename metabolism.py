#!/usr/bin/env python
# encoding: utf-8
"""
pyMetabolism.py

Some notes:
- Reactions and Compounds are regarded immutable

Created by Nikolaus Sonnenschein on 2009-09-10.
Copyright (c) 2009 . All rights reserved.
"""

import sys
import os
import networkx
from numpy import array, hstack, vstack, zeros

class Metabolism(object):
    """Implements the representation of a metabolic system.
    Parameter:
    'object' = ###list###
    """
    def __init__(self, reactions=None, name=''):
        super(Metabolism, self).__init__()
        self.name = name
        if reactions:
            self.reactions = list(reactions)
        self.compounds = set()
        for rxn in self.reactions:
            self.compounds.update(rxn.get_compounds()) 
        # self.reactions_dictionary = dict([(reac.id, reac) for reac in reactions])
        # self.metabolite_centric = MetaboliteCentricNetwork(self.reactions, name=self.name)
    
    def __str__(self):
        """Provides some statistics about the system e.g. no. of reactions"""
        info = """System name: %s
Number of Reactions: %i
Number of Metabolites: %i""" 
        return info % (self.name, len(self), len(self.get_compounds()))
    
    def __len__(self):
        """System size.
        Returns the number of reactions."""
        return len(self.reactions)
    
    def __contains__(self, reaction):
        if isinstance(reaction, Reaction):
            if reaction in self.reactions:
                return True
            else:
                False
        else:
            if reaction in self.reactions_dictionary:
                return True
            else:
                return False
    
    def __getitem__(self, indexkey):
        """Provides item access
        
        Both access via keys (dict) and indices (list, tuple, ...) is allowed.
        Indices have to be integers
        Keys have to be something else
        """
        if type(indexkey) == int:
            return self.reactions[indexkey]
        else:
            try:
                return self.reactions_dictionary[indexkey]
            except KeyError, msg:
                print " ".join(("Reaction",indexkey,"is not in the system!"))
                raise KeyError(msg)
    
    def get_compounds(self):
        return list(self.compounds())
        
    def get_stoichiometry_matrix(self):
        """"""
        # return StoichiometricMatrix()
        pass


class StoichiometricMatrix(object):
    """A class representing a stoichiometric matrix.
    
    Columns represent reactions
    Rows represent compounds
    Coefficients ...
    """
    def __init__(self):
        super(StoichiometricMatrix, self).__init__()
        self.matrix = None
        self.compound_map = dict()
        self.reaction_map = dict()
        
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
                self.matrix = hstack((self.matrix, zeros((self.matrix.shape[0], 1))))
                for compound in reaction:
                    if compound not in self.compound_map:
                        self.compound_map[compound] = i
                        i += 1
                        self.matrix = vstack((self.matrix, zeros((1, self.matrix.shape[1]))))
                    self.matrix[self.compound_map[compound]][self.reaction_map[reaction]] = -reaction.get_stoich_coeff(compound)
    
    def _init_matrix(self, reaction):
        num = len(reaction)
        self.matrix = zeros((num, 1))
        for i, compound in enumerate(reaction):
            self.compound_map[compound] = i
            self.matrix[i][0] = reaction.get_stoich_coeff(compound)
        self.reaction_map[reaction] = 0
        return num
        
    
class Compound(object):
    """A class modeling a chemical compound.
    
    Parameters:
        id = compound id
        ...
    """
    def __init__(self, id, compartment=None, synonyms=None, formula=None, InChI=None, InCheyKey=None, SMILES=None, Charge=None, Mass=None):
        super(Compound, self).__init__()
        self.id = id
        self.compartment = compartment
        self.synonyms = synonyms
        self.formula = formula
        self.InChI = InChI
        self.InCheyKey = InCheyKey
        self.SMILES = SMILES
        self.Charge = Charge
        self.Mass = Mass
        
    def __str__(self):
        return str(id)
        
    def __contains__(self, element):
        """Checks for atomic element in compound."""
        raise NotImplemented
    
    def __hash__(self):
        return hash(self.id)
    
    def __cmp__(self, other):
        return hash(self.id) - hash(other.id)


class Reaction(object):
    """A class modeling a chemical reaction.
    
    id: The identifier of the reaction (integers are converted to strings!)
    substrates: A tuple of the substrates of type Compound
    products: A tuplle of the products of type Compound
    stoichiometry: A tuple of the stoichiometric factors (floats and integers)
        e.g. 2 A + 4.3 B -> 1 C => [2, 4.3, 1]
    reversibleQ = If the reaction is reversible [False]
    """
    def __init__(self, id, substrates, products, stoichiometry, reversibleQ=False):
        super(Reaction, self).__init__()
        try:
            hash(id)
        except TypeError, msg:
            print "Sorry, the provided reaction id cannot be used as a dictionary key!"
            raise TypeError(msg)
        if type(id) == int:
            self.id = str(id)
        else:
            self.id = id
        self.substrates = tuple(substrates)
        self.products = tuple(products)
        self.stoichiometry = tuple(stoichiometry)
        self.stoichiometry_dict = dict(zip(list(self.substrates) + list(products), self.stoichiometry))
        self.reversibleQ = reversibleQ
        self.__consistency_check()
    
    def __consistency_check(self):
        """Makes some consistency checks.
        
        1. Equal No. of substrates + products and stoichiometric factors
        2. Check for stoichiometric balancing (if enough meta information e.g. elemental composition, mass etc is available)
        ...
        """
        assert len(self.products) + len(self.substrates) == len(self.stoichiometry), "The number of stoichimetric factory does not match the number of compounds"
        assert (set(self.products) & set(self.substrates)) == set([])
    
    def __iter__(self):
        return (elem for elem in self.substrates + self.products)
    
    def __str__(self):
        """Print the reaction
        e.g. 2 A + 4.3 B -> 1 C or 2 A + 4.3 B <=> 1 C for reversible reactions.
        """
        def util(compound_list):
            reactionStringList = list()
            for compound in compound_list:
                reactionStringList.append(str(stoichTmp.pop(0)))
                reactionStringList.append(compound.id)
                if not compound == compound_list[-1]:
                    reactionStringList.append('+')
            return reactionStringList
        stoichTmp = list(self.stoichiometry)
        reactionStringList = list()
        reactionStringList += util(self.substrates)
        if self.reversibleQ:
            reactionStringList.append('<=>')
        else:
            reactionStringList.append('->')
        reactionStringList += util(self.products)
        return ' '.join(reactionStringList)
    
    def __contains__(self, compound):
        if isinstance(compound, Compound):
            if compound in self.substrates or compound in self.products:
                return True
            else:
                False
        else:
            flag = False
            for item in self.products:
                if compound == item.id:
                    flag = True
            for item in self.substrates:
                if compound == item.id:
                    flag = True
            return flag
    
    def __len__(self):
        return len(self.substrates) + len(self.products)
    
    def get_compounds(self):
        """Return both """
        return tuple(list(self.substrates) + list(self.products))
    
    def get_stoich_coeff(self, compound):
        """docstring for getStoichiometricFactor"""
        try:
            if compound in self.substrates:
                return - self.stoichiometry_dict[compound]
            elif compound in self.products:
                return self.stoichiometry_dict[compound]
        except KeyError, msg:
            print ' '.join((compound.id,"is not participating in reaction", reaction.id))
            raise KeyError(msg)


# class BipartiteMetabolicNetwork(networkx.DiGraph):
#     """docstring for Bipartite"""
#     def __init__(self, reactions=None, name='', weighted=True):
#         super(BipartiteMetabolicNetwork, self).__init__(name=name, weighted=weighted)
#         if reactions:
#             self.__populate_graph_on_init(reactions)
#     
#     def __populate_graph_on_init(self, reactions):
#         for reac in reactions:
#             self.add_reaction(reac)
#     
#     def add_reaction(self, reaction):
#         """docstring for add_reaction"""
#         subEdges = [(elem, reaction, reaction.getStoichiometricFactor(elem)) for elem in reaction.substrates]
#         prodEdges = [(reaction, elem, reaction.getStoichiometricFactor(elem)) for elem in reaction.products]
#         if reaction.reversibleQ:
#             return subEdges + prodEdges + [(j, i, d) for i, j, d in subEdges] + [(j, i, d) for i, j, d in prodEdges]
#         else:
#             return subEdges + prodEdges

class MetaboliteCentricNetwork(networkx.MultiDiGraph):
    """"""
    def __init__(self, reactions=None, name='', weighted=True):
        super(MetaboliteCentricNetwork, self).__init__(name=name, weighted=weighted)
        if reactions:
            self.__populate_graph_on_init(reactions)
    
    def __populate_graph_on_init(self, reactions):
        for reac in reactions:
            self.add_reaction(reac)
    
    def add_reaction(self, reaction):
        """docstring for add_reaction"""
        edges = [(substr, prod, reaction) for prod in reaction.products for substr in reaction.substrates]
        if reaction.reversibleQ:
            self.add_edges_from(edges + [(j, i, d) for i, j, d in edges])
        else:
            self.add_edges_from(edges)


if __name__ == '__main__':
    # import io
    # rxns = io.read_reactions_from_csv(open('./test_data/iAF1260_cytsolicNet_CurrencyFree.csv', 'rU'))
    # met_system = MetabolicSystem(rxns, name='test_system')
    # print met_system

    print Compound('atp') == Compound('atp')
    print hash(Compound('atp'))
    print hash(Compound('atp'))
    rxn1 = Reaction('ATPhyrdolysis', (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (1,1,1,1), reversibleQ=True)
    rxn2 = Reaction('GTPhyrdolysis', (Compound('gtp'), Compound('h2o')), (Compound('gdp'), Compound('pi')), (1,1,1,1), True)
    rxn3 = Reaction('CTPhyrdolysis', (Compound('ctp'), Compound('h2o')), (Compound('cdp'), Compound('pi')), (1,1,1,1))    
    sys = Metabolism((rxn1, rxn2, rxn3))
    s = StoichiometricMatrix()
    s.add_stoichiometry_from(sys)
    print s.matrix
    print s.compound_map
    print s.reaction_map
    # print rxn3 in sys
    # print 'ATPhyrdolysis' in sys
    # 
    # testGraph = MetaboliteCentricNetwork((rxn1, rxn2, rxn3), name='testGraph')
    
