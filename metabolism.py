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
    _metabolism_memory = dict()
    _metabolism_counter = 0
    
    def __new__(cls, reactions=None, name='', *args, **kwargs):
        if name:
            name = str(name)
        if name in cls._metabolism_memory:
            return cls._metabolism_memory[name]
        else:
            instance = super(Metabolism, cls).__new__(cls, *args, **kwargs)
            return instance
    
    def __init__(self, reactions=None, name=None, *args, **kwargs):
        if name in self.__class__._metabolism_memory:
            return None
        super(Metabolism, self).__init__(*args, **kwargs)
        if name:
            self.name = name
        else:
            self.name = "Metabolism-%d" % self.__class__._metabolism_counter
            self.__class__._metabolism_counter += 1
        if reactions:
            self.reactions = list(reactions)
        self.compounds = set()
        for rxn in self.reactions:
            self.compounds.update(rxn.get_compounds()) 
        # self.reactions_dictionary = dict([(reac.id, reac) for reac in reactions])
        # self.metabolite_centric = MetaboliteCentricNetwork(self.reactions, name=self.name)
        self.__class__._metabolism_memory[self.name] = self
    
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
            reaction = str(reaction)
            for item in self.reactions:
                if reaction == item.identifier:
                    return True
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
    
    def __cmp__(self, other):
        return cmp(id(self), id(other))
    
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
    def __init__(self, *args, **kwargs):
        super(StoichiometricMatrix, self).__init__(*args, **kwargs)
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
                self.matrix = hstack((self.matrix, zeros((self.matrix.shape[0],
                    1))))
                for compound in reaction:
                    if compound not in self.compound_map:
                        self.compound_map[compound] = i
                        i += 1
                        self.matrix = vstack((self.matrix, zeros((1,
                            self.matrix.shape[1]))))
                    self.matrix[self.compound_map[compound]]\
                        [self.reaction_map[reaction]] =\
                        -reaction.get_stoich_coeff(compound)
    
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
    _compound_memory = dict()
    
    def __new__(cls, identifier, compartment=None, synonyms=None,
        formula=None, in_chi=None, in_chey_key=None, smiles=None, charge=None,
        mass=None, *args, **kwargs):
        identifier = str(identifier)
        if identifier in cls._compound_memory:
            return cls._compound_memory[identifier]
        else:
            instance = super(Compound, cls).__new__(cls, *args, **kwargs)
            return instance
    
    def __init__(self, identifier, compartment=None, synonyms=None,
        formula=None, in_chi=None, in_chey_key=None, smiles=None, charge=None,
        mass=None, *args, **kwargs):
        if identifier in self.__class__._compound_memory:
            return None
        super(Compound, self).__init__(*args, **kwargs)
        self.identifier = identifier
        self.compartment = compartment
        self.synonyms = synonyms
        self.formula = formula
        self.in_chi = in_chi
        self.in_chey_key = in_chey_key
        self.smiles = smiles
        if charge:
            self.charge = int(charge)
        else:
            self.charge = None
        if mass:
            self.mass = float(mass)
        else:
            self.mass = None
        self.__class__._compound_memory[self.identifier] = self
        
    def __str__(self):
        return self.identifier
        
    def __contains__(self, element):
        """Checks for atomic element in compound."""
        raise NotImplemented
    
    def __hash__(self):
        return hash(self.identifier)
    
    def __cmp__(self, other):
        return cmp(id(self), id(other))


class Reaction(object):
    """A class modeling a chemical reaction.
    
    id: The identifier of the reaction (integers are converted to strings!)
    substrates: A tuple of the substrates of type Compound
    products: A tuplle of the products of type Compound
    stoichiometry: A tuple of the stoichiometric factors (floats and integers)
        e.g. 2 A + 4.3 B -> 1 C => [2, 4.3, 1]
    reversibleQ = If the reaction is reversible [False]
    """
    _reaction_memory = dict()
    
    def __new__(cls, identifier, substrates, products, stoichiometry,
        reversible=False, synonyms=None, *args, **kwargs):
        identifier = str(identifier)
        if identifier in cls._reaction_memory:
            return cls._reaction_memory[identifier]
        else:
            instance = super(Reaction, cls).__new__(cls, *args, **kwargs)
            return instance
    
    def __init__(self, identifier, substrates, products, stoichiometry,
        reversible=False, synonyms=None, *args, **kwargs):
        if identifier in self.__class__._reaction_memory:
            return None
        super(Reaction, self).__init__(*args, **kwargs)
        self.identifier = identifier
        try:
            hash(self.identifier)
        except TypeError, msg:
            print "Sorry, the provided reaction id cannot be used as a dictionary key!"
            raise TypeError(msg)
        self.substrates = tuple(substrates)
        self.products = tuple(products)
        self.stoichiometry = tuple([int(coeff) for coeff in stoichiometry])
        self.stoichiometry_dict = dict(zip(list(self.substrates)
            + list(self.products), self.stoichiometry))
        self.reversible = bool(reversible)
        self._consistency_check()
        self.__class__._reaction_memory[self.identifier] = self
    
    def _consistency_check(self):
        """Makes some consistency checks.
        
        1. Equal No. of substrates + products and stoichiometric factors
        2. Check for stoichiometric balancing (if enough meta information e.g.
        elemental composition, mass etc is available)
        ...
        """
        assert (len(self.products) + len(self.substrates)) ==\
            len(self.stoichiometry), "The number of stoichimetric factory does"\
            " not match the number of compounds"
        assert (set(self.products) & set(self.substrates)) == set([])
    
    def __iter__(self):
        return (elem for elem in self.substrates + self.products)
    
    def __str__(self):
        """Print the reaction
        e.g. 2 A + 4.3 B -> 1 C or 2 A + 4.3 B <=> 1 C for reversible reactions.
        """
        def util(compound_list):
            reaction_string = list()
            for compound in compound_list:
                reaction_string.append(str(self.stoichiometry_dict[compound]))
                reaction_string.append(compound.identifier)
                if not (compound == compound_list[-1]):
                    reaction_string.append('+')
            return reaction_string
        reaction_string = list()
        reaction_string.extend(util(self.substrates))
        if self.reversible:
            reaction_string.append('<=>')
        else:
            reaction_string.append('->')
        reaction_string.extend(util(self.products))
        return ' '.join(reaction_string)
    
    def __contains__(self, compound):
        if isinstance(compound, Compound):
            if compound in self.substrates or compound in self.products:
                return True
        else:
            compound = str(compound)
            for item in self.products:
                if compound == item.identifier:
                    return True
            for item in self.substrates:
                if compound == item.identifier:
                    return True
        return False
    
    def __len__(self):
        return len(self.substrates) + len(self.products)
    
    def __cmp__(self, other):
        return cmp(id(self), id(other))
    
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
            print ' '.join((compound.idenifier, "is not participating in"\
                " reaction", reaction.identifier))
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

class MetaboliteCentricNetwork(networkx.DiGraph):
    """"""
    def __init__(self, reactions=None, name='', weighted=True, *args, **kwargs):
        super(MetaboliteCentricNetwork, self).__init__(name=name, weighted=weighted, *args, **kwargs)
        if reactions:
            self.__populate_graph_on_init(reactions)
    
    def __populate_graph_on_init(self, reactions):
        for reac in reactions:
            self.add_reaction(reac)
    
    def add_reaction(self, reaction):
        """docstring for add_reaction"""
        edges = [(substr, prod, reaction) for prod in reaction.products for substr in reaction.substrates]
        if reaction.reversible:
            self.add_edges_from(edges + [(j, i, d) for i, j, d in edges])
        else:
            self.add_edges_from(edges)


if __name__ == '__main__':
    # import io
    # rxns = io.read_reactions_from_csv(open('./test_data/iAF1260_cytsolicNet_CurrencyFree.csv', 'rU'))
    # met_system = MetabolicSystem(rxns, name='test_system')
    # print met_system
    rxn1 = Reaction('ATPhyrdolysis', (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (1,1,1,1), reversible=True)
    rxn2 = Reaction('GTPhyrdolysis', (Compound('gtp'), Compound('h2o')), (Compound('gdp'), Compound('pi')), (1,1,1,1), True)
    rxn3 = Reaction('CTPhyrdolysis', (Compound('ctp'), Compound('h2o')), (Compound('cdp'), Compound('pi')), (1,1,1,1))    
    sys = Metabolism((rxn1, rxn2, rxn3))
    s = StoichiometricMatrix()
    s.add_stoichiometry_from(sys)
    print s.matrix
    print rxn3 in sys
    print 'ATPhyrdolysis' in sys
    # 
    test_graph = MetaboliteCentricNetwork((rxn1, rxn2, rxn3), name='testGraph')
    print test_graph.edges()
    
