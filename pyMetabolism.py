#!/usr/bin/env python
# encoding: utf-8
"""
metabolicNetwork.py

Created by Nikolaus Sonnenschein on 2009-09-10.
Copyright (c) 2009 . All rights reserved.
"""

import sys
import os
import unittest
import networkx
import copy


class Compound(object):
    """docstring for Metabolite"""
    def __init__(self, id, synonyms=None, formula=None, InChI=None, InCheyKey=None, SMILES=None, Charge=None, Mass=None):
        super(Compound, self).__init__()
        self.id = id
        self.synonyms = synonyms
        self.formula = formula
        self.InChI = InChI
        self.InCheyKey = InCheyKey
        self.SMILES = SMILES
        self.Charge = Charge
        self.Mass = Mass

class Reaction(object):
    """A class modeling a chemical reaction.
    
    id: The identifier of the reaction
    substrates: A tuple of the substrates of type Compound
    products: A tuplle of the products of type Compound
    stoichiometry: A tuple of the stoichiometric factors (floats and integers)
        e.g. 2 A + 4.3 B -> 1 C => [2, 4.3, 1]
    reversibleQ = If the reaction is reversible [False]
    """
    def __init__(self, id, substrates, products, stoichiometry, reversibleQ=False):
        super(Reaction, self).__init__()
        self.id = id
        self.substrates = substrates
        self.products = products
        self.stoichiometry = stoichiometry
        self.reversibleQ = reversibleQ
        self.__consistencyCheck()
    
    def __consistencyCheck(self):
        """Makes some consistency checks."""
        assert len(self.products) + len(self.substrates) == len(self.stoichiometry), "The number of stoichimetric factory does not match the number of compounds"
    
    def __str__(self):
        """Print the reaction
        e.g. 2 A + 4.3 B -> 1 C or 2 A + 4.3 B <=> 1 C for reversible reactions.
        """
        stoichTmp = list(self.stoichiometry)
        reactionStringList = list()
        for compound in self.substrates:
            reactionStringList.append(str(stoichTmp.pop(0)))
            reactionStringList.append(compound.id)
        if self.reversibleQ:
            reactionStringList.append('<=>')
        else:
            reactionStringList.append('->')
        for compound in self.products:
            reactionStringList.append(str(stoichTmp.pop(0)))
            reactionStringList.append(compound.id)
        return ' '.join(reactionStringList)

class ReactionTests(unittest.TestCase):
    def setUp(self):
        self.reaction = Reaction('ATPhyrdolysis', (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (1,1,1,1))
    
    def test__str__(self):
        """docstring for test__str__"""
        self.assertEqual(self.reaction.__str__(), '1 atp 1 h2o -> 1 adp 1 pi')


class Bipartite(networkx.DiGraph):
    """docstring for Bipartite"""
    def __init__(self, name='', weighted=True):
        super(Bipartite, self).__init__(name=name, weighted=weighted)



class metabolicNetworkTests(unittest.TestCase):
    def setUp(self):
        pass


if __name__ == '__main__':
    rxn = Reaction('ATPhyrdolysis', (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (1,1,1,1))
    print rxn
    unittest.main()