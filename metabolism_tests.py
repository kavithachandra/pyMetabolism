#!/usr/bin/env python
# encoding: utf-8
"""
pyMetabolism_Tests.py

Created by Nikolaus Sonnenschein on 2009-09-10.
Copyright (c) 2009 . All rights reserved.
"""

import unittest
from pyMetabolism.metabolism import Metabolism, Compound, Reaction, Compartment, CompartCompound

class Item5MetabolismTests(unittest.TestCase):
    def setUp(self):
        self.r1 = Reaction('ATPhyrdolysis', (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (-1,-1,1,1))
        self.r2 = Reaction('GTPhyrdolysis', (Compound('gtp'), Compound('h2o')), (Compound('gdp'), Compound('pi')), (-1,-1,1,1))
        self.r3 = Reaction('CTPhyrdolysis', (Compound('ctp'), Compound('h2o')), (Compound('cdp'), Compound('pi')), (-1,-1,1,1))
        self.metabolic_system = Metabolism((self.r1, self.r2, self.r3))
        self.smaller_metabolic_system = Metabolism((self.r1, self.r2))
    
    # def test__str__(self):
    #     """Tests if the __str__ methods provides the correct statistics"""
    #     self.assertTrue(type(self.metabolic_system.__str__()) == str)
        
    def test__getitem__(self):
        """Test if metabolic_system[y] provides the correct x"""
        self.assertEqual(self.metabolic_system['ATPhyrdolysis'], self.r1)
        self.assertEqual(self.metabolic_system[0], self.r1)
        self.assertEqual(self.metabolic_system['GTPhyrdolysis'], self.r2)
        self.assertEqual(self.metabolic_system[1], self.r2)

    def test__contains__(self):
        """Tests if __contains__ works correctly"""
        self.assertTrue(self.r1 in self.metabolic_system)
        self.assertTrue(self.r2 in self.metabolic_system)
        self.assertTrue(self.r3 in self.metabolic_system)
        self.assertFalse(self.r3 in self.smaller_metabolic_system)
        
        self.assertTrue('ATPhyrdolysis' in self.metabolic_system)
        self.assertTrue('GTPhyrdolysis' in self.metabolic_system)
        self.assertTrue('CTPhyrdolysis' in self.metabolic_system)
        self.assertFalse('CTPhyrdolysis' in self.smaller_metabolic_system)
        
    
    
    # def test__add__(self):
    #     """Tests if it is possible to merge metabolic systems"""
    #     pass
    # 
    # def test__del__(self):
    #     """docstring for test__del__"""
    #     pass
    # 
    # def test__cmp__(self):
    #     """docstring for test__"""
    #     pass



class Item4ReactionTests(unittest.TestCase):
    def setUp(self):
        self.reaction = Reaction('ATPhyrdolysis', (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (1,1,1,1))
    
    # def test__init__(self):
    #     """Tests if wrong initializations are caught by the right Exceptions"""
    #     self.assertRaises(TypeError, self.reaction.__init__, ('asdf', (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (1,1,1,1)))
    #     self.assertRaises(TypeError, self.reaction.__init__, ([1,2], (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (1,1,1,1)))
    
    def test__str__(self):
        """Tests if the __str__ method works correctly"""
        self.assertEqual(self.reaction.__str__(), '1 atp + 1 h2o -> 1 adp + 1 pi')

    def test_contains_string_input(self):
        self.assertTrue('atp' in self.reaction)
        
    def test_contains_string_input(self):
        self.assertTrue(Compound('atp') in self.reaction)


class Item2CompartmentTests(unittest.TestCase):
    def setUp(self):
        args = {"suffix": "_c", "spatial_dimensions": 3, "size": 1., "units": "ml"}
        self.comp = Compartment("Cytosol", True, args)

    def test_options(self):
        def utility():
            self.comp.options = None
        self.assertRaises(AttributeError, utility())


class Item3CompartCompoundTests(unittest.TestCase):
    pass


class Item1CompoundTests(unittest.TestCase):
    pass


if __name__ == '__main__':
    tests = list()
    tests.append(unittest.TestLoader().loadTestsFromTestCase(Item1CompoundTests))
    tests.append(unittest.TestLoader().loadTestsFromTestCase(Item2CompartmentTests))
    tests.append(unittest.TestLoader().loadTestsFromTestCase(Item3CompartCompoundTests))
    tests.append(unittest.TestLoader().loadTestsFromTestCase(Item4ReactionTests))
    tests.append(unittest.TestLoader().loadTestsFromTestCase(Item5MetabolismTests))
    suite = unittest.TestSuite(tests)
    unittest.TextTestRunner(verbosity=4).run(suite)