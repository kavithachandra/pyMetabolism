#!/usr/bin/env python
# encoding: utf-8
"""
io_Tests.py

Created by Nikolaus Sonnenschein on 2009-09-11.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import unittest
from pyMetabolism.metabolism import *
from pyMetabolism.io.tabular import *
from pyMetabolism.io.sbml import *


class ioTests(unittest.TestCase):
    def setUp(self):
        pass

    def test_read_reactions_from_csv(self):
        # rxns = read_reactions_from_csv(open('./test_data/iAF1260_cytsolicNet_CurrencyFree.csv', 'rU'))
        # met_system = MetabolicSystem(rxns, name='test_system')

        self.assertRaises(Exception, tabular.parse_reaction_from_string, '2.3 atp + 44 ctp <=> 5 gtp')
        # self.assertRaises(Exception, parse_reaction_from_string, ('2.3 atp hasdhf + 44 ctp <=> 5 gtp',))
        

class sbmlTests(unittest.TestCase):
    """docstring for sbmlTests"""
    def setUp(self):
        self.smallModel = '../test_data/Ec_core_flux1.xml'
        self.bigModel = '../test_data/iAF1260.xml'
        self.parser_ecoli_core = SBMLParser('../test_data/Ec_core_flux1.xml')
        self.parser_iAF1260 = SBMLParser('../test_data/iAF1260.xml')
        
    def test_sbml_parser_get_compounds(self):
        """Tests if compounds in SBML models are correctly read."""
        compounds = self.parser_ecoli_core.get_compounds()
        self.assertTrue(isinstance(compounds[0], Compound))
        self.assertTrue(isinstance(compounds[-1], Compound))
        self.assertEqual(compounds[0].get_id(), 'M_13dpg_c')
        self.assertEqual(compounds[-1].get_id(), 'M_succ_b')
        self.assertEqual(len(compounds), 77)
        compounds = self.parser_iAF1260.get_compounds()
        self.assertTrue(isinstance(compounds[0], Compound))
        self.assertTrue(isinstance(compounds[-1], Compound))
        self.assertEqual(compounds[0].get_id(), 'M_10fthf_c')
        self.assertEqual(compounds[-1].get_id(), 'M_zn2_b')
        self.assertEqual(len(compounds), 1972)
        
    def test_sbml_get_reactions(self):
        """Tests if reactions in SBML models are correctly read."""
        reactions = self.parser_ecoli_core.get_reactions()
        self.assertTrue(isinstance(reactions[0], Reaction))
        self.assertTrue(isinstance(reactions[-1], Reaction))
        self.assertEqual(reactions[0].get_id(), 'R_ACKr')
        self.assertEqual(reactions[-1].get_id(), 'R_TPI')
        self.assertEqual(len(reactions), 77)
        reactions = self.parser_iAF1260.get_reactions()
        self.assertTrue(isinstance(reactions[0], Reaction))
        self.assertTrue(isinstance(reactions[-1], Reaction))
        self.assertEqual(reactions[0].get_id(), 'R_GSPMDS')
        self.assertEqual(reactions[-1].get_id(), 'R_FDMO6')
        self.assertEqual(len(reactions), 2381)
        
    def test_get_metabolic_systme(self):
        """Tests SBML models are correctly imported"""
        metbol = self.parser_ecoli_core.get_metabolic_system()
        self.assertEqual(metbol.name, 'Ec_core')
        metbol = self.parser_iAF1260.get_metabolic_system()
        self.assertEqual(metbol.name, 'E. coli iAF1260')
        
        
if __name__ == '__main__':
    # suite = unittest.TestLoader().loadTestsFromTestCase(ioTests)
    # unittest.TextTestRunner(verbosity=4).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(sbmlTests)
    unittest.TextTestRunner(verbosity=10).run(suite)
