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
from pyMetabolism import *
from io import *


class ioTests(unittest.TestCase):
    def setUp(self):
        pass

    def test_read_reactions_from_csv(self):
        # rxns = read_reactions_from_csv(open('./test_data/iAF1260_cytsolicNet_CurrencyFree.csv', 'rU'))
        # met_system = MetabolicSystem(rxns, name='test_system')

        self.assertRaises(Exception, parse_reaction_from_string, '2.3 atp + 44 ctp <=> 5 gtp')
        # self.assertRaises(Exception, parse_reaction_from_string, ('2.3 atp hasdhf + 44 ctp <=> 5 gtp',))

if __name__ == '__main__':
    print parse_reaction_from_string('2.3 atp + 44 ctp <=> 5 gtp')
    suite = unittest.TestLoader().loadTestsFromTestCase(ioTests)
    unittest.TextTestRunner(verbosity=4).run(suite)
