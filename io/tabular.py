#!/usr/bin/env python
# encoding: utf-8
"""
tabular.py

Created by Nikolaus Sonnenschein on 2010-02-02.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import csv
import re
from pyMetabolism.metabolism import Reaction, Compound

test_reac = '1.3 Mcrnc + 4 Mctbtcoac <=> Mcrncoac'

def parse_reaction_from_string(reaction_string, id='',rev='<=>', irev='->'):
    """docstring for parse_reacs_from_string"""
    def get_compounds_and_stoich(string):
        reacs_plus_stoich = dict()
        subElems = re.split('\s*\+\s*', string)
        for elem in subElems:
            stoich = re.findall('[\s*]*(\d+\.*\d*)[\s*]+', elem)
            if stoich == []:
                stoich.append('1.')
            # elem2 = elem.replace(stoich[0], '')
            elem2 = re.sub('[\s*]*'+stoich[0]+'[\s*]+', '', elem)
            compounds = re.findall('\s*(\w+)\s*', elem2)
            if len(stoich) > 1 or len(compounds) > 1 or len(compounds) < 1:
                raise Exception, 'The reaction string "' + string +'" is probably malformed!'
            reacs_plus_stoich[compounds[0]] = abs(float(stoich[0]))
        return reacs_plus_stoich
    if re.search(rev, reaction_string):
        revQ = True
    else:
        revQ = False
    if revQ:
        (substr_part, prod_part) = reaction_string.split(rev)
    else:
        (substr_part, prod_part) = reaction_string.split(irev)
    tmp = [get_compounds_and_stoich(string) for string in (substr_part, prod_part)]
    return Reaction(id, [Compound(elem) for elem in tmp[0].keys()], [Compound(elem) for elem in tmp[1].keys()], tmp[0].values() + tmp[1].values() ,reversibleQ=revQ)

def read_reactions_from_csv(filehandle, id_column=0, reac_column=1, rev='<=>', irev='->', delimiter=';'):
    """returns a list of Reactions"""
    content = csv.reader(filehandle, delimiter=delimiter)
    reaction_list = list()
    content = list(content)
    print content[1]
    for row in content[1:-2]:
        reaction_list.append(parse_reaction_from_string(row[reac_column], id=row[id_column], rev=rev, irev=irev))
    return reaction_list

if __name__ == '__main__':
    print parse_reaction_from_string('2.3 atp + 44 ctp <=> 5 gtp')
    for elem in read_reactions_from_csv(open('./test_data/iAF1260_cytsolicNet_CurrencyFree.csv', 'rU')):
        print elem
    print parse_reaction_from_string(test_reac, id='testReac')

