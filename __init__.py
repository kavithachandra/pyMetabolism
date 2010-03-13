"""
A python package to easily represent metabolism as an ensemble of reactions
which describe transformations of compounds. It aims to provide intuitive
interfaces using L{Metabolism} objects in flux balance analysis (FBA),
statistics (random), graph theory (graph_view), and in connecting metabolism
with transcription data (PDA).

@author: Nikolaus Sonnenschein
@author: Moritz Beber
@contact: niko(at)foo.org
@contact: moritz(at)foo.org
@copyright: Jacobs University Bremen. All rights reserved.
@since: 2009-09-11
@version: 0.1
@status: development
@todo: fba interface, pda module

"""


import subprocess
import sys


class OptionsManager(object):
    """
    This class unifies some global options, like certain name pre- and suffixes,
    number of cpus to use for certain parallel tasks, etc.
    """
    def __init__(self):
        self.n_cpus = 1
        self.metb_prefix = "M_"
        self.rxn_prefix = "R_"
        self.rev_rxn_suffix = "_Rev"
        self._find_cpus()
        self.logger_main = "pyMetabolism"

    def _find_cpus(self):
        pass

# eof