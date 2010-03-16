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

import logging
import os
from pyMetabolism.metabolism_logging import NullHandler
import subprocess


def new_property(func):
    """
    A decorator function for easy property creation.
    @author: runsun pan
    @copyright: MIT License
    @contact: http://code.activestate.com/recipes/576742/
    """
    ops = func() or dict()
    name = ops.get("prefix", '_') + func.__name__ # property name
    fget = ops.get("fget", lambda self: getattr(self, name))
    fset = ops.get("fset", lambda self, value: setattr(self, name, value))
#    fdel = ops.get("fdel", lambda self: delattr(self, name))
#    return property(fget=fget, fset=fset, fdel=fdel, doc=ops.get('doc',''))
    return property(fget, fset)


class OptionsManager(object):
    """
    This class unifies some global options, like certain name pre- and suffixes,
    number of cpus to use for certain parallel tasks, etc. This class is modelled
    as a singleton, i.e., only one instance of this class may exist.
    @todo: Consider thread safety, though I would defer that to using this class.
    """
    _singleton = None

    def __new__(cls, * args, ** kwargs):
        if not cls._singleton:
            cls._singleton = super(OptionsManager, cls).__new__(cls, * args, ** kwargs)
        return cls._singleton

    def __init__(self, * args, ** kwargs):
        super(OptionsManager, self).__init__(*args, ** kwargs)
        self._metb_prefix = "M_"
        self._rxn_prefix = "R_"
        self._rev_rxn_suffix = "_Rev"
        self._main_logger_name = "pyMetabolism"
        self._logger = logging.getLogger(self.main_logger_name)
        self.logger.addHandler(NullHandler)
        self._n_cpus = self._find_num_cpus()

    def _find_num_cpus(self):
        """
        Detects the number of effective CPUs in the system. Adapted from Parallel
        Python.
        @author: Vitalii Vanovschi
        @copyright: (c) 2005-2010. All rights reserved.
        @contact: http://www.parallelpython.com
        """
        num_cpus = None
        # for Linux, Unix and MacOS
        if hasattr(os, "sysconf"):
            if "SC_NPROCESSORS_ONLN" in os.sysconf_names:
                #Linux and Unix
                try:
                    num_cpus = int(os.sysconf("SC_NPROCESSORS_ONLN"))
                except:
                    pass
                else:
                    if num_cpus > 0:
                        return num_cpus
        # Linux, Unix failsafe
        cmd = ["grep", "-c", "'model name'", "'/proc/cpuinfo'"]
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        if p.returncode == 0:
            try:
                num_cpus = int(stdout)
            except:
                pass
            else:
                if num_cpus > 0:
                    return num_cpus
        # MacOS X
        cmd = ["sysctl", "-n", "hw.ncpu"]
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        if p.returncode == 0:
            try:
                num_cpus = int(stdout)
            except:
                pass
            else:
                if num_cpus > 0:
                    return num_cpus
        #for Windows
        if "NUMBER_OF_PROCESSORS" in os.environ:
            try:
                num_cpus = int(os.environ["NUMBER_OF_PROCESSORS"])
            except:
                pass
            else:
                if num_cpus > 0:
                    return num_cpus
        return 1

    @new_property
    def metb_prefix():
        pass

    @new_property
    def rxn_prefix():
        pass

    @new_property
    def rev_rxn_suffix():
        pass

    @new_property
    def logger():
        pass

    @new_property
    def main_logger_name():
        pass

    @new_property
    def n_cpus():
        pass
