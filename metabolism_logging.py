"""
Logging of information from the pyMetabolism package is organised in a
hierarchical fashion. pyMetabolism uses L{logging} for this pupose. If you want
to use the logging facilities for your own purposes simply create your own top-
level logger called "pyMetabolism" and then attach a handler of your choosing to
it. All package internal loggers inherit from "pyMetabolism" so make sure not to
over-ride your logger instance's attribute perpetuate to C{False}. A simple code
example of using logging with pyMetabolism is below.

>>> import logging
>>> logger = logging.getLogger("pyMetabolism")
>>> handler = logging.StreamHandler()
>>> handler.setLevel(logging.DEBUG)
>>> logger.addHandler(handler)

Now all messages from within the pyMetabolism modules are propagated up to this
logger and then handled by your handler. If you want more fine-grained control
you can access the sub-loggers and attach other handlers to them. The logger
inheritance path is always as follows:

pyMetabolism.Class.instance

That means if you want all messages occurring in reactions sent to you by e-mail
because they're so pretty then you can get the logger easily:

>>> reaction_logger = logging.getLogger("pyMetabolism.Reaction")

Just add your preferred handler to this logger and you're good to go.

@author: Nikolaus Sonnenschein
@author: Moritz Beber
@copyright: Jacobs University Bremen. All rights reserved.
@since: 2010-02-12
"""


import logging


class NullHandler(logging.Handler):
    """
    Stub handler that does not emit any messages. Used in library development.
    Any application should add its own handler to the library logger or root
    logger.
    """
    def emit(self, record):
        pass
