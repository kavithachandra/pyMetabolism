# encoding: utf-8


"""
General exceptions used in the whole pyMetabolism package.

@author: Nikolaus Sonnenschein
@author: Moritz Beber
@contact: niko(at)foo.org
@contact: moritz(at)foo.org
@copyright: Jacobs University Bremen. All rights reserved.
@since: 2010-04-11
@todo: complete docs
"""


class PyMetabolismError(StandardError):
    def __init__(self, msg="", *args, **kwargs):
        super(PyMetabolismError, self).__init__(*args, **kwargs)
#        self.errno
        self.strerror = msg

    def __str__(self):
        print self.strerror