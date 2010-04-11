# encoding: utf-8


"""
Some common functions.

@author: Nikolaus Sonnenschein
@author: Moritz Beber
@contact: niko(at)foo.org
@contact: moritz(at)foo.org
@copyright: Jacobs University Bremen. All rights reserved.
@since: 2010-04-11
@todo: complete docs
"""


class OutputEater(object):
     def write(self, string):
         pass


def gcd(a, b):
    """
    Greatest commond divisor
    """
    while a:
        (a, b) = (b % a, a)
    return b

def lcm(a, b):
    """
    Least common multiple
    """
    return (a * b / gcd(a, b))

def fxrange(*args):
    """
    @author: dwhall256
    @see: http://code.activestate.com/recipes/66472-frange-a-range-function-with-float-increments/
    """
    start = 0.0
    step = 1.0

    l = len(args)
    if l == 1:
        end = args[0]
    elif l == 2:
        start, end = args
    elif l == 3:
        start, end, step = args
        if step == 0.0:
            raise ValueError, "step must not be zero"
    else:
        raise TypeError, "frange expects 1-3 arguments, got %d" % l

    v = start
    while True:
        if (step > 0 and v >= end) or (step < 0 and v <= end):
            raise StopIteration
        yield v
        v += step


if __name__ == "__main__":
    print lcm(4, 2)
    print lcm(50, 30)
    print list(fxrange(0.1, 1.0, 0.1))