#!/usr/bin/env python3


p = [ 5.59407094e-01, 2.65860512e+00,-3.41368735e+01, 6.77689501e+02]

def f(x, A, B, C, D):
    return A - B/x -C/x/x -D/x/x/x

def g(x):
    return f(x, *p)

import numpy

for x in numpy.linspace(1, 800, 30):
    print(x, g(x))
