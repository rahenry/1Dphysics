#!/usr/bin/env python3

import numpy
import scipy

def eigenvalue_test(m, x):
    y = numpy.matmul(m, x)
    r = y[0] / x[0]
    n = numpy.inner(numpy.conjugate(x), y) / r
    n2 = numpy.inner(x, y) / r
    return r, n, n2
