#!/usr/bin/env python3

from exact_fp import *
import numpy


N = 3
flipZ = 0
def ratio(gamma):
    lam = gamma / (1.-gamma)
    e, x, pfs = exact_parafermions_matrix(100, N, lam, flipZ)
    scale = 1./(1.-gamma)
    pfs = [p / scale for p in pfs]
    return pfs[0] / pfs[-1]

def f1(gammas, target):
    return [ratio(gamma) - target for gamma in gammas]

from scipy.optimize import root
r = root(f1, [0.5], args=(0.5))
print(r)
print(r.x)
print(ratio(r.x[0]))
