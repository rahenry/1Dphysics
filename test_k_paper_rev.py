#!/usr/bin/env python3

from fp_k import *
from System import System
from solve import *

L = 4
M = 2*L+1
lamb = 1
phi = 0.1
ind = 0
for a in range(1,2*L+1):
    k = numpy.pi * a / M
    lf = numpy.exp(2.j*numpy.pi/3*phi)
    u = f_kj(k, lf, L, 3)
    eps = eps_k(k, lf, 3)
    # print(a, k, u)
    s = {"L" : L, "N":3, "lambda" : lamb, "phi":phi, "mu":-1}
    s = System(s)
    s.prepare()
    solve_matrix(s)

    if abs(u) < 1E-10:
        print(a, eps, s['pfs'][-ind], eps)
        ind += 1
