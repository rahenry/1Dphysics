#!/usr/bin/env python3

from System import *
from SystemSet import *
from args import Q
from exact_fp import *
from rah_utility import *
from output_figs import *

import numpy, scipy


def M(x, y, z, x0, y0, z0):
    M = numpy.full((2, 2), 1.j, dtype=numpy.complex)
    M[0, 0] = z + 1.0j * z0
    M[0, 1] = x - 1.0j * y + 1.0j * x0 + y0
    M[1, 0] = x + 1.0j * y + 1.0j * x0 - y0
    M[1, 1] = -z - 1.0j * z0
    return M


r = (1, 0, 0)
r0 = (0, 1, 0.1)
M1 = M(*r, *r0)
M1d = numpy.conj(M1.T)
M1d = M1.T
print(M1)
print(M1d)
print()

solR = numpy.linalg.eig(M1)
indsR = numpy.argsort(solR[0])
print(indsR, solR[0])
solL = numpy.linalg.eig(M1d)
indsL = numpy.argsort(solL[0])
print(indsL, solL[0])

print()
for I in range(len(solR)):
    for J in range(len(solL)):
        if I < J: continue
        i = indsR[I]
        j = indsL[J]
        eR = solR[0][i]
        vR = solR[1][:,i]
        eL = solL[0][j]
        vL = solL[1][:,j]
        print(I, J, i, j)
        print(eR)
        print(eL)
        if I > J:
            print(f'overlap of right states = {numpy.dot(solR[1][:,indsR[J]], vR)}')
            print(f'{i}, {indsR[J]}')
        # print(vR)
        # print(vL)
        # vL = vL.conj()
        p = numpy.dot(vL, vR)
        print('overlap = ', abs(p), p)
        # vL = vL.conj()
        # p = numpy.dot(vL, vR)
        # print('overlap = ', abs(p), p)
        print()
