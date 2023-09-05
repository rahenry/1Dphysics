#!/usr/bin/env python3

import sympy
import numpy
import scipy
# _, indexes = sympy.Matrix(matrix).T.rref()  # T is for transpose
def normalise(v):
    return v / numpy.vdot(v, v)
def red(u0, vs):
    vs = [normalise(x) for x in vs]
    # u = normalise(u0)
    # vs = numpy.array(vs).T
    # z = scipy.linalg.orth(vs)
    # for i in range(len(vs.T)):
    #     v = z[:,i]
    #     # print(v)
    #     u = u - numpy.vdot(u, v) * v
    # return numpy.vdot(u, u)
    u = normalise(u0)
    vs = numpy.array(vs+[u]).T
    z = scipy.linalg.orth(vs)
    # print('...')
    # print(u0)
    for i in range(len(vs)+1):
        if i >= z.shape[-1]: continue
        k = numpy.vdot(u, z[:,i])
        # print(k)
    return numpy.vdot(u, u)

def analyse_zeros(S):
    print(S)
    s0 = S.systems[0]
    for bid, b in s0.ed.blocks.items():
        if len(b.zeros) == 0: continue
        print(f'{b}, nzeros={len(b.zeros)}')
        for s in S.systems:
            print(s)
            vecs = s.ed.blocks[bid].zero_evecs
            for v in vecs:
                print(numpy.vdot(v, vecs[0]))
            for v in s.ed.blocks[bid].zero_evecs:
                # vecsys = numpy.array([x for x in b.zero_evecs]+[v])
                # print(vecsys)
                # red, indexes = sympy.Matrix(vecsys).rref()
                # print(indexes)
                # print(red)
                z = red(v, b.zero_evecs)
                # print(z)
