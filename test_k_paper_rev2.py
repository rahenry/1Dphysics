#!/usr/bin/env python3

from SystemSet import SystemSet
from fp_k import *

S = SystemSet("test_k_paper_rev")
S.prepare()
S.run()

def f_kj(k, lamb, L, N):
    return numpy.sin((L+1.) * k) + (lamb ** (-N/2.)) * numpy.sin(L*k)

def eps_to_k(e, g, N):
    if g < 0:
        g = complex(g)
    a = 1
    res = (e ** N) - 1. - (g ** (a*N))
    res /= 2. * (g ** ((a*N)/2.))
    # res = -res
    # res = abs(res)
    r1 = res
    res = numpy.arccos(res)
    return res
for s in S.systems:
    pfs = s['pfs']
    lamb = s['lambda_full']
    N = s['N']
    L = s['L']
    print(s, lamb)
    for p in pfs:
        print("...")
        k = eps_to_k(p, lamb, N)
        u = f_kj(k, lamb, L, N)
        # print(p, k/numpy.pi, u)
        print(p, k, u)
        p2 = eps_k(k, lamb, N)
        print(abs(p), abs(p2))
