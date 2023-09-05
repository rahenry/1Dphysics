#!/usr/bin/env python3

from exact_fp import *
import matplotlib.pyplot as plt
import matplotlib
import os

N = 3
omega = numpy.exp(2.j * numpy.pi / N)
Ls = [20, 200]
Ls = [20]
U = 101
from rah_utility import silent_remove
d = 'graphs/test_gap2'
silent_remove(d)
from rah_utility import mkdir
mkdir(d)


cm = matplotlib.cm.get_cmap('hsv')
M = 2
def find_smallest_gap(x):
    s = 1E6
    for i in range(N-1):
        g = (1. - omega ** (i+1)) * x
        if abs(g) < s:
            s = g
    return s

for L in Ls:
    gammas = numpy.linspace(0.01, 0.99, U)
    dg = 0.11
    gammas = numpy.linspace(0.5-dg, 0.5+dg, U)
    dp = 1.5 / L
    phis = numpy.linspace(0.5-dp, 0.5+dp, U)
    for m in range(M):
        i = 0
        data = numpy.zeros((len(gammas), len(phis)))
        datag = numpy.zeros((len(gammas), len(phis)))
        for g in gammas:
            j = 0
            for p in phis:
                phi = p
                gamma = g
                lamb = gamma / (1.-gamma)
                lambsc = 1./(1.-gamma)
                s = {'N' : N,
                    'lambda' : lamb,
                    'phi'  : phi,
                    'L' : L,
                    'lambda_scale' : lambsc}

                # pfs = pfs_matrix_limited(s, M)
                pfs, e0 = pfs_matrix(s)
                # gs = [abs(find_smallest_gap(x)) for x in pfs]
                # gs.sort()
                pfs = [abs(x) for x in pfs]
                pfs.sort()
                data[i][j] = abs(pfs[m]) * L
                # datag[i][j] = gs[m] * L
                j += 1
            print(i)
            i += 1

        P = plt.imshow(data, vmin=0)
        print(P)
        plt.colorbar(P)
        f = os.path.join(d, f'L={L}_m={m}.png')
        plt.savefig(f)
        plt.clf()

        # P = plt.imshow(datag)
        # plt.colorbar(P)
        # f = os.path.join(d, f'L={L}_m={m}_g.png')
        # plt.savefig(f)
        # plt.clf()
