#!/usr/bin/env python3

import math
import numpy
import scipy
from exact_fp import *
from decimal import *
from scipy.special import hyp2f1
from scipy.optimize import root
from scipy.optimize import brentq
from scipy.optimize import fsolve
from rah_utility import mkdir, silent_remove, rel_diff
from rah_numpy import eigenvalue_test
import os
mkdir('test_graphs')
numpy.set_printoptions(threshold=numpy.inf)

def exact_inf(N, gamma):
    return -hyp2f1(-1./N, -1./N, 1, gamma**N)

def make_m1(L, N, gamma, flipZ):
    if flipZ: gamma *= -1
    g = gamma
    g = g ** (-0.5*N)
    g2 = g
    g1 = numpy.conj(g)
    g1 = g
    # g2 = g
    z = 1 + g1*g2
    # z = 1+abs(g)**2
    rows = [0]
    cols = [0]
    data = [1]

    for i in range(L-1):
        rows.append(i+1)
        cols.append(i)
        data.append(g1)
        rows.append(i)
        cols.append(i+1)
        data.append(g2)
        rows.append(i+1)
        cols.append(i+1)
        data.append(z)
    m = numpy.zeros((L, L), dtype=complex)
    for (i,j,d) in zip(rows, cols, data):
        if flipZ: d *= -1
        m[i, j] = d
    return m


def make_m2(L, N, gamma, flipZ):
    if flipZ: gamma *= -1
    g = gamma
    g = g ** N
    L0 = N*L
    m = numpy.zeros((L0, L0), dtype=complex)
    m = numpy.zeros((L0, L0))
    b = 0
    count = 0
    while (True):
        if b+1 >= L0: break
        m[b+1][b] = 1
        b += 1
    s = 0
    while (True):
        x = s*N
        y = s*N+N-1
        if x >= L0 or y >= L0:
            break
        else:
            count += 1
            m[x][y] = g
        s += 1
    s = 0
    while (True):
        x = s*N+1
        y = s*N+N
        if x >= L0 or y >= L0:
            break
        else:
            count += 1
            m[x][y] = 1
        s += 1
    return m

def find_closest(e, candidates):
    best = 1000
    best_choice = None
    for c in candidates:
        d = abs((c-e)/(c+e))
        if d < best:
            best = d
            best_choice = c
    return best_choice, d

    1

def inf_test(N, L, gamma, flipZ):
    print(N, L, gamma, flipZ)
    1
    m1 = make_m1(L, N, gamma, flipZ)
    sol1 = scipy.linalg.eig(m1)
    eigs1 = [e ** (1./N) * gamma / L for e in sol1[0]]
    eigs1.sort()
    e0 = 0
    for i in range(len(sol1[0])):
        e0 -= eigs1[i]
    g = gamma
    if flipZ: g = g*-1
    einf = exact_inf(N, g)
    return e0


def inf_test1():

    Ls = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    Ls = [4, 5, 10, 11]
    Ls = [12, 40, 80]
    N = 3
    gammas = [0.5, 1.0, 1.5]
    gammas = [0.8]
    flipZs = [0]
    Ls = numpy.arange(100, 500, 50)
    e0s = []
    z = ''
    from scipy.optimize import curve_fit
    def fitfun(L, einf, Lcoef):
        return einf + Lcoef / L
    xdata = []
    ydata = []
    for L in Ls:
        for gamma in gammas:
            for flipZ in flipZs:
                # t = tester(L, N, gamma, flipZ)
                # z += f'{t:10.8f} {exact_inf(N, gamma):10.8f}\n'
                x= inf_test(N, L, gamma, flipZ)
                xdata.append(L)
                ydata.append(x)

    params = curve_fit(fitfun, xdata, ydata)
    print(params)
    print(params[0][0])
    print(exact_inf(N, -gammas[0]))
    print(exact_inf(N, gammas[0]))


def tester(L, N, gamma, flipZ):
    m1 = make_m1(L, N, gamma, flipZ)
    m2 = make_m2(L, N, gamma, flipZ)
    sol1 = scipy.linalg.eig(m1)
    sol2 = scipy.linalg.eig(m2)
    eigs1 = [e ** (1./N) * gamma / L for e in sol1[0]]
    eigs2 = sol2[0]
    eigs1.sort()
    eigs2 = sorted(eigs2, key=lambda x: -abs(x.imag))
    eigs2 = [e / L for e in eigs2]
    sign = ''
    if flipZ: sign = '-'
    print(f'L={L}, N={N}, lambda={sign}{gamma}')
    e0 = 0
    for i in range(len(sol1[0])):
        e0 -= eigs1[i]
    print(f'Ground state = {e0:10.7f}')
    print(f'Ground state inf = {exact_inf(N, gamma)}')
    print(f'{"Real":>10} {"Imag":>10}')
    for i in range(len(sol1[0])):
        e1 = eigs1[i]
        # e2 = eigs2[i]
        e2, d = find_closest(e1, eigs2)
        w = 10
        print(f'{e1.real:10.7f} {e1.imag:10.7f} {e2.real:10.7f} {e2.imag:10.7f}')
        # print(f'{e1.real:10.7f} {e1.imag:10.7f}')
    print('\n')
    # for i in range(len(sol2[0])):
    #     e1 = eigs2[i]
    #     print(e1)

    for i in range(len(sol1[0])):
        continue
        e = sol1[0][i]
        v = sol1[1][:,i]
        r, n, n2 = eigenvalue_test(m1, v)
        n = abs(n)
        n2 = abs(n2)
        n0 = numpy.inner(numpy.conjugate(v), v)
        print(n0, e)
        print(abs(e/r), n, n2)
    return e0

Ls = [2, 3, 4, 5, 6, 7, 8, 9, 10]
Ls = [12, 40, 80]
N = 3
gammas = [0.5, 1.0, 1.5]
gammas = [0.8]
flipZs = [0]
Ls = numpy.arange(100, 500, 50)
Ls = [4, 5, 10, 11]
flipZs = [1]
e0s = []
z = ''
for L in Ls:
    for gamma in gammas:
        for flipZ in flipZs:
            t = tester(L, N, gamma, flipZ)
