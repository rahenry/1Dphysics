#!/usr/bin/env python3

import math
import numpy
import scipy
import random
from decimal import *
from scipy.special import hyp2f1
from scipy.optimize import root
from scipy.optimize import brentq
from scipy.optimize import fsolve

def fix_pfs(pfs, N):
    pfs_new = []
    for p in pfs:
        best = -100 + 1.j
        for i in range(N):
            ptest = p * numpy.exp(i * 2.j * numpy.pi / N)
            if ptest.real > best.real:
                best = ptest
        pfs_new.append(best)
    return pfs_new

# def fix_pfs(pfs, N):
#     return pfs
#     pfs_new = []
#     a0 = 2. * numpy.pi / N
#     for p in pfs:
#         a = numpy.angle(p) / a0
#         f = numpy.floor_divide(a, a0)
#         # print(p, a, f)
#         anew = (a - f) * a0
#         # print(anew)
#         pnew = abs(p) * numpy.exp(1.j * anew)
#         pfs_new.append(pnew)
#     return pfs_new


def gap_inf(s):
    N = s['N']
    lamb = s['lambda']
    omega = numpy.exp(2.j * numpy.pi / N)
    phi = s['phi']
    g = lamb * numpy.exp(phi * 2.j * numpy.pi / N)
    u = 1./ s['lambda_scale']
    u *= numpy.exp(-phi * 1.j * numpy.pi / N)
    res = (1. - g ** (N / 2.)) ** (2. / N)
    best = 10000
    return res * (1.-omega) * u

    for i in range(N-1):
        z = (1. - (omega ** (i+1))) * res * u
        # if abs() < best.real:
        if z.real < best.real:
            best = z
    return best


def exact_inf(N, lamb):
    return -hyp2f1(-1./N, -1./N, 1, lamb**(N))

def einf(s):
    N = s['N']
    lamb = s['lambda']
    # C = -1
    # q = s['lambda'] * numpy.exp(C*numpy.pi * 2.j * s['phi']/N)
    q = s['lambda_full']

    if lamb < 1:
        v = exact_inf(N, q)
        return v
    else:
        v = exact_inf(N, 1./q) * q
        return v

def rhoX(N, z):
    if z< 1:
        return hyp2f1(-1./N, (N-1.)/N, 1, z**N)
    else:
        return 1./N*(z**(1.-N))*hyp2f1((N-1.)/N, (N-1.)/N, 2, z**(-N))

def epsilon(k, N, gamma):
    return (1. +gamma**N + 2.*gamma**(0.5*N) * numpy.cos(k))**(1./N)

def make_matrix_data(L, N, gamma, flipZ=False):
    rows = [0]
    cols = [0]
    data = [1]
    # gamma = (1.0-gamma)/gamma
    # if flipZ:
    #     g = -g
    # if flipZ:
    #     gamma = -gamma
    if gamma < 0:
        gamma = numpy.complex64(gamma)
    g = gamma ** (-0.5*N)
    # g = numpy.pow(gamma, 0.5*N)
    z = 1 + g*g
    g2 = g
    g1 = g
    # z = 1 + g*g

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
    if gamma < 0:
        data = [-x for x in data]
    return rows, cols, data

def exact_parafermions_matrix(L, N, gamma, flipZ=False):
    rows, cols, data = make_matrix_data(L, N, gamma, flipZ)
    m = numpy.zeros((L, L), dtype=complex)
    for (i,j,d) in zip(rows, cols, data):
        m[i, j] = d
    sol = scipy.linalg.eig(m)
    e = 0
    pfs = []
    for x in sol[0]:
        p = x ** (1./N)
        p = p/L*gamma
        if gamma < 0:
            p = -p
        pfs.append(p)
    pfs.sort()
    e = -numpy.sum(pfs)
    # e = e / L
    return (e, sol[0], pfs)

def some_parafermions(L, N, gamma, n):
    rows, cols, data = make_matrix_data(L, N, gamma)
    m = scipy.sparse.coo_matrix((data, (rows, cols)), (L, L))
    sol = scipy.sparse.linalg.eigs(m, k=n, which='LM')
    return sol


def make_pf_spectrum(g, pfs, N, nstates=None):
    pfs = sorted(pfs)
    n_pfs = 1
    all_energies = [g]
    while(True):
        energies = [g]
        for i in range(n_pfs):
            pf = pfs[i]
            energies_new = []
            for e in energies:
                for i in range(1, N):
                    enew = e + pf - numpy.exp(2.j*numpy.pi*i/N) * pf
                    energies_new.append(enew)
            energies = energies_new
        all_energies += energies
        n_pfs += 1
        if n_pfs > len(pfs):
            break
        if nstates and len(all_energies) > nstates:
            break
    print(len(all_energies), N ** len(pfs))
    exit()
    return all_energies

def make_coords(L, n):
    res = [[x] for x in range(L)]
    for i in range(n-1):
        res_new = []
        for x in res:
            for j in range(L):
                if j in x:
                    break
                res_new.append([*x, j])
        res = res_new
    return res


def make_pf_spectrum2(g, pfs, N, depth=None):
    pfs = sorted(pfs)
    n_pfs = 0
    all_energies = [g]
    if depth is None:
        depth = len(pfs)+1
    while(True):
        n_pfs += 1
        # if n_pfs > depth and not n_pfs > len(pfs) - depth:
        if n_pfs > len(pfs):
            break
        if n_pfs > depth:
            continue

        # print(f'{n_pfs} / {depth}')
        coords = make_coords(len(pfs), n_pfs)
        for c in coords:
            energies = [g]
            for x in c:
                pf = pfs[x]
                energies_new = []
                for e in energies:
                    for i in range(1, N):
                        enew = e + pf - numpy.exp(2.j*numpy.pi*i/N) * pf
                        energies_new.append(enew)
                energies = energies_new
            all_energies += energies
            # print(n_pfs, len(coords), len(energies), len(all_energies))



    all_energies.sort()
    return all_energies

def make_pf_spectrum3(g, pfs, N, depth=None):
    pfs = sorted(pfs)
    n_pfs = 0
    all_energies = [(g,0)]
    if depth is None:
        depth = len(pfs)
    while(True):
        n_pfs += 1
        # if n_pfs > depth and not n_pfs > len(pfs) - depth:
        if n_pfs > len(pfs):
            break
        if n_pfs > depth:
            continue

        # print(f'{n_pfs} / {depth}')
        coords = make_coords(len(pfs), n_pfs)
        for c in coords:
            energies = [(g,0)]
            for x in c:
                pf = pfs[x]
                energies_new = []
                for e, q in energies:
                    for i in range(1, N):
                        enew = e + pf - numpy.exp(2.j*numpy.pi*i/N) * pf
                        qnew = (q + i) % N
                        energies_new.append((enew, qnew))
                energies = energies_new
            all_energies += energies
            # print(n_pfs, len(coords), len(energies), len(all_energies))


    return all_energies

def make_pf_spectrum_rand(g, pfs, N, depth=None):
    pfs = sorted(pfs)
    n_pfs = 0
    all_energies = [(g,0)]
    if depth is None:
        depth = len(pfs)
    while(True):
        n_pfs += 1
        # if n_pfs > depth and not n_pfs > len(pfs) - depth:
        if n_pfs > len(pfs):
            break
        if n_pfs > depth:
            continue

        # print(f'{n_pfs} / {depth}')
        coords = make_coords(len(pfs), n_pfs)
        for c in coords:
            energies = [(g,0)]
            for x in c:
                pf = pfs[x]
                energies_new = []
                for e, q in energies:
                    i = random.randrange(1, N)
                    enew = e + pf - numpy.exp(2.j*numpy.pi*i/N) * pf
                    qnew = (q + i) % N
                    energies_new.append((enew, qnew))
                energies = energies_new
            all_energies += energies
            # print(n_pfs, len(coords), len(energies), len(all_energies))


    return all_energies

def make_m1(s):
    L = s["L"]
    N = s["N"]
    lamb = s["lambda"]
    phi = s["phi"]
    g = s['lambda_full']
    g = g ** (0.5*N)
    g1 = g
    g2 = g
    z = 1 + g1*g2
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
        m[i, j] = d
    return m



def pfs_matrix(s):
    L = s["L"]
    N = s["N"]
    lamb = s["lambda"]
    phi = s["phi"]
    lf = s["lambda_full"]
    m1 = make_m1(s)
    sol = scipy.linalg.eig(m1)
    g = lamb * numpy.exp(phi * 2.j * numpy.pi / N)
    # eigs = sol[0]
    # eigs = [e ** (1./N) * g for e in sol[0]]
    eigs = [e ** (1./N) for e in sol[0]]
    if s['fix_pfs']:
        eigs = fix_pfs(eigs, s['N'])
    # eigs.sort(key = lambda x : abs(x))
    eigs.sort(key = lambda x : x.real)
    e0 = 0
    for p in eigs:
        e0 -= p
    return eigs, e0

def k_inf(j, L, N, lamb):
    return numpy.pi * j / L - numpy.pi * j / L / L / (1. + lamb**(0.5*N))

def k_to_eps(k, lamb, N):
    res = 1. + lamb**(N) + 2.*(lamb**(N/2.))*numpy.cos(k)
    return res ** (1./N)

def eps_to_k(e, g, N):
    res = (e ** N) - 1. - (g ** N)
    res /= 2. * (g ** (N/2.))
    # res = -res
    # res = abs(res)
    r1 = res
    res = numpy.arccos(res)
    # print(e, g, N, r1, res)
    return res

def kinf2(j, L, N, lamb):
    res = g ** (N/2.) + numpy.cos(numpy.pi * j / L)
    res = res / numpy.sin(numpy.pi * j / L)
    res = numpy.arctan(-1./res)
    res = numpy.pi * j - res
    res = res / L
    return res
Ls = [1000]
Ns = [3]
lambdas = [0.2]
phis = [0.2]
def form(a):
    # return '%0.4f+%0.4fi' % (a.real, a.imag)
    return "{num.real:+0.04f}{num.imag:=+0.08f}j".format(num=a)
def formr(a):
    return "{num.real:0.04f}".format(num=a)
for L in Ls:
    for N in Ns:
        for phi in phis:
            for lamb in lambdas:
                1
                print(123)
                s = {}
                s['lambda'] = lamb
                s['phi'] = phi
                s['L'] = L
                s['N'] = N
                g = lamb * numpy.exp(2.j * numpy.pi * s['phi'] / 2. / N)
                s['lambda_full'] = g
                s['fix_pfs'] = 0
                eigs, e0 = pfs_matrix(s)
                ei = einf(s) * L
                print(e0, ei)
                kinfs = []
                print('...')
                j = 0
                for e in eigs:
                    # ki = k_inf(j+1, L, N, g)
                    ki = kinf2(j+1, L, N, g)
                    kinfs.append(ki)
                    j+= 1
                kinfs.sort(key=lambda x: -x)


                alpha = -(g ** (N/2.))
                def f1(x):
                    res = numpy.sin((L+1)*x)/numpy.sin(L*x) -alpha
                    res = abs(res)
                    return res

                def f2(x, j, L):
                    k = x
                    v = x - numpy.pi * (L-j+1) / L
                    v = v*L
                    res = numpy.sin(v) * numpy.cos(k) + numpy.cos(v) * numpy.sin(k)
                    res = res / numpy.sin(v) - alpha
                    res = abs(res)
                    return res

                j = 0
                for e in eigs:
                    ki = kinfs[j]
                    k = eps_to_k(e, g, N)
                    pi = k_to_eps(ki, lamb, N)

                    p2 = k_to_eps(k, lamb, N)
                    # print(e, k, p2, pi)
                    a1 = formr(f2(k, j, L))
                    a2 = formr(f2(ki, j, L))
                    # b1 = formr(f2(k, j, 100000))
                    # b2 = formr(f2(ki, j, 100000))
                    d = k - ki
                    a0 = numpy.pi * (L-j) / L
                    print(formr(a0), form(k), form(ki), form(d), formr(f1(k)), formr(f1(ki)), a1, a2)
                    # print(k, ki)
                    # print(f1(k), f1(ki))
                    # print('..')
                    j += 1


exit()
Ls = [1000]
Ns = [3]
lambdas = [0.2]
phis = [0.2]
for L in Ls:
    break
    for N in Ns:
        for phi in phis:
            for lamb in lambdas:
                s = {}
                s['lambda'] = lamb
                s['phi'] = phi
                s['L'] = L
                s['N'] = N
                g = lamb * numpy.exp(2.j * numpy.pi * s['phi'] / 2. / N)
                s['lambda_full'] = g
                s['fix_pfs'] = 0
                alpha = -(g ** (N/2.))
                for j in range(L):
                    ki = k_inf(j+1, L, N, g)
                    a0 = numpy.pi * (j+1) / L
                    c = numpy.arccos(alpha/2)
                    print(formr(a0), form(ki), c)


Ls = numpy.arange(100, 1000, 50)
Ns = [3]
lambdas = [0.2]
phis = [0.3]
for L in Ls:
    for N in Ns:
        for phi in phis:
            for lamb in lambdas:
                s = {}
                s['lambda'] = lamb
                s['phi'] = phi
                s['L'] = L
                s['N'] = N
                g = lamb * numpy.exp(2.j * numpy.pi * s['phi'] / 2. / N)
                s['lambda_full'] = g
                s['fix_pfs'] = 0
                alpha = -(g ** (N/2.))

                eigs, e0 = pfs_matrix(s)
                ei = einf(s) * L
                eigs = sorted(eigs, key=lambda x: x.real)
                kinfs = []
                j = 0
                for e in eigs:
                    ki = k_inf(j+1, L, N, g)
                    kinfs.append(ki)
                    j+= 1
                kinfs.sort(key=lambda x: x.real)


                alpha = -(g ** (N/2.))

                j = L // 2
                e = eigs[j]
                ki = kinfs[j]
                k = eps_to_k(e, g, N)
                pi = k_to_eps(ki, lamb, N)

                p2 = k_to_eps(k, lamb, N)
                a0 = numpy.pi * (L-j) / L
                j += 1
                print(formr(a0), form(k), form(ki))
