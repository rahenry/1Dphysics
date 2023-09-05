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

# This code fixes a set of quasienergies contained in the list "pfs" so they
# all have complex arguments in the interval (-pi/N, pi/N]. This is
# achieved by trying all the equivalent rotations of each pf, found by
# multiplying it by a power of omega = exp(2pi*i/N), and choosing the one
# with the largest real part (choosing the largest real part is equivalent
# to selecting an argument in the required interval, as it chooses the
# result furthest to the right in the complex plane).
def fix_pfs(pfs, N):
    pfs_new = []
    for p in pfs:
        best = -1E50 + 1.j
        for i in range(N):
            ptest = p * numpy.exp(i * 2.j * numpy.pi / N)
            if ptest.real > best.real:
                best = ptest
        pfs_new.append(best)
    return pfs_new

# This code fixes a set of quasienergies contained in the list "pfs" so they
# all have complex arguments in the interval (-pi/N, pi/N]. This is
# achieved by trying all the equivalent rotations of each pf, found by
# multiplying it by a power of omega = exp(2pi*i/N), and choosing the one
# with the largest real part (choosing the largest real part is equivalent
# to selecting an argument in the required interval, as it chooses the
# result furthest to the right in the complex plane).
def fix_pfs(pfs, N):
    pfs_new = []
    for p in pfs:
        best = p
        # This could be range(1, N) since we started with p * 1, but it
        # doesn't really matter
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
    return -hyp2f1(-1./N, -1./N, 1, lamb**(-N))
def f12(N, z):
    return -hyp2f1(-1./N, -1./N, 1., z)

def einf2(s):
    N = s['N']
    lamb = s['lambda']
    phi = s['phi']
    # C = -1
    # q = s['lambda'] * numpy.exp(C*numpy.pi * 2.j * s['phi']/N)
    q = s['lambda_full']
    if lamb <= 1:
        v = -hyp2f1(-1./N, -1./N, 1, q**N)
        return v
    else:
        # v = -hyp2f1(-1./N, -1./N, 1, q**(-N)) * abs(q)
        g = q * numpy.exp(2.j * numpy.pi / N * -phi)
        v = -hyp2f1(-1./N, -1./N, 1, q**(-N)) * q
        return v

def einf(s):
    N = s['N']
    lamb = s['lambda']
    # C = -1
    # q = s['lambda'] * numpy.exp(C*numpy.pi * 2.j * s['phi']/N)
    q = s['lambda_full']

    if lamb < 1:
        # v = f12(N,  c
        v = exact_inf(N, q) * q * q
        a = -1.j * numpy.pi * s['phi'] / N
        # v = v * numpy.exp(a)
        return v
    else:
        v = exact_inf(N, 1./q) * q
        a = -1.j * numpy.pi * s['phi'] / N
        # v = v * numpy.exp(a)
        # v = numpy.conj(v)
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

if __name__ == "__main__":
    L0 = 100
    L1 = 700
    nL = 5
    N = 3
    g0 = 0.499
    g1 = 0.501
    ng = 30
    npfs = 50
    from rah_utility import mkdir
    mkdir("test_graphs")
    import os
    biggest = 0
    for L in numpy.linspace(L0, L1, nL):
        xdata = []
        ydata = []
        L = int(L)
        name = f'L={L}.png'
        name = os.path.join('test_graphs',name)
        for gamma in numpy.linspace(g0, g1, ng):
            e, x, pfs = exact_parafermions_matrix(L, N, gamma)
            maxpf = max(pfs)
            if maxpf > biggest:
                biggest = maxpf
            pfs = sorted(pfs)
            xdata.append(gamma)
            ydata.append(pfs[0:npfs])
            ydata[-1].append(pfs[-1])
        import matplotlib.pyplot as plt
        plt.clf()
        plt.plot(xdata, ydata)
        # plt.ylim([0, biggest])
        plt.savefig(name)



    exit()
    plt.clf()
    gamma = 0.5
    xdata = []
    ydata = []
    L0 = 100
    L1 = 1000
    dL = 100
    for L in numpy.arange(L0, L1, dL):
        e, x, pfs = exact_parafermions_matrix(L, N, gamma)
        pfs.sort()
        pfs = pfs[:npfs]
        xdata.append(L)
        ydata.append(pfs)
    plt.plot(xdata, ydata)
    plt.show()


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
    if s['swap_lambda']:
        g = g ** (0.5*N)
    else:
        g = g ** (-0.5*N)
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
    m1 = make_m1(s)
    sol = scipy.linalg.eig(m1)
    g = s['lambda_full']
    if s['swap_lambda'] == 1:
        eigs = [e ** (1./N) for e in sol[0]]
    else:
        eigs = [e ** (1./N) * g for e in sol[0]]
    u = 1
    if s.get('scale', 'lambda') == 'gamma':
        u = 1. / s['lambda_scale']
    v = s.get('angle_style', 1)
    mu = s.get('mu', 0)
    u *= numpy.exp(mu * phi * 2.j * numpy.pi / N)
    eigs = [e * u for e in eigs]
    e0 = 0
    if s['fix_pfs']:
        eigs = fix_pfs(eigs, s['N'])
    eigs.sort(key = lambda x : x.real)
    for p in eigs:
        e0 -= p
    a = max([abs(x) for x in eigs])
    s['pf_scale'] = a
    return eigs, e0

def spectrum_matrix(s):
    pfs, e0 = pfs_matrix(s)
    spect = make_pf_spectrum2(e0, pfs, s['N'])
    spect.sort()
    pfs.sort()
    return spect, pfs

def make_m1_sparse(s):
    L = s["L"]
    N = s["N"]
    lamb = s["lambda"]
    phi = s["phi"]
    g = lamb * numpy.exp(phi * 2.j * numpy.pi / N)
    g = g ** (-0.5*N)
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

def pfs_matrix_limited(s, M):
    L = s["L"]
    N = s["N"]
    lamb = s["lambda"]
    phi = s["phi"]
    g = lamb * numpy.exp(phi * 2.j * numpy.pi / N)

    rows, cols, data = make_matrix_data(L, N, g, 0)
    m = scipy.sparse.coo_matrix((data, (rows, cols)), (L, L))

    sol = scipy.sparse.linalg.eigs(m, k=M, which='LM')

    eigs = [e ** (1./N) * g / L for e in sol[0]]
    u = 1./ s['lambda_scale']
    u *= numpy.exp(-phi * 1.j * numpy.pi / N)
    eigs = [e * u for e in eigs]
    e0 = 0
    eigs = fix_pfs(eigs, s['N'])
    eigs.sort()
    return eigs

def pfs_matrix_xyh(s):
    L = s["L"]
    N = s["N"]
    phi = s["phi"]
    m1 = make_m1_xyh(s)
    sol = scipy.linalg.eig(m1)
    g = s['gamma_full']
    eigs = [e ** (1./N) * 0.5 for e in sol[0]]
    e0 = 0
    if s['fix_pfs']:
        eigs = fix_pfs(eigs, s['N'])
    eigs.sort(key = lambda x : x.real)
    for p in eigs:
        e0 -= p
    a = max([abs(x) for x in eigs])
    s['pf_scale'] = a
    return eigs, e0

def make_m1_xyh(s):
    L = s["L"]
    N = s["N"]
    g = s['gamma_full']
    bc = s["bc"]
    a = 0.5
    b = g * 0.5

    rows = []
    cols = []
    data = []
    for i in range(L-1):
        rows.append(i+1)
        cols.append(i)
        data.append(a)
        rows.append(i)
        cols.append(i+1)
        data.append(a)
    if s["bc"] == 1:
        rows.append(L-1)
        cols.append(0)
        data.append(a)
        rows.append(0)
        cols.append(L-1)
        data.append(a)

    A = numpy.zeros((L, L), dtype=complex)
    for (i,j,d) in zip(rows, cols, data):
        A[i, j] = d

    rows = []
    cols = []
    data = []
    for i in range(L-1):
        rows.append(i+1)
        cols.append(i)
        data.append(b)
        rows.append(i)
        cols.append(i+1)
        data.append(-b)
    if s["bc"] == 1:
        rows.append(L-1)
        cols.append(0)
        data.append(-b)
        rows.append(0)
        cols.append(L-1)
        data.append(b)
    B = numpy.zeros((L, L), dtype=complex)
    for (i,j,d) in zip(rows, cols, data):
        B[i, j] = d


    res = (A-B)@(A+B)
    return res

def pfs_matrix_xyh2(s):
    L = s["L"]
    N = s["N"]
    phi = s["phi"]
    m1 = make_m1_xyh2(s)
    sol = scipy.linalg.eig(m1)
    g = s['gamma_full']
    eigs = [e ** (1./N) * 0.5 * 4. / (1. + s["gamma_full"]) for e in sol[0]]
    e0 = 0
    if s['fix_pfs']:
        eigs = fix_pfs(eigs, s['N'])
    eigs.sort(key = lambda x : x.real)
    for p in eigs:
        e0 -= p
    a = max([abs(x) for x in eigs])
    s['pf_scale'] = a
    return eigs, e0

def make_m1_xyh2(s):
    L = s["L"]
    N = s["N"]
    g = s['gamma_full']
    bc = s["bc"]
    a = 0.5
    b = g * 0.5

    rows = []
    cols = []
    data = []
    for i in range(L-1):
        rows.append(i+1)
        cols.append(i)
        data.append(a)
        rows.append(i)
        cols.append(i+1)
        data.append(a)
    if s["bc"] == 1:
        rows.append(L-1)
        cols.append(0)
        data.append(a)
        rows.append(0)
        cols.append(L-1)
        data.append(a)

    A = numpy.zeros((L, L), dtype=complex)
    for (i,j,d) in zip(rows, cols, data):
        A[i, j] = d

    rows = []
    cols = []
    data = []
    for i in range(L-1):
        rows.append(i+1)
        cols.append(i)
        data.append(b)
        rows.append(i)
        cols.append(i+1)
        data.append(-b)
    if s["bc"] == 1:
        rows.append(L-1)
        cols.append(0)
        data.append(-b)
        rows.append(0)
        cols.append(L-1)
        data.append(b)
    B = numpy.zeros((L, L), dtype=complex)
    for (i,j,d) in zip(rows, cols, data):
        B[i, j] = d


    res = (A-B)@(A+B)
    return res


def pfs_matrix_xyh2(s):
    L = s["L"]
    N = s["N"]
    phi = s["phi"]
    m1 = make_m1_xyh2(s)
    sol = scipy.linalg.eig(m1)
    g = s['gamma_full']
    # eigs = [e ** (1./N) / (1. + s["eta_full"]) for e in sol[0]]
    eigs = [e ** (1./N)  for e in sol[0]]
    e0 = 0
    if s['fix_pfs']:
        eigs = fix_pfs(eigs, s['N'])
    eigs.sort(key = lambda x : x.real)
    for p in eigs:
        e0 -= p
    a = max([abs(x) for x in eigs])
    s['pf_scale'] = a
    return eigs, e0

def make_m1_xyh2(s):
    L = s["L"]
    N = s["N"]
    g = s['eta_full']
    bc = s["bc"]
    a = g
    b = g

    rows = []
    cols = []
    data = []
    for i in range(L-1):
        rows.append(i+1)
        cols.append(i)
        data.append(a)
        rows.append(i)
        cols.append(i+1)
        data.append(1)
    if s["bc"] == 1:
        rows.append(L-1)
        cols.append(0)
        data.append(a)
        rows.append(0)
        cols.append(L-1)
        data.append(1)

    A = numpy.zeros((L, L), dtype=complex)
    for (i,j,d) in zip(rows, cols, data):
        A[i, j] = d

    rows = []
    cols = []
    data = []
    for i in range(L-1):
        rows.append(i+1)
        cols.append(i)
        data.append(1)
        rows.append(i)
        cols.append(i+1)
        data.append(b)
    if s["bc"] == 1:
        rows.append(L-1)
        cols.append(0)
        data.append(b)
        rows.append(0)
        cols.append(L-1)
        data.append(1)
    B = numpy.zeros((L, L), dtype=complex)
    for (i,j,d) in zip(rows, cols, data):
        B[i, j] = d


    # res = (A-B)@(A+B)
    res = A@B
    res = B@A
    return res
