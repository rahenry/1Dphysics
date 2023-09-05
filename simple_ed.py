#!/usr/bin/env python3
import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, os, shutil
from rah_numpy import eigenvalue_test
from rah_utility import rel_diff

def translate(u, L, N):
    t = 0
    for j in range(L):
        s = sigma(N, u, j)
        t = index(N, t, (j+1) % L, s)
    return t

def get_period(z):
    for i in range(1, len(z)+1):
        zp = z[i:] + z[:i]
        if zp == z:
            return i

def stringify(u, L, N):
    s = ''
    for j in range(L):
        s += str(sigma(N ,u, j))
    return s

def charge(u, L, N):
    q = 0
    for j in range(L):
        q += sigma(N, u, j)
    return q % N

def testt(s):
    N = s['N']
    L = s['L']
    n = N ** L
    basis = {}
    for q in range(N):
        for i in range(L):
            basis[(q, i)] = []
    for u in range(n):
        s = stringify(u, L, N)
        t = translate(u, L, N)
        p = get_period(s)
        Q = charge(u, L, N)
        for j in range(p):
            q = u
            v = numpy.full((n), 0.j)
            k = j * (L // p)
            for i in range(p):
                v[q] = 1 * numpy.exp(2.j * numpy.pi / N * i * k)
                q = translate(q, L, N)
            basis[(Q, k)].append(v)

    for p, vs in basis.items():
        print(p, len(vs))
        for v in vs:
            print(stringify(v, L, N))
    return basis

def make_extra_matrices(s):
    N = s['N']
    L = s['L']
    n = N ** L
    cZ = s.model.cZ
    cX = s.model.cX

    T = numpy.full((n, n), 0.j)
    for u in range(n):
        t = 0
        for j in range(L):
            s = sigma(N, u, j)
            t = index(N, t, (j+1) % L, s)
        T[u][t] = 1

    Q = numpy.full((n, n), 0.j)
    for u in range(n):
        t = 0
        for j in range(L):
            s = sigma(N, u, j)
            t = index(N, u, j, (s+1)%N)
        Q[u][t] = 1

    return T, Q

def get_sectors(s):
    T, Q = make_extra_matrices(s)
    sol = scipy.linalg.eig(T, Q)
    sectors = {}
    for i, e in enumerate(sol[0]):
        found = None
        for s in sectors:
            if rel_diff(e, s) < 1E-10:
                found = s
                sectors[found].append(sol[1][:,1])
        if not found:
            sectors[e] = [sol[1][:, i]]
    # for x, y in sectors.items():
    #     print(x, len(y))


def solve_simple_ed(s):
    N = s['N']
    # print(s)
    L = s['L']
    n = N ** L
    # testt(s)
    # return
    # get_sectors(s)
    # exit()
    # T, Q = make_extra_matrices(s)

    # sol = scipy.linalg.eig(T, Q)
    # sol = scipy.linalg.eig(T)
    # print(sol[0])
    # for i in range(n):
    #     v = sol[1][:,i]
    #     a, b, c = eigenvalue_test(T, v)
    #     print(a, b, c)

    # exit()
    m = make_H(s)
    # print(cZ, cX)
    sol = scipy.linalg.eig(m)
    es = sol[0]
    # es = sorted(es, key=lambda x: abs(x))
    inds = numpy.argsort(sol[0])
    s['energies'] = sol[0][inds]
    s.eigenvectors = sol[1][:,inds]
    s['e_0'] = sol[0][0]

    count = 0
    i = 0
    # for i in inds:
    #     # e = s['energies'][i]
    #     e = sol[0][i]
    #     v = sol[1][:,i]
    #     # v = s.eigenvectors[i]
    #     if abs(e) < 1E-6:
    #         print('...')
    #         count += 1
    #         c = dict(s.config_base)
    #         c['gamma'] = 0
    #         h1 = make_H(c, 0, 1)
    #         h2 = make_H(c, 1, 0)
    #         a = numpy.dot(v.conj(), v)
    #         b = numpy.dot((h1@v).conj(), h1@v)
    #         c = numpy.dot(v.conj(), h2@v)
    #         print(abs(e))
    #         print(abs(a), abs(b), abs(c))
    #         print(v)
    #         print(h1@v)
    #         print('_')
    #         continue
    #         for j in range(len(v)):
    #             if abs(v[j]) > 1E-5:
    #                 print(j, v[j])
    #         print(v)
    #         print(m@v)
    #         print(h1@v)

    #     i += 1
    # exit()
    # print(sol[0][inds])

def sigma(N, state, site):
    j = site
    i = state
    N = N
    return (i // (N ** j)) % N

def index(N, initial, j_modified, sigma_j_f):
    sigma_j_i = sigma(N, initial, j_modified)
    return initial + ((N ** j_modified) * (sigma_j_f - sigma_j_i))

def make_H(s, cX=None, cZ=None):
    N = s['N']
    L = s['L']
    n = N ** L
    m = numpy.full((n, n), 0.j)
    if cZ is None:
        cZ = s.model.cZ
    if cX is None:
        cX = s.model.cX

    for i in range(L-1):
        for u in range(n):
            delta = sigma(N, u, i) - sigma(N, u, i+1)
            m[u][u] += cZ * numpy.exp(2.j * numpy.pi / N * delta)

    for i in range(L):
        for u in range(n):
            uP = index(N, u, i, (sigma(N, u, i) + 1) % N)
            m[u][uP] += cX
    if s['bc'] == 1:
        for u in range(n):
            delta = sigma(N, u, L-1) - sigma(N, u, 0)
            m[u][u] += cZ * numpy.exp(2.j * numpy.pi / N * delta)

    # for u in range(n):
    #     print(m[u][u])


    return m
