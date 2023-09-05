#!/usr/bin/env python3


from System import *
from SystemSet import *
from args import Q
from exact_fp import *
from rah_utility import *
from output_figs import *

import numpy, scipy

def get_pfs(N, L, gamma, phi):
    f = os.path.join('pfdata', f'{N}_{L}_{phi}_{gamma}')
    pfs = None
    if os.path.exists(f) and not Q.rerun_measurements:
        with open(f, 'rb') as f1:
            pfs = dill.load(f1)
    else:
        lamb = gamma / (1.-gamma)
        lsc = 1. / (1.-gamma)
        s = {
            'L' : L,
            'gamma' : gamma,
            'phi':  phi,
            'N' : N,
            'lambda' : lamb,
            'lambda_scale' : lsc,
            'fix_pfs' : 1,
            }
        pfs, e0 = pfs_matrix(s)
        pfs = numpy.array(pfs)
        with open(f, 'wb+') as f1:
            dill.dump(pfs, f1)
            # print("saving")
    return pfs

def eps_k(k, lamb, N):
    res = 1. + lamb**N + 2.*(lamb**(N/2.))*numpy.cos(k)
    return res ** (1./N)

def add_roots(f, x0, x1, res, delta=1E-10):
    new_root = brentq(f, x0, x1)
    return res + [new_root] + add_roots(f, x0, new_root-delta) + add_roots(f, new_root, x1)
def find_kjs(gamma, L, N):
    res = []
    def f(k):
        return kj_roots(k, g, L, N)

    while (True):
        scipy.optimize.root
        new_root = brentq(f, minimum, maximum)
        1


    return res

def eps_to_k(e, g, N):
    # print(e, g, N)
    res = (e ** N) - 1. - (g ** N)
    res /= 2. * (g ** (N/2.))
    # res = -res
    # res = abs(res)
    res = numpy.arccos(res)
    return res

def get_pfs_lamb(N, L, lamb, phi):
    d1 = dict(L=L, N=N, phi=phi)
    d1['lambda'] = lamb
    m1 = make_m1(d1)
    sol = scipy.linalg.eig(m1)
    g = lamb * numpy.exp(phi * 2.j * numpy.pi / N)
    eigs = [e ** (1./N) for e in sol[0]]
    eigs = [e*g for e in eigs]
    eigs = fix_pfs(eigs, N)
    eigs.sort()
    return eigs

def f_kj(k, lamb, L, N):
    return numpy.sin((L+1.) * k) + (lamb ** (-N/2.)) * numpy.sin(L*k)

def f_kj1(k, lamb, L, N):
    return (L+1) * numpy.cos((L+1.) * k) + (lamb ** (-N/2.)) * L * numpy.cos(L*k)

def f_kj2(k, lamb, L, N):
    return - (L+1.) * (L+1.) * numpy.sin((L+1.) * k) - (lamb ** (-N/2.)) * L*L *  numpy.sin(L*k)

Ns = [3]
lambs = [0.1, 0.9]
lambs = [2.0]
phis = [0.1]
# m = [0.41421356, 0.34740645]
Ls = [3]
m = [0.41421356, 0.65259355]
v= 1.095893964469992
u = -0.3715713274619571j

# Ls = [4]
# m = [0.42715452,  0.5       ]
# u = -0.32666823811850276j

g1 = m[0]
lambs = [g1 / (1.-g1)]
phis = [m[1]]
for L in Ls:
    for N in Ns:
        for lamb in lambs:
            for phi in phis:
                print(L, N, lamb, phi)
                pfs = get_pfs_lamb(N, L, lamb, phi)
                g = lamb * numpy.exp(phi * 2.j * numpy.pi / N)

                e0 = 0
                pfs.sort()
                for p in pfs:
                    # p = g * p
                    e0 -= p
                    k = eps_to_k(p, g, N)
                    d = f_kj(k, g, L, N)
                    d = abs(d)
                    # print(p, k, eps_k(k, g, N), d)
                    print(p, k, d)

                def f1(k):
                    return f_kj(k+u, g, L, N)
                def f2(k):
                    return f_kj1(k+u, g, L, N)
                def f3(k):
                    return f_kj2(k+u, g, L, N)
                # z = scipy.optimize.root(f1, 1.1)
                xdata = numpy.linspace(0, numpy.pi, 100)
                xdata = numpy.linspace(v - 0.1 * v, v + 0.1*v, 100)
                ydata = f1(xdata)
                p1 = plt.scatter(xdata, numpy.real(ydata))
                p2 = plt.scatter(xdata, numpy.imag(ydata))
                ydata = f2(xdata)
                p3 = plt.scatter(xdata, numpy.real(ydata), marker='x')
                p4 = plt.scatter(xdata, numpy.imag(ydata), marker='x')
                # ydata = f3(xdata)
                # p5 = plt.scatter(xdata, numpy.real(ydata), marker='v')
                # p6 = plt.scatter(xdata, numpy.imag(ydata), marker='v')
                # plt.legend((p1, p2, p3, p4, p5, p6),
                #            (
                #            )
                plt.show()

                # print(z.x, f1(z.x), eps_k(z.x, lamb, N))



                gamma = lamb / (1.+lamb)
                lsc = 1. / (1.-gamma)
                d1 = dict(L=L, N=N, phi=phi, gamma=gamma, lambda_scale=lsc, fix_pfs=0)
                d1['lambda'] = lamb
                pfs2, e02 = pfs_matrix(d1)
                pfs2.sort()

                q = g**N
                d1['lambdas'] = [q, 1]
                # d1['lambdas'] = [1, q]
                d1['M'] = 2 * L - 1
                d1['p'] = 1
                pfs_poly = sympy_eps2(d1)
                pfs_poly = [p * L for p in pfs_poly]
                print(pfs_poly)
                # print(e0)
                # print(e02*lsc*L)
