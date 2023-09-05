#!/usr/bin/env python3

import math
import os
from exact_fp import *
import matplotlib.pyplot as plt
from rah_utility import rel_diff, mkdir, silent_remove
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt
import sympy

def find_eps_poly(s):
    M = s.model.M
    p = s.model.p
    N = s['N']
    lambdas = s.model.lambdas
    roots, poly = find_roots_poly(M, p, lambdas)
    pfs = [numpy.complex64(r)**(-1./s['N']) for r in roots]
    return pfs

def find_roots_poly(M, p, lambdas):
    P0 = Polynomial([1])
    polys = [P0]
    i = 0
    # Mbar = numpy.floor((M+p)/(p+1))
    while(True):
        q = lambdas[i % len(lambdas)]
        i += 1
        Pnew = polys[-1].copy()
        if len(polys) > p:
            Pnew -= q * Polynomial([0, 1]) * polys[-1-p]
        else:
            Pnew -= q * Polynomial([0, 1])
        polys.append(Pnew)
        if len(polys) == M+1:
            break

    roots = polys[-1].roots()
    return roots, polys[-1]

if __name__ == '__main__':
    f = os.path.join('graphs', 'fp_polynomial2')
    # silent_remove(f)
    ps = [1, 2, 3]
    ps = [1]
    Ms = numpy.arange(10, 30, 3)
    N = 3
    lambdas = [-0.01, 0.5, 0.5]
    for p in ps:
        continue
        for M in Ms:
            roots = find_roots_poly(M, p, lambdas)
            pfs = [numpy.complex64(r)**(-1./N) for r in roots]
            spectrum = make_pf_spectrum2(0, pfs, N)
            xdata = [x.real for x in spectrum]
            ydata = [x.imag for x in spectrum]
            plt.plot(xdata, ydata, linestyle='', marker='o', markersize=1.5)
            f = os.path.join('graphs', 'fp_polynomial2')
            mkdir(f)
            f = os.path.join(f, f'{p}_{M}.png')
            plt.savefig(f)
            plt.close()

    # Ms = numpy.arange(3, 35, 1)
    f = os.path.join('graphs', 'fp_polynomial2')
    # silent_remove(f)
    mkdir(f)
    Ms = numpy.arange(3, 35, 1)
    N = 3
    gammas = [0.2, 0.5, 0.8]
    p = 1
    Ms = [25]
    gammas = [0.1, 1.0, 5]
    gammas = [1, 3/2, 2]
    gammas = [2]
    ax_ind = 0
    n = 2*len(gammas)
    if 1 in gammas: n = n-1
    # fig, axs = plt.subplots(n, 1)
    for gamma in gammas:
        # continue
        # for reverse in [0, 1]:
        for reverse in [1]:
            if gamma == 1 and reverse == 1:
                continue
            lambdas = [gamma, 1]
            if reverse:
                lambdas.reverse()
            lambdas = [l**N for l in lambdas]
            print(f'lambdas = {lambdas}')
            xdata = []
            ydata = []
            xdata2 = []
            ydata2 = []
            fail_data = []
            for M in Ms:
                roots_old, poly = find_roots_poly(M, p, lambdas)
                u = sympy.Symbol('y')
                npoly = 0
                j = 0
                for x in poly:
                    npoly += x * (u**j)
                    j += 1
                npoly = sympy.Poly(npoly)
                roots = []
                try:
                    roots = npoly.nroots(maxsteps = 500, n=500)
                except:
                    roots = roots_old
                    print(1230123)
                    exit()

                roots_accepted = []
                roots_rejected = []
                for r in roots:
                    q = poly(r)
                    if abs(q) > 1E-12:
                        # print(gamma, M, r, q)
                        roots_rejected.append(r)
                    else:
                        roots_accepted.append(r)
                if M == 25:
                    print(lambdas, npoly)
                    for r in roots:
                        q = npoly(r)
                        print(q)
                        print(r)
                        print('.')
                        if abs(q) > 1E-12:
                            pass
                    exit()
                roots = roots_accepted
                # print(f'{len(roots_rejected)} roots rejected')
                fail_data.append(len(roots_rejected))

                pfs = [numpy.complex64(r)**(-1./N) for r in roots]
                for pf in pfs:
                    xdata.append(M)
                    ydata.append(abs(pf))


                M = int(M+0.01)
                if M % 2 == 1:
                    L = (M+1)//2
                    print(M, L)
                    # lamb = gamma / (1.-gamma)
                    lamb = gamma
                    if reverse:
                        lamb = 1./gamma
                    e, q, pfs = exact_parafermions_matrix(L, N, lamb)
                    # sc = 1./(1.-gamma)
                    sc = 1.
                    if reverse:
                        sc = gamma
                    xdata2 += [M for pf in pfs]
                    ydata2 += [abs(pf)*L*sc for pf in pfs]

            l = gamma
            if reverse: l = 1./gamma
            plt.close()
            plt.figure(figsize=(10, 10))
            fig, ax = plt.subplots()
            ax.set_title(f'lambda={l}')
            ax_ind += 1
            ax.plot(xdata, ydata, linestyle='', marker='_', markersize=10.0, color='red')
            ax.plot(xdata2, ydata2, linestyle='', marker='o', markersize=2.5, color='blue')
            ax.set_xlabel("M\n[no. failed roots]")
            ax.set_ylim(bottom=0)
            xtick_coords = Ms
            labs = []
            for i in range(len(Ms)):
                lab = f'{Ms[i]}\n[{fail_data[i]}]'
                labs.append(lab)

            ax.set_xticks(xtick_coords, labs)
            f = os.path.join('graphs', 'fp_polynomial2')
            f = os.path.join(f, f'lambda={l}.png')
            plt.tight_layout()
            plt.savefig(f)
    # plt.tight_layout()
    # plt.show()

    etas = [0.5, 1, 2]
    kappas = [0.5, 1, 2]
    coefs = [0.1, 0.5, 1, 2, 10]
    Ms = numpy.arange(3, 18, 1)
    for eta in coefs:
        for kappa in coefs:
            lambdas = [1, eta, kappa]
            lambdas = [l**N for l in lambdas]
            xdata = []
            ydata = []
            xdata2 = []
            ydata2 = []
            fail_data = []
            for M in Ms:
                print(M)
                roots_old, poly = find_roots_poly(M, p, lambdas)
                u = sympy.Symbol('y')
                npoly = 0
                j = 0
                for x in poly:
                    npoly += x * (u**j)
                    j += 1
                npoly = sympy.Poly(npoly)
                roots = []
                try:
                    roots = npoly.nroots(maxsteps = 500, n=200)
                except:
                    roots = roots_old

                roots_accepted = []
                roots_rejected = []
                for r in roots:
                    q = poly(r)
                    if abs(q) > 1E-12:
                        # print(gamma, M, r, q)
                        roots_rejected.append(r)
                    else:
                        roots_accepted.append(r)
                roots = roots_accepted
                print(f'{len(roots_rejected)} roots rejected')
                fail_data.append(len(roots_rejected))

                pfs = [numpy.complex64(r)**(-1./N) for r in roots]
                for pf in pfs:
                    xdata.append(M)
                    ydata.append(abs(pf))


            l = gamma
            plt.close()
            plt.figure(figsize=(10, 10))
            fig, ax = plt.subplots()
            T = f'cX=1, cY={eta}, cW={kappa}'
            ax.set_title(T)
            ax_ind += 1
            ax.plot(xdata, ydata, linestyle='', marker='o', markersize=2.0, color='blue')
            ax.set_xlabel("M\n[no. failed roots]")
            ax.set_ylim(bottom=0)
            xtick_coords = Ms
            labs = []
            for i in range(len(Ms)):
                lab = f'{Ms[i]}\n[{fail_data[i]}]'
                labs.append(lab)

            ax.set_xticks(xtick_coords, labs)
            f = os.path.join('graphs', 'fp_polynomial2')
            f = os.path.join(f, f'{T.replace(" ","")}.png')
            plt.tight_layout()
            plt.savefig(f)
    # plt.tight_layout()
    # plt.show()
