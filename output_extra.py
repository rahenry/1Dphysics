#!/usr/bin/env python3

from lmfit import minimize, Parameters
import matplotlib.pyplot as plt
import os
from Slicer import SystemSlicer
from pf_tester import pf_test
from tenpy.linalg.np_conserved import tensordot, inner, norm
from rah_utility import mkdir, matplotlib_set_aspect, extrema
import graph_settings
import matplotlib
import numpy
from figures import *
from args import Q
from data_set import DataSet, plot_multiple
import cProfile, pstats, io
from pstats import SortKey
from exact_fp import exact_parafermions_matrix, make_pf_spectrum2, rhoX, make_coords, fix_pfs, einf
from fp_more_matrices import *
from fp_polynomial import *
from matplotlib.offsetbox import AnchoredText
import matplotlib.gridspec as gridspec
from polytest import sympy_eps2
from point_generator import PointGenerator
from rah_mpl import draw_label
from labels import *
import scipy

def oe_1(S):
    ex = ['method']
    for x, y in S.slicer.get_sets(exclude=ex).items():
        print(x)
        plt.clf()
        for s in y:
            if s['method'] == 'dmrg':
                pltargs = dict(linestyle='', color='red', marker='x', zorder=10, markersize=10)
                e0 = s['e_0']
                xdata = [e0.real]
                ydata = [e0.imag]
                plt.plot(xdata, ydata, **pltargs)
            else:
                pltargs = dict(linestyle='', color='blue', marker='o', zorder=5, markersize=3)
                energies = make_pf_spectrum2(s['e_0'], s['pfs'], s['N'])
                xdata = [x.real for x in energies]
                ydata = [x.imag for x in energies]
                plt.plot(xdata, ydata, **pltargs)

                pltargs = dict(linestyle='', color='green', marker='+', zorder=11, markersize=15, linewidth=2)
                e0 = s['e_0']
                print(e0)
                xdata = [e0.real]
                ydata = [e0.imag]
                plt.plot(xdata, ydata, **pltargs)

                pltargs = dict(linestyle='', color='orange', marker='+', zorder=11, markersize=15, linewidth=2)
                e0 = s['e_inf']
                print(e0)
                xdata = [e0.real]
                ydata = [e0.imag]
                plt.plot(xdata, ydata, **pltargs)
        f = f'{x}.png'
        f = os.path.join(S.graph_dir, f)
        plt.savefig(f)

def residual(p, x, data):
    e0 = p['e0_real'] + p['e0_imag'] * 1.j
    A = p['A_real'] + p['A_imag'] * 1.j
    # B = p['B_real'] + p['B_imag'] * 1.j
    # model = e0 + A / x + B / x / x
    model = e0 + A / x
    z = data - model
    z = z.view('float')
    return z

def lmtest(data):
    params = Parameters()
    params.add('e0_real', value=-1)
    params.add('e0_imag', value=0)
    params.add('A_real', value=1)
    params.add('A_imag', value=0)
    # params.add('B_real', value=-1)
    # params.add('B_imag', value=0)
    xdata = data['xdata']
    ydata = data['ydata']
    out = minimize(residual, params, args=(xdata, ydata))
    p = out.params
    e0 = p['e0_real'] + p['e0_imag'] * 1.j
    A = p['A_real'] + p['A_imag'] * 1.j
    # B = p['B_real'] + p['B_imag'] * 1.j
    B = 0
    return e0, A, B, out

def gap_test(S):
    d = DataSet(S, xkey='L', ykey='gap')
    def fitf(L, einf, Lcoef): return -einf - Lcoef/L

    d.plot(linestyle='')
    ex = ['L']
    for x, y in S.slicer.get_sets(exclude=ex).items():
        print(x)
        d = {}
        d['xdata'] = [s['L'] for s in y]
        d['ydata'] = [s['gap'] for s in y]
        s = y[0]
        N = s['N']
        # s.init_model()
        v = s['gap_inf']
        print(f'gap_exact = {v}')
        e0, A, B, out = lmtest(d)
        print(f'gap_interpolated = {e0}')
        # print(f'gap = {s["gap"]}')
        print(A, B)

def gap_test2(S):
    d = DataSet(S, xkey='L', ykey='gap')
    def fitf(L, einf, Lcoef): return -einf - Lcoef/L

    i = 0
    for s in S.systems:
        i += 1
        z = (1. - s['omega'])
        pfs = numpy.array(s['pfs']) * s['L'] * (1. - s['omega'])
        pfs /= z
        xdata = numpy.real(pfs)
        ydata = numpy.imag(pfs)
        plt.clf()
        plt.plot(xdata, ydata, linestyle='', marker='o', markersize=2)
        g = s['gap_inf'] / z
        g = fix_pfs([g], s['N'])[0]
        plt.plot([g.real], [g.imag], linestyle='', marker='x', markersize=15, linewidth=5)
        f = f"{i:0>3}_{s['gamma']}_{s['phi']}.png"
        f = os.path.join(S.graph_dir, f)
        plt.savefig(f)

import numpy as np

def to_raster(X, y):
    X = numpy.array(X)
    y = numpy.array(y)
    def deduce_raster_params():
        """
        Computes raster dimensions based on min/max coordinates in X
        sample step computed from 2nd - smallest coordinate values
        """
        unique_sorted = np.vstack((np.unique(v) for v in X.T)).T
        d_min = unique_sorted[0] # x min, y min
        d_max = unique_sorted[-1] # x max, y max
        d_step = unique_sorted[1]-unique_sorted[0] # x, y step
        nsamples = (np.round((d_max - d_min) / d_step) + 1).astype(int)
        return d_min, d_max, d_step, nsamples

    d_min, d_max, d_step, nsamples = deduce_raster_params()
    # Allocate matrix / tensor for raster. Allow y to be vector (e.g. RGB triplets)
    A = np.full((*nsamples, 1 if y.ndim==1 else y.shape[-1]), np.NaN)
    # Compute index for each point in X
    ind = np.round((X - d_min) / d_step).T.astype(int)
    # Scalar/vector values assigned over outer dimension
    print(d_min, d_max, d_step, nsamples)
    print(X.shape)
    print(y.shape)
    print(A.shape)
    print(ind.shape)
    ind = ind.T
    A[list(ind)] = y  # cell id
    # Prepare extent in imshow format
    extent = np.vstack((d_min, d_max)).T.ravel()
    return A, extent

def gap_test3(S):
    M = 2
    gammas = S.slicer.get_possible_values('gamma')
    phis = S.slicer.get_possible_values('phi')
    for i in range(M):
        X = []
        y = []
        for s in S.systems:
            X.append(numpy.array([s['gamma'], s['phi']]))
            y.append(abs(s['pfs'][i]))

        X = numpy.array(X)
        d_min, d_max, d_step, nsamples = deduce_raster_params(X)
        ind = np.round((X - d_min) / d_step).T.astype(int)
        ind = ind.T
        z = numpy.zeros((len(gammas), len(phis)))
        print(z.shape)
        for i in range(len(X)):
            a, b = ind[i]
            z[a][b] = y[i]


        extent = np.vstack((d_min, d_max)).T.ravel()
        z = numpy.array(z)
        plt.imshow(z, extent=extent)
        plt.show()


def transform1(s, N):
    eigs = s['eigs']
    lamb = s['lambda']
    phi = s['phi']
    L = s['L']
    g = lamb * numpy.exp(phi * 2.j * numpy.pi / N)
    eigs = [e ** (1./N) * g for e in eigs]
    u = 1./ s['lambda_scale']
    u *= numpy.exp(-phi * 1.j * numpy.pi / N)
    eigs = [e * u for e in eigs]
    e0 = 0
    eigs = fix_pfs(eigs, N)
    # eigs = [abs(x) for x in eigs]
    eigs.sort()
    for p in eigs:
        e0 -= p
    return eigs, e0

def mtest1(S):
    for x, y in S.slicer.get_sets(exclude=['gamma', 'phi']).items():
        print(x)
        for m in range(S['M']):
            for N in S['Ns']:
                print(N, m)
                X = []
                Y = []
                for s in y:
                    eigs, e0 = transform1(s, N)
                    X.append(numpy.array([s['phi'], s['gamma']]))
                    Y.append(abs(eigs[m]))
                z, ex = rastify(X, Y)
                plt.clf()
                aspect_ratio = 1
                aspect = abs(ex[1]-ex[0])/abs(ex[3]-ex[2]) * aspect_ratio
                P = plt.imshow(z, vmin=0, extent=ex, aspect=aspect)
                plt.colorbar(P)
                lab = make_ascii_label(x)
                f = f'{lab}_N={N}_m={m}.png'
                f = os.path.join(S.graph_dir, f)
                plt.xlabel(r'$\phi$')
                plt.ylabel(r'$\gamma$')
                plt.savefig(f)

def mtest2(S):
    for m in range(S['M']):
        for N in S['Ns']:
            print(N, m)
            X = []
            Y = []
            Y_e0 = []
            Y_error = []
            biggest_error = (0, 0, 0)
            for x, y in S.slicer.get_sets(exclude=['L']).items():
                xdata = []
                ydata = []
                e0_ydata = []
                for s in y:
                    eigs, e0 = transform1(s, N)
                    xdata.append(s['L'])
                    ydata.append(eigs[m])
                    e0_ydata.append(e0 / s['L'])
                d = {'xdata':xdata, 'ydata':ydata}
                fit, A, B, out = lmtest(d)
                s = y[0]
                X.append(numpy.array([s['phi'], s['gamma']]))
                Y.append(abs(fit))

                d = {'xdata':xdata, 'ydata':e0_ydata}
                fit, A, B, out = lmtest(d)
                Y_e0.append(abs(fit))
                D = rel_diff(fit, einf(s))
                Y_error.append(D)
                if D > biggest_error[0]:
                    biggest_error = (D, s['gamma'], s['phi'])
                    print(fit, einf(s))


                # omega = numpy.exp(2.j * numpy.pi / N)
                # gexact = gap_inf(s) / (1.-omega)
                # D = rel_diff(fit, gexact)
                # q, v=  transform1(y[-1], N)
                # if D > 1E-4:
                #     print(fit, gexact, D)
                # print(fit, gexact, D)
            print(f'be = {biggest_error}')
            for name, Z in [
                    ('gap', Y),
                    ('e0', Y_e0),
                    ('error', Y_error),
                    ]:
                z, ex = rastify(X, Z)
                plt.clf()
                aspect_ratio = 1
                aspect = abs(ex[1]-ex[0])/abs(ex[3]-ex[2]) * aspect_ratio
                P = plt.imshow(z, vmin=0, extent=ex, aspect=aspect, origin='lower')
                plt.colorbar(P)
                f = f'{name}_N={N}_m={m}.png'
                f = os.path.join(S.graph_dir, f)
                plt.xlabel(r'$\phi$')
                plt.ylabel(r'$\gamma$')
                plt.savefig(f)


def proggy(i, N):
    q = 100
    def f(u):
        return numpy.floor(q*u/N)
    if not f(i) == f(i+1):
        print(f'Progress = {f(i)/q}')

def mtest3(S):
    for m in range(S['M']):
        X = []
        Y = []
        Y_e0 = []
        Y_error = []
        biggest_error = (0, 0, 0)
        sets = S.slicer.get_sets(exclude=['L'])
        i = 0
        for x, y in sets.items():
            i += 1
            proggy(i, len(sets))
            xdata = []
            ydata = []
            e0_ydata = []
            s0 = y[0]
            for s in y:
                eigs = s['pfs']
                e0 = s['e_0']
                xdata.append(s['L'])
                ydata.append(eigs[m] * s['L'])
                e0_ydata.append(e0)

            d = {'xdata':xdata, 'ydata':ydata}
            fit, A, B, out = lmtest(d)
            X.append(numpy.array([s['phi'], s['gamma']]))
            Y.append(abs(fit))

            d = {'xdata':xdata, 'ydata':e0_ydata}
            fit, A, B, out = lmtest(d)
            Y_e0.append(abs(fit))
            D = rel_diff(fit, einf(s))
            Y_error.append(D)
            if D > biggest_error[0]:
                biggest_error = (D, s['gamma'], s['phi'])
                print(fit, einf(s))


            # omega = numpy.exp(2.j * numpy.pi / N)
            # gexact = gap_inf(s) / (1.-omega)
            # D = rel_diff(fit, gexact)
            # q, v=  transform1(y[-1], N)
            # if D > 1E-4:
            #     print(fit, gexact, D)
            # print(fit, gexact, D)
        print(f'be = {biggest_error}')
        for name, Z in [
                ('gap', Y),
                ('e0', Y_e0),
                ('error', Y_error),
                ]:
            z, ex = rastify(X, Z)
            plt.clf()
            aspect_ratio = 1
            aspect = abs(ex[1]-ex[0])/abs(ex[3]-ex[2]) * aspect_ratio
            P = plt.imshow(z, vmin=0, extent=ex, aspect=aspect, origin='lower')
            plt.colorbar(P)
            f = f'{name}_N={s["N"]}_m={m}.png'
            f = os.path.join(S.graph_dir, f)
            plt.xlabel(r'$\phi$')
            plt.ylabel(r'$\gamma$')
            plt.savefig(f)


def oe_2(S):
    ex = []
    for x, y in S.slicer.get_sets(exclude=ex).items():
        print(x)
        plt.clf()
        for s in y:
            pltargs = dict(linestyle='', color='blue', marker='o', zorder=5)
            energies = make_pf_spectrum2(s['e_0'], s['pfs'], s['N'])
            xdata = [x.real for x in energies]
            ydata = [x.imag for x in energies]
            plt.plot(xdata, ydata, **pltargs)

            pltargs = dict(linestyle='', color='green', marker='+', zorder=11, markersize=15, linewidth=2)
            e0 = s['e_0']
            print(e0)
            xdata = [e0.real]
            ydata = [e0.imag]
            plt.plot(xdata, ydata, **pltargs)

            pltargs = dict(linestyle='', color='red', marker='+', zorder=11, markersize=15, linewidth=2)
            e0 = s['e_inf']
            print(e0)
            xdata = [e0.real]
            ydata = [e0.imag]
            plt.plot(xdata, ydata, **pltargs)
        f = f'{x}.png'
        f = os.path.join(S.graph_dir, f)
        plt.savefig(f)


def mtest3(S):
    for m in range(S['M']):
        X = []
        Y = []
        Y_e0 = []
        Y_error = []
        biggest_error = (0, 0, 0)
        sets = S.slicer.get_sets(exclude=['L'])
        i = 0
        for x, y in sets.items():
            i += 1
            proggy(i, len(sets))
            xdata = []
            ydata = []
            e0_ydata = []
            s0 = y[0]
            for s in y:
                eigs = s['pfs']
                e0 = s['e_0']
                xdata.append(s['L'])
                ydata.append(eigs[m] * s['L'])
                e0_ydata.append(e0)

            d = {'xdata':xdata, 'ydata':ydata}
            fit, A, B, out = lmtest(d)
            X.append(numpy.array([s['phi'], s['gamma']]))
            Y.append(abs(fit))

            d = {'xdata':xdata, 'ydata':e0_ydata}
            fit, A, B, out = lmtest(d)
            Y_e0.append(abs(fit))
            D = rel_diff(fit, einf(s))
            Y_error.append(D)
            if D > biggest_error[0]:
                biggest_error = (D, s['gamma'], s['phi'])
                print(fit, einf(s))


            # omega = numpy.exp(2.j * numpy.pi / N)
            # gexact = gap_inf(s) / (1.-omega)
            # D = rel_diff(fit, gexact)
            # q, v=  transform1(y[-1], N)
            # if D > 1E-4:
            #     print(fit, gexact, D)
            # print(fit, gexact, D)
        print(f'be = {biggest_error}')
        for name, Z in [
                ('gap', Y),
                ('e0', Y_e0),
                ('error', Y_error),
                ]:
            z, ex = rastify(X, Z)
            plt.clf()
            aspect_ratio = 1
            aspect = abs(ex[1]-ex[0])/abs(ex[3]-ex[2]) * aspect_ratio
            P = plt.imshow(z, vmin=0, extent=ex, aspect=aspect, origin='lower')
            plt.colorbar(P)
            f = f'{name}_N={s["N"]}_m={m}.png'
            f = os.path.join(S.graph_dir, f)
            plt.xlabel(r'$\phi$')
            plt.ylabel(r'$\gamma$')
            plt.savefig(f)

def p1(X, Y, name, **kwargs):
    z, ex = rastify(X, Y)
    plt.clf()
    aspect_ratio = 1
    aspect = abs(ex[1]-ex[0])/abs(ex[3]-ex[2]) * aspect_ratio
    P = plt.imshow(z, extent=ex, aspect=aspect, origin='lower', **kwargs)
    plt.colorbar(P)
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\gamma$')
    plt.savefig(name)

def get_excitations(s):
    exs_a = []
    exs_r = []
    for p in s['pfs']:
        a = 1000
        r = 1000
        for i in range(s['N']-1):
            u = 1. - s['omega'] ** (i+1)
            ex = u * p * s['L']
            if abs(ex) < abs(a):
                a = ex
            if ex.real < r.real:
                r = ex
        exs_a.append(a)
        exs_r.append(r)
    exs_a = sorted(exs_a, key = lambda x: abs(x))
    exs_r = sorted(exs_r, key = lambda x: x.real)
    s['exs_a'] = exs_a
    s['exs_r'] = exs_r

def mtest4(S):
    sets = S.slicer.get_sets(include=['L'])
    for s in S.systems:
        get_excitations(s)
    for a, b in sets.items():
        for x, y in b.get_sets(exclude = ['gamma', 'phi']).items():
            nsys = len(y)
            X = numpy.full((nsys, 2), 0.)
            for i in range(nsys):
                s = y[i]
                X[i] = s['phi'], s['gamma']

            keys = ['e_inf', 'gap_inf']
            for k in keys:
                Y = numpy.full((nsys), 0.)
                for i in range(nsys):
                    s = y[i]
                    Y[i] = abs(s[k])
                name = os.path.join(S.graph_dir, f'{k}.png')
                p1(X, Y, name)
        break

    sets = S.slicer.get_sets(exclude=['gamma','phi'])
    ind = 0
    for x, y in sets.items():
        ind += 1
        lab = make_ascii_label(x)
        print(lab)
        f = f'{ind:0>3}_{lab}.png'

        # pfs 1 2 3
        # e0
        # EP test
        nsys = len(y)
        X = numpy.full((nsys, 2), 0.)
        for i in range(nsys):
            s = y[i]
            X[i] = s['phi'], s['gamma']

        for j in range(3):
            for k in ['pfs', 'exs_a', 'exs_r']:
                Y = numpy.full((nsys), 0.)
                for i in range(nsys):
                    s = y[i]
                    Y[i] = abs(s[k][j])
                name = os.path.join(S.graph_dir, f'{k}{j}_{f}.png')
                # p1(X, Y, name, vmin=0)
                p1(X, Y, name)

        M = 2
        Mm = 2
        for d in range(Mm):
            for m in range(M):
                Y = numpy.full((nsys), 0.)
                for i in range(nsys):
                    s = y[i]
                    a = s['exs_a'][m]
                    b = s['exs_a'][m+d+1]
                    a = s['pfs'][m] * s['L']
                    b = s['pfs'][m+1+d] * s['L']
                    Y[i] = abs(a-b)
                    # best = 100
                    # for k in range(s['N']-1):
                    #     z = abs(a - b * (s['omega']**k))
                    #     if z < best:
                    #         best = z
                    # Y[i] = best
                name = os.path.join(S.graph_dir, f'pfdelta{m}_{d}_{f}.png')
                p1(X, Y, name)
        # p1(X, Y, name, vmin=0)


        Y = numpy.full((nsys), 0.)
        for i in range(nsys):
            s = y[i]
            g = s['gamma']
            if g < 0.5:
                Y[i] = s['pfs'][1]
            else:
                Y[i] = s['pfs'][0]
        name = os.path.join(S.graph_dir, f'gap_{f}.png')
        p1(X, Y, name)

        # keys = ['e_0', 'smallest_pf_diff']
        keys = ['e_0']
        for k in keys:
            Y = numpy.full((nsys), 0.)
            for i in range(nsys):
                s = y[i]
                Y[i] = abs(s[k])
            name = os.path.join(S.graph_dir, f'{k}_{f}.png')
            p1(X, Y, name)
