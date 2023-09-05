#!/usr/bin/env python3

import matplotlib.pyplot as plt
import os
from Slicer import SystemSlicer
from pf_tester import pf_test
from tenpy.linalg.np_conserved import tensordot, inner, norm
from rah_utility import mkdir, matplotlib_set_aspect, extrema, rastify
import graph_settings
import matplotlib
import numpy
from figures import *
from args import Q
from data_set import DataSet, plot_multiple
import cProfile, pstats, io
from pstats import SortKey
from exact_fp import exact_parafermions_matrix, make_pf_spectrum2, rhoX, make_coords
from fp_more_matrices import *
from fp_k import *
from fp_polynomial import *
from matplotlib.offsetbox import AnchoredText
import matplotlib.gridspec as gridspec
from polytest import sympy_eps2
from point_generator import PointGenerator
from rah_mpl import draw_label
import labels
import matplotlib.transforms as mtransforms
import string
import matplotlib.gridspec as gridspec
def cl_eig_info(S):
    basename = S.graph_dir
    # for x, y in S.slicer.get_sets(exclude=['lambda', 'phi']).items():
    for a, b in S.slicer.get_sets().items():
        s = b[0]

        f0 = labels.make_ascii_label(a)
        keys = ['deltaLR', 'overlapLRnc', 'overlapRR']

        for k in keys:
            plt.clf()
            vals = s[k]
            f = os.path.join(S.graph_dir, f'{f0}_{k}.png')
            p = plt.imshow(vals, origin='lower', vmin=0)
            plt.gca().xaxis.set_major_locator(plt.NullLocator())
            plt.gca().yaxis.set_major_locator(plt.NullLocator())
            plt.savefig(f)

def cl_eig_plot_set(S):
    targs = S['output_ep_targets']
    keys = ['none', 'deltaLR', 'overlapLRnc', 'overlapRR']
    m = len(keys) + 1
    n = len(targs) // 3
    gs1 = gridspec.GridSpec(n, m)

    scale = 1.5
    fig = plt.figure(figsize = (m*scale,n*scale))
    gs1 = gridspec.GridSpec(n, m)
    gs1.update(wspace=0.125, hspace=0.05) # set the spacing between axes.

    labels = {
        'EP' : r'EP',
        'EP_PT' : r'EP, $PT$-symmetric',
        'EP_trivial' : r'EP, trivial',
        }


    for j in range(len(targs) // 3):
        N, L, ind = targs[3*j:3*j+3]
        cid = (L, N)
        eps_data = S.eps_data[cid]['minima']
        eps = []
        for lamb, phi in eps_data:
            if phi < 0 or phi > 0.51:
                continue
            eps.append((lamb, phi))
        eps.sort(key=lambda x: x[-1])
        lamb, phi = eps[ind]
        d = {'N':N, 'L':L, 'lambda':lamb, 'phi':phi}
        y = S.slicer.find(d)
        s = y[0]

        for i, k in enumerate(keys):
            ind_gs = i + m * j
            ax = plt.subplot(gs1[ind_gs])
            if k == 'none':
                ax.axis('off')
                lab = S['output_ep_labels'][j]
                lab = labels[lab]
                text = fr'{lab} \\ $N={N},\; L={L}$\\ $\lambda={lamb:.4f}, $\\ $\phi={phi:.4f}$'
                # text = r'a \\ b'
                text= r'\begin{tabular}{@{}l@{}}' + text +r'\end{tabular}'
                ax.text(0.0, 0.5, text, va='center', ha='left')
                continue
            print(i, j)
            vals = s[k]
            p = ax.imshow(vals, origin='lower', vmin=0)

            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())

            trans = mtransforms.ScaledTranslation(-0/72, 7/72, fig.dpi_scale_trans)
            # label = labels[i]
            # if (i % 3) == 0:
            #     label = f'({labs[ind]})'
            #     ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            #             fontsize='medium', va='bottom', fontfamily='serif')
            #     ind += 1

    f = os.path.join(S.graph_dir, 'eig_plot_set.png')
    plt.savefig(f)
    f = os.path.join(S.pgf_dir, 'eig_plot_set.pgf')
    plt.savefig(f)


def cl_plot_k(S):
    for cid in S.eps_data:
        d = S.eps_data[cid]
        L, N = cid
        print(cid)
        z = d['z']
        z = numpy.power(z, S['zpower'])
        ex = d['ex']
        plt.clf()
        aspect_ratio = 1
        aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio
        plt.imshow(z, extent=d['ex'], aspect=aspect)
        f = os.path.join(S.graph_dir, f'k_{cid}.png')
        plt.savefig(f)

def get_q(pfs, N):
    L = len(pfs)
    pfs = sorted(pfs, key=lambda x: abs(x))
    d0 = 10000
    for j in range(L-1):
        # pfs = s['pfs']
        p0 = pfs[0]
        p1 = pfs[j+1]
        for i in range(N):
            p1a = p1 * numpy.exp(2.j * numpy.pi * i / N)
            d = abs(p0 - p1a)
            if d < d0:
                d0 = d
                yt1 = float(i)

    return d0
    p0 = pfs[0]
    q = 10000
    for p in pfs[1:]:
        a = abs(p0 - p)
        if a < q:
            q = a

def rotate_minima(minima, N):
    res = []
    for m in minima:
        a = m[0] * numpy.exp(2.j * m[1] * numpy.pi / N)
        for i in range(N):
            b = a * numpy.exp(2.j * numpy.pi * i / N)
            failed = 0
            for r in res:
                if rel_diff(r, b) < 1E-5:
                    failed = 1
            if not failed:
                res.append(b)
    return res

def cl_plot_grid(S):
    plotargs = {
        'marker' : 'x',
        'markersize' : 4,
        'color' : 'white',
        'linewidth' : 0.2,
        'alpha' : 0.20,
        }
    imargs = {
        'vmin' : 0,
        'origin' : 'lower',
        }

    for cid in S.eps_data:
        L, N = cid
        d = S.eps_data[cid]
        y = S.slicer.find(dict(L=L, N=N))
        n = S['grid_npoints']
        X = []
        Y = numpy.full((n*n), 0.)
        i = 0
        for s in y:
            x = numpy.array((s['lambda_real'], s['lambda_imag']))
            X.append(x)
            pfs = s['pfs']
            q = get_q(pfs, N)

            Y[i] = q
            i += 1

        z, ex = rastify(X, Y)
        f = os.path.join(S.graph_dir, f'grid_{cid}.png')
        aspect_ratio = 1
        aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio

        plt.clf()
        plt.xlabel(r'$\Re(\lambda)$')
        plt.ylabel(r'$\Im(\lambda)$')
        p = plt.imshow(z, extent=ex, aspect=aspect, **imargs)
        plt.colorbar(p)

        if S['draw_minima']:
            for m in rotate_minima(d['minima'], N):
                a = m
                # a = m[0] * numpy.exp(2.j * m[1] * numpy.pi / N)
                plt.plot([a.real], [a.imag], **plotargs)

        f = os.path.join(S.graph_dir, f'grid_{cid}.png')
        plt.savefig(f)
        f = os.path.join(S.pgf_dir, f'grid_{cid}.pgf')
        plt.savefig(f)

    for cid in S.eps_data:
        L, N = cid
        d = S.eps_data[cid]
        z = d['z']
        ex = d['ex']
        aspect_ratio = 1
        aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio

        plt.clf()
        plt.xlabel(r'$\Re(k)$')
        plt.ylabel(r'$\Im(k)$')
        z = numpy.power(z, s['zpow_k'])
        p = plt.imshow(z, extent=ex, aspect=aspect, **imargs)
        plt.colorbar(p)

        if S['draw_minima']:
            for m in d['minima_k']:
                a = m
                plt.plot([a.real], [a.imag], **plotargs)

        f = os.path.join(S.graph_dir, f'grid_k_{cid}.png')
        plt.savefig(f)
        f = os.path.join(S.pgf_dir, f'grid_k_{cid}.pgf')
        plt.savefig(f)
        print(cid)

def cl_plot_grid2(S):
    plotargs = {
        'marker' : 'x',
        'markersize' : 4,
        'color' : 'white',
        'linewidth' : 0.2,
        'alpha' : 0.20,
        }
    imargs = {
        'vmin' : 0,
        'origin' : 'lower',
        }

    for x, y in S.slicer.get_sets(exclude=['lambda_real', 'lambda_imag']).items():
        print(x, len(y))
        X = []
        Y = []
        for s in y:
            N = s['N']
            X.append((s['lambda_real'], s['lambda_imag']))
            pfs = s['pfs']
            q = get_q(pfs, N)
            Y.append(q)

        z, ex = rastify(X, Y)
        aspect_ratio = 1
        aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio

        plt.clf()
        plt.xlabel(r'$\Re(\lambda)$')
        plt.ylabel(r'$\Im(\lambda)$')
        p = plt.imshow(z, extent=ex, aspect=aspect, **imargs)
        plt.colorbar(p)

        if S['draw_minima']:
            for m in rotate_minima(d['minima'], N):
                a = m
                # a = m[0] * numpy.exp(2.j * m[1] * numpy.pi / N)
                plt.plot([a.real], [a.imag], **plotargs)

        f = os.path.join(S.graph_dir, f'grid_{x}.png')
        plt.savefig(f)
        f = os.path.join(S.pgf_dir, f'grid_{x}.pgf')
        plt.savefig(f)


def k_test(S):
    print(S)
    # return
    for s in S.systems:
        pfs = s['pfs']
        ks = s['ks']
        L = s['L']
        lamb = s['lambda_full']
        N = s['N']

        print(s)
        ls = s['lambda_scale']
        print(ls)
        print(s['E_0']/L)
        print(s['e_inf'])

        inf_data = []
        for i in range(len(pfs)):
            p = pfs[i]
            p = p * ls

            ki = k_inf(i+1, L, N, lamb)
            pi = eps_k(ki, lamb, N)
            inf_data.append((ki, pi))
        inf_data.sort(key=lambda x: abs(x[1]))
        i = 0
        nprint = -1
        for (ki, pi) in inf_data:
            p = pfs[i]
            p = p * ls
            k = eps_to_k(p, lamb, N)
            # print(k, ki, p, pi)
            i += 1
            if nprint > 0 and i > nprint:
                break
        print(123)
        kL = numpy.pi + 1.j*numpy.log(lamb ** (-N/2))
        kLi = eps_k(kL, lamb, N)
        print(kL, kLi)
        kL = numpy.pi - 1.j*numpy.log(lamb ** (-N/2))
        kLi = eps_k(kL, lamb, N)
        print(kL, kLi)

        print('------')
        alpha = -(lamb ** (N/2.))
        def f1(x):
            res = numpy.sin((L+1)*x)/numpy.sin(L*x) -alpha
            return res
        for i in range(len(pfs)):
            p = pfs[i]
            z = 1
            z = lamb
            # z = s['lambda'] / s['lambda_full']
            # z = 1./z
            # z = 1./lamb
            p = p * ls * z

            k = eps_to_k(p, lamb, N)
            # print(k, numpy.exp(1.j*k), alpha)
            print(k, f1(k))
