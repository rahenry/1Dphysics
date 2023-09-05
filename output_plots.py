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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from fp_k import k_inf

label_data = {
    'lambda_real' : r'$\Re(\lambda)$',
    'lambda_imag' : r'$\Im(\lambda)$',
    'k_real' : r'$\Re(k)$',
    'k_imag' : r'$\Im(k)$',
    }
def interpret_plot(S, i):
    p = S[f'plot{i}']
    if not p: return

    res = {}
    pnew = []
    for x in p:
        if '$' in x:
            pnew += S[x.replace('$', '')]
        else:
            pnew.append(x)

    p = pnew
    for x in p:
        if '=' in x:
            a, b = x.split('=')
            res[a] = b
        else:
            res[x] = 1

    res['target'] = {}
    res['target']['sset_name'] = S.name
    res['find_index'] = 0
    for x in S[f'plot{i}_target']:
        a, b = x.split('=')
        if a == 'from':
            res['target']['sset_name'] = S.added_sets[int(b)].name
            res['find_index'] = int(b)
            # print(b)
            # print(res)
    ind = res['find_index']
    for x in S[f'plot{i}_target']:
        a, b = x.split('=')
        if a == 'from':
            continue
        if '!' in b:
            b = b.replace('!', '')
            b = int(b)
            b = S.added_sets[ind].config_systems[a][b]
        res['target'][a] = b
    return res

def getty(k, pdata, S):
    if k not in pdata['target']:
        return S.config_systems[k][0]
    v = pdata['target'][k]
    t = type(S.config_systems[k][0])
    return t(v)

def get_eps_2d_data(S, pdata):
    S = S.added_sets[pdata['find_index']]
    N = getty('N', pdata, S)
    L = getty('L', pdata, S)
    cid = (L, N)
    res = S.eps_data[cid]
    return res['z'], res['ex']

def get_2d_data(S, pdata):
    if pdata['zkey'] == 'eps':
        return get_eps_2d_data(S, pdata)
    z = S.slicer.find(pdata['target'])
    if not z:
        print(pdata['target'])
        print('Nothing found')
        exit()
    s0 = z.systems[0]
    X = []
    Y = []
    for s in z.systems:
        x = s[pdata['xkey']]
        y = s[pdata['ykey']]
        z = s[pdata['zkey']]
        X.append((x,y))
        Y.append(z)
    z, ex = rastify(X, Y)
    return z, ex

def insert_cbar(fig, ax, p):
    ticks = numpy.linspace(p.norm._vmin, p.norm._vmax, 3)
    ticks[0] = 0
    axins1 = inset_axes(ax,
                    width="4%",  # width = 50% of parent_bbox width
                    height="40%",  # height : 5%
                    loc='lower right')


    fig.colorbar(p, cax=axins1, orientation="vertical", ticks=ticks)
    axins1.yaxis.set_ticks_position("left")
    # c_ax.yaxis.set_ticks_position('left')

imargs_default = {
    'origin' : 'lower',
    }

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

def plot_eps(S, ax, d, ex):
    S = S.added_sets[d['find_index']]
    plotargs = {
        'marker' : 'x',
        'markersize' : 4,
        'color' : 'white',
        'linewidth' : 0.2,
        'alpha' : 0.20,
        }
    N = getty('N', d, S)
    L = getty('L', d, S)
    cid = (L, N)
    eps = S.eps_data[cid]['minima_k']
    if not d['zkey'] == 'eps':
        eps = S.eps_data[cid]['minima']
        # a = e[0] * numpy.exp(e[1] * 2.j * numpy.pi / N)
        eps = rotate_minima(eps, N)
    for a in eps:
        x = a.real
        y = a.imag
        if x > ex[1] or x < ex[0] or y > ex[3] or y < ex[2]:
            continue
        print(a)
        ax.plot([a.real], [a.imag], **plotargs)

def grid_plot(S):
    i = 0
    pdata = []
    while (True):
        i += 1
        p = interpret_plot(S, i)
        if not p:
            break
        pdata.append(p)

    use_cbars = 0
    w2 = 0
    for d in pdata:
        if 'cbar' in d:
            use_cbars = 1
            w2 = 0.4

    n = S['plot_cols']
    m = len(pdata) // n + len(pdata) % n
    wr = [1, w2] * n
    hr = [1] * m
    # gs = gridspec.GridSpec(m, n, width_ratios=wr, height_ratios=hr)
    gs = gridspec.GridSpec(m, n)

    scale = 4
    plt.clf()
    fig = plt.figure(figsize = (n*scale,m*scale))
    # gs.update(wspace=0.025, hspace=0.15)
    # gs.update(left=0.1, right=0.95)
    # gs.tight_layout(fig)

    for (ind_gs, d) in enumerate(pdata):
        ax = plt.subplot(gs[ind_gs])
        print(ind_gs, d)

        if d['type'] == 'image':
            # print(z[0,0])
            z, ex = get_2d_data(S, d)
            if 'zpow' in d:
                z = numpy.power(z, float(d['zpow']))
            aspect_ratio = 1
            aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio

            imargs = dict(imargs_default)
            if 'vmin' in d:
                imargs['vmin'] = d['vmin']
            p = ax.imshow(z, extent=ex, aspect=aspect, **imargs)
            if 'cbar' in d:
                insert_cbar(fig, ax, p)
            if d['xkey'] in label_data:
                ax.set_xlabel(label_data[d['xkey']])
            if d['ykey'] in label_data:
                ax.set_ylabel(label_data[d['ykey']])
            if d['show_eps']:
                plot_eps(S, ax, d, ex)

    plt.tight_layout()
    f = 'grid_plot'
    savefig(S, f)

def plot_trajectories(S):
    default_x = 6.4
    default_y = 4.8
    do_rotations = S.get('do_rotations', 0)
    scale = 1
    plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)
    zkeys = S['trajectory_zkeys']
    xkeys = S['trajectory_xkeys']
    npoints = S['trajectory_npoints']
    plot_kinf = 0
    if not npoints or npoints < 1:
        npoints = None
    for k in zkeys:
        for xkey in xkeys:
            for x, y in S.slicer.get_sets(exclude=xkey).items():
                plt.clf()
                cm = matplotlib.cm.get_cmap('cool')
                ind = 0

                for s in y:
                    z = s[k]
                    xdata = numpy.real(z)
                    ydata = numpy.imag(z)
                    cind = ind / len(y)
                    c = cm(cind)
                    if npoints:
                        c = cm(ind / npoints)
                    plt.plot(xdata, ydata, linestyle='', marker='o', markersize=1.3, color=c, zorder=10)

                    ind += 1
                    if npoints and ind > npoints:
                        break

                    if plot_kinf:
                        kinf = k_inf(ind, s['L'], s['N'], s['lambda_full'])
                        plt.plot([kinf.real], [kinf.imag], linestyle='', marker='x', color='red', markersize=5, zorder=0)

                # plt.colorbar()
                m = matplotlib.cm.ScalarMappable(cmap=cm)
                plt.colorbar(m)
        # fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
                f = os.path.join(S.graph_dir, f'{k}_{xkey}_{labels.make_ascii_label(x)}.png')
                plt.xlabel(r'$\text{Re}(\varepsilon)$')
                plt.ylabel(r'$\text{Im}(\varepsilon)$')
                # add_radial_lines(plt.gca(), y[0], zorder=1)
                plt.savefig(f)
                plt.clf()
