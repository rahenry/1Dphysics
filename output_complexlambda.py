#!/usr/bin/env python3

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
from fp_k import eps_to_k, eps_k

def to_bins(Nbins, points, scale=0.5, offset=0):
    biggest = 0
    for p in points:
        m = max(p.real, p.imag)
        if m > biggest: biggest = m
    biggest *= scale
    bins = numpy.zeros((Nbins, Nbins))
    for p in points:
        z = numpy.array([p.real, p.imag])
        z = z / biggest * Nbins/2 + Nbins/2
        z = [round(x) for x in z]
        fail = 0
        for y in z:
            if y >= Nbins or y < 0:
                fail = 1
        if fail: continue
        bins[z[1]][z[0]] += 1
    m = numpy.max(bins)
    for i in range(Nbins):
        for j in range(Nbins):
            if bins[i][j] > 0:
                bins[i][j] += offset * m
    return bins, biggest

def fin_complexlambda1(S):
    # fig, axs = plt.subplots(2, 3)
    # fig, axs = plt.subplots(2, 3, sharex=True, sharey=True)
    # plt.figure(figsize=(40,60), dpi=100)
    nrows = 2
    ncols = 3
    gs1 = gridspec.GridSpec(nrows, ncols)
    gs1.update(wspace=0, hspace=0)
    ax_ind = 0
    all_axs = []
    # for z in axs: all_axs += [z0 for z0 in z]
    for x, y in S.slicer.get_sets().items():
        I = ax_ind // ncols
        J = ax_ind % ncols
        s = y[0]
        spect2, pfs2 = make_spectrum(s, mtype='m2')
        pfs_poly = find_eps_poly(s)
        pfs_poly = [p/s['L'] for p in pfs_poly]
        e0 = 0
        for p in pfs_poly:
            e0 -= p
        spect_poly = make_pf_spectrum2(e0, pfs_poly, s['N'])
        # continue
        xdata_diag = []
        ydata_diag = []
        xdata_pf = []
        ydata_pf = []
        for i in range(len(spect_poly)):
            a = spect_poly[i]
            xdata_pf.append(a.real)
            ydata_pf.append(a.imag)
            b = s['energies'][i]
            xdata_diag.append(b.real)
            ydata_diag.append(b.imag)


        # ax = all_axs[ax_ind]
        ax = plt.subplot(gs1[ax_ind])
        # plt.axis('on')
        ax.plot(xdata_pf, ydata_pf, linestyle='', marker='x', color='red', markersize=7)
        ax.plot(xdata_diag, ydata_diag, linestyle='', marker='o', color='blue', markersize=2.0)
        ax.set_xlim([-0.6, 0.6])
        ax.set_ylim([-0.6, 0.6])
        # matplotlib_set_aspect(ax, 1)
        if not I == 1:
            ax.set_xticklabels([])
            ax.set_xticks([])
        else:
            ax.set_xticks(numpy.linspace(-0.5, 0.5, 3))
            ax.set_xlabel(r"$\Re(E)$")
        if not J == 0:
            ax.set_yticks([])
            ax.set_yticklabels([])
        else:
            ax.set_yticks(numpy.linspace(-0.5, 0.5, 3))
            ax.set_ylabel(r"$\Im(E)$")

        # ax.set_title(rf'$\phi = {s["phi"]}$')
        u = f'{s["phi"]:.3f}'
        u = u.rstrip('0').rstrip('.')
        lab = rf'$\phi = {u}$'
        loc = "upper right"
        text_box = AnchoredText(lab, frameon=True, loc=loc, pad=0)
        text_box.patch.set_boxstyle("square,pad=0.15")
        plt.setp(text_box.patch, linestyle='', facecolor='lightgrey')
        ax.add_artist(text_box)
        ax_ind += 1
    f = os.path.join(S.graph_dir, f'spectra.png')
    # fig.tight_layout()
    # fig.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(f)


def fin_complexlambda2(S):
    nrows = 2
    ncols = 3
    gs1 = gridspec.GridSpec(nrows, ncols)
    gs1.update(wspace=0, hspace=0)
    ax_ind = 0
    all_axs = []
    # for z in axs: all_axs += [z0 for z0 in z]
    for x, y in S.slicer.get_sets().items():
        I = ax_ind // ncols
        J = ax_ind % ncols
        s = y[0]
        pfs_poly = sympy_eps2(s)
        e0 = 0
        for p in pfs_poly:
            e0 -= p
        spect_poly = make_pf_spectrum2(e0, pfs_poly, s['N'])
        xdata_diag = []
        ydata_diag = []
        xdata_pf = []
        ydata_pf = []
        for i in range(len(spect_poly)):
            a = spect_poly[i]
            xdata_pf.append(a.real)
            ydata_pf.append(a.imag)


        # ax = all_axs[ax_ind]
        ax = plt.subplot(gs1[ax_ind])
        # plt.axis('on')
        ax.plot(xdata_pf, ydata_pf, linestyle='', marker='o', color='blue', markersize=1)
        # ax.set_xlim([-0.6, 0.6])
        # ax.set_ylim([-0.6, 0.6])
        # matplotlib_set_aspect(ax, 1)
        if not I == 1:
            ax.set_xticklabels([])
            ax.set_xticks([])
        else:
            ax.set_xticks(numpy.linspace(-0.5, 0.5, 3))
            ax.set_xlabel(r"$\Re(E)$")
        if not J == 0:
            ax.set_yticks([])
            ax.set_yticklabels([])
        else:
            ax.set_yticks(numpy.linspace(-0.5, 0.5, 3))
            ax.set_ylabel(r"$\Im(E)$")

        # ax.set_title(rf'$\phi = {s["phi"]}$')
        u = f'{s["phi"]:.3f}'
        u = u.rstrip('0').rstrip('.')
        lab = rf'$\phi = {u}$'
        loc = "upper right"
        text_box = AnchoredText(lab, frameon=True, loc=loc, pad=0)
        text_box.patch.set_boxstyle("square,pad=0.15")
        plt.setp(text_box.patch, linestyle='', facecolor='lightgrey')
        ax.add_artist(text_box)
        ax_ind += 1
    f = os.path.join(S.graph_dir, f'spectra.png')
    # fig.tight_layout()
    # fig.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(f)


def fin_complexlambda3(S):
    # fig, axs = plt.subplots(2, 3)
    # fig, axs = plt.subplots(2, 3, sharex=True, sharey=True)
    # plt.figure(figsize=(40,60), dpi=100)
    nrows = 2
    ncols = 3
    # gs1 = gridspec.GridSpec(nrows, ncols)
    # gs1.update(wspace=0, hspace=0)
    # ax_ind = 0
    # all_axs = []
    # for z in axs: all_axs += [z0 for z0 in z]
    fig_index = 0
    for x, y in S.slicer.get_sets().items():
        fig_index += 1
        fig_ident = f'0000{fig_index}'[-3:]
        plt.clf()
        s = y[0]
        pfs_poly = s['pfs']
        xdata_pf = []
        ydata_pf = []
        for a in pfs_poly:
            xdata_pf.append(a.real)
            ydata_pf.append(a.imag)


        label = fr'$L = {s["L"]} \\ \gamma = {s["gamma"]:.2f} \\ \lambda = {s["lambda"]:.2f} \\ \phi = {s["phi"]:.2f}$'
        label = fr'$N = {s["N"]} \\ L = {s["L"]} \\ |\lambda| = {s["lambda"]:.2f} \\ \phi = {s["phi"]:.2f}$'

        # epsilon plot
        plt.plot(xdata_pf, ydata_pf, linestyle='', marker='o', color='blue', markersize=3)
        f = os.path.join(S.graph_dir, f'{fig_ident}_{x}_pfs.png')
        plt.xlabel(r'$\text{Re}(\epsilon)$')
        plt.ylabel(r'$\text{Im}(\epsilon)$')
        draw_label(plt.gca(), label)
        plt.savefig(f)



        # get point data
        P = PointGenerator(s)
        points, xdata, ydata = P.get_data()

        nbins = s["n_bins"]

        # full plot of energies
        bins, biggest = to_bins(nbins, points, scale=1.25, offset=0.1)
        plt.clf()
        ex = [-biggest, biggest, -biggest, biggest]
        plt.imshow(bins, extent=ex)
        plt.xlabel(r'$\text{Re}(E)$')
        plt.ylabel(r'$\text{Im}(E)$')
        f = os.path.join(S.graph_dir, f'{fig_ident}_{x}_random_bins.png')
        draw_label(plt.gca(), label)
        plt.savefig(f)

        # zoomed-in plot of energies
        bins, biggest = to_bins(nbins, points, scale=0.2, offset=0)
        plt.clf()
        ex = [-biggest, biggest, -biggest, biggest]
        plt.imshow(bins, extent=ex)
        plt.xlabel(r'$\text{Re}(E)$')
        plt.ylabel(r'$\text{Im}(E)$')
        draw_label(plt.gca(), label)
        f = os.path.join(S.graph_dir, f'{fig_ident}_{x}_random_zoomed.png')
        plt.savefig(f)



def fin_complexlambda_multi(S):
    ncols = S['multiplot_ncols']
    nrows = len(S.systems) // ncols + 1
    gs1 = gridspec.GridSpec(nrows, ncols)
    gs1.update(wspace=0, hspace=0)
    all_axs = []
    # for z in axs: all_axs += [z0 for z0 in z]
    default_x = 6.4
    default_y = 4.8
    scale = 0.3
    plt.figure(figsize=(default_x*scale*ncols, default_y*scale*nrows), dpi=100)

    limbase = [None, None]
    limbase = numpy.array([numpy.nan, numpy.nan])
    xlims = [numpy.array(limbase) for i in range(ncols)]
    ylims = [numpy.array(limbase) for i in range(nrows)]
    def update_lims(lims, i, data):
        a = min(data)
        b = max(data)
        l = lims[i]
        if numpy.isnan(l[0]) or a < l[0]:
            l[0] = a
        if numpy.isnan(l[1]) or b > l[1]:
            l[1] = b


    ax_ind = 0
    for x, y in S.slicer.get_sets().items():
        s = y[0]
        I = ax_ind // ncols
        J = ax_ind % ncols
        pfs = s['pfs']
        xdata = [x.real for x in pfs]
        print(len(xdata))
        ydata = [x.imag for x in pfs]
        # xdata, ydata = s['pfs_real'], s['pfs_imag']
        if not xdata: continue
        update_lims(xlims, J, xdata)
        update_lims(ylims, I, ydata)
        ax_ind += 1

    # print("Axis limits:")
    # for i in range(ncols):
    #     print(xlims[i])
    # for i in range(nrows):
    #     print(ylims[i])

    ax_ind = 0
    ex = S['xkey']
    for x, y in S.slicer.get_sets(exclude=ex).items():
    # for x, y in S.slicer.get_sets().items():
        for s in y:
            I = ax_ind // ncols
            J = ax_ind % ncols
            s = y[0]
            P = PointGenerator(s)
            # points, xdata, ydata = P.get_data()
            xdata, ydata = s['pfs_real'], s['pfs_imag']


            ax = plt.subplot(gs1[ax_ind])
            ax.plot(xdata, ydata, linestyle='', marker='o', color='blue')
            x_nticks = 5
            y_nticks = 5
            tick_scaleoff = 0.8
            mscale = 1.1
            ax.set_xlim(xlims[J])
            ax.set_ylim(ylims[I])
            xlims_ax = tick_scaleoff * xlims[J]
            ylims_ax = tick_scaleoff * ylims[I]
            xlims_ax = numpy.linspace(*xlims_ax, x_nticks)
            ylims_ax = numpy.linspace(*ylims_ax, y_nticks)
            if not I == ncols-1:
                ax.set_xticklabels([])
                ax.set_xticks([])
            else:
                # ax.set_xticks(numpy
                ax.set_xticks(xlims_ax)
                ax.set_xlabel(r"$\Re(E)$")
            if not J == 0:
                ax.set_yticks([])
                ax.set_yticklabels([])
            else:
                ax.set_yticks(ylims_ax)
                ax.set_ylabel(r"$\Im(E)$")

            u = f'{s["phi"]:.3f}'
            u = u.rstrip('0').rstrip('.')
            lab = labels.make_latex_label(x)
            # lab = rf'$\phi = {u}$'
            loc = "upper right"
            text_box = AnchoredText(lab, frameon=True, loc=loc, pad=0)
            text_box.patch.set_boxstyle("square,pad=0.15")
            plt.setp(text_box.patch, linestyle='', facecolor='lightgrey')
            ax.add_artist(text_box)
            ax_ind += 1
        for x, y in S.slicer.get_sets().items():
            I = ax_ind // ncols
            J = ax_ind % ncols
            s = y[0]
            P = PointGenerator(s)
            xdata, ydata = s['pfs_real'], s['pfs_imag']
            ax = plt.subplot(gs1[-1])
            ax.plot(xdata, ydata, linestyle='', marker='o', color='blue')
    f = os.path.join(S.graph_dir, f'spectra.png')
    # fig.tight_layout()
    # fig.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(f)

    scale = 1
    plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)
    ex = S['xkey']
    for x, y in S.slicer.get_sets(exclude=ex).items():
        plt.clf()
        cm = matplotlib.cm.get_cmap('cool')
        ind = 0
        for s in y:
            P = PointGenerator(s)
            xdata, ydata = s['pfs_real'], s['pfs_imag']
            ind += 1
            c = cm(ind / len(y))
            plt.plot(xdata, ydata, linestyle='', marker='o', markersize=2, color=c)
        f = os.path.join(S.graph_dir, f'pfs_{labels.make_ascii_label(x)}.png')
        plt.savefig(f)

def fin_cl_pf_trajectories(S):
    default_x = 6.4
    default_y = 4.8
    do_rotations = S.get('do_rotations', 0)
    scale = 1
    plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)
    exclude_keys = []
    if S['xkeys']:
        exclude_keys = S['xkeys']
    else:
        exclude_keys = [S['xkey']]
    for k in exclude_keys:
        for x, y in S.slicer.get_sets(exclude=k).items():
            print(x)
            plt.clf()
            cm = matplotlib.cm.get_cmap('cool')
            # cm = matplotlib.cm.get_cmap('nipy_spectral')
            ind = 0
            n_plotted = 10000
            for s in y:
                for j in range(s['N']):
                    pfs = numpy.array(s['pfs']) * s['L']
                    m = numpy.max(abs(pfs))
                    lim = 0.2
                    xdata = numpy.full(len(pfs), 0.)
                    ydata = numpy.full(len(pfs), 0.)
                    for i in range(len(pfs)):
                        if i > n_plotted:
                            continue
                        pf = pfs[i]
                        # if abs(pf)/m > lim:
                        #     continue

                        z = pf * numpy.exp(2.j * numpy.pi * j / s['N'])
                        xdata[i] = z.real
                        ydata[i] = z.imag
                    # xdata, ydata = s['pfs_real'], s['pfs_imag']
                    c = cm(ind / len(y))
                    plt.plot(xdata, ydata, linestyle='', marker='o', markersize=1.3, color=c, zorder=10)

                    if not do_rotations:
                        break
                ind += 1

            # plt.colorbar()
            m = matplotlib.cm.ScalarMappable(cmap=cm)
            plt.colorbar(m)
    # fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            f = os.path.join(S.graph_dir, f'pfs_key={k}_{labels.make_ascii_label(x)}.png')
            plt.xlabel(r'$\text{Re}(\varepsilon)$')
            plt.ylabel(r'$\text{Im}(\varepsilon)$')
            add_radial_lines(plt.gca(), y[0], zorder=1)
            plt.savefig(f)
            plt.clf()

def fin_cl_spectrum_test(S):
    # test energies with pfs and ED
    for x, y in S.slicer.get_sets(exclude=[]).items():
        s = y[0]
        default_x = 6.4
        default_y = 4.8
        scale = 1
        plt.clf()
        plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)


        xdata2 = []
        ydata2 = []
        pfs_poly = sympy_eps2(s)
        e0 = 0
        for p in pfs_poly:
            e0 -= p
        spect_poly = make_pf_spectrum2(e0, pfs_poly, s['N'])
        for e in spect_poly:
            xdata2.append(e.real)
            ydata2.append(e.imag)
        plt.plot(xdata2, ydata2, linestyle='', marker='x', markersize=5, color='red')

        xdata = []
        ydata = []
        for e in s['energies']:
            xdata.append(e.real)
            ydata.append(e.imag)
        plt.plot(xdata, ydata, linestyle='', marker='o', markersize=2)

        f = os.path.join(S.graph_dir, f'spectrum_test_{labels.make_ascii_label(x)}.png')
        plt.xlabel(r'$\text{Re}(E)$')
        plt.ylabel(r'$\text{Im}(E)$')
        plt.savefig(f)
        plt.clf()

def remove_small_imag(ydata):
    def f(y):
        if abs(y) < 1E-10:
            return 0
        return y
    return [f(y) for y in ydata]

def add_zero_to_lims(xlims, ylims):
    xlims = numpy.array(xlims)
    ylims = numpy.array(ylims)
    if xlims[0] > 0: xlims[0] = 0
    if xlims[1] < 0: xlims[1] = 0
    if ylims[0] > 0: ylims[0] = 0
    if ylims[1] < 0: ylims[1] = 0
    return xlims, ylims

def get_multiplot_data(s, target):
    xdata, ydata, plotargs = [], [], {}
    if target == 'pfs':
        xdata, ydata = s['pfs_real'], s['pfs_imag']
    if target == 'energies':
        es = s['energies']
        print('Point count:', len(es), s['N'] ** s['L'])
        xdata = [e.real for e in es]
        ydata = [e.imag for e in es]
        plotargs['markersize'] = 0.33333
        plotargs['color'] = 'blue'
    ydata = remove_small_imag(ydata)
    return xdata, ydata, plotargs

def add_radial_lines(ax, s, **plotargs):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    xlim = numpy.array(xlim)
    ylim = numpy.array(ylim)
    for i in range(s['N']):
        z = numpy.exp(i * 2.j*numpy.pi / s['N'])
        u = max(max(abs(xlim)), max(abs(xlim))) * 3.
        xdata = [0, z.real*u]
        ydata = [0, z.imag*u]
        pa = plotargs
        pa.update({
            # 'color' : 'mediumseagreen',
            'color' : 'darkgrey',
            'marker' : '',
            # 'linestyle' : (0, (1, 5)),
            'linestyle' : '--',
            'linewidth' : 1.0,
            })
        ax.plot(xdata, ydata, **pa)

        z = numpy.exp((i+0.5) * 2.j * numpy.pi / s['N'])
        u = max(max(abs(xlim)), max(abs(xlim))) * 3.
        xdata = [0, z.real*u]
        ydata = [0, z.imag*u]
        pa.update({
            # 'color' : 'indianred',
            })
        ax.plot(xdata, ydata, **pa)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def fin_cl1_multiplot(S):
    ts = S['multiplot_targets']
    if ts == []:
        ts = ['pfs', 'energies']

    for t in ts:
        ax_ind = 0
        ncols = S.get('multiplot_ncols', 3)
        nrows = len(S.systems) // ncols
        if ncols * nrows < len(S.systems):
            nrows += 1
        gs1 = gridspec.GridSpec(nrows, ncols)
        gs1.update(wspace=0, hspace=0)
        all_axs = []
        default_x = 6.4
        default_y = 4.8
        scale = 0.5
        figsize=(default_x*scale*ncols, default_y*scale*nrows)
        fig, axs = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows == 1:
            all_axs = axs
        else:
            for z in axs: all_axs += [z0 for z0 in z]

        for x, y in S.slicer.get_sets().items():
            for s in y:
                I = ax_ind // ncols
                J = ax_ind % ncols
                ax = all_axs[ax_ind]

                # plot
                xdata, ydata, plotargs2 = get_multiplot_data(s, t)
                plotargs = dict(marker='o', color='blue', linestyle='')
                plotargs.update(plotargs2)
                ax.plot(xdata, ydata, **plotargs)

                # axis labels
                if I == nrows-1:
                    ax.set_xlabel(r"$\Re(\varepsilon)$")
                if J == 0:
                    ax.set_ylabel(r"$\Im(\varepsilon)$")

                # labels
                lab = labels.make_latex_label(x)
                loc = "upper right"
                if s['label_loc']:
                    loc = ' '.join(s['label_loc'])
                text_box = AnchoredText(lab, frameon=True, loc=loc, pad=0)
                text_box.patch.set_boxstyle("square,pad=0.15")
                plt.setp(text_box.patch, linestyle='', facecolor='grey', alpha=0.3)
                ax.add_artist(text_box)

                # include zero
                xlim, ylim = add_zero_to_lims(ax.get_xlim(), ax.get_ylim())
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)

                # radial lines
                if t == 'pfs':
                    add_radial_lines(ax, s)

                ax_ind += 1
        fig.tight_layout()
        f = os.path.join(S.graph_dir, f'{t}.png')
        plt.savefig(f)
        f = os.path.join(S.pgf_dir, f'{t}.pgf')
        plt.savefig(f)

def fin_traject_multiplot(S):
    print(S['multiplot_targets'])
    default_x = 6.4
    default_y = 4.8
    do_rotations = S.get('do_rotations', 0)
    scale = 1
    plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)
    exclude_keys = []
    if S['xkeys']:
        exclude_keys = S['xkeys']
    else:
        exclude_keys = [S['xkey']]
    for k in exclude_keys:
        ax_ind = 0
        ncols = S['multiplot_ncols']
        slice = S.slicer.get_sets(exclude=k)
        nrows = len(slice) // ncols + len(slice) % ncols
        plt.clf()
        all_axs = []
        default_x = 6.4
        default_y = 4.8
        scale = 0.5
        figsize=(default_x*scale*ncols, default_y*scale*nrows)
        fig, axs = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows == 1:
            all_axs = axs
        else:
            for z in axs: all_axs += [z0 for z0 in z]

        for x, y in S.slicer.get_sets(exclude=k).items():
            I = ax_ind // ncols
            J = ax_ind % ncols
            ax = all_axs[ax_ind]

            # plot
            cm = matplotlib.cm.get_cmap('cool')
            ind = 0
            for s in y:

                for i in range(s['N']):
                    P = PointGenerator(s)
                    pfs = s['pfs']
                    xdata = []
                    ydata = []
                    for pf in pfs:
                        z = pf * numpy.exp(2.j * numpy.pi * i / s['N']) / s['pf_scale']
                        a = z.real
                        b = z.imag
                        xdata.append(a)
                        ydata.append(b)
                    c = cm(ind / len(y))
                    plotargs = dict(linestyle='', marker='o', markersize=2, color=c, zorder=3)
                    if rel_diff(s[k], 0.5) < 1E-10:
                        plotargs.update({
                            'color' : 'black',
                            'markersize' : 3,
                            'zorder' : 4,
                            })
                    ax.plot(xdata, ydata, **plotargs, rasterized=True)
                    # ax.set_rasterized(True)

                    if not do_rotations:
                        break
                ind += 1

            # axis labels
            if I == nrows-1:
                ax.set_xlabel(r"$\Re(\varepsilon)$")
            if J == 0:
                ax.set_ylabel(r"$\Im(\varepsilon)$")

            # labels
            lab = labels.make_latex_label(x)
            loc = "upper right"
            if s['label_loc']:
                loc = ' '.join(s['label_loc'])
            text_box = AnchoredText(lab, frameon=True, loc=loc, pad=0)
            text_box.patch.set_boxstyle("square,pad=0.15")
            plt.setp(text_box.patch, linestyle='', facecolor='lightgray', alpha=1.0, zorder=8)
            ax.add_artist(text_box)

            # include zero
            xlim, ylim = add_zero_to_lims(ax.get_xlim(), ax.get_ylim())
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            # radial lines
            add_radial_lines(ax, s)

            ax_ind += 1

        t = 'test'
        fig.tight_layout()
        f = os.path.join(S.graph_dir, f'{t}.png')
        plt.savefig(f)
        f = os.path.join(S.pgf_dir, f'{t}.pgf')
        plt.savefig(f)

    # fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)

from lmfit import minimize, Parameters


def residual(p, x, data):
    e0 = p['e0_real'] + p['e0_imag'] * 1.j
    A = p['A_real'] + p['A_imag'] * 1.j
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
    xdata = data['xdata']
    ydata = data['ydata']
    out = minimize(residual, params, args=(xdata, ydata))
    p = out.params
    e0 = p['e0_real'] + p['e0_imag'] * 1.j
    A = p['A_real'] + p['A_imag'] * 1.j
    return e0, A, out

def gse_test(S):
    d = DataSet(S, xkey='L', ykey='e_0')
    def fitf(L, einf, Lcoef): return -einf - Lcoef/L

    # d.plot(linestyle='')
    for key, data in d.data.items():
        s = S.slicer.find(key)[0]
        print(s)
        N = s['N']
        # s.init_model()
        v = s['e_inf']
        print(f'einf_exact = {v}')
        e0, A, out = lmtest(data)
        print(f'einf_interpolated = {e0}')
        print(f'e0 = {s["e_0"]}')
        # print(s['energies'][:5])

def cl_basic_plots(S):
    for x, y in S.slicer.get_sets(exclude=[]).items():
        s = y[0]
        default_x = 6.4
        default_y = 4.8
        scale = 1
        plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)

        xdata = []
        ydata = []
        for e in s['energies']:
            xdata.append(e.real)
            ydata.append(e.imag)
        plt.plot(xdata, ydata, linestyle='', marker='o', markersize=2)

        u = s['e_inf']
        plt.plot([u.real], [u.imag], linestyle='', marker='x', c='red', markersize=6)
        u = s['energies'][0]
        plt.plot([u.real], [u.imag], linestyle='', marker='x', c='green', markersize=6)


        f = os.path.join(S.graph_dir, f'spectrum_test_{labels.make_ascii_label(x)}.png')
        plt.xlabel(r'$\text{Re}(E)$')
        plt.ylabel(r'$\text{Im}(E)$')

        plt.savefig(f)
        plt.clf()


def gse1(S):
    def fitf(L, einf, Lcoef): return -einf - Lcoef/L
    d_real = DataSet(S, xkey='L', ykey='e_real')
    d_imag = DataSet(S, xkey='L', ykey='e_imag')
    d = DataSet(S, xkey='L', ykey='e_0')

    cell_text = numpy.empty(shape=(len(d_real.data), 8), dtype=object)
    col_labels = [
        r'$\gamma$',
        r'$\phi$',
        r'$\Re{(e_\infty^{\text{fit}})}$',
        r'$\Im{(e_\infty^{\text{fit}})}$',
        r'$\Re{(e_\infty^{\text{exact}})}$',
        r'$\Im{(e_\infty^{\text{exact}})}$',
        r'Rel.\ Diff.',
        r'$\chi^2_\nu$',
        ]
    col_colours = ['lightgray'] * len(col_labels)
    ind = 0
    for key, data_real in d_real.data.items():
        keydict = {}
        for a, b in key:
            keydict[a] = b

        data_imag = d_imag[key]
        phi, gamma = keydict['phi'], keydict['gamma']
        e0_fit, A, out = lmtest(d[key])
        u = ind
        s = S.slicer.find(key)[0]
        einf = s['e_inf']
        cell_text[u][0] = f'{gamma:1.2f}'
        cell_text[u][1] = f'{phi:1.2f}'
        F1 = numpy.format_float_positional
        Fargs = dict(precision=6)
        cell_text[u][2] = f'{F1(e0_fit.real, **Fargs)}'
        cell_text[u][3] = f'{F1(e0_fit.imag, **Fargs)}'
        cell_text[u][4] = f'{F1(einf.real, **Fargs)}'
        cell_text[u][5] = f'{F1(einf.imag, **Fargs)}'
        Delta = rel_diff(e0_fit, einf)
        F1 = numpy.format_float_scientific
        Fargs = dict(precision=4)
        cell_text[u][6] = f'{F1(Delta, **Fargs)}'
        cell_text[u][7] = f'{F1(out.redchi, **Fargs)}'
        ind += 1
    fig, ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    table = ax.table(cell_text, loc='center', colLabels= col_labels, colColours=col_colours)
    table.auto_set_font_size(False)
    fig.tight_layout()
    print(1)
    f = os.path.join(S.graph_dir, f'gse_table.pdf')
    plt.savefig(f)
    f = os.path.join(S.pgf_dir, f'gse_table.pgf')
    plt.savefig(f)
    print(1)

def gap1(S):
    d_real = DataSet(S, xkey='L', ykey='gap_real')
    d_imag = DataSet(S, xkey='L', ykey='gap_imag')
    d = DataSet(S, xkey='L', ykey='gap')

    cell_text = numpy.empty(shape=(len(d_real.data), 8), dtype=object)
    col_labels = [
        r'$\gamma$',
        r'$\phi$',
        r'$\Re{(\Delta^{\text{fit}})}$',
        r'$\Im{(\Delta^{\text{fit}})}$',
        r'$\Re{(\Delta^{\text{exact}})}$',
        r'$\Im{(\Delta^{\text{exact}})}$',
        r'Rel.\ Diff.',
        r'$\chi^2_\nu$',
        ]
    col_colours = ['lightgray'] * len(col_labels)
    ind = 0
    for key, data_real in d_real.data.items():
        keydict = {}
        for a, b in key:
            keydict[a] = b

        data_imag = d_imag[key]
        phi, gamma = keydict['phi'], keydict['gamma']
        e0_fit, A, out = lmtest(d[key])
        u = ind
        s = S.slicer.find(key)[0]
        einf = s['gap_inf']
        cell_text[u][0] = f'{gamma:1.2f}'
        cell_text[u][1] = f'{phi:1.2f}'
        F1 = numpy.format_float_positional
        Fargs = dict(precision=6)
        cell_text[u][2] = f'{F1(e0_fit.real, **Fargs)}'
        cell_text[u][3] = f'{F1(e0_fit.imag, **Fargs)}'
        cell_text[u][4] = f'{F1(einf.real, **Fargs)}'
        cell_text[u][5] = f'{F1(einf.imag, **Fargs)}'
        Delta = rel_diff(e0_fit, einf)
        F1 = numpy.format_float_scientific
        Fargs = dict(precision=4)
        cell_text[u][6] = f'{F1(Delta, **Fargs)}'
        cell_text[u][7] = f'{F1(out.redchi, **Fargs)}'
        ind += 1
    fig, ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    table = ax.table(cell_text, loc='center', colLabels= col_labels, colColours=col_colours)
    table.auto_set_font_size(False)
    fig.tight_layout()
    print(2)
    f = os.path.join(S.graph_dir, f'gap_table.pdf')
    plt.savefig(f)
    print(2)
    f = os.path.join(S.pgf_dir, f'gap_table.pgf')
    plt.savefig(f)
    print(2)


def cl_spectrum_trajectory(S):
    default_x = 6.4
    default_y = 4.8
    do_rotations = S.get('do_rotations', 0)
    scale = 1
    plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)
    exclude_keys = []
    if S['xkeys']:
        exclude_keys = S['xkeys']
    else:
        exclude_keys = [S['xkey']]
    for k in exclude_keys:
        for x, y in S.slicer.get_sets(exclude=k).items():
            print(x)
            plt.clf()
            cm = matplotlib.cm.get_cmap('cool')
            # cm = matplotlib.cm.get_cmap('nipy_spectral')
            ind = 0
            n_plotted = 1E10
            for s in y:
                for j in range(s['N']):
                    pfs = numpy.array(s['pfs'])
                    # pfs = numpy.array(s['energies']) * s['L']
                    m = numpy.max(abs(pfs))
                    lim = 0.2
                    xdata = numpy.full(len(pfs), 0.)
                    ydata = numpy.full(len(pfs), 0.)
                    for i in range(len(pfs)):
                        if i > n_plotted:
                            continue
                        pf = pfs[i]
                        # if abs(pf)/m > lim:
                        #     continue

                        z = pf * numpy.exp(2.j * numpy.pi * j / s['N'])
                        xdata[i] = z.real
                        ydata[i] = z.imag
                    # xdata, ydata = s['pfs_real'], s['pfs_imag']
                    c = cm(ind / len(y))
                    plt.plot(xdata, ydata, linestyle='', marker='o', markersize=1, color=c, zorder=10)

                    if not do_rotations:
                        break
                ind += 1

            # plt.colorbar()
            m = matplotlib.cm.ScalarMappable(cmap=cm)
            plt.colorbar(m)
    # fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            f = os.path.join(S.graph_dir, f'energies_key={k}_{labels.make_ascii_label(x)}.png')
            plt.xlabel(r'$\text{Re}(\varepsilon)$')
            plt.ylabel(r'$\text{Im}(\varepsilon)$')
            add_radial_lines(plt.gca(), y[0], zorder=1)
            plt.savefig(f)
            plt.clf()


def spectrum_test_new(S):
    1


def fin_traject_multiplot_k(S):
    print(S['multiplot_targets'])
    default_x = 6.4
    default_y = 4.8
    do_rotations = S.get('do_rotations', 0)
    scale = 1
    plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)
    exclude_keys = []
    if S['xkeys']:
        exclude_keys = S['xkeys']
    else:
        exclude_keys = [S['xkey']]
    for k in exclude_keys:
        ax_ind = 0
        ncols = S['multiplot_ncols']
        slice = S.slicer.get_sets(exclude=k)
        nrows = len(slice) // ncols + len(slice) % ncols
        plt.clf()
        all_axs = []
        default_x = 6.4
        default_y = 4.8
        scale = 0.5
        figsize=(default_x*scale*ncols, default_y*scale*nrows)
        fig, axs = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows == 1:
            all_axs = axs
        else:
            for z in axs: all_axs += [z0 for z0 in z]

        for x, y in S.slicer.get_sets(exclude=k).items():
            I = ax_ind // ncols
            J = ax_ind % ncols
            ax = all_axs[ax_ind]

            # plot
            cm = matplotlib.cm.get_cmap('cool')
            ind = 0
            for s in y:

                for i in range(s['N']):
                    P = PointGenerator(s)
                    pfs = s['pfs']
                    xdata = []
                    ydata = []
                    for pf in pfs:
                        z = pf * numpy.exp(2.j * numpy.pi * i / s['N']) / s['pf_scale']
                        a = z.real
                        b = z.imag
                        xdata.append(a)
                        ydata.append(b)
                    c = cm(ind / len(y))
                    plotargs = dict(linestyle='', marker='o', markersize=2, color=c, zorder=3)
                    if rel_diff(s[k], 0.5) < 1E-10:
                        plotargs.update({
                            'color' : 'black',
                            'markersize' : 3,
                            'zorder' : 4,
                            })
                    ax.plot(xdata, ydata, **plotargs, rasterized=True)
                    # ax.set_rasterized(True)

                    if not do_rotations:
                        break
                ind += 1

            # axis labels
            if I == nrows-1:
                ax.set_xlabel(r"$\Re(\varepsilon)$")
            if J == 0:
                ax.set_ylabel(r"$\Im(\varepsilon)$")

            # labels
            lab = labels.make_latex_label(x)
            loc = "upper right"
            if s['label_loc']:
                loc = ' '.join(s['label_loc'])
            text_box = AnchoredText(lab, frameon=True, loc=loc, pad=0)
            text_box.patch.set_boxstyle("square,pad=0.15")
            plt.setp(text_box.patch, linestyle='', facecolor='lightgray', alpha=1.0, zorder=8)
            ax.add_artist(text_box)

            # include zero
            xlim, ylim = add_zero_to_lims(ax.get_xlim(), ax.get_ylim())
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            # radial lines
            add_radial_lines(ax, s)

            ax_ind += 1

        t = 'test'
        fig.tight_layout()
        f = os.path.join(S.graph_dir, f'{t}.png')
        plt.savefig(f)
        f = os.path.join(S.pgf_dir, f'{t}.pgf')
        plt.savefig(f)

    # fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)


def fin_cl_pf_trajectories_k(S):
    default_x = 6.4
    default_y = 4.8
    do_rotations = S.get('do_rotations', 0)
    scale = 1
    plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)
    exclude_keys = []
    if S['xkeys']:
        exclude_keys = S['xkeys']
    else:
        exclude_keys = [S['xkey']]
    for k in exclude_keys:
        for x, y in S.slicer.get_sets(exclude=k).items():
            print(x)
            plt.clf()
            cm = matplotlib.cm.get_cmap('cool')
            # cm = matplotlib.cm.get_cmap('nipy_spectral')
            ind = 0
            n_plotted = 10000
            for s in y:
                pfs = numpy.array(s['pfs'])
                # m = numpy.max(abs(pfs))
                lim = 0.2
                xdata = numpy.full(len(pfs), 0.)
                ydata = numpy.full(len(pfs), 0.)
                for i in range(len(pfs)):
                    if i > n_plotted:
                        continue
                    pf = pfs[i]
                    # if abs(pf)/m > lim:
                    #     continue

                    # z = pf * numpy.exp(2.j * numpy.pi * j / s['N'])
                    z = eps_to_k(pf, s['lambda'], s['N'])
                    q = eps_k(z, s['lambda'], s['N'])
                    d = abs(abs(q) - abs(pf))
                    # print(pf, z, q, d)
                    print(s['lambda'], abs(z))
                    xdata[i] = z.real
                    ydata[i] = z.imag
                    c = cm(ind / len(y))
                    plt.plot(xdata, ydata, linestyle='', marker='o', markersize=1.3, color=c, zorder=10)

                    # if not do_rotations:
                    #     break
                ind += 1

            # plt.colorbar()
            m = matplotlib.cm.ScalarMappable(cmap=cm)
            plt.colorbar(m)
    # fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            f = os.path.join(S.graph_dir, f'ks_key={k}_{labels.make_ascii_label(x)}.png')
            plt.xlabel(r'$\text{Re}(\varepsilon)$')
            plt.ylabel(r'$\text{Im}(\varepsilon)$')
            add_radial_lines(plt.gca(), y[0], zorder=1)
            plt.savefig(f)
            plt.clf()


def trivial_test(S):
    default_x = 6.4
    default_y = 4.8
    do_rotations = S.get('do_rotations', 0)
    scale = 1
    plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)
    exclude_keys = []
    if S['xkeys']:
        exclude_keys = S['xkeys']
    else:
        exclude_keys = [S['xkey']]
    for k in exclude_keys:
        for x, y in S.slicer.get_sets(exclude=k).items():

            plt.clf()
            cm = matplotlib.cm.get_cmap('cool')
            # cm = matplotlib.cm.get_cmap('nipy_spectral')
            ind = 0
            n_plotted = 10000
            s0 = y[0]
            for ind in range(len(s0['pfs'])):
                xdata = []
                ydata = []
                ydatai = []
                for s in y:
                    p = s['pfs'][ind]
                    p = eps_to_k(p, s['lambda'], s['N'])

                    # xdata.append(s['lambda'])
                    # ydata.append(numpy.real(ys))
                    # ydatai.append(numpy.imag(ys))
                    xdata.append(s['lambda'])
                    ydata.append(numpy.real(p))
                    ydatai.append(numpy.imag(p))

                    # z = eps_to_k(z, s['lambda'], s['N'])
                c = cm(ind / len(s0['pfs']))
                plt.plot(xdata, ydata, linestyle='', marker='o', markersize=1.3, c=c)
                plt.plot(xdata, ydatai, linestyle='', marker='x', c=c)

                    # if not do_rotations:
                    #     break

            # plt.colorbar()
    # fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            f = f'ttest_{x}.png'
            f = os.path.join(S.graph_dir, f)
            print(f)
            plt.savefig(f)
            plt.clf()
