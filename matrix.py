#!/usr/bin/env python3

import numpy
from tenpy.linalg.np_conserved import tensordot, inner, norm, outer
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as mtransforms
import string
import matplotlib.gridspec as gridspec

def sort1(E):
    n = len(E)
    res = []
    inds = numpy.argsort(E)
    res.append(inds[0])
    available = set([inds[x] for x in range(1,n)])
    while len(available) > 0:
        d0 = 100000
        e = E[res[-1]]
        j = None
        for i in available:
            d = abs(E[i] - e)
            if d < d0:
                d0 = d
                j = i
        res.append(j)
        available.remove(j)

    # for i in range(n):
    return res



def analyse_eigs(ER, EL, VR, VL):
    n = len(ER)
    # indsR = numpy.argsort([abs(x) for x in ER])
    indsR = numpy.argsort([abs(x) for x in ER])
    indsR = sort1(ER)
    # indsL = []
    # errors = []
    # for j in range(n):
    #     eL = EL[j]
    #     d0 = 1E50
    #     k = None
    #     for i in range(n):
    #         eR = ER[i]
    #         d = abs(eL-eR)
    #         if d < d0:
    #             d0 = d
    #             k = i
    #     indsL.append(k)
    #     errors.append(d0)
    indsL = numpy.full((n), 0)
    errors = numpy.full((n), 0)
    for i in range(n):
        I = indsR[i]
        eR = ER[I]
        d0 = 1E50
        k = None
        for j in range(n):
            eL = EL[j]
            d = abs(eL-eR)
            if d < d0:
                d0 = d
                k = j
        indsL[i] = k
        errors[i] = d0

    for i in range(n):
        eL = EL[indsL[i]]
        eR = ER[indsR[i]]



    fields = ['overlapLR', 'overlapLRnc', 'overlapRR', 'overlapRRnc', 'deltaLR', 'deltaRR']
    min_overlap = 1000
    min_overlap2 = 1000
    res = {}
    def inf(*args, **kwargs):
        return inner(*args, **kwargs, axes='range')
    for f in fields:
        res[f] = numpy.full((n, n), 0.)
    for I in range(n):
        i = indsR[I]
        for J in range(n):
            jR = indsR[J]
            j = indsL[J]
            eR = ER[i]
            vR = VR[:,i]
            eL = EL[j]
            vL = VL[:,j]
            eRj = ER[jR]
            vRj = VR[:,jR]
            o = abs(inf(vL, vR, do_conj=1))
            res['overlapLR'][I][J] = o
            o = abs(inf(vL, vR, do_conj=0))
            if I == J:
                if o < min_overlap:
                    min_overlap = o

            res['overlapLRnc'][I][J] = o
            o = abs(inf(vRj, vR, do_conj=1))

            if I == J:
                if o < min_overlap2:
                    min_overlap2 = o

            res['overlapRR'][I][J] = o

            o = abs(inf(vRj, vR, do_conj=0))
            res['overlapRRnc'][I][J] = o
            d = abs(eR-eL)
            res['deltaLR'][I][J] = d
            d = abs(eR-eRj)
            res['deltaRR'][I][J] = d

    fields = ['overlapLR', 'overlapLRnc', 'overlapRR', 'overlapRRnc', 'deltaLR', 'deltaRR']
    for f in fields:
        res[f'{f}_diag'] = numpy.full((n), 0.)
        for I in range(n):
            res[f'{f}_diag'][I] = res[f][I][I]

    res['min_overlap'] = min_overlap
    res['min_overlap2'] = min_overlap2
    res['overlap0'] = res['overlapLRnc'][0][0]
    return res

def plot_eig_info(data, basename):
    for key, vals in data.items():
        continue
        f = f'{basename}_{key}.png'
        plt.clf()
        n = vals.shape[0]
        ex = [0.5, n+0.5, 0.5,n+0.5]
        p = plt.imshow(vals, origin='lower', vmin=0, extent=ex)
        plt.colorbar(p)
        # plt.xlabel(r'$i$')
        # plt.ylabel(r'$j$')
        # ticks = [x+1 for x in range(n)]
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        # plt.xticks(ticks)
        # plt.yticks(ticks)
        # plt.axis('off')
        plt.savefig(f)

    keys = ['deltaLR', 'overlapLRnc', 'overlapRR']
    plt.clf()
    fig, axs = plt.subplots(1, 3)
    labels = ['(a)', '(b)', '(c)']
    f = f'{basename}.png'
    for i, k in enumerate(keys):
        ax = axs[i]
        vals = data[k]
        p = ax.imshow(vals, origin='lower', vmin=0)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.15)
        fig.colorbar(p, cax=cax, orientation='vertical')
        # ax.colorbar(p)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())

        trans = mtransforms.ScaledTranslation(-0/72, 7/72, fig.dpi_scale_trans)
        label = labels[i]
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
                fontsize='medium', va='bottom', fontfamily='serif')

    fig.tight_layout()
    plt.savefig(f)


def plot_eig_sets(datasets, basename):
    keys = ['deltaLR', 'overlapLRnc', 'overlapRR']
    plt.clf()
    n = len(datasets)
    scale = 1.5
    fig = plt.figure(figsize = (3*scale,n*scale))
    gs1 = gridspec.GridSpec(n, 3)
    gs1.update(wspace=0.125, hspace=0.05) # set the spacing between axes.
    # fig, axs = plt.subplots(len(datasets), 3)
    labs = string.ascii_lowercase
    ind = 0
    ind_gs = 0
    for j, data in enumerate(datasets):
        for i, k in enumerate(keys):
            # ax = axs[j,i]
            ax = plt.subplot(gs1[ind_gs])
            ind_gs += 1
            vals = data[k]
            p = ax.imshow(vals, origin='lower', vmin=0)

            # divider = make_axes_locatable(ax)
            # cax = divider.append_axes('right', size='5%', pad=0.15)
            # fig.colorbar(p, cax=cax, orientation='vertical')

            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())

            trans = mtransforms.ScaledTranslation(-0/72, 7/72, fig.dpi_scale_trans)
            # label = labels[i]
            if (i % 3) == 0:
                label = f'({labs[ind]})'
                ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
                        fontsize='medium', va='bottom', fontfamily='serif')
                ind += 1

    # plt.subplots_adjust(wspace=0, hspace=0)
    f = f'{basename}.png'
    plt.savefig(f)
    # f = f'{basename}.pgf'
    # plt.savefig(f)
