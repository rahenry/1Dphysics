#!/usr/bin/env python3


from fp_polynomial import *
from rah_utility import rel_diff
import graph_settings
from scipy.optimize import fsolve
import matplotlib.transforms as transforms
from matplotlib.offsetbox import AnchoredText

def poly_to_sympy(poly):
    u = sympy.Symbol('y')
    npoly = 0
    j = 0
    for x in poly:
        npoly += x * (u**j)
        j += 1
    npoly = sympy.Poly(npoly)
    return npoly

def converter(x):
    x = x+0.001
    sgn = 1
    if x.real < 0: sgn = -1
    x = abs(x) * sgn
    if rel_diff(x, int(x)) < 1E-2:
        return int(x)

do_convert = 0
do_rescale = 0
def sympy_roots(s, maxsteps=1000, n=1000):
    lambdas = s['lambdas']
    l = min(lambdas)
    lambdas = [x / l for x in lambdas]
    if do_rescale:
        a = max(lambdas)
        b = min(lambdas)
        avg = numpy.sqrt(a*b)
        lambdas = [x/avg for x in lambdas]
    if do_convert:
        lambdas = [converter(x) for x in lambdas]

    roots_old, poly_base = find_roots_poly(s['M'], s['p'], lambdas)
    poly = poly_to_sympy(poly_base)
    roots = []
    try:
        roots = poly.nroots(maxsteps=maxsteps, n=n)
    except:
        roots = []
    return roots, poly, l, poly_base

def sympy_eps(s):
    ms = 500
    pp = 500
    N = s['N']
    roots, poly, l, poly_base = sympy_roots(s, ms, pp)
    pfs = [numpy.complex64(r/l)**(-1./N) for r in roots]
    pfs = [p/s['L'] for p in pfs]
    roots_accepted = []
    roots_rejected = []
    for r in roots:
        q = poly(r)
        if abs(q) > 1E-3:
            roots_rejected.append(r)
        else:
            roots_accepted.append(r)
    if len(roots_rejected) > 0:
        print(f'{s} has {len(roots_rejected)} bad solutions')
        for x in roots_rejected:
            print(x, poly(x))
    return pfs


def sympy_eps2(s):
    ms = 2000
    pp = 500
    N = s['N']
    roots, poly, l, poly_base = sympy_roots(s, ms, pp)
    if len(roots) == 0:
        print(roots)
        print(s)
        exit()
    pfs = [numpy.complex64(r/l)**(-1./N) for r in roots]
    pfs = [p/s['L'] for p in pfs]
    roots_accepted = []
    roots_rejected = []
    for r in roots:
        q = poly(r)
        if abs(q) > 1E-3:
            roots_rejected.append(r)
        else:
            roots_accepted.append(r)
    if len(roots_rejected) > 0:

        if len(roots_rejected) == 1 and abs(poly(roots_rejected[0])) > 1E20:
            1
        else:
            print(f'{s} has {len(roots_rejected)} bads, with vals:')
            for x in roots_rejected:
                print(poly(x))
    return pfs


def polytest2(S):
    base_dir = S.graph_dir
    fail_data = {}
    for x, y in S.slicer.get_sets(exclude=['L', 'Mplus']).items():
        fail_data[x] = []
        xtick_coords = []
        for s in y:
            s.xdata = []
            s.ydata = []
            s.xdata2 = []
            s.ydata2 = []
            N = s['N']
            N = s['N']
            M = s['M']
            xtick_coords.append(M)
            L = s['L']
            roots, poly, l, s.poly_base = sympy_roots(s, S['poly_maxsteps'][0], S['poly_precision'][0])
            roots_accepted = []
            roots_rejected = []
            for r in roots:
                q = poly(r)
                if abs(q) > 1E-3:
                    roots_rejected.append(r)
                else:
                    roots_accepted.append(r)
            roots = roots_accepted
            roots = [r/l for r in roots]
            s.roots = roots
            s.poly = poly
            s.missing = len(roots_rejected)
            fail_data[x].append(len(roots_rejected))

            pfs = [numpy.complex64(r)**(-1./N) for r in roots]
            s.pfs_poly = pfs
            for pf in pfs:
                s.xdata.append(M)
                s.ydata.append(abs(pf))

            if s['M'] % 2 == 1 and s['model_type'] == 'fp':
                pfs2 = [p*L for p in s.pfs]
                s.xdata2 = [s['M'] for pf in pfs2]
                s.ydata2 = [abs(pf) for pf in pfs2]

            # ax.set_xlabel("M\n[no. failed roots]")

    # try for more roots

    # individual plots
    labs = []
    for i in range(len(xtick_coords)):
        lab = f'{xtick_coords[i]}\n[{fail_data[x][i]}]'
        lab = f'{xtick_coords[i]}'
        labs.append(lab)
    for x, y in S.slicer.get_sets(exclude=['L', 'Mplus']).items():
        plt.close()
        for s in y:
            if len(s.xdata2) > 0:
                plt.plot(s.xdata2, s.ydata2, linestyle='', marker='x', markersize=10.0, color='red')
            plt.plot(s.xdata, s.ydata, linestyle='', marker='o', markersize=2.5, color='blue')
        plt.xticks(xtick_coords, labs)
        # plt.ylabel("$\varepsilon$")
        plt.ylabel(r"$\varepsilon$", rotation=0)
        plt.xlabel("M\n[no. missing roots]")
        plt.ylim(bottom=0)

        f = os.path.join(base_dir, f'{x}_polytest.png')
        plt.tight_layout()
        plt.savefig(f)
        f = f.replace('.png', '.pgf')
        plt.savefig(f)

    # composite plot

    ncols = S['columns'][0]
    sets = S.slicer.get_sets(exclude=['L', 'Mplus']).items()
    print(len(sets), ncols)
    nrows = len(sets) // ncols
    print(ncols, nrows)
    plt.close()
    aspect = 0.15
    # plt.figure(figsize=(20, 35), dpi=300)
    a4_size = numpy.array([8.25, 11.75*0.75])
    scale = 0.80
    figsize = scale * a4_size
    fig, axs = plt.subplots(nrows, ncols, sharex=True, sharey=False, figsize=figsize, dpi=200)
    fig, axs = plt.subplots(nrows, ncols, sharex=True, sharey=False, figsize=figsize)
    # fig, axs = plt.subplots(nrows, ncols)
    ax_ind = 0
    for x, y in sets:
        i = ax_ind // ncols
        j = ax_ind % ncols
        ax = None
        print(x, i, j)
        if ncols == 1:
            ax = axs[i]
        else:
            ax = axs[i][j]
        for s in y:
            pass
        ax_ind += 1

        plot_zero = False
        ylims = [100, 0.1]
        for s in y:
            yvals = s.ydata + s.ydata2
            m = max(yvals)
            if m > ylims[1]:
                ylims[1] = m
            m = min(yvals)
            if m < ylims[0]:
                ylims[0] = m
        d = ylims[1]-ylims[0]
        ylims = [ylims[0] - d*0.08, ylims[1] +d*0.08]
        ax.set_ylim(*ylims)
        if ylims[0] < 0:
            plot_zero = True
            ylims[0] = 0
        else:
            ylims[0] += 0.10*d
        yticks = numpy.linspace(ylims[0], ylims[1], 3, endpoint=False)

        ylabs = [f'{x:1.2f}' for x in yticks]
        if abs(yticks[0] < 1E-10):
            ylabs[0] = '0'

        if plot_zero:
            xdata = ax.get_xlim()
            ydata = [0, 0]
            ax.axhline(color='lightgray', linewidth=1, linestyle='dashed')
        ax.set_yticks(yticks, ylabs)

        for s in y:
            ax.plot(s.xdata, s.ydata, linestyle='', marker='_', markersize=7, color='blue', linewidth=1)
            # ax.plot(s.xdata, s.ydata, linestyle='', marker='o', markersize=2.5, color='blue', linewidth=1)
            if len(s.xdata2) > 0:
                ax.plot(s.xdata2, s.ydata2, linestyle='', marker='o', markersize=2.5, color='red')
                print(s.ydata2)
            # ax.plot(s.xdata, s.ydata, linestyle='', marker='_', markersize=10, color='green')
            # if s.missing > 0:
            if True:
                trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                ax.text(s.xdata[0]+0.5, 0.03, f'[{s.missing}]', transform=trans, color='red', size='xx-small', ha='center', va='bottom')

        # ax.set_yticks(yticks)
        ax.set_box_aspect()
        ax.set_xticks(xtick_coords, labs)
        ax.set_xticks(xtick_coords)
        ax.set_ylabel(r"$\varepsilon$", fontsize='large', rotation=0)
        label = ''
        for l in S['labels']:
            z = s[l]
            label += rf"$\{l} = {z:1.2f}$, "
        label = label.strip(r' ')
        label = label.strip(r',')
        # props = dict(boxstyle='square', facecolor='lightgrey', alpha=1.00, linestyle='', pad=0)
        # ax.text(1, 1, label, transform=ax.transAxes, ha='right', va='center', bbox=props)
        # ax.text(1, 0.99, label, transform=ax.transAxes, ha='right', va='top', bbox=props)
        # text_box = AnchoredText(label, frameon=True, loc='upper right', pad=0)
        loc = S['text_loc'][0].replace('_', ' ')
        text_box = AnchoredText(label, frameon=True, loc=loc, pad=0)
        # text_box.patch.set_boxstyle("facecolor='lightgrey',linestyle=''")
        text_box.patch.set_boxstyle("square,pad=0.15")
        plt.setp(text_box.patch, linestyle='', facecolor='lightgrey')
        # ax.setp(text_box.patch, facecolor='white', alpha=0.5)
        ax.add_artist(text_box)

    f = os.path.join(base_dir, f'polytest.png')
    # plt.show()
    ax.set_xlabel("M")
    plt.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(f)
    f = f.replace('.png', '.pgf')
    plt.savefig(f)

if __name__ == '__main__':
    f = os.path.join('graphs', 'polytest3')
    silent_remove(f)
    mkdir(f)
