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
from labels import *

def plot_spectrum(S):
    for x, y in S.slicer.get_sets().items():
        s = y[0]
        plt.clf()
        xdata = []
        ydata = []
        es = s['energies']
        N = s['N']
        L = s['L']
        print(len(es), N**L)
        for e in s['energies']:
            # print(e)
            xdata.append(e.real/L)
            ydata.append(e.imag/L)
        plt.plot(xdata, ydata, linestyle='', marker='o', markersize=2)
        e = s['e_0']
        plt.plot([e.real], [e.imag], linestyle='', marker='x', markersize=5, c='red')
        e = s['e_inf2']
        plt.plot([e.real], [e.imag], linestyle='', marker='+', markersize=5, c='green')

        f = make_ascii_label(x)
        f = f'{f}.png'
        f = os.path.join(S.graph_dir, f)
        plt.savefig(f)
    1
