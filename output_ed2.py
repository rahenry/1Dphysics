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

def plot_ed2(S):
    for s in S.systems:
        cm = matplotlib.cm.get_cmap('hsv')
        i = 0
        plt.clf()
        default_x = 6.4
        default_y = 4.8
        scale = 2
        plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)
        for (q, k), B in s.ed_system.blocks.items():
            c = cm(i / len(s.ed_system.blocks))
            i += 1
            xdata = [x.real for x in B.energies]
            ydata = [x.imag for x in B.energies]
            plt.plot(xdata, ydata, marker='o', markersize=2.5, color=c, linestyle='', alpha=0.3)


        f = f'{s}_spectrum.png'
        f = os.path.join(S.graph_dir, f)
        plt.savefig(f)
        print(f)

    for s in S.systems:
        cm = matplotlib.cm.get_cmap('hsv')
        i = 0
        for (q, k), B in s.ed_system.blocks.items():
            default_x = 6.4
            default_y = 4.8
            scale = 2
            plt.clf()
            plt.figure(figsize=(default_x*scale, default_y*scale), dpi=100)
            c = cm(i / len(s.ed_system.blocks))
            i += 1
            xdata = [x.real for x in B.energies]
            ydata = [x.imag for x in B.energies]
            plt.plot(xdata, ydata, marker='o', markersize=2.5, color=c, linestyle='')


            f = f'{s}_{q}_{k}_spectrum.png'
            f = os.path.join(S.graph_dir, f)
            plt.savefig(f)
        print(f)
