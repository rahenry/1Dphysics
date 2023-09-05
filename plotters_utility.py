#!/usr/bin/env python3

import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import os
import matplotlib
import numpy
from rah_utility import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from labels import *
import random
from cmap import truncate_colormap
import plotly.graph_objects as go
option_types = [int, float, str]

def insert_cbar(fig, ax, p, plotter):
                # include_zero=0, pos='lower right'):
    include_zero = plotter.get('include_zero', 0)
    pos = plotter.get('cbar_pos', 'lower left')
    side = plotter.get('cbar_side', 'left')
    ticks = numpy.linspace(p.norm._vmin, p.norm._vmax, 3)
    if include_zero:
        ticks[0] = 0
    axins1 = inset_axes(ax,
                    width="4%",  # width = 50% of parent_bbox width
                    height="40%",  # height : 5%
                    loc=pos)


    fig.colorbar(p, cax=axins1, orientation="vertical", ticks=ticks)
    axins1.yaxis.set_ticks_position(side)

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

def deduce_raster_params(X):
    # unique_sorted = numpy.vstack((numpy.unique(v) for v in X.T)).T
    things = [numpy.unique(v) for v in X.T]
    unique_sorted = numpy.vstack(things).T
    d_min = unique_sorted[0] # x min, y min
    d_max = unique_sorted[-1] # x max, y max
    d_step = unique_sorted[1]-unique_sorted[0] # x, y step
    nsamples = (numpy.round((d_max - d_min) / d_step) + 1).astype(int)
    return d_min, d_max, d_step, nsamples

def rastify(X, y):
    X = numpy.array(X)
    d_min, d_max, d_step, nsamples = deduce_raster_params(X)
    ind = numpy.round((X - d_min) / d_step).T.astype(int)
    ind = ind.T
    z = numpy.zeros(nsamples)
    for i in range(len(X)):
        a, b = ind[i]
        z[b][a] = y[i]
    extent = numpy.vstack((d_min, d_max)).T.ravel()
    z = numpy.array(z)
    return z, extent

plotargs_default = {
    'linestyle' : '',
    'marker' : 'o',
    'markersize' : 3,
    'zorder' : 10,
    }

def interpret_option(x):
    a, b = x.split('=')
    b = b.replace('__', ' ')
    b = b.replace('EQ', '=')
    for t in option_types:
        try:
            b = t(b)
            break
        except:
            pass
    return a, b


def matplotlib_set_aspect(ax, ratio):
    # ratio = 1.0
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)

def eliminate_close(v):
    v = numpy.sort(v)
    res = []
    tol = 1E-10 * numpy.max(v) + 1E-16
    for x in v:
        success = 1
        for y in res:
            if rel_diff(x, y) < tol:
                success = 0
        if success:
            res.append(x)

    return res

import matplotlib.transforms as transforms
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt

def draw_label(ax, label, position='upper right', scale=1, c='lightgrey'):
    text_box = AnchoredText(label, frameon=True, loc=position, pad=0, prop=dict(size=scale))
    # text_box.patch.set_boxstyle("facecolor='lightgrey',linestyle=''")
    text_box.patch.set_boxstyle(f"square,pad=0.15")
    plt.setp(text_box.patch, linestyle='', facecolor=c)
    # ax.setp(text_box.patch, facecolor='white', alpha=0.5)
    ax.add_artist(text_box)


def transform_data(data, key):
    if isinstance(data, list):
        data = numpy.array(data)
    res = data
    if key == 'log':
        res = numpy.log(data)
    if key == 'log10':
        res = numpy.log10(data)
    return res

def labelify(x, ex=['sset_name']):
    if isinstance(x, dict):
        x = x.items()
        print(x)
    res = ''
    for a, b in x:
        if a in ex:
            continue
        res += f'${a} = {b}$, '
    res = res.strip()
    res = res.strip(',')
    return res

def axis_off(ax):
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ax.axis('off')
    ax.axis('tight')

def axis_light(ax):
    # print(ax.get_xlim(), ax.get_ylim())
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # ax.axis('off')
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    # self.ax.axis('off')
    # self.ax.axis('tight')
    for s in ax.spines:
        sp = ax.spines[s]
        sp.set_linestyle((0,(1,5)))
        # sp.set_linestyle(':')
        sp.set_color('gray')
        # sp.set_alpha(0.5)

def rescale_cmap(cmap, a, b, z):
    res = []
    npoints = 100
    for i in range(npoints):
        cind = a + (b-a) * (i / npoints)
        c = cmap(cind)
        # c = cmap((depth-d-1) / depth)
        c = [x * z for x in c]
        res.append(c)

    return matplotlib.colors.LinearSegmentedColormap.from_list("map1", res)
