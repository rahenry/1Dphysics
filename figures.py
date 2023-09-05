#!/usr/bin/env python3


import matplotlib.pyplot as plt
from Slicer import SystemSlicer
import numpy
import os
from rah_utility import mkdir
import matplotlib
# matplotlib.use('pgf')

def graph1(sset, xcoord, y, exclude=[], rescale_to_max = False, ylims=[None, None]):
    s = sset.systems[0]
    print(s)
    if not y in s.measurement_data:
        print(f"Can't find measurement {y} to graph")
        # for k in s.measurement_data:
        #     print(k)
        return
    xdatas = []
    ydatas = []
    labels = []
    ymin = None
    ymax = None
    if rescale_to_max:
        for key, systems in sset.slicer.get_sets(exclude=exclude+[xcoord]).items():
            for s in systems:
                y = s[y]
                if ymin is None:
                    ymin = y
                elif y < ymin:
                    ymin = y
                if ymax is None:
                    ymax = y
                elif y > ymax:
                    ymax = y

    for key, systems in sset.slicer.get_sets(exclude=exclude+[xcoord]).items():
        xdata = []
        ydata = []
        for s in systems:
            xdata.append(s[xcoord])
            ydata.append(s[y])
        xdatas.append(xdata)
        ydatas.append(ydata)
        labels.append(key)
    ydatas = numpy.swapaxes(ydatas, 0, 1)
    plt.clf()
    plt.plot(xdatas[0], ydatas, marker='o')
    plt.ylim(*ylims)
    plt.legend(labels)
    plt.xlabel(xcoord)
    plt.ylabel(y)

    f = f'{y}_{xcoord}'
    if exclude:
        f += f'_{exclude}'
    f = os.path.join(sset.graph_dir, f'{f}.png')
    plt.savefig(f)


def graph_all(sset, xcoord, y, exclude=[], rescale_to_max = False, ylims=[None, None]):
    for key, systems in sset.slicer.get_sets(exclude=exclude+[xcoord]).items():
        xdata = []
        ydata = []
        for s in systems:
            xdata.append(s[xcoord])
            ydata.append(s[y])
        plt.clf()
        plt.plot(xdata, ydata, marker='o')
        plt.xlabel(xcoord)
        plt.ylabel(y)
        base_dir = os.path.join(sset.graph_dir, y)
        mkdir(base_dir)
        f = os.path.join(base_dir, f'{key}.png')
        plt.savefig(f)

def graph2(sset, xcoord, y, exclude=[], rescale_to_max = False, ylims=[None, None], apply_func=None, title=None,
           xlab=None, ylab=None, **kwargs):
    s = sset.systems[0]
    if ylims is None:
        ylims = [None, None]
    # print(s)
    if not y in s.measurement_data:
        print(f"Can't find measurement {y} to graph")
        # for k in s.measurement_data:
        #     print(k)
        return
    xdatas = []
    ydatas = []
    labels = []
    ymin = None
    ymax = None
    if rescale_to_max:
        for key, systems in sset.slicer.get_sets(exclude=exclude+[xcoord]).items():
            for s in systems:
                y = s[y]
                if ymin is None:
                    ymin = y
                elif y < ymin:
                    ymin = y
                if ymax is None:
                    ymax = y
                elif y > ymax:
                    ymax = y

    for key, systems in sset.slicer.get_sets(exclude=exclude+[xcoord]).items():
        xdata = []
        ydata = []
        for s in systems:
            xdata.append(s[xcoord])
            ydata.append(s[y])
        xdatas.append(xdata)
        if apply_func:
            ydata = [apply_func(y) for y in ydata]
        ydatas.append(ydata)
        labels.append(key)
    # print(ydatas)
    # print(ydatas.shape)
    plt.clf()
    fig, (ax1) = plt.subplots(1)
    defaults = {'markersize':5, 'marker':'o'}
    defaults.update(kwargs)
    for i in range(len(xdatas)):
        X = xdatas[i]
        Y = ydatas[i]
        ax1.plot(X, Y, **defaults)
    # ydatas = numpy.swapaxes(ydatas, 0, 1)
    ax1.set_ylim(*ylims)
    ax1.legend(labels, bbox_to_anchor=(1.04, 1), loc="upper left")
    ax1.set_xlabel(xcoord)
    ax1.set_ylabel(y)
    if xlab:
        ax1.set_xlabel(xlab)
    if ylab:
        ax1.set_ylabel(ylab)
    fig.tight_layout()

    f = f'{y}_{xcoord}'
    if exclude:
        f += f'_{exclude}'
    if title:
        f = f'{title}'
    plt.savefig(os.path.join(sset.graph_dir, f'{f}.png'))
    f = os.path.join(sset.graph_dir, f'{f}.pgf')
    plt.savefig(f)

def graph_sets(sset, data, title, **kwargs):
    labels = [k for k in data]
    plt.clf()
    fig, (ax1) = plt.subplots(1)
    defaults = {'markersize':5, 'marker':'o'}
    defaults.update(kwargs)
    for key, vals in data.items():
        print(key)
        X = vals['xdata']
        Y = vals['ydata']
        ax1.plot(X, Y, **defaults)
    f = os.path.join(sset.graph_dir, f'{title}.png')
    ax1.legend(labels, bbox_to_anchor=(1.04, 1), loc="upper left")
    fig.tight_layout()
    plt.savefig(f)


import numpy as np
def deduce_raster_params(X):
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

def rastify(X, y):
    X = numpy.array(X)
    d_min, d_max, d_step, nsamples = deduce_raster_params(X)
    ind = np.round((X - d_min) / d_step).T.astype(int)
    ind = ind.T
    z = numpy.zeros(nsamples)
    for i in range(len(X)):
        a, b = ind[i]
        z[b][a] = y[i]


    extent = np.vstack((d_min, d_max)).T.ravel()
    z = numpy.array(z)
    return z, extent

def savefig(S, name):
    f = os.path.join(S.graph_dir, f'{name}.png')
    plt.savefig(f)
    f = os.path.join(S.pgf_dir, f'{name}.pgf')
    plt.savefig(f)
