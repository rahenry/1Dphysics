#!/usr/bin/env python3

import matplotlib.pyplot as plt
from rah_utility import mkdir, matplotlib_set_aspect
from Slicer import SystemSlicer
import numpy
import os
from rah_utility import mkdir
import matplotlib
from scipy.optimize import curve_fit
from ansatz import generate_ansatzs
from args import Q
# matplotlib.use('pgf')

labels = {
    'fsL' : r'$\chi_s / L$',
    # 'x' : r'$\langle X_1 X_L^\dagger \rangle / L)$',
    # 'z' : r'$\langle Z_1 Z_L^\dagger \rangle / L)$',
    'x' : r'$x_{1L}$',
    'z' : r'$z_{1L}$',
    'z1L' : r'$z_{1L}$',
    'gamma' : r'$\gamma',
    'o' : r'$\langle L | R \rangle$',
    }

def savefig(name):
    plt.savefig(name, bbox_inches='tight')

def make_label(x):
    if x in labels:
        return labels[x]
    return x

class DataSet:
    def __init__(self, sset=None, ykey=None, from_data=None,
                 ansatz=None, mode=None, xkey=None, separate=[], name=None,
                 **plot_defaults):
        self.xkey = xkey
        self.ykey = ykey
        self.sset = sset
        self.data = None
        # if name is None and sset is not None:
        #     name = sset.name
        self.name = name
        self.plot_defaults=plot_defaults
        self.labels = []
        if from_data is not None:
            self.data = from_data
        elif sset is not None:
            if self.xkey is None:
                self.xkey = sset['xkey']
            exclude = [self.xkey]
            if mode == "in_system":
                self.data = sset.gds(ykey, include=separate, exclude=exclude)
            else:
                self.data = sset.gd(ykey, include=separate, exclude=exclude, xcoord=self.xkey)
            # for key, vals in self.data.items():
            #     print(key, len(vals))

        else:
            print(f'Cannot initialise data_set')
            exit()

        for key, vals in self.data.items():
            self.labels.append(key)
            vals['ansatz_data'] = {}

        self.output_dir = sset.graph_dir
        self.pgf_output_dir = os.path.join(sset.graph_dir, 'pgf')
        mkdir(self.pgf_output_dir)
        self.available_ansatzs = generate_ansatzs(['default', self.sset.name])
        if ansatz:
            self.make_ansatz(ansatz)

    def __getitem__(self, x):
        return self.data[x]

    def plot(self, aspect=0.5, use_legend=True, **kwargs):
        kwargs.update(self.plot_defaults)
        # w, h = figaspect(aspect)
        # fig = Figure(figsize=(w, h))
        plt.clf()
        fig = plt.figure()
        ax1 = plt.subplot()
        self._plot(ax1, **kwargs)
        # save
        if use_legend:
            ax1.legend(self.labels, bbox_to_anchor=(1.04, 1), loc="upper left")
        matplotlib_set_aspect(ax1, aspect)
        fig.tight_layout()
        f0 = self.ykey
        if 'name'in kwargs: f0 = kwargs['name']
        elif self.name:
            f0 = self.name
        f = os.path.join(self.output_dir, f'{f0}.png')
        savefig(f)
        # plt.savefig(f)

        if not Q.skip_pgf:
            f = os.path.join(self.pgf_output_dir, f'{f0}.pgf')
            savefig(f)
        # plt.savefig(f)

    def _plot(self, ax, ylims=None, xlims=None, name=None, xlab=None, ylab=None, apply_func=None, **kwargs):
        if ylims is None: ylims = [None, None]
        if xlims is None: xlims = [None, None]
        defaults = {'markersize':5, 'marker':'o'}
        defaults.update(kwargs)

        for key, vals in self.data.items():
            X = vals['xdata']
            Y = vals['ydata']
            if apply_func:
                try:
                    Y = apply_func(Y)
                except:
                    print(vals)
                    print(Y)
                    print(type(Y))
                    print("Failed to apply function")
                    exit()

            ax.plot(X, Y, **defaults)

        n_ansatz = 0
        for key, vals in self.data.items():
            n_ansatz += len(vals['ansatz_data'])
        cm_ansatz = matplotlib.cm.get_cmap('cool')

        ind = 0
        for key, vals in self.data.items():
            for ansatz_name, d in vals['ansatz_data'].items():
                c = cm_ansatz(ind / n_ansatz)
                ind += 1
                ax.plot(d['X'], d['Y'], marker='', linestyle='--', color=c)

        ax.set_ylim(*ylims)
        ax.set_xlim(*xlims)

        # extras
        if xlab is None:
            xlab = self.xkey
        xlab = make_label(xlab)
        if ylab is None:
            ylab = self.ykey
        ylab = make_label(ylab)

        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)

    def make_ansatz(self, name, Xfraction=None, label=None):
        funcstr = self.available_ansatzs[name]
        # scipy.optimize.curve_fit(f, xdata, ydata,
        c = f'def {funcstr}'
        exec(c, globals())
        for key, data in self.data.items():
            data['ansatz_data'][name] = {}
            d = data['ansatz_data'][name]
            X = data['xdata']
            Y = data['ydata']
            if Xfraction:
                n = int(len(X) * Xfraction)
                X = X[:n]
                Y = Y[:n]

            p, pcov = curve_fit(f, X, Y)
            d['X'] = numpy.linspace(X[0], X[-1], 100)
            d['Y'] = [f(x, *p) for x in d['X']]
            d['p'] = p
            d['pcov'] = pcov


        for key, data in self.data.items():
            d = data['ansatz_data'][name]
        if label is None:
            label = name
        self.labels.append(label)
        return d['p'], d['pcov']

    def print_ansatz(self):
        for key, data in self.data.items():
            for name, d in data['ansatz_data'].items():
                print(name, d['p'])


def plot_multiple(ds, name=None, aspect=1, **kwargs):
    n = len(ds)
    plt.clf()
    fig, axs = plt.subplots(1, n)
    name_default = ''
    for i in range(n):
        d = ds[i]
        p = d._plot(axs[i], **kwargs)
        name_default += d.ykey + '_'
        matplotlib_set_aspect(axs[i], aspect)

    axs[-1].legend(ds[0].labels, bbox_to_anchor=(1.04, 1), loc="upper left", fontsize='x-small', markerscale=0.7)

    fig.tight_layout()
    f0 = name_default
    if name:
        f0 = name
    s1 = ds[0]
    f = os.path.join(s1.output_dir, f'{f0}.png')
    savefig(f)

    if not Q.skip_pgf:
        f = os.path.join(s1.pgf_output_dir, f'{f0}.pgf')
        savefig(f)


def make_ansatz_data():
    1
