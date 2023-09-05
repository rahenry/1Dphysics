#!/usr/bin/env python3

import plotters_linefuncs

from plotters import *

class PlotterMomentum(LinePlot):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.s0 = S.systems[0]
        self.nplots = len(self.s0['block_energies'])
        print(self.nplots)
        for b in self.s0['block_energies']:
            print(b)

    def set_limits(self):
        if self['xinclude'] is not None:
            self.xlims = self.update_lims(self.xlims, [self['xinclude']])
        if self['yinclude'] is not None:
            self.ylims = self.update_lims(self.ylims, [self['yinclude']])
        d = self.get('lim_buffer', 0.05)
        dx = d * (self.xlims[1] - self.xlims[0])
        dy = d * (self.xlims[1] - self.xlims[0])
        xmin = self.get('xmin', self.xlims[0]-dx)
        xmax = self.get('xmax', self.xlims[1]+dx)
        ymin = self.get('ymin', self.ylims[0]-dx)
        ymax = self.get('ymax', self.ylims[1]+dx)
        # self.xlims = [xmin, xmax]
        # self.ylims = [ymin, ymax]
        # self.register_data(xmin, ymin)
        # self.register_data(xmax, ymax)
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)

    def draw(self):
        i = 0
        ex = []
        for a in self.s0['block_energies']:
            self.ax = self.grid.axs[i+self.grid.gs_offset]
            # self.ax.axis('off')
            # self.ax.axis('tight')
            axis_light(self.ax)
            i += 1
            for x, y in self.slicer.get_sets(exclude=ex).items():
                b0 = self.s0['block_energies'][a]
                d = len(b0)
                xdata = numpy.empty(d, dtype=float)
                ydata = numpy.empty(d, dtype=float)
                j = 0
                s = y[j]
                e = s['block_energies'][a]
                xdata = numpy.real(e)
                ydata = numpy.imag(e)
                xdata, ydata = self.register_data(xdata, ydata)

        i = 0
        ex = []
        for a in self.s0['block_energies']:
            self.ax = self.grid.axs[i+self.grid.gs_offset]
            bname = f'$Q = {a[0]}$, $k = {a[1]}$'
            i += 1
            for x, y in self.slicer.get_sets(exclude=ex).items():
                b0 = self.s0['block_energies'][a]
                d = len(b0)
                xdata = numpy.empty(d, dtype=float)
                ydata = numpy.empty(d, dtype=float)
                j = 0
                s = y[j]
                e = s['block_energies'][a]
                xdata = numpy.real(e)
                ydata = numpy.imag(e)
                # xdata, ydata = self.register_data(xdata, ydata)

                self.set_limits()
                c = self.get('lab_col', 'lightgrey')
                pos = self.get('pos', 'upper left')
                draw_label(self.ax, bname, 'upper left', self.fs_labs, c=c)
                self.ax.plot(xdata, ydata, **self.pltargs)
                self.process_axis()

class PlotterMomentumMatrix(PlotterMomentum):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.s0 = S.systems[0]
        i = 0
        for x, M in self.s0['block_measurements'].items():
            if self['key'] in M:
                i += 1
        self.nplots = i
        print(self.nplots)

    def draw(self):
        ex = []
        i = 0
        for a, M in self.s0['block_measurements'].items():
            if self['key'] not in M:
                continue
            self.ax = self.grid.axs[i+self.grid.gs_offset]
            bname = f'$Q = {a[0]}$, $k = {a[1]}$'
            i += 1
            for x, y in self.slicer.get_sets(exclude=ex).items():
                b0 = self.s0['block_measurements'][a]
                d = len(b0)
                z = b0[self['key']]
                self.ax.imshow(z)

                c = self.get('lab_col', 'lightgrey')
                pos = self.get('pos', 'upper left')
                lab = bname
                if self['key'] == 'zero_info':
                    lab = 'hi'
                    lab = f'{bname}, $n_0 = {b0["n_zeros"]}$, $D = {b0["zero_det"]}$'
                    lab = f'{bname}, $n_0 = {b0["n_zeros"]}$'
                draw_label(self.ax, lab, 'upper left', self.fs_labs, c=c)
                self.process_axis()
