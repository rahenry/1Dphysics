#!/usr/bin/env python3
#
from plotters_utility import *
import plotters_linefuncs

class Plotter:
    list_info = []
    def process_info(self):
        for x in self.list_info:
            self.parameters[x] = []

        for l in self.info:
            sp = l.split()
            if sp[0]=='filter':
                for x in sp[1:]:
                    a, b = interpret_option(x)
                    self.system_filter[a] = b
            elif '=' in sp[0]:
                for x in sp:
                    if '=' in x:
                        a, b = interpret_option(x)
                        if a in self.parameters:
                            if not isinstance(self.parameters[a], list):
                                self.parameters[a] = [self.parameters[a]]
                            self.parameters[a].append(b)
                        else:
                            self.parameters[a] = b
            elif sp[0] in self.list_info:
                self.parameters[sp[0]].append(sp[1:])
            else:
                self.parameters[sp[0]] = ' '.join(sp[1:])
        1
    def __init__(self, S, grid, info):
        self.S = S
        self.grid = grid
        self.info = info
        self.parameters = {}
        self.nplots = 1
        self.system_filter = {}
        self.process_info()
        self.f12 = self.grid.font12
        self.fscale = self.grid.ncols / 2. * self.grid.scale / 4.
        self.fscale = self.grid.scale / 4.
        self.fscale=1 * self.grid.fontscale

        self.nticks = self.get('nticks', 3)
        self.fs_ticks = self.get('fs_ticks', 10) * self.fscale
        self.fs_axislabs = self.get('fs_axislabs', 12) * self.fscale
        self.fs_labs = self.get('fs_labs', 9) * self.fscale
        if not self.nticks:
            self.nticks = 3


        # if "sset_name" in self.system_filter:
        #     S1 = self.S.added_sets[self.system_filter["sset_name"]]
        #     self.S = S1
        self.slicer = S.slicer.find(self.system_filter)

    def process_axis(self):
        ax = self.ax
        if self.grid['noaxis']:
            axis_light(self.ax)
            # ax.xaxis.set_major_locator(plt.NullLocator())
            # ax.yaxis.set_major_locator(plt.NullLocator())
            # ax.set_xticklabels([])
            # ax.set_yticklabels([])
            # # ax.axis('off')
            # ax.set_xlabel(None)
            # ax.set_ylabel(None)
            # # self.ax.axis('off')
            # # self.ax.axis('tight')
            # for s in ax.spines:
            #     sp = ax.spines[s]
            #     sp.set_linestyle((0,(1,5)))
            #     # sp.set_linestyle(':')
            #     sp.set_color('gray')
            #     sp.set_alpha(0.5)

    def register_data(self, xdata, ydata):
        if self['transform_x'] or self['transform']:
            t = self['transform_x']
            if t is None: t = self['transform']
            xdata = self.transform_data(xdata, t)
        if self['transform_y']:
            t = self['transform_y']
            ydata = self.transform_data(ydata, t)
        return xdata, ydata

    def transform_data(self, data, t):
        res = [eval(f"{t}(d)") for d in data]
        return res

    def __len__(self):
        return self.nplots

    def draw(self):
        return

    def __getitem__(self, key):
        if not key in self.parameters:
            return None
        return self.parameters[key]

    def get(self, key, default, i=0):
        res = self[key]
        if res is None:
            res = default
        if isinstance(res, list):
            if i >= len(res):
                return res[-1]
            return res[i]
        return res

    def set_ticks(self, ax, **kwargs):
        xticks = numpy.linspace(*self.xlims, self.nticks)
        yticks = numpy.linspace(*self.ylims, self.nticks)
        yticks = eliminate_close(yticks)
        def mklb(x):
            return f'{x:.2f}'
        xlabs = [mklb(x) for x in xticks]
        ax.set_xticks(xticks, xlabs, **kwargs)
        ylabs = [mklb(x) for x in yticks]
        ax.set_yticks(yticks, ylabs, **kwargs)

    def test_inside(self, x, y):
        return x > self.xlims[0] and x < self.xlims[1] and y > self.ylims[0] and y < self.ylims[1]

    def make_label(self, x=None, y=None):
        lab = self["lab"]
        print(lab)
        if lab == "ETAS":
            if y is not None:
                s = y[0]
                e1 = s["eta_real"]
                e2 = s["eta_imag"]
                lab = f"$\eta = {e1:2.2f}+{e2:2.2f}i$"
                print(lab)
                return lab
        return lab

class Plotter1D(Plotter):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.xlims = [None, None]
        self.ylims = [None, None]

    def update_lims(self, lims, data):
        data = numpy.array(data)
        a, b = numpy.min(data), numpy.max(data)
        res = numpy.array(lims)
        if lims[0] is None or a < lims[0]:
            res[0] = a
        if lims[1] is None or b > lims[1]:
            res[1] = b
        return res

    def register_data(self, xdata=None, ydata=None):
        xdata, ydata = super().register_data(xdata, ydata)
        if xdata is not None:
            self.xlims = self.update_lims(self.xlims, xdata)
        if ydata is not None:
            self.ylims = self.update_lims(self.ylims, ydata)
        return xdata, ydata

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
        self.xlims = [xmin, xmax]
        self.ylims = [ymin, ymax]
        # self.register_data(xmin, ymin)
        # self.register_data(xmax, ymax)
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)

    def add_line(self, func):
        xdata = numpy.linspace(*self.xlims, 100)
        # exec(f'f = plotters_linefuncs.{func}', globals(), locals())
        a = f'plotters_linefuncs.{func}'
        def q(x):
            b = f'{a}(x)'
            z = eval(b)
            return eval(b)
        ydata = [q(x) for x in xdata]
        self.ax.plot(xdata, ydata)



class Trajectory(Plotter1D):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        if 'xkey' not in self.parameters:
            self.parameters['xkey'] = 'lambda'

        self.nplots = 0
        for x, y in self.slicer.get_sets(exclude=[self['xkey']]).items():
            self.nplots += 1


    def draw(self):
        i = 0
        cm = matplotlib.cm.get_cmap('gnuplot')
        for x, y in self.slicer.get_sets(exclude=[self['xkey']]).items():
            ax = self.grid.axs[i+self.grid.gs_offset]

            ind_s = 0
            for s in y:
                z = s[self['zkey']]
                xdata = numpy.real(z)
                ydata = numpy.imag(z)
                c = cm(ind_s / len(y))
                if self['colour']:
                    c = self['colour']
                ax.plot(xdata, ydata, linestyle='', marker='o', markersize=1.3, color=c)
                ind_s += 1

            i += 1

        super().draw()

class MatrixInfo(Plotter):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.nplots = 1
    def draw(self):
        s = self.slicer[0]
        ev = 'fr"'+self['message']+'"'
        text = fr'{eval(ev)}'
        ax = self.grid.axs[self.grid.gs_offset]
        ax.axis('off')
        text= r'\begin{tabular}{@{}l@{}}' + text +r'\end{tabular}'
        ax.text(0.0, 0.5, text, va='center', ha='left', size=13*self.fscale)

class MatrixPlot(Plotter):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.nplots = 1
    def draw(self):
        ax = self.grid.axs[self.grid.gs_offset]
        vals = self.slicer[0][self.parameters['key']]
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        p = ax.imshow(vals, origin='lower', vmin=0, rasterized=True)


class ETest(Plotter1D):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.nplots = 1

    def draw(self):
        ax = self.grid.axs[self.grid.gs_offset]
        s = self.slicer[0]
        vals = s['energies']
        xdata = [x.real for x in vals]
        ydata = [x.imag for x in vals]
        pltargs = dict(plotargs_default)
        ax.plot(xdata, ydata, **pltargs)

        vals = s['energies_matrix']
        xdata = [x.real for x in vals]
        ydata = [x.imag for x in vals]
        xdata, ydata = self.register_data(xdata, ydata)
        pltargs = dict(plotargs_default)
        pltargs.update({
            'marker':'x',
            'c' : 'red',
            'markersize' : 8,
            'zorder' : 5,
            })
        ax.plot(xdata, ydata, **pltargs, rasterized=True)

        ax.set_xlabel(r'$\Re(E)$', size=self.fs_axislabs)
        ax.set_ylabel(r'$\Im(E)$', size=self.fs_axislabs)

        self.set_ticks(ax, size=self.fs_ticks)



class LinePlot(Plotter1D):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.nplots = 1
        defaults = {
            'xkey' : 'lambda',
            'key' : 'entropy',
        }
        pltargs_default = {
            'linestyle' : '',
            'linewidth' : 0,
            'marker' : 'o',
            'markeredgewidth' : 0,
            'markersize' : 2,
            'alpha' : 1,
            }
        self.pltargs = dict(pltargs_default)
        for x, y in self.parameters.items():
            if x in pltargs_default:
                self.pltargs[x] = y

    def draw(self):
        ax = self.grid.axs[self.grid.gs_offset]
        self.ax = ax
        xdata = []
        ydata = None
        i = 0
        for s in self.slicer:
            xdata.append(s[self['xkey']])
            Y = s[self['key']]
            Y = numpy.array(Y)
            lY = 1
            try:
                lY = len(Y)
                Y = numpy.sort(Y)
            except:
                1
            if ydata is None:
                ydata = numpy.full((len(self.slicer), lY), 0.)

            try:
                iter(Y)
                # u = ydata[i,0:len(Y)]
                ydata[i,0:len(Y)] = Y
            except:
                try:
                    ydata[i,0] = Y
                except TypeError:
                    ydata[i,0] = numpy.sign(Y) * abs(Y)
            i += 1

        xdata, ydata = self.register_data(xdata, ydata)
        if self['add_line']:
            self.add_line(self['add_line'])
        if not self['auto_lims']:
            self.set_limits()
        # self.set_ticks(ax, size=self.fs_ticks)
        pltargs = dict(plotargs_default)
        pltargs.update(self.pltargs)
        # if self['sequence']:
        #     for i in range(len(xdata)):
        #         p = ax.plot(xdata[i], ydata[i], **pltargs, rasterized=True)
        # else:
        p = ax.plot(xdata, ydata, **pltargs, rasterized=True)
        xlab = self.get('xlab', r'$\lambda$')
        ax.set_xlabel(xlab, size=self.fs_axislabs)
        ylab = self.get('ylab', '?')
        ax.set_ylabel(ylab, size=self.fs_axislabs)

        if self['lab']:
            pos = self.get('pos', 'upper left')
            c = self.get('lab_col', 'lightgrey')
            draw_label(ax, self['lab'], pos, self.fs_labs, c=c)


class Scatter(Plotter1D):
    list_info = ['scatter', 'zkey']
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.nplots = 0
        # if not self['xkey'] and not self['zkey']:
        #     self.parameters['zkey'] = ['lambda']
        self.ex = self['zkey']
        if self['xkey']:
            self.ex.append(self['xkey'])
        for x, y in self.slicer.get_sets(exclude=self.ex).items():
            self.nplots += 1

    def draw(self):
        defaults = {
            'cmap' : 'cool',
            'colour' : None,
            'key' : 'pfs',
            'transy' : None,
            'transx' : None,
            }
        pltargs_default = {
            'linestyle' : '',
            'marker' : 'o',
            'markersize' : 3,
            }

        i = 0
        inc = []
        if self['order']:
            order = self['order'].split(',')
            order.reverse()
            inc = order
        for x, y in self.slicer.get_sets(exclude=self.ex, include=inc).items():
            ax = self.grid.axs[self.grid.gs_offset+i]
            self.ax = ax

            for data in self['scatter']:
                args = dict(defaults)
                pltargs = dict(pltargs_default)
                for q in data:
                    if '=' in q:
                        a, b = interpret_option(q)
                        if a in defaults:
                            args[a] = b
                        else:
                            pltargs[a] = b
                cm = matplotlib.cm.get_cmap(args['cmap'])

                if self['xkey'] is None:
                    ind_s = 0
                    for s in y:
                        self.xlims = [None, None]
                        self.ylims = [None, None]
                        z = s[args['key']]

                        xdata = numpy.real(z)
                        ydata = numpy.imag(z)
                        xdata, ydata = self.register_data(xdata, ydata)

                        if args['colour']:
                            col = args['colour']
                            if col == 'rand':
                                for k in range(len(xdata)):
                                    c = cm(random.random()*0.7)
                                    ax.plot([xdata[k]], [ydata[k]], **pltargs, color=c, rasterized=True)

                            else:
                                c = args['colour']
                                ax.plot(xdata, ydata, **pltargs, color=c, rasterized=True)
                        else:
                            c = cm(ind_s / len(y))
                            ax.plot(xdata, ydata, **pltargs, color=c, rasterized=True)
                        ind_s += 1
                else:
                    sets = y.get_sets(exclude=[self['xkey']])
                    ind_s = 0
                    for x1, y1 in sets.items():
                        xdata = []
                        ydata = []
                        for s in y1:
                            z = s[args['key']]
                            z = numpy.real(z)
                            xdata.append(s[self['xkey']])
                            ydata.append(z)
                        c = cm(ind_s / len(sets))
                        if args['colour']:
                            c = args['colour']
                        if 'transy' in args:
                            ydata = transform_data(ydata, args['transy'])
                        xdata, ydata = self.register_data(xdata, ydata)
                        ax.plot(xdata, ydata, **pltargs, color=c, rasterized=True)
                        ind_s += 1

            xlab = self.get('xlab', r'$\lambda$')
            if not xlab == 'none':
                ax.set_xlabel(xlab, size=self.fs_axislabs)
            ylab = self.get('ylab', '?')
            if not ylab == 'none':
                ax.set_ylabel(ylab, size=self.fs_axislabs)

            # if self['xlab']:
            #     ax.set_xlabel(self['xlab'])
            # if self['ylab']:
            #     ax.set_ylabel(self['ylab'])

            i += 1
            self.set_limits()
            self.set_ticks(ax, size=self.fs_ticks)
            matplotlib_set_aspect(ax, 1)
            if self['lab']:
                lab = self.make_label(x, y)
                pos = self.get('pos', 'upper left')
                c = self.get('lab_col', 'lightgrey')
                draw_label(ax, lab, pos, self.fs_labs, c=c)
            elif not len(x) == 0:
                lab = make_latex_label(x)
                pos = self.get('label_pos', 'upper right', i)
                draw_label(ax, lab, pos, self.fs_labs)

            self.process_axis()

        super().draw()

class LinePlot2(Plotter1D):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.nplots = 1
        self.pltargs = {}
        pltargs = self['pltargs']
        if pltargs:
            for x in pltargs.split()[0:]:
                a, b = interpret_option(x)
                self.pltargs[a] = b


    def draw(self):
        ax = self.grid.axs[self.grid.gs_offset]
        self.ax = ax
        xdata = []
        ydata = None
        pltargs = dict(plotargs_default)
        pltargs.update(self.pltargs)
        i = 0
        for x, y in self.slicer.get_sets(exclude=[self['xkey']]).items():
            xdata = numpy.full(len(y), 0.)
            ydata = numpy.full(len(y), 0.)
            for i in range(len(y)):
                s = y[i]
                xdata[i] = s[self['xkey']]
                ydata[i] = s[self['key']]

            xdata, ydata = self.register_data(xdata, ydata)
            p = ax.plot(xdata, ydata, **pltargs, rasterized=True)

        xlab = self.get('xlab', r'$\lambda$')
        if xlab is not 'none':
            ax.set_xlabel(xlab, size=self.fs_axislabs)
        ylab = self.get('ylab', '?')
        if ylab is not 'none':
            ax.set_ylabel(ylab, size=self.fs_axislabs)
        if self['lab']:
            pos = self.get('pos', 'upper left')
            c = self.get('lab_col', 'lightgrey')
            draw_label(ax, self['lab'], pos, self.fs_labs, c=c)

        self.set_limits()
        self.set_ticks(ax, size=self.fs_ticks)


class Plot2D(Plotter):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        xkey = self['xkey']
        ykey = self['ykey']
        self.nplots = len(self.slicer.get_sets(exclude=[xkey, ykey]))
        if "sset_name" in self.system_filter:
            self.eps_data = self.S.added_sets[self.system_filter["sset_name"]].eps_data

    def get_2d_data(self):
        xkey = self['xkey']
        ykey = self['ykey']
        self.data = {}
        for a, b in self.slicer.get_sets(exclude=[xkey, ykey]).items():
            xvals = b.get_possible_values(xkey)
            yvals = b.get_possible_values(ykey)
            X = []
            Z = []
            for s in b:
                x = s[xkey]
                y = s[ykey]
                z = s[self['key']]
                if z == []:
                    z = 0
                Z.append(z)
                X.append((x, y))
            z, ex = rastify(X, Z)
            self.data[a] = (z, ex)


    def draw(self):
        self.get_2d_data()
        i = 0
        for a, (z, ex) in self.data.items():
            ax = self.grid.axs[self.grid.gs_offset+i]

            aspect_ratio = 1
            aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio
            # exit()
            ext1 = [self.get('extent_left', ex[0]),
            self.get('extent_right', ex[1]),
            self.get('extent_lower', ex[2]),
            self.get('extent_upper', ex[3])]
            a1 = self.get('cmin', 0)
            a2 = self.get('cmax', 1)
            cmap = self.get('cmap', 'viridis')
            cmap = plt.get_cmap(cmap)
            cmap = truncate_colormap(cmap, a1, a2, n=1000)
            p = ax.imshow(z, origin='lower', rasterized=True, extent=ext1, aspect=aspect, cmap=cmap)
            insert_cbar(self.grid.fig, ax, p, self)
            # pos=self.get('cbar_pos', 'lower right'))
            i += 1

            xlab = self.get('xlab', r'$\lambda$')
            ax.set_xlabel(xlab, size=self.fs_axislabs)
            ylab = self.get('ylab', '?')
            ax.set_ylabel(ylab, size=self.fs_axislabs)
            if self['lab']:
                pos = self.get('pos', 'upper left')
                c = self.get('lab_col', 'lightgrey')
                draw_label(ax, self['lab'], pos, self.fs_labs, c=c)


