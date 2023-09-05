#!/usr/bin/env python3

from plotters import *

class EP_Plot(Plotter):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.nplots = 1
        # if self['draw_eps']:
        self.eps_data = self.S.eps_data
        if "sset_name" in self.system_filter:
            self.eps_data = self.S.added_sets[self.system_filter["sset_name"]].eps_data

    def get_2d_data(self):
        X = []
        Y = []
        xkey = self.get("xkey", "lambda")
        print(self['key'])
        for s in self.slicer:
            x = s[f'{xkey}_real']
            y = s[f'{xkey}_imag']
            z = s[self['key']]
            X.append((x, y))
            Y.append(z)

        z, ex = rastify(X, Y)
        return z, ex

    def draw(self):
        z, ex = None, None
        s = self.slicer[0]
        L, N = s['L'], s['N']
        cid = (L, N)
        # labkey = r'\lambda'
        labkey = self.get("labkey", r"\lambda")
        if self['key'] == 'epk':
            d = self.eps_data[cid]
            z = d['z']
            ex = d['ex']
            labkey = 'k'
        else:
            z, ex = self.get_2d_data()


        aspect_ratio = 1
        aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio
        ax = self.grid.axs[self.grid.gs_offset]
        self.cm = matplotlib.colormaps['viridis']
        if True:
            cm = matplotlib.colormaps['plasma']
            self.cm = rescale_cmap(cm, 0.15, 1.0, 1.0)
        imargs = {"vmin":0,
                  "cmap":self.cm,
                  }
        # if self["two_colour_EPs"]:
        #     imargs.update(dict(cmap='gray', vmin=0))
        p = ax.imshow(z, origin='lower', extent=ex, aspect=aspect, rasterized=True, **imargs)
        ax.set_xlabel(fr'$\Re({labkey})$', size=self.fs_axislabs)
        ax.set_ylabel(fr'$\Im({labkey})$', size=self.fs_axislabs)


        pos = self.get('label_pos', 'upper right')
        lab = self.get('label', None)
        if lab:
            draw_label(ax, lab, pos, self.fs_labs)

        self.xlims = (ex[0], ex[1])
        self.ylims = (ex[2], ex[3])
        self.set_ticks(ax, size=self.fs_ticks)

        # insert_cbar(self.grid.fig, ax, p)
        insert_cbar(self.grid.fig, ax, p, self)

        pltEPargs = {
            'marker' : 'x',
            'markersize' : 7,
            'color' : 'white',
            'linewidth' : 1.5,
            # 'alpha' : 0.80,
            }
        pargs_plus = {
            "marker" : "+",
            'markersize' : 8,
            'linewidth' : 2.0,
            }
        if self['draw_eps']:
            skip_trivial = self.get("skip_trivial", 1)
            d = self.eps_data[cid]
            if self['key'] == 'epk':
                for k in d['minima_k']:
                    pargs = dict(pltEPargs)
                    if self["two_colour_EPs"] and k.imag < 0:
                        pargs.update(pargs_plus)
                    if skip_trivial and abs(k) < 1E-6:
                        continue
                    ax.plot([k.real], [k.imag], **pargs)
            else:
                if self["no_rotate"]:
                    xkey = self.get("xkey", "lambda")
                    minima = d["minima"]
                    for i, a in enumerate(d["minima"]):
                        if "gamma" in xkey:
                            a = (1 - a) / (1 + a)
                        pargs = dict(pltEPargs)
                        if self["two_colour_EPs"]:
                            j = i % (len(minima)//2)
                            k = d["minima"][j]
                            if i >= len(minima) // 2:
                                pargs.update({
                                    'color' : 'black',
                                })
                            else:
                                pargs.update({
                                    # 'color' : 'green',
                                    'color' : 'white',
                                })
                            if k.imag < 0:
                                pargs.update(pargs_plus)
                        if not self.test_inside(a.real, a.imag):
                            continue
                        ax.plot([a.real], [a.imag], **pargs)
                else:
                    minima = rotate_minima(d['minima'], N)
                    # pltEPargs["color"] = "red"
                    for i, a in enumerate(minima):
                        # print(i, len(minima) // 2, a)
                        if skip_trivial:
                            ang = numpy.angle(a)
                            ang = ang / 2 / numpy.pi * N
                            skip = 0
                            for nn in range(-N-1, N):
                                if rel_diff(ang, nn) < 1E-6:
                                    skip = 1
                            if skip: continue

                        if not self.test_inside(a.real, a.imag):
                            continue
                        ax.plot([a.real], [a.imag], **pltEPargs)

        if self['lab']:
            pos = self.get('pos', 'upper left')
            c = self.get('lab_col', 'lightgrey')
            draw_label(ax, self['lab'], pos, self.fs_labs, c=c)
