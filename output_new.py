#!/usr/bin/env python3
#
from plotters import *
from EP_Plot import *
from plotters_table import *
from plotters_pbc import *

def process_output(S):
    print("Processing output")
    spec = S.output_spec

    grids = []
    plot_specs = []

    for line in spec:
        line = line.strip()
        if line == '': continue
        if line[0] == '#': continue
        if line[0] == '$': continue

        sp = line.split()
        if sp[0] == 'new':
            if not plot_specs == []:
                grids[-1].process(plot_specs)
            grids.append(PlotGrid(S, sp[1:]))
            plot_specs = []
        elif sp[0] == 'plot':
            ptype = sp[1]
            plot_specs.append((ptype, []))
        else:
            plot_specs[-1][-1].append(line)

    if not plot_specs == []:
        grids[-1].process(plot_specs)

    for g in grids:
        g.init_image()
        g.draw()


options_default = {
    'ncols' : 3,
    }
class PlotGrid:
    def __init__(self, S, options=[]):
        self.S = S
        self.nplots = 0
        self.plotters = []
        self.flags = []
        self.options = dict(options_default)
        for x in options:
            if '=' in x:
                a, b = interpret_option(x)
                self.options[a] = b
            else:
                self.flags.append(x)
        self.ncols = self['ncols']
        scale = 4
        if 'scale' in self.options: scale = self['scale']
        self.scale = scale
        self.fontscale = self.scale / 2.
        self.fontscale = 1
        if 'fontscale' in self.options:
            self.fontscale *= self['fontscale']
        self.font12 = int(12 * self.fontscale)

    def __getitem__(self, key):
        if key not in self.options:
            return None
        return self.options[key]

    def __setitem__(self, key, val):
        self.options[key] = val

    def process(self, plot_specs):
        for ptype, info in plot_specs:
            ev = f'{ptype}(self.S, self, info)'
            # q = Trajectory(self.S, self, info)
            try:
                eval(ptype)
            except:
                print(f'Unknown ptype {ptype}')
                continue
            try:
                pnew = eval(ev)
                self.plotters.append(pnew)
            except:
                raise

        for p in self.plotters:
            self.nplots += p.nplots
        # self.nrows = self.nplots // self.ncols + self.nplots % self.ncols
        # self.nrows = self.nplots // self.ncols + 1
        self.nrows = (self.nplots + self.ncols - 1) // self.ncols

    def init_image(self):
        dpi = self['dpi']
        if self['nosize']:
            self.fig = plt.figure(dpi=dpi)
        else:
            self.fig = plt.figure(figsize=(self.ncols*self.scale, self.nrows*self.scale), dpi=dpi)
            print(self.ncols, self.nrows, self.nplots)

        self.gs = gridspec.GridSpec(self.nrows, self.ncols)
        self.axs = []
        for i in range(self.nplots):
            if self['transpose']:
                j = (i % self.nrows) * self.ncols + i // self.nrows
                i = j
            ax = plt.subplot(self.gs[i])
            self.axs.append(ax)

    def draw(self):
        self.gs_offset = 0
        print("Drawing")
        for p in self.plotters:
            p.draw()
            self.gs_offset += p.nplots

        if 'letter_labels' in self.flags:
            for i in range(self.nplots):
                ax = self.axs[i]
                loc = 'upper left'
                letter = chr(ord('a') + i)
                text = f'({letter})'
                text_box = AnchoredText(text, frameon=True, loc=loc, pad=0)
                text_box.patch.set_boxstyle(f"square,pad=0.15")
                plt.setp(text_box.patch, linestyle='', facecolor='lightgray', alpha=1.0, zorder=8)
                ax.add_artist(text_box)
        name = self.options["name"]
        space = self['spacing']
        if self['invis']:
            self.fig.patch.set_visible(False)
        if not space:
            plt.tight_layout()
        else:
            print(space)
            self.gs.update(wspace=space, hspace=space)
        plt.savefig(os.path.join(self.S.graph_dir, f'{name}.png'))
        plt.savefig(os.path.join(self.S.pgf_dir, f'{name}.pgf'))
