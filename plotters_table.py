#!/usr/bin/env python3
#
from plotters import *
import itertools

class Table(Plotter):
    def get_2d_data(self):
        self.xvals = self.slicer.get_possible_values(self['xkey'])
        self.yvals = self.slicer.get_possible_values(self['ykey'])
        self.xlabs = [f'${self["xkey"]}={x}$' for x in self.xvals]
        self.ylabs = [f'${self["ykey"]}={y}$' for y in self.yvals]
        table_data = {}
        self.table_data = numpy.ndarray((len(self.xvals), len(self.yvals)), dtype=object)
        self.table_data.fill('-')
        tdd = {}
        for (x, y) in itertools.product(self.xvals, self.yvals):
            tdd[(x, y)] = '-'

        for s in self.slicer:
            x = s[self['xkey']]
            y = s[self['ykey']]
            tdd[(x, y)] = s[self['key']]

        for i in range(len(self.xvals)):
            for j in range(len(self.yvals)):
                x = self.xvals[i]
                y = self.yvals[j]
                self.table_data[i, j] = tdd[(x, y)]
    def stringy(self, x):
        exit()

    def draw(self):
        self.get_2d_data()
        # ax = self.grid.axs[self.grid.gs_offset]
        # ax.axis('off')
        # ax.axis('tight')
        # self.get_2d_data()
        # print(self.table_data)
        # table = ax.table(self.table_data, rowLabels=self.xlabs, colLabels=self.ylabs)

        shp = self.table_data.T.shape
        k1 = 3.2
        k2 = 4.5
        siz = (3+shp[0]/k1, shp[1]/k2+0.5)
        fig, ax = plt.subplots(figsize=siz)
        ax.axis('off')
        ax.axis('tight')
        c1 = 'lightgray'
        self.ny = self.table_data.shape[1]
        self.nx = self.table_data.shape[0]
        tabargs = {
            'colColours' : [c1]*self.ny,
            'rowColours' : [c1]*self.nx,
            'fontsize' : self.get('fsize', 8),
            }
        table = ax.table(self.table_data, rowLabels=self.xlabs, colLabels=self.ylabs, loc='center', **tabargs)
        S = self.slicer[0].parent
        f = os.path.join(S.graph_dir, f'table_{self["key"]}.pdf')
        fig.tight_layout()
        fig.savefig(f, bbox_inches='tight')


class Table2(Plotter):
    def get_2d_data(self):
        self.xvals = self.slicer.get_possible_values(self['xkey'])
        self.yvals = self.slicer.get_possible_values(self['ykey'])
        self.xlabs = [f'${self["xkey"]}={x}$' for x in self.xvals]
        self.ylabs = [f'${self["ykey"]}={y}$' for y in self.yvals]
        table_data = {}
        self.table_data = numpy.ndarray((len(self.xvals), len(self.yvals)), dtype=object)
        self.table_data.fill('-')
        tdd = {}
        for (x, y) in itertools.product(self.xvals, self.yvals):
            tdd[(x, y)] = '-'

        for s in self.slicer:
            x = s[self['xkey']]
            y = s[self['ykey']]
            g, a = s['n_zeros_alg'], s['n_zeros_geo']
            tdd[(x, y)] = f'{g}, {a}'

        for i in range(len(self.xvals)):
            for j in range(len(self.yvals)):
                x = self.xvals[i]
                y = self.yvals[j]
                self.table_data[i, j] = tdd[(x, y)]
    def stringy(self, x):
        exit()

    def post_table_process(self):
        1
    def draw(self):
        self.get_2d_data()
        shp = self.table_data.T.shape
        k1 = 3.2
        k2 = 4.5
        siz = (3+shp[0]/k1, shp[1]/k2+0.5)
        fig, ax = plt.subplots(figsize=siz)
        ax.axis('off')
        ax.axis('tight')
        c1 = 'lightgray'
        self.ny = self.table_data.shape[1]
        self.nx = self.table_data.shape[0]
        tabargs = {
            'colColours' : [c1]*self.ny,
            'rowColours' : [c1]*self.nx,
            'fontsize' : self.get('fsize', 8),
            }
        self.table = ax.table(self.table_data, rowLabels=self.xlabs, colLabels=self.ylabs, loc='center', **tabargs)
        self.post_table_process()
        fig.tight_layout()
        S = self.slicer[0].parent
        f = os.path.join(S.graph_dir, f'table2.pdf')
        fig.savefig(f, bbox_inches='tight')
        # table = ax.table(self.table_data, rowLabels=self.xlabs, colLabels=self.ylabs)
        # table.auto_set_font_size(False)
        # table.set_fontsize(5)
class ZeroTable(Plotter):
    def __init__(self, S, grid, info):
        self.rotate = 1
        super().__init__(S, grid, info)
    def post_table_process(self):
        1

    def get_2d_data(self):
        m = 0
        blocks = []
        ylabs = []
        for x, y in self.slicer.get_sets().items():
            head = labelify(x)
            s = y[0]
            for b, meas in s['block_measurements'].items():
                # if b not in blocks:
                #     blocks.append(b)
                if meas['n_zeros_alg'] > 0 and b not in blocks:
                    blocks.append(b)

        def sortfun(x):
            return x[1] + 1000 * x[0]
        blocks = sorted(blocks, key=sortfun)
        ny = len(blocks)
        nx = len(self.slicer.get_sets())
        self.ny = ny
        self.nx = nx
        I = 0
        t = numpy.empty((ny, nx), dtype=object)
        t.fill('-')
        xlabs = []
        for x0, y0 in self.slicer.get_sets(exclude=['L']).items():
            for x, y in y0.get_sets(exclude=['N']).items():
                # h = dict(x0)
                # h.update(dict(x))
                # head = labelify(h)
                s = y[0]
                head = f'{s["N"]}, {s["L"]}'
                xlabs.append(head)
                J = 0
                for b in blocks:
                    if not b in s['block_measurements']:
                        J += 1
                        continue
                    meas = s['block_measurements'][b]
                    ng = meas['n_zeros_geo']
                    na = meas['n_zeros_alg']
                    z = f'$n_a = {na}$, $n_g = {ng}$'
                    z = f'${na}$, ${ng}$'
                    t[J, I] = z
                    J += 1
                I += 1
        self.xlabs = xlabs
        for b in blocks:
            ylabs.append(f'{b}')
        self.ylabs = ylabs
        self.table_data = t

    def draw(self):
        self.get_2d_data()
        if self.rotate:
            self.table_data = self.table_data.T
        else:
            temp = self.ylabs
            self.ylabs = self.xlabs
            self.xlabs = temp
            temp = self.nx
            self.nx = self.ny
            self.ny = temp
        shp = self.table_data.shape
        k1 = 1.2
        k2 = 0.7
        A = 0.40
        siz = (0.1+shp[0]*k1, shp[1]*k2+0.1)
        siz = [x*A for x in siz]
        fig, ax = plt.subplots(figsize=siz)
        fig.subplots_adjust(0.1,0,1,0.9)
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')
        c1 = 'lightgray'
        tabargs = {
            'colColours' : [c1]*self.ny,
            'rowColours' : [c1]*self.nx,
            'fontsize' : self.get('fsize', 8),
            }
        table = ax.table(self.table_data, rowLabels=self.xlabs, colLabels=self.ylabs, loc='center', **tabargs)
        self.table = table
        self.post_table_process()
        table.auto_set_font_size(False)
        S = self.slicer[0].parent
        f = os.path.join(S.graph_dir, 'degen_table.pdf')
        # ax.set_title('Sector (charge, momentum)')
        ax.set_title('N, L')
        fig.suptitle('Sector (charge, momentum)', x=0, y=0.5, rotation=90, va='center', ha='center')
        # fig.suptitle('N, L', x=0.5, y=1.0, va='top', ha='center')

        # fig.tight_layout()
        fig.savefig(f, bbox_inches='tight')



# Rotated
class ZeroTableR(ZeroTable):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.rotate = 0

    def post_table_process(self):
        return
        for cell in self.table._cells:
            if cell[0] ==0:
                self.table._cells[cell].get_text().set_rotation(90)
                self.table._cells[cell].set_height(0.17)


class ZeroTablePlotly(ZeroTable):
    def __init__(self, S, grid, info):
        super().__init__(S, grid, info)
        self.rotate = 0

    def post_table_process(self):
        return
        for cell in self.table._cells:
            if cell[0] ==0:
                self.table._cells[cell].get_text().set_rotation(90)
                self.table._cells[cell].set_height(0.17)

    def draw(self):
        self.get_2d_data()
        print(self.table_data.shape)
        t = self.table_data.tolist()
        print(t)
        T = go.Table(cells=dict(values=self.table_data))
        fig = go.Figure(data=[T])
        fig.show()
        exit()
