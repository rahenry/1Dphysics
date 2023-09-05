#!/usr/bin/env python3
import scipy.spatial
import numpy
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.markers
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection
from scipy.spatial import ConvexHull

def random_marker():
    angle = numpy.random.uniform(0,360)
    t = matplotlib.markers.MarkerStyle(marker='+')
    t._transform = t.get_transform().rotate_deg(angle)
    return t

known_coords = set()
class Polygon:
    def __init__(self, N, points=None, coords=None):
        self.N = N
        self.Q = 0
        self.d = 0
        if coords is not None:
            self.coords = coords
        else:
            self.coords = numpy.zeros(N)
        if points is not None:
            self.points = list(points)
        else:
            self.points = [numpy.array([0,0])]

    def expand(self, minimum, mode='ground'):
        N = self.N
        res = []
        for i in range(0,N):
            if i == 0 and mode == 'ground': continue
            coords = numpy.array(self.coords)
            coords[i] += 1
            sc = str(coords)
            if sc in known_coords:
                continue
            else:
                known_coords.add(sc)
            new = []
            unit = numpy.exp(2.j * numpy.pi * i / N)
            if mode == 'ground':
                unit = 1-unit
            # unit = 1-unit
            for p in self.points:
                z = numpy.array([unit.real, unit.imag])
                new.append(p + z*minimum)
                new.append(p + z)
            P = Polygon(self.N, new, coords)
            P.Q = (self.Q + i) % N
            P.d = self.d + 1

            P.hull()
            res.append(P)
        return res

    def __len__(self):
        return len(self.points)
    def hull(self):
        test = len([x for x in self.coords if x > 0])
        if test == 1:
            self.points = [self.points[0], self.points[-1]]
            return
        try:
            if len(self.points) < 3:
                return
            verts = ConvexHull(self.points).vertices
            self.points = [self.points[x] for x in verts]
        except scipy.spatial.qhull.QhullError:
            print(12312)
            pass

    def reorder(self):
        v = -1

        xmax = max([p[0] for p in self.points])
        ymax = max([p[1] for p in self.points])
        xmin = min([p[0] for p in self.points])
        ymin = min([p[1] for p in self.points])
        center = 0.5*numpy.array([xmax+xmin, ymax+ymin])

        angles = []
        for p in self.points:
            z = p - center
            a = z[0] + 1.j*z[1]
            theta = numpy.angle(a)
            angles.append([theta, abs(a)])
        angles.sort(key = lambda x: x[0])
        R = 10
        angle_sets = {}
        for theta, a in angles:
            if a < 1E-10 and len(angles) > 1:
                continue
            theta = numpy.around(theta, R)
            if not theta in angle_sets:
                angle_sets[theta] = [a]
            else:
                angle_sets[theta].append(a)
        res = []
        index = 0
        for theta, mags in angle_sets.items():
            x, y = center
            m = max(mags)
            X = m * numpy.cos(theta) + x
            Y = m * numpy.sin(theta) + y
            res.append(numpy.array([X, Y]))
            index += 1
        self.points = res

class Generator:
    def __init__(self, N, m, depth, Q=None, mode='ground', dexact=False, alpha=0.15):
        self.N = N
        self.m = m
        self.depth = depth
        self.mode = mode
        self.polygons = [Polygon(self.N)]
        old_polygons = self.polygons
        for i in range(depth):
            print(f'Depth = {i+1} / {depth}')
            new_polygons = []
            for o in old_polygons:
                new_polygons += o.expand(self.m, self.mode)
            self.polygons += new_polygons
            old_polygons = new_polygons

        # ptcs += [
        #     Wedge((.3, .7), .1, 0, 360),             # Full circle
        #     Wedge((.7, .8), .2, 0, 360, width=0.05),  # Full ring
        #     Wedge((.8, .3), .2, 0, 45),              # Full sector
        #     Wedge((.8, .3), .2, 45, 90, width=0.10),  # Ring sector
        # ]
        fig, ax = plt.subplots()
        for p in self.polygons:
            if Q is not None:
                if not p.Q == Q:
                    continue
            if dexact:
                if not p.d == depth:
                    continue
            xdata = [x[0] for x in p.points]
            ydata = [x[1] for x in p.points]
            # for x in p.points:
            #     ax.plot((x[0]), (x[1]), marker=random_marker(), markersize=10)
            if len(p) == 2:
                plt.plot(xdata, ydata, marker='+', markersize=1)
            elif len(p) == 1:
                plt.plot(xdata, ydata, marker='.', markersize=5)
            # else:
            #     plt.plot(xdata, ydata, marker='.', markersize=1, linestyle='')
        ptcs = []
        for p in self.polygons:
            if Q is not None:
                if not p.Q == Q:
                    continue
            if dexact:
                if not p.d == depth:
                    continue
            # for x in p.points:
            ptcs.append(matplotlib.patches.Polygon(p.points, True))
            # ptcs.append(matplotlib.patches.Polygon(numpy.random.rand(3, 2), True))

        colors = 100 * numpy.random.rand(len(ptcs))
        p = PatchCollection(ptcs, alpha=alpha)
        p.set_array(colors)
        ax.add_collection(p)
        # fig.colorbar(p, ax=ax)
        ax.autoscale_view()
        plt.show()


Q = 0
Q = 1
Q = 0
Q = 0
Q = None
Q = 2
N = 3
m = 0.5233
m = 0.433
m = 0.2
m = 0.3
d = 10
Q = None
Q = None
m = 0.36
ex = False
Q = None
md = 'ground'
alpha = 0.5
Generator(N, m, d, Q, md, False, alpha)
# Generator(N, m, d, Q, 'ground')

# z = numpy.random.rand(5, 2)
