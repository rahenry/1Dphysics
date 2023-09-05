#!/usr/bin/env python3

import random
import matplotlib.pyplot as plt
import math
import numpy
from exact_fp import *
from rah_utility import mkdir, matplotlib_set_aspect
base_dir = 'graphs/fp_generator'
mkdir(base_dir)

def to_bins(Nbins, points, scale=0.5):
    biggest = 0
    for p in points:
        m = p[0]**2+p[1]**2
        m = numpy.sqrt(m)
        m = max(p[0], p[1])
        if m > biggest: biggest = m
    biggest *= scale
    bins = numpy.zeros((Nbins, Nbins))
    for p in points:
        z = numpy.array([p[0], p[1]])
        # a = math.sqrt(z[0]**2+z[1]**2)
        # print(a, biggest)
        # if a > biggest: continue
        z = z / biggest * Nbins/2 + Nbins/2
        z = [round(x) for x in z]
        fail = 0
        for y in z:
            if y >= Nbins or y < 0:
                fail = 1
        if fail: continue
        bins[z[1]][z[0]] += 1
    return bins

class Generator:
    1
    def jump(self, origin):
        mag = random.uniform(self.min, self.max)
        up = random.choice([-1, 1])
        x = origin[0] + 0.5 * mag + mag
        y = origin[1] + up * math.sqrt(3) * 0.5 * mag
        q = origin[2] + up
        return (x, y, q)

    def add_path(self):
        p = (0, 0, 0)
        while (True):
            pnew = self.jump(p)
            self.points.append(pnew)
            x, y, q = pnew
            mag = x*x + y*y
            if x > self.maglim:
                break
            p = pnew

    def plot(self):
        for q in range(3):
            q = 2 - q
            points1 = [p for p in self.points if p[2] % 3 == q]
            xdata = [p[0] for p in points1]
            ydata = [p[1] for p in points1]
            plt.plot(xdata, ydata, linestyle='', marker='.', markersize=1)
        plt.show()

    def __init__(self, x0, x1, maglim=5):
        self.min = x0
        self.max = x1
        self.maglim = maglim
        self.points = []

        for i in range(10000):
            self.add_path()
        self.plot()


class Generator2:
    def jump(self, origin):
        mag = random.uniform(self.min, self.max)
        q = random.choice([0, 1, 2])
        x = origin[0] - mag * math.cos(q*2.*math.pi/3.)
        y = origin[1] - mag * math.sin(q*2.*math.pi/3.)
        q = origin[2] + q
        return (x, y, q)

    def add_path(self):
        p = (0, 0, 0)
        j = 0
        while (True):
            pnew = self.jump(p)
            # self.points.append(pnew)
            x, y, q = pnew
            mag = x*x + y*y
            j += 1
            if j > self.maglim:
                self.points.append(pnew)
                break
            p = pnew

    def plot(self):
        biggest = 0
        for p in self.points:
            m = p[0]**2+p[1]**2
            m = numpy.sqrt(m)
            if m > biggest: biggest = m
        scale = 0.5
        biggest *= scale

        # for q in range(3):
        #     points1 = [p for p in self.points if p[2] % 3 == q]
        #     points1 = [p for p in self.points]
        #     Nbins = 200
        #     bins = numpy.zeros((Nbins, Nbins))
        #     for p in points1:
        #         z = numpy.array([p[0], p[1]])
        #         # a = math.sqrt(z[0]**2+z[1]**2)
        #         # print(a, biggest)
        #         # if a > biggest: continue
        #         z = z / biggest * Nbins/2 + Nbins/2
        #         z = [round(x) for x in z]
        #         fail = 0
        #         for y in z:
        #             if y >= Nbins or y < 0:
        #                 fail = 1
        #         if fail: continue
        #         bins[z[1]][z[0]] += 1
        #     # bins = numpy.log(bins+1)
        #     # bins = numpy.log(bins+1)
        #     plt.imshow(bins)
        #     plt.show()
        #     exit()
        # exit()
        for q in range(3):
            q = 2 - q
            points1 = [p for p in self.points if p[2] % 3 == q]
            xdata = [p[0] for p in points1]
            ydata = [p[1] for p in points1]
            plt.plot(xdata, ydata, linestyle='', marker='.', markersize=1)
        plt.show()

    def __init__(self, x0, x1, maglim=5, npoints=10000):
        self.min = x0
        self.max = x1
        self.maglim = maglim
        self.points = []

        for i in range(npoints):
            self.add_path()
        self.plot()



# g = Generator(0.7, 1, 5)

class Generator3:
    def __init__(self, N, gamma, L):
        e, x, pfs = exact_parafermions_matrix(L, N, gamma)
        self.energies = [0]
        for p in pfs:
            energies_new = []
            for E in self.energies:
                for i in range(N):
                    e = E + numpy.exp(2.j*numpy.pi/N*i) * p
                    energies_new.append(e)
            self.energies = energies_new
        print(len(self.energies))
        xdata = [e.real for e in self.energies]
        ydata = [e.imag for e in self.energies]
        plt.plot(xdata, ydata, linestyle='', marker='.', markersize=1)
        plt.show()


class Generator4:
    def generate_point(self):
        res = 0
        for p in self.pfs:
            z = random.choice(range(self.N))
            res += p * numpy.exp(2.j*numpy.pi/self.N*z)
        return res

    def __init__(self, N, gamma, L, Npoints=1000):
        self.N = N
        self.L = L
        e0, x, self.pfs = exact_parafermions_matrix(L, N, gamma)
        self.points = []
        for p in range(Npoints):
            self.points.append(self.generate_point())
        xdata = [x.real for x in self.points]
        ydata = [y.imag for y in self.points]
        plt.plot(xdata, ydata, linestyle='', marker='.', markersize=1)
        plt.show()
        # pts = [[p.real, p.imag] for p in self.points]
        # Nbins = 100
        # bins = to_bins(Nbins, pts, scale=0.01)
        # plt.imshow(bins)
        # plt.show()

# g = Generator4(3, 0.58, 100, 100000)

class Generator5:
    def generate_point(self):
        res = 0
        for p in self.pfs:
            z = random.choice(range(self.N))
            res += p * numpy.exp(2.j*numpy.pi/self.N*z)
        return res

    def __init__(self, N, gamma, L, Npoints=1000, lim=0.1):
        self.N = N
        self.L = L
        e0, x, self.pfs = exact_parafermions_matrix(L, N, gamma)
        self.points = []
        for p in range(Npoints):
            self.points.append(self.generate_point())
        xdata = [x.real for x in self.points]
        ydata = [y.imag for y in self.points]
        plt.plot(xdata, ydata, linestyle='', marker='.', markersize=1)
        ax = plt.gca()
        l = [-lim,lim]
        ax.set_xlim(l)
        ax.set_ylim(l)
        plt.show()
        # pts = [[p.real, p.imag] for p in self.points]
        # Nbins = 100
        # bins = to_bins(Nbins, pts, scale=0.01)
        # plt.imshow(bins)
        # plt.show()

# g = Generator5(3, 0.88, 100, 100000, lim=3.1)

# g = Generator2(0.28, 1, 5, 100000)

class Generator7:
    def jump(self, origin):
        mag = random.uniform(self.min, self.max)
        q = random.choice([0, 1, 2])
        x = origin[0] - mag * math.cos(q*2.*math.pi/3.)
        y = origin[1] - mag * math.sin(q*2.*math.pi/3.)
        q = origin[2] + q
        return (x, y, q)

    def add_path(self):
        p = (0, 0, 0)
        j = 0
        while (True):
            pnew = self.jump(p)
            # self.points.append(pnew)
            x, y, q = pnew
            mag = x*x + y*y
            j += 1
            if j > self.maglim:
                self.points.append(pnew)
                break
            p = pnew

    def plot(self):
        biggest = 0
        for p in self.points:
            m = p[0]**2+p[1]**2
            m = numpy.sqrt(m)
            m = max(p[0], p[1])
            if m > biggest: biggest = m
        scale = 0.5
        biggest *= scale

        # for q in range(3):
        #     points1 = [p for p in self.points if p[2] % 3 == q]
        #     points1 = [p for p in self.points]
        #     Nbins = 200
        #     bins = numpy.zeros((Nbins, Nbins))
        #     for p in points1:
        #         z = numpy.array([p[0], p[1]])
        #         # a = math.sqrt(z[0]**2+z[1]**2)
        #         # print(a, biggest)
        #         # if a > biggest: continue
        #         z = z / biggest * Nbins/2 + Nbins/2
        #         z = [round(x) for x in z]
        #         fail = 0
        #         for y in z:
        #             if y >= Nbins or y < 0:
        #                 fail = 1
        #         if fail: continue
        #         bins[z[1]][z[0]] += 1
        #     # bins = numpy.log(bins+1)
        #     # bins = numpy.log(bins+1)
        #     plt.imshow(bins)
        #     plt.show()
        #     exit()
        # exit()
        for q in range(3):
            q = 2 - q
            q = 1
            points1 = [p for p in self.points if p[2] % 3 == q]
            xdata = [p[0] for p in points1]
            ydata = [p[1] for p in points1]
            # plt.plot(xdata, ydata, linestyle='', marker='.', markersize=1)
            Nbins = 100
            bins = to_bins(Nbins, points1)
            plt.imshow(bins)
            plt.show()
            exit()
        plt.show()

    def __init__(self, x0, x1, maglim=5, npoints=10000):
        self.min = x0
        self.max = x1
        self.maglim = maglim
        self.points = []

        for i in range(npoints):
            self.add_path()
        self.plot()

# g = Generator7(0.38, 1, 10, 100000)

class Generator8:
    def __init__(self, N, L, gamma, npoints=None, flipZ=0):
        self.N = N
        self.L = L
        self.gamma = gamma
        e, x, pfs = exact_parafermions_matrix(L, N, gamma, flipZ)
        # if flipZ:
        #     e = -e
        #     pfs = [-x for x in pfs]

        energies = make_pf_spectrum2(e, pfs, N, npoints)


        xdata = [x.real for x in energies]
        ydata = [x.imag for x in energies]
        plt.plot(xdata, ydata, linestyle='', marker='.')
        plt.show()
# g = Generator8(3, 100, 0.38, 3, 0)


class Generator9:
    def make_point(self):
        path = []
        i = 0
        while (True):
            x = random.randrange(0, self.L)
            if x in path:
                continue
            path.append(x)
            i += 1
            if i >= self.depth:
                break

        # print(path)
        pfs = [self.pfs[x] for x in path]
        # res = make_pf_spectrum3(self.g, pfs, self.N, None)
        res = make_pf_spectrum_rand(self.g, pfs, self.N, None)
        self.points += res

    def __init__(self, N, L, gamma, npaths=100, depth=3, flipZ=0, extend=0, setQ=None):
        args = locals()
        self.N = N
        self.L = L
        self.gamma = gamma
        self.depth = depth
        print(0)
        e, x, pfs = exact_parafermions_matrix(L, N, gamma, flipZ)
        pfs.sort()
        pfs.reverse()
        print(2)
        self.pfs = pfs
        self.g = e
        self.npaths = 0
        self.points = []
        if depth > 0:
            for i in range(npaths):
                self.make_point()

        print(1)
        energies = self.points
        print(1)

        if extend:
            energies_extended = []
            for i in range(0, N):
                energies_extended += [e * numpy.exp(2.j*numpy.pi*i/N) for e in energies]
            energies = energies_extended
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.set_xlabel("Re(E)")
        ax1.set_ylabel("Im(E)")
        ax2.set_xlabel("Re(E)")
        ax2.set_ylabel("Im(E)")
        ax1.set_title("Parafermion energies")
        ax2.set_title("Eigenstate energies (sampled)")
        print(1)
        for Q in range(N):
            if setQ is not None: Q = setQ
            xdata = [x[0].real for x in energies if x[1] == Q]
            ydata = [x[0].imag for x in energies if x[1] == Q]
            ax2.plot(xdata, ydata, linestyle='', marker='o', markersize=0.5)
            if setQ is not None: break

        xdata = [p.real for p in pfs]
        ydata = [p.imag for p in pfs]
        ax1.plot(xdata, ydata, linestyle='', marker='o', markersize=0.5)
        name = ''
        skip = 'extend'
        figtext = ''
        for x, y in args.items():
            if x == 'self': continue
            name += f'{x}={y},'
            figtext += f'{x} = {y}\n'
        lamd = (1.0-gamma)/gamma
        if flipZ == 1: lamd *= -1
        figtext += f'lambda = {lamd}'
        fig.text(0.01, 0.99, figtext, va='top')

        name = name.strip(',')
        matplotlib_set_aspect(ax1, 1)
        matplotlib_set_aspect(ax2, 1)
        plt.savefig("/home/rah/graph.png")
        fig.tight_layout()
        plt.show()

gamma = 0.46
gamma = 0.7066
gamma = 0.5636
flipZ = 0
setQ = None
npoints = 10000
L = 250
depth = 6
N = 3
setQ = 0
g = Generator9(N, L, gamma, npoints, depth, flipZ, 0, setQ)

e, x, pfs = exact_parafermions_matrix(100, N, gamma, flipZ)
pfs.sort()
print(pfs[0] / pfs[-1])
