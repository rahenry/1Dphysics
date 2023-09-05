#!/usr/bin/env python3


import numpy
import random
from exact_fp import exact_parafermions_matrix, make_pf_spectrum2, rhoX, make_coords

class PointGenerator:
    def generate_point(self):
        res = 0
        for p in self.s['pfs']:
            z = random.choice(range(self.N))
            res -= p * numpy.exp(2.j*numpy.pi/self.N*z)
        return res

    def generate_point_gs(self, x):
        1

    def generate_points_gs(self):
        n = self.s['n_points_generated']
        if n == 'all':
            return
        depth = self.s['generate_gs_depth']
        for i in range(depth+1):
            points = make_coords(self.L, i)
            for p in points:
                # print(f'p = {p}')
                to_add = make_coords(self.N, len(p))
                for x in to_add:
                    e = self.e0
                    for k in range(len(x)):
                        y = x[k]
                        # print(x, k, y)
                        z = (1. - numpy.exp(-2.j * numpy.pi / self.N * y))
                        e += z * self.s['pfs'][k]
                    self.add_point(e)

    def __init__(self, s):
        self.s = s
        self.N = s['N']
        self.L = s['L']
        self.points = set([])
        if s['energies_generated']:
            self.points = s['energies_generated']
        self.e0 = 0
        for p in s['pfs']:
            self.e0 -= p

    def add_point(self, p):
        for j in range(self.N):
            u = p * numpy.exp(numpy.pi * -2.j * j / self.N)
            self.points.add(u)

    def get_data(self):
        xdata = [x.real for x in self.points]
        ydata = [y.imag for y in self.points]
        return self.points, xdata, ydata

    def make_points(self):
        n = self.s['n_points_generated']
        if n == 'all':
            self.points = set([])
            spect_poly = make_pf_spectrum2(self.e0, self.s['pfs'], self.s['N'])
            for p in spect_poly:
                self.points.add(p)
        else:
            for i in range(n):
                p = self.generate_point()
                self.add_point(p)
