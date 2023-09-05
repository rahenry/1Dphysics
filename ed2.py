#!/usr/bin/env python3
import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, os, shutil
from rah_numpy import eigenvalue_test
from rah_utility import rel_diff
import matplotlib.pyplot as plt
import matplotlib


class block:
    def __init__(self, parent, q, k):
        self.parent = parent
        self.L = parent.L
        self.N = parent.N
        self.q = q
        self.k = k
        self.basis = []

    def solve2(self):
        n = len(self)
        b = numpy.array(self.basis).T
        self.m = b.conj().T @ self.parent.H @ b
        sol = scipy.linalg.eig(self.m.T)
        self.energies = sol[0]
        self.zeros = []
        self.zero_evecs = []
        for i, e in enumerate(sol[0]):
            if abs(e) < 1E-6:
                v = sol[1][:,i]
                self.zeros.append(e)
                self.zero_evecs.append(v)
        self.evecs_basis = sol[1]
    # def solve(self, i):
    #     n = len(self)
    #     self.m = numpy.full((n, n), 0.j)
    #     for j in range(self.n):
    #         1

    #     M = self.parent.THT
    #     self.m = M[i : i + n, i : i + n]
    #     sol = scipy.linalg.eig(self.m.T)
    #     self.energies = sol[0]
    #     # print(sol[0])
    #     self.evecs_basis = sol[1]

    # def __hash__(self):
    #     return hash(f'{self.L}_{self.N}_{self.q}_{self.k}')

    # def __eq__(self, other):
    #     return self.__hash__() == other.__hash__()
    #
    #
    def __str__(self):
        return f"Block q={self.q}, k={self.k}, nstates={len(self.basis)}"

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.basis)

    def add_state(self, v):
        self.basis.append(v)
        # print(f'{self} added {v}')


class ed_system:
    def __init__(self, s):
        self.L = s["L"]
        self.N = s["N"]
        self.n = self.N**self.L
        self.sys = s
        s.ed_system = self

        self.blocks = {}
        self.make_H()
        self.make_H_swapped()
        self.make_basis()
        # print(self.H)
        # print(self.T)
        self.solve()
        # self.test1()
        #
    def solve(self):
        i = 0
        # self.energies = numpy.array([])
        es = []
        for (q, k), B in self.blocks.items():
            # print(i)
            # B.solve(i)
            B.solve2()
            i += len(B)
            # self.energies.concatenate(B.energies)
            for x in B.energies:
                es.append(x)
        es = sorted(es, key=lambda x:x.real)
        self.sys['energies'] = es
        self.sys['E_0'] = es[0]
        self.sys['e_0'] = es[0] / self.L
        for B in self.blocks:
            print(B)
        exit()

    def make_H_swapped(self):
        self.H = numpy.full((self.n, self.n), 0.0j)
        cX = self.sys.model.cZ
        cZ = self.sys.model.cX

        for i in range(self.L):
            for u in range(self.n):
                self.H[u, u] += cX * numpy.exp(
                    2.0j * numpy.pi / self.N * self.get_site(u, i)
                )

        lim = self.L - 1
        if self.sys['bc'] == 1:
            lim = self.L
        for i in range(lim):
            for u in range(self.n):
                uP = self.mod_site(u, i, (self.get_site(u, i) + 1) % self.N)
                j = (i + 1) % self.L
                uP = self.mod_site(uP, j, (self.get_site(uP, j) - 1) % self.N)
                self.H[uP, u] += cZ
        q = abs(self.H)
        sol = scipy.linalg.eig(self.H)
        es = sorted(sol[0], key=lambda x:x.real)
        # for e in es:
        #     print(e / self.L)
        # print(q)

    def make_H(self):
        self.H = numpy.full((self.n, self.n), 0.0j)
        cX = self.sys.model.cZ
        cZ = self.sys.model.cX

        for i in range(self.L - 1):
            for u in range(self.n):
                delta = self.get_site(u, i) - self.get_site(u, i + 1)
                self.H[u][u] += cZ * numpy.exp(2.0j * numpy.pi / self.N * delta)

        for i in range(self.L):
            for u in range(self.n):
                uP = self.mod_site(u, i, (self.get_site(u, i) + 1) % self.N)
                self.H[u][uP] += cX
        if self.sys["bc"] == 1:
            for u in range(self.n):
                delta = self.get_site(u, self.L - 1) - self.get_site(u, 0)
                self.H[u][u] += cZ * numpy.exp(2.0j * numpy.pi / self.N * delta)
        q = abs(self.H)
        sol = scipy.linalg.eig(self.H)
        es = sorted(sol[0], key=lambda x:x.real)
        # print(cX, cZ)
        # for e in es:
        #     print(e / self.L)
        # print(q)

    def make_basis(self):
        for q in range(self.N):
            for k in range(self.L):
                self.blocks[(q, k)] = block(self, q, k)

        n = self.n
        used = set()
        for u in range(n):
            if u in used:
                continue
            p = self.period(u)
            q = self.charge(u)
            for j in range(p):
                v = numpy.full((n), 0.0j)
                k = j * (self.L // p)
                x = u
                for i in range(p):
                    v[x] = (
                        1 * numpy.exp(2.0j * numpy.pi / self.L * i * k) / numpy.sqrt(p)
                    )
                    x = self.translate(x)
                    used.add(x)

                self.blocks[(q, k)].add_state(v)

        self.T = numpy.full((n, n), 0.0j)
        i = 0

        o = numpy.full((n, n), 0.0)

        for (q, k), B in self.blocks.items():
            for v in B.basis:
                self.T[i] = v

                j = 0
                for (q2, k2), B2 in self.blocks.items():
                    for v2 in B2.basis:
                        z = numpy.dot(v.conj(), v2)
                        if abs(z) > 1e-6 and not i == j:
                            print(v)
                            print(v2)
                            print(z)
                            print("...")
                        o[i][j] = numpy.dot(v.conj(), v2)
                        j += 1

                i += 1



        # self.THT = self.T @ self.H @ (self.T.T.conj())


        # m = numpy.abs(o) > 1E-5
        # m = numpy.where(m, 1, 0)
        # self.THT = self.T.T @ self.H @ (self.T)
        # self.THT = self.T @ self.H @ self.T.T
        # self.THT = numpy.conj(self.T.T) @ self.H @ (self.T)
        # def thing(x):
        #     return abs(x) > 0
        # m = numpy.abs(self.THT) > 1E-5
        # m = numpy.where(m, 1, 0)
        # print(m)

        # print(numpy.abs(self.THT))

    def get_site(self, state, site):
        j = site
        i = state
        return (i // (self.N**j)) % self.N

    def mod_site(self, initial, j_modified, sigma_j_f):
        sigma_j_i = self.get_site(initial, j_modified)
        return initial + ((self.N**j_modified) * (sigma_j_f - sigma_j_i))

    def translate(self, u):
        t = 0
        for j in range(self.L):
            s = self.get_site(u, j)
            t = self.mod_site(t, (j + 1) % self.L, s)
        return t

    def period(self, z):
        z = self.stringify(z)
        for i in range(1, len(z) + 1):
            zp = z[i:] + z[:i]
            if zp == z:
                return i

    def stringify(self, u):
        s = ""
        for j in range(self.L):
            s += str(self.get_site(u, j))
        return s

    def charge(self, u):
        q = 0
        for j in range(self.L):
            q += self.get_site(u, j)
        return q % self.N

    def walls(self, u):
        q = 0
        for j in range(self.L):
            a = self.get_site(u, j)
            b = self.get_site(u, (j + 1) % self.L)
            q += a - b
        return q % self.N

    def test1(self):
        print(123)
        print(self.N, self.L)
        u = 5
        for u in range(self.n):
            for k in range(self.L):
                s0 = self.get_site(u, k)
                s1 = (s0-1)%self.N
                uP = self.mod_site(u, k, s1)
                print(u, k, s0, s1, self.get_site(uP, k))
        exit()
def test_energies(s):
    c = dict(s.config_base)
    c['method'] = 'simple_ed'
    from System import System
    stest = System(c, just_go=1)
    es = sorted(stest['energies'], key=lambda x:x.real)
    es2 = sorted(s['energies'], key=lambda x:x.real)
    for i, e in enumerate(es):
        print(i, e, es2[i])

def solve_ed2(s):
    s.ed = ed_system(s)
    # test_energies(s)
