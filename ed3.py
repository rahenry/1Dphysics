#!/usr/bin/env python3
import numpy, scipy, math, cmath, random, scipy.sparse, scipy.sparse.linalg, scipy.special, os, shutil
from rah_numpy import eigenvalue_test
from rah_utility import rel_diff
import matplotlib.pyplot as plt
import matplotlib
import cProfile, pstats, io
from pstats import SortKey

tola = 1E-5
tolg = 1E-5
class block:
    def __init__(self, parent, Q, K):
        self.parent = parent
        self.L = parent.L
        self.N = parent.N
        self.K = K
        self.Q = Q
        self.basis = []
        self.measurements = {}

    def degen_test(self):
        db = {}
        tol = tola
        # tol = numpy.sqrt(tol1)
        etest = 1E100
        itest = None
        self.zero_block = None
        self.zero_block_left = None
        for i in range(len(self)):
            e = self.es[i]
            if rel_diff(e, etest) < tol:
                db[itest].append(i)
            elif (abs(e) < tol and abs(etest) < tol):
                db[itest].append(i)
            else:
                db[i] = [i]
                etest = e
                itest = i
                if abs(e) < tol:
                    self.zero_block = i

            eL = self.esL[i]
            if abs(eL) < tol and self.zero_block_left is None:
                self.zero_block_left = i


        self.degen_blocks = db

        non_zero = []
        nz = 0
        for i, states in db.items():
            if len(states) > 1:
                e = self.es[i]
                psis = self.vs[:,states]
                orth = scipy.linalg.orth(psis, tolg)
                north = orth.shape[1]
                if abs(e) > 1E-6:
                    nz += north
                # print(psis.shape, orth.shape)
                # print(i, self.zero_block, len(states), north)
                if len(states) == 2:
                    o = numpy.vdot(self.vs[:,i], self.vs[:,i+1])
                    # print(o)
                    if north == 1:
                        if abs(o) - 1 < 1E-6:
                            non_zero.append(e)
                        else:
                            print(f'Problem with state {i} {e}')

        if self.zero_block is not None:
            zstates = db[self.zero_block]
            psis = self.vs[:,zstates]
            orth = scipy.linalg.orth(psis, tolg)
            north = orth.shape[1]
            self.nzeros = north
        else:
            self.nzeros = 0


        self.measurements['non_zero_degens'] = non_zero
        self.measurements['n_degens_nonzero'] = nz

                # for s1 in states:
                #     for s2 in states:
                #         psiL = self.vs[:,s1]
                #         o = numpy.dot(psiL, psiR)

        #         print(states)
        #         print(self.es[states])

    def measure(self):

        if self.parent.sys['measure_matrices']:
            sref = self.parent.sys.parent.ref_system
            edref = sref.ed3
            bref = edref.blocks[(self.Q, self.K)]
            matrix_self = self.matrix_measure(self)
            self.measurements.update(matrix_self)
            matrix_ref = self.matrix_measure(bref, '_ref')
            self.measurements.update(matrix_ref)
            bref.degen_test()

        # if self.zero_block:
        #     m = len(self.degen_blocks[self.zero_block])

        # nzd = self.measurements['non_zero_degens']
        # if len(nzd) > 0:
        #     print(self, nzd)
        if self.zero_block is not None:
            # a = len(self.degen_blocks[self.zero_block])
            a = self.nzeros
            self.measurements['n_zeros_alg'] = a
            g = scipy.linalg.null_space(self.H, tolg).shape[1]

            # ns = scipy.linalg.null_space(self.H, tolg)
            if len(self) == a:
                g = a
            u, s, v = scipy.linalg.svd(self.H)
            count = 0
            for s in s:
                if abs(s) < tolg:
                    count += 1

            self.measurements['n_zeros_geo'] = count
            # for i in self.degen_blocks[self.zero_block]:
            #     v = self.vs[:,i]
            #     u = self.H @ v

        else:
            self.measurements['n_zeros_alg'] = 0
            self.measurements['n_zeros_geo'] = 0
        #     z = numpy.full((m, m), 0.)
        #     zd = numpy.full((m, m), 0.j)
        #     i = 0
        #     for i in self.degen_blocks[self.zero_block]:
        #         v = self.vs[:,i]
        #         print(i, numpy.vdot(v, v))
        #         vP = bref.H @ v
        #         print(numpy.vdot(vP, vP))

        # self.block_det(bref)
        # nnull = self.null.shape[1]
        # orth = scipy.linalg.orth(self.H, tol1)
        # north = orth.shape[1]
        # n = len(self)
        # nnull2 = 0
        # if nnull > 0:
        #     print(self, nnull)
        # if self.null.shape[1] > 0:
        #     null2 = numpy.concatenate([self.null, bref.null], axis=1)
        #     null2 = scipy.linalg.orth(null2, tol1)
        #     nnull2 = null2.shape[1]
        #     print(len(self.degen_blocks[self.zero_block]), nnull, north, n, nnull2)

        # print(nnull, north, n)
        # if self.null.shape[1] > 0:
        #     M = numpy.full((n, n), 0.j)
        #     M2 = numpy.full((nnull, nnull), 0.j)
        #     for iL in range(nnull):
        #         for iR in range(nnull):
        #             vL = self.nullL[:,iL]
        #             vR = self.null[:,iR]
        #             o = numpy.dot(vL, vR)
        #             M2[iL, iR] = abs(o)
                    # print(iL, iR)
                    # print(o)

            # print(M2)
            # for iL in range(nnull):
            #     vL = self.nullL[:,iL]
            #     vR = self.null[:,iL]
            #     M += numpy.outer(vR, vL)

            # d = numpy.linalg.det(M)
            # d = abs(d)
            # print(d, self)

        # self.measurements['n_zeros'] = 0
        # if self.zero_block:
        #     m = len(self.degen_blocks[self.zero_block])
        #     z = numpy.full((m, m), 0.)
        #     zd = numpy.full((m, m), 0.j)
        #     i = 0
        #     for i in self.degen_blocks[self.zero_block]:
        #         J = 0
        #         for j in bref.degen_blocks[bref.zero_block]:
        #             o = numpy.vdot(bref.vs[:,j], self.vs[:, i])
        #             z[I, J] = abs(o)
        #             zd[I, J] = o
        #             J += 1
        #         I += 1
        #     self.measurements['zero_info'] = z
        #     D = numpy.linalg.det(zd)
        #     # print(D, self)
        #     # print(zd)
        #     self.measurements['zero_det'] = numpy.linalg.det(zd)
        #     self.measurements['n_zeros'] = m

            # m = len(self.degen_blocks[self.zero_block])
            # zd = numpy.full((m, m), 0.j)
            # zdR = numpy.full((m, m), 0.j)
            # for I in range(m):
            #     for J in range(m):
            #         iR = self.degen_blocks[self.zero_block][0] + I
            #         iL = self.zero_block_left + J
            #         iLR = bref.zero_block_left + J
            #         o = numpy.vdot(bref.vsL[:,iLR], self.vs[:, iR])
            #         zdR[I, J] = o
            #         o = numpy.dot(self.vsL[:,iL], self.vs[:, iR])
            #         zd[I, J] = o
            #         J += 1
            #     I += 1
            # print(self)
            # print(zd)
            # D = numpy.linalg.det(zd)
            # print(D)
            # print(zd)
            # D = numpy.linalg.det(zd)
            # print(D)

    def matrix_measure(self, other, label=''):
        res = {}

        n = len(self)
        overlaps = numpy.empty((n, n), dtype=float)
        for i in range(n):
            for j in range(n):
                psiL = other.vsL[:,i]
                psiR = self.vs[:,j]
                o = numpy.vdot(psiL, psiR)
                overlaps[i, j] = abs(o)
        res['overlapLR'+label] = overlaps

        overlaps = numpy.empty((n, n), dtype=float)
        for i in range(n):
            for j in range(n):
                psiL = other.vs[:,i]
                psiR = self.vs[:,j]
                o = numpy.vdot(psiL, psiR)
                overlaps[i, j] = abs(o)
        res['overlapRR'+label] = overlaps
        return res


    def solve(self):
        # sol = scipy.linalg.eig(self.H)
        sol = numpy.linalg.eig(self.H)
        inds = numpy.argsort(sol[0])
        self.es = sol[0][inds]
        self.vs = sol[1][:,inds]

        sol = scipy.linalg.eig(numpy.conj(self.H.T))
        inds = numpy.argsort(sol[0])
        self.esL = sol[0][inds]
        self.vsL = sol[1][:,inds]

        # nulltol = tolg
        # self.null = scipy.linalg.null_space(self.H, nulltol)
        # self.nullL = scipy.linalg.null_space(self.H.T, nulltol)

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

    def __str__(self):
        return f"Block q={self.Q}, k={self.K}, nstates={len(self.basis)}"

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.basis)

    def add_state(self, v):
        self.basis.append(v)

    def make_H(self):
        n = len(self.basis)
        self.H = numpy.full((n, n), 0.j)
        for i in range(n):
            for j in range(n):
                # self.H[i,j] = self.Hbasis(i, j)
                self.H[i,j] = self.Hbasis2(i, j)


    def Helem(self, u, v):
        res = 0
        cX = self.parent.sys.model.cZ
        cZ = self.parent.sys.model.cX
        cX = 0
        # cZ = 0
        if u == v:
            for i in range(self.L):
                # res += cX * self.parent.Z[i, u]
                1
                res += cX * numpy.exp (2.j * numpy.pi / self.N * self.get_site(u, i))

        for i in range(self.L):
            uP = self.mod_site(u, i, (self.get_site(u, i) + 1) % self.N)
            j = (i + 1) % self.L
            uP = self.mod_site(uP, j, (self.get_site(uP, j) - 1) % self.N)
            # uP = self.parent.XX[i, u]
            if uP == v:
                res += cZ
                1

        return res

    def Hbasis(self, i, j):
        s1, sc1 = self.basis[i]
        s2, sc2 = self.basis[j]
        res = 0
        for I in range(s1.shape[0]):
            for J in range(s2.shape[0]):
                x1, a1 = s1[I], sc1[I]
                x2, a2 = s2[J], sc2[J]
                a1 = numpy.conj(a1)
                h = self.Helem(x1, x2)
                res += a1 * a2 * h

        return res

    def Hbasis2(self, i, j):
        res = 0
        s1, sc1 = self.basis[i]
        s2, sc2 = self.basis[j]
        cX = self.parent.sys.model.cZ
        cZ = self.parent.sys.model.cX
        # if s1.shape[0] == 1 or s2.shape[0] == 1:
        #     return 0
        if i == j:
            Z = self.parent.Z[:,s1]
            a = Z * sc1[None,:]
            b = numpy.einsum('j,ij', sc2.conj(), a)
            b = numpy.sum(b)
            res += cX * b
        # X = numpy.subtract.outer(s2, s1)
        # X = numpy.isin(X, self.parent.Xmask)
        # X = sc2[:,None].conj() * X * sc1
        # # X = sc2.conj() * X * sc1[:,None]
        # X = numpy.sum(X) * cZ
        # res += X
        # res = 0
        # res += self.Hbasis(i, j)
        X = self.parent.XX[:,s1]
        X = numpy.equal(s2[:,None,None], X)
        X = sc2.conj()[:,None,None] * X * sc1
        res += cZ * numpy.sum(X)
        return res

class ed3(block):
    def __str__(self):
        return "ed3"
    def __init__(self, s):
        self.L = s['L']
        self.sys = s
        self.N = s['N']
        self.n = self.N**self.L
        s.ed_system = self

        self.precalc()
        self.make_blocks()
        self.make_H()
        self.solve()

    def precalc(self):
        self.Z = numpy.full((self.L, self.n), 0.j)
        self.X = numpy.full((self.L, self.n), 0)
        self.XX = numpy.full((self.L, self.n), 0)
        self.Xmask = numpy.array([(self.N-1) * (self.N ** x) for x in range(1,self.L)]+[1-self.N**(self.L-1)])
        for i in range(self.L):
            for u in range(self.n):
                uP = self.mod_site(u, i, (self.get_site(u, i) + 1) % self.N)
                self.X[i, u] = uP
                j = (i + 1) % self.L
                uP = self.mod_site(uP, j, (self.get_site(uP, j) - 1) % self.N)
                self.XX[i, u] = uP
                self.Z[i, u] = numpy.exp(2.j * numpy.pi / self.N * self.get_site(u, i))

    def in_Xmask(self, u, v):
        a = u - v
        return numpy.isin(a, self.Xmask)

    def make_blocks(self):
        used = set()
        self.blocks = {}
        for Q in range(self.N):
            for K in range(self.L):
                b = block(self, Q, K)
                self.blocks[(Q, K)] = b
        for u in range(self.n):
            if u in used:
                continue
            p = self.period(u)
            q = self.charge(u)
            # for each momentum quantum number
            for j in range(p):
                s = numpy.full(p, 0)
                sc = numpy.full(p, 0.j)
                k = j * (self.L // p)
                x = u
                # add each component of the state and its coefficient
                for i in range(p):
                    a = numpy.exp(2.0j * numpy.pi / self.L * i * k) / numpy.sqrt(p)
                    s[i] = x
                    sc[i] = a
                    x = self.translate(x)
                    used.add(x)

                self.blocks[(q, k)].add_state((s, sc))

    def make_H(self):
        t = 0
        for (Q, K), b in self.blocks.items():
            t += len(b.basis)
            b.make_H()


    def solve(self):
        for (q, k), b in self.blocks.items():
            b.solve()
        es = numpy.concatenate([b.es for (q, k), b in self.blocks.items()])
        es = es.flatten()
        print(es.shape)

        s = self.sys
        s['energies'] = es

        z = numpy.empty((self.N, self.L), dtype=object)
        for (q, k), b in self.blocks.items():
            z[q,k] = numpy.array(b.es)
        z = {}
        for (q, k), b in self.blocks.items():
            z[(q,k)] = numpy.array(b.es)
        s.measurement_data['block_energies'] = z

    def measure(self):
        z = {}
        for (q, k), b in self.blocks.items():
            b.degen_test()
        for (q, k), b in self.blocks.items():
            b.measure()
            z[(q,k)] = b.measurements
        self.sys.measurement_data['block_measurements'] = z
        a = 0
        g = 0
        t = 0
        o = 0
        for (q, k), b in self.blocks.items():
            a += b.measurements['n_zeros_alg']
            g += b.measurements['n_zeros_geo']
            o += b.measurements['n_degens_nonzero']
            # t += b.measurements['n_zeros_alg']
            t += b.measurements['n_zeros_geo']
            t += b.measurements['n_degens_nonzero']
        self.sys.measurement_data['n_zeros_alg'] = a
        self.sys.measurement_data['n_zeros_geo'] = g
        self.sys.measurement_data['n_degens_nonzero'] = o
        self.sys.measurement_data['n_degens_total'] = t


def solve_ed3(s):
    pro = 0
    if pro:
        pr = cProfile.Profile()
        pr.enable()
        ed3(s)
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        # ps = pstats.Stats(pr, stream=s).sort_stats(sortby).reverse_order()
        ps.print_stats()
        print(s.getvalue())
        exit()

    s.ed3 = ed3(s)
    s.measurement_data['converged_0'] = 1
    s.save_measurements()
