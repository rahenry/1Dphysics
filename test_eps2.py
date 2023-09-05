#!/usr/bin/env python3


from System import *
from SystemSet import *
from args import Q
from exact_fp import *
from rah_utility import *
from output_figs import *

import numpy, scipy

def find_minima(z):
    X, Y = z.shape
    res = []
    for i in range(X-2):
        x = i+1
        for j in range(Y-2):
            y = j+1
            u = z[x][y+1]
            d = z[x][y-1]
            l = z[x-1][y]
            r = z[x+1][y]
            t = z[x][y]
            things = [-1, 0, 1]
            accepted = 1
            for ox in things:
                for oy in things:
                    if oy == ox == 0:
                        continue
                    if t > z[x+ox][y+oy]:
                        accepted = 0

            if accepted:
                res.append((x, y))
            # if t < u and t < d and t < l and t < r:
            #     res.append((x, y))

    return res
z = numpy.random.rand(5,5)
print(z)
ex = find_minima(z)
print(ex)


Ls = [3, 4, 5, 6]
Ns = [2, 3]
Ns = [3]
g0 = 0.3
g1 = 0.5
p0 = 0.
p1 = 1.
npoints = 101

Ls = [51]
g0 = 0.3
g1 = 0.6
p0 = 0.
p1 = 1.0
npoints = 151


Ls = [5]
npoints = 51

Ls = [100]
npoints = 100

gammas = numpy.linspace(g0, g1, npoints)
phis = numpy.linspace(p0, p1, npoints)

def get_pfs(N, L, gamma, phi):
    f = os.path.join('pfdata', f'{N}_{L}_{phi}_{gamma}')
    pfs = None
    if os.path.exists(f) and not Q.rerun_measurements:
        with open(f, 'rb') as f1:
            pfs = dill.load(f1)
    else:
        lamb = gamma / (1.-gamma)
        lsc = 1. / (1.-gamma)
        s = {
            'L' : L,
            'gamma' : gamma,
            'phi':  phi,
            'N' : N,
            'lambda' : lamb,
            'lambda_scale' : lsc,
            'fix_pfs' : 1,
            'angle_style' : 0,
            }
        pfs, e0 = pfs_matrix(s)
        pfs = numpy.array(pfs)
        with open(f, 'wb+') as f1:
            dill.dump(pfs, f1)
            # print("saving")
    return pfs

def min_func(N, L, gamma, phi):
    pfs = get_pfs(N, L, gamma, phi)
    D = abs(pfs[0] - pfs[1])
    D = abs(pfs[1])
    return D

def min_func2(N, L, gamma, phi):
    pfs = get_pfs(N, L, gamma, phi)
    D = abs(pfs[0] - pfs[1])
    # D = abs(pfs[0]) - abs(pfs[1])
    return D

mkdir('pfdata')
fbase = 'graphs/test_eps'
silent_remove(fbase)
mkdir('graphs/test_eps')
for N in Ns:
    for L in Ls:
        print('............')
        print(N, L)
        M = numpy.full((len(gammas), len(phis), L-1), 0.)
        for i in range(npoints):
            for j in range(npoints):
                phi = phis[j]
                gamma = gammas[i]
                pfs = get_pfs(N, L, gamma, phi)
                for k in range(L-1):
                    # M[i,j,k] = abs(pfs[k+1]-pfs[k])
                    M[i,j,k] = abs(pfs[k])

        m = M[:, :, 1]
        f = os.path.join('graphs/test_eps', f'{N}_{L}_{k}.png')
        ex = [p0, p1, g0, g1]
        aspect = abs(ex[1]-ex[0])/abs(ex[3]-ex[2]) * 1
        plt.clf()
        p = plt.imshow(m, extent=ex, aspect=aspect, origin='lower')
        plt.colorbar(p)

        # plt.imshow(m)
        plt.ylabel(r'$\gamma$')
        plt.xlabel(r'$\phi$')

        plt.savefig(f)
        ms = find_minima(m)

        def mf(x):
            return min_func(N, L, x[0], x[1])


        ms2 = []
        out = []
        for m in ms:
            phi = phis[m[1]]
            gamma = gammas[m[0]]
            out.append(f'm1 = {phi} {gamma}')
            plt.scatter(phi, gamma, marker='x', c='red')

            # bounds = ((g0+0.01, g1-0.01), (p0+0.01, p1-0.01))
            #
            bounds = ((g0, g1), (p0, p1))
            dphi = 0.1*phi
            bounds = ((g0, g1), (phi-dphi, phi+dphi))
            q = scipy.optimize.minimize(mf, (gamma, phi), method="Powell", bounds=bounds)
            plt.scatter(q.x[1], q.x[0], marker='o', c='green')
            ms2.append(q.x)
        plt.savefig(f)

        plt.clf()

        M2 = numpy.abs(M[:,:,1] - M[:,:,0])
        f = os.path.join('graphs/test_eps', f'M2_{N}_{L}_{k}.png')
        ex = [p0, p1, g0, g1]
        aspect = abs(ex[1]-ex[0])/abs(ex[3]-ex[2]) * 1

        p = plt.imshow(M2, extent=ex, aspect=aspect, origin='lower')

        plt.colorbar(p)
        plt.ylabel(r'$\gamma$')
        plt.xlabel(r'$\phi$')

        def mf(x):
            return min_func2(N, L, x[0], x[1])
        ms3 = []
        for m in ms2:
            gamma, phi = m
            bounds = ((g0, g1), (p0, p1))
            dphi = 0.1*phi
            bounds = ((g0, g1), (phi-dphi, phi+dphi))
            # ops = {'xatol' : 1E-10, 'fatol' : 1E-10}
            ops = {'eps' : 1E-10}
            method = "L-BFGS-B"
            q = scipy.optimize.minimize(mf, (gamma, phi), bounds=bounds, tol=1E-10, options=ops, method=method)
            plt.scatter(q.x[1], q.x[0], marker='x', c='green')
            # print(gamma, phi, q.x)
            ms3.append(q.x)
        plt.savefig(f)


        sols = []
        ms3 = [(0.1, 0)] + ms3

        for m in ms3:
            out.append(f'm3 = {m}, {mf(m)}')
            data = open('data/eps_base').readlines()
            data = [d.strip('\n') for d in data]
            data.append(f'phi {m[1]}')
            data.append(f'gamma {m[0]}')
            data.append(f'n_states {N**L}')
            data.append(f'conserve None')
            data.append(f'N {N}')
            data.append(f'L {L}')
            open('data/eps_test', 'w+').write('\n'.join(data))
            S = SystemSet("eps_test")
            S.run()
            # S.purge_graphs()
            # S.extras()
            s = S.systems[0]
            pfs, e0 = pfs_matrix(s)
            # sols.append(m, s)
            print(e0)
            psiL = s.eng.Vleft[:,s.eng.inds[0]]
            psiR = s.eng.V[:,s.eng.indsleft[0]]
            o2 = inner(psiL, psiR, do_conj=0)
            print(f'overlap = {o2}')
            out.append(f'overlaps = {o2}')
            print(norm(psiL), norm(psiR))

            for i in range(N**L):
                continue
                j = s.eng.inds[i]
                k = s.eng.indsleft[i]
                psiL = s.eng.Vleft[:,j]
                psiR = s.eng.V[:,k]
                o2 = inner(psiL, psiR, do_conj=0)
                o2 = abs(o2)
                # print(i, o2, s.eng.E[j], s.eng.Eleft[j])
                print(i, o2)
                # print(i, o2, s.eigs[j][0], s.eigs_left[i][0])

        print('\n'.join(out))
            # print(m, s)
            # print(e0)
            # for i in range(N**L):
            #     for j in range(N**L):
            #         psiL = s.eigs_left[i][1]
            #         psiR = s.eigs[j][1]
            #         o = inner(psiL, psiR, do_conj=1)
            #         o = abs(o)
            #         o2 = inner(psiL, psiR, do_conj=0)
            #         o2 = abs(o2)
            #         print(i, j, o, o2, s.eigs[j][0], s.eigs_left[i][0])

                    # for c in [0, 1]:
                    #     o = inner(psiL, psiR, do_conj=c)
                    #     print(i, j, c, o)

            # for i in range(L):
            #     biggest = (None, 0)
            #     for j in range(L):
            #         psiL = s.eigs_left[i][1]
            #         psiR = s.eigs[j][1]
            #         o = inner(psiL, psiR, do_conj=1)
            #         o = abs(o)
            #         o2 = inner(psiL, psiR, do_conj=0)
            #         o2 = abs(o2)
            #         if o2 > biggest[1]:
            #             biggest = (j, o2)
            #     # print(i, j, o, o2, s.eigs[j][0], s.eigs_left[i][0])
            #     print(i, *biggest)

