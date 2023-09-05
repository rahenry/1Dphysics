#!/usr/bin/env python3

import math
import numpy
import scipy
from exact_fp import *
from decimal import *
from scipy.special import hyp2f1
from scipy.optimize import root
from scipy.optimize import brentq
from scipy.optimize import fsolve
from rah_utility import mkdir, silent_remove
from rah_numpy import eigenvalue_test
import os
mkdir('test_graphs')
numpy.set_printoptions(threshold=numpy.inf)
# import matplotlib.pyplot as plt
# matplotlib.use("pgf")
# plt.rcParams.update(
#     {
#         "pgf.texsystem": "pdflatex",
#         "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}\usepackage{amsfonts}",
#         # "text.latex.preamble": r'\usepackage{amsfonts}\usepa',
#         # "font.family": "serif",
#         "text.usetex": True,
#         # "pgf.rcfonts": False,
#     }
# )
# xdata = [1,2,3]
# ydata = [2,3,5]
# # plt.plot(xdata, ydata)
# # plt.xlabel('hi $\mu \leqq$')
# fig, axs = plt.subplots(2,2)
# axs[0][0].plot(xdata, ydata)
# a = axs[0][0]
# a.set_xlabel(r"$\leq\leqq$")
# plt.savefig("junk.pdf")
# plt.savefig("junk.pgf")
# plt.close()

import matplotlib as mpl
mpl.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "pgf.preamble": "\n".join([
         r"\usepackage{url}",            # load additional packages
         r"\usepackage{unicode-math}",   # unicode math setup
         r"\setmainfont{DejaVu Serif}",  # serif font via preamble
    ])
})

fig, ax = plt.subplots(figsize=(4.5, 2.5))

ax.plot(range(5))

ax.set_xlabel("unicode text: я, ψ, €, ü")
ax.set_ylabel(r"\url{https://matplotlib.org}")
ax.legend(["unicode math: $λ=∑_i^∞ μ_i^2$"])

fig.tight_layout(pad=.5)

fig.savefig("pgf_preamble.pdf")
fig.savefig("pgf_preamble.png")
fig.savefig("pgf_preamble.pgf")
exit()

def f1():
    L0 = 100
    L1 = 700
    nL = 5
    N = 2
    g0 = 0.499
    g1 = 0.501
    ng = 30
    npfs = 50
    base = 'test_graphs/f1'
    silent_remove(base)
    mkdir(base)
    biggest = 0
    for L in numpy.linspace(L0, L1, nL):
        xdata = []
        ydata = []
        L = int(L)
        name = f'L={L}.png'
        name = os.path.join(base,name)
        for gamma in numpy.linspace(g0, g1, ng):
            e, x, pfs = exact_parafermions_matrix(L, N, gamma)
            maxpf = max(pfs)
            if maxpf > biggest:
                biggest = maxpf
            pfs = sorted(pfs)
            xdata.append(gamma)
            ydata.append(pfs[0:npfs])
            ydata[-1].append(pfs[-1])
        import matplotlib.pyplot as plt
        plt.clf()
        plt.plot(xdata, ydata)
        # plt.ylim([0, biggest])
        plt.savefig(name)



    exit()
    plt.clf()
    gamma = 0.5
    xdata = []
    ydata = []
    L0 = 100
    L1 = 1000
    dL = 100
    for L in numpy.arange(L0, L1, dL):
        print(L)
        e, x, pfs = exact_parafermions_matrix(L, N, gamma)
        pfs.sort()
        pfs = pfs[:npfs]
        xdata.append(L)
        ydata.append(pfs)
    plt.plot(xdata, ydata)
    plt.show()

def f2():
    L0 = 100
    L1 = 700
    nL = 5
    N = 3
    g0 = 0.3
    g1 = 1.0 - g0
    ng = 100
    base = 'test_graphs/f2'
    silent_remove(base)
    mkdir(base)
    biggest = 0
    for L in numpy.linspace(L0, L1, nL):
        print(L)
        xdata = []
        ydata = []
        L = int(L)
        name = f'L={L}.png'
        name = os.path.join(base,name)
        for gamma in numpy.linspace(g0, g1, ng):
            e, x, pfs = exact_parafermions_matrix(L, N, gamma)
            pfs = sorted(pfs)
            ratios = [pfs[i]/pfs[-1] for i in range(len(pfs)-1)]
            xdata.append(gamma)
            ydata.append(ratios)
        import matplotlib.pyplot as plt
        plt.clf()
        plt.plot(xdata, ydata)
        # plt.ylim([0, biggest])
        plt.savefig(name)
    1
# f1()
# f2()

def f3():
    L0 = 100
    L1 = 400
    nL = 5
    N = 3
    g0 = 0.43
    g1 = 0.45
    ng = 10
    base = 'test_graphs/f3'
    silent_remove(base)
    mkdir(base)
    biggest = 0
    for L in numpy.linspace(L0, L1, nL):
        print(L)
        xdata = []
        ydata = []
        L = int(L)
        name = f'L={L}.png'
        name = os.path.join(base,name)
        for gamma in numpy.linspace(g0, g1, ng):
            e, x, pfs = exact_parafermions_matrix(L, N, gamma)
            pfs = sorted(pfs)
            ratios = [pfs[i]/pfs[-1] for i in range(len(pfs)-1)]
            xdata.append(gamma)
            ydata.append(ratios[1:5])
        import matplotlib.pyplot as plt
        plt.clf()
        plt.plot(xdata, ydata)
        # plt.ylim([0, biggest])
        plt.savefig(name)


def make_matrix_data(L, N, gamma):
    g = gamma
    g = g ** (0.5*N)
    rows = [0]
    cols = [0]
    data = [1]

    for i in range(L-1):
        rows.append(i+1)
        cols.append(i)
        data.append(g)
        rows.append(i)
        cols.append(i+1)
        data.append(g)
        rows.append(i+1)
        cols.append(i+1)
        data.append(1+g**2)
    m = numpy.zeros((L, L), dtype=complex)
    for (i,j,d) in zip(rows, cols, data):
        m[i, j] = d
    return m

def gfun(i, gamma, phi):
    res = gamma
    #res = (gamma * (phi ** i)) ** 3
    print(res)
    return res ** 3

def make_matrix_data_alt(L, N, gamma, phi=1):
    g = gamma
    g = g ** N
    rows = [0]
    cols = [0]
    data = [1]

    L0 = 3*L
    m = numpy.zeros((L0, L0), dtype=complex)
    for i in range(3*L-1):
        m[i+1][i] = 1
    for i in range(3*L-2):
        if i % 3 == 2:
            pass
        elif i % 3 == 1:
            m[i][i+2] = gfun(i, gamma, phi)
            m[i][i+2] = g
        else:
            m[i][i+2] = 1
    return m

def make_matrix_data_alt2(L, N, gamma, phi=1):
    g = gamma
    g = g ** N
    L0 = N*L
    m = numpy.zeros((L0, L0), dtype=complex)
    m = numpy.zeros((L0, L0))
    b = 0
    count = 0
    while (True):
        if b+1 >= L0: break
        m[b+1][b] = 1
        b += 1
    s = 0
    while (True):
        x = s*N
        y = s*N+N-1
        if x >= L0 or y >= L0:
            break
        else:
            count += 1
            m[x][y] = 1
        s += 1
    s = 0
    while (True):
        x = s*N+1
        y = s*N+N
        if x >= L0 or y >= L0:
            break
        else:
            count += 1
            m[x][y] = g
        s += 1
    print(count, 2*L-1)
    return m

def mred(L, N, gamma, phi=1):
    m = make_matrix_data_alt2(L, N, gamma, phi=1)
    rows = []
    for i in range(N*L):
        row = m[i]
        if i % N == 1:
            rows.append(row)
    m = numpy.array(rows)
    m = numpy.transpose(m)
    cols = []
    for i in range(N*L):
        col = m[i]
        if i % N == 0:
            cols.append(col)
    m = numpy.array(cols)
    m = numpy.transpose(m)
    return m

def mred2(L, e, gamma, phi=1):
    m = make_matrix_data_alt2(L, 2, gamma, phi=1)
    m2 = numpy.matmul(m, m)
    m2 = m
    B = numpy.zeros((L, L))
    i = 0
    cols = []
    d = 0
    while (True):
        I = 2*i + d
        if I >= len(m2):
            break
        j = 0
        row = []
        while (True):
            J = 2*j + d
            if J >= len(m2):
                break
            row.append(m2[I, J])
            j += 1
        cols.append(row)
        i += 1
    B = numpy.array(cols)
    return B
def make_matrix_data_ff(L, N, gamma, phi=1):
    L0 = 2*L
    m = numpy.zeros((L0, L0), dtype=complex)
    m = numpy.zeros((L0, L0))
    b = 0
    count = 0
    while (True):
        if b+1 >= L0: break
        z = 1
        if b % 2 == 1:
            z = gamma**2
        m[b+1][b] = 1
        m[b][b+1] = z
        b += 1
    return m

def f3():
    L = 40
    N = 3
    gamma = 10.j + 10
    gamma = -1.j
    a = 3
    c = numpy.exp(2.j*numpy.pi/a)
    b = 1
    gamma = b*c
    gamma = 0.1*1.j
    gamma = -1000
    gamma = 1.1 + 0.1 * 5.j
    gamma = 0.30
    gamma = 0.01
    # gamma = 1.j
    phi = numpy.exp(2.j*numpy.pi/2)
    phi = 1
    M1 = make_matrix_data(L, N, gamma)
    # M2 = make_matrix_data_alt(L, N, gamma, phi)
    M3 = make_matrix_data_alt2(L, N, gamma, phi)
    cubed = numpy.matmul(M3, M3)
    cubed = numpy.matmul(M3, cubed)
    M4 = make_matrix_data_ff(L, N, gamma, phi)
    which = 'SM'
    m1 = scipy.sparse.coo_matrix(M1)
    m2 = scipy.sparse.coo_matrix(M3)
    mff = scipy.sparse.coo_matrix(M4)
    mr = mred2(L, N, gamma, phi)
    mr = scipy.sparse.coo_matrix(mr)
    # print(M3)
    # print(M4)
    # print(M3-M4)
    # print(M2)
    # print('...')
    # print(M3)
    # print('...')
    # print(M3-M2)

    # sol1 = scipy.linalg.eig(M1)
    # sol2 = scipy.linalg.eig(M3)
    # sol0 = scipy.linalg.eig(M3, left=False, right=True)
    # sol0 = numpy.linalg.eig(M3)
    # m0 = M3 @ M3 @ M3
    # sol0 = scipy.sparse.linalg.eigs(m0, k=10, which='SM', maxiter=L*1000)
    # for i in range(len(sol0[0])):
    #     e = sol0[0][i]
    #     v = sol0[1][:,i]
    #     r, n, n2 = eigenvalue_test(m0, v)
    #     n = abs(n)
    #     n2 = abs(n2)
    #     n0 = numpy.inner(numpy.conjugate(v), v)
    #     # print(n0, e)
    #     # print(abs(e/r), n, n2)
    #     print(e)
    # exit()
    q = int(L / 3.)
    sol1 = scipy.sparse.linalg.eigs(m1, k=q, which=which)
    # sol2 = scipy.sparse.linalg.eigs(mr, k=q, which=which)
    sol2 = scipy.sparse.linalg.eigs(cubed, k=3*q, which=which)
    solff = scipy.sparse.linalg.eigs(mff, k=2*q, which=which)
    e1s = sorted(sol1[0])
    things = []
    for e in sol2[0]:
        e = e ** (1./N)
        a = numpy.angle(e, deg=1)
        things.append((a, abs(e)))
        z = f'{a:>5.2f} {e:4.5f}'

    e1s.sort()
    e1s.reverse()
    things.sort(key=lambda x:x[1])
    things.reverse()
    effs = solff[0]
    effs.sort()
    effs = [x.real for x in effs if x > 0]
    effs.reverse()

    for i in range(len(e1s)):
        e1 = e1s[i].real
        e1 = e1 ** (1./N)
        e2 = things[i*N][1]
        # e2 = things[i][1]
        delta = e2-e1
        # arg = things[i][0]
        arg = things[i*N][0]
        eff = effs[i]
        # eff = eff ** (N/2.)
        # eff = eff ** (2./N)
        print(f'{e1:10f} {e2:10f} {delta:10f} {arg:10f} {eff:10f}')
        # print(f'{e1:10f}  {eff:10f}')

    # for t in things:
        # print(t)
        # print(z)
        # if abs(e.imag) < 1E-6:
        #     e = e.real
        #     print(e)
            # print(e**(1./N))

f3()
exit()

from scipy.special import hyp2f1
from scipy.optimize import root
from scipy.optimize import brentq
from scipy.optimize import fsolve

x = exact_parafermions_matrix(20, 3, 1.0, flipZ=True)
x = exact_parafermions_matrix(50, 3, 0.4, flipZ=1)
print(x[0])
pfs = x[2]
pfs.sort()
print(pfs)


#trying to do PT and find left eigenstate

        # def t2(psi):
        #     psi = psi.copy()
        #     for i in range(psi.L):
        #         site = psi.sites[i]
        #         print(site)
        #         sites = psi.sites
        #         print(sites)
        #         q = site.Z
        #         v = site.plocal
        #         o = site.get_op('Z')
        #         o = site.get_op('plocal')

        #         1
        # def transmogrify(psi):
        #     psi = psi.copy()
        #     psi = psi.spatial_inversion()
        #     Bs = []
        #     for q in range(psi.L):
        #         b = psi.get_B(q)
        #         a = b.to_ndarray()[:,[0,2,1],:]
        #         legs = b.legs
        #         legs[0].qconj *= -1
        #         legs[2].qconj *= -1
        #         bnew = npc.Array.from_ndarray(a, legs, labels=b._labels, qtotal=b.qtotal)
        #         Bs.append(bnew)
        #     snew = MPS(s.model.lat.mps_sites(), Bs, psi._S)
        #     return snew
        # s = y[1][0]
        # s.load_states()
        # s1 = s.states[0].psi
        # s2 = s.left_system.states[0].psi
        # s3 = t2(s2)
        # exit()
        # s3 = transmogrify(s2)
        # o1 = s1.overlap(s1)
        # o12 = s1.overlap(s2)
        # o3 = s3.overlap(s3)
        # o13 = s1.overlap(s3)
        # # print(o1)
        # # print(o12)
        # # print(o3)
        # # print(o13)
        # s3 = s2.spatial_inversion()
        # # s3 = transmogrify(s2, sref)
        # q = 2
        # b0 = s3.get_B(q);
        # bref = s1.get_B(q);
        # b = b0.to_ndarray()
        # # b = b[:,[2,1,0],:]
        # # b = b[:,:,[0,2,1]]
        # print(b0.shape)
        # print(b.shape)
        # print(bref.shape)
        # for i in range(b0.shape[0]):
        #     print('========================')
        #     for j in range(b0.shape[1]):
        #         print('---')
        #         for k in range(b0.shape[2]):
        #             x1 = b[i][j][k]
        #             x2 = bref[i][j][k]
        #             print(f'{x1:50} {x2:50}')

        # exit()
        # # b = npc.Array.from_ndarray(b, bref.legs, labels=bref._labels, qtotal=bref.qtotal)
        # c = npc.Array.from_ndarray(b, b0.legs, labels=b0._labels, qtotal=b0.qtotal)
        # print(c)
        # exit()
        # # shitty sprectrum test
        # # s = y[1][0]
        # # e, x, pfs = exact_parafermions_matrix(s['L'], s['N'], s['gamma'], s['flipZ'])
        # # spectrum = make_pf_spectrum(e, pfs, s['N'])
        # # for x in spectrum[:10]:
        # #     print(x)
