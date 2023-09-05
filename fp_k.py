#!/usr/bin/env python3


import numpy, scipy
from exact_fp import make_m1, fix_pfs
from rah_utility import rastify, cartesian_product, find_minima, mkdir, silent_remove
from tenpy.linalg.np_conserved import tensordot, inner, norm, outer
import matplotlib.pyplot as plt
import labels
import os
from matrix import analyse_eigs, plot_eig_info

def k_inf(j, L, N, lamb):
    return numpy.pi * j / L - numpy.pi * j / L / L / (1. + lamb**(-0.5*N))

def eps_k(k, lamb, N):
    if lamb < 0:
        lamb = complex(lamb)
    res = 1. + lamb**(N) + 2.*(lamb**(N/2.))*numpy.cos(k)
    return res ** (1./N)

def eps_to_k(e, g, N):
    if g < 0:
        g = complex(g)
    res = (e ** N) - 1. - (g ** -N)
    res /= 2. * (g ** (-N/2.))
    # res = -res
    # res = abs(res)
    r1 = res
    res = numpy.arccos(res)
    return res

def get_pfs_lamb(N, L, lamb, phi):
    d1 = dict(L=L, N=N, phi=phi)
    d1['lambda'] = lamb
    m1 = make_m1(d1)
    sol = scipy.linalg.eig(m1)
    g = numpy.exp(phi * 2.j * numpy.pi / N) * lamb
    eigs = [e ** (1./N) for e in sol[0]]
    eigs = [e*g for e in eigs]
    eigs = fix_pfs(eigs, N)
    eigs.sort(key=lambda x: x.real)
    e0 = 0
    for e in eigs:
        e0 -= e
    return eigs, e0

def f_kj(k, lamb, L, N):
    return numpy.sin((L+1.) * k) + (lamb ** (-N/2.)) * numpy.sin(L*k)

def f_kj1(k, lamb, L, N):
    return (L+1) * numpy.cos((L+1.) * k) + (lamb ** (-N/2.)) * L * numpy.cos(L*k)

def f_kj2(k, lamb, L, N):
    return - (L+1.) * (L+1.) * numpy.sin((L+1.) * k) - (lamb ** (-N/2.)) * L*L *  numpy.sin(L*k)

def h(k, L, D=1):
    # res = (L+1) * numpy.cos((L+1)*k) * numpy.sin(L*k) - L * numpy.sin((L+1)*k)*numpy.cos(L*k)
    return (L+D) * numpy.cos((L+D)*k) * numpy.sin(L*k) - L * numpy.sin((L+D)*k)*numpy.cos(L*k)
def h2(k, L):
    return numpy.sin((2*L+1)* k) - (2*L+1) * numpy.sin(k)

def k_to_lambda(k, L, N, D=1):
    res = -(L+D)/L
    if abs(k) > 1E-16:
        res = -numpy.sin((L+D)*k) / numpy.sin(L*k)
    res = res ** (-2./N)
    q = 2*numpy.pi/N
    a = numpy.angle(res)
    v = q
    if a < 0: v = -q
    while(abs(a) + 1E-6 > q):
        a -= v
    phi = a * N / 2. / numpy.pi

    return abs(res), phi

def make_bounds(a, b, ks, qs, ms):
    xs = (None, None)
    if len(ms) > 1:
        l = int(len(ms) / 2)
        if l+1 >= len(ms): l = 0
        dx = 0.5 * abs(ms[l+0][0] - ms[l+1][0])
        xs = (a - dx, a + dx)
    return (xs, (qs[0], qs[-1]))



cbase = """
N 3
method exact
quiet 1
phi 0 0.1
angle_style 0
fix_pfs 1
conserve None
"""
def ep_test2(L, N, lamb, phi):
    # phi = abs(phi)
    c = cbase + "\n".join(
        [
            f"phi {phi}",
            f"lambda {lamb}",
            f"L {L}",
            f"N {N}",
        ]
    )
    from SystemSet import SystemSet
    S = SystemSet.from_str(c)
    S.run()
    s = S.systems[0]
    if s["method"] == "matrix":
        p0 = s["pfs"][0]
        p1 = s["pfs"][1]
        d = p1 - p0
    res = analyse_eigs(s.eng.E, s.eng.Eleft, s.eng.V, s.eng.Vleft)
    c = c.replace('exact', 'matrix')
    c += '\n calc_spectrum 1'
    S = SystemSet.from_str(c)
    S.run()
    s = S.systems[0]
    return res

def ep_test(L, N, lamb, phi):
    # phi = abs(phi)
    c = cbase + "\n".join(
        [
            f"phi {phi}",
            f"lambda {lamb}",
            f"L {L}",
            f"N {N}",
        ]
    )
    from SystemSet import SystemSet
    S = SystemSet.from_str(c)
    S.run()
    s = S.systems[0]
    if s["method"] == "matrix":
        p0 = s["pfs"][0]
        p1 = s["pfs"][1]
        d = p1 - p0

    nstates = len(s.eng.E)
    overlap = numpy.full((nstates, nstates), 0.)
    overlap2 = numpy.full((nstates, nstates), 0.)
    delta = numpy.full((nstates, nstates), 0.)
    for i in range(nstates):
        for j in range(nstates):
            I = s.eng.inds[i]
            J = s.eng.indsleft[j]
            if i > j:
                continue
            eR = s.eng.E[I]
            eL = s.eng.Eleft[J]
            vR = s.eng.V[:,I]
            vL = s.eng.Vleft[:,J]

            d = abs(eR - eL)
            delta[i][j] = d
            delta[j][i] = d

            o = inner(vL, vR, do_conj=0)
            o = abs(o)
            overlap[i][j] = o
            overlap[j][i] = o

            o = inner(s.eng.V[:, s.eng.inds[j]], vR, do_conj=1)
            o = abs(o)
            overlap2[i][j] = o
            overlap2[j][i] = o
    return delta, overlap, overlap2

def ep_test_large(L, N, lamb, phi):
    phi = abs(phi)
    c = cbase + "\n".join(
        [
            f"phi {phi}",
            f"lambda {lamb}",
            f"L {L}",
            f"N {N}",
            f"method matrix",
        ]
    )
    from SystemSet import SystemSet
    S = SystemSet.from_str(c)
    S.run()
    s = S.systems[0]
    max_ang = 0
    for i in range(4):
        p = s['pfs'][i] * L
        ang = numpy.angle(p) / numpy.pi
        if ang > max_ang:
            max_ang = ang

    if max_ang > 0:
        print(N, L, max_ang)
        exit()

def plot1(z, name, **kwargs):
    plt.clf()
    if 'extent' in kwargs:
        ex = kwargs['extent']
        aspect_ratio = 1
        aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio
        kwargs['aspect'] = aspect
    # P = plt.imshow(z, vmin=0, vmax=1, extent=ex, aspect=aspect, origin='lower')
    P = plt.imshow(
        z, origin="lower", **kwargs
    )
    plt.colorbar(P)
    plt.savefig(name)

def get_eps(L, N, k0, k1, q0, q1, npoints=501, D=1):
    ks = numpy.linspace(k0, k1, npoints)
    qs = numpy.linspace(q0/L, q1/L, npoints)

    X = cartesian_product(ks, qs)
    ybase = numpy.array([a[0] + a[1] * 1.j for a in X])
    func = h

    y = func(ybase, L, D)
    y = abs(y)
    z, ex = rastify(X, y)
    res = {}
    res['z'] = z
    res['ex'] = ex

    ms = find_minima(z)
    ms_coord = [(ks[m[1]], qs[m[0]]) for m in ms]
    ms_coord = sorted(ms_coord, key=lambda x: x[0])
    res['minima_bins'] = ms_coord

    def mf(x):
        return abs(func(x[0] + 1.j * x[1], L, D))
    params = []
    aprev = None
    res['minima'] = []
    res['minima_k'] = []
    res['evals'] = []
    for a, b in ms_coord:
        bounds = make_bounds(a, b, ks, qs, ms_coord)
        meth = "Nelder-Mead"
        q = scipy.optimize.minimize(mf, (a, b), method=meth, bounds=bounds, tol=1E-26)
        k = q.x[0] + 1.j * q.x[1]
        lamb, phi = k_to_lambda(k, L, N, D)
        # phi = numpy.angle(lamb) * N / 2. / numpy.pi
        # lamb = abs(lamb)
        res['minima'].append((lamb, phi))
        res['minima_k'].append(k)
        res['evals'].append(mf(q.x))
    return res

def get_ep_lambdas_grid(L, N, data, npoints=50, n_eps=-1, eps_ind_offset=0, lims=[[None, None], [None, None]], style='square'):
    print(f"Making lambda grid, style = {style}")
    minima = data["minima"]

    minima_filtered = []
    for x, y in minima:
        if y < 0:
            continue
        minima_filtered.append((x, y))
    minima = minima_filtered
    minima.sort(key = lambda x: x[0])
    minima = minima[eps_ind_offset:n_eps+eps_ind_offset]

    xlims = [None, None]
    ylims = [None, None]
    # for y, x in minima:
    #     if style == 'circle':
    #         1
    #     if not xlims[0] or x < xlims[0]:
    #         xlims[0] = x
    #     if not xlims[1] or x > xlims[1]:
    #         xlims[1] = x
    #     if not ylims[0] or y < ylims[0]:
    #         ylims[0] = y
    #     if not ylims[1] or y > ylims[1]:
    #         ylims[1] = y
    # dx = xlims[1] - xlims[0]
    # dy = xlims[1] - xlims[0]
    # dx = 0.1 * dx
    # dy = 0.1 * dy
    # if len(minima) == 1:
    #     k = 0.02
    #     dx = k * abs(xlims[0]) + 0.001
    #     dy = k * abs(ylims[0]) + 0.001
    # xlims = [xlims[0] - dx, xlims[1] + dx]
    # ylims = [ylims[0] - dy, ylims[1] + dy]

    if lims[0][0]: xlims[0] = lims[0][0]
    if lims[0][1]: xlims[1] = lims[0][1]
    if lims[1][0]: ylims[0] = lims[1][0]
    if lims[1][1]: ylims[1] = lims[1][1]

    lambstr = f"lambda {ylims[0]} {ylims[1]} {npoints}"
    phistr = f"phi {xlims[0]} {xlims[1]} {npoints}"
    if style == 'circle':
        lambstr = f"lambda_real {xlims[0]} {xlims[1]} {npoints}"
        phistr = f"lambda_imag {ylims[0]} {ylims[1]} {npoints}"


    z = numpy.full((npoints, npoints), 0.)
    c = "\n".join(
        [
            phistr,
            lambstr,
            f"L {L}",
            f"N {N}",
            f"method matrix",
            # f"method matrix_limited",
            "quiet' 0",
            "skip_model 1",
            "fix_pfs 1",
        ]
    )
    from SystemSet import SystemSet
    S = SystemSet.from_str(c)
    S.run()
    X = []
    y = []
    ytest1 = []
    for s in S.systems:
        d0 = 10000
        yt1 = 0
        pfs = sorted(s['pfs'], key=lambda x: abs(x))
        for j in range(L-1):
            # pfs = s['pfs']
            p0 = pfs[0]
            p1 = pfs[j+1]
            for i in range(N):
                p1a = p1 * numpy.exp(2.j * numpy.pi * i / N)
                d = abs(p0 - p1a)
                if d < d0:
                    d0 = d
                    yt1 = float(i)
        # for j in range(L):
        #     for k in range(L):
        #         if j >= k: continue
        #         p0 = pfs[j]
        #         p1 = pfs[k]
        #         for i in range(N):
        #             p1a = p1 * numpy.exp(2.j * numpy.pi * i / N)
        #             d = abs(p0 - p1a)
        #             if d < d0:
        #                 d0 = d
        #                 yt1 = float(i)
        if style == 'circle':
            X.append((s['lambda_real'], s['lambda_imag']))
        else:
            X.append((s['phi'], s['lambda']))
        # X.append((s['phi'], s['lambda']))
        y.append(d0)
        ytest1.append(yt1)


    z, ex = rastify(X, y)
    res = {}
    res['z'] = z
    res['ex'] = ex
    return res

def solve_find_eps(s):
    lab = s.identifier
    L = s['L']
    N = s['N']
    q1 = 3.5
    q0 = -0.1
    k0 = -3
    k1 = -k0
    res = get_eps(L, N, k0, k1, q0, q1, 501)
    f = os.path.join(s.parent.graph_dir, f'{lab}_ks.png')
    z = res['z']
    plot1(z, f, extent=res['ex'])
    for i in range(len(res['minima'])):
        k = res['minima_k'][i]
        lamb, phi = res['minima'][i]
        e = res['evals'][i]
        # print(k, lamb, phi, e)


    # lims = [[s['xlim0'], s['xlim1']], [s['ylim0'], s['ylim1']]]
    # res2 = get_ep_lambdas_grid(L, N, res, npoints=s['lambda_grid_npoints'], n_eps=10, eps_ind_offset=5, lims=lims)
    # f = os.path.join(s.parent.graph_dir, f'{lab}_eps.png')
    # ex = res2['ex']
    # z = res2['z']
    # aspect_ratio = 1
    # aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio
    # plt.clf()
    # p = plt.imshow(z, extent=ex, origin='lower', aspect=aspect)
    # plt.colorbar(p)
    # plt.xlabel(r'$\phi$')
    # plt.ylabel(r'$\lambda$')
    # plt.savefig(f)

    lims = [[-s['clim0'], s['clim0']], [-s['clim0'], s['clim0']]]
    res2 = get_ep_lambdas_grid(L, N, res, npoints=s['lambda_grid_npoints'], n_eps=10, eps_ind_offset=5, lims=lims, style='circle')
    f = os.path.join(s.parent.graph_dir, f'{lab}_eps_circ.png')
    ex = res2['ex']
    z = res2['z']
    aspect_ratio = 1
    aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio
    plt.clf()
    p = plt.imshow(z, extent=ex, origin='lower', aspect=aspect)
    plt.colorbar(p)
    plt.xlabel(r'$\Re(\lambda)$')
    plt.ylabel(r'$\Im(\lambda)$')

    for i in range(len(res['minima'])):
        k = res['minima_k'][i]
        lamb, phi = res['minima'][i]
        e = res['evals'][i]
        a = lamb * numpy.exp(2.j * numpy.pi / s['N'] * phi)
        plt.plot([a.real], [a.imag], marker='x', color='red')
    plt.savefig(f)


def eps_to_k_xy(eps, gamma):
    eps *= 2
    q = 1. - eps**2
    q /= (1. - gamma**2)
    q = numpy.arcsin(numpy.sqrt(q))
    return q

def eps_to_k_xy2(eps, gamma):
    eps *= 2
    eps /= 4. / (1. + gamma)
    q = 1. - eps**2
    q /= (1. - gamma**2)
    q = numpy.arcsin(numpy.sqrt(q))
    return q
