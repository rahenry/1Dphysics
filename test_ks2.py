#!/usr/bin/env python3


from System import *
from SystemSet import *
from args import Q
from exact_fp import *
from rah_utility import *
from output_figs import *
from fp_k import *
from matrix import *

import numpy, scipy

basename = "fin1"
Ls = [3, 4, 5]
# Ls = [100, 150]
Ns = [2, 3, 4, 5, 6]
q1 = 3.5
q0 = -0.1
zpow = 0.2
vmax = None
vmin = 0
Ls = [50]
Ns = [3]
Ns = [2, 3]
k0 = -0.1
k1 = 3.0
Ls = [3]
Ls = [3, 4]


# Ls = [3, 5, 20, 51, 100, 200, 300]
# q1 = 0
# k0 = k1 / 5
# k1 = 1.2 * k0

basename = 'graphs/ep_test'
silent_remove(basename)
mkdir(basename)
cbase = f"""
lambda 1
method exact matrix
quiet 1
phi 0 0.1
angle_style 0
fix_pfs 1
conserve None
"""
npoints = 1001
ks = numpy.linspace(k0, k1, npoints)
basedir = 'graphs/test_ks2'
silent_remove(basedir)
mkdir(basedir)

basename = 'graphs/ep_test'
silent_remove(basename)
mkdir(basename)
ind0 = 0
to_plot = [1,2,3]
datasets = []
for N in Ns:
    for L in Ls:
        # basename = os.path.join('graphs/test_ks2', f'{N}_{L}')
        res = get_eps(L, N, k0, k1, q0, q1, 501)
        # print(res['evals'])
        z = numpy.power(res['z'], 0.2)
        f = os.path.join(basename, f'{N}_{L}')
        plot1(z, f+'.png', extent=res['ex'])
        # res["minima"] = []
        if L > 10:
            for i, (lamb, phi) in enumerate(res["minima"]):
                ep_test_large(L, N, lamb, phi)

        else:
            # res["minima"].append((1, 0))
            # res["minima"].append((0.1, 0))
            # print(res['evals'])
            for i, (lamb, phi) in enumerate(res["minima"]):
                print(lamb, phi)
                res = ep_test2(L, N, lamb, phi)
                datasets.append(res)
                f = os.path.join(basename, f'{N}_{L}_{lamb}_{phi}_{i}_{ind0}')
                plot_eig_info(res, f)
                ind0 += 1
                # delta, overlap, overlap2 = ep_test(L, N, lamb, phi)
                # name = basename + f'_{lamb}_{phi}_delta.png'
                # plot1(delta, name)
                # name = basename + f'_{lamb}_{phi}_overlap.png'
                # plot1(overlap, name)
                # name = basename + f'_{lamb}_{phi}_overlap2.png'
                # plot1(overlap2, name)

        # get_ep_lambdas_grid(L, N, res)

f = os.path.join(basename, f'{N}_{L}')
plot_eig_sets(datasets, f)
f = os.path.join(basename, 'plotA')
plot_eig_sets([datasets[q] for q in to_plot], f)

exit()

c = """
lambda 1
lambda 0.76
L 5
N 3
method exact matrix
quiet 1
phi 0 0.1 0.769
angle_style 0
fix_pfs 1
calc_spectrum 1
conserve None
"""
S = SystemSet.from_str(c)
S.run()
for s in S.systems:
    print(s, s["gamma"], s["L"])
    print(s["e_0"])
    # pfs, e0 = get_pfs_lamb(s['N'], s['L'], s['lambda'], s['phi'])

    # print(e0 / s['lambda_scale'] / s['L'])
    f = os.path.join(S.graph_dir, f"{s}.png")
    xdata = [e.real for e in s["energies"]]
    ydata = [e.imag for e in s["energies"]]
    plt.clf()
    plt.plot(xdata, ydata, linestyle="", marker="o")
    plt.savefig(f)
    print(f)
exit()

c = c + "\n gamma 0.5"
S = SystemSet.from_str(c)
S.run()
for s in S.systems:
    print(s)
    print(s["e_0"])

d = f"graphs/test_ks2_{basename}"
silent_remove(d)
mkdir(d)
for L in Ls:
    print("_" * 20)
    print(f"L = {L}")
    qs = numpy.linspace(q0 / L, q1 / L, npoints)
    # X = numpy.column_stack((ks, qs))
    X = cartesian_product(ks, qs)
    ind = 0
    ybase = numpy.array([a[0] + a[1] * 1.0j for a in X])
    # for func in [h, h2]:
    for func in [h]:
        ind += 1
        f = f"{L}_test_{func.__name__}.png"
        f = os.path.join(d, f)
        y = func(ybase, L)
        y = abs(y)
        z, ex = rastify(X, y)
        zplot = numpy.power(z, zpow)
        plt.clf()
        aspect_ratio = 1
        aspect = abs(ex[1] - ex[0]) / abs(ex[3] - ex[2]) * aspect_ratio
        # P = plt.imshow(z, vmin=0, vmax=1, extent=ex, aspect=aspect, origin='lower')
        P = plt.imshow(
            zplot, extent=ex, aspect=aspect, origin="lower", vmax=vmax, vmin=vmin
        )
        plt.colorbar(P)

        ms = find_minima(z)
        ms_coord = [(ks[m[1]], qs[m[0]]) for m in ms]
        ms_coord = sorted(ms_coord, key=lambda x: x[0])
        print(len(ms))

        def mf(x):
            return abs(func(x[0] + 1.0j * x[1], L))

        params = []
        aprev = None
        for a, b in ms_coord:
            plt.scatter(a, b, marker="+", c="red")
            bounds = make_bounds(a, b, ks, qs, Ls, ms_coord)

            meth = "Nelder-Mead"
            q = scipy.optimize.minimize(
                mf, (a, b), method=meth, bounds=bounds, tol=1e-26
            )
            if aprev:
                dx = q.x[0] - aprev
                # print(f'dx = {dx}')
            aprev = q.x[0]
            print(q)
            print((a, b), q.x, abs(mf((a, b))), abs(mf(q.x)))
            # print(abs(mf((a,b))), abs(mf(q.x)))
            plt.plot(q.x[0], q.x[1], marker="o", c="white", markersize=2)

            k = q.x[0] + 1.0j * q.x[1]
            lamb = k_to_lambda(k, L, N)
            phi = numpy.angle(lamb) * N / 2.0 / numpy.pi
            lamb = abs(lamb)

            params.append((lamb, phi))
        plt.savefig(f)

        for lamb, phi in params:
            print("@#*$&(@")
            print(lamb, phi)
            pfs, e0 = get_pfs_lamb(N, L, lamb, phi)
            spect = make_pf_spectrum2(e0, pfs, N)
            print(pfs[0], pfs[1])
            print(".")
            # for p in pfs:
            #     print(p)
            # print('.')
            if L < 20:
                count = 0
                for i in range(len(spect)):
                    for j in range(len(spect)):
                        if i >= j:
                            continue
                        e1 = spect[i]
                        e2 = spect[j]
                        if rel_diff(e1, e2) < 1e-6:
                            count += 1
                print(f"count = {count}")
            print(".")

            nprint = 3
            for x in spect[:nprint]:
                print(x)

        if L < 6:
            for lamb, phi in params:
                ep_test(N, L, lamb, phi)

        # xdata = ks
        # ydata = func(xdata, L)
        # plt.clf()
        # plt.scatter(xdata, ydata)
        # f = f'{L}_scatter_{func.__name__}.png'
        # f = os.path.join(d, f)
        # plt.savefig(f)

