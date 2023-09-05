#!/usr/bin/env python3

import matplotlib.pyplot as plt
import os
from Slicer import SystemSlicer
from pf_tester import pf_test
from tenpy.linalg.np_conserved import tensordot, inner, norm
from rah_utility import mkdir, matplotlib_set_aspect
import graph_settings
import matplotlib
import numpy
from figures import *
from args import Q
from data_set import DataSet, plot_multiple
import cProfile, pstats, io
from pstats import SortKey
from exact_fp import exact_parafermions_matrix, make_pf_spectrum2, rhoX
from output_complexlambda import *
from output_extra import *
from output_eptest import *
from output_plots import *
from output_periodic import *
from output_ed2 import *


def standard_output(S):
    xkey = S["xkey"]
    if xkey == []:
        xkey = "gamma"
    for s in S.systems:
        print(s)
        print(s['fsL'])
        # print(s['e_0'], s['e_exact'])
        # print(s.auxillary_systems[0]['e_0'])
        saux = s.auxillary_systems[0]
        s.load_states()
        saux.load_states()
        for E in saux['solve_energies']:
            print(E/s['L'])
        for E in s['solve_energies']:
            print(E/s['L'])
        e1 = s['e_0']
        e2 = saux['e_0']
        e3 = saux['solve_energies'][1] / s['L']
        # print(e1, e2, e3)
        p1 = s.states[0].psi
        p2 = saux.states[0].psi
        def getE(p):
            return s.model.H_MPO.expectation_value(p)

        print(getE(p1), getE(p2))

        f = p1.overlap(p2) * p2.overlap(p1)
        d = s['fidelity_delta']
        fsl = (1.-f)/d/d
        print(f)
        print(fsl)
    for s in S.systems:
        print(s)
        print(s['fsL'])
        # print

    exit()
    for x, y in S.slicer.get_sets(exclude=[xkey]).items():
        print(x)
        slicer = SystemSlicer(y)
        # z = slicer.tabulate(["e_0", "e_exact", "o", "time", "sweep"])
        # print(z)
        # z = slicer.tabulate(
        #     ["e_0", "e_exact", "norm_0", "o2", "o3", "z", "z2", "x", "zc", "z2c", "xc"]
        # )
        # z = slicer.tabulate(["o", "o2", "o3", "z", "z2", "x", "zc", "z2c", "xc"])

        # z = slicer.tabulate(["o", "oL", "X1L", "x", "Z1L", "z", "smallest_pf_rel"])
        # print(z)
        # z = slicer.tabulate([xkey, "fidelity", "fsusceptibility"])
        # print(z)
        # for i in range(1, y[1][0]['n_states']):
        # z = slicer.tabulate([f'e_{i}_real', f'e_{i}_imag', f'e{i}_exact_real', f'e{i}_exact_imag'])
        #     print(z)
        #
        # z = slicer.tabulate([xkey, "e_0", "e_exact", "o"])
        # print(z)
        # z = slicer.tabulate(["fidelity", "fsL"])
        # print(z)
        z = slicer.tabulate([xkey, "e_0_real", "e_0_imag", "e_exact"])
        print(z)
        z = slicer.tabulate([xkey, "gs1", "gs1_new", "gs2", "gs3"])
        print(z)
        z = slicer.tabulate(['lambda', 'h', "fidelity", "fsL"])
        print(z)
        # z = slicer.tabulate(['o'])
        # print(z)
        # z = slicer.tabulate(["smallest_pf_rel"])
        # print(z)
        # z = slicer.tabulate(["x", "z"])
        # print(z)
        # z = slicer.tabulate(["sweep", "time"])
        # print(z)
        # # z = slicer.tabulate(["rhoX", "rhoX_exact"])
        # # print(z)
        # z = slicer.tabulate(["rhoX", "rhoX_exact", "rhoXtest", "rhoXtest2", "x1Ltest"])
        # print(z)
        # z = slicer.tabulate(["rhoZ", "rhoZ_exact", "rhoZtest", "rhoZtest2"])
        # print(z)
    return


def fit_z1L(S):
    d = {}
    xdata = []
    ydata = []
    for x, y in S.slicer.get_sets(exclude=S["xkey"]).items():
        xdata.append(y[0]["L"])
        m = 1000000
        gmin = None
        for s in y:
            z = s["z1L"]
            if z < m:
                gmin = s["gamma"]
                m = z
        ydata.append(gmin)

    d["xdata"] = xdata
    d["ydata"] = ydata
    print(d)
    dbase = {"all": d}

    d = DataSet(S, from_data=dbase)
    d1 = d
    p, pcov = d.make_ansatz("z1Lc", label="test")
    print(p)
    perr = numpy.sqrt(numpy.diag(pcov))
    print(perr)
    d.plot(linestyle="", xlab="L", ylab=r"$\gamma_{min}", use_legend=0, name="fit_L")
    gamma = p[-1]
    l = gamma / (1.0 - gamma)
    e, x, pfs = exact_parafermions_matrix(1000, 3, l, 0)
    print(pfs[0] / pfs[-1])

    d = DataSet(S, from_data=dbase)
    d2 = d
    p, pcov = d.make_ansatz("z1Le", label="test2")
    print(p)
    d.plot(linestyle="", xlab="L", ylab=r"$\gamma_{min}", use_legend=0, name="fit_Lpow")
    print(d)
    perr = numpy.sqrt(numpy.diag(pcov))
    print(perr)
    gamma = p[-1]
    l = gamma / (1.0 - gamma)
    e, x, pfs = exact_parafermions_matrix(1000, 3, l, 0)
    print(pfs[0] / pfs[-1])

    plt.clf()
    fig, axs = plt.subplots(1, 2)
    aspect = 1
    d1._plot(axs[0], xlab="L", ylab=r"$\gamma_\text{min}$", linestyle="")
    d2._plot(axs[1], xlab="L", ylab=r"$\gamma_\text{min}$", linestyle="")
    matplotlib_set_aspect(axs[0], aspect)
    matplotlib_set_aspect(axs[1], aspect)
    fig.tight_layout()
    f = os.path.join(S.graph_dir, "fits.png")
    plt.savefig(f)
    f = os.path.join(S.graph_dir, "pgf/fits.pgf")
    plt.savefig(f)
    # plt.show()


def standard_graphs(S):
    print("Generating standard graphs")
    xkey = S["xkey"]
    if not xkey:
        xkey = "gamma"

    def rescale1(y):
        if abs(y) > 50:
            y = -abs(y)
        return y

    def rescale_log(y):
        return numpy.log10(y)

    to_exec = []

    # DataSet(S, 'o', name='overlap_log', apply_func=rescale_log).plot()
    # DataSet(S, 'x', name='x_log', apply_func=rescale_log).plot()
    # DataSet(S, 'z', name='z_log', apply_func=rescale_log).plot()
    # DataSet(S, 'rhoZ', name='rhoZ').plot()
    # DataSet(S, 'rhoZtest').plot()
    # DataSet(S, 'rhoZtest').plot()
    # keyz = ['rhoXtest2', 'x1Ltest', 'z1Ltest', 'rhoX', 'xhalf', 'zhalf']
    keyz = S["graph_keys"]
    if keyz == []:
        keyz = ["z1L", "x1L", "rhoX", "xhalf", "zhalf"]
    for k in keyz:
        DataSet(S, k).plot()
    if os.path.exists(S.graph_info_file):

        # pr = cProfile.Profile()
        # pr.enable()
        for l in open(S.graph_info_file).readlines():
            if l[0] == "@":
                c = f"DataSet(S, {l[1:]}).plot()"
                print(c)
                exec(c)
            if l[0] == "$":
                l = l[1:]
                l = l.split()
                k1, k2 = l[0:2]
                args = " ".join(l[2:])
                c = f"plot_multiple([DataSet(S, '{k1}'), DataSet(S, '{k2}')], {args})"
                print(c)
                exec(c)
        # pr.disable()
        # s = io.StringIO()
        # sortby = SortKey.CUMULATIVE
        # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        # ps.print_stats()
        # print(s.getvalue())
        # exit()
    basic_keys = ["x", "z"]
    # for k in basic_keys:
    #     DataSet(S, k).plot()
    DataSet(S, "fsL", ylims=[-10, 10], name="fsL_small").plot()
    DataSet(S, "fsL", ylims=None, name="fsL_all").plot()
    DataSet(S, "e0_error").plot()
    DataSet(S, "o").plot()
    DataSet(S, "o", name="overlap_log", apply_func=rescale_log).plot()

    # DataSet(S, 'fsL', ylims = [1.0, 2.5], name='fsL_3').plot()


def exact_test1(S):
    for s in S.systems:
        print(s)
        xdata = [e.real for e in s.eng.E]
        ydata = [e.imag for e in s.eng.E]
        scale = abs(s.eng.E[0])
        xdata = [e.real / scale for e in s.eng.E]
        ydata = [e.imag / scale for e in s.eng.E]
        plt.clf()
        plt.plot(xdata, ydata, linestyle="", marker="o", markersize=0.5)
        plt.xlabel("Re(E)")
        plt.ylabel("Im(E)")
        g = s["gamma"]
        if s["flipZ"]:
            g = -g
        txt = f"L = {s['L']}, gamma = {g}, eta = {s['eta']}"
        plt.title(txt)
        f = f"{s}_spectrum.png"
        f = os.path.join(S.graph_dir, f)
        # plt.xlim([-0.5, -0.3])
        # plt.ylim([-0.1, 0.1])
        # plt.xlim([-0.42, -0.38])
        # plt.ylim([-0.01, 0.01])
        plt.tight_layout()
        plt.savefig(f, dpi=200)

        pf_test(s)

    1


def fp_test2(S):
    # for x, y in S.slicer.get_sets(exclude=['gamma']).items():
    labels = ["gamma", "L", "e0_ED_real", "e0_ED_imag", "e0_pfs_real", "e0_pfs_imag"]
    keys = ["gamma", "L", "e_0_real", "e_0_imag", "e_exact_real", "e_exact_imag"]
    widths = {"gamma": 4}
    z = S.slicer.tabulate(keys, labels, widths)
    print(z)
    print("")
    label = ["L", "gamma"]
    for s in S.systems:
        print(s.make_label(label))
        print(f"Ground state = {s['e_0_real']} + {s['e_0_imag']}i")
        print("Parafermion energies:")
        print(f'{"Real":>10} {"Imag":>10}')
        for x in s.pfs:
            e1 = x / (1.0 - s["gamma"])
            print(f"{e1.real:10.7f} {e1.imag:10.7f}")
        print("")


def fp_test1(S):
    # for x, y in S.slicer.get_sets(exclude=['gamma']).items():
    labels = [
        "lambda",
        "L",
        "e0_ED_lambda_real",
        "e0_ED_lambda_imag",
        "e0_pfs_lambda_real",
        "e0_pfs_lambda_imag",
    ]
    keys = [
        "lambda",
        "L",
        "e_0_lambda_real",
        "e_0_lambda_imag",
        "e_exact_lambda_real",
        "e_exact_lambda_imag",
    ]
    widths = {"lambda": 4}
    z = S.slicer.tabulate(keys, labels, widths)
    print(z)
    print("")
    label = ["L", "lambda"]
    for s in S.systems:
        print(s.make_label(label))
        print(f"Ground state = {s['e_0_lambda_real']} + {s['e_0_lambda_imag']}i")
        print("Parafermion energies:")
        print(f'{"Real":>10} {"Imag":>10}')
        for x in s.pfs:
            e1 = x / (1.0 - s["gamma"])
            print(f"{e1.real:10.7f} {e1.imag:10.7f}")
        print("")


def XYW_test(S):
    # if S['method'][0] == 'dmrg':
    if True:
        xkey = S["xkey"]
        for x, y in S.slicer.get_sets(exclude=[xkey]).items():
            print(x)
            slicer = SystemSlicer(y)
            # z = slicer.tabulate([xkey, "e_0", "e_exact", "o", "time", "sweep"])
            # print(z)
            # z = slicer.tabulate([xkey, "fidelity", "fsusceptibility"])
            # print(z)
            z = slicer.tabulate([xkey, "e_0", "e_exact", "o"])
            print(z)
            z = slicer.tabulate(["fidelity", "fsL"])
            print(z)
            z = slicer.tabulate(["x", "z"])
            print(z)
        return
    else:
        for s1 in S.systems:
            print(s1)
            srt = sorted(s1.eng.E, key=lambda x: x.real)
            print(srt[0] / s1["L"], s1["e_exact"], s1.left_system.eigs[0][0] / s1["L"])
            print(s1["fidelity"], s1["fsusceptibility"])
            print(s1["o"])
            print(s1["eta"], s1.fidelity_system["eta"])
            z = s1.fidelity_system
            q = inner(s1.eigs[0][1], z.eigs[0][1], do_conj=True)
            print(q)
        return
    for s1 in S.systems:
        for s2 in S.systems:
            print(s1, s2)
            for i in range(len(s1.eng.E)):
                best = 0
                ind = None
                for j in range(len(s1.eng.E)):
                    c = False
                    v1 = s1.eng.V[i]
                    v1 /= norm(v1)
                    v2 = s2.eng.V[j]
                    v2 /= norm(v2)
                    # v1 = v1.conj()
                    o = inner(v1, v2, do_conj=c)
                    o = abs(o)
                    if o > best:
                        ind = j
                        best = o
                print(i, ind, best)
            pass


def XYWp_test(S):
    for s in S.systems:
        print(s)
        print(s["e_exact"])
        # for e, v in s.eigs[:1]:
        #     print(e / s["L"])
    s1 = S.slicer.find({"bc": 1, "L": 7})[0]
    s2 = S.slicer.find({"bc": 0, "L": 8})[0]
    for i in range(len(s1["eigs"])):
        e1, v1 = s1["eigs"][i]
        e2, v2 = s1["eigs"][i]
        # print(e1 / s1["L"], e2 / s2["L"])
        print(e1, e2)


def basic_graphs(S):
    print("Making basic graphs")
    xkey = S["xkey"]
    if not xkey:
        xkey = "gamma"
    else:
        xkey = xkey[0]
    graph_keys = ["e_0", "e_exact", "e0_error"]
    if S["fidelity_parameter"]:
        graph_keys += ["fidelity", "fsusceptibility"]
    for x in graph_keys:
        print(x)
        graph1(S, xkey, x)


def show_spectrum_pf(S):
    for x, y in S.slicer.get_sets().items():
        for s in y:
            L = s["L"]
            xdata = [L * x.real for x in s.spectrum]
            ydata = [L * x.imag for x in s.spectrum]
            plt.plot(
                xdata,
                ydata,
                marker="o",
                markersize=1.1,
                linestyle="",
                color="blue",
            )
            plt.axis("equal")
            plt.xlabel("Re(E)")
            plt.ylabel("Im(E)")
            f = os.path.join(S.graph_dir, f"{x}.png")
            # plt.show()
            plt.savefig(f)
            plt.close()


def show_spectrum(S):
    for x, y in S.slicer.get_sets().items():
        for s in y:
            s.eigs = s["eigs"]
            L = s["L"]
            xdata_sim = [L * x[0].real for x in s.eigs]
            ydata_sim = [L * x[0].imag for x in s.eigs]
            plt.plot(
                xdata_sim,
                ydata_sim,
                marker="o",
                markersize=1.1,
                linestyle="",
                color="blue",
            )
            plt.axis("equal")
            plt.xlabel("Re(E)")
            plt.ylabel("Im(E)")
            f = os.path.join(S.graph_dir, f"{x}.png")
            # plt.show()
            title = rf"$L={s['L']},\;\gamma={s['gamma']},\;\eta={s['eta']},\;\kappa={s['kappa']}$"
            cZ = abs(s.model.cZ)
            cY = abs(s.model.cY)
            cX = abs(s.model.cX)
            cW = abs(s.model.cW)
            title += (
                rf"\\$c_Z={cZ:1.3f},\; c_X={cX:1.3f},\; c_Y={cY:1.3f},\; c_W={cW:1.3f}$"
            )
            plt.title(title)
            plt.savefig(f)
            plt.close()


def pf_tester(S):
    for x, y in S.slicer.get_sets().items():
        for s in y:
            L = s["L"]
            xdata_sim = [L * x[0].real for x in s.eigs]
            ydata_sim = [L * x[0].imag for x in s.eigs]
            xdata_exact = [L * x.real for x in s.spectrum]
            ydata_exact = [L * x.imag for x in s.spectrum]

            plt.clf()
            plot_args = dict(linestyle="", linewidth=10, markeredgewidth=1)
            plt.plot(
                xdata_exact,
                ydata_exact,
                marker="x",
                markersize=8,
                color="red",
                **plot_args,
            )
            plt.plot(
                xdata_sim,
                ydata_sim,
                marker="o",
                markersize=2,
                color="blue",
                **plot_args,
            )
            plt.xlabel("Re(E)")
            plt.ylabel("Im(E)")

            f = f"{x}_pftest.png"
            f = os.path.join(S.graph_dir, f)
            plt.savefig(f)
            # plt.show()
            plt.close()


def fidelity_full(S):
    xkey = S["fidelity_parameter"]
    for x, y in S.slicer.get_sets(exclude=[xkey]).items():
        xvals = y.get_possible_values(xkey)
        n = len(xvals)
        zdata = numpy.zeros((n, n))
        for i in range(n):
            s1 = y[i]
            X = s1[xkey]
            for j in range(n):
                s2 = y[j]
                Y = s2[xkey]
                psi1 = s1.states[0].psi
                psi2 = s2.states[0].psi
                f = psi1.overlap(psi2)
                zdata[i][j] = abs(f)
                psi1 = s1.left_system.states[0].psi
                f = psi1.overlap(psi2) / psi1.overlap(psi1) / psi2.overlap(psi2)
                # f = s1.states[0].psi.overlap(s2.states[0].psi)
        plt.imshow(zdata)
        plt.show()
        plt.clf()


def fidelity(S):
    xkey = S["fidelity_parameter"]

    print("Graphing fidelity")
    if Q.nprocs > 1:
        return
    for x, y in S.slicer.get_sets(exclude=[xkey]).items():
        f = f"fidelity_{x}.png"
        f = os.path.join(S.graph_dir, f)
        zdata = S.fidelity_data["nonh"][x]
        plt.imshow(zdata)
        plt.savefig(f)
        plt.clf()


def conv_output(S):
    xkey = S["xkey"]
    if xkey == []:
        xkey = "gamma"
    for x, y in S.slicer.get_sets(exclude=[xkey]).items():
        print(x)
        slicer = SystemSlicer(y)
        z = slicer.tabulate(["e_0", "e_exact", "o"])
        print(z)
        z = slicer.tabulate(["fidelity", "fsL"])
        print(z)
    return



from fp_more_matrices import *
from fp_polynomial import *


def fptest1(S):
    for s in S.systems:
        print("...")
        print(s["phi"], s["flipZ"])
        print(s["phi"])
        # pf_test(s)
        spect, pfs = make_spectrum(s)
        spect2, pfs2 = make_spectrum(s, mtype="m2")
        # print(pfs)
        # print(s['N'] ** s['L'], len(spect), len(s['energies']))
        pfs_poly = find_eps_poly(s)
        pfs_poly = [p / s["L"] for p in pfs_poly]
        e0 = 0
        for p in pfs_poly:
            e0 -= p
        spect_poly = make_pf_spectrum2(e0, pfs_poly, s["N"])
        for i in range(len(pfs_poly)):
            print(pfs_poly[i], pfs2[i])
        for i in range(len(spect)):
            b = spect2[i]
            b = s["energies"][i]
            # print(spect[i], spect2[i], s['energies'][i])
            # print(spect2[i], s['energies'][i])
            print(s["energies"][i], spect_poly[i])


def complex_plot(S):
    fig, axs = plt.subplots(1,2)
    print(axs)
    labels = []
    for x, y in S.slicer.get_sets(exclude="L").items():
        xdata1 = []
        ydata1 = []
        xdata2 = []
        ydata2 = []
        labels.append(x)
        for s in y:
            L = s["L"]
            xdata1.append(L)
            xdata2.append(L)
            ydata1.append(s["e_0_real"])
            ydata2.append(s["e_0_imag"])
        pltsty = dict(
            linestyle='-',
            marker='o')
        axs[0].plot(xdata1, ydata1, **pltsty)
        axs[1].plot(xdata2, ydata2, **pltsty)

        print(x)

    plt.legend(labels)
    plt.show()


def sep06_test(S):
    for s in S.systems:
        print(s['swap_lambda'], s['lambda'], s['phi'], s['lambda_full'])
        print(s['e_0'], s['e_inf2'])
        # q = s['e_inf2']
        # print(s['lambda'], s['phi'], s['e_0'], s['e_inf2'], abs(q))


def test5(S):
    for s in S.systems:
        print(s, s['e_exact'])
        for i in range(s['n_states']):
            print(s[f'e_{i}'])


    fp = S.systems[0]['fidelity_parameter']
    for x, y in S.slicer.get_sets(exclude=["fidelity_delta"]).items():
        print(x)
        for s in y:
            sf = s.auxillary_systems[0]
            s.load_states()
            for i in range(s['n_states']):
                # print(i)
                psi1 = s.states[i].psi
                fp = s['fidelity_parameter']
                # print(s[fp], sf[fp])
                psi2 = sf.states[i].psi

                o1 = psi1.overlap(psi2)
                o2 = psi2.overlap(psi1)
                f =  o1 * o2
                d = s['fidelity_delta']
                chi = (1.-f)/d/d

                E = s.model.H_MPO.expectation_value(psi1)
                E2 = sf.model.H_MPO.expectation_value(psi2)
                print(s['fidelity_delta'], f, chi)

def print_energies(S):
    print(S)
    for s in S.systems:
        print(s)
        es = sorted(s['solve_energies'])
        for e in es:
            print(e)

import free_test


def free_test1(S):
    for s in S.systems:
        print(s)
        # f = free_test.free_test(s, s["N"], 0)
        # print(f)
        print("_____")
        f2 = free_test.free_test_newnew(s["energies"], s["N"], 0)
        for z in f2:
            print(z)
        print("_____")
        print("_____")
        continue
        # for q in s.spectrum:
        #     print(q)
        u = 0
        for z in s.pfs:
            z = z * 2
            u += z
            print(z)
        print(u)
    print(2342)
    exit()
    1
