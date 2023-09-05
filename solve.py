#!/usr/bin/env python3

from tenpy.algorithms.dmrg import SingleSiteDMRGEngine, TwoSiteDMRGEngine
import numpy as np
from rah_utility import rel_diff
from tenpy.tools.params import Config
from tenpy.models.model import NearestNeighborModel
from tenpy.algorithms import tebd
from exact_diag import ExactDiagNonH
# from exact_fp import spectrum_matrix, pfs_matrix, fix_pfs, pfs_matrix_limited, make_m1, make_pf_spectrum2
from exact_fp import *
from NonHDMRG import NonHermitianDMRG, SSNonHermitianDMRG
from tenpy.linalg.np_conserved import tensordot, inner, norm
from State import State
from polytest import sympy_eps2
from point_generator import PointGenerator
from fp_k import *
from matrix import *
from simple_ed import *
from ed2 import solve_ed2
from ed3 import solve_ed3
import scipy
import logging

config_info = {
    "none":{
        },
    "sparse": {
        },
    "exact": {
        },
    "tebd": {
        "trunc_params" : ["chi_max", "svd_min"],
        },
    "dmrg": {
        "trunc_params": [
            "chi_max",
            "svd_min",
            "degeneracy_tol",
            "chi_min",
            "trunc_cut",
            "log_error_target",
        ],
        "mixer_params": ["decay", "amplitude", "disable_after"],
        "error_params": ["svd_min", "max_E_err", "max_S_err"],
        "lanczos_params": ["P_tol", "reortho", "N_max", "verbose", "N_min"],
    }
}
convert_to_bool = ['mixer']

def solve_matrix_limited(s):
    pfs = pfs_matrix_limited(s, 3)
    s['pfs'] = pfs
    s['pfs_real'] = [x.real for x in pfs]
    s['pfs_imag'] = [x.imag for x in pfs]
    s.other_data[f'converged_0'] = 1
    lowest = pfs[0]
    if s['gamma'] < 0.5:
        lowest = pfs[1]
    lowest *= s['L']
    lowest *= (1. - s['omega'])
    s['gap'] = lowest
    s['gap_real'] = lowest.real
    s['gap_imag'] = lowest.imag
    if s['calc_spectrum'] or s['method'] == 'exact':
        spect = make_pf_spectrum2(e0, pfs, s['N'])
        s['energies'] = spect
    if s['depth']:
        s['energies'] = make_pf_spectrum2(e0, pfs, s['N'], depth=s['depth'])

def solve_matrix(s):
    pfs, e0 = pfs_matrix(s)
    s['pfs'] = pfs
    s['pfs_real'] = [x.real for x in pfs]
    s['pfs_imag'] = [x.imag for x in pfs]
    s['e_0'] = e0 / s['L']
    s['E_0'] = e0
    s.other_data[f'converged_0'] = 1
    lowest = pfs[0]
    if s['gamma'] < 0.5:
        lowest = pfs[1]
    lowest *= s['L']
    lowest *= (1. - s['omega'])
    s['gap'] = lowest
    s['gap_real'] = lowest.real
    s['gap_imag'] = lowest.imag
    if s['depth']:
        s['energies'] = make_pf_spectrum2(e0, pfs, s['N'], depth=s['depth'])
    elif s['calc_spectrum'] or s['method'] == 'exact':
        spect = make_pf_spectrum2(e0, pfs, s['N'])
        s['energies'] = spect


def solve_exact(s, c, index):
    if s['model_type'] == 'fp':
        solve_matrix(s)
        s['energies_matrix'] = list(s['energies'])
        s['pfs_matrix'] = list(s['pfs'])
    elif s['model_type'] == 'xyh':
        solve_matrix_xyh(s)
        s['energies_matrix'] = list(s['energies'])
        s['pfs_matrix'] = list(s['pfs'])
    elif s['model_type'] == 'xyh2':
        solve_matrix_xyh2(s)
        s['energies_matrix'] = list(s['energies'])
        s['pfs_matrix'] = list(s['pfs'])
    charge = None
    if s['conserve'] == 'Q':
        charge = [0]
    eng = ExactDiagNonH(s.model, max_size=1E10, charge_sector=charge)
    eng.build_full_H_from_mpo()
    eng.full_diagonalization_nonH()
    # eng.E.sort()
    # eng.E = sorted(eng.E, key=lambda x:x.real)
    s.eng = eng
    s.eigs = []
    s.eigs_left = []
    s.states = []
    # s['energies'] = s.eng.E[s.eng.inds]
    s['E_0'] = s.eng.E[s.eng.inds[0]]
    s['e_0'] = s.eng.E[s.eng.inds[0]] / s['L']
    # s['e_0_real'] = (s['e_0']+0j).real
    s['energies'] = s.eng.E
    s['es'] = [x/s['L'] for x in s.eng.E]
    es = []
    for i in range(len(s.eng.E)):
        continue
        j = s.eng.inds[i]
        v = s.eng.V[:,j]
        # g = inner(s.eng.Vleft[:,j], v)
        st = State(j, s)
        s.states.append(st)
    for i in range(len(s.eng.E)):
        continue
        j = s.eng.inds[i]
        v = s.eng.V[:,j]
        # n = norm(v)
        # else:
        #     v = v / norm(v)
        e = s.eng.E[j]/s['L']
        s.eigs.append((e, v))
        es.append(e)
        st = State(j, s)
        st.psi = s.eng.full_to_mps(v)
        s.states.append(st)

        j = s.eng.indsleft[i]
        g = inner(v, s.eng.Vleft[:,j])
        g = inner(s.eng.Vleft[:,j], v)
        v = s.eng.Vleft[:,j]
        # v = v / norm(v)
        st = State(j, s)
        st.psi = s.eng.full_to_mps(v)
        s.states_left.append(st)
        e = s.eng.Eleft[j] / s['L']
        s.eigs_left.append((e, v))
    # s.eigs = sorted(s.eigs, key = lambda x:x[0].real)
    # s['energies'] = sorted(es, key = lambda x:x.real)
    # s['energies'] = sorted(es, key = lambda x: -abs(x))
    # s['energies'] = es
    # s['eigs'] = s.eigs
    # s['e_0'] = s.eigs[0][0]

    # s.eigs_left = []
    # for i in range(len(s.eng.Eleft)):
    #     v = s.eng.Vleft[:,i]
    #     n = norm(v)
    #     else:
    #         v = v / norm(v)
    #     s.eigs_left.append((s.eng.Eleft[i]/s['L'], v))
    # s.eigs_left = sorted(s.eigs_left, key = lambda x:x[0].real)
    # s['eigs_left'] = s.eigs_left
    # s['e_0_left'] = s.eigs_left[0][0]
    # for i in range(len(eng.E)):
    #     e1 = s.eng.E[i]
    #     e2 = s.spectrum[i]
    s.other_data[f'converged_0'] = 1
    if s['analyse_eigs']:
        res = analyse_eigs(s.eng.E, s.eng.Eleft, s.eng.V, s.eng.Vleft)
        m = {}
        for x, y in res.items():
            try:
                y = y.tolist()
                m[x] = y
            except:
                m[x] = y
        # print(res)
        # exit()
        s.measurement_data.update(m)
        s.save_measurements()
    s.save_other()

    if s['free_test2']:
        pfs1 = s["pfs"]
        pfs2 = s["pfs_matrix"]
        es1 = numpy.sort(s["energies"])
        es2 = numpy.sort(s["energies_matrix"])
        # for i in range(len(es1)):
        #     print(es1[i], es2[i])
        def f_kj_xy(k, g, L, u):
            return numpy.sin((L+2.) * k) + (((1.-g)/(1.+g)) ** u) * numpy.sin(L*k)
        prev = None
        for i in range(len(es1)):
            print(es1[i], es2[i])
        for p in pfs1:
            print(p)
        for p in pfs2:
            print(p)
        for p in pfs1:
            g = s["gamma_full"]
            k = eps_to_k_xy2(p, g)
            # print(p, k)
            # if abs(f) < 1E-6:
            #     if prev:
            #         print(k, prev - k)
            #     prev = k
            f = f_kj_xy(k, g, s["L"],-1)
            f = f_kj_xy(k, g, s["L"],+1)
            print(p, k)

        # print(pfs1)
        # print(pfs2)
        exit()
        for p in pfs2:
            print(p)
        import free_test
        f2 = free_test.free_test_newnew(s["energies"], s["N"], 0)
        for i in range(len(pfs1)):
            # print(pfs1[i], pfs2[i])
            print(pfs1[i])
        print("...")
        es1 = numpy.sort(s["energies"])
        es2 = numpy.sort(s["energies_matrix"])
        # for i in range(len(es1)):
        #     print(es1[i], es2[i])
        for x in f2:
            if len(x[1]) == 1:
                print(x[0]/2)
        exit()
        for x in f2:
            if len(x[1]) == 1:
                print(x[0] * 2)
        f2 = free_test.free_test_newnew(s["energies_matrix"], s["N"], 0)
        for x in f2:
            if len(x[1]) == 1:
                print(x[0])
        exit()

    if s['free_test']:
        import free_test
        f2 = free_test.free_test_newnew(s["energies"], s["N"], 0)
        pfs = []
        largest = 0
        for x, y in f2:
            if len(y) > largest:
                largest = len(y)
            if len(y) == 1:
                pfs.append(x*1.0)

        for x in f2:
            print(x)
        # if largest == s['L']:
        #     s['pfs'] = pfs
        #     pfm = s['pfs_matrix']
        #     print(pfs, pfm)
        #     for i in range(len(pfs)):
        #         print(pfs[i], pfm[i])
        # else:
        #     print(f"{s} doesn't have a free spectrum")
        #     s['pfs'] = None

        # def eps_to_k_xy(e, g):
        #     res = (1.-e**2) / (1.-g**2)
        #     res = numpy.arcsin(numpy.sqrt(res))
        #     return res

        def f_kj_xy(k, g, L, u):
            return numpy.sin((L+2.) * k) + (((1.-g)/(1.+g)) ** u) * numpy.sin(L*k)


        j = 1
        L = s["L"]
        pfs = s["pfs"]
        for p in pfs:
            print("_______")
            g = s['gamma_full']
            k = eps_to_k_xy(p, g)
            print(g)
            q = f_kj_xy(k, g, s['L'], 1)
            q2 = f_kj_xy(k, g, s['L'], -1)
            print(p)
            print(p*p)
            print(k)
            print(q, q2)

        exit()


def solve_dmrg(s, c, index):
    # logging.getLogger("tenpy.algorithms.dmrg").addHandler(s.handler)
    m = s['dmrg_method']
    psi = s.states[index].psi
    ortho = [x.psi for x in s.states[:index]]
    c['orthogonal_to'] = ortho
    eng = None
    hc = s['hc']
    if not hc:
        if s['active_sites'] == 1:
            eng = SSNonHermitianDMRG(psi, s.model, c)
        else:
            eng = NonHermitianDMRG(psi, s.model, c)
    else:
        if s['active_sites'] == 1:
            eng = SingleSiteDMRGEngine(psi, s.model, c)
        else:
            eng = TwoSiteDMRGEngine(psi, s.model, c)

    if m == 'ramp_up':
        chi = c.get('chi0', 5)
        chi_max = c.get('chi_max', 20)
        chi_mult = c.get('chi_mult', 1.5)
        while (True):
            if chi > chi_max:
                break
            c.options['trunc_params']['chi_max'] = 10
            c['max_sweeps'] = 1
            eng = NonHermitianDMRG(psi, s.model, c)
            eng.run()
            E = s.model.H_MPO.expectation_value(psi)
            sstats = eng.sweep_stats
            # h = s.model.H_MPO
            # p = psi.copy()
            # o = dict(c.options)
            # o['compression_method'] = 'SVD'
            # h.apply(p, o)
            # q = p.overlap(psi)

            # chi *= chi_mult
            # chi_new = int(chi*chi_mult)
            # if chi == chi_new: chi_new += 1
            # chi = chi_new
    else:
        a = eng.run()
        E = a[0]
        s['solve_energies'][index] = E
        sstats = eng.sweep_stats
        shelve = eng.shelve
        save_keys = ['time', 'sweep']
        for k in save_keys:
            s[k] = sstats[k][-1]
    # q = 5
    # for thing in [eng.update_stats, eng.sweep_stats]:
    #     for x, y in thing.items():
    # res = dmrg.run(psi, s.model, c)
    # E = res['E']
    s.other_data[f'converged_{index}'] = 1
    s.save_other()
    logging.getLogger("tenpy.algorithms.dmrg").removeHandler(s.handler)

def solve(s, index):
    c = Config(dict(s.config.options), "solver_config")
    for x in convert_to_bool:
        if x in c.options:
            c.options[x] = bool(c.options[x])
    method = s.get("method", "dmrg")
    if method in config_info:
        for key, vals in config_info[method].items():
            c.options[key] = {}
            for v in vals:
                if v in c.options:
                    c.options[key][v] = c.options[v]

    if method == "dmrg":
        psi = s.states[index].psi
        solve_dmrg(s, c, index)
    elif method == 'tebd':
        psi = s.states[index].psi
        if index > 0:
            return
        if c['bc'] == 1:
            return
        d = s['delta_tau_min']
        l = -round(np.log10(d))
        dtau_list = [10**(-x) for x in range(1, l+1)]
        c['order'] = s['tebd_order']
        c['delta_tau_list'] = dtau_list
        NNM = NearestNeighborModel.from_MPOModel(s.model)
        eng  = tebd.TEBDEngine(psi, NNM, c)
        z = eng.run_GS()
        E = s.model.H_MPO.expectation_value(psi)
        E = complex(E)
        if not s['bc'] == 'infinite':
            E = E / s['L']

        # s[f'e_{index}'] = E
        s.other_data[f'converged_{index}'] = 1
        s.save_other()

    elif method == "exact":
        solve_exact(s, c, index)
    elif method == "sparse":
        eng = ExactDiagNonH(s.model)
        eng.build_full_H_from_mpo()
        eng.partial_diag()
    elif method == "simple_ed":
        solve_matrix(s)
        s['energies_matrix'] = list(s['energies'])
        s['pfs_matrix'] = list(s['pfs'])
        solve_simple_ed(s)
    elif method == "ed2":
        solve_ed2(s)
    elif method == "ed3":
        solve_ed3(s)
    elif method == "none":
        # exact values only
        pass
    elif method == "poly":
        pfs_poly = sympy_eps2(s)
        pfs_poly = fix_pfs(pfs_poly, s['N'])
        s['pfs'] = pfs_poly
        s['pfs_real'] = [x.real for x in pfs_poly]
        s['pfs_imag'] = [x.imag for x in pfs_poly]
        P = PointGenerator(s)
        P.make_points()
        P.generate_points_gs()
        s['energies_generated'] = P.points
        s['energies'] = sorted([e for e in P.points])
        s['converged_0'] = 1
        s['e_0'] = s['energies'][0]
        s['E_0'] = s['energies'][0] * s['L']
    elif method == "matrix":
        if s['model_type'] == 'xyh':
            solve_matrix_xyh(s)
        elif s['model_type'] == 'xyh2':
            solve_matrix_xyh2(s)
        else:
            solve_matrix(s)
    elif method == 'matrix_limited':
        solve_matrix_limited(s)
    elif method == 'matrix_reusable':
        L = s["L"]
        N = s["N"]
        lamb = s["lambda"]
        phi = s["phi"]
        m1 = make_m1(s)
        sol = scipy.linalg.eig(m1)
        s['eigs'] = sorted(sol[0])
        s['converged_0'] = 1
    elif method == 'find_eps':
        solve_find_eps(s)

    else:
        print(f"Unknown method {method}")
        exit()

def solve_matrix_xyh(s):
    pfs, e0 = pfs_matrix_xyh(s)
    # s['pfs'] = pfs
    s['pfs'] = pfs
    s['pfs_real'] = [x.real for x in pfs]
    s['pfs_imag'] = [x.imag for x in pfs]
    s['e_0'] = e0 / s['L']
    s['E_0'] = e0
    s.other_data['converged_0'] = 1

    if s['depth']:
        s['energies'] = make_pf_spectrum2(e0, pfs, s['N'], depth=s['depth'])
    elif s['calc_spectrum'] or s['method'] == 'exact':
        spect = make_pf_spectrum2(e0, pfs, s['N'])
        s['energies'] = spect
    # spect = make_pf_spectrum2(e0, pfs, s['N'])
    # s['energies'] = spect

def solve_matrix_xyh2(s):
    pfs, e0 = pfs_matrix_xyh2(s)
    # s['pfs'] = pfs
    s['pfs'] = pfs
    s['pfs_real'] = [x.real for x in pfs]
    s['pfs_imag'] = [x.imag for x in pfs]
    s['e_0'] = e0 / s['L']
    s['E_0'] = e0
    s.other_data['converged_0'] = 1

    if s['depth']:
        s['energies'] = make_pf_spectrum2(e0, pfs, s['N'], depth=s['depth'])
    elif s['calc_spectrum'] or s['method'] == 'exact':
        spect = make_pf_spectrum2(e0, pfs, s['N'])
        s['energies'] = spect
    # spect = make_pf_spectrum2(e0, pfs, s['N'])
    # s['energies'] = spect
