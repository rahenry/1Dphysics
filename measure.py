#!/usr/bin/env python3
import numpy

from measure_funcs import *
from tenpy.linalg.np_conserved import tensordot, inner, norm, outer

def measure_sxx(s):
    for i in range(s['n_states']):
        index = s['eig_indices'][i]
        psi = s.states[index].psi
        E = s.model.H_MPO.expectation_value(psi)
        s[f'E_{i}'] = E
        if not s['bc'] == 'infinite':
            E = E / s['L']

        s[f'e_{i}'] = E
        if i > 0:
            s[f'gap_{i}'] = E - s['e_0']


    i = 0
    index = s['eig_indices'][i]
    s['e_error'] = rel_diff(s['e_0'], s['e_exact'])
    psi = s.states[index].psi
    mag_z = numpy.sum(psi.expectation_value("Z"))
    mag_x = numpy.sum(psi.expectation_value("X"))
    mag_y = numpy.sum(psi.expectation_value("Y"))
    s['mag_z'] = mag_z
    s['mag_x'] = mag_x
    s['mag_y'] = mag_y
    L = s['L']
    con_z = psi.correlation_function("Z", "Z", [L//2], [L//2+1])[0]
    s['con_z'] = con_z
    con_x = psi.correlation_function("X", "X", [L//2], [L//2+1])[0]
    s['con_x'] = con_x
    con_y = psi.correlation_function("Y", "Y", [L//2], [L//2+1])[0]
    s['con_y'] = con_y
    E = s.model.H_MPO.expectation_value(psi)
    s['entropy'] = psi.entanglement_entropy(1, [s['L']//2])[0]
    if s['e_exact']:
        s['e0_error'] = rel_diff(s['e_0'], s['e_exact'])

    if s['max_hours'] and s['time']:
        s['time_relative'] = s['time'] / s['max_hours'] / 3600.

    fp = s['fidelity_parameter']
    if fp or s.parent["fidelity_parameters"]:
        fs = s.auxillary_systems[0]
        psi2 = fs.states[fs['eig_indices'][0]].psi
        E = s.model.H_MPO.expectation_value(psi)
        E2 = fs.model.H_MPO.expectation_value(psi2)

        f = psi2.overlap(psi) * psi.overlap(psi2)
        d = s['fidelity_delta']
        z1 = psi.overlap(psi)
        z2 = psi2.overlap(psi2)
        # f /= numpy.sqrt(z1*z2)
        f = f / z1 / z2
        chi = (1.-f) / d/d
        s['fidelity'] = f
        s['fsusceptibility'] = chi
        s['fsL'] = chi / s['L']
        Ef = s.model.H_MPO.expectation_value(psi2) / s['L']
        s['fidelity_delta_e'] = s['e_0'] - Ef
    # for index in range(s['n_states']):
    return
def measure_MPS(s):
    if s['model_type'] == 'sxx':
        measure_sxx(s)
        return
    s.drop_charges()
    L = s['L']
    N = s['N']
    for index in range(s['n_states']):
        psiR = s.chargeless_states[index]
        # psiL = s.chargeless_states_left[index]
        psiL = s.make_left_state(psiR)
        psiR.canonical_form()
        psiL.canonical_form()
        psiR.norm = 1
        psiL.norm = 1
        E = s.chargeless_model.H_MPO.expectation_value(psiR)
        s[f'E_{index}'] = E
        if not s['bc'] == 'infinite':
            E = E / s['L']

        s[f'e_{index}'] = E
        E2 = s.model.H_MPO.expectation_value(s.states[index].psi) / s['L']
        if s['e_exact']:
            s['e0_error'] = rel_diff(s['e_0'], s['e_exact'])
        if s['max_hours'] and s['time']:
            s['time_relative'] = s['time'] / s['max_hours'] / 3600.

        o = psiL.overlap(psiR)
        s['o'] = abs(o)
        s[f'norm_{index}'] = abs(o)
        s['oL'] = abs(o) / L
        s[f'o{index}'] = abs(o)

        fp = s['fidelity_parameter']
        if fp or s.parent["fidelity_parameters"]:
            fs = s.auxillary_systems[0]
            psiR2 = fs.chargeless_states[index]
            psiL2 = fs.make_left_state(psiR2)

            f = psiL2.overlap(psiR) * psiL.overlap(psiR2)
            d = s['fidelity_delta']
            z1 = psiL.overlap(psiR)
            z2 = psiL2.overlap(psiR2)
            # f /= numpy.sqrt(z1*z2)
            f = f / z1 / z2
            chi = (1.-f) / d/d
            s['fidelity'] = f
            s['fsusceptibility'] = chi
            s['fsL'] = chi / s['L']

        if not s['N'] == 3:
            return
        st = L // 2
        omega = numpy.exp(2.j * numpy.pi / N)
        z = numpy.zeros([N, N], dtype=complex)
        z[0][0] = 1
        z[1][1] = omega
        z[2][2] = omega**2

        phi = psiR.copy()
        z = s.chargeless_model.lat.mps_sites()[st].get_op("Z")
        phi.apply_local_op(st, z)
        u = psiL.overlap(phi)
        s['rhoZtest2'] = abs(u/o)

        phi = psiR.copy()
        x = s.chargeless_model.lat.mps_sites()[st].get_op("X")
        phi.apply_local_op(st, x)
        x = s.chargeless_model.lat.mps_sites()[st+1].get_op("XP")
        phi.apply_local_op(st+1, x)
        u = psiL.overlap(phi)
        s['rhoXtest2'] = abs(u/o)

        phi = psiR.copy()
        x = s.chargeless_model.lat.mps_sites()[0].get_op("X")
        phi.apply_local_op(0, x)
        x = s.chargeless_model.lat.mps_sites()[L-1].get_op("XP")
        phi.apply_local_op(L-1, x)
        u = psiL.overlap(phi)
        s['x1L'] = abs(u/o)
        s['x'] = abs(u/o)/L

        phi = psiR.copy()
        x = s.chargeless_model.lat.mps_sites()[0].get_op("X")
        phi.apply_local_op(0, x)
        x = s.chargeless_model.lat.mps_sites()[st].get_op("XP")
        phi.apply_local_op(st, x)
        u = psiL.overlap(phi)
        s['xhalf'] = abs(u/o)

        phi = psiR.copy()
        z = s.chargeless_model.lat.mps_sites()[0].get_op("Z")
        phi.apply_local_op(0, z)
        z = s.chargeless_model.lat.mps_sites()[st].get_op("ZP")
        phi.apply_local_op(st, z)
        u = psiL.overlap(phi)
        s['zhalf'] = abs(u/o)

        phi = psiR.copy()
        z = s.chargeless_model.lat.mps_sites()[0].get_op("Z")
        phi.apply_local_op(0, z)
        z = s.chargeless_model.lat.mps_sites()[L-1].get_op("ZP")
        phi.apply_local_op(L-1, z)
        u = psiL.overlap(phi)
        s['z1L'] = abs(u/o)
        s['z'] = abs(u/o) / L


        phi = psiR.copy()
        x = s.chargeless_model.lat.mps_sites()[0].get_op("X")
        phi.apply_local_op(0, x)
        x = s.chargeless_model.lat.mps_sites()[L-1].get_op("XP")
        phi.apply_local_op(L-1, x)
        u = psiL.overlap(phi)
        s['x1Ltest'] = abs(u/o)

        st = s['L'] // 2
        op1 = s.chargeless_model.lat.mps_sites()[st].get_op("Z")
        env = MPSEnvironment(psiL, psiR)
        ex = env.expectation_value([op1], sites=[st])
        s['rhoZ'] = abs(ex[0]/o)

        o1 = s.model.lat.mps_sites()[st].get_op("X").drop_charge()
        o2 = s.model.lat.mps_sites()[st+1].get_op("XP").drop_charge()
        o1 = o1.replace_labels(['p', 'p*'], ['p0', 'p0*'])
        o2 = o2.replace_labels(['p', 'p*'], ['p1', 'p1*'])
        o3 = outer(o1, o2)
        env = MPSEnvironment(psiL, psiR)
        s['rhoX'] = abs(env.expectation_value([o3], sites=[st])[0]/o)


        m = s['n_correlation_sites']
        if not m == []:
            if m > L: m = L
            if m < 0: m = L
            env = MPSEnvironment(psiL, psiR)
            sites = [int(i) for i in numpy.linspace(0, L-1, m)]
            for i in sites:
                Z = s.chargeless_model.lat.mps_sites()[i].get_op("Z")
                ex = env.expectation_value([Z], sites=[i])[0]
                s[f'Z{i}'] = abs(ex/o)
                # s[f'phiZ{i}'] = numpy.angle(ex)/o

            sites = [int(i) for i in numpy.linspace(0, L-2, m)]
            for i in sites:
                st = i
                o1 = s.model.lat.mps_sites()[st].get_op("X").drop_charge()
                o2 = s.model.lat.mps_sites()[st+1].get_op("XP").drop_charge()
                o1 = o1.replace_labels(['p', 'p*'], ['p0', 'p0*'])
                o2 = o2.replace_labels(['p', 'p*'], ['p1', 'p1*'])
                o3 = outer(o1, o2)
                ex = env.expectation_value([o3], sites=[st])[0]
                s[f'XX{i}'] = abs(ex/o)
                # s[f'phiXX{i}'] = numpy.angle(ex)/o


def lambda_scale(s):
    i = 0
    while (True):
        e = s[f'e_{i}_real']
        if not e: break
        s[f'e_{i}_lambda_real'] = s[f'e_{i}_real'] * s['lambda_scale']
        s[f'e_{i}_lambda_imag'] = s[f'e_{i}_imag'] * s['lambda_scale']
        i += 1

def invert(psi):
    psi = psi.copy()
    # psi.sites = psi.sites[::-1]
    psi.form = [(f if f is None else (f[1], f[0])) for f in psi.form[::-1]]
    psi._B = [
        B.replace_labels(['vL', 'vR'], ['vR', 'vL']).transpose(psi._B_labels)
        for B in psi._B[::-1]
    ]
    psi._S = psi._S[::-1]
    psi.canonical_form()
    return psi

def conjify(psi):
    psi = psi.copy()
    # psi.sites = psi.sites[::-1]
    psi._B = [B.conj().conj(complex_conj=False) for B in psi._B]
    psi.canonical_form()
    return psi

def make_left(psi, s):
    N = s['N']
    z = np.zeros([N, N])
    z[0][0] = 1
    z[2][1] = 1
    z[1][2] = 1
    o = npc.Array.from_ndarray_trivial(z, labels=['p', 'p*'])
    ops = [o for i in range(s["L"])]
    psi = psi.copy()
    psi.apply_product_op(ops)
    psi = invert(psi)
    return psi


def measure_exact(s):
    s.eigs = s['eigs']
    for i in range(s['n_states']):
        s[f'e_{i}'] = s.eigs[i][0] / s['L']
        if s.left_system:
            s[f'e_{i}_left'] = s.left_system.eigs[i][0] / s['L']
    i = 0
    fp = s['fidelity_parameter']
    index = 0
    psiR = s.eigs[index][1]
    psiL = psiR
    c = False
    c = 0
    if s.left_system:
        psiL = s.left_system.eigs[index][1]
    s['o'] = inner(psiL, psiR, do_conj=c)
    if fp:
        psiR2 = s.fidelity_system.eigs[index][1]
        psiL2 = psiR2
        if s.left_system:
            psiL2 = s.fidelity_system.left_system.eigs[index][1]
        f = inner(psiL2, psiR, do_conj=c) * inner(psiL, psiR2, do_conj=c)
        z1 = inner(psiL, psiR, do_conj=c)
        z2 = inner(psiL2, psiR2, do_conj=c)
        f = f / z1 / z2
        d = s['fidelity_delta']
        chi = (1.-f) / d/d
        s['fidelity'] = f
        s['fsusceptibility'] = chi

    z = s.model.lat.mps_sites()[2].get_op("Z")
    exit()

def get_min_pf_dist0n(s):
    pfs = s['pfs']
    if len(pfs) == 1:
        return 1
    N = s['N']
    L = s['L']
    pfs = sorted(pfs, key=lambda x: abs(x))
    pfs_abs = [abs(x) for x in pfs]
    scale = max(pfs_abs)
    d0 = 10000
    best = None
    # for i in range(L-1):
    for i in range(len(pfs)-1):
        p0 = pfs[0]
        p1 = pfs[i+1]
        for k in range(N):
            p1a = p1 * numpy.exp(2.j * numpy.pi * k / N)
            d = abs(p0 - p1a)
            if d < d0:
                d0 = d
                yt1 = float(k)
                best = (p0, p1)
    return d0/scale

def get_min_pf_dist0(s):
    pfs = s['pfs']
    if len(pfs) == 1:
        return 1
    N = s['N']
    L = len(pfs)
    L = s['L']
    pfs = sorted(pfs, key=lambda x: abs(x))
    pfs_abs = [abs(x) for x in pfs]
    scale = max(pfs_abs)
    scale = 1
    offset = 0
    # if abs(pfs[0]) < 1E-6:
    #     offset = 1


    d0 = 10000
    best = None
    # for i in range(L-1):
    for i in range(len(pfs)-1-offset):
        p0 = pfs[0+offset]
        p1 = pfs[i+1+offset]
        for k in range(N):
            p1a = p1 * numpy.exp(2.j * numpy.pi * k / N)
            d = abs(p0 - p1a)
            if d < d0:
                d0 = d
                yt1 = float(k)
                best = (p0, p1)
    return d0

def get_min_pf_dist_minus(s):
    pfs = s['pfs']
    if len(pfs) == 1:
        return 1
    N = s['N']
    L = len(pfs)
    L = s['L']
    pfs = sorted(pfs, key=lambda x: abs(x))
    pfs_abs = [abs(x) for x in pfs]
    scale = max(pfs_abs)
    scale = 1

    d0 = 10000
    d1 = 0
    best = None
    # for i in range(L):
    for i in range(len(pfs)):
        # for j in range(L):
        for j in range(len(pfs)):
            if i == j:
                continue
            p0 = pfs[i]
            p1 = pfs[j]
            for k in range(N):
                p1a = p1 * numpy.exp(2.j * numpy.pi * k / N)
                d = abs(p0 - p1a) / scale
                if d < d0:
                    d0 = d
                    yt1 = float(k)
                    best = (p0, p1)
                if d > d1:
                    d1 = d
    return d0

def get_min_pf_dist_all(s):
    pfs = s['pfs']
    N = s['N']
    L = s['L']

    best = 10000
    if len(pfs) == 1:
        return 1
    for i in range(len(pfs)):
        for j in range(len(pfs)):
            if i == j:
                continue
            p0 = pfs[i]
            p1 = pfs[j]
            for k in range(N):
                p1a = p1 * numpy.exp(2.j * numpy.pi * k / N)
                d = abs(p0 - p1a)
                if d < best and d > 1E-10:
                    best = d
    return best


def get_min_pf_dist(s):
    pfs = s['pfs']
    N = s['N']
    L = s['L']
    pfs = sorted(pfs, key=lambda x: abs(x))
    pfs_abs = [abs(x) for x in pfs]
    scale = max(pfs_abs)
    scale = 1


    d0 = 10000
    d1 = 0
    best = None
    # for i in range(L):
    if len(pfs) == 1:
        return 1
    for i in range(len(pfs)):
        # for j in range(L):
        for j in range(len(pfs)):
            if i == j:
                continue
        # pfs = s['pfs']
            # if i >= j:
            #     continue
            p0 = pfs[i]
            p1 = pfs[j]
            for k in range(N):
                p1a = p1 * numpy.exp(2.j * numpy.pi * k / N)
                d = abs(p0 - p1a) / scale
                if d < d0:
                    d0 = d
                    yt1 = float(k)
                    best = (p0, p1)
                if d > d1:
                    d1 = d
    return d0 / d1


def get_min_E_dist(s):
    # if not s['calc_spectrum']:
    #     return 0
    d = 100000
    es = s['energies']
    if es is None or es == []:
        return 0
    for i in range(len(s['energies'])):
        for j in range(len(s['energies'])):
            if i == j: continue
            e1 = s['energies'][i]
            e2 = s['energies'][j]
            a = abs(e1-e2)
            if a < d:
                d = a
    return d

def get_q(pfs, N):
    L = len(pfs)
    pfs = sorted(pfs, key=lambda x: abs(x))
    d0 = 10000
    best = (0, 0)
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
                best = (p0, p1)

    return d0
    p0 = pfs[0]
    q = 10000
    for p in pfs[1:]:
        a = abs(p0 - p)
        if a < q:
            q = a

def measure_simple_ed(s):
    s['min_E_distance'] = get_min_E_dist(s)
    N = s['N']
    # s['min_pf_distance'] = get_q(pfs, N)
    L = s['L']
    return
    for i in range(N ** L):
        v1 = s.eigenvectors[i]
        v2
        1
    exit()

def measure_pfs(s):
    pfs = s['pfs']
    if len(pfs) == 0:
        return
    N = s['N']
    # s['min_pf_distance'] = get_q(pfs, N)
    s['min_pf_distance0'] = get_min_pf_dist0(s)
    if s["measure_all_dist"]:
        s['min_pf_distance'] = get_min_pf_dist(s)
        s['min_pf_distance0n'] = get_min_pf_dist0n(s)
        s['min_pf_distance_all'] = get_min_pf_dist_all(s)
        s['min_pf_distance_minus'] = get_min_pf_dist_minus(s)
        s['min_pf_distance_test'] = get_min_pf_dist_all(s) - get_min_pf_dist0(s)
        L = s['L']
        if (L < 10):
            s['min_E_distance'] = get_min_E_dist(s)
        else:
            s['min_E_distance'] = 1
    if s['model_type'] == 'fp':
        ks = [eps_to_k(p, s['lambda_full'], s['N']) for p in s['pfs']]
        ks = sorted(ks)
        s['ks'] = ks



def measure_XYW(s):
    pfs = s['pfs']
    N = s['N']
    L = s['L']
    # s['min_pf_distance'] = get_q(pfs, N)
    s['min_pf_distance'] = get_min_pf_dist(s)
    s['min_pf_distance0'] = get_min_pf_dist0(s)
    s['min_pf_distance_minus'] = get_min_pf_dist_minus(s)
    L = s['L']
    if (L < 10):
        s['min_E_distance'] = get_min_E_dist(s)
    else:
        s['min_E_distance'] = 1
    ks = [eps_to_k(p, s['lambda_full'], s['N']) for p in s['pfs']]
    ks = sorted(ks)
    s['ks'] = ks

def measure_SICP(s):
    L = s['L']
    if (L < 10):
        s['min_E_distance'] = get_min_E_dist(s)
    else:
        s['min_E_distance'] = 1

def measure(s):
    method = s['method']
    if method == 'poly' or method == 'matrix' or method == 'exact':
        if s['model_type'] == 'fpXYW':
            measure_XYW(s)
        elif s['model_type'] == 'SICP':
            measure_SICP(s)
        else:
            measure_pfs(s)

    if method == 'simple_ed':
        measure_simple_ed(s)
    # if method == 'exact' or method == 'sparse':
    #     measure_exact(s)
    # else:
    #     measure_MPS(s)
    if method in ['dmrg']:
        measure_MPS(s)
        lambda_scale(s)

    if method == 'ed3':
        s.ed3.measure()
    s.other_data['measured'] = True
    s.save_other()
    s.save_measurements()
