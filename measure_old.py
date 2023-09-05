#!/usr/bin/env python3
import numpy

from measure_funcs import *
from tenpy.linalg.np_conserved import tensordot, inner, norm, outer

def measure_MPS(s):
    s.drop_charges()
    L = s['L']
    N = s['N']
    for index in range(s['n_states']):
        psiR = s.chargeless_states[index]
        psiL = s.make_left_state(psiR)
        psiR.canonical_form()
        psiL.canonical_form()
        psiR.norm = 1
        psiL.norm = 1
        E = s.chargeless_model.H_MPO.expectation_value(psiR)
        if not s['bc'] == 'infinite':
            E = E / s['L']

        s[f'e_{index}'] = E
        if s['e_exact']:
            s['e0_error'] = rel_diff(s['e_0'], s['e_exact'])
        if s['max_hours']:
            s['time_relative'] = s['time'] / s['max_hours'] / 3600.

        o = psiL.overlap(psiR)
        s['o'] = abs(o)
        s[f'norm_{index}'] = abs(o)
        s['oL'] = abs(o) / L

        Z1L = correlation_function(psiL, psiR, "Z", "ZP", [0], [s['L']-1], do_conj=1)[0][0]
        s['Z1L'] = abs(Z1L / o)
        X1L = correlation_function(psiL, psiR, "X", "XP", [0], [s['L']-1], do_conj=1)[0][0]
        s['X1L'] = abs(X1L / o)
        s['z'] = s['Z1L'] / s['L']
        s['x'] = s['X1L'] / s['L']

        Z2L = correlation_function(psiL, psiR, "Z", "ZP", [0], [s['L']//2], do_conj=1)[0][0]
        s['Z2L'] = abs(Z2L / o)
        X2L = correlation_function(psiL, psiR, "X", "XP", [0], [s['L']//2], do_conj=1)[0][0]
        s['X2L'] = abs(X2L / o)
        s['z2'] = s['Z2L'] / s['L']
        s['x2'] = s['X2L'] / s['L']

        st = L // 2
        Ztest = correlation_function(psiL, psiR, "Z", "ZP", [st], [st+1], do_conj=1)[0][0]
        s['rhoZtest'] = abs(Ztest / o)
        Xtest = correlation_function(psiL, psiR, "X", "XP", [st], [st+1], do_conj=1)[0][0]
        s['rhoXtest'] = abs(Xtest / o)
        # s['rhoZtest'] = s['rhoZtest'] / s['L']
        # s['rhoXtest'] = s['rhoXtest'] / s['L']
        omega = numpy.exp(2.j * numpy.pi / N)
        z = numpy.zeros([N, N])
        z[0][0] = 1
        z[1][1] = omega
        z[2][2] = omega**2

        phi = psiR.copy()
        z = s.chargeless_model.lat.mps_sites()[st].get_op("Z")
        phi.apply_local_op(st, z)
        # z = s.chargeless_model.lat.mps_sites()[st+1].get_op("ZP")
        # phi.apply_local_op(st+1, z)
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
        s['x1Ltest'] = abs(u/o)

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
        s['z1Ltest'] = abs(u/o)


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
        # print(o3)
        env = MPSEnvironment(psiL, psiR)
        s['rhoX'] = abs(env.expectation_value([o3], sites=[st])[0]/o)

        fp = s['fidelity_parameter']
        if fp:
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

        m = s['n_correlation_sites']
        if m is not None:
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
                # print(o3)
                ex = env.expectation_value([o3], sites=[st])[0]
                s[f'XX{i}'] = abs(ex/o)
                # s[f'phiXX{i}'] = numpy.angle(ex)/o


def lambda_scale(s):
    i = 0
    print(i, type(i))
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
        print(f, z1, z2)
        f = f / z1 / z2
        d = s['fidelity_delta']
        chi = (1.-f) / d/d
        s['fidelity'] = f
        s['fsusceptibility'] = chi

    z = s.model.lat.mps_sites()[2].get_op("Z")
    print(z)
    exit()

def measure(s):

    method = s['method']
    if method == 'exact' or method == 'sparse':
        measure_exact(s)
    else:
        measure_MPS(s)
    lambda_scale(s)

    s.other_data['measured'] = True
    s.save_other()
    s.save_measurements()
