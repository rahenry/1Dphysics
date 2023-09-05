#!/usr/bin/env python3
import numpy, scipy
import warnings
import random
from rah_utility import rel_diff
from fp_k import eps_to_k

from tenpy.linalg import np_conserved as npc
from tenpy.linalg import sparse
# from .site import GroupedSite, group_sites
from tenpy.tools.misc import to_iterable, to_array, get_recursive
from tenpy.tools.math import lcm, speigs, entropy
from tenpy.tools.params import asConfig
from tenpy.tools.cache import DictCache
from tenpy.tools import hdf5_io
from tenpy.algorithms.truncation import TruncationError, svd_theta
from tenpy.networks.mps import MPSEnvironment
from tenpy.linalg.np_conserved import tensordot, inner, norm
from System import System

def ortho_test(s):
    vs = []
    for z1 in s.zero_states:
        psi1 = s.eng.V[:,z1].to_ndarray()
        vs.append(psi1)
        for z2 in s.zero_states:
            psi2 = s.eng.V[:,z2].to_ndarray()
            a = numpy.inner(psi1.conj(), psi2)
            # print(z1, z2, abs(a))

    vs = numpy.array(vs).T
    vs, R = numpy.linalg.qr(vs)
    print('...')
    print(vs.shape)
    for z1 in s.zero_states:
        psi1 = s.eng.V[:,z1].to_ndarray()
        print(psi1.shape)
        a = 0
        for v in range(vs.shape[-1]):
            psi2 = vs[:,v]
            a += abs(numpy.vdot(psi1, psi2))**2
        print(z1, a)


    
def get_zero_states(s):
    s.zero_states = []
    i = 0
    tol = 1E-4
    # print(s)
    for e in s['energies']:
        # print(abs(e))
        # if abs(e) < 1E-2:
        #     print(abs(e))
        if abs(e) < tol:
            s.zero_states.append(i)
        i += 1

    tol2 = 1E-12
    count2 = 0
    for e in s['energies']:
        if abs(e) < tol2:
            count2 += 1
    print(s, len(s.zero_states), count2)

def get_fidel(phi, srefs, srefs_conj=[]):
    phi = phi.to_ndarray()
    res = 0
    vs = []
    for sref in srefs:
        for x in sref.zero_states:
            psi = sref.eng.V[:,x]
            psi = psi.to_ndarray()
            vs.append(psi)
        vs_test = scipy.linalg.orth(numpy.array(vs).T)
        print(vs_test.shape)
    for sref in srefs_conj:
        for x in sref.zero_states:
            psi = sref.eng.V[:,x]
            psi = psi.to_ndarray()
            psi = numpy.conj(psi)
            vs.append(psi)
        vs_test = scipy.linalg.orth(numpy.array(vs).T)
        print(vs_test.shape)
    vs = numpy.array(vs).T
    vs = scipy.linalg.orth(vs, rcond=1E-10)
    print(vs.shape)
    for i in range(vs.shape[-1]):
    # for psi in vs:
        # x = sref.zero_states[i]
        # psi = sref.eng.V[:,x]
        psi = vs[:,i]
        # a = inner(phi, psi, do_conj=1, axes='range')
        a = numpy.inner(phi.conj(), psi)

        # print(x, a)
        res += abs(a) ** 2
    return res

def count_zeros(S):
    for s in S.slicer:
        get_zero_states(s)
        s['nzeros'] = len(s.zero_states)

def periodic_lambda_test(S):
    s0 = S.systems[0]
    if not s0.will_run:
        return
    # search = {'lambda':0.5}
    lambdas = S.slicer.get_possible_values('lambda')
    search = {'lambda' : lambdas[-1]}
    sref = S.slicer.find(search)[0]
    for s in S.slicer:
        get_zero_states(s)
    # ortho_test(sref)
    for s in S.slicer:
        # c = dict(s.config_base)
        # c.update({'conj':1, 'lambda':numpy.conj(s['lambda'])})
        # print(s, c)
        # sleft = System(c, just_go=1)
        # get_zero_states(sleft)

        fids = []
        for z in sref.zero_states:
            # phi = s.eng.V[:,z]
            # f = get_fidel(phi, sref)
            # fids.append(f)
            phi = sref.eng.V[:,z]
            # f = get_fidel(phi, [s], [sleft])
            f = get_fidel2(phi, [s])
            fids.append(f)
        s.measurement_data['zero_fidelities'] = fids
        s.save_measurements()

def get_fidel2(phi, srefs, srefs_conj=[]):
    phi = phi.to_ndarray()
    res = 0
    vs = []
    for sref in srefs:
        for x in sref.zero_states:
            psi = sref.eng.V[:,x]
            psi = psi.to_ndarray()
            vs.append(psi)
        # vs_test, R = numpy,linalg.qr.(numpy.array(vs).T)
        # print(vs_test.shape)
    vs = numpy.array(vs).T
    vs, R = numpy.linalg.qr(vs)
    a, b, c = numpy.linalg.svd(vs)
    print(len(b), b)

        # vs_test, R = numpy,linalg.qr.(numpy.array(vs).T)
    # print(vs.shape)
    for i in range(vs.shape[-1]):
    # for psi in vs:
        # x = sref.zero_states[i]
        # psi = sref.eng.V[:,x]
        psi = vs[:,i]
        # a = inner(phi, psi, do_conj=1, axes='range')
        # a = numpy.inner(phi.conj(), psi)
        a = numpy.vdot(phi, psi)

        # print(x, a)
        res += abs(a) ** 2
    return res

def simple_zero_test(S):
    for s in S.systems:
        get_zero_states(s)
        print(s, len(s['zero_states']))

    print(S)
    exit()
