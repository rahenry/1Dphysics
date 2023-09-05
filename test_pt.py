#!/usr/bin/env python3

import numpy
import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=FutureWarning)
warnings.filterwarnings(action='ignore', category=numpy.ComplexWarning)
import math
import scipy
from measure_funcs import *
from SystemSet import SystemSet
from exact_fp import *
from decimal import *
from scipy.special import hyp2f1
from scipy.optimize import root
from scipy.optimize import brentq
from scipy.optimize import fsolve
from rah_utility import mkdir, silent_remove
from rah_numpy import eigenvalue_test
import os
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.linalg.np_conserved import tensordot, inner, norm, outer

# trying to do PT and find left eigenstate

def t2(psi):
    psi = psi.copy()
    for i in range(psi.L):
        site = psi.sites[i]
        print(site)
        sites = psi.sites
        print(sites)
        q = site.Z
        # v = site.plocal
        o = site.get_op('Z')
        # o = site.get_op('plocal')

        1
def transmogrify(psi):
    psi = psi.copy()
    psi = psi.spatial_inversion()
    return psi
    Bs = []
    for q in range(psi.L):
        b = psi.get_B(q)
        a = b.to_ndarray()[:,[0,2,1],:]
        legs = b.legs
        legs[0].qconj *= -1
        legs[2].qconj *= -1
        bnew = npc.Array.from_ndarray(a, legs, labels=b._labels, qtotal=b.qtotal)
        Bs.append(bnew)
    svs = list(psi._S)
    svs.reverse()
    print(svs)
    print(type(svs))
    snew = MPS(s.model.lat.mps_sites(), Bs, svs)
    return snew

def remove_charges():
    psi = psi.copy()
    # psi.sites = psi.sites[::-1]
    psi._B = [B.drop_charge() for B in psi._B]
    psi.canonical_form()
    return psi

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
    #

    Bs = []
    for q in range(psi.L):
        # b = psi.get_B(q)
        b = psi._B[q]
        Bs.append(b)
        continue

        # b = psi.get_B(q)
        b = psi._B[q]
        if (N==3):
            a = b.to_ndarray()[:,[0,2,1],:]
        else:
            a = b.to_ndarray()
        legs = b.legs
        print(legs)
        legs[0].qconj *= -1
        legs[2].qconj *= -1
        bnew = npc.Array.from_ndarray(a, legs, labels=b._labels, qtotal=b.qtotal)
        Bs.append(bnew)
    svs = []
    for q in range(psi.L+1):
        sv = psi._S[q]
        # if len(sv) > 1:
        #     print(sv)
        #     if (N==3):
        #         sv = sv[[0,2,1]]
        svs.append(sv)

    snew = MPS(psi.sites, Bs, svs)
    # psi.test_sanity()
    return psi

def conjify(psi):
    psi = psi.copy()
    # psi.sites = psi.sites[::-1]
    psi._B = [B.conj().conj(complex_conj=False) for B in psi._B]
    psi.canonical_form()
    return psi
S = SystemSet("pttest1")
for s in S.systems:
    print(s)
    N = s['N']
    s.run()
    s.measure()
    s.init_model()
    s.load_states()
    s.drop_charges()
    s1 = s.chargeless_states[0]
    s1.canonical_form()
    s2 = s.make_left_state(s1)
    # s2 = s.left_system.chargeless_states[0]
    s2.canonical_form()
    st = s['L'] // 2
    print(s['rhoZ_exact'])
    print(s['rhoX_exact'])
    def f1(a, b):
        things = [['X', 'XP'],
                  ['XP', 'X'],
                  ['Z', 'ZP'],
                  ['ZP', 'Z'],
                  ]
        for t in things:
            print(t)
            d = expectation_value_multi_sites(a, b, t, st)
            c = d / a.overlap(b)
            print(abs(d), abs(c))
    # f1(s2, s1)

    psiL = s2
    psiR = s1
    c = psiL.overlap(psiR)
    e = MPSEnvironment(psiL, psiR)
    # print(c, e.full_contraction(0))
    c = e.full_contraction(0)
    st = s['L'] //2
    o1 = s.model.lat.mps_sites()[st].get_op("ZP").drop_charge()
    ops = [o1]
    ex = e.expectation_value(ops, sites=[st])
    # print(abs(ex[0]/c))
    o1 = s.model.lat.mps_sites()[st].get_op("Z").drop_charge()
    ops = [o1]
    ex = e.expectation_value(ops, sites=[st])
    print(abs(ex[0]/c))
    o1 = s.model.lat.mps_sites()[st].get_op("X").drop_charge()
    o2 = s.model.lat.mps_sites()[st+1].get_op("XP").drop_charge()
    o1 = o1.replace_labels(['p', 'p*'], ['p0', 'p0*'])
    o2 = o2.replace_labels(['p', 'p*'], ['p1', 'p1*'])
    o3 = outer(o1, o2)
    # print(o3)
    ex = e.expectation_value([o3], sites=[st])
    print(abs(ex[0]/c))
exit()
print(abs(s1.overlap(s2)))
s2 = s.left_system.chargeless_states[0]
l = s2
s4 = s1.copy()
# ops = ["P1", "P2"]
# ops = [f'P' for i in range(s["L"])]
# s4.apply_product_op(ops)
o = s1.sites[1].get_op("Z")
print(o)
print(type(o))
print(o)
N = s['N']
z = np.zeros([N, N])
z[0][0] = 1
z[2][1] = 1
z[1][2] = 1
o = npc.Array.from_ndarray_trivial(z, labels=['p', 'p*'])
ops = [o for i in range(s["L"])]
s4.apply_product_op(ops)
c = s1.overlap(conjify(s1))
c = s1.overlap(conjify(conjify(s1)))
print(123)
a = abs(s1.overlap(conjify(s1)))
print(a)
a = abs(l.overlap(s1))
print(a)
a = abs(l.overlap(conjify(s1)))
print(a)
a = abs(l.overlap(s4))
print(a)
st = s['L'] // 2
sz = [st, st+1]
b = expectation_value_multi_sites(s2, s1, ['X', 'XP'], st)
print(b, abs(b))
b = expectation_value_multi_sites(s2, s1, ['X', 'XP'], st)
print(b, abs(b))
b = expectation_value_multi_sites(s2, s1, ['X', 'XP'], st)
print(b, abs(b))
b = expectation_value_multi_sites(s2, s1, ['XP', 'X'], st)
print(b, abs(b))
exit()
# s4 = invert(s4)
o12 = s2.overlap(s1)

print(s['rhoX'], s['rhoX_exact'])
def rhoX(psiL, psiR):
    psiL.canonical_form()
    psiR.canonical_form()
    # rhoX = correlation_function(psiL, psiR, "X", "XP", [s['L']//2], [s['L']//2+1], do_conj=1)[0][0]
    a = s['L']//2
    b = a+1
    print(a, b)
    # rhoX = correlation_function(psiL, psiR, "X", "XP", [0], [s['L']-2], do_conj=1)[0][0]
    rhoX = correlation_function(psiL, psiR, "X", "XP", [a], [b], do_conj=1)[0][0]
    print(rhoX)
    return abs(rhoX / psiL.overlap(psiR))
print(rhoX(s4, s1))
print(rhoX(s1, s1))
print('..')
print(rhoX(s2, s1))
print(rhoX(s2, s1))
print(rhoX(s2, s1))
print('..')
print(rhoX(s1, s1))
print(rhoX(s1, s1))
print(rhoX(s1, s1))
print(abs(o12))
print(s4.overlap(s1))
print(abs(s4.overlap(s2)))
exit()
# s3 = s2.spatial_inversion()
# s3 = transmogrify(s2, sref)
q = 1
b0 = s3.get_B(q);
bref = s1.get_B(q);
# b = b0.to_ndarray()
b = b0
# b = b[:,[2,1,0],:]
# b = b[:,:,[0,2,1]]
print(s['e_0'], s['e_exact'])
# exit()
print(b0.shape)
print(b.shape)
print(bref.shape)
L = s['L']
for i in range(s['L']):
    print('========================')
    print(i)
    B1 = s1.get_B(i)
    B2 = s3.get_B(L-i-1, 'B')
    print(B1)
    print(B2)
exit()
for i in range(s['L']):
    print(i)
    print('========================')
    x1 = s3._S[i]
    x2 = s1._S[i]
    print(x1)
    print(x2)
    for j in range(b0.shape[1]):
        print('---')
        for k in range(b0.shape[2]):
            x1 = b[i][j][k]
            x2 = bref[i][j][k]
            print(f'{x1:50} {x2:50}')
