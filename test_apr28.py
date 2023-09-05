#!/usr/bin/env python3

import math
import numpy
import scipy
from measure_funcs import *
from measure import make_left
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
from tenpy.networks.mps import MPS, MPSEnvironment
from tenpy.linalg.np_conserved import tensordot, inner, norm
from Slicer import SystemSlicer

S = SystemSet("pttestexact1")
S.run()

for x, y in S.slicer.get_sets(exclude=['method']).items():
    print(f'~~{x}~~')
    for s in y:
        if s['method'] == 'dmrg':
            continue
        print('...')
        s.force_run()
        k = {'method':'dmrg'}

        sdmrg = SystemSlicer(y).find(k)[0]
        sdmrg.load_states()
        print(s)
        print(sdmrg)

        gR = s['eigs'][0][1]
        gL = s['eigs_left'][0][1]
        psiL = s.eng.full_to_mps(gL)
        psiR = s.eng.full_to_mps(gR)

        sC = s.copy({'conj':1})
        sC.force_run()
        psiLC = sC.eng.full_to_mps(sC['eigs_left'][0][1])
        psiRC = sC.eng.full_to_mps(sC['eigs'][0][1])

        psiLdmrg = sdmrg.left_system.states[0].psi
        psiRdmrg = sdmrg.states[0].psi
        psiRdmrg = sdmrg.chargeless_states[0]

        psiL.canonical_form()
        psiR.canonical_form()
        psiLdmrg.canonical_form()
        psiRdmrg.canonical_form()

        print(s['e_0'], sdmrg['e_0'])
        E = sdmrg.model.H_MPO.expectation_value(psiRdmrg) / s['L']
        print(E)
        E = sdmrg.model.H_MPO.expectation_value(psiR) / s['L']
        print(E)
        print(abs(psiRdmrg.overlap(psiR)))
        print(abs(psiLdmrg.overlap(psiL)))
        print(abs(psiL.overlap(psiR)))
        print(abs(psiLdmrg.overlap(psiRdmrg)))
        def f1(a, b):
            print(abs(a.overlap(b)))
        f1(psiLC, psiL)
        f1(psiRC, psiR)
        f1(psiLC, psiR)
        f1(psiRC, psiL)
        f1(psiR, psiRdmrg)
        f1(psiL, s.make_left_state(psiRdmrg))
        f1(psiL, s.make_left_state(psiR))
        f1(psiL, psiR)

