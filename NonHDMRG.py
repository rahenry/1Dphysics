#!/usr/bin/env python3

import numpy as np
import scipy

from tenpy.networks.mps import MPS
from tenpy.models.tf_ising import TFIChain
from tenpy.models.spins import SpinModel
from tenpy.algorithms import dmrg
import tenpy.linalg.np_conserved as npc
from tenpy.linalg.sparse import FlatLinearOperator
from tenpy.tools.math import speigs
from tenpy.algorithms import mps_common
from tenpy.linalg import np_conserved as npc
from rah_utility import rel_diff

def diag_arpack4(H,  psi, options={}, move_right=True):
    H_flat, psi_flat = FlatLinearOperator.from_guess_with_pipe(H.matvec, psi, dtype=np.complex128)
    tol = options.get('P_tol', 1.e-8)
    N_min = options.get('N_min', None)
    z = 'SR'
    try:
        Es, Vs = speigs(H_flat, k=1, which=z, v0=psi_flat, tol=tol, ncv=N_min, maxiter=1E10)
    except scipy.sparse.linalg.ArpackNoConvergence:
        # simply try again with larger "k", that often helps
        new_k = min(6, H_flat.shape[1])
        if new_k <= 1:
            raise
        Es, Vs = speigs(H_flat, k=new_k, which=z, v0=psi_flat, tol=tol, ncv=N_min)
    v = Vs[:,0]
    v1 = H_flat.matvec(v)
    u = np.inner(v1.conj(), v)
    error = rel_diff(u, Es[0])
    # if error > 1E-10:
    #     print(error)

    psi0 = H_flat.flat_to_npc(Vs[:, 0]).split_legs(0)
    psi0.itranspose(psi.get_leg_labels())

    return Es[0], psi0

def diag_arpack3(m, H_flat, psi, psi_flat, options={}, move_right=True):
    # if not move_right:
    #     H = H.adjoint()
    # H = H.adjoint()
    tol = options.get('P_tol', 1.e-8)
    N_min = options.get('N_min', None)
    z = 'SR'
    try:
        Es, Vs = speigs(m, k=1, which=z, v0=psi_flat, tol=tol, ncv=N_min, maxiter=1E10)
    except scipy.sparse.linalg.ArpackNoConvergence:
        # simply try again with larger "k", that often helps
        new_k = min(6, m.shape[1])
        if new_k <= 1:
            raise
        Es, Vs = speigs(m, k=new_k, which=z, v0=psi_flat, tol=tol, ncv=N_min)
    psi0 = H_flat.flat_to_npc(Vs[:, 0]).split_legs(0)
    psi0.itranspose(psi.get_leg_labels())

    return Es[0], psi0

def diag_arpack(H,  psi, options={}, move_right=True):
    # if not move_right:
    #     H = H.adjoint()
    # H = H.adjoint()
    H_flat, psi_flat = FlatLinearOperator.from_guess_with_pipe(H.matvec, psi, dtype=np.complex128)
    tol = options.get('P_tol', 1.e-8)
    N_min = options.get('N_min', None)
    z = 'SR'
    try:
        Es, Vs = speigs(H_flat, k=1, which=z, v0=psi_flat, tol=tol, ncv=N_min, maxiter=1E10)
    except scipy.sparse.linalg.ArpackNoConvergence:
        # simply try again with larger "k", that often helps
        new_k = min(6, H_flat.shape[1])
        if new_k <= 1:
            raise
        Es, Vs = speigs(H_flat, k=new_k, which=z, v0=psi_flat, tol=tol, ncv=N_min)
    psi0 = H_flat.flat_to_npc(Vs[:, 0]).split_legs(0)
    psi0.itranspose(psi.get_leg_labels())

    return Es[0], psi0

def diag_arpack2(H, psi, options={}, move_right=True):
    # if not move_right:
    #     H = H.adjoint()
    H_flat, psi_flat = FlatLinearOperator.from_guess_with_pipe(H.matvec, psi, dtype=np.complex128)
    tol = options.get('P_tol', 1.e-6)
    N_min = options.get('N_min', None)
    k = 10
    best = None
    while (True):
        Es, Vs = speigs(H_flat, k=k, which='SR', v0=psi_flat, tol=tol, ncv=N_min, maxiter=None)
        for i in range(len(Es)):
            e = Es[i]
            if e.imag > 0:
                continue
            if best is None or abs(e.imag) < abs(best.imag):
                best = i
        if best is None:
            k *= 2
        else:
            break
            # if abs(e.imag) < abs(Es[best].imag):
            #     best = i
    # print(best, Es[best].imag)

    psi0 = H_flat.flat_to_npc(Vs[:, best]).split_legs(0)
    psi0.itranspose(psi.get_leg_labels())

    return Es[best], psi0

class NonHermitianDMRG(dmrg.TwoSiteDMRGEngine):
    def diag(self, theta_guess):
        if self.diag_method == 'arpack_eigs':
            E, theta = diag_arpack4(self.eff_H, theta_guess, self.lanczos_params, self.move_right)
        else:
            return super().diag(theta_guess)
        ov_change = 1. - abs(npc.inner(theta_guess, theta, 'labels', do_conj=True))
        return E, theta, -1, ov_change


class SSNonHermitianDMRG(dmrg.SingleSiteDMRGEngine):
    def diag(self, theta_guess):
        if self.diag_method == 'arpack_eigs':
            E, theta = diag_arpack(self.eff_H, theta_guess, self.lanczos_params, self.move_right)
        else:
            return super().diag(theta_guess)
        ov_change = 1. - abs(npc.inner(theta_guess, theta, 'labels', do_conj=True))
        return E, theta, -1, ov_change
