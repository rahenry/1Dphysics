#!/usr/bin/env python3

import warnings

print(1)
from tenpy.linalg.sparse import FlatLinearOperator
from tenpy.tools.events import EventHandler
from tenpy.tools.params import asConfig
from tenpy.algorithms.mps_common import TwoSiteH
from NonHDMRG import diag_arpack3
from ZnModel import ZnModel

class AlgorithmNH:
    def __init__(self, psiL, psiR, modelL, modelR, options, *, resume_data=None):
        self.options = asConfig(options, self.__class__.__name__)
        self.trunc_params = self.options.subconfig('trunc_params')
        self.psiL = psiL
        self.psiR = psiR
        self.modelL = modelL
        self.modelR = modelR
        if resume_data is None:
            resume_data = {}
        self.resume_data = resume_data
        self.checkpoint = EventHandler("algorithm")
        self._resume_psiL = None
        self._resume_psiR = None

    @property
    def verbose(self):
        warnings.warn(
            "verbose is deprecated, we're using logging now! \n"
            "See https://tenpy.readthedocs.io/en/latest/intro/logging.html", FutureWarning, 2)
        return self.options.get('verbose', 1.)

    def run(self):
        raise NotImplementedError("Sublcasses should implement this.")

    def resume_run(self):
        return

print(2)
import numpy as np
import time
import warnings
import copy
import logging
logger = logging.getLogger(__name__)

from tenpy.linalg import np_conserved as npc
from tenpy.algorithms.truncation import svd_theta, TruncationError
from tenpy.networks.mps import MPSEnvironment, MPS
from tenpy.networks.mpo import MPOEnvironment
from tenpy.linalg.sparse import NpcLinearOperator, SumNpcLinearOperator, OrthogonalNpcLinearOperator


class SweepNH(AlgorithmNH):
    def __init__(self, psiL, psiR, modelL, modelR, options, *, orthogonal_to=None, **kwargs):
        if not hasattr(self, "EffectiveH"):
            raise NotImplementedError("Subclass needs to set EffectiveH")
        super().__init__(psiL, psiR, modelL, modelR, options, **kwargs)
        options = self.options

        self.combine = options.get('combine', False)
        self.finite = self.psiL.finite
        self.S_inv_cutoff = 1.e-15
        self.lanczos_params = options.subconfig('lanczos_params')

        self.env = None
        self.ortho_to_envs = []
        self.init_env(modelL, modelR, resume_data=self.resume_data, orthogonal_to=orthogonal_to)
        self.i0 = 0
        self.move_right = True
        self.update_LP_RP = (True, False)
        self.time0 = time.time()

    @property
    def engine_params(self):
        warnings.warn("renamed self.engine_params -> self.options", FutureWarning, stacklevel=2)
        return self.options

    @property
    def _all_envs(self):
        return [self.envL, self.envR, self.envIL, self.envIR] + self.ortho_to_envs

    @property
    def n_optimize(self):
        return self.EffectiveH.length

    def init_env(self, modelL=None, modelR=None, resume_data=None, orthogonal_to=None):
        HL = modelL.H_MPO
        HR = modelR.H_MPO

        c = self.modelL.model_params
        c['model_type'] = 'ident'
        modelI = ZnModel(c)
        IR = modelI.H_MPO
        IL = IR

        # actually initialize the environment
        self.envR = MPOEnvironment(self.psiR, HR, self.psiL)
        self.envL = MPOEnvironment(self.psiL, HL, self.psiR)
        self.envIR = MPOEnvironment(self.psiR, IR, self.psiL)
        self.envIL = MPOEnvironment(self.psiL, IL, self.psiR)

    def sweep(self, optimize=True):
        self.E_trunc_list = []
        self.trunc_err_list = []
        schedule = self.get_sweep_schedule()

        # the actual sweep
        for i0, move_right, update_LP_RP in schedule:
            self.i0 = i0
            self.move_right = move_right
            self.update_LP_RP = update_LP_RP
            logger.debug("in sweep: i0 =%d", i0)
            # --------- the main work --------------
            thetaL, thetaR = self.prepare_update()
            update_data = self.update_local(thetaL, thetaR, optimize=optimize)
            self.update_env(**update_data)
            self.post_update_local(**update_data)

        if optimize:  # count optimization sweeps
            self.sweeps += 1
        return np.max(self.trunc_err_list)

    def get_sweep_schedule(self):
        L = self.psiR.L
        n = self.EffectiveH.length
        i0s = list(range(0, L - n)) + list(range(L - n, 0, -1))
        move_right = [True] * (L - n) + [False] * (L - n)
        update_LP_RP = [[True, False]] * (L - n) + [[False, True]] * (L - n)
        return zip(i0s, move_right, update_LP_RP)


    def prepare_update(self):
        self.make_eff_H()  # self.eff_H represents tensors LP, W0, RP
        # make theta
        thetaL = self.psiL.get_theta(self.i0, n=self.n_optimize, cutoff=self.S_inv_cutoff)
        thetaR = self.psiR.get_theta(self.i0, n=self.n_optimize, cutoff=self.S_inv_cutoff)
        thetaL = self.eff_HL.combine_theta(thetaL)
        thetaR = self.eff_HR.combine_theta(thetaR)
        return thetaL, thetaR

    def make_eff_H(self):
        """Create new instance of `self.EffectiveH` at `self.i0` and set it to `self.eff_H`."""
        self.eff_HL = self.EffectiveH(self.envL, self.i0, self.combine, self.move_right)
        self.eff_HR = self.EffectiveH(self.envR, self.i0, self.combine, self.move_right)
        self.eff_IL = self.EffectiveH(self.envIL, self.i0, self.combine, self.move_right)
        self.eff_IR = self.EffectiveH(self.envIR, self.i0, self.combine, self.move_right)

    def update_local(self, thetaL, thetaR, **kwargs):
        raise NotImplementedError("needs to be overridden by subclass")

    def update_env(self, **update_data):
        i_L, i_R = self._update_env_inds()  # left and right updated sites
        self.envL.del_LP(i_R)
        self.envR.del_LP(i_R)
        self.envL.del_RP(i_L)
        self.envR.del_RP(i_L)
        self.envIL.del_LP(i_R)
        self.envIR.del_LP(i_R)
        self.envIL.del_RP(i_L)
        self.envIR.del_RP(i_L)
        update_LP, update_RP = self.update_LP_RP
        combine = self.combine
        if update_LP:
            self.eff_HL.update_LP(self.envL, i_R, update_data['UL'])  # possibly optimized
            self.eff_HR.update_LP(self.envR, i_R, update_data['UR'])  # possibly optimized
            self.eff_IL.update_LP(self.envIL, i_R, update_data['UL'])  # possibly optimized
            self.eff_IR.update_LP(self.envIR, i_R, update_data['UR'])  # possibly optimized
        if update_RP:
            self.eff_HL.update_RP(self.envL, i_L, update_data['VHL'])  # possibly optimized
            self.eff_HR.update_RP(self.envR, i_L, update_data['VHR'])  # possibly optimized
            self.eff_IL.update_RP(self.envIL, i_L, update_data['VHL'])  # possibly optimized
            self.eff_IR.update_RP(self.envIR, i_L, update_data['VHR'])  # possibly optimized

    def _update_env_inds(self):
        n = self.n_optimize  # = 1 or 2
        move_right = self.move_right
        if n == 2 or move_right:
            i_L = self.i0
            i_R = self.i0 + 1
        else:  # n == 1 and left moving
            i_L = self.i0 - 1
            i_R = self.i0
        return i_L, i_R

    def post_update_local(self, err, **update_data):
        """Algorithm-specific actions to be taken after local update.

        An example would be to collect statistics.
        """
        self.trunc_err_list.append(err.eps)

print(3)
import numpy as np
import time
import warnings
import logging
logger = logging.getLogger(__name__)

from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPSEnvironment
from tenpy.linalg.lanczos import lanczos, lanczos_arpack
from tenpy.algorithms.truncation import truncate, svd_theta
from tenpy.tools.params import asConfig
from tenpy.tools.params import Config
from tenpy.tools.misc import find_subclass
from tenpy.tools.math import matvec_to_array
from tenpy.tools.process import memory_usage
from tenpy.algorithms.mps_common import Sweep, OneSiteH, TwoSiteH


class DMRGEngineNH(SweepNH):
    EffectiveH = None
    DefaultMixer = None

    def __init__(self, psiL, psiR, modelL, modelR, options, **kwargs):
        options = asConfig(options, self.__class__.__name__)
        self.mixer = None
        self.diag_method = options.get('diag_method', 'default')
        super().__init__(psiL, psiR, modelL, modelR, options, **kwargs)

    @property
    def DMRG_params(self):
        warnings.warn("renamed self.DMRG_params -> self.options", FutureWarning, stacklevel=2)
        return self.options

    def run(self):
        self.E0 = 1000
        options = self.options
        start_time = self.time0
        self.shelve = False
        # parameters for lanczos
        p_tol_to_trunc = options.get('P_tol_to_trunc', 0.05)
        if p_tol_to_trunc is not None:
            svd_min = self.trunc_params.silent_get('svd_min', 0.)
            svd_min = 0. if svd_min is None else svd_min
            trunc_cut = self.trunc_params.silent_get('trunc_cut', 0.)
            trunc_cut = 0. if trunc_cut is None else trunc_cut
            p_tol_min = max(1.e-30, svd_min**2 * p_tol_to_trunc, trunc_cut**2 * p_tol_to_trunc)
            p_tol_min = options.get('P_tol_min', p_tol_min)
            p_tol_max = options.get('P_tol_max', 1.e-4)
        e_tol_to_trunc = options.get('E_tol_to_trunc', None)
        if e_tol_to_trunc is not None:
            e_tol_min = options.get('E_tol_min', 5.e-16)
            e_tol_max = options.get('E_tol_max', 1.e-4)

        # parameters for DMRG convergence criteria
        N_sweeps_check = options.get('N_sweeps_check', 1)
        min_sweeps = int(1.5 * N_sweeps_check)
        min_sweeps = options.get('min_sweeps', min_sweeps)
        max_sweeps = options.get('max_sweeps', 1000)
        max_E_err = options.get('max_E_err', 1.e-8)
        max_S_err = options.get('max_S_err', 1.e-5)
        max_seconds = 3600 * options.get('max_hours', 24 * 365)
        if not self.finite:
            update_env = options.get('update_env', N_sweeps_check // 2)
        E_old, S_old = np.nan, np.mean(self.psiR.entanglement_entropy())  # initial dummy values
        E, Delta_E, Delta_S = 1., 1., 1.
        self.diag_method = options['diag_method']

        is_first_sweep = True
        # loop over sweeps
        self.sweeps = 0
        while True:
            loop_start_time = time.time()
            # check convergence criteria
            if self.sweeps >= max_sweeps:
                break
            if (self.sweeps > min_sweeps and -Delta_E < max_E_err * max(abs(E), 1.)
                    and abs(Delta_S) < max_S_err):
                break
            if loop_start_time - start_time > max_seconds:
                self.shelve = True
                logger.warning("DMRG: maximum time limit reached. Shelve simulation.")
                break
            if not is_first_sweep:
                self.checkpoint.emit(self)
            # --------- the main work --------------
            logger.info('Running sweep with optimization')
            for i in range(N_sweeps_check - 1):
                self.sweep(meas_E_trunc=False)
            max_trunc_err = self.sweep(meas_E_trunc=True)
            max_E_trunc = np.max(self.E_trunc_list)
            # --------------------------------------
            # update lancos_params depending on truncation error(s)
            if p_tol_to_trunc is not None and max_trunc_err > p_tol_min:
                P_tol = max(p_tol_min, min(p_tol_max, max_trunc_err * p_tol_to_trunc))
                self.lanczos_params['P_tol'] = P_tol
                logger.debug("set lanczos_params['P_tol'] = %.2e", P_tol)
            if e_tol_to_trunc is not None and max_E_trunc > e_tol_min:
                E_tol = max(e_tol_min, min(e_tol_max, max_E_trunc * e_tol_to_trunc))
                self.lanczos_params['E_tol'] = E_tol
                logger.debug("set lanczos_params['E_tol'] = %.2e", E_tol)
            # update environment
            if not self.finite:
                self.environment_sweeps(update_env)

            # update values for checking the convergence
            try:
                S = self.psiR.entanglement_entropy()
                max_S = max(S)
                S = np.mean(S)
                Delta_S = (S - S_old) / N_sweeps_check
            except ValueError:
                # with a mixer, psi._S can be 2D arrays s.t. entanglement_entropy() fails
                S = np.nan
                max_S = np.nan
                Delta_S = 0.
            S_old = S
            E = self.E0
            Delta_E = (E - E_old) / N_sweeps_check
            E_old = E
            norm_err = np.linalg.norm(self.psiR.norm_test())

            # update statistics
            # self.sweep_stats['sweep'].append(self.sweeps)
            # self.sweep_stats['N_updates'].append(len(self.update_stats['i0']))
            # self.sweep_stats['E'].append(E)
            # self.sweep_stats['S'].append(S)
            # self.sweep_stats['time'].append(time.time() - start_time)
            # self.sweep_stats['max_trunc_err'].append(max_trunc_err)
            # self.sweep_stats['max_E_trunc'].append(max_E_trunc)
            # self.sweep_stats['max_chi'].append(np.max(self.psi.chi))
            # self.sweep_stats['norm_err'].append(norm_err)

            # status update
            logger.info(
                "checkpoint after sweep %(sweeps)d\n"
                "energy=%(E).16f, max S=%(S).16f, norm_err=%(norm_err).1e\n"
                "Current memory usage %(mem).1fMB, wall time: %(wall_time).1fs\n"
                "Delta E = %(dE).4e, Delta S = %(dS).4e (per sweep)\n"
                "max trunc_err = %(trunc_err).4e, max E_trunc = %(E_trunc).4e\n"
                "chi: %(chi)s\n"
                "%(sep)s", {
                    'sweeps': self.sweeps,
                    'E': E,
                    'S': max_S,
                    # 'age': self.update_stats['age'][-1],
                    'norm_err': norm_err,
                    'mem': memory_usage(),
                    'wall_time': time.time() - loop_start_time,
                    'dE': Delta_E,
                    'dS': Delta_S,
                    'trunc_err': max_trunc_err,
                    'E_trunc': max_E_trunc,
                    'chi': self.psiR.chi if self.psiR.L < 40 else max(self.psiR.chi),
                    'sep': "=" * 80,
                })
            is_first_sweep = False


        self._canonicalize(True)
        logger.info("DMRG finished after %d sweeps, max chi=%d", self.sweeps, max(self.psiR.chi))
        return E, self.psiL, self.psiR

    def _canonicalize(self, warn=False):
        # update environment until norm_tol is reached
        norm_err = np.linalg.norm(self.psiR.norm_test())
        norm_tol = self.options.get('norm_tol', 1.e-5)
        if norm_tol is None or norm_err < norm_tol:
            return
        if warn:
            logger.warning(
                "final DMRG state not in canonical form up to "
                "norm_tol=%.2e: norm_err=%.2e", norm_tol, norm_err)
        self.psiL.canonical_form()
        self.psiR.canonical_form()

    def reset_stats(self, resume_data=None):
        """Reset the statistics, useful if you want to start a new sweep run.

        .. cfg:configoptions :: DMRGEngine

            chi_list : dict | None
                A dictionary to gradually increase the `chi_max` parameter of
                `trunc_params`. The key defines starting from which sweep
                `chi_max` is set to the value, e.g. ``{0: 50, 20: 100}`` uses
                ``chi_max=50`` for the first 20 sweeps and ``chi_max=100``
                afterwards. Overwrites `trunc_params['chi_list']``.
                By default (``None``) this feature is disabled.
            sweep_0 : int
                The number of sweeps already performed. (Useful for re-start).
        """
        super().reset_stats(resume_data)
        self.update_stats = {
            'i0': [],
            'age': [],
            'E_total': [],
            'N_lanczos': [],
            'time': [],
            'err': [],
            'E_trunc': [],
            'ov_change': []
        }
        self.sweep_stats = {
            'sweep': [],
            'N_updates': [],
            'E': [],
            'S': [],
            'time': [],
            'max_trunc_err': [],
            'max_E_trunc': [],
            'max_chi': [],
            'norm_err': []
        }

    def sweep(self, optimize=True, meas_E_trunc=False):
        self._meas_E_trunc = meas_E_trunc
        res = super().sweep(optimize)
        return res

    def update_local(self, thetaL, thetaR, optimize=True):
        i0 = self.i0
        n_opt = self.n_optimize
        age = self.envL.get_LP_age(i0) + n_opt + self.envL.get_RP_age(i0 + n_opt - 1)
        if optimize:
            E0, thetaL, thetaR, N, ov_change = self.diag(thetaL, thetaR)
        else:
            E0, N, ov_change = None, 0, 0.

        thetaR = self.prepare_svd(thetaR)
        UR, SR, VHR, errR = self.mixed_svd(thetaR, self.envR)
        self.set_BR(UR, SR, VHR)


        thetaL = self.prepare_svd(thetaL)
        i0 = self.i0
        tp = self.trunc_params
        qtotal_i0 = self.envL.bra.get_B(i0, form=None).qtotal
        UL, SL, VHL, errL, _ = svd_theta(thetaL,
                                         tp,
                                        qtotal_LR=[qtotal_i0, None],
                                        inner_labels=['vR', 'vL'])
        self.set_BL(UL, SL, VHL)
        print(len(SL), len(SR))
        update_data = {
            'E0': E0,
            'err': errR,
            'N': N,
            'age': age,
            'UR': UR,
            'VHR': VHR,
            'UL': UL,
            'VHL': VHL,
            'ov_change': ov_change
        }
        return update_data

    def post_update_local(self, E0, age, N, ov_change, err, **update_data):
        self.E0 = E0
        i0 = self.i0
        E_trunc = None
        if self._meas_E_trunc or E0 is None:
            E_trunc = self.envR.full_contraction(i0).real  # uses updated LP/RP (if calculated)
            if E0 is None:
                E0 = E_trunc
            E_trunc = E_trunc - E0

        # collect statistics
        # self.update_stats['i0'].append(i0)
        # self.update_stats['age'].append(age)
        # self.update_stats['E_total'].append(E0)
        # self.update_stats['E_trunc'].append(E_trunc)
        # self.update_stats['N_lanczos'].append(N)
        # self.update_stats['ov_change'].append(ov_change)
        # self.update_stats['err'].append(err)
        # self.update_stats['time'].append(time.time() - self.time0)
        self.trunc_err_list.append(err.eps)
        self.E_trunc_list.append(E_trunc)

    def diag(self, theta_guessL, theta_guessR):
        # print(type(self.eff_HL))
        # print(type(self.envL))
        # e = MPSEnvironment(self.psiL, self.psiR)
        # print(e)
        # I, z = FlatLinearOperator.from_guess_with_pipe(e.matvec, theta_guessL)
        IL, z = FlatLinearOperator.from_guess_with_pipe(self.eff_IL.matvec, theta_guessL)
        IR, z = FlatLinearOperator.from_guess_with_pipe(self.eff_IR.matvec, theta_guessR)

        N = -1  # (unknown)

        H_flatL, psi_flatL = FlatLinearOperator.from_guess_with_pipe(self.eff_HL.matvec, theta_guessL, dtype=np.complex128)
        H_flatR, psi_flatR = FlatLinearOperator.from_guess_with_pipe(self.eff_HR.matvec, theta_guessR, dtype=np.complex128)

        HL = matvec_to_array(H_flatL)
        HR = matvec_to_array(H_flatR)
        IL = matvec_to_array(IL)
        IR = matvec_to_array(IR)

        def ct(x):
            return np.conjugate(np.transpose(x))
        a = numpy.matmul(IL, psi_flatL)
        b = numpy.matmul(ct(psi_flatL), HR)
        m1 = numpy.outer(numpy.matmul(IL, psi_flatL), numpy.matmul(ct(psi_flatL), HR))
        m2 = numpy.outer(numpy.matmul(IR, psi_flatR), numpy.matmul(ct(psi_flatR), HL))
        print(m1.shape)
        print(HR.shape)
        print(m1)
        print(HR)
        exit()
        # m1 = HR
        # m2 = HL
        # m1 = numpy.matmul(mL, mR)
        # m2 = numpy.matmul(mR, mL)
        # z = inner(theta_guessL, theta_guessR, do_conj=True)
        # print(f'z={z}')
        # print(type(H_flatL))
        # print(H_flatL)
        # print(H_flatL.shape)
        # exit()
        # EL, thetaL = diag_arpack3(H_flatL, theta_guessL, psi_flatL, self.lanczos_params, self.move_right)
        # ER, thetaR = diag_arpack3(H_flatR, theta_guessR, psi_flatR, self.lanczos_params, self.move_right)
        EL, thetaL = diag_arpack3(m2, H_flatL, theta_guessL, psi_flatL, self.lanczos_params, self.move_right)
        ER, thetaR = diag_arpack3(m1, H_flatR, theta_guessR, psi_flatR, self.lanczos_params, self.move_right)
        print(ER, EL)
        # ER, thetaR = diag_arpack(self.eff_HR, theta_guessR, self.lanczos_params, self.move_right)
        # EL = EL / o2
        # ER = ER / o1

        ov_changeL = 1. - abs(npc.inner(theta_guessL, thetaL, 'labels', do_conj=True))
        ov_changeR = 1. - abs(npc.inner(theta_guessR, thetaR, 'labels', do_conj=True))
        o = inner(thetaL, thetaR, do_conj=True)
        o = numpy.sqrt(o)
        o = 1
        # thetaL /= o
        # thetaR /= o
        o = inner(thetaL, thetaR, do_conj=True)
        print(o, abs(o))
        return ER, thetaL, thetaR, -1, ov_changeR


class TwoSiteDMRGEngineNH(DMRGEngineNH):
    EffectiveH = TwoSiteH

    def prepare_svd(self, theta):
        """Transform theta into matrix for svd."""
        if self.combine:
            return theta  # Theta is already combined.
        else:
            return theta.combine_legs([['vL', 'p0'], ['p1', 'vR']],
                                      new_axes=[0, 1],
                                      qconj=[+1, -1])

    def mixed_svd(self, theta, env):
        i0 = self.i0
        qtotal_i0 = env.bra.get_B(i0, form=None).qtotal
        U, S, VH, err, _ = svd_theta(theta,
                                        self.trunc_params,
                                        qtotal_LR=[qtotal_i0, None],
                                        inner_labels=['vR', 'vL'])
        return U, S, VH, err

    def set_BL(self, U, S, VH):
        B0 = U.split_legs(['(vL.p0)']).replace_label('p0', 'p')
        B1 = VH.split_legs(['(p1.vR)']).replace_label('p1', 'p')
        i0 = self.i0
        self.psiL.set_B(i0, B0, form='A')  # left-canonical
        self.psiL.set_B(i0 + 1, B1, form='B')  # right-canonical
        self.psiL.set_SR(i0, S)


    def set_BR(self, U, S, VH):
        B0 = U.split_legs(['(vL.p0)']).replace_label('p0', 'p')
        B1 = VH.split_legs(['(p1.vR)']).replace_label('p1', 'p')
        i0 = self.i0
        self.psiR.set_B(i0, B0, form='A')  # left-canonical
        self.psiR.set_B(i0 + 1, B1, form='B')  # right-canonical
        self.psiR.set_SR(i0, S)


import numpy
import warnings
import math
import scipy
from measure_funcs import *
from SystemSet import SystemSet
from System import System
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

S = SystemSet("nonHdmrgtest")

for s in S.systems:
    L = s['L']
    s.init_model()
    s.init_states()
    c = Config(dict(s.config.options), "solver_config")
    psiL = s.states[0].psi
    psiR = psiL.copy()

    # eng = TwoSiteDMRGEngineNH(psiR, psiR, s.model, s.model, c)
    # E, psiL, psiR = eng.run()
    # print(E/L, s['e_exact'])
    # exit()

    d = dict(s.config_base)
    d['conj'] = 1
    c['trunc_params'] = {
        'chi_max': 500,
        'svd_min': 0,
        'trunc_cut':0,
        }
    sleft = System(d)
    eng = TwoSiteDMRGEngineNH(psiL, psiR, s.model, sleft.model, c)
    E, psiL, psiR = eng.run()
    print(psiL.norm, psiR.norm)
    print(234234)
    o = psiL.overlap(psiR)
    E = s.model.H_MPO.expectation_value(psiR)
    print(o)
    print(E, s['e_exact'])
    print(abs(E/s['L']))
    print(abs(E/o)/s['L'])
