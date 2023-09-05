#!/usr/bin/env python3

import numpy as np
import warnings
import random
from rah_utility import rel_diff
from functools import reduce
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

def correlation_function(
    MPS1,
    MPS2,
    ops1,
    ops2,
    sites1=None,
    sites2=None,
    opstr=None,
    str_on_first=True,
    hermitian=False,
    autoJW=True,
        do_conj=1,
):
    MPS1 = MPS1.copy();
    MPS2 = MPS2.copy();
    if opstr is not None:
        autoJW = False
    ops1, ops2, sites1, sites2, opstr = MPS1._correlation_function_args(
        ops1, ops2, sites1, sites2, opstr
    )
    if (len(sites1) > 2 * len(sites2) and min(sites2) > max(sites1) - len(sites2)) or (
        len(sites2) > 2 * len(sites1) and min(sites1) > max(sites2) - len(sites1)
    ):
        warnings.warn(
            "Inefficent evaluation of MPS.correlation_function(), "
            "it's probably faster to use MPS.term_correlation_function_left()",
            stracklevel=2,
        )
    if autoJW and not all([isinstance(op1, str) for op1 in ops1]):
        warnings.warn(
            "Non-string operator: can't auto-determine Jordan-Wigner!", stacklevel=2
        )
        autoJW = False
    if autoJW:
        need_JW = []
        for i in sites1:
            need_JW.append(MPS1.sites[i % MPS1.L].op_needs_JW(ops1[i % len(ops1)]))
        for j in sites2:
            need_JW.append(MPS2.sites[j % MPS2.L].op_needs_JW(ops1[j % len(ops1)]))
        if any(need_JW):
            if not all(need_JW):
                raise ValueError("Some, but not any operators need 'JW' string!")
            if not str_on_first:
                raise ValueError("Need Jordan Wigner string, but `str_on_first`=False`")
            opstr = ["JW"]
    if hermitian and np.any(sites1 != sites2):
        warnings.warn(
            "MPS correlation function can't use the hermitian flag", stacklevel=2
        )
        hermitian = False
    C = np.empty((len(sites1), len(sites2)), dtype=complex)
    for x, i in enumerate(sites1):
        # j > i
        j_gtr = sites2[sites2 > i]
        if len(j_gtr) > 0:
            C_gtr = corr_up_diag(MPS1, MPS2, ops1, ops2, i, j_gtr, opstr, str_on_first, True, do_conj=do_conj)
            C[x, (sites2 > i)] = C_gtr
            if hermitian:
                C[x + 1 :, x] = np.conj(C_gtr)
        # j == i
        j_eq = sites2[sites2 == i]
        if len(j_eq) > 0:
            # on-site correlation function
            op12 = npc.tensordot(
                MPS1.get_op(ops1, i), MPS1.get_op(ops2, i), axes=["p*", "p"]
            )
            C[x, (sites2 == i)] = expectation_value(MPS1, MPS2, op12, i, [["p"], ["p*"]], do_conj=do_conj)
    if not hermitian:
        #  j < i
        for y, j in enumerate(sites2):
            i_gtr = sites1[sites1 > j]
            if len(i_gtr) > 0:
                C[(sites1 > j), y] = corr_up_diag(
                    MPS1, MPS2, ops2, ops1, j, i_gtr, opstr, str_on_first, False, do_conj=do_conj
                )
                # exchange ops1 and ops2 : they commute on different sites,
                # but we apply opstr after op1 (using the last argument = False)
    return np.real_if_close(C)

def corr_up_diag(MPS1, MPS2, ops1, ops2, i, j_gtr, opstr, str_on_first, apply_opstr_first, do_conj=1):
    """correlation function above the diagonal: for fixed i and all j in j_gtr, j > i."""
    op1 = MPS1.get_op(ops1, i)
    opstr1 = MPS1.get_op(opstr, i)
    if opstr1 is not None and str_on_first:
        axes = ['p*', 'p'] if apply_opstr_first else ['p', 'p*']
        op1 = npc.tensordot(op1, opstr1, axes=axes)
    theta1 = MPS1.get_theta(i, n=1)
    theta2 = MPS2.get_theta(i, n=1).conj(complex_conj=do_conj)
    C = npc.tensordot(op1, theta1, axes=['p*', 'p0'])
    C = npc.tensordot(theta2, C, axes=[['p0*', 'vL*'], ['p', 'vL']])
    # C has legs 'vR*', 'vR'
    js = list(j_gtr[::-1])  # stack of j, sorted *descending*
    res = []
    for r in range(i + 1, js[0] + 1):  # js[0] is the maximum
        B1 = MPS1.get_B(r, form='B')
        B2 = MPS1.get_B(r, form='B').conj(complex_conj=do_conj)
        C = npc.tensordot(C, B1, axes=['vR', 'vL'])
        if r == js[-1]:
            Cij = npc.tensordot(MPS1.get_op(ops2, r), C, axes=['p*', 'p'])
            Cij = npc.inner(B2, Cij, axes=[['vL*', 'p*', 'vR*'], ['vR*', 'p', 'vR']])
            res.append(Cij)
            js.pop()
        if len(js) > 0:
            op = MPS1.get_op(opstr, r)
            if op is not None:
                C = npc.tensordot(op, C, axes=['p*', 'p'])
            C = npc.tensordot(B2, C, axes=[['vL*', 'p*'], ['vR*', 'p']])
    return res

def expectation_value(MPS1, MPS2, ops, sites=None, axes=None, do_conj=1):
        ops, sites, n, (op_ax_p, op_ax_pstar) = MPS1._expectation_value_args(ops, sites, axes)
        ax_p = ['p' + str(k) for k in range(n)]
        ax_pstar = ['p' + str(k) + '*' for k in range(n)]
        E = []
        for i in sites:
            op = MPS1.get_op(ops, i)
            op = op.replace_labels(op_ax_p + op_ax_pstar, ax_p + ax_pstar)
            theta1 = MPS1.get_theta(i, n)
            theta2 = MPS2.get_theta(i, n).conj(complex_conj=do_conj)
            C = npc.tensordot(op, theta1, axes=[ax_pstar, ax_p])  # C has same labels as theta
            E.append(npc.inner(theta2, C, axes='labels', do_conj=False))
        return np.real_if_close(np.array(E))

def corr_ops_LP(MPS1, MPS2, operators, i0):
    """Contract the left part of a correlation function.

    Same as :meth:`expectation_value_multi_sites`, but with the right-most legs left open,
    with labels ``'vR*', 'vR'``.
    """
    op = operators[0]
    if (isinstance(op, str)):
        op = MPS2.sites[MPS2._to_valid_index(i0)].get_op(op)
    theta2 = MPS2.get_B(i0, 'Th')
    C = npc.tensordot(op, theta2, axes=['p*', 'p'])
    axes = [['vL*'] + MPS2._get_p_label('*'), ['vL'] + MPS2._p_label]
    theta1 = MPS1.get_B(i0, 'Th')
    C = npc.tensordot(theta1.conj(), C, axes=axes)
    axes[1][0] = 'vR*'
    for j in range(1, len(operators)):
        op = operators[j]  # the operator
        is_str = isinstance(op, str)
        i = i0 + j  # the site it acts on
        B = MPS2.get_B(i, form='B')
        C = npc.tensordot(C, B, axes=['vR', 'vL'])
        if not (is_str and op == 'Id'):
            if is_str:
                op = MPS2.sites[MPS2._to_valid_index(i)].get_op(op)
            C = npc.tensordot(op, C, axes=['p*', 'p'])
        C = npc.tensordot(MPS1.get_B(i, form='B').conj(), C, axes=axes)
    return C

def expectation_value_multi_sites(MPS1, MPS2, operators, i0):
    C = corr_ops_LP(MPS1, MPS2, operators, i0)
    exp_val = npc.trace(C, 'vR*', 'vR')
    return np.real_if_close(exp_val)

class MPSEnvironment2:
    def __init__(self, bra, ket, cache=None, do_conj=True, **init_env_data):
        self.do_conj = do_conj
        if ket is None:
            ket = bra
        if ket is not bra:
            ket._gauge_compatible_vL_vR(bra)  # ensure matching charges
        self.bra = bra
        self.ket = ket
        self.dtype = np.find_common_type([bra.dtype, ket.dtype], [])
        self.L = L = lcm(bra.L, ket.L)
        if hasattr(self, 'H'):
            self.L = L = lcm(self.H.L, L)
        self._finite = self.ket.finite  # just for _to_valid_index
        self._LP_keys = ['LP_{0:d}'.format(i) for i in range(L)]
        self._RP_keys = ['RP_{0:d}'.format(i) for i in range(L)]
        self._LP_age = [None] * L
        self._RP_age = [None] * L
        if cache is None:
            cache = DictCache.trivial()
        self.cache = cache
        if not self.cache.long_term_storage.trivial and L < 8:
            warnings.warn("non-trivial cache for short-length environment: "
                          "Much overhead for a little RAM saving. Necessary?")
        self.init_first_LP_last_RP(**init_env_data)
        self.test_sanity()

    def init_first_LP_last_RP(self,
                              init_LP=None,
                              init_RP=None,
                              age_LP=0,
                              age_RP=0,
                              start_env_sites=0):
        """(Re)initialize first LP and last RP from the given data.

        Parameters
        ----------
        init_LP : ``None`` | :class:`~tenpy.linalg.np_conserved.Array`
            Initial very left part ``LP``. If ``None``, build one with :meth`init_LP`.
        init_RP : ``None`` | :class:`~tenpy.linalg.np_conserved.Array`
            Initial very right part ``RP``. If ``None``, build one with :meth:`init_RP`.
        age_LP : int
            The number of physical sites involved into the contraction of `init_LP`.
        age_RP : int
            The number of physical sites involved into the contraction of `init_RP`.
        start_env_sites : int
            If `init_LP` and `init_RP` are not specified, contract each `start_env_sites` for them.
        """
        vL_ket, vR_ket = self.ket._outer_virtual_legs()
        vL_bra, vR_bra = self.bra._outer_virtual_legs()
        ket_U, ket_V = self.ket.segment_boundaries
        bra_U, bra_V = self.bra.segment_boundaries
        if init_LP is not None:
            compatible = (init_LP.get_leg('vR') == vL_ket.conj()
                          and init_LP.get_leg('vR*') == vL_bra)
            if not compatible:
                logger.warning("dropping `init_LP` with incompatible MPS legs")
                init_LP = None
        if init_RP is not None:
            compatible = (init_RP.get_leg('vL') == vR_ket.conj()
                          and init_RP.get_leg('vL*') == vR_bra)
            if not compatible:
                logger.warning("dropping `init_RP` with incompatible MPS legs")
                init_RP = None
        if init_LP is None:
            init_LP = self.init_LP(0, start_env_sites)
            age_LP = start_env_sites
        else:
            if ket_U is not None:
                init_LP = npc.tensordot(init_LP, ket_U, axes=['vR', 'vL'])
            if bra_U is not None:
                init_LP = npc.tensordot(bra_U.conj(complex_conjugate=self.do_conj), init_LP, axes=['vL*', 'vR*'])
        if init_RP is None:
            init_RP = self.init_RP(self.L - 1, start_env_sites)
            age_RP = start_env_sites
        else:
            if ket_V is not None:
                init_RP = npc.tensordot(ket_V, init_RP, axes=['vR', 'vL'])
            if bra_V is not None:
                init_RP = npc.tensordot(init_RP, bra_V.conj(complex_conjugate=self.do_conj), axes=['vL*', 'vR*'])
        self.set_LP(0, init_LP, age=age_LP)
        self.set_RP(self.L - 1, init_RP, age=age_RP)

    def test_sanity(self):
        """Sanity check, raises ValueErrors, if something is wrong."""
        assert (self.bra.finite == self.ket.finite == self._finite)
        assert any(key in self.cache for key in self._LP_keys)
        assert any(key in self.cache for key in self._RP_keys)

    def init_LP(self, i, start_env_sites=0):
        """Build initial left part ``LP``.

        If `bra` and `ket` are the same and in left canonical form, this is the environment
        you get contracting he overlaps from the left infinity up to bond left of site `i`.

        For segment MPS, the :attr:`~tenpy.networks.mps.MPS.segment_boundaries` are read out
        (if set).

        Parameters
        ----------
        i : int
            Build ``LP`` left of site `i`.
        start_env_sites : int
            How many sites to contract to converge the `init_LP`; the initial `age_LP`.

        Returns
        -------
        init_LP : :class:`~tenpy.linalg.np_conserved.Array`
            Identity contractible with the `vL` leg of ``ket.get_B(i)``, labels ``'vR*', 'vR'``.
        """
        if self.ket.bc == "segment" and self.bra is not self.ket:
            U_bra, V_bra = self.bra.segment_boundaries
            U_ket, V_ket = self.ket.segment_boundaries
            if U_bra is not None or U_ket is not None:
                if U_bra is not None and U_ket is not None:
                    init_LP = npc.tensordot(U_bra.conj(complex_conjugate=self.do_conj), U_ket, axes=['vL*', 'vL'])
                elif U_bra is not None:
                    init_LP = U_bra.conj(complex_conjugate=self.do_conj).ireplace_label('vL*', 'vR')
                else:
                    init_LP = U_ket.replace_label('vL', 'vR*')
                return init_LP
        leg_ket = self.ket.get_B(i - start_env_sites, None).get_leg('vL')
        leg_bra = self.bra.get_B(i - start_env_sites, None).get_leg('vL')
        leg_ket.test_equal(leg_bra)
        init_LP = npc.diag(1., leg_ket, dtype=self.dtype, labels=['vR*', 'vR'])
        for j in range(i - start_env_sites, i):
            init_LP = self._contract_LP(j, init_LP)
        return init_LP

    def init_RP(self, i, start_env_sites=0):
        """Build initial right part ``RP`` for an MPS/MPOEnvironment.

        If `bra` and `ket` are the same and in right canonical form, this is the environment
        you get contracting from the right infinity up to bond right of site `i`.

        For segment MPS, the :attr:`~tenpy.networks.mps.MPS.segment_boundaries` are read out
        (if set).

        Parameters
        ----------
        i : int
            Build ``RP`` right of site `i`.
        start_env_sites : int
            How many sites to contract to converge the `init_RP`; the initial `age_RP`.

        Returns
        -------
        init_RP : :class:`~tenpy.linalg.np_conserved.Array`
            Identity contractible with the `vR` leg of ``ket.get_B(i)``, labels ``'vL*', 'vL'``.
        """
        if self.ket.bc == "segment" and self.bra is not self.ket:
            U_bra, V_bra = self.bra.segment_boundaries
            U_ket, V_ket = self.ket.segment_boundaries
            if V_bra is not None or V_ket is not None:
                if V_bra is not None and V_ket is not None:
                    init_RP = npc.tensordot(V_bra.conj(complex_conjugate=self.do_conj), V_ket, axes=['vR*', 'vR'])
                elif V_bra is not None:
                    init_RP = V_bra.conj(complex_conjugate=self.do_conj).ireplace_label('vR*', 'vL')
                else:
                    init_RP = V_ket.replace_label('vR', 'vL*')
                return init_RP
        leg_ket = self.ket.get_B(i + start_env_sites, None).get_leg('vR')
        leg_bra = self.bra.get_B(i + start_env_sites, None).get_leg('vR')
        leg_ket.test_equal(leg_bra)
        init_RP = npc.diag(1., leg_ket, dtype=self.dtype, labels=['vL*', 'vL'])
        for j in range(i + start_env_sites, i, -1):
            init_RP = self._contract_RP(j, init_RP)
        return init_RP

    def get_LP(self, i, store=True):
        """Calculate LP at given site from nearest available one.

        The returned ``LP_i`` corresponds to the following contraction,
        where the M's and the N's are in the 'A' form::

            |     .-------M[0]--- ... --M[i-1]--->-   'vR'
            |     |       |             |
            |     LP[0]   |             |
            |     |       |             |
            |     .-------N[0]*-- ... --N[i-1]*--<-   'vR*'

        Parameters
        ----------
        i : int
            The returned `LP` will contain the contraction *strictly* left of site `i`.
        store : bool
            Wheter to store the calculated `LP` in `self` (``True``) or discard them (``False``).

        Returns
        -------
        LP_i : :class:`~tenpy.linalg.np_conserved.Array`
            Contraction of everything left of site `i`,
            with labels ``'vR*', 'vR'`` for `bra`, `ket`.
        """
        # find nearest available LP to the left.
        for i0 in range(i, i - self.L, -1):
            key = self._LP_keys[self._to_valid_index(i0)]
            LP = self.cache.get(key, None)
            if LP is not None:
                break
            # (for finite, LP[0] should always be set, so we should abort at latest with i0=0)
        else:  # no break called
            raise ValueError("No left part in the system???")
        age = self.get_LP_age(i0)
        for j in range(i0, i):
            LP = self._contract_LP(j, LP)
            age = age + 1
            if store:
                self.set_LP(j + 1, LP, age=age)
        return LP

    def get_RP(self, i, store=True):
        """Calculate RP at given site from nearest available one.

        The returned ``RP_i`` corresponds to the following contraction,
        where the M's and the N's are in the 'B' form::

            |     'vL'  ->---M[i+1]-- ... --M[L-1]----.
            |                |              |         |
            |                |              |         RP[-1]
            |                |              |         |
            |     'vL*' -<---N[i+1]*- ... --N[L-1]*---.


        Parameters
        ----------
        i : int
            The returned `RP` will contain the contraction *strictly* right of site `i`.
        store : bool
            Wheter to store the calculated `RP` in `self` (``True``) or discard them (``False``).

        Returns
        -------
        RP_i : :class:`~tenpy.linalg.np_conserved.Array`
            Contraction of everything left of site `i`,
            with labels ``'vL', 'vL*'`` for `ket`, `bra`.
        """
        # find nearest available RP to the right.
        for i0 in range(i, i + self.L):
            key = self._RP_keys[self._to_valid_index(i0)]
            RP = self.cache.get(key, None)
            if RP is not None:
                break
            # (for finite, RP[-1] should always be set, so we should abort at latest with i0=L-1)
        else:  # no break called
            raise ValueError("No right part in the system???")
        age = self.get_RP_age(i0)
        for j in range(i0, i, -1):
            RP = self._contract_RP(j, RP)
            age = age + 1
            if store:
                self.set_RP(j - 1, RP, age=age)
        return RP

    def get_LP_age(self, i):
        """Return number of physical sites in the contractions of get_LP(i).

        Might be ``None``.
        """
        return self._LP_age[self._to_valid_index(i)]

    def get_RP_age(self, i):
        """Return number of physical sites in the contractions of get_RP(i).

        Might be ``None``.
        """
        return self._RP_age[self._to_valid_index(i)]

    def set_LP(self, i, LP, age):
        """Store part to the left of site `i`."""
        i = self._to_valid_index(i)
        self.cache[self._LP_keys[i]] = LP
        self._LP_age[i] = age

    def set_RP(self, i, RP, age):
        """Store part to the right of site `i`."""
        i = self._to_valid_index(i)
        self.cache[self._RP_keys[i]] = RP
        self._RP_age[i] = age

    def del_LP(self, i):
        """Delete stored part strictly to the left of site `i`."""
        i = self._to_valid_index(i)
        del self.cache[self._LP_keys[i]]
        self._LP_age[i] = None

    def del_RP(self, i):
        """Delete stored part scrictly to the right of site `i`."""
        i = self._to_valid_index(i)
        del self.cache[self._RP_keys[i]]
        self._RP_age[i] = None

    def clear(self):
        """Delete all partial contractions except the left-most `LP` and right-most `RP`."""
        for key in self._LP_keys[1:] + self._RP_keys[:-1]:
            if key in self.cache:
                del self.cache[key]
        self._LP_age[1:] = [None] * (self.L - 1)
        self._RP_age[:-1] = [None] * (self.L - 1)

    def has_LP(self, i):
        """Return True if `LP` left of site `i` is stored."""
        return self._LP_keys[self._to_valid_index(i)] in self.cache

    def has_RP(self, i):
        """Return True if `RP` right of site `i` is stored."""
        return self._RP_keys[self._to_valid_index(i)] in self.cache

    def cache_optimize(self, short_term_LP=[], short_term_RP=[], preload_LP=None, preload_RP=None):
        """Update `short_term_keys` for the cache and possibly preload tensors.

        Parameters
        ----------
        short_term_LP, short_term_RP : list of int
            `i` indices for :meth:`get_LP` and :meth:`get_RP`, respectively, for which a repeated
            look-up could happen, i.e., for which tensors should be kept in RAM until the next
            call to this function.
        preload_LP, preload_RP : int | None
            If not None, preload the tensors for the corrsponding :meth:`get_LP` and :meth:`get_RP`
            call, respectively, from disk.
        """
        LP_keys = self._LP_keys
        RP_keys = self._RP_keys
        preload = []
        if preload_LP is not None:
            preload.append(LP_keys[self._to_valid_index(preload_LP)])
        if preload_RP is not None:
            preload.append(RP_keys[self._to_valid_index(preload_RP)])
        self.cache.set_short_term_keys(*(LP_keys[self._to_valid_index(i)] for i in short_term_LP),
                                       *(RP_keys[self._to_valid_index(i)] for i in short_term_RP),
                                       *preload)
        self.cache.preload(*preload)

    def get_initialization_data(self, first=0, last=None):
        """Return data for (re-)initialization of the environment.

        Parameters
        ----------
        first, last : int
            The first and last site, to the left and right of which we should return the
            environments.  Defaults to 0 and :attr:`L` - 1.

        Returns
        -------
        init_env_data : dict
            A dictionary with the following entries.

            init_LP, init_RP : :class:`~tenpy.linalg.np_conserved.Array`
                `LP` on the left of site `first` and `RP` on the right of site `last`, which can be
                used as `init_LP` and `init_RP` for the initialization of a new environment.
            age_LP, age_RP : int
                The number of physical sites involved into the contraction yielding `init_LP` and
                `init_RP`, respectively.
        """
        L = self.L
        if last is None:
            last = self.L - 1
        data = {'init_LP': self.get_LP(first, True), 'init_RP': self.get_RP(last, True)}
        data['age_LP'] = self.get_LP_age(first)
        data['age_RP'] = self.get_RP_age(last)
        return data

    def full_contraction(self, i0):
        """Calculate the overlap by a full contraction of the network.

        The full contraction of the environments gives the overlap ``<bra|ket>``,
        taking into account :attr:`MPS.norm` of both `bra` and `ket`.
        For this purpose, this function contracts
        ``get_LP(i0+1, store=False)`` and ``get_RP(i0, store=False)`` with appropriate singular
        values in between.

        Parameters
        ----------
        i0 : int
            Site index.
        """
        if self.ket.finite and i0 + 1 == self.L:
            # special case to handle `_to_valid_index` correctly:
            # get_LP(L) is not valid for finite b.c, so we use need to calculate it explicitly.
            LP = self.get_LP(i0, store=False)
            LP = self._contract_LP(i0, LP)
        else:
            LP = self.get_LP(i0 + 1, store=False)
        # multiply with `S`: a bit of a hack: use 'private' MPS._scale_axis_B
        S_bra = self.bra.get_SR(i0)
        if self.do_conj:
            S_bra = S_bra.conj()
        LP = self.bra._scale_axis_B(LP, S_bra, form_diff=1., axis_B='vR*', cutoff=0.)
        # cutoff is not used for form_diff = 1
        S_ket = self.ket.get_SR(i0)
        LP = self.bra._scale_axis_B(LP, S_ket, form_diff=1., axis_B='vR', cutoff=0.)
        RP = self.get_RP(i0, store=False)
        contr = npc.inner(LP, RP, axes=[['vR*', 'vR'], ['vL*', 'vL']], do_conj=False)
        return contr * self.bra.norm * self.ket.norm

    def expectation_value(self, ops, sites=None, axes=None):
        """Expectation value ``<bra|ops|ket>`` of (n-site) operator(s).

        Calculates n-site expectation values of operators sandwiched between bra and ket.
        For examples the contraction for a two-site operator on site `i` would look like::

            |          .--S--B[i]--B[i+1]--.
            |          |     |     |       |
            |          |     |-----|       |
            |          LP[i] | op  |       RP[i+1]
            |          |     |-----|       |
            |          |     |     |       |
            |          .--S--B*[i]-B*[i+1]-.

        Here, the `B` are taken from `ket`, the `B*` from `bra`.
        The call structure is the same as for :meth:`MPS.expectation_value`.

        .. warning ::

            In contrast to :meth:`MPS.expectation_value`, this funciton does not normalize,
            thus it also takes into account :attr:`MPS.norm` of both `bra` and `ket`.

        Parameters
        ----------
        ops : (list of) { :class:`~tenpy.linalg.np_conserved.Array` | str }
            The operators, for wich the expectation value should be taken,
            All operators should all have the same number of legs (namely `2 n`).
            If less than ``len(sites)`` operators are given, we repeat them periodically.
            Strings (like ``'Id', 'Sz'``) are translated into single-site operators defined by
            :attr:`sites`.
        sites : list
            List of site indices. Expectation values are evaluated there.
            If ``None`` (default), the entire chain is taken (clipping for finite b.c.)
        axes : None | (list of str, list of str)
            Two lists of each `n` leg labels giving the physical legs of the operator used for
            contraction. The first `n` legs are contracted with conjugated `B`,
            the second `n` legs with the non-conjugated `B`.
            ``None`` defaults to ``(['p'], ['p*'])`` for single site (n=1), or
            ``(['p0', 'p1', ... 'p{n-1}'], ['p0*', 'p1*', .... 'p{n-1}*'])`` for `n` > 1.

        Returns
        -------
        exp_vals : 1D ndarray
            Expectation values, ``exp_vals[i] = <bra|ops[i]|ket>``, where ``ops[i]`` acts on
            site(s) ``j, j+1, ..., j+{n-1}`` with ``j=sites[i]``.

        """
        ops, sites, n, (op_ax_p, op_ax_pstar) = self.ket._expectation_value_args(ops, sites, axes)
        ax_p = ['p' + str(k) for k in range(n)]
        ax_pstar = ['p' + str(k) + '*' for k in range(n)]
        E = []
        for i in sites:
            LP = self.get_LP(i, store=True)
            RP = self.get_RP(i + n - 1, store=True)
            op = self.ket.get_op(ops, i)
            op = op.replace_labels(op_ax_p + op_ax_pstar, ax_p + ax_pstar)
            C = self.ket.get_theta(i, n)
            C = npc.tensordot(op, C, axes=[ax_pstar, ax_p])  # same labels
            C = npc.tensordot(LP, C, axes=['vR', 'vL'])  # axes_p + (vR*, vR)
            C = npc.tensordot(C, RP, axes=['vR', 'vL'])  # axes_p + (vR*, vL*)
            C.ireplace_labels(['vR*', 'vL*'], ['vL', 'vR'])  # back to original theta labels
            theta_bra = self.bra.get_theta(i, n)
            E.append(npc.inner(theta_bra, C, axes='labels', do_conj=self.do_conj))
        return np.real_if_close(np.array(E)) * self.bra.norm * self.ket.norm

    def _contract_LP(self, i, LP):
        """Contract LP with the tensors on site `i` to form ``self.get_LP(i+1)``"""
        LP = npc.tensordot(LP, self.ket.get_B(i, form='A'), axes=('vR', 'vL'))
        axes = (self.ket._get_p_label('*') + ['vL*'], self.ket._p_label + ['vR*'])
        # for a ususal MPS, axes = (['p*', 'vL*'], ['p', 'vR*'])
        LP = npc.tensordot(self.bra.get_B(i, form='A').conj(complex_conj=self.do_conj), LP, axes=axes)
        return LP  # labels 'vR*', 'vR'

    def _contract_RP(self, i, RP):
        """Contract RP with the tensors on site `i` to form ``self.get_RP(i-1)``"""
        RP = npc.tensordot(self.ket.get_B(i, form='B'), RP, axes=('vR', 'vL'))
        axes = (self.ket._p_label + ['vL*'], self.ket._get_p_label('*') + ['vR*'])
        # for a ususal MPS, axes = (['p', 'vL*'], ['p*', 'vR*'])
        RP = npc.tensordot(RP, self.bra.get_B(i, form='B').conj(complex_conj=self.do_conj), axes=axes)
        return RP  # labels 'vL', 'vL*'

    def _to_valid_index(self, i):
        """Make sure `i` is a valid index (depending on `finite`)."""
        if not self._finite:
            return i % self.L
        if i < 0:
            i += self.L
        if i >= self.L or i < 0:
            raise KeyError("i = {0:d} out of bounds for MPSEnvironment".format(i))
        return i

    def expectation_value_multi_sites(self, operators, i0):
        r"""Expectation value  ``<psi|op0_{i0}op1_{i0+1}...opN_{i0+N}|psi>/<psi|psi>``.

        Calculates the expectation value of a tensor product of single-site operators
        acting on different sites next to each other.
        In other words, evaluate the expectation value of a term
        ``op0_i0 op1_{i0+1} op2_{i0+2} ...``, looking like this (with `op` short for `operators`,
        for ``len(operators)=3``):

            |          .--S--B[i0]---B[i0+1]--B[i0+2]--B[i0+3]--.
            |          |     |       |        |        |        |
            |          |     op[0]   op[1]    op[2]    op[3]    |
            |          |     |       |        |        |        |
            |          .--S--B*[i0]--B*[i0+1]-B*[i0+2]-B*[i0+3]-.


        .. warning ::
            This function does *not* automatically add Jordan-Wigner strings!
            For correct handling of fermions, use :meth:`expectation_value_term` instead.

        Parameters
        ----------
        operators : List of { :class:`~tenpy.linalg.np_conserved.Array` | str }
            List of one-site operators. This method calculates the
            expectation value of the n-sites operator given by their tensor
            product.
        i0 : int
            The left most index on which an operator acts, i.e.,
            ``operators[i]`` acts on site ``i + i0``.

        Returns
        -------
        exp_val : float/complex
            The expectation value of the tensorproduct of the given onsite operators,
            ``<psi|operators[0]_{i0} operators[1]_{i0+1} ... |psi>/<psi|psi>``,
            where ``|psi>`` is the represented MPS.
        """
        C = self._corr_ops_LP(operators, i0)
        exp_val = npc.trace(C, 'vR*', 'vR')
        return np.real_if_close(exp_val)

