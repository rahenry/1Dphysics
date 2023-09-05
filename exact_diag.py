#!/usr/bin/env python3

from tenpy.algorithms.exact_diag import ExactDiag
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.linalg.np_conserved import tensordot, inner, norm
import numpy
from rah_utility import rel_diff


class ExactDiagNonH(ExactDiag):
    def full_diagonalization_nonH(self, *args, **kwargs):
        if self.full_H is None:
            raise ValueError("You need to call one of `build_full_H_*` first!")
        E, V = npc.eig(self.full_H, *args, **kwargs)
        inds = numpy.argsort(E)
        # V.iset_leg_labels(['ps', 'ps*'])
        self.inds = inds
        self.E = E
        self.V = V

        Hleft = self.full_H.transpose()
        # Hleft = self.full_H.conj()
        Eleft, Vleft = npc.eig(Hleft, *args, **kwargs)
        indsleft = numpy.argsort(Eleft)
        # indsleft = []
        # for i in range(len(E)):
        #     k = inds[i]
        #     eR = E[k]
        #     for j in range(len(E)):
        #         eL = Eleft[j]
        #         if rel_diff(eL, eR) < 1E-10:
        #             indsleft.append(j)
        #             break
        # Eleft = Eleft[indsleft]
        self.indsleft = indsleft
        self.Eleft = Eleft
        self.Vleft = Vleft
        # for i in range(len(E)):
        #     j = inds[i]
        #     z = inner(Vleft[:,j], V[:,j])
        #     print(i, j, z, E[j])

    def partial_diag(self, *args, **kwargs):
        if self.full_H is None:
            raise ValueError("You need to call one of `build_full_H_*` first!")
        1
        print(420)
        exit()

    def full_to_mps(self, psi, canonical_form='B'):
        if not isinstance(psi.legs[0], npc.LegPipe):
            # projected onto charge_sector: need to restore the LegPipe.
            full_psi = npc.zeros([self._pipe], psi.dtype, psi.qtotal)
            full_psi[self._mask] = psi
            psi = full_psi
        psi.iset_leg_labels(['(' + '.'.join(self._labels_p) + ')'])
        psi = psi.split_legs([0])  # split the combined leg into the physical legs of the sites
        return MPS.from_full(self._sites, psi, form=canonical_form)
