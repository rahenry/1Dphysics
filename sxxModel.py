#!/usr/bin/env python3

#!/usr/bin/env python3

from tenpy.networks.site import Site
from tenpy.models.model import CouplingMPOModel
from tenpy.linalg import np_conserved as npc
import numpy as np
import math

class sxxSite(Site):
    def __init__(self, conserve=None):
        if conserve not in ["Q", "parity", None]:
            raise ValueError("invalid `conserve`: " + repr(conserve))

        Z = np.zeros([2, 2])
        Z[0][0] = 1
        Z[1][1] = -1
        ZP = Z
        X = np.zeros([2, 2])
        X[0][1] = 1
        X[1][0] = 1
        XP = X
        Y = np.zeros([2, 2], dtype=np.complex64)
        Y[0][1] = -1.0j
        Y[1][0] = 1.0j
        YP = Y

        ops = dict(Z=Z, ZP=ZP, X=X, XP=XP, Y=Y, YP=YP)
        Q = [1, -1]
        if conserve == "Q":
            chinfo = npc.ChargeInfo([2], ["Q"])
            leg = npc.LegCharge.from_qflat(chinfo, Q)
        else:
            if conserve == "parity":
                chinfo = npc.ChargeInfo([2], ["parity_Q"])
                leg = npc.LegCharge.from_qflat(chinfo, Q)
            else:
                leg = npc.LegCharge.from_trivial(2)
        self.conserve = conserve
        names = [str(i) for i in np.arange(2)]
        Site.__init__(self, leg, names, **ops)

    def __repr__(self):
        return f"sxx_site"

class sxxModel(CouplingMPOModel):
    def init_sites(self, model_params):
        self.model_params = model_params
        m = model_params
        self.L = m.get('L', 3)
        self.conserve = m.get('conserve', None)
        if self.conserve == 'None':
            self.conserve = None
        self.lamb = m.get('lambda', 1)
        self.h = m.get('h', 0)
        self.bc = m.get('bc', 0)
        # self.angle_style = m.get('angle_style', 1)
        exports = 'c1, c2, ch'
        self.export = exports.split()
        return sxxSite(self.conserve)

    def init_terms(self, model_params):
        self.c1 = -self.lamb
        self.c2 = -1.
        self.ch = -self.h
        L = self.L
        c1 = self.c1
        c2 = self.c2
        ch = self.ch

        g = model_params.get("gamma", 0)
        g = 0
        k = model_params.get("lambda", 1)
        h = model_params.get("h", 0)
        a = model_params.get("a", 1)
        a = 0
        for u in range(len(self.lat.unit_cell)):
            self.add_onsite(ch, u, "Z", plus_hc=False)

        for u1, u2, dx in self.lat.pairs["nearest_neighbors"]:
            wx1, wx2, wy1, wy2 = [1 + g, k * (1 + g), (1 - g), k * a * (1 - g)]
            wx1, wx2, wy1, wy2 = [c1, c2, c1, c2]
            i = 0
            # wxs = [wx2, wx1] * L
            # wxs = wxs[:L-1]
            # wys = [wy2, wy1] * L
            # wys = wys[:L-1]
            wxs = [wx2, wx1]
            wys = [wy2, wy1]
            self.add_coupling(wxs, u1, "X", u2, "X", dx, plus_hc=False)
            self.add_coupling(wys, u1, "Y", u2, "Y", dx, plus_hc=False)
        # if model_params['bc'] == 1:
        # # if model_params['bc'] == 1 or model_params['bc'] == 'infinite':
        #     c = c1
        #     if L % 2 == 0:
        #         c = c2
        #     self.add_coupling_term(c, 0, L-1, "X", "X", plus_hc=False)
        #     self.add_coupling_term(c, 0, L-1, "Y", "Y", plus_hc=False)
        #     1

        # for u in range(len(self.lat.unit_cell)):
        #     self.add_onsite(ch, u, "Z", plus_hc=False)

        # for u1, u2, dx in self.lat.pairs["nearest_neighbors"]:
        #     ws = [c1, c2] * L
        #     ws = ws[:L-1]
        #     self.add_coupling(ws, u1, "X", u2, "X", dx, plus_hc=False)
        #     # ws = [-c1, -c2] * L
        #     # ws = ws[:L-1]
        #     self.add_coupling(ws, u1, "Y2", u2, "Y2", dx, plus_hc=False)

        # if self.bc == 1:
        #     c = c1
        #     if L % 2 == 0:
        #         c = c2
        #     self.add_coupling_term(c, 0, L-1, "X", "X", plus_hc=False)
        #     self.add_coupling_term(c, 0, L-1, "Y", "Y", plus_hc=False)

