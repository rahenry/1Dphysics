#!/usr/bin/env python3

#!/usr/bin/env python3

from ZnSite import ZnSite
from tenpy.models.model import CouplingMPOModel
import numpy as np

class ZnModel(CouplingMPOModel):
    def init_sites(self, model_params):
        self.model_params = model_params
        m = model_params
        self.N = m.get('N', 3)
        self.L = m.get('L', 3)
        self.conserve = m.get('conserve', 'Q')
        if self.conserve == 'None':
            self.conserve = None
        self.alpha = m.get('alpha', 1)
        self.gamma = m.get('gamma', 0.5)
        self.eta = m.get('eta', 0)
        self.k = m.get('k', 1)
        self.thetaX = m.get('thetaX', 0)
        self.thetaY = m.get('thetaY', 0)
        self.phi = m.get('phi', 0)
        self.hc = m.get('hc', 1)
        self.conj = m.get('conj', 0)
        self.swapXZ = m.get('swapXZ', 0)
        self.flipZ = m.get('flipZ', 0)
        self.switchConjTerm = m.get('switchConjTerm', 1)
        return ZnSite(self.N, self.conserve)

    def init_terms(self, model_params):
        g2 = self.gamma
        g1 = (1.-self.gamma)
        cX = -self.alpha * g1 * (1.0 - self.eta) * np.exp(-np.pi * 1.0j * self.thetaX)
        cY = -self.alpha * g1 * self.eta * np.exp(-np.pi * 1.0j * self.thetaY)
        cZ = -self.alpha * g2 * np.exp(-np.pi * 1.0j * self.phi)
        if self.flipZ:
            cZ *= -1
        out = f"Couplings:\ncX = {cX}\ncY = {cY}\ncZ = {cZ}"

        o1 = "Z"
        o3 = "X"
        if self.swapXZ:
            o1 = "X"
            o3 = "Z"
        o2 = o1
        o4 = o3
        if self.conj:
            o1 += "P"
            if self.switchConjTerm:
                o3 += "P"
            else:
                o4 += "P"
        else:
            o2 += "P"
            if self.switchConjTerm:
                o4 += "P"
            else:
                o3 += "P"
        for u in range(len(self.lat.unit_cell)):
            # Z
            self.add_onsite(cZ, u, o1, plus_hc=self.hc)
            # self.add_onsite(cZ, u, o2, plus_hc=self.hc)

        for u1, u2, dx in self.lat.pairs["nearest_neighbors"]:
            # XP X
            # self.add_coupling(cX, u1, o3, u2, o4, dx, plus_hc=self.hc)
            self.add_coupling(cX, u1, o3, u2, o3, dx, plus_hc=self.hc)
            # self.add_coupling(cY, u1, "YP", u2, "Y", dx, plus_hc=self.hc)
            # self.add_coupling(cY, u1, "Y", u2, "YP", dx, plus_hc=self.hc)
            # self.add_coupling(cY, u1, "X", u2, "X", dx, plus_hc=self.hc)


if __name__ == "__main__":
    z = ZnModel({})
    print(z)
