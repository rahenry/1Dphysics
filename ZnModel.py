#!/usr/bin/env python3

#!/usr/bin/env python3

from ZnSite import ZnSite
from tenpy.models.model import CouplingMPOModel
import numpy as np
import math

def SICP_factor(r, N):
    u = np.exp(math.pi * 1.J * (2.*r-N)/2./N)/np.sin(math.pi*r/N)
    return u

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
        self.gamma_full = m.get('gamma_full', 0.5)
        self.eta_full = m.get('eta_full', 0.5)
        self.eta = m.get('eta', 0)
        self.kappa = m.get('kappa', 0)
        self.mu = m.get('mu', -1)
        self.k = m.get('k', 1)
        self.thetaX = m.get('thetaX', 0)
        self.thetaY = m.get('thetaY', 0)
        self.thetaW = m.get('thetaW', 0)
        self.phi = m.get('phi', 0)
        self.hc = m.get('hc', 0)
        self.conj = m.get('conj', 0)
        self.swapXZ = m.get('swapXZ', 0)
        self.flipZ = m.get('flipZ', 0)
        self.switchConjTerm = m.get('switchConjTerm', 0)
        self.model_type = m.get('model_type', 'fp')
        self.bc = m.get('bc', 0)
        self.Mplus = m.get('Mplus', 0)
        self.h = m.get('h', 0.)
        # self.angle_style = m.get('angle_style', 1)
        exports = 'cX cY cW cZ p M lambdas'
        self.export = exports.split()
        self.scale = m.get('scale', 'lambda')
        return ZnSite(self.N, self.conserve)


    def init_terms(self, model_params):
        # a1 = np.exp(-np.pi*1.0j*self.phi/self.N)
        # a2 = np.exp(np.pi*1.0j*self.phi/self.N)
        # if self.angle_style == 0:
        #     a1 = 1
        #     a2 = np.exp(np.pi*2.0j*self.phi/self.N)
        # if self.angle_style == 2:
        #     a1 = np.exp(np.pi*2.0j*self.phi/self.N)
        #     a2 = 1
        g1 = (1.-self.gamma)
        g2 = self.gamma  * np.exp(np.pi * 2.j * self.phi / self.N)
        sc = 1.
        if self.scale == 'lambda' and not self.gamma == 1:
            sc = 1. / (1.-self.gamma)
        sc *= np.exp(self.mu * self.phi * 2.j * np.pi / self.N)

        cX = -self.alpha * g1 * (1.0 - self.eta) * (1.0 - self.kappa) * np.exp(np.pi * 2.0j * self.thetaX) * sc
        cY = -self.alpha * g1 * self.eta * (1.0 - self.kappa) * np.exp(np.pi * 2.0j * self.thetaY) * sc
        cW = -self.alpha * g1 * self.eta * self.kappa * np.exp(np.pi * 2.0j * self.thetaW) * sc
        cZ = -self.alpha * g2 * sc
        N = self.N
        self.cX = cX
        self.cY = cY
        self.cZ = cZ
        self.cW = cW
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
                o4 += "P"
            else:
                o3 += "P"
        else:
            o2 += "P"
            if self.switchConjTerm:
                o3 += "P"
            else:
                o4 += "P"

        c = self.model_type
        if c == 'fp':
            self.lambdas = [(-cZ)**N, (-cX)**N]
            for u in range(len(self.lat.unit_cell)):
                # Z
                self.add_onsite(cZ, u, o1, plus_hc=self.hc)

            for u1, u2, dx in self.lat.pairs["nearest_neighbors"]:
                # XP X
                self.add_coupling(cX, u1, o3, u2, o4, dx, plus_hc=self.hc)
                # self.add_coupling(cX, u1, o4, u2, o3, dx, plus_hc=self.hc)
            # if self.bc == 1:
            #     self.add_coupling_term(cX, 0, self.L-1, o4, o3, plus_hc=self.hc)
                # self.add_coupling_term(cX, 0, o4, self.L-1, o3, plus_hc=self.hc)
        elif c == 'fpXYW':
            self.M = self.L - 1 + self.Mplus
            if self.bc == 1:
                self.M += 1
            self.p = 1
            self.lambdas = [(-cX)**N, (-cY)**N, (-cW)**N]
            for u in range(len(self.lat.unit_cell)):
                o1 = "Z"
                # Z
                self.add_onsite(cZ, u, o1, plus_hc=self.hc)
            for i in range(self.L-1):
                if i % self.N == 0:
                    self.add_coupling_term(cX, i, i+1, 'XP', 'X')
                if i % self.N == 1:
                    self.add_coupling_term(cY, i, i+1, 'YP', 'Y')
                if i % self.N == 2:
                    self.add_coupling_term(cW, i, i+1, 'WP', 'W')
                # self.add_coupling_term(cZ, i, i+1, "XP", "X")
            if self.bc == 1:
                i = self.L-1
                if i % self.N == 0:
                    self.add_coupling_term(cX, 0, self.L-1, 'X', 'XP')
                if i % self.N == 1:
                    self.add_coupling_term(cY, 0, self.L-1, 'Y', 'YP')
                if i % self.N == 2:
                    self.add_coupling_term(cW, 0, self.L-1, 'W', 'WP')
                # self.add_coupling_term(cX, 0, self.L-1, 'X', 'XP')
        elif c == 'fptest':
            for u in range(len(self.lat.unit_cell)):
                # self.add_onsite([cZ, 0], u, 'Z', plus_hc=self.hc)
                # self.add_onsite([0, cZ], u, 'W', plus_hc=self.hc)
                self.add_onsite(cZ, u, 'Z', plus_hc=self.hc)
            for u1, u2, dx in self.lat.pairs["nearest_neighbors"]:
                # XP X
                a = 1.
                self.add_coupling(a*cX, u1, "XP", u2, "X", dx, plus_hc=self.hc)
                # self.add_coupling(a*cY, u1, "YP", u2, "Y", dx, plus_hc=self.hc)
                # self.add_coupling(a*cW, u1, "W", u2, "WP", dx, plus_hc=self.hc)
            if self.bc == 1:
                self.add_coupling_term(cX, 0, self.L-1, 'X', 'XP')
        elif c == 'ident':
            for u in range(len(self.lat.unit_cell)):
                self.add_onsite(1, u, 'I', plus_hc=0)
        elif c == 'fptest2':
            for u in range(len(self.lat.unit_cell)):
                self.add_onsite(cZ, u, f'{o1}', plus_hc=self.hc)
            for u1, u2, dx in self.lat.pairs["nearest_neighbors"]:
                # XP X
                self.add_coupling(cX, u1, f'{o3}', u2, f'{o4}', dx, plus_hc=self.hc)
                # self.add_coupling(cX, u1, o4, u2, o3, dx, plus_hc=self.hc)
            if self.bc == 1:
                self.add_coupling_term(cX, 0, self.L-1, f'{o4}', f'{o3}', plus_hc=self.hc)
        elif c == 'SICP':
            for i in range(1, N):
                u = SICP_factor(i, N)
                a = -g1 * u * sc
                b = -g2 * u * sc
                for u in range(len(self.lat.unit_cell)):
                    # Z
                    self.add_onsite(b, u, f'{o1}{i}', plus_hc=self.hc)
                for u1, u2, dx in self.lat.pairs["nearest_neighbors"]:
                    # XP X
                    self.add_coupling(a, u1, f'{o3}{i}', u2, f'{o4}{i}', dx, plus_hc=self.hc)
                    # self.add_coupling(cX, u1, o4, u2, o3, dx, plus_hc=self.hc)
                if self.bc == 1:
                    self.add_coupling_term(b, 0, self.L-1, f'{o4}{i}', f'{o3}{i}', plus_hc=self.hc)

        elif c == 'fpzeros':
            for u in range(len(self.lat.unit_cell)):
                self.add_onsite(1, u, o1, plus_hc=self.hc)
        elif c == 'xyh':
            # g = self.gamma * np.exp(2.j*np.pi * self.phi / self.N)
            g = self.gamma_full
            aX = -0.25 * (1.+g)
            aY = -0.25 * (1.-g)
            aZ = self.h
            for u1, u2, dx in self.lat.pairs["nearest_neighbors"]:
                self.add_coupling(aX, u1, "XP", u2, "X", dx, plus_hc=self.hc)
                self.add_coupling(aY, u1, "YP", u2, "Y", dx, plus_hc=self.hc)
            for u in range(len(self.lat.unit_cell)):
                self.add_onsite(-self.h, u, 'Z', plus_hc=self.hc)
            if self.bc == -1:
                self.add_coupling_term(-aX, 0, self.L-1, "X", "XP", plus_hc=self.hc)
                self.add_coupling_term(-aY, 0, self.L-1, "Y", "YP", plus_hc=self.hc)
            if self.bc == -2:
                self.add_coupling_term(aX, 0, self.L-1, "X", "XP", plus_hc=self.hc)
                self.add_coupling_term(aY, 0, self.L-1, "Y", "YP", plus_hc=self.hc)
        elif c == 'xyh2':
            # g = self.gamma * np.exp(2.j*np.pi * self.phi / self.N)
            g = self.eta_full
            aX = 1.
            aY = g
            aZ = self.h
            for u1, u2, dx in self.lat.pairs["nearest_neighbors"]:
                self.add_coupling(aX, u1, "XP", u2, "X", dx, plus_hc=self.hc)
                self.add_coupling(aY, u1, "YP", u2, "Y", dx, plus_hc=self.hc)
            for u in range(len(self.lat.unit_cell)):
                self.add_onsite(self.h, u, 'Z', plus_hc=self.hc)
            # if self.bc == -1:
            #     self.add_coupling_term(-aX, 0, self.L-1, "X", "XP", plus_hc=self.hc)
            #     self.add_coupling_term(-aY, 0, self.L-1, "Y", "YP", plus_hc=self.hc)
            # if self.bc == -2:
            #     self.add_coupling_term(aX, 0, self.L-1, "X", "XP", plus_hc=self.hc)
            #     self.add_coupling_term(aY, 0, self.L-1, "Y", "YP", plus_hc=self.hc)


if __name__ == "__main__":
    z = ZnModel({})
    print(z)
