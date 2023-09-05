#!/usr/bin/env python3

from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
import numpy as np

class ZnSite(Site):
    def __init__(self, N=2, conserve="Q"):
        if conserve not in ["Q", "parity", None]:
            raise ValueError("invalid `conserve`: " + repr(conserve))
        self.N = N
        I = np.zeros([N, N], dtype=np.complex64)
        for i in range(N):
            I[i][i] = 1
        Z = np.zeros([N, N], dtype=np.complex64)
        for i in range(N):
            Z[i][i] = np.exp(2.0j * i * np.pi / N)
        ZP = np.zeros([N, N], dtype=np.complex64)
        for i in range(N):
            ZP[i][i] = np.exp(-2.0j * i * np.pi / N)
        X = np.zeros([N, N])
        XP = np.zeros([N, N])
        for i in range(N):
            iP = (i + 1) % N
            XP[i][iP] = 1
            X[iP][i] = 1
        Y = np.zeros([N, N], dtype=np.complex64)
        for i in range(N):
            iP = (i + 1) % N
            Y[iP][i] = np.exp(2.0j*i*np.pi/N)
        YP = Y.conj().T
        W = np.zeros([N, N], dtype=np.complex64)
        for i in range(N):
            iP = (i + 1) % N
            W[iP][i] = np.exp(-2.0j*i*np.pi/N)
        WP = W.conj().T

        Y2 = np.zeros([2, 2], dtype=np.complex64)
        Y2[0][1] = -1.0j
        Y2[1][0] = 1.0j
        # A = ZP * YP
        # B = Z * YP

        Q = np.arange(N)
        # ops = dict(Z=Z, ZP=ZP, X=X, XP=XP, Y=Y, YP=YP, W=W, WP=WP, A=A, B=B)
        ops = dict(Z=Z, ZP=ZP, X=X, XP=XP, Y=Y, YP=YP, W=W, WP=WP, I=I, Y2=Y2)
        for i in range(1, N):
            Zpow = np.linalg.matrix_power(Z, i)
            ops[f"Z{i}"] = Zpow
            ZPpow = np.linalg.matrix_power(ZP, i)
            ops[f"ZP{i}"] = ZPpow
            Xpow = np.linalg.matrix_power(X, i)
            ops[f"X{i}"] = Xpow
            XPpow = np.linalg.matrix_power(XP, i)
            ops[f"XP{i}"] = XPpow
            Ypow = np.linalg.matrix_power(Y, i)
            ops[f"Y{i}"] = Ypow
            YPpow = np.linalg.matrix_power(YP, i)
            ops[f"YP{i}"] = YPpow
            Wpow = np.linalg.matrix_power(W, i)
            ops[f"W{i}"] = Wpow
            WPpow = np.linalg.matrix_power(WP, i)
            ops[f"WP{i}"] = WPpow
        if conserve == "Q":
            chinfo = npc.ChargeInfo([N], ["Q"])
            leg = npc.LegCharge.from_qflat(chinfo, Q)
        else:
            if conserve == "parity":
                chinfo = npc.ChargeInfo([N], ["parity_Q"])
                leg = npc.LegCharge.from_qflat(chinfo, Q)
            else:
                leg = npc.LegCharge.from_trivial(N)
                if N == 3:
                    P = np.zeros([N, N])
                    P[0][0] = 1
                    P[1][2] = 1
                    P[2][1] = 1
                    ops["P"] = P

        self.conserve = conserve
        names = [str(i) for i in np.arange(N)]
        Site.__init__(self, leg, names, **ops)

    def __repr__(self):
        return f"ZnSite(N={self.N})"

if __name__ == "__main__":
    z = ZnSite(2, None)
    print(z)
    z2 = ZnSite(3, None)
    print(z2)
    z3 = ZnSite(3, "Q")
    print(z3)
