#!/usr/bin/env python3
import math
from exact_fp import *
import matplotlib.pyplot as plt
from rah_utility import rel_diff
from numpy.polynomial import Polynomial
from System import System


def test_present(e, exclude):
    for p in exclude:
        if rel_diff(e, p) < 1E-7:
            return True
    return False

def pf_test(s):

    # energies = s.eng.E
    energies = s['energies']
    energies = sorted(energies, key=lambda x: x.real)
    # deltas = [x - s.eng.E[0] for x in energies]
    deltas = [x - s['energies'][0] for x in energies]
    candidates = [x / (1.5 + 1.j * math.sqrt(3) / 2.) for x in deltas[1:]]
    # print(candidates[:10])

    pfs = []
    exclude = []
    for d in deltas[1:]:
        if d.imag < 0: continue
        p = d / (1.5 + 1.j * math.sqrt(3) / 2.)
        test = test_present(d, exclude)
        print(d, test)
        if not test:
            pfs.append(p)
            exclude = make_pf_spectrum2(0, pfs, s['N'])
        if len(pfs) == s['L']:
            break

    # pfs = energies[1:s['L']+1]
    # pfs = [x - energies[0] for x in pfs]
    # pfs = [x.real * 2./3. for x in pfs]

    # spectrum = make_pf_spectrum2(s.eng.E[0], pfs, s['N'])
    spectrum = make_pf_spectrum2(s['energies'][0], pfs, s['N'])
    spectrum = sorted(spectrum, key=lambda x: x.real)
    # print(s['e_exact'], s.eng.E[0] / s['L'])
    # for i in range(len(pfs)):
    #     print(pfs[i], s.pfs[i] * s['L'])


    pfs_poly = find_eps_poly(s)
    spectrum = make_pf_spectrum2(s['energies'][0], pfs_poly, s['N'])
    spectrum = sorted(spectrum, key=lambda x: x.real)
    plt.clf()
    spectrum_x = [x.real for x in spectrum]
    spectrum_y = [x.imag for x in spectrum]
    engs_x = [x.real for x in s['energies']]
    engs_y = [x.imag for x in s['energies']]
    plt.plot(spectrum_x, spectrum_y, linestyle='', marker='x', markersize=6, color='red')
    plt.plot(engs_x, engs_y, linestyle='', marker='o', markersize=1.5, color='blue')
    plt.savefig('pf_test.png')
    plt.show()



def find_eps_poly(s):
    eta = s['eta'] * 0.5
    P0 = Polynomial([1])
    polys = [P0]
    i = 0
    Mmax = s['L'] * 2
    N = s['N']
    L = s['L']
    M = 2*L - 2
    M = L - 1
    p = 1
    Mbar = numpy.floor((M+p)/(p+1))
    # print(Mmax, Mbar)
    cX = s.model.cX
    cY = s.model.cY
    cW = s.model.cW
    while(True):
        # print(polys[-1])
        i += 1
        kappa = cX
        if (i % 3) == 2:
            kappa = cY
        if (i % 3) == 0:
            kappa = cW
        kappa = (-kappa) ** N
        # if i%2==0:
        #     kappa = 1.-eta
        Pnew = polys[-1].copy()
        if len(polys) > p:
            Pnew -= kappa * Polynomial([0, 1]) * polys[-1-p]
        else:
            Pnew -= kappa * Polynomial([0, 1])
        polys.append(Pnew)
        if len(polys) == M+1:
            break

    # print(polys[-1])
    roots = polys[-1].roots()
    # roots = []
    # for p in polys:
    #     z = p.roots()
    #     for y in z:
    #         roots.append(y)
    # pfs = [r**(-1./s['N']) for r in roots]
    # print(roots)
    pfs = [numpy.complex64(r)**(-1./s['N']) for r in roots]
    # pfs2 = s.pfs
    # a, b, pfs2 = exact_parafermions_matrix(L, N, 1.0, flipZ=1)
    # pfs2 = [p * L * 0.5 for p in pfs2]
    # print(N, L, eta)
    # print(pfs2)
    # print(pfs)
    return pfs


if __name__ == "__main__":
    eta = 0.7
    kappa = 0.6
    N = 3
    L = 20
    c = {
        'N' : N,
        'L' : L,
        'eta' : eta,
        'kappa' : eta,
        }
    s = System(c)
    pfs = find_eps_poly(s)
    print(pfs)
    1

    pfs_poly = find_eps_poly(s)
    spectrum = make_pf_spectrum2(0, pfs_poly, s['N'])
    spectrum = sorted(spectrum, key=lambda x: x.real)
    spectrum_x = [x.real for x in spectrum]
    spectrum_y = [x.imag for x in spectrum]
    plt.plot(spectrum_x, spectrum_y, linestyle='', marker='o', markersize=1, color='red')
    plt.show()
