#!/usr/bin/env python3

from exact_fp import *
import matplotlib.pyplot as plt
import matplotlib

N = 3
gammas = numpy.linspace(0.01, 0.99, 100)
phis = numpy.linspace(0, 3*N, 100)

def gap_inf(s):
    N = s['N']
    lamb = s['lambda']
    omega = numpy.exp(2.j * numpy.pi / N)
    phi = s['phi']

    best = 1E20
    for i in range(N-1):
        g = lamb * numpy.exp(phi * 2.j * numpy.pi / N)
        u = 1./ s['lambda_scale']
        u *= numpy.exp(-phi * 1.j * numpy.pi / N)
        res = (1. - g ** (N / 2.)) ** (2. / N)
        res = (1. - omega**(i+1)) * res * u
        if abs(res) < best:
            best = res

    return best

cm = matplotlib.cm.get_cmap('hsv')
data = numpy.zeros((len(gammas), len(phis)))
i = 0
for g in gammas:
    j = 0
    for p in phis:
        phi = p
        gamma = g
        lamb = gamma / (1.-gamma)
        lambsc = 1./(1.-gamma)
        s = {'N' : N,
             'lambda' : lamb,
             'phi'  : phi,
             'lambda_scale' : lambsc}

        z = gap_inf(s)
        z = abs(z)
        data[i][j] = z
        j += 1
    i += 1

P = plt.imshow(data)
plt.colorbar(P)
plt.show()
