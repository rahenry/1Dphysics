#!/usr/bin/env python3


from System import *
from SystemSet import *
from args import Q
from exact_fp import *
from rah_utility import *
from output_figs import *

import numpy, scipy
from fp_k import ep_test, plot1

def eps_to_k(e, g, N):
    # print(e, g, N)
    res = (e ** N) - 1. - (g ** N)
    res /= 2. * (g ** (N/2.))
    # res = -res
    # res = abs(res)
    res = numpy.arccos(res)
    return res

def eps_k(k, lamb, N):
    res = 1. + lamb**N + 2.*(lamb**(N/2.))*numpy.cos(k)
    return res ** (1./N)

def f_kj(k, lamb, L, N):
    return numpy.sin((L+1.) * k) + (lamb ** (-N/2.)) * numpy.sin(L*k)

def f_kj1(k, lamb, L, N):
    return (L+1) * numpy.cos((L+1.) * k) + (lamb ** (-N/2.)) * L * numpy.cos(L*k)

def f_kj2(k, lamb, L, N):
    return - (L+1.) * (L+1.) * numpy.sin((L+1.) * k) - (lamb ** (-N/2.)) * L*L *  numpy.sin(L*k)


npoints = 1000
k0 = -2*numpy.pi
k0 = -0.5
k1 = -k0
k0 = numpy.pi - 0.5
k1 = numpy.pi + 0.5
ks = numpy.linspace(k0, k1, npoints)
lamb = 0.88554880
# lamb = 0.90
# lamb = 0.88
# lamb = -lamb
L = 3
N = 3
# print(ks)
lamb = (-L / (L+1)) ** (2./N)
print(abs(lamb))
print(numpy.angle(lamb)/numpy.pi)
print(lamb)
print(f_kj1(0, lamb, L, N))
print(eps_k(numpy.pi, lamb, N))
# lamb = numpy.conj(lamb)
y = f_kj(ks, lamb, L, N)

plt.plot(ks, y, linestyle='', marker='o', markersize=2)
y2 = f_kj1(ks, lamb, L, N)
plt.plot(ks, y2, linestyle='', marker='o', markersize=2)
# plt.ylim([min(y), max(y)])
plt.gca().axhline(0)
plt.show()
exit()

phi = 1
lamb = abs(lamb)
print(lamb, phi)
delta, overlap, overlap2 = ep_test(L, N, lamb, phi)
silent_remove('graphs/test_k3')
mkdir('graphs/test_k3')
basename = os.path.join('graphs/test_k3', f'{N}_{L}')
name = basename + f'_{lamb}_{phi}_delta.png'
plot1(delta, name)
name = basename + f'_{lamb}_{phi}_overlap.png'
plot1(overlap, name)
name = basename + f'_{lamb}_{phi}_overlap2.png'
plot1(overlap2, name)
