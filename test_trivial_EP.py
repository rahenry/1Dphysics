#!/usr/bin/env python3

from System import *
from SystemSet import *
from args import Q
from exact_fp import *
from rah_utility import *
from output_figs import *
from fp_k import *
from matrix import *

import numpy, scipy

basename = "fin1"

d = f"graphs/test_ks2_{basename}"
silent_remove(d)
mkdir(d)
Ls = [10]
Ns = [2, 3, 4, 5, 6]
q1 = 3.5
q0 = -0.1
zpow = 0.2
vmax = None
vmin = 0
k0 = -0.5
k1 = 2.5
npoints = 101
N = 2
ks = numpy.linspace(k0, k1, npoints)
Ls = [5]
lamb = 0.9384
lamb = 0.8
lamb = 1.1
lamb = 0.9384
lamb = 0.1


basedir = 'trivial_EP_test'
graph_dir = os.path.join('graphs', basedir)
silent_remove(graph_dir)
mkdir(graph_dir)
L = 10
l = ((L/(L+1.))**2)**(1./3.)
g = -(l ** (-3./2.))
kk = 1E-9
a = numpy.sin(kk*(L+1.)) + g * numpy.sin(L*kk)
b = (L+1.) * numpy.cos((L+1.)*kk) + g * L * numpy.cos(L*kk)
print(l, g)
print(a, b)
# exit()
# for L in Ls:
#     eps = get_eps(L, N, k0, k1, q0, q1)
#     for x, y in eps.items():
#         print(x)
#     # print(eps['minima_k'])
#     for i in range(len(eps['minima'])):
#         e = eps['minima'][i]
#         k = eps['minima_k'][i]
#         print(e, k)
for L in Ls:
    plt.clf()

    # ydata = f_kj(ks, lamb, L, N)
    ydata = [f_kj(k, lamb, L, N) for k in ks]
    plt.plot(ks, ydata)

    # ydata = abs(ydata)
    ydata = [1./L*f_kj1(k, lamb, L, N) for k in ks]
    plt.plot(ks, ydata)

    ydata = h(ks, L)
    ydata = abs(ydata)
    plt.plot(ks, ydata)

    ydata = h2(ks, L)
    ydata = abs(ydata)
    plt.plot(ks, ydata)

    y = 0 * ks
    plt.plot(ks, y, linestyle='--', c='gray', alpha=0.5)

    # y = []
    # for k in ks:
    #     z = numpy.sin((L+1)*k) + numpy.sin(L*k)
    #     y.append(z)
    # plt.plot(ks, y, c='red')

    f = os.path.join(graph_dir, f'{L}.png')
    plt.savefig(f)
