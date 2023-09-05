#!/usr/bin/env python3

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

L = 10
N = 3
lamb = (L+1) / L
lamb = lamb ** (-N/2.)

print(lamb)

e = eps_k(numpy.pi, lamb, N)
print(e)
for a in numpy.linspace(0.21, 0.23, 20):
    k = numpy.pi + a * 1.j
    e = eps_k(k, lamb, N)
    print(e, k)
