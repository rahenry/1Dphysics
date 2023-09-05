#!/usr/bin/env python3


from System import *
from SystemSet import *
from args import Q

import numpy, scipy

def find_minima(z):
    X, Y = z.shape
    res = []
    for i in range(X-2):
        x = i+1
        for j in range(Y-2):
            y = j+1
            u = z[x][y+1]
            d = z[x][y-1]
            l = z[x-1][y]
            r = z[x-1][y]
            t = z[x][y]
            if t < u and t < d and t < l and t < r:
                res.append((x, y))

    return res
z = numpy.random.rand(5,5)
print(z)
ex = find_minima(z)
print(ex)

def make_M():
    S = SystemSet("eps_test1")
    S.run()
    S.purge_graphs()
    S.extras()

    gammas = S.slicer.get_possible_values('gamma')
    phis = S.slicer.get_possible_values('phi')
    L = S['L']
    s = S.systems[0]
    nsys = len(S.systems)
    Y = numpy.full((nsys, L+2), 0.)
    for i in range(nsys):
        s = S.systems[i]
        Y[i][0] = s['gamma']
        Y[i][1] = s['phi']
        for j in range(L):
            Y[i][j+2] = s['pfs'][j]
    return Y

f = 'M1'
X = None
Y = None
if os.path.exists(f) and not Q.rerun_measurements:
    print("Loading")
    with open(f, 'rb') as f1:
        M = dill.load(f1)
else:
    M = make_M()
    with open(f, 'wb+') as f1:
        dill.dump(M, f1)
        print("saving")
