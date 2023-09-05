#!/usr/bin/env python3

import numpy

omega = numpy.exp(2.*numpy.pi*1.j/3.)
X = numpy.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
Z = numpy.array([[1, 0, 0], [0, omega, 0], [0, 0, omega*omega]])
Y = numpy.array([[0, 0, 1], [omega,0,0],[0,omega*omega,0]])
W = numpy.array([[0, 0, 1], [omega*omega,0,0],[0,omega,0]])
# W = numpy.array([[0, 1, 0], [0,0,omega],[omega*omega,0,0]])
# Y = X@Z - Z@X
# Y2 = Y / numpy.sqrt(3) / 1.j
# Y = Y2
# Y = Y*Y
# Y = Y.conj().T
def com(a, b):
    return a*b - b*a
def acom(a, b):
    return a*b + b*a

c1 = com(X, Z)



eigs = numpy.linalg.eig(W)
def process_elem(x, evec):
    if abs(x) < 1E-10: return (0, 0)
    base = evec[0]
    if abs(base) < 1E-10: base = 1
    a = abs(x)
    a = numpy.around(a, 5)
    b = numpy.angle(x/base)
    # b = numpy.angle(x)
    b = numpy.around(b/2./numpy.pi, 5)
    return a,b
def printy(things, eigs):

    res = ''
    for i in range(len(things[0])):
        b = eigs[i]
        res += f'{b.real:>5.2f}+{b.imag:>5.2f}i   '
    res += '\n'
    for i in range(len(things[0])):
        for v in things:
            a, b = v[i]
            # res += f'{int(a)} {b:>5.2f}\n'
            res += f'{b:>12.2f}   '
        res += '\n'
    print(res)
    return res

# Y = omega*Y
# Y = omega*Y
ops = [Z, X, Y, W]
ops = {
    'Z':Z,
    'X':X,
    'Y':Y,
    'W':W,
    }
# ops = {
#     'Y':Y,
#     'Q':Q,
#     }
# ops = {
#     'Z':Z,
#     'X':X,
#     }
for name, o in ops.items():
    eigs, evecs = numpy.linalg.eig(o)
    print(eigs)
for name, o in ops.items():
    print(name)
    eigs, evecs = numpy.linalg.eig(o)
    things = []
    for i in range(len(eigs)):
        v = numpy.array(evecs[:,i])
        v = [process_elem(x, v) for x in v]
        things.append(v)
    printy(things, eigs)
    print(name+'D')
    eigs, evecs = numpy.linalg.eig(o.conj().T)
    things = []
    for i in range(len(eigs)):
        v = numpy.array(evecs[:,i])
        v = [process_elem(x, v) for x in v]
        things.append(v)
    printy(things, eigs)


# W = W.conj().T
print((X@Z) / (Z@X))
print((X@Y) / (Y@X))
print((Y@Z) / (Z@Y))
print((Z@W) / (W@Z))
print((X@W) / (W@X))
print((Y@W) / (W@Y))
def commute_test2(A, B):
    a = commute_test(A, B)
    b = commute_test(B, A)
    print(a, b)
def commute_test(A, B):
    p1 = (A@B) / (B@A)
    rep = None
    for x in p1:
        found = False
        for y in x:
            if not numpy.isnan(y):
                rep = y
                found = True
                break
        if found: break
    return rep
commute_test2(X, Z)
commute_test2(X, Y)
commute_test2(Y, Z)
commute_test2(W, Z)
commute_test2(W, X)
commute_test2(Y, W)
print(X)
print(Z)
print(Y)
print(W)
print(Z@Y)
print(Z@W)
print(Y@W)
