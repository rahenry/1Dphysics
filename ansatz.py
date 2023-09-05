#!/usr/bin/env python3
import os

def interpret(x):
    name = x.split('(')[0]
    return name, x.replace(name, 'f')

def generate_ansatzs(files):
    res = {}
    for f in files:
        f = os.path.join('data', f'{f}.ansatz')
        if os.path.exists(f):
            data = open(f).readlines()
            for l in open(f).readlines():
                z = interpret(l)
                res[z[0]] = z[1]

    return res
class Ansatz:
    1
