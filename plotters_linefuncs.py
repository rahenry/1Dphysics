#!/usr/bin/env python3

from plotters_utility import *
from scipy.special import hyp2f1

def lambda_to_gamma(l):
    return l / (1.+l)

def gamma_to_lambda(g):
    return g / (1.-g)

def rho(N, z):
    if z< 1:
        return hyp2f1(-1./N, (N-1.)/N, 1, z**N)
    else:
        return 1./N*(z**(1.-N))*hyp2f1((N-1.)/N, (N-1.)/N, 2, z**(-N))

def rhoX(l, N=3):
    return rho(N, l)

def rhoZ(l, N=3):
    return rho(N, 1./l)

def rhoXg(l, N=3):
    l = gamma_to_lambda(l)
    return rho(N, l)

def rhoZg(l, N=3):
    l = gamma_to_lambda(l)
    return rho(N, 1./l)
