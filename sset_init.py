#!/usr/bin/env python3

import os
from args import Q
import dill
import numpy
import time
import matplotlib.pyplot as plt
from rah_utility import (
    # dict_product,
    # proc,
    mkdir,
    read_parameter_lines,
    dict_product,
    silent_remove,
)
from System import System
from Slicer import SystemSlicer
import pathos
# import multiprocessing
import multiprocess as mp
from report import report
import logging
import SystemSet_run
logging.basicConfig(filename="tenpy.log", filemode="w+", level=logging.INFO)

import cProfile, pstats, io
from pstats import SortKey
from data_set import DataSet
from fp_k import get_eps

import labels

class edata(dict):
    def __init__(self, d, h):
        self.hash = h
        super().__init__(d)

def find_eps(S):
    # print(S.config)
    # print(S.hash)
    self = S
    # systems_base = [System(x, self) for x in dict_product(self.config_systems)]
    self.systems = []
    self.eps_data = edata({}, self.hash)
    self.eps_data_file = os.path.join(self.result_dir, 'eps_data')
    if os.path.exists(self.eps_data_file) and not Q.rerun_eps:
        with open(self.eps_data_file, 'rb') as f:
            self.eps_data = dill.load(f)
    if not self.eps_data.hash == self.hash:
        self.eps_data = edata({}, self.hash)
    for c in dict_product(self.config_systems):
        L = c['L']
        N = c['N']
        cid = (L, N)

    for c in dict_product(self.config_systems):
        L = c['L']
        N = c['N']
        cid = (L, N)
        res = self.eps_data[cid]
        for lamb, phi in res['minima']:
            cnew = dict(c)
            cnew['lambda'] = lamb
            cnew['phi'] = phi
            s = System(cnew, self)
            self.systems.append(s)

def grid_eps(S):
    find_eps(S)
    S.systems = []
    S.grid_systems = {}
    for c in dict_product(S.config_systems):
        N, L = c['N'], c['L']
        cid = (L, N)
        S.grid_systems[cid] = []
        minima_base = S.eps_data[cid]['minima']

        minima = []
        for m in minima_base:
            if S['grid_style'] == 'circle':
                a = m[0] * numpy.exp(2.j * m[1] * numpy.pi / N)
                minima.append((a.real, a.imag))
        X = [m[0] for m in minima]
        Y = [m[1] for m in minima]
        xlims = numpy.array((min(X), max(X)))
        ylims = numpy.array((min(Y), max(Y)))
        if S['grid_lim']:
            g = S['grid_lim']
            xlims = [-g, g]
            ylims = [-g, g]
        n = S['grid_npoints']
        # print(xlims, ylims)
        for x in numpy.linspace(*xlims, n):
            for y in numpy.linspace(*ylims, n):
                cnew = dict(c)
                if S['grid_style'] == 'circle':
                    cnew['lambda_real'] = x
                    cnew['lambda_imag'] = y
                s = System(cnew, S)
                S.systems.append(s)
                S.grid_systems[cid].append(s)

    # print(len(S.systems))

def find_eps2(S):
    self = S
    systems_base = [System(x, self) for x in dict_product(self.config_systems)]
    s0 = systems_base[0]
    self.systems = []
    self.eps_data = edata({}, self.hash)
    self.eps_data_file = os.path.join(self.result_dir, 'eps_data')
    if os.path.exists(self.eps_data_file) and not Q.rerun_eps:
        with open(self.eps_data_file, 'rb') as f:
            self.eps_data = dill.load(f)
    if not self.eps_data.hash == self.hash:
        self.eps_data = edata({}, self.hash)
    conf = self.config
    confs = {}
    for k in ['L', 'N']:
        confs[k] = self.config_systems[k]
    for u in dict_product(confs):
        L = u['L']
        N = u['N']
        cid = (L, N)
        if cid in self.eps_data:
            continue
        print(f'{cid} running')
        c = self
        D = S.get("L_offset", 1)
        res = get_eps(L, N, c['k0'], c['k1'], c['q0'], c['q1'], 501, D=D)
        if s0['model_type'] == 'xyh' or s0["model_type"] == "xyh2":
            def cond1(k):
                return k.real >= numpy.pi / 2 or k.real < 0
                return k.real >= numpy.pi / 2 or k.real < 0 or k.imag < 0
            def cond2(k):
                if S["skip_trivial_EPs"]:
                    if numpy.isclose(abs(k), numpy.pi/2, atol=1E-2):
                        return False
                    if numpy.isclose(abs(k), 0, atol=1E-2):
                        return False
                return True

            def f2(k):
                if not cond2(k):
                    return None
                if numpy.isclose(abs(k), numpy.pi/2):
                    return None
                res = -numpy.sin((L+D)*k) / numpy.sin(L*k)
                return res

                res = 1
                if numpy.isclose(abs(k), numpy.pi/2):
                    res = -1
                    return 1E100
                else:
                    res = -numpy.sin((L+D)*k) / numpy.sin(L*k)
                if cond1(k):
                    return None

                # res  = (1 - res) / (1 + res)
                a = numpy.angle(res)
                phi = a * N / 2. / numpy.pi

                return res
                return abs(res), phi

            def f1(k):
                if not cond2(k):
                    return None
                if numpy.isclose(abs(k), numpy.pi/2):
                    return None
                res = -numpy.sin((L+D)*k) / numpy.sin(L*k)
                return 1./res
                res = 1
                if numpy.isclose(abs(k), numpy.pi/2):
                    res = -1
                    return 0
                    return 0, 0
                else:
                    res = -numpy.sin((L+D)*k) / numpy.sin(L*k)
                if cond1(k):
                    return None
                res = 1. / res
                # res  = (1 - res) / (1 + res)
                a = numpy.angle(res)
                phi = a * N / 2. / numpy.pi
                return res
                return abs(res), phi
            res['minima_k'] = [k for k in res['minima_k'] if cond2(k)]
            res['minima'] = [f1(x) for x in res['minima_k']] + [f2(x) for x in res['minima_k']]
            res["minima"] = [x for x in res["minima"] if x is not None]
        self.eps_data[cid] = res
    with open(self.eps_data_file, "wb+") as f:
        # print(f)
        dill.dump(self.eps_data, f)
    if S['add_regular_systems']:
        self.systems = [System(x, self) for x in dict_product(self.config_systems)]
    if S['add_ep_systems']:
        add_ep_systems(S)


def add_ep_systems(S):
    config_base = dict_product(S.config_systems)[0]
    S.ep_systems = {}
    for cid in S.eps_data:
        # print('...')
        # print(cid)
        S.ep_systems[cid] = []
        d = S.eps_data[cid]
        cbase = dict(config_base)
        cbase.update({'L' : cid[0], 'N' : cid[1]})

        eps = sorted(d['minima'], key=lambda x: x[1])
        i = 0
        for lamb, phi in eps:
            # print(lamb, phi, i)
            c = dict(cbase)
            c['lambda'] = lamb
            c['phi'] = phi
            c['ep_ind'] = i
            c['ep_ind_offset'] = i - cid[0] // 2
            s = System(c, S)
            S.systems.append(s)
            i += 1


def add_ep_line(S):
    # print(S)
    for cid in S.eps_data:
        1
        # print(cid)
        d = S.eps_data[cid]

    exit()

    1
