#!/usr/bin/env python3
import os
from args import Q
import dill
import numpy
import time
import matplotlib.pyplot as plt
import json
from rah_utility import (
    # dict_product,
    # proc,
    mkdir,
    read_parameter_lines,
    read_parameter_lines_str,
    dict_product,
    silent_remove,
    simple_hash,
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
from pbc_zeros import analyse_zeros
from periodic import *

import labels
from sset_init import *

class SystemSet:
    def process_input(self, input_name):
        input_file = os.path.join("data", input_name)
        if not os.path.exists(input_file):
            print(f'{input_file} not found')
            exit()

        config = {}
        config_systems = {}

        for l in open(input_file).readlines():
            l = l.split()
            if not len(l) == 2:
                continue
            x, y = l
            if x == 'input' and not y==input_name:
                self.process_input(y)

        lines = []
        for l in open(input_file).readlines():
            l = l.strip('\n')
            if l == '':
                continue
            if l[0] == '#':
                continue
            sp = l.split()
            if '$NAME' in l:
                l = l.replace('$NAME', input_name)
            if sp[0][0] == '$' and sp[0] not in self.input_vars:
                self.input_vars[sp[0]] = ' '.join(sp[1:])
                continue
            lines.append(l)
        lines_new = []
        for l in lines:
            for v in self.input_vars:
                if v in l:
                    l = l.replace(v, self.input_vars[v])
            lines_new.append(l)
        lines = lines_new

        lines_new = []
        lines_output = []
        self.subsets = []
        output_found = 0
        in_subset = 0
        for l in lines:
            if l == 'OUTPUT':
                output_found = 1
                continue
            if output_found:
                lines_output.append(l)
            else:
                if l == "SUBSET":
                    self.subsets.append([])
                    in_subset = 1
                    continue
                if in_subset:
                    if l == "END":
                        in_subset = 0
                        continue
                    else:
                        self.subsets[-1].append(l)
                else:
                    lines_new.append(l)
        lines = lines_new
        self.lines_raw = lines
        self.output_spec = lines_output


        # for x, y in read_parameter_lines(input_file).items():
        for x, y in read_parameter_lines_str(lines).items():
            if x[-1] == "'":
                self.config[x.strip("'")] = y[0]
            elif x[-1] == '"':
                self.config[x.strip('"')] = y
            else:
                self.config[x] = y[0]
                self.config_systems[x] = y

        self.config.update(config)
        self.hash = simple_hash(json.dumps(self.config))
        self.config_systems.update(config_systems)

    @classmethod
    def from_str(cls, c):
        open('data/tempsys', 'w+').write(c)
        return cls('tempsys')

    def __getitem__(self, key):
        if key in self.config:
            return self.config[key]
        return []
        # s = self.systems[0]
        # return s[key]

    def renew_names(self):
        # *** PROFILING CODE ***
        # pr = cProfile.Profile()
        # pr.enable()
        # self.slicer = SystemSlicer(self)
        # pr.disable()
        # s = io.StringIO()
        # sortby = SortKey.CUMULATIVE
        # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        # ps.print_stats()
        # print(s.getvalue())
        #

        self.output("Renewing names")
        self.slicer = SystemSlicer(self)
        u = self.slicer.relevant_identifiers
        for s in self.systems:
            l = []
            for x in u:
                z = (x, s[x])
                if s[x] == []: z = (x, "None")
                l.append(z)
            l = labels.make_ascii_label(l)
            s.rename(l)

    def output(self, s):
        if not s:
            return
        if not self['quiet']:
            print(s)

    def __init__(self, name, top_level=True):
        self.name = name
        self.inherited = []
        self.systems = []
        # self.added_sets = [self]
        self.added_sets = {name : self}
        self.graph_info_file = os.path.join('data', f'{self.name}.graphs')
        self.log_dir = os.path.join(Q.storage_prefix, f'logs/{self.name}')
        silent_remove(self.log_dir)
        mkdir(self.log_dir)
        self.log_file = os.path.join(self.log_dir, f'{self.name}.log')
        self.log = logging.getLogger(f'{self}')
        self.log.setLevel(logging.INFO)
        self.result_dir = os.path.join(Q.storage_prefix, f'results/{self.name}')
        f = logging.FileHandler(self.log_file, mode='w+')
        # self.log.addHandler(f)
        self.log.propagate = False
        self.log.info(f"Initialised log for {self}")
        self.input_file = os.path.join("data", self.name)
        self.input_vars = {}

        data = {}
        self.config = {}
        self.config_systems = {}
        self.process_input('default')
        self.process_input(self.name)
        self.output(f"Initialised system set {name}")

        self.retrieval_data = {}
        self.used_retrievals = set()

        for key in ['add', 'input', 'do']:
            if key in self.config_systems:
                del self.config_systems[key]
        self.config_systems['sset_name'] = [self.name]

        # *** PROFILING CODE ***
        # pr = cProfile.Profile()
        # pr.enable()
        # for x in dict_product(self.config_systems):
        #     s = System(x)
        #     break
        # pr.disable()
        # s = io.StringIO()
        # sortby = SortKey.CUMULATIVE
        # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        # ps.print_stats()
        # print(s.getvalue())
        # exit()
        #
        if self['refine']:
            r = self['refine']
            S1 = None
            while True:
                S_base = SystemSet(r)
                if S1 is None:
                    S1 = S_base
                S_base.run()
                self.inherited.append(S_base)
                r = S_base['refine']
                if not r:
                    break

            # print(S_base)
            # print(self['refine'])
            # print(self['banana'])
            v = S_base['prog_var']
            # for x, y in S_base.slicer.get_sets(exclude=S_base['prog_var']).items():
            for x, y in S1.slicer.get_sets(exclude=S_base['prog_var']).items():
                m = 100000
                imin = 0
                for i in range(len(y)):
                    v = y[i]['z1L']
                    if v < m:
                        imin = i
                        m = v

                d = S_base['prog_delta']
                d = y[imin+1]['gamma'] - y[imin]['gamma']
                x = y[imin]['gamma']
                x0 = x + d
                x1 = x - d
                points = numpy.linspace(x0, x1, int(self['refine_npoints'])+2)[1:-1]
                for p in points:
                    c = dict(y[0].config_base)
                    c['sset_name'] = self.name
                    c['gamma'] = p
                    snew = System(c, self)
                    self.systems.append(snew)
        elif self['expand']:
            def z1Le(x, A, B, C): return -A*numpy.power(x, -B)+C
            p = [1.41560118, 0.87778943, 0.56040851]
            def q(x):
                return z1Le(x, *p)

            u = ' '.join(self['expand_fun'])

            for s in dict_product(self.config_systems):
                v = eval(u)
                L = s['L']
                a = q(L)
                D = self['expand_delta']
                points = numpy.linspace(a * (1.-D), a*(1.+D), self['expand_npoints'])
                for pt in points:
                    d0 = dict(s)
                    d0['gamma'] = pt
                    S = System(d0, self)
                    self.systems.append(S)
                # print(v)
                1

        elif self['find_eps']:
            find_eps(self)
        elif self['find_eps2']:
            find_eps2(self)
        elif self['grid_eps']:
            grid_eps(self)
        else:
            # *** PROFILING CODE ***
            # pr = cProfile.Profile()
            # pr.enable()
            # for x in dict_product(self.config_systems):
            #     s = System(x)
            # pr.disable()
            # s = io.StringIO()
            # sortby = SortKey.CUMULATIVE
            # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
            # ps.print_stats()
            # print(s.getvalue())
            # exit()
            self.systems = [System(x, self) for x in dict_product(self.config_systems)]
        if self['make_systems']:
            self.systems = [System(x, self) for x in dict_product(self.config_systems)]
        if self['empty']:
            self.systems = []

        for x in self['add']:
            if not top_level:
                break
            S = SystemSet(x, top_level=False)
            self.systems += S.systems
            self.added_sets[x] = S
        for (i, subset) in enumerate(self.subsets):
            content = list(self.lines_raw) + subset
            content = '\n'.join(content)
            f = self.input_file + f'_gen{i}'
            a = os.path.basename(f)
            open(f, 'w+').write(content)
            S = SystemSet(a, top_level=False)
            self.systems += S.systems
            self.added_sets[a] = S
        if top_level:
            for s in self.systems:
                s.parent = self

        self.ref_system = None
        for s in self.systems:
            if s['reference']:
                self.ref_system = s
    def prepare(self):
        self.renew_names()
        self.output("Preparing systems")

        profile = 0
        if profile:
            pr = cProfile.Profile()
            pr.enable()
            for s in self.systems:
                s.prepare()
            if self['retrieve_from']:
                exit()
            pr.disable()
            s = io.StringIO()
            sortby = SortKey.CUMULATIVE
            ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
            ps.print_stats()
            print(s.getvalue())
            exit()

        for s in self.systems:
            s.prepare()
        if self['retrieve_from']:
            print(123)
            exit()


        self.graph_dir = os.path.join(Q.storage_prefix, os.path.join('graphs', self.name))
        mkdir(self.graph_dir)
        self.pgf_dir = os.path.join(Q.storage_prefix, os.path.join(Q.storage_prefix, os.path.join(self.graph_dir, 'pgf')))
        mkdir(self.pgf_dir)

    def get(self, key, default=None):
        return self.config.get(key, default)

    def purge_graphs(self):
        if self['purge_graphs'] == 0: return
        silent_remove(self.graph_dir)
        mkdir(self.graph_dir)
        mkdir(self.pgf_dir)

    def run(self):
        SystemSet_run.run(self)

    def fidelity(self):
        self.fidelity_file = os.path.join(Q.storage_prefix, os.path.join(self.result_dir, 'fidelity.dill'))
        self.fidelity_data = {}
        if os.path.exists(self.fidelity_file) and not Q.rerun_measurements:
            self.output("fidelity done already")
            try:
                with open(self.fidelity_file, 'rb') as f:
                    self.fidelity_data = dill.load(f)
                return
            except:
                self.output("Loading fidelity data fialed")

        self.fidelity_data["nonh"] = {}
        S = self
        xkey = S['fidelity_parameter']

        for x, y in S.slicer.get_sets(exclude=[xkey]).items():
            xvals = y.get_possible_values(xkey)
            n = len(xvals)
            zdata = numpy.zeros((n, n))
            for i in range(n):
                s1 = y[i]
                s1.drop_charges()
                s1.load_states()
                X = s1[xkey]
                for j in range(n):
                    s2 = y[j]
                    s2.drop_charges()
                    s2.load_states()
                    Y = s2[xkey]
                    psi1 = s1.chargeless_states[0]
                    psi2 = s2.chargeless_states[0]
                    if not s1['hc']:
                        psi1 = s1.make_left_state(psi1)
                    # f = psi1.overlap(psi2)
                    # psi1 = s1.left_system.states[0].psi
                    f = psi1.overlap(psi2) / psi1.overlap(psi1) / psi2.overlap(psi2)
                    zdata[i][j] = abs(f)
                    # f = s1.states[0].psi.overlap(s2.states[0].psi)
            self.fidelity_data["nonh"][x] = zdata
        with open(self.fidelity_file, 'wb+') as f:
            dill.dump(self.fidelity_data, f)


    def extras(self):
        things = self['do']
        for t in things:
            if t == 'fidelity':
                self.fidelity()
            else:
                ex = f'{t}(self)'
                exec(ex)
            # if t == 'analyse_zeros':
            #     analyse_zeros(self)
            # if t == 'periodic_lambda_test':
            #     periodic_lambda_test(self)
        1
    def find(self, data):
        return self.slicer.find(data)

    def gd(self, y, transform=None, xcoord=None, exclude=[], include=[]):
        if xcoord is None: xcoord = self['xkey']
        # exclude = exclude + [xcoord]
        res = {}
        slc = SystemSlicer([self, *self.inherited])
        exclude = exclude + ['sset_name']
        # for key, systems in self.slicer.get_sets(exclude=exclude, include=include).items():
        for key, systems in slc.get_sets(exclude=exclude, include=include).items():
            res[key] = {}
            X = []
            Y = []
            systems = sorted(systems, key=lambda x: x[xcoord])
            for s in systems:
                a = s[xcoord]
                b = s[y]
                if b is None or b == []:
                    continue
                X.append(a)
                Y.append(b)
            X = numpy.array(X)
            Y = numpy.array(Y)
            res[key]['xdata'] = X
            res[key]['ydata'] = Y
        return res

    def gds(self, ykey, transform=None, xcoord='L', exclude=[], include=[]):
        res = {}
        slc = SystemSlicer([self, *self.inherited])
        exclude = exclude + ['sset_name']
        for key, systems in slc.get_sets(exclude=exclude, include=include).items():
        # for key, systems in self.slicer.get_sets(exclude=exclude, include=include).items():
            res[key] = {}
            systems = sorted(systems, key=lambda x: x[xcoord])
            X = []
            Y = []
            s = systems[0]
            for i in range(s['L']):
                k = f'{ykey}{i}'
                y = s[f'{ykey}{i}']
                if y == []:
                    continue
                X.append(i)
                Y.append(y)
            res[key]['xdata'] = numpy.array(X)
            res[key]['ydata'] = numpy.array(Y)
        return res

if __name__ == "__main__":
    # A = SystemSet("test1")
    A = SystemSet("test_slow")

    # sets = A.slicer.get_sets()
    A.run()
