#!/usr/bin/env python3
from args import Q
from report import report
from System import System
import multiprocess as mp
import time
import cProfile, pstats, io
from pstats import SortKey
import numpy

def proggy(i, N):
    q = 10
    def f(u):
        return numpy.floor(q*u/N)
    if not f(i) == f(i+1):
        return f'Progress = {f(i)/q}'

def run_prog_parallel(S):
    for s in S.systems:
        s.run_test()

    def runfun(s):
        r = 0
        if not s.run_all_finished:
            r = 1
            S.log.info(f'{s} started')
        s.run()
        s.measure()
        if r:
            S.log.info(f'{s} finished')
        return 1

    active = []
    finished = []
    ind = 0
    for s in S.systems:
        s.load_measurements()


    for a in S.systems:
        s = a
        for x in a.auxillary_systems:
            s = x
            for y in x.auxillary_systems:
                s = y
    # pr = cProfile.Profile()
    # pr.enable()
    while True:
        if ind == len(S.systems):
            break
        s = S.systems[ind]
        if not s.run_all_finished:
            ind += 1
            continue
        if s.prog_tested:
            ind += 1
            continue
        s.prog_tested = 1
        if not is_prog_finished(S, s):
            add_prog_system(S, s)
            ind = 0
        else:
            ind += 1
    # exit()
    # pr.disable()
    # s = io.StringIO()
    # sortby = SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # exit()
    to_run = [s for s in S.systems if not s.run_all_finished]
    S.log.info(f'{len(to_run)} base systems will run')

    mold = ''
    while True:
        while len(active) < Q.nprocs and len(to_run) > 0:
            s = to_run.pop()
            P = mp.Process(target=runfun, args=(s,))
            P.start()
            active.append((P, s))
        ind = 0
        while True:
            if ind == len(active):
                break
            P, s = active[ind]
            if not P.is_alive():
                del active[ind]
                s.load_measurements()
                finished.append((P, s))
                ind = 0
                if not is_prog_finished(S, s):
                    snew = add_prog_system(S, s)
                    if not snew.run_all_finished:
                        to_run.append(snew)
            else:
                ind += 1

        if len(to_run) == 0 and len(active) == 0:
            S.log.info('all finished')
            break

        m = f'{len(active)} procs active, {len(finished)} procs finished'
        if not m == mold:
            mold = m
            S.log.info(f'{time.time()}: {len(active)} procs active, {len(finished)} procs finished')
        time.sleep(3)

    for f in finished:
        S.log.info(f)
def run_standard(S):
    self = S
    self.log.info(f'Running systemset {self}')
    if Q.nprocs > 1:
        self.log.info('Submitting run')
        def runfun(s):
            if s.will_run:
                self.log.info(f'{s} started')
            s.run()
            s.measure()
            if s.will_run:
                self.log.info(f'{s} finished')
            return 1

        for s in self.systems:
            s.run_test()
        report(self)
        # with mp.get_context('spawn').Pool(Q.nprocs) as pool:
        #     pool.map(runfun, self.systems)
        active = []
        finished = []
        ind = 0
        to_run = [s for s in self.systems if s.will_run or s.will_measure]
        N = len(to_run)
        self.log.info(f'{N}/{len(self.systems)} systems will run')
        while (True):
            while len(active) < Q.nprocs:
                if ind >= len(to_run):
                    break
                s = to_run[ind]
                P = mp.Process(target=runfun, args=(s,))
                P.start()
                active.append((ind, P, s))
                ind += 1

            self.log.info(f'{len(active)} procs active, {len(finished)} procs finished, ind={ind}/{len(to_run)}')

            to_delete = []
            for j in range(len(active)):
                i1, P1, s1 = active[j]
                if not P1.is_alive():
                    to_delete.append(j)
                    finished.append(active[j])
            to_delete.sort()
            to_delete.reverse()
            for d in to_delete:
                del active[d]

            if ind == N and len(active) == 0:
                self.log.info("Finished!")
                break

            time.sleep(3)


        for f in finished:
            self.log.info(f)

        # with pathos.pools._ProcessPool(Q.nprocs) as pool:
        #     pool.map(runfun, self.systems, chunksize=1)
        for s in self.systems:
            s.load_measurements()
    else:
        for s in self.systems:
            s.run_test()
        S.output(report(self))
        n = len(self.systems)
        i = 0
        for s in self.systems:
            if not s.will_run:
                continue
            if n < 50:
                S.output(f'Running {s}')
            S.output(proggy(i, len(S.systems)))
            i += 1
            s.run()
            if not S['delay_measurements']:
                s.measure()
        if S['delay_measurements']:
            i = 0
            for s in self.systems:
                if n < 50:
                    S.output(f'Measuring {s}')
                S.output(proggy(i, len(S.systems)))
                i += 1
                s.measure()
        for s in self.systems:
            s.load_measurements()

def add_prog_system(S, s):
    prog_var = S['prog_var']
    prog_type = S['prog_type']
    c = dict(s.config_base)
    p = s[prog_var]
    p = p * S['prog_mult'] + S['prog_delta']
    if prog_type == 'int':
        p = int(p)
    S.log.info(f'prog_test failed; creating a new system from {s} with {prog_var} = {p}')
    c[prog_var] = p
    snew = System(c)
    s.next = snew
    snew.previous= s
    S.systems.append(snew)
    S.renew_names()
    return snew


def is_prog_finished(S, s):
    c = S['prog_cond']
    prog_var = S['prog_var']
    p = s[prog_var]
    if s[prog_var] >= S['prog_max']:
        S.log.info(f'system {s} has {p} exceeding max; stopping chain')
        return True
    if s[prog_var] <= S['prog_min']:
        S.log.info(f'system {s} has {p} less than min; stopping chain')
        return True
    if c == 'fsL1':
        # if not s.previous:
        #     return False
        # x1 = s.previous['fsL']
        x2 = s['fsL']
        if x2 is None:
            S.log.info(f"system {s} hasn't been measured")
            S.log.info(f"{s['e_0']}, {s['e_exact']}")
            return False

        if x2 < 0:
            S.log.info(f'system {s} has fsL = {x2}; stopping chain')
            return True
        else:
            S.log.info(f'system {s} has fsL = {x2}; continuing')
    elif c == 'z1L':
        if not s.previous:
            return False
        if not s.previous.previous:
            return False
        x0 = s.previous.previous['z1L']
        if not x0:
            S.log.info(f"system {S} has no x0")
            return False
        x1 = s.previous['z1L']
        x2 = s['z1L']
        S.log.info(f"system {s} has {x0}, {x1}, {x2}")
        if x2 is None:
            S.log.info(f"{s['e_0']}, {s['e_exact']}")
            S.log.info(f"system {s} hasn't been measured")
            return False

        if x2 > x1 and x2 > x0:
            S.log.info(f'system {s} has z1L = {x2}; x0={x0}, x1={x1}; stopping chain')
            return True
        else:
            S.log.info(f'system {s} has z1L = {x2}; continuing')
        1


def run_prog(S):
    prog_var = S['prog_var']
    prog_type = S['prog_type']
    for s in S.systems:
        s.run_test()
    to_run = S.systems
    while (True):
        any_failed = 0
        new_systems = []
        for s in to_run:
            s.run()
            s.measure()
            if not is_prog_finished(S, s):
                c = dict(s.config_base)
                p = s[prog_var]
                p = p * S['prog_mult'] + S['prog_delta']
                if prog_type == 'int':
                    p = int(p)
                print(f'prog_test failed; creating a new system from {s} with {prog_var} = {p}')
                c[prog_var] = p
                snew = System(c)
                s.next = snew
                snew.pervious = s
                snew.run_test()
                new_systems.append(snew)
                any_failed = 1
        if any_failed:
            to_run = new_systems
            S.systems += new_systems
            S.renew_names()
        else:
            break


    1

def run(S):
    pr = None
    if Q.profile_run:
        print("Profiling")
        pr = cProfile.Profile()
        pr.enable()
    m = S['run_mode']
    if not m:
        m = 'standard'
    S.output(f'run_mode = {m}')
    if m == 'standard':
        run_standard(S)
    if m == 'prog':
        if Q.nprocs > 1:
            run_prog_parallel(S)
        else:
            run_prog(S)

    if Q.profile_run:
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        z = s.getvalue()
        open("profile_data", "w+").write(z)
