#!/usr/bin/env python3
from rah_utility import newliner

def report(S):
    res = newliner()
    n = len(S.systems)
    res += f'System set {S} length {len(S.systems)}'
    runs = 0
    for s in S.systems:
        if s.will_run:
            runs += 1
    res += f'{runs} / {n} will run'
    # print(res)
    S.log.info(res)
    return res
