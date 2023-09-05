#!/usr/bin/env python3

from System import System
from rah_utility import dict_product, is_numeric, rel_diff
import numpy
def is_float(d):
    return isinstance(d, float) or isinstance(d, numpy.floating)

class SystemSlicer:
    def __init__(self, ssets):
        if not isinstance(ssets, list):
            ssets = [ssets]
        if ssets and isinstance(ssets[0], System):
            self.systems = ssets
        else:
            self.ssets = ssets
            self.systems = []
            for S in ssets:
                self.systems += S.systems
        self.float_width = 14
        self.sep_width = 2
        self.sep = ' ' * self.sep_width
        all_identifiers = set()
        for s in self.systems:
            all_identifiers.update([x for x in s.config_base])

        all_identifiers = sorted(list(all_identifiers))
        s0 = self.systems[0]
        R = {}
        for x in all_identifiers:
            z = s0.get_base(x)
            R[x] = set()
            for s in self.systems:
                z = s.get_base(x)
                if z is None: continue
                if z == []: continue
                R[x].add(z)
            if R[x] == set():
                del R[x]
        self.all_identifiers = R

        for x, y in self.all_identifiers.items():
            y = list(y)
            y.sort()
            self.all_identifiers[x] = y
        S = {}
        for x, data in R.items():
            if len(data) == 1:
                continue
            S[x] = data
        self.relevant_identifiers = S

        W = {}
        for x, data in R.items():
            W[x] = len(str(x))
            for d in data:
                z = len(str(d))
                if is_float(d):
                    z = self.float_width
                if (z > W[x]):
                    W[x] = z

        self.width_info = W

    def find(self, data):
        if isinstance(data, frozenset):
            data = dict(data)
        res = []
        for s in self.systems:
            match = True
            for x, y in data.items():
                z = s.get(x, None)
                if z is None or z is []: continue
                tol = 1E-10
                if is_numeric(z):
                    if not rel_diff(z, y) < tol:
                        match = False
                        break
                elif not z == y:
                    match = False
                    break
            if match:
                res.append(s)
        if not res:
            return []
        res = SystemSlicer(res)
        return res

    def get_sets(self, exclude=[], include=[]):
        combs = {}
        if not include == []:
            for inc in include:
                combs[inc] = self.relevant_identifiers[inc]
        else:
            for x, data in self.relevant_identifiers.items():
                if x in exclude:
                    continue
                combs[x] = data
        combs = dict_product(combs)
        res = {}
        for c in combs:
            matches = self.find(c)
            z = frozenset(c.items())
            if len(matches) == 0:
                continue
            res[z] = matches
        return res

    def print_ident(self, data):
        res = ''
        for x, y in data.items():
            if is_float(y):
                if abs(y) < 1E-4:
                    res += f'{x}={y}, '
                else:
                    q = f'{y:2.8f}'
                    res += f'{x}={q}, '
            elif is_numeric(y):
                res += f'{x}={y}, '
            else:
                res += f'{x}={y}, '
            continue
            if is_float(y):
                q = f'{y:010f}'
                res += f'{x}={q}, '
                # if abs(y) > 1E-6:
                #     q = f'{y:010f}'
                #     res += f'{x}={q}, '
                # else:
                #     res += f'{x}={y}, '
            elif is_numeric(y):
                q = f'{y:10}'
                q = q.replace(' ', '0')
                res += f'{x}={q}, '
            else:
                res += f'{x}={y}, '
        if res == '':
            res = 'all'
        res = res.strip()
        res = res.strip(',')
        return res

    def printify_short(self, s):
        res = ''
        for x in self.relevant_identifiers:
            res += f'{x}={s[x]},'
        res = res.strip(',')
        return res

    def printify(self, s):
        res = f'System {s["model"]} {s["method"]}, '
        for x in self.relevant_identifiers:
            res += f'{x}={s[x]}, '
        res = res.strip()
        res = res.strip(',')
        return res

    def tabulate(self, targets=[], relabel=None, specify_widths={}):
        header = ''
        # for x in self.relevant_identifiers:
        #     w = self.width_info[x]
        #     header += f'{x:>{w}}' + self.sep

        labels = {}
        i = 0
        for t in targets:
            if relabel and i < len(relabel):
                labels[t] = relabel[i]
            else:
                labels[t] = t
            i += 1

        widths = {}
        for t in targets:
            wmax = 0
            for s in self.systems:
                w = str(s[t])
                if isinstance(s[t], float):
                    w = w.rstrip('0')
                w = len(w)
                if t in specify_widths:
                    w = specify_widths[t]
                if w > wmax:
                    wmax = w
            if wmax > self.float_width: wmax = self.float_width
            if len(labels[t]) > wmax:
                wmax = len(labels[t])
            widths[t] = wmax

        for t in targets:
            w = widths[t]
            header += f'{labels[t]:>{w}}' + self.sep
        lines = [header]

        for s in self.systems:
            line = ''
            # for x in self.relevant_identifiers:
            #     y = s[x]
            #     w = self.width_info[x]
            #     form = f'>{w}'
            #     if y is None:
            #         y = "n/a"
            #         form = f'>{w}s'
            #     elif is_float(y):
            #         form = f'>2.{w-2}f'
            #     line += f'{y:{form}}' + self.sep
            for x in targets:
                y = s[x]
                # w = max(self.float_width, len(str(x)))
                w = widths[x]
                if x in specify_widths:
                    w = specify_widths[x]
                q = ''
                if y:
                    if isinstance(y, str):
                        form = f'>{w}'
                        q = f'{s[x]:{form}}'
                    elif isinstance(y, int):
                        form = f'>{w}'
                        q = f'{s[x]:{form}}'
                    else:
                        z = w-3
                        if z < 0:
                            z = 0
                        if z == 0:
                            z = 1
                        form = f'>{w}.{z}f'
                        q = f'{s[x]:{form}}'
                else:
                    form = f'>{w}'
                    q = f'{str(y):{form}}'
                # line += q
                line += f'{q:>{widths[x]}}' + self.sep
            lines.append(line)

        res = '\n'.join(lines)
        return res

    def get_possible_values(self, coord, filter=None):
        res = set()
        for s in self.systems:
            match = True
            if filter:
                for x, y in filter.items():
                    if not s[x]==y:
                        match = False
                        break
            if not match: continue
            res.add(s[coord])
        return sorted(list(res))

    def __getitem__(self, x):
        return self.systems[x]

    def __iter__(self):
        self.position = 0
        return self

    def __len__(self):
        return len(self.systems)

    def __next__(self):
        if self.position < len(self.systems):
            res = self.systems[self.position]
            self.position += 1
            return res
        else:
            raise StopIteration

