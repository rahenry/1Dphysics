#!/usr/bin/env python3

import numpy
from rah_utility import dict_product, is_numeric, rel_diff

latex_label_data = {
    'phi' : r'\phi',
    'gamma' : r'\gamma',
}
def is_float(d):
    return isinstance(d, float) or isinstance(d, numpy.floating)

def make_ascii_label(data):
    res = r''
    if isinstance(data, dict):
        x = frozenset(data.items())
    data = sorted(data)
    for x, y in data:
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

def is_float(d):
    return isinstance(d, float) or isinstance(d, numpy.floating)

def make_latex_label(x):
    res = r''
    if isinstance(x, dict):
        x = frozenset(x.items())
    x = sorted(x)
    for a, b in x:
        if a in latex_label_data:
            a = latex_label_data[a]
        # res += rf'{a} = {b}, '
        if is_float(b):
            b = f'{b:.3f}'
        res += rf'${a} = {b}$\\ '

# \pgftext{\begin{tabular}{@{}l@{}}can haz line\\break plz?\end{tabular}}

    res = res[:-3]
    # res = r'\begin{minipage}{1cm}$' + res + r'$\end{minipage}'
    # res = r'\begin{minipage}{5cm}can haz line \\ break plz?\end{minipage}'
    res = r'\begin{tabular}{@{}l@{}}' + res + '\end{tabular}'
    return res
