#!/usr/bin/env python3

import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import os
import matplotlib
import numpy
from rah_utility import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def make_pos(z):
    return f'({z.real}, {z.imag})'

def tf_diag1():
    Npoints = 4
    len1 = 1.8
    isep = 0.5
    res = newliner()
    res += r'\begin{tikzpicture}'
    lw = '1.0pt'
    w1 = '0.80cm'

    nodestyle = fr'[draw,circle,fill=blue!0,inner sep={isep}, minimum width={w1}]'
    linestyle = fr'[line width={lw}]'
    for j in range(Npoints):
        x = j * len1
        y = -j * len1
        res += fr'\node{nodestyle} (j{j}) at ({x}, {y}) {{$j_{j+1}$}};'
        x = len1 * j
        y = -len1 * (j+1)
        res += fr"\node{nodestyle} (jp{j}) at ({x}, {y}) {{$j'_{j+1}$}};"
    for j in range(Npoints):
        res += fr'\draw{linestyle} (j{j}) -- (jp{j});'
    for j in range(Npoints-1):
        res += fr'\draw{linestyle} (jp{j}) -- (j{j+1});'
    # res += fr'\draw(jp{Npoints-1}) arc[radius={Npoints*len1}, start angle=0, delta angle=90] -- (j0);'
    # a = len1
    # z1 = a + 1.j * a + a * (2.j)
    # z2 = Npoints * a - Npoints*a * 1.j + a * (2. + 2.j)
    # res += fr'''\draw{linestyle}(jp{Npoints-1})
    # .. controls {make_pos(z2)} and {make_pos(z1)} ..
    # (j0);'''
    l2 = len1 * 0.7
    lstyle2 = fr'[line width={lw}, dash pattern=on 2pt off 1pt]'
    res += fr'\draw{lstyle2}(j0) -- ++ ({-l2}, 0);'
    res += fr'\draw{lstyle2}(jp{Npoints-1}) -- ++ ({l2}, 0);'


    f = '/home/rah/s/phd2/tikz_transfer1.tex'
    res += r'\end{tikzpicture}'
    res = str(res)
    open(f, 'w+').write(res)

def build():
    c = ['plx', 'tikz_transfer']
    p = proc(' '.join(c), cwd='/home/rah/s/phd2', shell=True)
    print(p)

tf_diag1()
build()
