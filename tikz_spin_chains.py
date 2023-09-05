#!/usr/bin/env python3

import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import os
import matplotlib
import numpy
from rah_utility import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import random

def make_pos(z):
    return f'({z.real}, {z.imag})'

def tf_diag1():
    Npoints = 4
    len1 = 1.8
    isep = 0.5
    res = newliner()
    res += r'\usetikzlibrary{backgrounds}'
    # res += r'\begin{tikzpicture}[background rectangle/.style={fill=gray!10}, show background rectangle]'
    res += r'\begin{tikzpicture}'
    res += r'\usetikzlibrary {arrows.meta}'
    lw = '1.0pt'
    w1 = '0.80cm'
    L = 10

    dx = 1.4
    dy = 1.4

    # atoms
    nodestyle1 = fr'circle, fill=red!49, minimum size=0.9cm'

    # arrows
    alen = 0.25
    alen2 = 0.6
    A1 = fr'{{Latex[length=2.0mm]}}'
    A2 = fr'line width=1pt,'

    ketstr = ''
    for i in range(L):
        atype = i % 2
        atype = random.choice([0, 1])
        z = 0.j + i * dx
        res += fr'\node[{nodestyle1}] (A{i}) at {make_pos(z)} {{}};'
        if atype:
            z1 = z - alen * 1.j * alen2
            z2 = z + alen * 1.j
            res += fr'\draw [{A2}, -{A1}] {make_pos(z1)} -- {make_pos(z2)};'
        else:
            z1 = z - alen * 1.j
            z2 = z + alen * 1.j * alen2
            res += fr'\draw [{A2}, {A1}-] {make_pos(z1)} -- {make_pos(z2)};'

        z = 0.j + i * dx - 1.j * dy + 0.1
        res += fr'\node[inner sep=0pt] (B_{i}) at {make_pos(z)} {{\Large$|{atype}\rangle_{i}$}};'

        if i < L-1:
            z = 0.j + (i+0.5) * dx - 1.j * dy + 0.1
            res += fr'\node[inner sep=0pt] (B_{i}) at {make_pos(z)} {{$\otimes$}};'

        if i > 0:
            res += fr'\draw (A{i-1}) -- (A{i});'

        ketstr += f'{atype}'

    z = 0.j - 1.j * dy - 1.0
    res += fr'\node[] (C) at {make_pos(z)} {{\Large$\sim$}};'
    z = 0.j - 2.j * dy + 5.1
    res += fr'\node[] (C) at {make_pos(z)} {{\Huge$\sim\;$\Huge $|{ketstr}\rangle$}};'
    f = '/home/rah/s/phd2/tikz_spin_chains.tex'
    res += r'\end{tikzpicture}'
    res = str(res)
    open(f, 'w+').write(res)

def build():
    c = ['plx', 'tikz_gen_test']
    p = proc(' '.join(c), cwd='/home/rah/s/phd2', shell=True)
    print(p)

tf_diag1()
build()
