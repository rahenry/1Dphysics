#!/usr/bin/env python3

import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import os
import matplotlib
import numpy
from rah_utility import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.colors as mcolors

class FF_Point:
    def __init__(self, z, parent, q):
        self.z = z
        self.parent = parent
        if parent:
            self.name = parent.name + f'_{q}'
        else:
            self.name = '0'

def make_pos(z):
    return f'({z.real}, {z.imag})'

def make_tikz_col(name, c):
    c = mcolors.CSS4_COLORS[c]

    c = mcolors.to_rgb(c)
    cstring = f'{{{c[0]},{c[1]},{c[2]}}}'
    res = r'\definecolor{'+name+r'}{rgb}'+cstring
    return res

def tf_diag1():
    res = newliner()
    res += r'\begin{tikzpicture}'

    # res += fr'\node[circle,fill=black] (a) at (0, 0) {{}};'
    epsilons = [2, 1.5, 1.19908, 0.83234, 0.455, 0.3123,]
    height = 0.8
    linescale = 1.
    lwidth=0.8
    pattern = ''
    pattern = 'on 3pt off 1pt'

    size1=3.5
    size2=5
    c1 = 'royalblue'
    c2 = 'firebrick'
    res += make_tikz_col('c1', c1)
    res += make_tikz_col('c2', c2)
    nodestyle1 = fr'[circle,fill opacity=0.7,fill=c1,minimum size={size1}pt,inner sep=0pt]'
    nodestyle2 = fr'[circle,fill opacity=0.7,fill=c2,minimum size={size2}pt,inner sep=0pt]'
    linestyle = fr'[dash pattern={pattern}, line width={lwidth}pt,color=gray]'

    axiswidth = sum(epsilons)*1.1
    axisheight = 0.7
    zaxis = 0. - 1.j * len(epsilons) * height - 1.j * axisheight
    arrow_axis = '{-Latex[scale=1.5]}'
    res += fr'\draw[arrows={arrow_axis}] {make_pos(zaxis-axiswidth)} -- {make_pos(zaxis+axiswidth)};'
    l1 = -0.4
    zlabel = zaxis + l1 * 1.j
    res += fr'\draw {make_pos(zaxis)} -- {make_pos(zlabel)};'
    res += fr'\node[circle,fill=black!0, inner sep=1pt] (zero) at {make_pos(zlabel)} {{0}};'
    zlabel = zaxis + axiswidth + l1 * 1.j
    res += fr'\node[circle,fill=black!0, inner sep=0pt] (zero) at {make_pos(zlabel)} {{$E$}};'


    # epsilon labels
    epsilon_label_height = 0.2
    linestyle_epsilon_label = '[<->]'
    linestyle_epsilon_label = '[Straight Barb-Straight Barb,color=black!80]'

    points = {0 : [FF_Point(0.j, None, None)]}
    depth = 0
    for ep in epsilons:
        depth += 1
        points[depth] = []
        for p in points[depth-1]:
            for q in [-1, 1]:
                znew = p.z - linescale*(1.j*height + q*ep)
                new_point = FF_Point(znew, p, q)
                points[depth].append(new_point)

    for depth in points:
        for p in points[depth]:
            if depth != len(epsilons):
                # res += fr'\coordinate ({p.name}) at {make_pos(p.z)} {{}};'
                res += fr'\node{nodestyle1} ({p.name}) at {make_pos(p.z)} {{}};'
            else:
                res += fr'\node{nodestyle2} ({p.name}) at {make_pos(p.z)} {{}};'
            if p.parent:
                res += fr'\draw{linestyle} ({p.name}) -- ({p.parent.name});'

    for depth in points:
        if depth == 0: continue
        h = epsilon_label_height
        p1 = points[depth-1][0]
        p2 = points[depth][0]
        print(p1.name, p2.name)
        z1 = p1.z + h * 1.j
        z2 = p2.z.real + p1.z.imag * 1.j + h * 1.j
        zlabel = 0.5 * (z1.real + z2.real) + 1.j * z1.imag
        res += fr'\draw {linestyle_epsilon_label} {make_pos(z1)} -- {make_pos(z2)};'
        res += fr'\node [label=above:{{$\epsilon_{depth}$}}] ({depth}_epslabel) at {make_pos(zlabel)} {{}};'

    f = '/home/rah/s/phd2/tikz_free_fermion_1.tex'
    res += r'\end{tikzpicture}'
    res = str(res)
    print(res)
    open(f, 'w+').write(res)

def build():
    c = ['plx', 'tikz_gen_test']
    p = proc(' '.join(c), cwd='/home/rah/s/phd2', shell=True)
    print(p)

tf_diag1()
build()
