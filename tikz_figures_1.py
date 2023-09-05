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

class FP_Point:
    def __init__(self, z, parent, name, size):
        self.z = z
        self.parent = parent
        self.name = name
        self.size = size
        self.coord = f'({z.real}, {z.imag})'
        self.q = [int(x) for x in self.name.split('_')[1:]]

    def plus_coord(self, other):
        z = self.z + other.z
        return f'({z.real}, {z.imag})'

def fp_coords4():
    res = newliner()
    # how deep
    depth = 4
    fsize=12

    # define epsilon colours
    #
    epsilon_colours = []
    cmap = matplotlib.cm.rainbow
    cmap = matplotlib.cm.plasma
    cmap = matplotlib.cm.brg
    cmap = matplotlib.cm.inferno
    light = 0.8
    offset_left = 0.25
    offset_right = 0.15
    reverse = 0
    for d in range(depth):
        if reverse:
            d = depth - d - 1
        cind = offset_left + (1.-offset_left-offset_right) * (d/(depth-1))
        c = cmap(cind)
        # c = cmap((depth-d-1) / depth)
        c = [x * light for x in c]
        cstring = f'{{{c[0]},{c[1]},{c[2]}}}'
        colname = f'col{d}'
        res += r'\definecolor{'+colname+r'}{rgb}'+cstring
        epsilon_colours.append(colname)

    # epsilon line styles
    use_patterns = 0
    epsilon_patterns = [
        'on 11pt off 3pt',
        'on 7pt off 2pt',
        'on 5pt off 1.5pt',
        'on 2pt off 0.8pt',
        'on 2pt off 0.8pt',
        ]
    Epointsize=5
    base_length = 3.5
    line_lengths = [1, 0.5, 0.25, 0.12, 0.05]
    line_lengths = [base_length * x for x in line_lengths]
    line_widths= [2.5, 2.25, 2.0, 1.8, 1.5]
    if not use_patterns:
        line_widths= [3.0, 2.50, 2.0, 1.5, 1.2]

    epsilon_linestyles = []
    for d in range(depth):
        s = fr'dash pattern={epsilon_patterns[d]}, color={epsilon_colours[d]}, line width={line_widths[d]}pt'
        if not use_patterns:
            s = fr'color={epsilon_colours[d]}, line width={line_widths[d]}pt'
        epsilon_linestyles.append(s)

    # legend
    legend_scale = 1.2

    # label line pattern
    pattern_2 = 'on 2pt off 2pt'
    line_width_label = 1
    # label line params
    lblank = 0.15
    labels = {
        '0_0_0_0_0' : [+0.2, 2.0, 90],
        '0_1_0_0_1' : [-0.35, 1.5, -90],
        # '0_1_0_0_1' : [-0.40, 1.2, -90],
        '0_2_0_0_2' : [+0.30, 1.5, 90],
        }
    # labels = {
    #     '0_0_0_0_0_0' : [+0.2, 2.0, 90],
    #     '0_0_1_1_1_0' : [-0.30, 1.5, -90],
    #     '0_2_0_0_2_0' : [+0.30, 1.5, 90],
    #     }
    # unused?
    base_size = 30
    size_ratio = 0.6


    N = 3


    points = {
        -1 : [],
        0 : [FP_Point(0.j, None, '0', base_size)],
    }
    for d in range(1, depth+1):
        old_points = points[d-1]
        points[d] = []
        for p in old_points:
            for i in range(N):
                # znew = p.z - length * (ratio ** d) * numpy.exp(2.j * numpy.pi / N * i)
                znew = p.z - line_lengths[d-1] * numpy.exp(2.j * numpy.pi / N * i)
                pnew = FP_Point(znew, p, f'{p.name}_{i}', p.size * size_ratio)
                points[d].append(pnew)

    pltargs = {'linestyle' : '',
               'marker' : 'o',
               }
    # for d, ps in points.items():
    #     print(d, len(ps))
    #     for p in ps:
    #         plt.plot([p.real], [p.imag], **pltargs)
    # plt.show()

    res += r'\begin{tikzpicture}'
    for d in points:
        ps = points[d]
        ind = 0
        for p in ps:
            x, y = p.z.real, p.z.imag
            # res += fr'\node[circle,fill=blue!20,minimum size={p.size}pt, inner sep=0pt] ({p.name}) at ({x}, {y}) {{}};'
            ind += 1

            if p.parent:
                # res += fr'\draw ({p.name}) -- ({p.parent.name});'
                # res += fr'\draw[dash pattern={pattern_1}, color={colname}, line width={line_widths[d-1]}pt] {p.coord} -- {p.parent.coord};'
                res += fr'\draw[{epsilon_linestyles[d-1]}] {p.coord} -- {p.parent.coord};'


    for p in points[depth]:
        res += fr'\node[circle,fill=black,minimum size={Epointsize}pt, inner sep=0pt] ({p.name}) at {p.coord} {{}};'
        if p.name in labels:
            # res += fr'\node[circle,fill=red,minimum size={p.size}pt, inner sep=0pt] ({p.name}) at {p.coord} {{}};'
            a, l, labpos = labels[p.name]
            zp = p.z + numpy.exp(a*2.j*numpy.pi) * l
            zp1 = p.z + lblank * (zp - p.z) / abs(zp-p.z)

            lab = ''
            for i in range(depth):
                x = p.q[i]
                omg = fr'\omega^{x}'
                if x == 0:
                    omg = ''
                if x == 1:
                    omg = '\omega'
                lab += fr'-{omg}\epsilon_{i+1}'
            lab = fr'${lab.strip("+")}$'
            res += fr'\draw[dash pattern={pattern_2}, line width={line_width_label}pt] {make_pos(zp1)} -- {make_pos(zp)};'
            nodestyle = fr'circle,fill=red!0,minimum size={p.size}pt, inner sep=0pt'
            nodestyle += fr',label={{[scale=1.0]{labpos}:{lab}}}'
            res += fr'\node [{nodestyle}] () at {make_pos(zp)} {{}};'


    # LEGEND
    legend_style = 1
    if legend_style == 1:
        q = 0
        for d in range(depth):
            q += line_lengths[d]
        x, y = q/2. + 1.5, q/2.*numpy.sqrt(3)
        # res += fr'\node[circle,fill=red,minimum size=25pt, inner sep=0pt] () at ({x}, {y}) {{}};'

        base = (x, y)
        x = base[0]
        y = base[1]
        for d in range(depth):
            l = line_lengths[d]
            ps = points[d]
            p = ps[0]
            c = cmap(d / (len(points)-1))
            c = [x * light for x in c]
            colname = f'col{d+1}'
            res += fr'\node[circle, anchor=east, fill=red!0, minimum size=15pt, inner sep=0pt, scale={legend_scale}] (N{d}) at ({x}, {y}) {{$\epsilon_{d+1}$}};'
            res += fr'\draw[{epsilon_linestyles[d]}] ({x}, {y}) -- ({x+l}, {y});'
            # res += fr'\draw[dash pattern={pattern_1}, color={colname}, line width={line_widths[d]}pt] ({x}, {y}) -- ({x+l}, {y});'
            # res += fr'\draw[dash pattern=on 3pt off 1pt, color={colname}, line width=1.2pt] (N{d}) -- +({l}, 0);'
            y -= 1
    else:
        q = 0
        for d in range(depth):
            q += line_lengths[d]
        x, y = -q, -q/2.*numpy.sqrt(3) +0.1
        # res += fr'\node[circle,fill=red,minimum size=25pt, inner sep=0pt] () at ({x}, {y}) {{}};'

        base = (x, y)
        x = base[0]
        y = base[1]
        for d in range(depth):
            l = line_lengths[d]
            ps = points[d]
            p = ps[0]
            c = cmap(d / (len(points)-1))
            c = [x * light for x in c]
            colname = f'col{d+1}'
            res += fr'\node[circle, anchor=east, fill=red!0, minimum size=15pt, inner sep=0pt, scale={legend_scale}] (N{d}) at ({x}, {y}) {{$\epsilon_{d+1}$}};'
            res += fr'\draw[{epsilon_linestyles[d]}] ({x}, {y}) -- ({x+l}, {y});'
            # res += fr'\draw[dash pattern={pattern_1}, color={colname}, line width={line_widths[d]}pt] ({x}, {y}) -- ({x+l}, {y});'
            # res += fr'\draw[dash pattern=on 3pt off 1pt, color={colname}, line width=1.2pt] (N{d}) -- +({l}, 0);'
            y += 1


    # center/origin
    res += r'\node[circle,fill=white, minimum size=13pt, inner sep=0pt, scale=1.5] () at (0, 0) {$0$};'

    res += r'\end{tikzpicture}'
    res = str(res)
    f = '/home/rah/s/phd2/tikz_spectrum4.tex'
    f = '/home/rah/s/paper_complex_lambda/tikz_spectrum4.tex'
    f = '/home/rah/s/paper_cl/tikz_spectrum4.tex'
    print(f)
    # print(res)
    open(f, 'w+').write(res)

