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

def fp_coords():
    depth = 4
    base_size = 30
    size_ratio = 0.6
    points = {
        -1 : [],
        0 : [FP_Point(0.j, None, '0', base_size)],
    }
    length = 5
    ratio = 0.5
    N = 3
    for d in range(1, depth+1):
        old_points = points[d-1]
        points[d] = []
        for p in old_points:
            for i in range(N):
                znew = p.z - length * (ratio ** d) * numpy.exp(2.j * numpy.pi / N * i)
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

    res = newliner()
    res += r'\begin{tikzpicture}'
    for d, ps in points.items():
        ind = 0
        for p in ps:
            x, y = p.z.real, p.z.imag
            res += fr'\node[circle,fill=blue!20,minimum size={p.size}pt, inner sep=0pt] ({p.name}) at ({x}, {y}) {{}};'
            ind += 1

            if p.parent:
                res += fr'\draw ({p.name}) -- ({p.parent.name});'

    res += r'\end{tikzpicture}'
    res = str(res)
    f = '/home/rah/s/phd2/tikz_review1.tex'
    # print(res)
    open(f, 'w+').write(res)

    c = ['plx', 'tikz_review']
    p = proc(' '.join(c), cwd='/home/rah/s/phd2', shell=True)
    # print(p)

def fp_coords2():
    res = newliner()
    # how deep
    depth = 4
    fsize=12

    # define epsilon colours
    epsilon_colours = [
        'red',
        'blue',
        'green',
        'orange',
        ]

    epsilon_colours = []
    cmap = matplotlib.cm.rainbow
    cmap = matplotlib.cm.plasma
    light = 0.8
    offset_left = 0.0
    offset_right = 0.3
    for d in range(depth):
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
    Epointsize=6
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
    f = '/home/rah/s/phd2/tikz_spectrum.tex'
    # print(res)
    open(f, 'w+').write(res)

def fp_arms1():
    ds1 = 1.
    ds2 = 1.8
    Ns = [2, 3, 4]
    l1 = 1.
    l2 = 2.
    r1 = 0.5 * l1
    node_distances = [l1, l1, l2, l1]
    line_styles = [0, 0, 1, 0]
    nodes_filled = {
        2 : [0, 0, 1, 0],
        3 : [1, 0, 0, 2],
        4 : [0, 1, 0, 0],
        }
    end_names = {
        2 : ["$1$", "$\omega$"],
        3 : ["$1$", "$\omega$", "$\omega^2$"],
        4 : ["$1$", "$\omega$", "$\omega^2$", "$\omega^3$"],
        }
    node_labels = [
        "$1$", '$2$', '$L\hspace{1cm}-2$', '$L$',
        ]

    lend = 0.8 * l1
    siz1=10
    siz0=5

    for N in Ns:
        res = newliner()
        res += r'\begin{tikzpicture}'
        res2 = newliner()
        res3 = newliner()

        nodestyle = f'circle, fill=black, draw=black, minimum size={siz0}pt, inner sep=0pt'
        res2 += fr'\node[{nodestyle}] () at (0, 0) {{}};'

        for i in range(N):
            zp = 0.j
            posp = '(0, 0)'
            a = numpy.exp(2.j * numpy.pi * i / N)
            for j in range(len(node_distances)):
                z = zp + a * node_distances[j]
                pos = f'({z.real}, {z.imag})'
                if line_styles[j] == 1:
                    linestyle1 = 'line width=1.2pt,'
                    linestyle2 = f'line width=1.2pt, dash pattern=on {ds1}pt off {ds1*ds2}pt'
                    q1 = zp + r1 * (z - zp) / abs(z-zp)
                    p1 = make_pos(q1)
                    q2 = z - r1 * (z - zp) / abs(z-zp)
                    p2 = make_pos(q2)
                    res += fr'\draw[{linestyle1}] {posp} -- {p1};'
                    res += fr'\draw[{linestyle2}] {p1} -- {p2};'
                    res += fr'\draw[{linestyle1}] {p2} -- {pos};'
                else:
                    linestyle = 'line width=1.2pt,'
                    res += fr'\draw[{linestyle}] {posp} -- {pos};'

                fill = 'red!0'
                if i == nodes_filled[N][j]:
                    fill = 'blue'
                nodestyle = f'circle, fill={fill}, draw=black, minimum size={siz1}pt, inner sep=0pt'
                if i == 0:
                    nodestyle += f',label=below:{node_labels[j]}'
                res2 += fr'\node[{nodestyle}] () at {pos} {{}};'
                zp = z
                posp = pos

            zend = z + a * lend
            res3 += fr'\node[circle] () at {make_pos(zend)} {{{end_names[N][i]}}};'



        res += res2
        res += res3
        res += r'\end{tikzpicture}'
        res = str(res)
        f = f'/home/rah/s/phd2/tikz_arm{N}.tex'
        # print(res)
        open(f, 'w+').write(res)

def fp_arms2():
    ds1 = 1.
    ds2 = 1.8
    Ns = [2, 3, 4]
    l1 = 0.5
    l2 = 2*l1
    r1 = 0.5 * l1
    node_distances = [l1, l1, l2, l1]
    line_styles = [0, 0, 1, 0]
    siz1=6
    siz0=0.5 * siz1
    ncols = 1
    lend = 0.8 * l1
    lextra = 0.5 * l1
    ltotal = sum(node_distances) + lend + lextra
    nodes_filled = {
        2 : [0, 0, 1, 0],
        3 : [1, 0, 0, 2],
        4 : [0, 1, 0, 0],
        }
    end_names = {
        2 : ["$1$", "$\omega$"],
        3 : ["$1$", "$\omega$", "$\omega^2$"],
        4 : ["$1$", "$\omega$", "$\omega^2$", "$\omega^3$"],
        }
    node_labels = [
        "$1$", '$2$', '$L{-}1$', '$L$',
        ]
    node_label_scale = 0.8

    base_positions = [0.j, 0.8, 1.1]

    res = newliner()
    res += r'\begin{tikzpicture}'
    res2 = newliner()
    res3 = newliner()

    z0s = [0.]
    for Nind in range(1, len(Ns)):
        b = base_positions[Nind]
        a1 = 2.0 * ltotal
        a2 = 2.j * ltotal
        z0 = z0s[-1] + b.real * a1 + b.imag * a2
        z0s.append(z0)

    for Nind in range(len(Ns)):
        N = Ns[Nind]
        z0 = z0s[Nind]

        nodestyle = f'circle, fill=black, draw=black, minimum size={siz0}pt, inner sep=0pt'
        res2 += fr'\node[{nodestyle}] () at {make_pos(z0)} {{}};'

        for i in range(N):
            zp = z0
            posp = make_pos(zp)
            a = numpy.exp(2.j * numpy.pi * i / N)
            for j in range(len(node_distances)):
                z = zp + a * node_distances[j]
                pos = f'({z.real}, {z.imag})'
                if line_styles[j] == 1:
                    linestyle1 = 'line width=1.2pt,'
                    linestyle2 = f'line width=1.2pt, dash pattern=on {ds1}pt off {ds1*ds2}pt'
                    q1 = zp + r1 * (z - zp) / abs(z-zp)
                    p1 = make_pos(q1)
                    q2 = z - r1 * (z - zp) / abs(z-zp)
                    p2 = make_pos(q2)
                    res += fr'\draw[{linestyle1}] {posp} -- {p1};'
                    res += fr'\draw[{linestyle2}] {p1} -- {p2};'
                    res += fr'\draw[{linestyle1}] {p2} -- {pos};'
                else:
                    linestyle = 'line width=1.2pt,'
                    res += fr'\draw[{linestyle}] {posp} -- {pos};'

                fill = 'red!0'
                if i == nodes_filled[N][j]:
                    fill = 'blue'
                nodestyle = f'circle, fill={fill}, draw=black, minimum size={siz1}pt, inner sep=0pt'
                if i == 0:
                    nodestyle += f',label={{[scale={node_label_scale}]below:{node_labels[j]}}}'
                res2 += fr'\node[{nodestyle}] () at {pos} {{}};'
                zp = z
                posp = pos

            zend = z + a * lend
            res3 += fr'\node[circle] () at {make_pos(zend)} {{{end_names[N][i]}}};'



    res += res2
    res += res3
    res += r'\end{tikzpicture}'
    res = str(res)
    f = f'/home/rah/s/phd2/tikz_arm.tex'
    # print(res)
    open(f, 'w+').write(res)

def fp_arms3():
    ds1 = 1.
    ds2 = 1.8
    Ns = [2, 3, 4, 5]
    l1 = 0.5
    l2 = 2*l1
    r1 = 0.5 * l1
    node_distances = [l1, l1, l2, l1]
    line_styles = [0, 0, 1, 0]
    siz1=6
    siz0=0.5 * siz1
    ncols = 2
    lend = 0.8 * l1
    lextra = 0.5 * l1
    ltotal = sum(node_distances) + lend + lextra
    nodes_filled = {
        2 : [1, 0, 1, 0],
        3 : [1, 0, 0, 2],
        4 : [2, 1, 3, 0],
        5 : [0, 1, 0, 4],
        }
    end_names = {
        2 : ["$1$", "$\omega$"],
        3 : ["$1$", "$\omega$", "$\omega^2$"],
        4 : ["$1$", "$\omega$", "$\omega^2$", "$\omega^3$"],
        5 : ["$1$", "$\omega$", "$\omega^2$", "$\omega^3$", "$\omega^4$"],
        }
    node_labels = [
        # "$1$", '$2$', '$L{-}1$', '$L$',
        "$L$", '$L{-}1$', '$2$', '$1$',
        ]
    node_label_scale = 0.7


    res = newliner()
    res += r'\begin{tikzpicture}'
    res2 = newliner()
    res3 = newliner()

    z0s = []
    for Nind in range(len(Ns)):
        irow = Nind % ncols
        icol = Nind // ncols
        a1 = 2.0 * ltotal
        a2 = -2.j * ltotal
        z0 = a1 * irow + a2 * icol
        z0s.append(z0)

    for Nind in range(len(Ns)):
        N = Ns[Nind]
        z0 = z0s[Nind]

        nodestyle = f'circle, fill=black, draw=black, minimum size={siz0}pt, inner sep=0pt'
        res2 += fr'\node[{nodestyle}] () at {make_pos(z0)} {{}};'

        for i in range(N):
            zp = z0
            posp = make_pos(zp)
            a = numpy.exp(2.j * numpy.pi * i / N)
            for j in range(len(node_distances)):
                z = zp + a * node_distances[j]
                pos = f'({z.real}, {z.imag})'
                if line_styles[j] == 1:
                    linestyle1 = 'line width=1.2pt,'
                    linestyle2 = f'line width=1.2pt, dash pattern=on {ds1}pt off {ds1*ds2}pt'
                    q1 = zp + r1 * (z - zp) / abs(z-zp)
                    p1 = make_pos(q1)
                    q2 = z - r1 * (z - zp) / abs(z-zp)
                    p2 = make_pos(q2)
                    res += fr'\draw[{linestyle1}] {posp} -- {p1};'
                    res += fr'\draw[{linestyle2}] {p1} -- {p2};'
                    res += fr'\draw[{linestyle1}] {p2} -- {pos};'
                else:
                    linestyle = 'line width=1.2pt,'
                    res += fr'\draw[{linestyle}] {posp} -- {pos};'

                fill = 'red!0'
                if i == nodes_filled[N][j]:
                    fill = 'blue!70'
                nodestyle = f'circle, fill={fill}, draw=black, minimum size={siz1}pt, inner sep=0pt'
                if i == 0:
                    nodestyle += f',label={{[scale={node_label_scale}]below:{node_labels[j]}}}'
                res2 += fr'\node[{nodestyle}] () at {pos} {{}};'
                zp = z
                posp = pos

            zend = z + a * lend
            res3 += fr'\node[circle] () at {make_pos(zend)} {{{end_names[N][i]}}};'



    res += res2
    res += res3
    res += r'\end{tikzpicture}'
    res = str(res)
    f = f'/home/rah/s/phd2/tikz_arm2.tex'
    # print(res)
    open(f, 'w+').write(res)

def fp_coords3():
    res = newliner()
    # how deep
    depth = 4
    fsize=12
    N = 3

    # define epsilon colours
    epsilon_colours = [
        'red',
        'blue',
        'green',
        'orange',
        ]

    epsilon_colours = []
    cmap = matplotlib.cm.rainbow
    cmap = matplotlib.cm.plasma
    light = 0.8
    offset_left = 0.0
    offset_right = 0.3
    for d in range(depth):
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
    Epointsize=6
    base_length = 3.5
    line_lengths = [1., 0.5, 0.25, 0.12]
    line_angles = [0.3, -0.11, -0.1, 0.2]
    line_lengths = [base_length * x for x in line_lengths]
    line_lengths = [line_lengths[i] * numpy.exp(2.j * numpy.pi * line_angles[i] / N) for i in range(depth)]
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
            q += abs(line_lengths[d])
        x, y = q/2. + 1.5, q/2.*numpy.sqrt(3)
        # res += fr'\node[circle,fill=red,minimum size=25pt, inner sep=0pt] () at ({x}, {y}) {{}};'

        base = (x, y)
        x = base[0]
        y = base[1]
        for d in range(depth):
            l = abs(line_lengths[d])
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
            q += abs(line_lengths[d])
        x, y = -q, -q/2.*numpy.sqrt(3) +0.1
        # res += fr'\node[circle,fill=red,minimum size=25pt, inner sep=0pt] () at ({x}, {y}) {{}};'

        base = (x, y)
        x = base[0]
        y = base[1]
        for d in range(depth):
            l = abs(line_lengths[d])
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
    f = '/home/rah/s/phd2/tikz_spectrum2.tex'
    # print(res)
    open(f, 'w+').write(res)

def build():
    c = ['plx', 'tikz_review']
    p = proc(' '.join(c), cwd='/home/rah/s/phd2', shell=True)
    print(p)

from tikz_figures_1 import *
fp_coords2()
fp_coords3()
fp_coords4()
fp_arms2()
fp_arms3()
build()
