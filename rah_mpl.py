#!/usr/bin/env python3

import matplotlib.transforms as transforms
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt

def draw_label(ax, label, position='upper right'):
    text_box = AnchoredText(label, frameon=True, loc=position, pad=0)
    # text_box.patch.set_boxstyle("facecolor='lightgrey',linestyle=''")
    text_box.patch.set_boxstyle("square,pad=0.15")
    plt.setp(text_box.patch, linestyle='', facecolor='lightgrey')
    # ax.setp(text_box.patch, facecolor='white', alpha=0.5)
    ax.add_artist(text_box)
