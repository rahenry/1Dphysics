#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt

preamb = [r"\usepackage{amsmath}", r"\usepackage{amssymb}", r"\usepackage{amsfonts}", r"\setlength{\parindent}{0pt}"]
preamb = "\n".join(preamb)
plt.rcParams.update(
    {
        "text.latex.preamble": preamb,
        "font.family": "serif",
        "text.usetex": True,
        "pgf.texsystem": "pdflatex",
        "pgf.rcfonts": False,
    }
)
plt.rcParams["pgf.preamble"] = preamb
