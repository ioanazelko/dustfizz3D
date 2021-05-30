
from __future__ import division
# to use division the way it is used in python 3
import os
import sys


import corner
import matplotlib as mpl
import matplotlib.backends.backend_pdf
from matplotlib.collections import LineCollection  # For plotting
import matplotlib.pyplot as plt
import numpy as np




def plot_parameters_corner(data,var, fiducial,bins,plot_countours):
    samples = data
    fontsize = 20
    ndim = len(fiducial)

    print("Fiducial is", list(zip(var,fiducial)))
    figure = corner.corner(samples, labels = var, label_kwargs={"fontsize":fontsize+10}, bins=bins,show_titles=True,scale_hist = True,title_fmt=".1e",truths=fiducial,plot_countours=plot_countours)
    # # Extract the axes
    axes = np.array(figure.axes).reshape((ndim, ndim))

    # Loop over the diagonal
    for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(fiducial[i], color="g")

    # Loop over the histograms
    for yi in range(ndim):
        for xi in range(yi):
            ax = axes[yi, xi]

            ax.axvline(fiducial[xi], color="g")
            ax.axhline(fiducial[yi], color="g")

    return figure
