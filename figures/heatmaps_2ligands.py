import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

import plotting_dictionaries as pd
import matplotlib.ticker as ticker
import equations2ligands as eqns2l
from decimal import Decimal


#from load_inputs import DATADICT
from settings import DIR_OUTPUT, DIR_INPUT
#from settings import COLOR_SCHEME as cs
from settings import KON, KP, T, KF, ALPHA, C1, C2, KOFF, KOFF2


plt.style.use('parameters.mplstyle')  # particularIMporting

# plot params
FS = 20
SHOW = True

# axes
POINTS_BETWEEN_TICKS = 10
LOG_START_C = -2
LOG_END_C = 2
TOTAL_POINTS_C = (LOG_END_C - LOG_START_C) * POINTS_BETWEEN_TICKS + 1
CRANGE = np.logspace(LOG_START_C, LOG_END_C, TOTAL_POINTS_C)
LOG_START_KOFF = -2
LOG_END_KOFF = 2
TOTAL_POINTS_KOFF = (LOG_END_KOFF - LOG_START_KOFF) * POINTS_BETWEEN_TICKS + 1
KOFFRANGE = np.logspace(LOG_START_KOFF, LOG_END_KOFF, TOTAL_POINTS_KOFF)

def plot_heatmap(arr, crange, koffrange, fname, xy_label, label, show=SHOW, save=True, log_norm=True,
                 xy_label_force=None, **kwargs):
    """
    crange: range of values for c
    koffrange: range of koff values
    fname: file saved as pdf and eps format with this name
    show: shows plots
    save: saves plots
    log_norm: heatmap is log, this is only used when plotting the posterior
    kwargs: a variety of keyword args are used below. They are mostly used for contour plot lines if I understand correctly. Using Duncan's default ones mostly, but we can talk about it.
    """
    # default parameters
    if 'levels' in kwargs.keys(): levels = kwargs['levels']
    else: levels = [1E0, 1E5, 1E10, 1E20]

    if 'vmin' in kwargs.keys(): vmin = kwargs['vmin']
    else: vmin = np.min(arr)

    if 'vmax' in kwargs.keys(): vmax = kwargs['vmax']
    else: vmax = np.max(arr)

    if 'contour_linestyle' in kwargs.keys(): contour_linestyle = kwargs['contour_linestyle']
    else: contour_linestyle = '-'

    if 'contour_color' in kwargs.keys(): contour_color = kwargs['contour_color']
    else: contour_color = 'k'

    if 'contour_linewidths' in kwargs.keys(): contour_lindewidths = kwargs['contour_linewidths']
    else: contour_lindewidths = 2

    if 'cmap_colour' in kwargs.keys(): cmap_colour = kwargs['cmap_colour']
    else: cmap_colour = 'YlGnBu'

    if 'fmt' in kwargs.keys(): fmt = kwargs['fmt']
    else: fmt = ticker.LogFormatterMathtext()

    print(vmin,vmax)

    if log_norm:
        imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax, 'norm': mpl.colors.LogNorm(vmin,vmax)}
    else:
        imshow_kw = {'cmap': 'viridis', 'aspect': None, 'vmin': vmin, 'vmax': vmax}

    # TODO change colour scheme, see https://matplotlib.org/examples/color/colormaps_reference.html
    # TODO fix ticks randomly disappearing on colourbar + flip colourbar minor ticks or remove?
    """
    Colours viridis, YlGnBu, terrain, plasma
    """
    #print 'arr limits:', np.min(arr), np.max(arr)
    # plot setup
    f = plt.figure()
    #im = plt.imshow(arr, interpolation='spline36', **imshow_kw)
    im = plt.imshow(arr, interpolation='hamming', **imshow_kw)
    #im = plt.imshow(arr, **imshow_kw)

    # axes setup
    fig = plt.gcf(); ax = plt.gca()

    # method 1
    if log_norm:
        ax.set_xticks([i for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS == 0])
        ax.set_yticks([i for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0])
        ax.set_xticklabels([r'$10^{%d}$' % np.log10(cval) for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
        ax.set_yticklabels([r'$10^{%d}$' % np.log10(kval) for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
    else:
        # Set ticks later because for some reason plotting contours messes up the ticks that get set here
        pass

    ax.invert_yaxis()
    ax.set_xlabel(xy_label[0], fontsize=FS); ax.set_ylabel(xy_label[1], fontsize=FS)

    # create colorbar
    cbar = fig.colorbar(im)
    #cbar.locator = ticker.LogLocator(base=10)
    cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=FS, labelpad=20); cbar.ax.tick_params(labelsize=FS)
    cbar.ax.minorticks_off();
    # UNCOMMENT THIS ONLY WHEN TICKS DON'T APPEAR
    #cbar.set_ticks([round(vmin,3)+0.001,round(vmax,3)-0.001])
    cbar.update_ticks()

    # TODO IDK why do ticks hide sometimes?
    CL = plt.contour(arr, levels=levels, linestyles=contour_linestyle, colors=contour_color, linewidths=contour_lindewidths)
    plt.clabel(CL, CL.levels, inline=True, fmt=fmt)

    # save
    if save == True:
        plt.savefig(DIR_OUTPUT + os.sep + fname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + fname + '.eps'); plt.savefig(DIR_OUTPUT + os.sep + fname + '.png')
    if show:
        plt.show()

    plt.close()
    return fig, ax

def sigmaEst(c1, koff, c2, koff2):
    # eqns2l are the equations imported from Mathematica and turned into matrices (function of c1,c2,koff,koff2)
    A = eqns2l.matrix_dmudthetaInv(c1, c2, koff, koff2)
    B = eqns2l.matrix_sigmadata(c1, c2, koff, koff2)
    Atrans = np.transpose(A)
    matrix = np.linalg.multi_dot([A, B, Atrans]) # matrix multiplication

    return matrix, np.linalg.det(matrix)


if __name__ == '__main__':


    value = [C1, KOFF, C2, KOFF2]
    label = [r'$c_1$', r'$k_{\text{off}}$', r'$c_2$', r'$k_{\text{off,2}}$']
    dimension = [CRANGE, KOFFRANGE, CRANGE, KOFFRANGE]

    # choose any of these 2 to be arrays
    dim = {'x' : 0, 'y' : 2, 'name' : 'c1c2DET2'}

    arrSigmaEst = np.zeros( (len(dimension[dim['x']]), len(dimension[dim['y']]), 4, 4) ) # heatmap for any element of sigmaEst(4x4 matrix)
    arrDetSigmaEst = np.zeros( (len(dimension[dim['x']]), len(dimension[dim['y']]) ) ) # heatmap for the determinant

    # making our heatmap,
    for i, xi in enumerate( dimension[dim['x']] ):
        value[dim['x']] = xi
        for j, yi in enumerate( dimension[dim['y']] ):
            value[dim['y']] = yi
            arrSigmaEst[i,j,:,:], arrDetSigmaEst[i,j] = sigmaEst(value[0], value[1], value[2], value[3])

    # for plotting the dedimensionalized values
    dedimension = [dimension[0]*KON/KP, dimension[1]/KP, dimension[3]*KON/KP, dimension[3]/KP]
    dedimension_label = [r'$k_{on}c_1/k_{p}$', r'$k_{\text{off}}/k_{p}$', r'$k_{on}c_2/k_{p}$', r'$k_{\text{off,2}}/k_{p}$']

    # Here I have selected the determinnant to plot. TODO make arr potentially about
    arr = arrDetSigmaEst[:,:]; label = r'Determinant $\Sigma_{est}$';

    plot_heatmap(arr, dedimension[dim['x']], dedimension[dim['y']], dim['name'], [ dedimension_label[dim['x']], dedimension_label[dim['y']]], label)
