import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
SHOW = False

# axes
POINTS_BETWEEN_TICKS = 10
LOG_START_C = 0
LOG_END_C = 4
TOTAL_POINTS_C = (LOG_END_C - LOG_START_C) * POINTS_BETWEEN_TICKS + 1
CRANGE = np.logspace(LOG_START_C, LOG_END_C, TOTAL_POINTS_C)
#CRANGE2 = np.logspace(LOG_START_C, LOG_END_C+1, 2*TOTAL_POINTS_C)
LOG_START_KOFF = 0
LOG_END_KOFF = 2
LOG_START_KOFF2 = 2
LOG_END_KOFF2 = 4
TOTAL_POINTS_KOFF = (LOG_END_KOFF - LOG_START_KOFF) * POINTS_BETWEEN_TICKS + 1
KOFFRANGE = np.logspace(LOG_START_KOFF, LOG_END_KOFF, TOTAL_POINTS_KOFF)
KOFFRANGE2 = np.logspace(LOG_START_KOFF2, LOG_END_KOFF2, TOTAL_POINTS_KOFF) + 0.1
#KOFFRANGE2 = KOFFRANGE + 0.1

def heatmap(ax, arr, xrange, yrange, xy_label, label, log_norm=True,
                 xy_label_force=None, **kwargs):
    """
    xrange: range of values for x
    yrange: range of y values
    fname: file saved as pdf and eps format with this name
    show: shows plots
    save: saves plot
    log_norm: heatmap is log, this is only used when plotting the posterior
    kwargs: a variety of keyword args are used below. They are mostly used for contour plot lines if I understand correctly. Using Duncan's default ones mostly, but we can talk about it.
    """
    # default parameters, this is relic ofthe previous way we did things. Still might be useful.
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

    if log_norm:
        vmin = np.min( arr[arr > 0.0] )
        imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax, 'norm': mpl.colors.LogNorm(vmin,vmax)}
        print("Logscale")
    else:
        imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax}

    # TODO change colour scheme, see https://matplotlib.org/examples/color/colormaps_reference.html
    # TODO fix ticks randomly disappearing on colourbar + flip colourbar minor ticks or remove?
    """
    Colours viridis, YlGnBu, terrain, plasma
    """
    #print 'arr limits:', np.min(arr), np.max(arr)

    #im = plt.imshow(arr, interpolation='spline36', **imshow_kw)
    arr = np.ma.array(arr, mask=( arr<=0.0 ) )
    #arr[arr <= ]

    im = ax.imshow(arr, interpolation='none', **imshow_kw)
    #im = plt.imshow(arr, **imshow_kw)

    # axes setup

    # axes log scaled
    ax.set_xticks([i for i, cval in enumerate(xrange) if i % POINTS_BETWEEN_TICKS == 0])
    ax.set_yticks([i for i, kval in enumerate(yrange) if i % POINTS_BETWEEN_TICKS == 0])
    ax.set_xticklabels([r'$10^{%d}$' % np.log10(xval) for i, xval in enumerate(xrange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
    ax.set_yticklabels([r'$10^{%d}$' % np.log10(yval) for i, yval in enumerate(yrange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)

    ax.invert_yaxis()
    ax.set_xlabel(xy_label[0], fontsize=FS); ax.set_ylabel(xy_label[1], fontsize=FS)

    # create colorbar
    #divider = make_axes_locatable(ax) # this does not work, proposed solution online for colorbar for multiplots
    #cax = divider.append_axes("right",size="5%",pad=0.05)
    #cbar = plt.colorbar(im, ax=cax)

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=FS, labelpad=20); cbar.ax.tick_params(labelsize=FS)
    cbar.ax.minorticks_off();
    #cbar.set_ticks([round(vmin,3)+0.001,round(vmax,3)-0.001]) # UNCOMMENT THIS ONLY WHEN TICKS DON'T APPEAR
    cbar.update_ticks()

    # TODO IDK why do ticks hide sometimes?
    CL = ax.contour(arr, levels=levels, linestyles=contour_linestyle, colors=contour_color, linewidths=contour_lindewidths)
    ax.clabel(CL, CL.levels, inline=True, fmt=fmt)

    plt.tight_layout(h_pad=1)

    return ax

def single_heatmap(arr, xrange, yrange, fname, xy_label, label, show=SHOW, save=True, log_norm=True,
                 xy_label_force=None, **kwargs):
    # Make only 1 heatmap at a time.

    f = plt.figure()
    # axes setup
    fig = plt.gcf(); ax = plt.gca()

    heatmap(ax, arr, xrange, yrange, xy_label, label, log_norm=log_norm, xy_label_force=None, **kwargs)

    # save
    if save == True:
        plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.png'); #plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.eps');
    if show:
        plt.show()

    plt.close()
    return fig, ax

def multiple_heatmaps(arrDetSigmaEst, arrRelErrorEst, array_x, array_y, fname, labels, log_select):
    # makes a figure with many subplots.

    fig = plt.figure(figsize=(24,10));
    ax0 = plt.subplot2grid((2,4), (0,0), colspan=2, rowspan=2);
    ax1 = plt.subplot2grid((2,4), (0,2), colspan=1, rowspan=1);
    ax2 = plt.subplot2grid((2,4), (0,3), colspan=1, rowspan=1);
    ax3 = plt.subplot2grid((2,4), (1,2), colspan=1, rowspan=1);
    ax4 = plt.subplot2grid((2,4), (1,3), colspan=1, rowspan=1);

    # Determinant plot
    heatmap(ax0, arrDetSigmaEst[:,:], array_x, array_y, labels, r'Det($\Sigma_{est}$)', log_norm=True)

    # each of the diagonals
    heatmap(ax1, arrRelErrorEst[:,:,0,0], array_x, array_y, labels, r'$\langle \delta {c_1}^2 \rangle / {c_1}^2$', log_norm=log_select)
    heatmap(ax2, arrRelErrorEst[:,:,1,1], array_x, array_y, labels, r'$\langle \delta {k_{off}}^2 \rangle / {k_{off}}^2$', log_norm=log_select)
    heatmap(ax3, arrRelErrorEst[:,:,2,2], array_x, array_y, labels, r'$\langle \delta {c_2}^2 \rangle / {c_2}^2$', log_norm=log_select)
    heatmap( ax4, arrRelErrorEst[:,:,3,3], array_x, array_y, labels, r'$\langle \delta {k_{off,2}}^2 \rangle / {k_{off,2}}^2$', log_norm=log_select)

    # save
    #fig.suptitle(r"IFN Crosstalk $k_{off}/k_{off,2}=100$", fontsize=int(1.5*FS)) # TODO customize the title
    fig.suptitle(r"T Cell antigen detection $c_2/c_1=1000$", fontsize=int(1.5*FS))
    plt.tight_layout(h_pad=1, rect=[0,0.03,1,0.95]) # need to change this to not have the title overlap, figuring it out still

    plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.png'); #plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.eps');

    plt.close()

    return fig


def sigmaEst(c1, koff, c2, koff2):
    # eqns2l are the equations imported from Mathematica and turned into matrices (function of c1,c2,koff,koff2)
    A = eqns2l.matrix_dmudthetaInv(c1, c2, koff, koff2)
    B = eqns2l.matrix_sigmadata(c1, c2, koff, koff2)
    Atrans = np.transpose(A)
    matrix = np.linalg.multi_dot([A, B, Atrans]) # matrix

    # use to make the entries relative estimates
    rel = np.array([ [c1*c1, c1*koff, c1*c2, c1*koff2], [koff*c1, koff*koff, koff*c2, koff*koff2], [c2*c1, c2*koff, c2*c2, c2*koff2], [koff2*c1, koff2*koff, koff2*c2, koff2*koff2] ])
    relErrorMatrix = np.divide(matrix,rel)

    return matrix, np.linalg.det(matrix), relErrorMatrix

def determinants_mathematica(c1, koff, c2, koff2):
    # Calculating the determinant of matrixes imported from Matehmatica
    return np.linalg.det( eqns2l.matrix_sigmadata(c1, c2, koff, koff2) ), np.linalg.det( eqns2l.matrix_dmudthetaInv(c1, c2, koff, koff2) )


if __name__ == '__main__':

    # choose any of these 2 to be arrays [0: c1, 1: koff, 2: c2, 3: koff2], they will be your axes
    dim = {'x' : 1, 'y' : 3}

    # axis we'd like to plot
    value = [C1, KOFF, C2, KOFF2]
    label_fig = ['c1', 'koff', 'c2', 'koff2']
    label = [r'$c_1$', r'$k_{off}$', r'$c_2$', r'$k_{off,2}$']
    dimension = [CRANGE, KOFFRANGE, CRANGE, KOFFRANGE2]
    axes = label_fig[dim['x']]+ '_' + label_fig[dim['y']]
    multi_fname = "multi_"+axes

    # for plotting the dedimensionalized values
    dedimension = [dimension[0]*KON/KP, dimension[1]/KP, dimension[2]*KON/KP, dimension[3]/KP]
    dedimension_label = [r'$k_{on}c_1/k_{p}$', r'$k_{off}/k_{p}$', r'$k_{on}c_2/k_{p}$', r'$k_{off,2}/k_{p}$']

    # heatmap for any element of sigmaEst(4x4 matrix)
    arrSigmaEst = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4) )
    arrRelErrorEst = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4) )

    # heatmap for the det estimate covariance
    arrDetSigmaEst = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )
    arrDetSigmaData = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )
    arrDetdmudthetaInv = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )

    # elements of certain matrices, less interested in those right now
    #element11SigmaData = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )
    #element31SigmaData = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )
    #element11dmudthetaInv = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )
    #element31dmudthetaInv = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )

    # making our heatmap,
    for i, xi in enumerate( dimension[dim['x']] ):
        value[dim['x']] = xi
        for j, yi in enumerate( dimension[dim['y']] ):
            value[dim['y']] = yi
            arrSigmaEst[j,i,:,:], arrDetSigmaEst[j,i], arrRelErrorEst[j,i,:,:] = sigmaEst(value[0], value[1], value[2], value[3])
            arrDetSigmaData[j,i], arrDetdmudthetaInv[j,i] = determinants_mathematica(value[0], value[1], value[2], value[3])
            #element11SigmaData[j,i] = eqns2l.sigmadata11(value[0], value[2], value[1], value[3])
            #element31SigmaData[j,i] = eqns2l.sigmadata31(value[0], value[2], value[1], value[3])
            #element11dmudthetaInv[j,i] = eqns2l.dmudthetaInv11(value[0], value[2], value[1], value[3])
            #element31dmudthetaInv[j,i] = eqns2l.dmudthetaInv31(value[0], value[2], value[1], value[3])


    # Here I have selected the plots I want, when we are making individual heatmaps. LOG_SELECT is kinda old, used to be for when we didn't want log colorbar, but we nearly always want it now
    """
    LOG_SELECT = True
    single_heatmap(arrDetSigmaEst[:,:], dedimension[dim['x']], dedimension[dim['y']], 'DetSigmaEst_'+axes, [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Det($\Sigma_{est}$)', log_norm=LOG_SELECT)
    single_heatmap(arrDetSigmaData[:,:], dedimension[dim['x']], dedimension[dim['y']], 'DetSigmaData_'+axes, [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Det($\Sigma_{data}$)', log_norm=LOG_SELECT)

    single_heatmap( arrRelErrorEst[:,:,0,0], dedimension[dim['x']], dedimension[dim['y']], 'SigmaEst_00_'+axes, [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'$\langle \delta {c_1}^2 \rangle / {c_1}^2$', log_norm=LOG_SELECT)
    single_heatmap( arrRelErrorEst[:,:,1,1], dedimension[dim['x']], dedimension[dim['y']], 'SigmaEst_11_'+axes, [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'$\langle \delta {k_{off}}^2 \rangle / {k_{off}}^2$', log_norm=LOG_SELECT)
    single_heatmap( arrRelErrorEst[:,:,2,2], dedimension[dim['x']], dedimension[dim['y']], 'SigmaEst_22_'+axes, [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'$\langle \delta {c_2}^2 \rangle / {c_2}^2$', log_norm=LOG_SELECT)
    single_heatmap( arrRelErrorEst[:,:,3,3], dedimension[dim['x']], dedimension[dim['y']], 'SigmaEst_33_'+axes, [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'$\langle \delta {k_{off,2}}^2 \rangle / {k_{off,2}}^2$', log_norm=LOG_SELECT)
    """

    # make a figure with multiple subplots
    LOG_SELECT = True
    multiple_heatmaps(arrDetSigmaEst, arrRelErrorEst,  dedimension[dim['x']],  dedimension[dim['y']], multi_fname, [ dedimension_label[dim['x']], dedimension_label[dim['y']]], LOG_SELECT)




    # old plots to check if conversion from matehmatica worked fine
    #plot_heatmap(arrDetdmudthetaInv[:,:], dedimension[dim['x']], dedimension[dim['y']], 'c1c2detdmudtheta', [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Determinant $d\mu /d \theta^{-1}$')
    #plot_heatmap(element11SigmaData[:,:], dedimension[dim['x']], dedimension[dim['y']], 'c1c2sigmadata11', [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Determinant $\Sigma_{data,11}$')
    #plot_heatmap(element31SigmaData[:,:], dedimension[dim['x']], dedimension[dim['y']], 'c1c2sigmadata31', [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Determinant $\Sigma_{data,31}$')
    #plot_heatmap(element11dmudthetaInv[:,:], dedimension[dim['x']], dedimension[dim['y']], 'c1c2detdmudtheta11', [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Determinant $d\mu /d \theta^{-1}_{11}$')
    #plot_heatmap(element31dmudthetaInv[:,:], dedimension[dim['x']], dedimension[dim['y']], 'c1c2detdmudtheta31', [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Determinant $d\mu /d \theta^{-1}_{31}$')
