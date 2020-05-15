import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import scipy as sp
import os

import plotting_dictionaries as pd
import matplotlib.ticker as ticker
import equations2ligands as eqns2l
import equations as eqns
from decimal import Decimal


#from load_inputs import DATADICT
from heatmaps import combined_error_c, combined_error_koff
from settings import DIR_OUTPUT, DIR_INPUT, KON, KP, T, KF, ALPHA, C1, C2, KOFF, KOFF2
#from settings import COLOR_SCHEME as cs



plt.style.use('parameters.mplstyle')  # particularIMporting

# plot params
FS = 22
SHOW = False

# axes
POINTS_BETWEEN_TICKS = 12
LOG_START_C = -9
LOG_END_C = -4
TOTAL_POINTS_C = (LOG_END_C - LOG_START_C) * POINTS_BETWEEN_TICKS + 1
CRANGE = np.logspace(LOG_START_C, LOG_END_C, TOTAL_POINTS_C)
LOG_START_KOFF = 0
LOG_END_KOFF = 2
LOG_START_KOFF2 = -1
LOG_END_KOFF2 = 3
TOTAL_POINTS_KOFF = (LOG_END_KOFF - LOG_START_KOFF) * POINTS_BETWEEN_TICKS + 1
TOTAL_POINTS_KOFF2 = (LOG_END_KOFF2 - LOG_START_KOFF2) * POINTS_BETWEEN_TICKS + 1
KOFFRANGE = np.logspace(LOG_START_KOFF, LOG_END_KOFF, TOTAL_POINTS_KOFF)
KOFFRANGE2 = np.logspace(LOG_START_KOFF2, LOG_END_KOFF2, TOTAL_POINTS_KOFF2)


class MidPointLogNorm(mpl.colors.LogNorm):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        mpl.colors.LogNorm.__init__(self,vmin=vmin, vmax=vmax, clip=clip)
        self.midpoint=midpoint
    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [np.log(self.vmin), np.log(self.midpoint), np.log(self.vmax)], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(np.log(value), x, y))


def heatmap(ax, arr, xrange, yrange, xy_label, label, log_norm=True, xy_label_force=None, **kwargs):
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
    else: levels = [1E0]#, 1E5, 1E10, 1E20]

    if 'vmin' in kwargs.keys(): vmin = kwargs['vmin']
    else: vmin = 1E-3#np.min(arr)#1E-3#

    if 'vmax' in kwargs.keys(): vmax = kwargs['vmax']
    else:
        if np.max(arr) > 1E10:
            vmax = 1E5
        else:
            vmax = 1E5#np.max(arr)#1E5#

    if 'contour_linestyle' in kwargs.keys(): contour_linestyle = kwargs['contour_linestyle']
    else: contour_linestyle = '-'

    if 'contour_color' in kwargs.keys(): contour_color = kwargs['contour_color']
    else: contour_color = 'k'

    if 'contour_linewidths' in kwargs.keys(): contour_lindewidths = kwargs['contour_linewidths']
    else: contour_lindewidths = 2

    if 'cmap_colour' in kwargs.keys(): cmap_colour = kwargs['cmap_colour']
    else: cmap_colour = 'RdBu_r'#'YlGnBu'

    if 'fmt' in kwargs.keys(): fmt = kwargs['fmt']
    else: fmt = ticker.LogFormatterMathtext()

    if log_norm:
        arr = np.ma.array(arr, mask=( arr<=0.0 ) )
        imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax, 'norm': mpl.colors.LogNorm(vmin,vmax)}
        imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax, 'norm': MidPointLogNorm(vmin,vmax,1.0)}
        print("Logscale")
    else:
        imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax}

    # TODO change colour scheme, see https://matplotlib.org/examples/color/colormaps_reference.html
    # TODO fix ticks randomly disappearing on colourbar + flip colourbar minor ticks or remove?
    """
    Colours viridis, YlGnBu, terrain, plasma
    """

    im = ax.imshow(arr, interpolation='none', **imshow_kw)
    #im = plt.imshow(arr, interpolation='spline36', **imshow_kw)

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
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04) # Fig 3 way

    cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=FS, labelpad=20); cbar.ax.tick_params(labelsize=FS)
    cbar.ax.minorticks_off();
    #cbar.set_ticks([round(vmin,3)+0.001,round(vmax,3)-0.001]) # UNCOMMENT THIS ONLY WHEN TICKS DON'T APPEAR
    cbar.update_ticks()

    # TODO IDK why do ticks hide sometimes?
    CL = ax.contour(arr, levels=levels, linestyles=contour_linestyle, colors=contour_color, linewidths=contour_lindewidths)
    ax.clabel(CL, CL.levels, inline=True, fmt=fmt, fontsize=FS-4)

    #Might be an issue
    #plt.tight_layout(h_pad=1)
    plt.tight_layout(h_pad=1,w_pad=1)

    return ax, cbar, im

def single_heatmap(arr, xrange, yrange, fname, xy_label, label, show=SHOW, save=True, log_norm=True,
                 xy_label_force=None, **kwargs):
    """
    arr: array to used
    xrange, yrange: the range for the axis
    fname: filename to save
    xy_label: the label of the x,y axis (need this to be an array of 2 values, both strings)
    show, save: whehther to show or save the plots
    log_norm: whether or not to have a log norm for the heatmap
    """
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

def multiple_heatmaps(arrRelDetSigmaEst, arrRelErrorEst, array_x, array_y, fname, labels, log_select):
    """
    arr*: various 2D arrays to used for heatmap
    array_x, array_y: the range for the axis
    fname: filename to save
    labels: the label of the x,y axis (need this to be an array of 2 values, both strings)
    log_select: whether or not to have a log norm for the heatmap
    """
    # makes a figure with many subplots.

    #fig = plt.figure(figsize=(24,11));
    fig = plt.figure(figsize=(23,10));
    #gs = fig.add_gridspec(2,4, hspace=0.05, wspace=0.05, width_ratios=[1.,1.,1.,1.], height_ratios=[1.,1.])

    #ax0 = fig.add_subplot(gs[:,:-2])
    det_on_left = False
    if det_on_left:
        gs = fig.add_gridspec(2, 6, hspace=-0.01, wspace=0.05, width_ratios=[1., 1., 0.08, 1., 1., 0.1], height_ratios=[1., 1.])

        ax0 = fig.add_subplot(gs[:,:2])  # the det
        ax1 = fig.add_subplot(gs[0,3])
        ax2 = fig.add_subplot(gs[0,4])
        ax3 = fig.add_subplot(gs[1,3])
        ax4 = fig.add_subplot(gs[1,4])
        ax5 = fig.add_subplot(gs[:,-1])   # the cbar

        # Determinant plot
        """
        #heatmap(ax0, arrDetSigmaEst[:,:], array_x, array_y, labels, r'Det($\Sigma_{est}$)', log_norm=True)
        ax0, cbar0, im0 = heatmap(ax0, arrRelDetSigmaEst[:,:], array_x, array_y, labels, r'Det($\Sigma_{est}$)/(${c_1}^2{c_1}^2{k_{\mathrm{off},1}}^2{k_{\mathrm{off},2}}^2$)', log_norm=True)
        cbar0.remove(); ax0.set_title(r'det($\Sigma_{est}$)/(${c_1}^2{c_2}^2{k_{\mathrm{off},1}}^2{k_{\mathrm{off},2}}^2$)', fontsize=FS)

        # each of the diagonals
        ax1, cbar1, _ = heatmap(ax1, arrRelErrorEst[:,:,0,0], array_x, array_y, labels, r'$\langle \delta {c_1}^2 \rangle / {c_1}^2$', log_norm=log_select)
        cbar1.remove(); ax1.tick_params(labelbottom=False); ax1.tick_params(labelleft=False); ax1.set_xticklabels([]); ax1.set_yticklabels([]); ax1.set_ylabel(''); ax1.xaxis.set_label_position('top'); ax1.set_xlabel(r'$\langle \delta {c_1}^2 \rangle / {c_1}^2$')

        a2, cbar2, _ = heatmap(ax2, arrRelErrorEst[:,:,1,1], array_x, array_y, labels, r'$\langle \delta {k_{off}}^2 \rangle / {k_{off}}^2$', log_norm=log_select)
        cbar2.remove(); ax2.tick_params(labelbottom=False); ax2.tick_params(labelleft=False); ax2.set_xticklabels([]); ax2.set_yticklabels([]); ax2.set_ylabel(''); ax2.xaxis.set_label_position('top'); ax2.set_xlabel(r'$\langle \delta {k_{\mathrm{off}}}^2 \rangle / {k_{\mathrm{off}}}^2$')

        a3, cbar3, _ = heatmap(ax3, arrRelErrorEst[:,:,2,2], array_x, array_y, labels, r'$\langle \delta {c_2}^2 \rangle / {c_2}^2$', log_norm=log_select)
        cbar3.remove(); ax3.tick_params(labelbottom=False); ax3.tick_params(labelleft=False); ax3.set_xticklabels([]); ax3.set_yticklabels([]); ax3.set_ylabel(''); ax3.xaxis.set_label_position('top'); ax3.set_xlabel( r'$\langle \delta {c_2}^2 \rangle / {c_2}^2$')

        a4, cbar4, _ = heatmap( ax4, arrRelErrorEst[:,:,3,3], array_x, array_y, labels, r'$\langle \delta {k_{off,2}}^2 \rangle / {k_{off,2}}^2$', log_norm=log_select)
        cbar4.remove(); ax4.tick_params(labelbottom=False); ax4.tick_params(labelleft=False); ax4.set_xticklabels([]); ax4.set_yticklabels([]); ax4.set_ylabel(''); ax4.xaxis.set_label_position('top'); ax4.set_xlabel( r'$\langle \delta {k_{\mathrm{off},2}}^2 \rangle / {k_{\mathrm{off},2}}^2$')

        cb = fig.colorbar(im0, cax=ax5); cb.ax.tick_params(labelsize=FS)
        """

    else:
        gs = fig.add_gridspec(2, 6, hspace=-0.12, wspace=0.05, width_ratios=[1., 1., 0.2, 1., 1., 0.1], height_ratios=[1., 1.])

        ax0 = fig.add_subplot(gs[:,3:5])  # the det
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])
        ax3 = fig.add_subplot(gs[1,0])
        ax4 = fig.add_subplot(gs[1,1])
        ax5 = fig.add_subplot(gs[:,-1])   # the cbar

    # Determinant plot
    #heatmap(ax0, arrDetSigmaEst[:,:], array_x, array_y, labels, r'Det($\Sigma_{est}$)', log_norm=True)
    ax0, cbar0, im0 = heatmap(ax0, arrRelDetSigmaEst[:,:], array_x, array_y, labels, r'Det($\Sigma_{est}$)/(${c_1}^2{c_1}^2{k_{\mathrm{off},1}}^2{k_{\mathrm{off},2}}^2$)', log_norm=True)
    cbar0.remove(); ax0.set_title(r'det($\Sigma_{est}$)/(${c_1}^2{c_2}^2{k_{\mathrm{off},1}}^2{k_{\mathrm{off},2}}^2$)', fontsize=FS)

    # each of the diagonals
    ax1, cbar1, _ = heatmap(ax1, arrRelErrorEst[:,:,0,0], array_x, array_y, labels, r'$\langle \delta {c_1}^2 \rangle / {c_1}^2$', log_norm=log_select)
    cbar1.remove(); ax1.tick_params(labelbottom=False); ax1.set_xticklabels([]);  ax1.set_xlabel('')
    ax1.set_title(r'$\langle \delta {c_1}^2 \rangle / {c_1}^2$', fontsize=FS)

    a2, cbar2, _ = heatmap(ax2, arrRelErrorEst[:,:,1,1], array_x, array_y, labels, r'$\langle \delta {k_{off}}^2 \rangle / {k_{off}}^2$', log_norm=log_select)
    cbar2.remove(); ax2.tick_params(labelbottom=False); ax2.tick_params(labelleft=False); ax2.set_xticklabels([]); ax2.set_yticklabels([]); ax2.set_ylabel(''); ax2.set_xlabel('')
    ax2.set_title(r'$\langle \delta {k_{\mathrm{off}}}^2 \rangle / {k_{\mathrm{off}}}^2$', fontsize=FS)

    a3, cbar3, _ = heatmap(ax3, arrRelErrorEst[:,:,2,2], array_x, array_y, labels, r'$\langle \delta {c_2}^2 \rangle / {c_2}^2$', log_norm=log_select)
    cbar3.remove()
    ax3.set_title(r'$\langle \delta {c_2}^2 \rangle / {c_2}^2$', fontsize=FS)

    a4, cbar4, _ = heatmap( ax4, arrRelErrorEst[:,:,3,3], array_x, array_y, labels, r'$\langle \delta {k_{off,2}}^2 \rangle / {k_{off,2}}^2$', log_norm=log_select)
    cbar4.remove(); ax4.tick_params(labelleft=False); ax4.set_yticklabels([]); ax4.set_ylabel('')
    ax4.set_title(r'$\langle \delta {k_{\mathrm{off},2}}^2 \rangle / {k_{\mathrm{off},2}}^2$', fontsize=FS)

    cb = fig.colorbar(im0, cax=ax5, fraction=0.9); cb.ax.tick_params(labelsize=FS)

    # save
    #plt.tight_layout(h_pad=0.05, w_pad=0.05, rect=[0.,0.,1.,1.]) # need to change this to not have the title overlap, figuring it out still

    # note tight layout seems incompatible with gridspec/subplots
    plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.pdf', bbox_inches='tight')
    plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.png', bbox_inches='tight')

    plt.close()

    return fig

def eigen_heatmaps(arrValues, arrVectors, array_x, array_y, fname, labels, direction_label):
    """
    arr*: various 2D arrays to used for heatmap
    array_x, array_y: the range for the axis
    fname: filename to save
    labels: the label of the x,y axis (need this to be an array of 2 values, both strings)
    log_select: whether or not to have a log norm for the heatmap
    """
    # makes a figure with many subplots, each row is for an eigenvalue
    flag_multiple = True

    if not flag_multiple:
        #fig = plt.figure(figsize=(24,11));
        fig = plt.figure(figsize=(30,24), constrained_layout=True);
        gs = fig.add_gridspec(4, 5, hspace=-0.12, wspace=0.05, width_ratios=[1., 1., 1., 1., 1.], height_ratios=[1., 1., 1., 1.])

        plot_projection_params = {'vmin' : -1., 'vmax' : 1. }

        for i in [0,1,2,3]:
            ax = fig.add_subplot(gs[i,0])
            # plotting eigenvalues
            plot_eigen_params = {'vmin' : np.min( arrValues[:,:,i]), 'vmax' : np.min( arrValues[:,:,i]), 'cmap_colour' : 'YlGnBu' }
            heatmap(ax, arrValues[:,:,i], array_x, array_y, labels, '', log_norm=False, **plot_eigen_params )
            ax.set_title( str(i)+r'$^{th}$ eigenvalue', fontsize=FS )

            for j in [0,1,2,3]:
                axj = fig.add_subplot(gs[i,1+j])
                # plotting components of the arrVectors
                heatmap(axj, arrVectors[:,:,j,i], array_x, array_y, labels, '', log_norm=False, **plot_projection_params)
                axj.set_title( r'projection on'+direction_label[j], fontsize=FS )

        plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.pdf', bbox_inches='tight')
        plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.png', bbox_inches='tight')

        plt.close()
    else:
        for i in [0,1,2,3]:
            #fig = plt.figure(figsize=(30,3), constrained_layout=True);
            fig = plt.figure(figsize=(30,3), tight_layout=True);
            fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.1)
            gs = fig.add_gridspec(1, 5, hspace=-0.12, wspace=0.05, width_ratios=[1., 1., 1., 1., 1.], height_ratios=[1.])

            plot_projection_params = {'vmin' : -1., 'vmax' : 1. }
            ax = fig.add_subplot(gs[0,0])
            # plotting eigenvalues
            plot_eigen_params = {'vmin' : np.min( arrValues[:,:,i]), 'vmax' : np.min( arrValues[:,:,i]), 'cmap_colour' : 'YlGnBu' }
            heatmap(ax, arrValues[:,:,i], array_x, array_y, labels, '', log_norm=False, **plot_eigen_params )
            ax.set_title( str(i)+r'$^{th}$ eigenvalue', fontsize=FS )

            for j in [0,1,2,3]:
                axj = fig.add_subplot(gs[0,1+j])
                # plotting components of the arrVectors
                heatmap(axj, arrVectors[:,:,j,i], array_x, array_y, labels, '', log_norm=False, **plot_projection_params)
                axj.set_title( r'projection on '+direction_label[j], fontsize=FS )

            plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + str(i) + '.pdf', bbox_inches='tight')
            plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + str(i) + '.png', bbox_inches='tight')

            plt.close()


    return fig


def fig23_heatmaps(arrDetSigmaEst, arrRelErrorEst, array_x, array_y, fname, labels, log_select):
    """
    arr*: various 2D arrays to used for heatmap
    array_x, array_y: the range for the axis
    fname: filename to save
    labels: the label of the x,y axis (need this to be an array of 2 values, both strings)
    log_select: whether or not to have a log norm for the heatmap
    """
    # makes a figure with many subplots.

    fig = plt.figure(figsize=(5,20));
    ax1 = plt.subplot2grid((4,1), (0,0), colspan=1, rowspan=1);
    ax2 = plt.subplot2grid((4,1), (1,0), colspan=1, rowspan=1);
    ax3 = plt.subplot2grid((4,1), (2,0), colspan=1, rowspan=1);
    ax4 = plt.subplot2grid((4,1), (3,0), colspan=1, rowspan=1);

    # each of the diagonals
    heatmap(ax1, arrRelErrorEst[:,:,0,0].astype(np.float64)*(KP*T), array_x, array_y, labels, r'$k_{p} t \langle \delta {c_1}^2 \rangle / {c_1}^2$', log_norm=log_select)
    heatmap(ax2, arrRelErrorEst[:,:,1,1].astype(np.float64)*(KP*T), array_x, array_y, labels, r'$k_{p} t \langle \delta {k_{off,1}}^2 \rangle / {k_{off,1}}^2$', log_norm=log_select)
    heatmap(ax3, arrRelErrorEst[:,:,2,2].astype(np.float64)*(KP*T), array_x, array_y, labels, r'$k_{p} t \langle \delta {c_2}^2 \rangle / {c_2}^2$', log_norm=log_select)
    heatmap( ax4, arrRelErrorEst[:,:,3,3].astype(np.float64)*(KP*T), array_x, array_y, labels, r'$k_{p} t \langle \delta {k_{off,2}}^2 \rangle / {k_{off,2}}^2$', log_norm=log_select)

    # save
    #fig.suptitle(r"IFN Crosstalk $k_{off}/k_{off,2}=100$", fontsize=int(1.5*FS)) # TODO customize the title
    #fig.suptitle(r"T Cell antigen detection $c_2/c_1=1000$", fontsize=int(1.5*FS))
    plt.tight_layout(h_pad=1, rect=[0,0.03,1,0.95]) # need to change this to not have the title overlap, figuring it out still

    plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.png'); #plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.eps');

    plt.close()

    return fig

def eigens_of_sigmaEst(sigmaEst):
    #value, vector = np.linalg.eig(sigmaEst)
    value, vector = sp.linalg.eig(sigmaEst)

    # the eigenvalues are not in any particular order: I will order them now
    idx = value.argsort()

    # returns ordered eigenvalues, smallest to largest
    return value[idx], vector[:,idx], np.prod(value)



def sigmaEst(c1, koff, c2, koff2, add_trace_term=True):
    # eqns2l are the equations imported from Mathematica and turned into matrices (function of c1,c2,koff,koff2)
    A = eqns2l.matrix_dmudthetaInv(c1, koff, c2, koff2)
    B = eqns2l.matrix_sigmadata(c1, koff, c2, koff2)
    Atrans = np.transpose(A)

    def build_dcov_dtheta_idx(theta_idx):

        label_data = {0: 'N1', 1:'M1', 2:'N2', 3:'M2'}
        label_theta = {0: 'c1', 1:'koff', 2:'c2', 3:'koff2'}

        def fetch_eqn(i, j):
            # TODO this is the text label is it the eqn label though
            # TODO rename sigmadataIJ_dK
            # old
            # new 'sigmadata%s%s_d%s' % (label_data[i], label_data[j], label_theta[theta_idx])
            if j > i:
                method_to_call = getattr(eqns2l,
                                         'Cov%s%sd%s' % (label_data[i], label_data[j], label_theta[theta_idx]))
            else:
                method_to_call = getattr(eqns2l,
                                         'Cov%s%sd%s' % (label_data[j], label_data[i], label_theta[theta_idx]))
            return method_to_call

        dcov_dtheta_idx = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                dcov_dtheta_idx[i, j] = fetch_eqn(i, j)(c1, koff, c2, koff2)

        return dcov_dtheta_idx

    def build_FI_trace_term():
        cov_matrix = eqns2l.matrix_sigmadata(c1, koff, c2, koff2)
        cov_matrix_inv = np.linalg.inv(cov_matrix)
        dcov_dtheta_idx_dict = {0: build_dcov_dtheta_idx(0),
                                1: build_dcov_dtheta_idx(1),
                                2: build_dcov_dtheta_idx(2),
                                3: build_dcov_dtheta_idx(3)}

        def compute_trace(n, m):
            arr = np.linalg.multi_dot([cov_matrix_inv,
                                       dcov_dtheta_idx_dict[n],
                                       cov_matrix_inv,
                                       dcov_dtheta_idx_dict[m]])
            return np.trace(arr)

        FI_trace_term = np.zeros((4, 4))
        for n in range(4):
            for m in range(4):
                FI_trace_term[n, m] = 0.5 * compute_trace(n, m)

        return FI_trace_term

    error_matrix = np.linalg.multi_dot([A, B, Atrans])
    if add_trace_term:
        # idea is to invert the matrix above, add the trace term, then invert the whole thing
        FI_term_base = np.linalg.inv(error_matrix)
        FI_term_trace = build_FI_trace_term()
        FI_full = FI_term_base + FI_term_trace
        error_matrix = np.linalg.inv(FI_full)

    # use to make the entries relative estimates
    rel = np.array([[c1 * c1, c1 * koff, c1 * c2, c1 * koff2], [koff * c1, koff * koff, koff * c2, koff * koff2],
                    [c2 * c1, c2 * koff, c2 * c2, c2 * koff2], [koff2 * c1, koff2 * koff, koff2 * c2, koff2 * koff2]])
    relErrorMatrix = np.divide(error_matrix, rel)

    return error_matrix, np.linalg.det(error_matrix), relErrorMatrix, np.linalg.det(error_matrix)/(c1**2 * c2**2 * koff**2 * koff2**2 )


def dmudthetaInv(c1, koff, c2, koff2):
    """ Create the matrix from exported mathematica expressions """
    return eqns2l.matrix_dmudthetaInv(c1, koff, c2, koff2), np.linalg.det( eqns2l.matrix_dmudthetaInv(c1, koff, c2, koff2) )

def sigmaData(c1, koff, c2, koff2):
    """ Create the matrix from exported mathematica expressions """
    return eqns2l.matrix_sigmadata(c1, koff, c2, koff2), np.linalg.det( eqns2l.matrix_sigmadata(c1, koff, c2, koff2) )


def fig_2ligands_vs_1ligand_diag(arrErrorRatio, array_x, array_y, fname, labels, log_select=True):

    # makes a figure with many subplots.
    fig = plt.figure(figsize=(5, 20))
    ax1 = plt.subplot2grid((4, 1), (0, 0), colspan=1, rowspan=1)
    ax2 = plt.subplot2grid((4, 1), (1, 0), colspan=1, rowspan=1)
    ax3 = plt.subplot2grid((4, 1), (2, 0), colspan=1, rowspan=1)
    ax4 = plt.subplot2grid((4, 1), (3, 0), colspan=1, rowspan=1)

    # each of the diagonals
    heatmap(ax1, arrErrorRatio[:, :, 0, 0].astype(np.float64), array_x, array_y, labels,
            r'$\langle \delta {c_1}^2 \rangle_B / \langle \delta {c}^2 \rangle_{A1}$', log_norm=log_select)
    heatmap(ax2, arrErrorRatio[:, :, 1, 1].astype(np.float64), array_x, array_y, labels,
            r'$\langle \delta {k_{off,1}}^2 \rangle_B / \langle \delta {k_{off}}^2 \rangle_{A1}$', log_norm=log_select)
    heatmap(ax3, arrErrorRatio[:, :, 2, 2].astype(np.float64), array_x, array_y, labels,
            r'$\langle \delta {c_2}^2 \rangle_B / \langle \delta {c}^2 \rangle_{A2}$', log_norm=log_select)
    heatmap(ax4, arrErrorRatio[:, :, 3, 3].astype(np.float64), array_x, array_y, labels,
            r'$\langle \delta {k_{off,2}}^2 \rangle_B / \langle \delta {k_{off}}^2 \rangle_{A2}$', log_norm=log_select)

    # save
    # fig.suptitle(r"IFN Crosstalk $k_{off}/k_{off,2}=100$", fontsize=int(1.5*FS)) # TODO customize the title
    # fig.suptitle(r"T Cell antigen detection $c_2/c_1=1000$", fontsize=int(1.5*FS))
    plt.tight_layout(h_pad=1, rect=[0, 0.03, 1,
                                    0.95])  # need to change this to not have the title overlap, figuring it out still

    plt.savefig(DIR_OUTPUT + os.sep + 'ligands2_over_ligands1' + os.sep + fname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + 'ligands2_over_ligands1' + os.sep + fname + '.png')

    plt.close()

    return fig


def fig_2ligands_vs_1ligand_det(arrDetErrorRatio, array_x, array_y, fname, labels, log_select=True):
    # makes a figure with many subplots.
    fig = plt.figure(figsize=(5, 5))
    ax = plt.gca()

    # each of the diagonals
    heatmap(ax, arrDetErrorRatio[:, :].astype(np.float64), array_x, array_y, labels,
            r'$ \frac{\mathrm{det} \Sigma_B}{\mathrm{det} \Sigma_{A1} \mathrm{det} \Sigma_{A2}} $', log_norm=log_select)#, vmax=10**5)

    # save
    # fig.suptitle(r"IFN Crosstalk $k_{off}/k_{off,2}=100$", fontsize=int(1.5*FS)) # TODO customize the title
    # fig.suptitle(r"T Cell antigen detection $c_2/c_1=1000$", fontsize=int(1.5*FS))
    #plt.tight_layout(h_pad=1, rect=[0, 0.03, 1,
    #                                0.95])  # need to change this to not have the title overlap, figuring it out still

    plt.savefig(DIR_OUTPUT + os.sep + 'ligands2_over_ligands1' + os.sep + fname + '_det' + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + 'ligands2_over_ligands1' + os.sep + fname + '_det' + '.png')

    plt.close()

    return fig


if __name__ == '__main__':

    for dirs in ['ligands2', 'ligands2_over_ligands1']:
        if not os.path.exists(DIR_OUTPUT + os.sep + dirs):
            os.makedirs(DIR_OUTPUT + os.sep + dirs)

    flag_general = True
    flag_compare_2ligands_vs_1ligand = False
    ADD_TRACE_TERM = True
    LOG_SELECT = True

    # choose any of these 2 to be arrays [0: c1, 1: koff, 2: c2, 3: koff2], they will be your axes
    dim = {'x' : 2, 'y' : 3}

    # axis we'd like to plot
    value = [C1, KOFF, C2, KOFF2]
    label_fig = ['c1', 'koff1', 'c2', 'koff2']
    label = [r'$c_1$', r'$k_{\mathrm{off},1}$', r'$c_2$', r'$k_{\mathrm{off},2}$']
    dimension = [CRANGE, KOFFRANGE, CRANGE, KOFFRANGE2]
    axes = label_fig[dim['x']]+ '_' + label_fig[dim['y']]

    # for plotting the dedimensionalized values
    dedimension = [dimension[0] * KON / KP, dimension[1] / KP, dimension[2] * KON / KP, dimension[3] / KP]
    dedimension_label = [r'$k_{\mathrm{on}}c_1/k_{p}$', r'$k_{\mathrm{off},1}/k_{p}$', r'$k_{\mathrm{on}}c_2/k_{p}$', r'$k_{\mathrm{off},2}/k_{p}$']
    # ratio
    # dedimension = [dimension[0]*KON/KP, dimension[1]/KP, dimension[2]*KON/KP, dimension[3]/dimension[1]]
    # dedimension_label = [r'$k_{on}c_1/k_{p}$', r'$k_{off,1}/k_{p}$', r'$k_{on}c_2/k_{p}$', r'$k_{off,2}/k_{off,1}$']
    # different
    # dedimension = [dimension[0]*KON/KP, dimension[1]/KP, dimension[2]*KON/KP, (dimension[3]-dimension[1])/KP]
    # dedimension_label = [r'$k_{on}c_1/k_{p}$', r'$k_{off,1}/k_{p}$', r'$k_{on}c_2/k_{p}$', r'$(k_{off,2}-k_{off,1})/k_{p}$']

    # heatmap for any element of sigmaEst(4x4 matrix)
    arrSigmaEst = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4) )
    arrRelErrorEst = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4) )
    arrSigmaData = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4) )
    arrDmudthetaInv = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4) )
    arrEigenvalues = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]), 4) )
    arrEigenvectors = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4) )

    # heatmap for the det estimate covariance
    arrDetSigmaEst = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )
    arrDetEigenvalue = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )
    arrDetSigmaData = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )
    arrDetDmudthetaInv = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )
    arrRelDetSigmaEst = np.zeros( (len(dimension[dim['y']]), len(dimension[dim['x']]) ) )

    # making our heatmap
    num_heatmap_points = len(dimension[dim['x']]) * len(dimension[dim['y']])
    pixel_num = 0.0
    progress_update_threshold = 5.0
    print("Computing matrices")
    for i, xi in enumerate(dimension[dim['x']]):
        value[dim['x']] = xi
        for j, yi in enumerate(dimension[dim['y']]):
            # progress update
            pixel_num += 1.0
            if 100*pixel_num/num_heatmap_points > progress_update_threshold:
                print(str(progress_update_threshold) + "%")
                progress_update_threshold += 5.0
            # matrix calculations
            value[dim['y']] = yi
            arrSigmaEst[j, i, :, :], arrDetSigmaEst[j, i], arrRelErrorEst[j, i, :, :], arrRelDetSigmaEst[j,i] = \
                sigmaEst(value[0], value[1], value[2], value[3], add_trace_term=ADD_TRACE_TERM)
            arrEigenvalues[j,i,:], arrEigenvectors[j,i,:,:], arrDetEigenvalue[j,i] = eigens_of_sigmaEst( arrSigmaEst[j,i,:,:] )
    print("Done computing matrices\n Plotting results")
    ver = "new_"
    multi_fname = "multi_" + ver + axes

    if flag_general:
        # Here I have selected the plots I want, when we are making individual heatmaps. LOG_SELECT is kinda old, used to be for when we didn't want log colorbar, but we nearly always want it now

        # make a Fig3 from paper
        """
        multiple_heatmaps(arrRelDetSigmaEst, arrRelErrorEst, dedimension[dim['x']], dedimension[dim['y']], multi_fname,
                          [dedimension_label[dim['x']], dedimension_label[dim['y']]], LOG_SELECT)
        """

        # Eigenvalue plots
        #eigen_heatmaps( arrEigenvalues, arrEigenvectors, dedimension[dim['x']], dedimension[dim['y']], 'eigen_plots', [dedimension_label[dim['x']], dedimension_label[dim['y']]], label)
        single_heatmap( arrDetEigenvalue, dedimension[dim['x']], dedimension[dim['y']], 'det_eigen_plots', [dedimension_label[dim['x']], dedimension_label[dim['y']]], 'determinant', log_norm=True)
        """
        # covariance
        multi_dcov_fname = "multi_jacob_"+ver+axes
        multiple_heatmaps(arrDetDmudthetaInv, arrDmudthetaInv,  dedimension[dim['x']],  dedimension[dim['y']], multi_dcov_fname, [ dedimension_label[dim['x']], dedimension_label[dim['y']]], LOG_SELECT)

        # Jacobian
        multi_jacob_fname = "multi_dcov_"+ver+axes
        multiple_heatmaps(arrDetSigmaData, arrSigmaData,  dedimension[dim['x']],  dedimension[dim['y']], multi_jacob_fname, [ dedimension_label[dim['x']], dedimension_label[dim['y']]], False)

        # Jacobian
        multi_jacob_fname = "multi_dcovlog_"+ver+axes
        multiple_heatmaps(arrDetSigmaData, arrSigmaData,  dedimension[dim['x']],  dedimension[dim['y']], multi_jacob_fname, [ dedimension_label[dim['x']], dedimension_label[dim['y']]], LOG_SELECT)
        """

        # old plots to check if conversion from matehmatica worked fine
        # plot_heatmap(arrDetdmudthetaInv[:,:], dedimension[dim['x']], dedimension[dim['y']], 'c1c2detdmudtheta', [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Determinant $d\mu /d \theta^{-1}$')
        # plot_heatmap(element11SigmaData[:,:], dedimension[dim['x']], dedimension[dim['y']], 'c1c2sigmadata11', [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Determinant $\Sigma_{data,11}$')
        # plot_heatmap(element31SigmaData[:,:], dedimension[dim['x']], dedimension[dim['y']], 'c1c2sigmadata31', [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Determinant $\Sigma_{data,31}$')
        # plot_heatmap(element11dmudthetaInv[:,:], dedimension[dim['x']], dedimension[dim['y']], 'c1c2detdmudtheta11', [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Determinant $d\mu /d \theta^{-1}_{11}$')
        # plot_heatmap(element31dmudthetaInv[:,:], dedimension[dim['x']], dedimension[dim['y']], 'c1c2detdmudtheta31', [ dedimension_label[dim['x']], dedimension_label[dim['y']]], r'Determinant $d\mu /d \theta^{-1}_{31}$')

    if flag_compare_2ligands_vs_1ligand:
        """
        From Jeremy pre-April 27 code we create the following arrays in main:
        size = (4 x 4 x axis1 x axis2)
            - arrSigmaEst
            - arrRelErrorEst
        size = (1 x axis1 x axis2)
            - arrDetSigmaEst
            - arrDetSigmaData
            - arrDetdmudthetaInv

        We want to compare these arrays to their 1 ligand counterparts
        """

        rel_err_model_A1_c = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))
        rel_err_model_A1_koff = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))
        rel_err_model_A2_c = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))
        rel_err_model_A2_koff = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))
        det_err_A1 = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))
        det_err_A2 = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))

        for i, xval in enumerate(dimension[dim['x']]):
            value[dim['x']] = xval
            for j, yval in enumerate(dimension[dim['y']]):
                value[dim['y']] = yval

                rel_err_model_A1_c[j, i] = eqns.RelErrC2NoTrace(value[0], value[1], kon=KON, T=T, KP=KP, KF=None)
                rel_err_model_A1_koff[j, i] = eqns.RelErrK2NoTrace(value[0], value[1], kon=KON, T=T, KP=KP, KF=None)
                rel_err_model_A2_c[j, i] = eqns.RelErrC2NoTrace(value[2], value[3], kon=KON, T=T, KP=KP, KF=None)
                rel_err_model_A2_koff[j, i] = eqns.RelErrK2NoTrace(value[2], value[3], kon=KON, T=T, KP=KP, KF=None)

                det_err_A1[j, i] = eqns.DetSigmacrlb2NoTrace(value[0], value[1], kon=KON, T=T, KP=KP, KF=None)
                det_err_A2[j, i] = eqns.DetSigmacrlb2NoTrace(value[2], value[3], kon=KON, T=T, KP=KP, KF=None)

        arrErrorRatio = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4))
        arrDetErrorRatio = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))

        for i, xval in enumerate(dimension[dim['x']]):
            value[dim['x']] = xval
            for j, yval in enumerate(dimension[dim['y']]):
                value[dim['y']] = yval
                arrErrorRatio[j, i, 0, 0] = arrRelErrorEst[j, i, 0, 0] / rel_err_model_A1_c[j, i]
                arrErrorRatio[j, i, 1, 1] = arrRelErrorEst[j, i, 1, 1] / rel_err_model_A1_koff[j, i]
                arrErrorRatio[j, i, 2, 2] = arrRelErrorEst[j, i, 2, 2] / rel_err_model_A2_c[j, i]
                arrErrorRatio[j, i, 3, 3] = arrRelErrorEst[j, i, 3, 3] / rel_err_model_A2_koff[j, i]
                arrDetErrorRatio[j, i] = arrDetSigmaEst[j, i] / (det_err_A1[j, i] * det_err_A2[j, i])

        axis_labels = [dedimension_label[dim['x']], dedimension_label[dim['y']]]
        fig_2ligands_vs_1ligand_diag(arrErrorRatio, dedimension[dim['x']], dedimension[dim['y']], multi_fname,
                                     axis_labels, log_select=True)
        fig_2ligands_vs_1ligand_det(arrDetErrorRatio, dedimension[dim['x']], dedimension[dim['y']], multi_fname,
                                    axis_labels, log_select=True)
