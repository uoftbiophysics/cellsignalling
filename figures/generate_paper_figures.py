import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
from decimal import Decimal

import matplotlib.ticker as ticker
import equations2ligands as eqns2l
import equations as eqns

from settings import DIR_OUTPUT, DIR_INPUT, KON, KP, T, KF, ALPHA, C1, C2, KOFF, KOFF2, N, COLOR_SCHEME

plt.style.use('parameters.mplstyle')


# plot params
FS = 10
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
        mpl.colors.LogNorm.__init__(self, vmin=vmin, vmax=vmax, clip=clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [np.log(self.vmin), np.log(self.midpoint), np.log(self.vmax)], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(np.log(value), x, y))


def __sigmaEst__(c1, koff, c2, koff2, add_trace_term=True):
    # eqns2l are the equations imported from Mathematica and turned into matrices (function of c1,c2,koff,koff2)
    A = eqns2l.matrix_dmudthetaInv(c1, koff, c2, koff2)
    B = eqns2l.matrix_sigmadata(c1, koff, c2, koff2)
    Atrans = np.transpose(A)

    def build_dcov_dtheta_idx(theta_idx):

        label_data = {0: 'N1', 1: 'M1', 2: 'N2', 3: 'M2'}
        label_theta = {0: 'c1', 1: 'koff', 2: 'c2', 3: 'koff2'}

        def fetch_eqn(i, j):
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

    return error_matrix, np.linalg.det(error_matrix), relErrorMatrix, np.linalg.det(error_matrix)/(c1**2 * c2**2 * koff**2 * koff2**2)


def __dmudthetaInv__(c1, koff, c2, koff2):
    """ Create the matrix from exported mathematica expressions """
    return eqns2l.matrix_dmudthetaInv(c1, koff, c2, koff2), np.linalg.det(eqns2l.matrix_dmudthetaInv(c1, koff, c2, koff2))


def __sigmaData__(c1, koff, c2, koff2):
    """ Create the matrix from exported mathematica expressions """
    return eqns2l.matrix_sigmadata(c1, koff, c2, koff2), np.linalg.det(eqns2l.matrix_sigmadata(c1, koff, c2, koff2))


def __heatmap__(ax, arr, xrange, yrange, xy_label, label, log_norm=True, less_xticks=False,
                skip_cbar=False, cbar_white_loc=1.0, **kwargs):
    """
    xrange: range of values for x
    yrange: range of y values
    xy_label: list with first element the x label (str) and second element the
    y label (str)
    fname: file saved as pdf and eps format with this name
    show: shows plots
    save: saves plot
    log_norm: heatmap is log, this is only used when plotting the posterior
    kwargs: a variety of keyword args are used below. They are mostly used for
    contour plot lines
    """
    if 'levels' in kwargs.keys():
        levels = kwargs['levels']
    else:
        levels = [1E0]  # , 1E5, 1E10, 1E20]

    if 'vmin' in kwargs.keys():
        vmin = kwargs['vmin']
    else:
        vmin = 1E-3  # np.min(arr)

    if 'vmax' in kwargs.keys():
        vmax = kwargs['vmax']
    else:
        if np.max(arr) > 1E10:
            vmax = 1E5
        else:
            vmax = 1E5  # np.max(arr)

    if 'contour_linestyle' in kwargs.keys():
        contour_linestyle = kwargs['contour_linestyle']
    else:
        contour_linestyle = '-'

    if 'contour_color' in kwargs.keys():
        contour_color = kwargs['contour_color']
    else:
        contour_color = 'k'

    if 'contour_linewidths' in kwargs.keys():
        contour_lindewidths = kwargs['contour_linewidths']
    else:
        contour_lindewidths = 1.2

    if 'cmap_colour' in kwargs.keys():
        cmap_colour = kwargs['cmap_colour']
    else:
        cmap_colour = 'RdBu_r'  # 'YlGnBu' or 'RdBu_r'

    if 'fmt' in kwargs.keys():
        fmt = kwargs['fmt']
    else:
        fmt = ticker.LogFormatterMathtext()

    if log_norm:
        arr = np.ma.array(arr, mask=(arr <= 0.0))
        if cmap_colour == 'RdBu_r':
            imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax, 'norm': MidPointLogNorm(vmin, vmax, cbar_white_loc)}
        else:
            imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax, 'norm': mpl.colors.LogNorm(vmin, vmax)}
        print("Logscale")
    else:
        imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax}

    """
    Colours viridis, YlGnBu, terrain, plasma
    """

    im = ax.imshow(arr, interpolation='none', **imshow_kw)
    # im = plt.imshow(arr, interpolation='spline36', **imshow_kw)

    # axes setup

    # axes log scaled
    if less_xticks:
        ax.set_xticks([i for i, xval in enumerate(xrange) if (i % POINTS_BETWEEN_TICKS == 0 and np.log10(xval) % 2 == 0)])
        ax.set_xticklabels([r'$10^{%d}$' % np.log10(xval) for i, xval in enumerate(xrange) if (i % POINTS_BETWEEN_TICKS == 0 and np.log10(xval) % 2 == 0)], fontsize=FS)
    else:
        ax.set_xticks([i for i, cval in enumerate(xrange) if i % POINTS_BETWEEN_TICKS == 0])
        ax.set_xticklabels([r'$10^{%d}$' % np.log10(xval) for i, xval in enumerate(xrange) if (i % POINTS_BETWEEN_TICKS == 0)], fontsize=FS)
    ax.set_yticks([i for i, kval in enumerate(yrange) if i % POINTS_BETWEEN_TICKS == 0])
    ax.set_yticklabels([r'$10^{%d}$' % np.log10(yval) for i, yval in enumerate(yrange) if i % POINTS_BETWEEN_TICKS == 0], fontsize=FS)

    ax.invert_yaxis()
    ax.set_xlabel(xy_label[0], fontsize=FS)
    ax.set_ylabel(xy_label[1], fontsize=FS)

    # create colorbar
    if skip_cbar:
        cbar = None
    else:
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.05)  # Fig 3 way pad 0.04

        cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=FS, labelpad=-20)
        cbar.ax.tick_params(labelsize=(FS), width=1, length=1)
        cbar.ax.minorticks_off()
        # cbar.set_ticks([round(vmin,3)+0.001,round(vmax,3)-0.001]) # UNCOMMENT THIS ONLY WHEN TICKS DON'T APPEAR
        cbar.update_ticks()

    CL = ax.contour(arr, levels=levels, linestyles=contour_linestyle, colors=contour_color, linewidths=contour_lindewidths)
    ax.clabel(CL, CL.levels, inline=True, fmt=fmt, fontsize=FS-2)

    plt.tight_layout(h_pad=1, w_pad=1)

    return ax, cbar, im


def __heatmap2__(arr, crange, koffrange, fname, label, show=SHOW, save=True, log_norm=True, log_norm_axes_only=False,
                 dedim=False, xy_label_force=None, less_xticks=True, return_cb=False, **kwargs):
    """
    crange: range of values for c
    koffrange: range of koff values
    fname: file saved as pdf and eps format with this name
    show: shows plots
    save: saves plots
    log_norm: heatmap is log, this is only used when plotting the posterior
    dedim: makes the axis dedimensionalized. scales them (see by how much below)
    kwargs: a variety of keyword args are used below. They are mostly used for contour plot lines
    """
    # default parameters
    if 'levels' in kwargs.keys():
        levels = kwargs['levels']
    else:
        levels = [1E0, 1E1, 1E2, 1E3]

    if 'vmin' in kwargs.keys():
        vmin = kwargs['vmin']
    else:
        vmin = np.min(arr)

    if 'vmax' in kwargs.keys():
        vmax = kwargs['vmax']
    else:
        vmax = np.max(arr)

    if 'contour_linestyle' in kwargs.keys():
        contour_linestyle = kwargs['contour_linestyle']
    else:
        contour_linestyle = '-'

    if 'contour_color' in kwargs.keys():
        contour_color = kwargs['contour_color']
    else:
        contour_color = 'k'

    if 'contour_linewidths' in kwargs.keys():
        contour_lindewidths = kwargs['contour_linewidths']
    else:
        contour_lindewidths = 2

    if 'cmap_colour' in kwargs.keys():
        cmap_colour = kwargs['cmap_colour']
    else:
        cmap_colour = 'YlGnBu'

    if 'fmt' in kwargs.keys():
        fmt = kwargs['fmt']
    else:
        fmt = ticker.LogFormatterMathtext()

    if 'fig_width' in kwargs.keys():
        fig_width = kwargs['fig_width']
    else:
        fig_width = 3.0

    if dedim:
        # if flag is true, this is how we scale the axis. Simple.
        crange = crange*KON/KP
        koffrange = koffrange/KP
        xy_label = [r'$k_{on}c/k_{p}$', r'$k_{off}/k_{p}$']
    else:
        xy_label = [r'${c}$', r'${k}_{off}$']
    if xy_label_force is not None:
        xy_label = xy_label_force

    # we want most of the functionality that comes with log_norm, but not the log colorbar
    if log_norm and not log_norm_axes_only:
        imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax, 'norm': mpl.colors.LogNorm(vmin, vmax)}
    else:
        imshow_kw = {'cmap': cmap_colour, 'aspect': None, 'vmin': vmin, 'vmax': vmax}

    """
    Colours viridis, YlGnBu, terrain, plasma
    """
    # plot setup
    f = plt.figure(figsize=(fig_width, fig_width/1.2))
    im = plt.imshow(arr, interpolation='spline36', **imshow_kw)

    # axes setup
    fig = plt.gcf()
    ax = plt.gca()

    # axes log scaled
    if less_xticks:
        ax.set_xticks([i for i, cval in enumerate(crange) if (i % POINTS_BETWEEN_TICKS == 0 and np.log10(cval) % 2 == 0)])
        ax.set_xticklabels([r'$10^{%d}$' % np.log10(cval) for i, cval in enumerate(crange) if (i % POINTS_BETWEEN_TICKS == 0 and np.log10(cval) % 2 == 0)], fontsize=FS)
    else:
        ax.set_xticks([i for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS == 0])
        ax.set_xticklabels([r'$10^{%d}$' % np.log10(cval) for i, cval in enumerate(crange) if (i % POINTS_BETWEEN_TICKS == 0)], fontsize=FS)
    ax.set_yticks([i for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0])
    ax.set_yticklabels([r'$10^{%d}$' % np.log10(yval) for i, yval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0], fontsize=FS)

    ax.invert_yaxis()
    ax.set_xlabel(xy_label[0], fontsize=FS)
    ax.set_ylabel(xy_label[1], fontsize=FS)

    # create colorbar
    cbar = fig.colorbar(im, fraction=0.0375, pad=0.04)
    cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=FS, labelpad=5)
    cbar.ax.tick_params(labelsize=FS)
    cbar.ax.minorticks_off()
    # UNCOMMENT THIS ONLY WHEN TICKS DON'T APPEAR
    # cbar.set_ticks([round(vmin,3)+0.001,round(vmax,3)-0.001])
    cbar.update_ticks()
    cbar.ax.minorticks_off()

    CL = plt.contour(arr, levels=levels, linestyles=contour_linestyle, colors=contour_color, linewidths=contour_lindewidths)
    plt.clabel(CL, CL.levels, inline=True, fmt=fmt)

    # now we can set the ticks without them getting messed up
    if not log_norm:
        ax = plt.gca()
        ax.set_xticks([i for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS == 0] + [len(crange)])
        ax.set_yticks([i for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0] + [len(koffrange)])

        # Exhausting process to make ticks look nice
        nice_ticks = ["{}".format(cval) for cval in crange if Decimal(str(cval)) % Decimal(str(crange[-1] / 5)) == 0]
        ctilde_ticks = []
        for i in range(len(nice_ticks) * 2):
            if i % 2 == 0:
                ctilde_ticks.append("")
            else:
                ctilde_ticks.append(nice_ticks.pop(0))
        ctilde_ticks = ["0"] + ctilde_ticks

        nice_ticks = ["{}".format(koffval) for koffval in koffrange if
                      Decimal(str(koffval)) % Decimal(str(koffrange[-1] / 5)) == 0]
        kofftilde_ticks = []
        for i in range(len(nice_ticks) * 2):
            if i % 2 == 0:
                kofftilde_ticks.append("")
            else:
                kofftilde_ticks.append(nice_ticks.pop(0))
        kofftilde_ticks = ["0"] + kofftilde_ticks

        ax.set_xticklabels(ctilde_ticks, fontsize=FS)
        ax.set_yticklabels(kofftilde_ticks, fontsize=FS)

    # save
    if save:
        plt.savefig(DIR_OUTPUT + os.sep + fname + '.pdf')
        plt.savefig(DIR_OUTPUT + os.sep + fname + '.eps')
        plt.savefig(DIR_OUTPUT + os.sep + fname + '.png')
    if show:
        plt.show()

    plt.close()
    if return_cb:
        return fig, ax, cbar
    else:
        return fig, ax


def __heatmap3__(equation, label, filename, log_norm, crange=CRANGE, koffrange=KOFFRANGE, dedim=False,
                 contour_args={'None': None}, verbose=False):
    if verbose:
        print("Plotting equation:", equation)
    # create array to plot
    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = equation(cval, koffval)+10**(-10)

    # call heatmap plotting
    __heatmap2__(arr, crange, koffrange, filename, label, dedim=dedim, log_norm=log_norm, **contour_args)
    return 0


def __heatmap_ratio__(eqn1, eqn2, label, filename, log_norm, crange=CRANGE, koffrange=KOFFRANGE, dedim=False, contour_args=None):
    # creating arrays to plot
    arr = np.zeros((len(koffrange), len(crange)))
    arr1 = np.zeros((len(koffrange), len(crange)))
    arr2 = np.zeros((len(koffrange), len(crange)))
    # filling arrays with ratio
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = eqn1(cval, koffval)/eqn2(cval, koffval)
            arr1[i, j] = eqn1(cval, koffval)  # these are unnecessary, however good for checks
            arr2[i, j] = eqn2(cval, koffval)

    __heatmap2__(arr, crange, koffrange, filename, label, log_norm=log_norm, dedim=dedim, **contour_args)


def figure_1_and_2_heatmaps(arrRelErrorEst1X, arrRelErrorEst2C, arrRelErrorEst2KOFF, array_x, array_y, fname1, fname2, fname3, labels, log_select):
    """
    Makes figures 1D, 2C,D for publication.
    """

    arr1 = arrRelErrorEst1X/N
    args1 = {'vmin': np.min(arr1), 'vmax': np.max(arr1), 'levels': [1E-1, 1E0, 1E1, 1E2], 'less_xticks': True}
    arr2 = arrRelErrorEst2C/N
    args2 = {'vmin': np.min(arr2), 'vmax': np.max(arr2), 'levels': [1E-1, 1E0, 1E1, 1E2], 'less_xticks': True}
    arr3 = arrRelErrorEst2KOFF/N
    args3 = {'vmin': np.min(arr3), 'vmax': np.max(arr3), 'levels': [1E-1, 1E0, 1E1, 1E2], 'less_xticks': True}

    fig0 = plt.figure(figsize=(2.8, 2.2))
    fig0 = plt.gcf()
    ax0 = plt.gca()
    # axes setup
    ax0, cbar0, im0 = __heatmap__(ax0, arr1, array_x, array_y, labels, '', log_norm=True, skip_cbar=False, cbar_white_loc=10, **args1)
    ax0.set_title(r'$\frac{k_p t}{N} \frac{\langle\delta x^{2}\rangle}{x^{2}}$', fontsize=FS)

    if not os.path.exists(DIR_OUTPUT + os.sep + 'ligand1'):
        os.mkdir(DIR_OUTPUT + os.sep + 'ligand1')
    plt.savefig(DIR_OUTPUT + os.sep + 'ligand1' + os.sep + fname1 + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + 'ligand1' + os.sep + fname1 + '.png')  #plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.eps');
    plt.close()

    fig1 = plt.figure(figsize=(2.8, 2.2))
    fig1 = plt.gcf()
    ax1 = plt.gca()
    # axes setup
    ax1, cbar1, im1 = __heatmap__(ax1, arr2, array_x, array_y, labels, r'', log_norm=True, skip_cbar=False, cbar_white_loc=10, **args2)
    ax1.set_title(r'$\frac{k_p t}{N} \frac{\langle\delta c^{2}\rangle}{c^{2}}$', fontsize=FS)
    plt.savefig(DIR_OUTPUT + os.sep + 'ligand1' + os.sep + fname2 + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + 'ligand1' + os.sep + fname2 + '.png')  # plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.eps');
    plt.close()

    fig2 = plt.figure(figsize=(2.8, 2.2))
    fig2 = plt.gcf()
    ax2 = plt.gca()
    # axes setup
    ax2, cbar2, im2 = __heatmap__(ax2, arr3, array_x, array_y, labels, r'', log_norm=True, skip_cbar=False, cbar_white_loc=10, **args3)
    ax2.set_title(r'$\frac{k_p t}{N} \frac{\langle{\delta k_{off}^2}\rangle}{k_{off}^{2}}$', fontsize=FS)
    plt.savefig(DIR_OUTPUT + os.sep + 'ligand1' + os.sep + fname3 + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + 'ligand1' + os.sep + fname3 + '.png')  # plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.eps');
    plt.close()

    return fig0, fig1, fig2


def figure_3_heatmaps(arrRelDetSigmaEst, arrRelErrorEst, array_x, array_y, fname, labels, log_select, rescale=True, axis_label=True):
    """
    arr*: various 2D arrays to used for heatmap
    array_x, array_y: the range for the axis
    fname: filename to save
    labels: the label of the x,y axis (need this to be an array of 2 values, both strings)
    log_select: whether or not to have a log norm for the heatmap
    """

    # makes a figure with many subplots.
    layout_horizontal = False
    det_on_left = False

    if not rescale:
        if layout_horizontal:
            fig = plt.figure(figsize=(7, 3.2))
        else:
            fig = plt.figure(figsize=(4.6, 8.4))

        if layout_horizontal and det_on_left:
            gs = fig.add_gridspec(2, 6, hspace=-0.14, wspace=0.05, width_ratios=[1., 1., 0.08, 1., 1., 0.1], height_ratios=[1., 1.])

            ax0 = fig.add_subplot(gs[:, :2])  # the det
            ax1 = fig.add_subplot(gs[0, 3])
            ax2 = fig.add_subplot(gs[0, 4])
            ax3 = fig.add_subplot(gs[1, 3])
            ax4 = fig.add_subplot(gs[1, 4])
            ax5 = fig.add_subplot(gs[:, -1])   # the cbar

        elif layout_horizontal:
            gs = fig.add_gridspec(2, 6, hspace=0.0, wspace=0.05, width_ratios=[1., 1., 0.4, 1., 1., 0.1], height_ratios=[1., 1.])

            ax0 = fig.add_subplot(gs[:, 3:5])  # the det
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[0, 1])
            ax3 = fig.add_subplot(gs[1, 0])
            ax4 = fig.add_subplot(gs[1, 1])
            ax5 = fig.add_subplot(gs[:, -1])   # the cbar

        else:
            gs = fig.add_gridspec(8, 8, hspace=0.2, wspace=-0.1,
                                  width_ratios=[0.4, 0.15, 0.4, 1.2, 1.2, 0.4, 0.15, 0.4],
                                  height_ratios=[1., 0.07, 1., 0.01, 0.1, 1., 1, 0.1])
            ax0 = fig.add_subplot(gs[4:, 1:5])  # the det
            ax1 = fig.add_subplot(gs[0, 0:4])
            ax2 = fig.add_subplot(gs[0, 4:8])
            ax3 = fig.add_subplot(gs[2, 0:4])
            ax4 = fig.add_subplot(gs[2, 4:8])
            ax5 = fig.add_subplot(gs[5:7, -2])  # the cbar (det)

        diag_cbar_white_loc = 1.0
        det_cbar_white_loc = 1.0
        diag_kw = {'vmin': 1E-3, 'vmax': 1E5, 'levels': [1E0]}
        det_kw = {'vmin': 1E-3, 'vmax': 1E5, 'levels': [1E0]}
        if KP * T < 100:
            print('Warning: these figure_3_heatmaps(...) settings chosen with kp*t=10^3 in mind')

    else:
        if not axis_label:
            fig = plt.figure(figsize=(4.6, 8.4))  # Layout without all axes in diagonal heatmaps
            gs = fig.add_gridspec(8, 9, hspace=0.2, wspace=-0.1,
                                  width_ratios=[0.15, 0.4, 0.4, 1.2, 1.2, 0.4, 0.4, 0.15, 0.15],
                                  height_ratios=[1., 0.07, 1., 0.01, 0.1, 1., 1., 0.1])

        else:
            fig = plt.figure(figsize=(5.5, 8.5))
            gs = fig.add_gridspec(8, 9, hspace=.5, wspace=2.0,
                                  width_ratios=[0.05, 0.4, 0.4, 1.2, 1.2, 0.4, 0.4, 0.2, 0.2],
                                  height_ratios=[1., 0.1, 1., -0.5, 0.4, 0.4, 1., 0.4])
        ax0 = fig.add_subplot(gs[4:, 2:6])  # the det
        ax1 = fig.add_subplot(gs[0, 0:4])
        ax2 = fig.add_subplot(gs[0, 4:8])
        ax3 = fig.add_subplot(gs[2, 0:4])
        ax4 = fig.add_subplot(gs[2, 4:8])
        ax5 = fig.add_subplot(gs[5:7, -2])  # the cbar (det) (5:7)
        ax6 = fig.add_subplot(gs[0:3, -1])  # the cbar (4)

        N = 1E2  # the number of the receptors
        scale_error_elements = KP * T / N
        arrRelDetSigmaEst = scale_error_elements ** 4 * arrRelDetSigmaEst
        arrRelErrorEst = scale_error_elements * arrRelErrorEst
        diag_min = np.min([np.min(arrRelErrorEst[:, :, idx, idx]) for idx in range(4)])
        diag_max = np.max([np.max(arrRelErrorEst[:, :, idx, idx]) for idx in range(4)])
        print([np.min(arrRelErrorEst[:, :, idx, idx]) for idx in range(4)])

        diag_cbar_white_loc = 10
        diag_max_override = 1e6
        det_cbar_white_loc = diag_cbar_white_loc ** 4
        det_max_override = 1e12  # np.max(arrRelDetSigmaEst[:,:])
        det_min_override = np.min(arrRelDetSigmaEst[:, :])  # 1e-2  # np.min(arrRelDetSigmaEst[:,:])

        diag_levels = [1E-1, 1E0, 1E1, 1E2]
        diag_kw = {'vmin': diag_min, 'vmax': diag_max_override, 'levels': diag_levels}
        det_kw = {'vmin': det_min_override, 'vmax': det_max_override, 'levels': [a**4 for a in diag_levels]}

    if rescale:
        scalelabel = r'$\frac{k_p t}{N}$'
    else:
        scalelabel = ''

    # determinant plot
    ax0, cbar0, im0 = __heatmap__(ax0, arrRelDetSigmaEst[:, :], array_x, array_y, labels, r'det($\Sigma_{est}$)', log_norm=True, skip_cbar=True, cbar_white_loc=det_cbar_white_loc, **det_kw)
    ax0.set_title(r'det(%s$\mathbf{\Sigma}_{est}$)' % scalelabel, fontsize=FS)

    # each of the diagonals
    ax1, cbar1, im1 = __heatmap__(ax1, arrRelErrorEst[:, :, 0, 0], array_x, array_y, labels, r'$\langle \delta {c_1}^2 \rangle / {c_1}^2$', less_xticks=True, log_norm=log_select, skip_cbar=True, cbar_white_loc=diag_cbar_white_loc, **diag_kw)
    ax1.set_title(r'%s$\frac{\langle \delta {c_1}^2 \rangle}{{c_1}^2}$' % scalelabel, fontsize=FS)

    a2, cbar2, _ = __heatmap__(ax2, arrRelErrorEst[:, :, 1, 1], array_x, array_y, labels, r'$\langle \delta {k_{off},1}^2 \rangle / {k_{off},1}^2$', less_xticks=True, log_norm=log_select, skip_cbar=True, cbar_white_loc=diag_cbar_white_loc, **diag_kw)
    ax2.set_title(r'%s$\frac{\langle \delta {k_{\mathrm{off},1}}^2 \rangle}{{k_{\mathrm{off},1}}^2} $' % scalelabel, fontsize=FS)

    a3, cbar3, _ = __heatmap__(ax3, arrRelErrorEst[:, :, 2, 2], array_x, array_y, labels, r'$\langle \delta {c_2}^2 \rangle / {c_2}^2$', less_xticks=True, log_norm=log_select, skip_cbar=True, cbar_white_loc=diag_cbar_white_loc, **diag_kw)
    ax3.set_title(r'%s$\frac{\langle \delta {c_2}^2 \rangle}{{c_2}^2}$' % scalelabel, fontsize=FS)

    a4, cbar4, _ = __heatmap__(ax4, arrRelErrorEst[:, :, 3, 3], array_x, array_y, labels, r'$\langle \delta {k_{off,2}}^2 \rangle / {k_{off,2}}^2$', less_xticks=True, log_norm=log_select, skip_cbar=True, cbar_white_loc=diag_cbar_white_loc, **diag_kw)
    ax4.set_title(r'%s$\frac{\langle \delta {k_{\mathrm{off},2}}^2 \rangle}{{k_{\mathrm{off},2}}^2}$' % scalelabel, fontsize=FS)

    cb_det = fig.colorbar(im0, cax=ax5, fraction=1.0)
    cb_det.ax.tick_params(labelsize=FS)

    if rescale:
        cb_diags = fig.colorbar(im1, cax=ax6, fraction=1.0)
        cb_diags.ax.tick_params(labelsize=FS)
        cb_diags.minorticks_off()

    # note tight layout seems incompatible with gridspec/subplots
    plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.pdf', bbox_inches='tight')
    plt.savefig(DIR_OUTPUT + os.sep + 'ligands2' + os.sep + fname + '.png', bbox_inches='tight')

    plt.close()

    return fig


def plot_1F():
    figname = 'heatmap_log_posterior_mode1'

    nobs = 0.1 * KP * T  # 10

    LOG_START_KOFF_1F = -1
    LOG_END_KOFF_1F = 3
    TOTAL_POINTS_KOFF_1F = (LOG_END_KOFF_1F - LOG_START_KOFF_1F) * POINTS_BETWEEN_TICKS + 1
    crange = CRANGE
    koffrange = np.logspace(LOG_START_KOFF_1F, LOG_END_KOFF_1F, TOTAL_POINTS_KOFF_1F)

    def log_posterior_x(c, koff, n):
        if c == 0:
            c = 1E-16
        if koff == 0:
            koff = 1E-16
        x = KON * c / koff
        mu = KP * T * x / (1 + x)
        var = (KP * T * x)/(1 + x) + (2 * KP**2 * T * x)/(koff * (1 + x)**3) * (1 + (np.exp(-T * koff * (1 + x)) - 1)/(T * koff * (1 + x)))
        post = -(n - mu)**2 / (2 * var) - 0.5 * np.log(2 * 3.14159 * var)
        return post

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = log_posterior_x(cval, koffval, nobs)

    arr = -1.0 * arr  # hack to get negative probability with log axes, linear colorbar
    label = r'ln $P(c, k_{\mathrm{off}}|n)$'
    linear_contours = []

    vmin_setting = 0
    vmax_setting = 100
    if np.all(arr) > 0.0:
        fig, ax, cb = __heatmap2__(arr, crange * KON / KP, koffrange / KP, figname, '', save=False, log_norm_axes_only=True,
                                   levels=linear_contours, vmin=vmin_setting, vmax=vmax_setting, contour_color='k', contour_linewidths=0.5,
                                   cmap_colour='viridis_r', fig_width=2.8, return_cb=True)
    else:
        print("Not all posterior points evaluated to non-zero probability")
        return 1

    # relabel colorbar with the negative sign
    if vmin_setting == 0:
        tick_space = (vmax_setting - vmin_setting)/(len(cb.ax.get_yticklabels())-1)
        new_cbar_ticks = [u'$0$'] + [u'-${}$'.format(str(i)) for i in np.linspace(tick_space, vmax_setting, len(cb.ax.get_yticklabels())-1, dtype=int)]
    else:
        new_cbar_ticks = [u'-' + item.get_text() for item in cb.ax.get_yticklabels()]

    cb.ax.set_yticklabels(new_cbar_ticks, fontsize=FS)
    cb.ax.invert_yaxis()  # flips the colorbar so that most negative number is on bottom

    ax.set_xlabel(r'$k_{\mathrm{on}}c/k_{p}$', fontsize=FS)
    ax.set_ylabel(r'$k_{\mathrm{off}}/k_{p}$', fontsize=FS)
    ax.set_title(label, fontsize=FS)

    # Superimpose heuristic estimate
    def heuristic_estimate(c, n):
        koff_est = ((KP * T - n) * KON * c) / n
        return koff_est
    koff_star = [heuristic_estimate(ctilde, nobs) for ctilde in crange]

    # rescale onto weird axes scale which do not run from min to max of koffrange and crange anymore
    max_log10koff = np.log10(np.max(koffrange))
    min_log10koff = np.log10(np.min(koffrange))
    max_koff_pixel_pos = len(koffrange)
    slope = max_koff_pixel_pos / (max_log10koff - min_log10koff)
    intercept = min_log10koff * max_koff_pixel_pos / (min_log10koff - max_log10koff)
    koff_star_rescaled = []
    c_pts = []
    for idx, kstar in enumerate(np.log10(koff_star)):
        if kstar > min_log10koff and kstar < max_log10koff:
            pixel = kstar * slope + intercept
            if pixel < max_koff_pixel_pos and pixel > 0:
                koff_star_rescaled.append(pixel)
                c_pts.append(idx)
    ax.plot(c_pts, koff_star_rescaled, 'k--')

    # save figure
    fig.savefig(DIR_OUTPUT + os.sep + 'Figure_1F' + '.pdf', transparent=True)
    fig.savefig(DIR_OUTPUT + os.sep + 'Figure_1F' + '.eps')


def plot_1E_and_2B():
    figname1 = 'Figure_1E'
    figname2 = 'Figure_2B'
    font_size = FS

    mpl.rcParams['xtick.labelsize'] = font_size
    mpl.rcParams['ytick.labelsize'] = font_size

    # Fig 1E
    koffrange = np.logspace(LOG_START_KOFF-1, LOG_END_KOFF, TOTAL_POINTS_KOFF)
    c = KP/KON
    vecErrorX1 = eqns.dedimRelErrorX1NoTrace(c, koffrange)

    # Fig 2B
    crange = CRANGE*10
    koff = 1.0
    vecRelErrorEst2C = eqns.dedimRelErrC2NoTrace(crange, koff)
    vecRelErrorEst2K = eqns.dedimRelErrK2NoTrace(crange, koff)

    # plot
    fig_width_1E = 2.8
    plt.figure(figsize=(fig_width_1E, fig_width_1E*0.8))
    ax = plt.gca()

    plt.plot(koffrange/KP, vecErrorX1/N, color='purple')
    ax.set_xlabel(r'$k_{\mathrm{off}}/k_{p}$', fontsize=font_size)
    ax.set_xscale('log')
    ax.set_ylabel(r'$\frac{k_{p}t}{N} \frac{\langle\delta x^{2}\rangle} {x^{2}}$', fontsize=font_size)
    ax.set_xlim([1E-1, 1E2])
    plt.ylim([0, 1])
    plt.savefig(DIR_OUTPUT + os.sep + figname1 + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname1 + '.eps')

    plt.close()

    # plot
    fig_width_2B = 2.2
    plt.figure(figsize=(fig_width_2B, fig_width_2B*0.8))
    ax = plt.gca()

    plt.plot(crange*KON/KP, vecRelErrorEst2C/N, color='purple', label=r'$c$')
    plt.plot(crange*KON/KP, vecRelErrorEst2K/N, color='orangered', label=r'$k_{\mathrm{off}}$')
    plt.legend()
    ax.set_xlabel(r'$k_{on}c/k_{p}$', fontsize=font_size)
    ax.set_xscale('log')
    ax.set_ylabel(r'$\frac{k_{p}t}{N} \frac{\langle\delta (\cdot)^{2}\rangle}{(\cdot)^{2}}$', fontsize=font_size)
    ax.set_xlim([5E-4, 1.5E2])
    plt.ylim([0, 50])
    plt.savefig(DIR_OUTPUT + os.sep + figname2 + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname2 + '.eps')

    plt.close()


def plot_S1():
    """
    'mode1_MLE_compare' figure, 2 curves (numeric vs heuristic)
    """
    figname = 'mode1_MLE_compare'
    curve1 = DATADICT[figname + '_numeric']
    curve2 = DATADICT[figname + '_heuristic']
    # plot
    plt.figure(figsize=(6, 4))
    c2_part1 = plt.plot(curve2['xpts'][0:400], curve2['ypts'][0:400], color=COLOR_SCHEME['heuristic'], label='heuristic')
    c1 = plt.plot(curve1['xpts'], curve1['ypts'], marker='o', linestyle='None', color=COLOR_SCHEME['numerical_fisher_sp'],
                  label='numeric')
    # plt.title('Mode 1 MLE: Numeric vs Heuristic')# ($k_p=10$, $t=1000$, $k_{off}=1$)')
    plt.xlabel(r'$n$')
    plt.ylabel(r'$x^{*}$')
    plt.legend()
    # save figure
    plt.gca().set_ylim([-5, max(curve1['ypts'])])
    plt.savefig(DIR_OUTPUT + os.sep + 'Figure_S1' + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + 'Figure_S1' + '.eps')


def plot_S2():
    """
    'combined_MLE_composite_compare' figure, 4 curves (numeric vs heuristic) for c and koff
    """
    figname = 'combined_MLE_composite_compare'
    curveL1 = DATADICT[figname + '_c_numeric']
    curveL2 = DATADICT[figname + '_c_heuristic']
    curveR1 = DATADICT[figname + '_koff_numeric']
    curveR2 = DATADICT[figname + '_koff_heuristic']

    # plot
    fig = plt.figure()
    gs1 = gridspec.GridSpec(1, 2)
    axarr = [fig.add_subplot(ss) for ss in gs1]
    fig.set_size_inches(8, 4)

    # left subplot
    cL2 = axarr[0].plot(curveL2['xpts'], curveL2['ypts'], color=COLOR_SCHEME['heuristic'], label='heuristic')
    cL1 = axarr[0].plot(curveL1['xpts'], curveL1['ypts'], marker='o', linestyle='None', color=COLOR_SCHEME['numerical_fisher_sp'], label='numeric')

    axarr[0].set_xlabel(r'$n$')
    axarr[0].set_ylabel(r'$c^*$')
    axarr[0].legend()

    # right subplot
    cR2 = axarr[1].plot(curveR2['xpts'], curveR2['ypts'], color=COLOR_SCHEME['heuristic'], label='heuristic')
    cR1 = axarr[1].plot(curveR1['xpts'][1:], curveR1['ypts'][1:], marker='o', linestyle='None', color=COLOR_SCHEME['numerical_fisher_sp'], label='numeric')

    axarr[1].set_xlabel(r'$n$')
    axarr[1].set_ylabel(r'$k_{off}^*$')
    axarr[1].legend()

    # save figure
    axarr[1].set_ylim([0, max(curveR1['ypts'][1:]) * 1.1])
    gs1.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
    plt.savefig(DIR_OUTPUT + os.sep + 'Figure_S2' + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + 'Figure_S2' + '.eps')


def plot_S3():
    dictionary_SI = {'subdir2': 'SI_ratios', 'log': True, 'plots': {
                            'ratioSigmaC2': {'num': eqns.SigmacrlbC2NoTrace, 'denom': eqns.SigmacrlbC2,
                                             'label': r'$\langle\delta c^{2}\rangle^A$/$\langle\delta c^{2}\rangle^B$'},
                            'ratioSigmaK2': {'num': eqns.SigmacrlbK2NoTrace, 'denom': eqns.SigmacrlbK2,
                                             'label': r'$\langle\delta k_{off}^{2}\rangle/^A\langle\delta {k_{off}^{2}}\rangle^B$'},
                            'ratioDetSigma2': {'num': eqns.DetSigmacrlb2NoTrace, 'denom': eqns.DetSigmacrlb2,
                                               'label': r'det$ (\Sigma^A)$/det$ (\Sigma^B)$'},
                        }
                     }
    subdir_2_use = 'ligand1'
    if not os.path.exists(DIR_OUTPUT + os.sep + subdir_2_use + os.sep + dictionary_SI['subdir2']):
        os.makedirs(DIR_OUTPUT + os.sep + subdir_2_use + os.sep + dictionary_SI['subdir2'])

    contour_args_SI = {'cmap_colour': 'PuBu', 'levels': [1.01, 1.1, 10.0, 100.0, 1000.0], 'fmt': '%.2f'}

    # Plot for choice of KP*T used in settings.py
    id = str(int(KP*T))
    for ratio in dictionary_SI['plots']:
        __heatmap_ratio__(dictionary_SI['plots'][ratio]['num'], dictionary_SI['plots'][ratio]['denom'],
                          dictionary_SI['plots'][ratio]['label'], subdir_2_use + os.sep + dictionary_SI['subdir2'] + os.sep + ratio + id,
                          dictionary_SI['log'], crange=np.logspace(LOG_START_C, LOG_END_C, TOTAL_POINTS_C),
                          koffrange=np.logspace(-2, 4.2, TOTAL_POINTS_C), dedim=True, contour_args=contour_args_SI)


def plot_S5():
    """
    'mode1_MLE_with_prior_compare' figure, 2 curves (numeric vs heuristic)
    """
    figname = 'mode1_MLE_with_Prior_compare'
    curve1 = DATADICT['mode1_MLE_with_Prior_compare_numeric']
    curve2 = DATADICT['mode1_MLE_compare_heuristic']
    # plot
    plt.figure(figsize=(6, 4))
    c2_part1 = plt.plot(curve2['xpts'][0:400], curve2['ypts'][0:400], color=COLOR_SCHEME['heuristic'], label='heuristic')
    c1 = plt.plot(curve1['xpts'], curve1['ypts'], marker='o', linestyle='None', color=COLOR_SCHEME['numerical_fisher_sp'],
                  label='numeric')
    # plt.title('Mode 1 MLE: Numeric vs Heuristic')# ($k_p=10$, $t=1000$, $k_{off}=1$)')
    plt.xlabel(r'$n$')
    plt.ylabel(r'$x^{*}$')
    plt.legend()
    # save figure
    plt.gca().set_ylim([-5, max(curve1['ypts'])])
    plt.savefig(DIR_OUTPUT + os.sep + 'Figure_S5' + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + 'Figure_S5' + '.eps')


if __name__ == '__main__':
    # ----------------------------------------------
    # BEGIN: USER CONTROLS
    # ----------------------------------------------
    # line plot control variables
    flag_Fig1F = False
    flag_Fig1E_and_2B = False
    # heatmap control variables
    flag_Fig1_and_Fig2 = False
    flag_Fig3 = False
    ADD_TRACE_TERM = False
    LOG_SELECT = True
    # supplementary figures control variables
    flag_S1 = False
    flag_S2 = False
    flag_S3 = False
    flag_S5 = True

    # choose any of these 2 to be arrays [0: c1, 1: koff, 2: c2, 3: koff2], they will be your axes
    dim = {'x': 2, 'y': 3}

    # axis we'd like to plot
    value = [C1, KOFF, C2, KOFF2]
    label_fig = ['c1', 'koff1', 'c2', 'koff2']
    label = [r'$c_1$', r'$k_{\mathrm{off},1}$', r'$c_2$', r'$k_{\mathrm{off},2}$']
    dimension = [CRANGE, KOFFRANGE, CRANGE, KOFFRANGE2]
    axes = label_fig[dim['x']] + '_' + label_fig[dim['y']]

    # ----------------------------------------------
    # END: USER CONTROLS
    # ----------------------------------------------
    if flag_Fig3 or flag_Fig1_and_Fig2:
        for dirs in ['ligand1', 'ligands2', 'ligands2_over_ligands1']:
            if not os.path.exists(DIR_OUTPUT + os.sep + dirs):
                os.makedirs(DIR_OUTPUT + os.sep + dirs)

        # for plotting with dedimensionalized variables
        dedimension = [dimension[0] * KON / KP, dimension[1] / KP, dimension[2] * KON / KP, dimension[3] / KP]
        dedimension_label = [r'$k_{\mathrm{on}}c_1/k_{p}$', r'$k_{\mathrm{off},1}/k_{p}$', r'$k_{\mathrm{on}}c_2/k_{p}$', r'$k_{\mathrm{off},2}/k_{p}$']

        # Figure 1 and 2 heatmaps
        arrRelErrorEst1X = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))
        arrRelErrorEst2C = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))
        arrRelErrorEst2KOFF = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))

        # heatmap for any element of sigmaEst(4x4 matrix)
        arrSigmaEst = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4))
        arrRelErrorEst = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4))
        arrSigmaData = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4))
        arrDmudthetaInv = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']]), 4, 4))

        # heatmap for the det estimate covariance
        arrDetSigmaEst = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))
        arrDetSigmaData = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))
        arrDetDmudthetaInv = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))
        arrRelDetSigmaEst = np.zeros((len(dimension[dim['y']]), len(dimension[dim['x']])))

        # making heatmaps
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
                arrSigmaEst[j, i, :, :], arrDetSigmaEst[j, i], arrRelErrorEst[j, i, :, :], arrRelDetSigmaEst[j, i] = \
                    __sigmaEst__(value[0], value[1], value[2], value[3], add_trace_term=ADD_TRACE_TERM)
                arrRelErrorEst1X[j, i] = eqns.dedimRelErrorX1NoTrace(value[2], value[3])
                arrRelErrorEst2C[j, i] = eqns.dedimRelErrC2NoTrace(value[2], value[3])
                arrRelErrorEst2KOFF[j, i] = eqns.dedimRelErrK2NoTrace(value[2], value[3])

        print("Done computing matrices\n Plotting results")
        ver = "new_"
        multi_fname = "multi_" + ver + axes

        if flag_Fig1_and_Fig2:
            figure_1_and_2_heatmaps(arrRelErrorEst1X, arrRelErrorEst2C, arrRelErrorEst2KOFF, dedimension[dim['x']],
                                    dedimension[dim['y']], 'Fig1D', 'Fig2C', 'Fig2D',
                                    [r'$k_{\mathrm{on}}c/k_{p}$', r'$k_{\mathrm{off}}/k_{p}$'], LOG_SELECT)

        if flag_Fig3:
            figure_3_heatmaps(arrRelDetSigmaEst, arrRelErrorEst, dedimension[dim['x']], dedimension[dim['y']],
                              multi_fname,
                              [dedimension_label[dim['x']], dedimension_label[dim['y']]], LOG_SELECT)

    if flag_Fig1E_and_2B:
        plot_1E_and_2B()

    if flag_Fig1F:
        plot_1F()

    if flag_S1 or flag_S2 or flag_S5:
        from load_inputs import gen_datadict
        DATADICT = gen_datadict(verbose=False)

    if flag_S1:
        plot_S1()

    if flag_S2:
        plot_S2()

    if flag_S3:
        plot_S3()

    if flag_S5:
        plot_S5()
