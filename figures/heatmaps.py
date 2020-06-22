import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

import plotting_dictionaries as pd
import matplotlib.ticker as ticker
import equations as eqns
from decimal import Decimal


#from load_inputs import DATADICT
from settings import DIR_OUTPUT, DIR_INPUT
from settings import COLOR_SCHEME as cs
from settings import KON, KP, T, KF, ALPHA, N


plt.style.use('parameters.mplstyle')  # particularIMporting

# de-dimensionalise params
c0 = KP/KON
alpha = KP*T

# plot params
FS = 10
SHOW = False

# axes
POINTS_BETWEEN_TICKS = 10
LOG_START_C = -9
LOG_END_C = -4
TOTAL_POINTS_C = (LOG_END_C - LOG_START_C) * POINTS_BETWEEN_TICKS + 1
CRANGE = np.logspace(LOG_START_C, LOG_END_C, TOTAL_POINTS_C)
LOG_START_KOFF = -2
LOG_END_KOFF = 4
TOTAL_POINTS_KOFF = (LOG_END_KOFF - LOG_START_KOFF) * POINTS_BETWEEN_TICKS + 1
KOFFRANGE = np.logspace(LOG_START_KOFF, LOG_END_KOFF, TOTAL_POINTS_KOFF)

CTILDERANGE = np.divide(list(CRANGE), c0)
ZRANGE = np.divide(list(KOFFRANGE), KP)

class MidPointLogNorm(mpl.colors.LogNorm):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        mpl.colors.LogNorm.__init__(self,vmin=vmin, vmax=vmax, clip=clip)
        self.midpoint=midpoint
    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [np.log(self.vmin), np.log(self.midpoint), np.log(self.vmax)], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(np.log(value), x, y))

def plot_heatmap(arr, crange, koffrange, fname, label, show=SHOW, save=True, log_norm=True, dedim=False,
                 xy_label_force=None, less_xticks=True, **kwargs):
    """
    crange: range of values for c
    koffrange: range of koff values
    fname: file saved as pdf and eps format with this name
    show: shows plots
    save: saves plots
    log_norm: heatmap is log, this is only used when plotting the posterior
    dedim: makes the axis dedimensionalized. scales them (see by how much below)
    kwargs: a variety of keyword args are used below. They are mostly used for contour plot lines if I understand correctly.
    Using Duncan's default ones mostly, but we can talk about it.
    """
    # default parameters
    if 'levels' in kwargs.keys(): levels = kwargs['levels']
    else: levels = [1E0, 1E1, 1E2, 1E3]

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

    if 'fig_width' in kwargs.keys(): fig_width = kwargs['fig_width']
    else: fig_width = 3.0

    if dedim:
        # if flag is true, this is how we scale the axis. Simple.
        crange = crange*KON/KP;
        koffrange = koffrange/KP;
        xy_label = [r'$k_{\mathrm{on}}c/k_{p}$', r'$k_{\mathrm{off}}/k_{p}$'];
    else:
        xy_label = [r'${c}$', r'${k}_{\mathrm{off}}$']
    if xy_label_force is not None:
        xy_label = xy_label_force

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
    f = plt.figure(figsize=(fig_width,fig_width/1.2))
    im = plt.imshow(arr, interpolation='spline36', **imshow_kw)

    # axes setup
    fig = plt.gcf(); ax = plt.gca()

    # method 1
    # axes log scaled
    # axes log scaled
    if less_xticks:
        ax.set_xticks([i for i, cval in enumerate(crange) if (i % POINTS_BETWEEN_TICKS==0 and np.log10(cval)%2 == 0)])
        ax.set_xticklabels([r'$10^{%d}$' % np.log10(cval) for i, cval in enumerate(crange) if (i % POINTS_BETWEEN_TICKS==0 and np.log10(cval)%2 == 0)], fontsize=FS)
    else:
        ax.set_xticks([i for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS == 0])
        ax.set_xticklabels([r'$10^{%d}$' % np.log10(cval) for i, cval in enumerate(crange) if (i % POINTS_BETWEEN_TICKS==0)], fontsize=FS)
    ax.set_yticks([i for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0])
    ax.set_yticklabels([r'$10^{%d}$' % np.log10(yval) for i, yval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
    """
    if log_norm:
        ax.set_xticks([i for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS == 0])
        ax.set_yticks([i for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0])
        ax.set_xticklabels([r'$10^{%d}$' % np.log10(cval) for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
        ax.set_yticklabels([r'$10^{%d}$' % np.log10(kval) for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
    else:
        # Set ticks later because for some reason plotting contours messes up the ticks that get set here
        pass
    """

    ax.invert_yaxis()
    ax.set_xlabel(xy_label[0], fontsize=FS); ax.set_ylabel(xy_label[1], fontsize=FS)

    # create colorbar
    cbar = fig.colorbar(im,fraction=0.0375, pad=0.04)
    #cbar.locator = ticker.LogLocator(base=10)
    cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=FS, labelpad=5); cbar.ax.tick_params(labelsize=FS)
    cbar.ax.minorticks_off();
    # UNCOMMENT THIS ONLY WHEN TICKS DON'T APPEAR
    #cbar.set_ticks([round(vmin,3)+0.001,round(vmax,3)-0.001])
    cbar.update_ticks()
    cbar.ax.minorticks_off();

    # TODO IDK why do ticks hide sometimes?
    CL = plt.contour(arr, levels=levels, linestyles=contour_linestyle, colors=contour_color, linewidths=contour_lindewidths)
    plt.clabel(CL, CL.levels, inline=True, fmt=fmt)


    """
    print(vmin, vmax)
    cbar.set_ticks([1.0, 10.0, 100.0])
    cbar.set_ticklabels([1.0, 10.0, 100.0])
    """

    # TODO IDK why do ticks hide sometimes?
    #for t in cbar.ax.get_yticklabels(): print(t.get_text())
    #plt.contour(arr, levels=levels, linestyles=contour_linestyle, colors=contour_color, linewidths=contour_lindewidths)
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
    if save == True:
        plt.savefig(DIR_OUTPUT + os.sep + fname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + fname + '.eps'); plt.savefig(DIR_OUTPUT + os.sep + fname + '.png')
    if show:
        plt.show()

    plt.close()
    return fig, ax


def heatmap_mode1_error_x(crange=CTILDERANGE, koffrange=ZRANGE, make_heatmap=True, make_panel=False,
                          scale_factor=alpha, divide_by_N=True, label_style=0):
    if make_heatmap == True:
        def mode1_error_x(ctilde, z, scale_factor=scale_factor):
            c = ctilde * c0
            koff = z * KP
            x = c * KON / koff
            val = (1 + x) / (KP * T * x) * ((1 + x) ** 2 + 2 * KP / koff)
            if divide_by_N:
                return scale_factor * val / N
            else:
                return scale_factor * val

        arr = np.zeros((len(koffrange), len(crange)))
        for i, koffval in enumerate(koffrange):
            for j, cval in enumerate(crange):
                arr[i, j] = mode1_error_x(cval, koffval)

        label = r'$\alpha \langle\delta x^{2}\rangle$/$x^{2}$'
        plot_heatmap(arr, crange, koffrange, 'heatmap_mode1_heuristic_error_x', label)

    if make_panel == True:
        def cross_section_mode1_error_c(crange=CTILDERANGE):
            def mode1_error_c(ctilde, z, scale_factor=scale_factor):
                c = ctilde * c0
                koff = z * KP
                x = c * KON / koff
                val = (1 + x) / (KP * T * x) * ((1 + x) ** 2 + 2 * KP / koff)
                if divide_by_N:
                    return scale_factor * val / N
                else:
                    return scale_factor * val

            arr = [mode1_error_c(cval, 1) for cval in crange]
            return dict({'xpts': crange, 'ypts':arr})

        figname = 'mode1_error_c_cross_section'
        curve1 = cross_section_mode1_error_c()
        # plot
        plt.figure(figsize=(3.25*0.8, 3*0.8))
        #plt.figure()
        plt.plot(curve1['xpts'], curve1['ypts'], color=cs['simple_fisher'], label='Simple Fisher', zorder=1)
        # axis
        if label_style == 0:
            #plt.title('Mode 1: MLE relative error comparison \n'+r'($\tilde{c}_0=10$, $\alpha=1 \times 10^4$, $k_{p}=10$)')
            if divide_by_N:
                plt.ylabel(r'$\frac{k_{p} t}{N} \langle\delta c^{2}\rangle$/$c^{2}$')
            else:
                plt.ylabel(r'$k_{p} t \langle\delta c^{2}\rangle$/$c^{2}$')
        elif label_style == 1:
            plt.title('Mode 1: MLE relative error comparison \n($k_p=10$, $t=100$, $k_{\mathrm{off}}=k_{\mathrm{on}}=1$)')
            plt.ylabel(r'$\langle\delta c^{2}\rangle$/$c^{2}$')

        plt.xlabel(r'$k_{\mathrm{on}}c/k_{p}$')

        plt.gca().set_xscale('log')
        plt.xlim([0.03, 10])
        plt.ylim([0, 100])
        plt.xscale('log')
        #plt.legend()
        # save figure
        plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
        plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    return


def combined_error_c(ctilde, z, scale_factor=alpha):
    c = ctilde * c0
    koff = z * KP
    x = KON * c / koff
    num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff**2 * T * x ** 2 * (1 + x) ** 3
    den = koff ** 2 * KP * T**2 * x * (1 + x) ** 2
    val = num / den
    return scale_factor * val


def heatmap_combined_error_c(crange=CTILDERANGE, koffrange=ZRANGE, scale_factor=alpha, label_style=0):
    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = combined_error_c(cval, koffval, scale_factor=scale_factor)

    if label_style == 0:
        label = r'$k_{p}t \langle\delta c^{2}\rangle$/$c^{2}$'
    elif label_style == 1:
        label = r'$\langle\delta c^{2}\rangle$/$c^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_combined_heuristic_error_c', label)
    return


def combined_error_koff(ctilde, z, scale_factor=alpha):
    c = ctilde * c0
    koff = z * KP
    x = c * KON / koff
    num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff**2 * T * (1 + x) ** 3
    den = koff ** 2 * KP * T**2 * x * (1 + x) ** 2
    val = num / den
    return scale_factor * val


def heatmap_combined_error_koff(crange=CRANGE, koffrange=KOFFRANGE, scale_factor=alpha, label_style=0):
    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = combined_error_koff(cval, koffval, scale_factor=scale_factor)

    if label_style == 0:
        label = r'$k_{p}t \langle\delta k_{\mathrm{off}}^{2}\rangle$/$k_{\mathrm{off}}^{2}$'
    elif label_style == 1:
        label = r'$\langle\delta k_{\mathrm{off}}^{2}\rangle$/$k_{\mathrm{off}}^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_combined_heuristic_error_koff', label)
    return


def figure_2_combined_cross_sections(crange=CRANGE, koffrange=KOFFRANGE,
                             scale_factor=alpha, label_style=0):
    """
    Produces cross sections of the c and koff relative error heatmaps.
    Both cross sections are for constant koff = 1 while varying c (ie. horizontal cross sections)
    """
    def combined_error_c(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        x = c * KON / koff
        num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff ** 2 * T * x ** 2 * (1 + x) ** 3
        den = koff ** 2 * KP * T ** 2 * x * (1 + x) ** 2
        val = num / den
        return scale_factor * val
    def combined_error_koff(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        x = c * KON / koff
        num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff**2 * T * (1 + x) ** 3
        den = koff ** 2 * KP * T**2 * x * (1 + x) ** 2
        val = num / den
        return scale_factor * val

    def cross_section_combined_error_c():
        arr = [combined_error_c(cval, 1) for cval in crange]
        return dict({'xpts': crange, 'ypts':arr})
    def cross_section_combined_error_koff():
        arr = [combined_error_koff(c, 1) for c in crange]
        return dict({'xpts': koffrange, 'ypts':arr})

    figname = 'combined_error_cross_sections'
    curve1 = cross_section_combined_error_c()
    curve2 = cross_section_combined_error_koff()
    # plot
    plt.figure(figsize=(3, 3))
    ax1 = plt.gca()
    ax2 = ax1.twiny()

    if label_style == 0:
        ln1 = ax1.plot(curve1['xpts'], curve1['ypts'], color=cs['simple_fisher'], label=r'$c$', zorder=1)
        ln2 = ax2.plot(curve2['xpts'], curve2['ypts'], color=cs['heuristic'], label=r'$k_{\mathrm{off}}$', zorder=1)
        #plt.title('Mode 2: MLE relative error comparison\n' + r'($\tilde{c}_0=10$, $\alpha=1 \times 10^4$, $k_{p}=10$)')

    elif label_style == 1:
        ln1 = ax1.plot(curve1['xpts'], curve1['ypts'], color=cs['simple_fisher'], label=r'$\delta c^{2}/c^{2}$', zorder=1)
        ln2 = ax2.plot(curve2['xpts'], curve2['ypts'], color=cs['heuristic'],label=r'$\delta k_{\mathrm{off}}^{2}/k_{\mathrm{off}}^{2}$', zorder=1)
        plt.title('Mode 2: MLE relative error comparison\n' + r'($k_p=10$, $t=100$, $k_{\mathrm{off}}=k_{\mathrm{on}}=1$)')
        plt.ylabel(r'$\langle\delta (\cdot)^{2}\rangle$/$(\cdot)^{2}$')

    # axis
    ax1.set_xlabel(r'$k_{\mathrm{on}}c/k_{p}$')
    ax1.set_ylabel(r'$k_{p}t \langle\delta (\cdot)^{2}\rangle$/$(\cdot)^{2}$')
    #ax2.set_xlabel(r'$k_{\mathrm{off}}$')

    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax1.set_xlim([1E-2, 1E2])
    ax2.set_xlim([1E-2, 1E2])
    plt.ylim([0, 0.01*alpha])

    lns = ln1 + ln2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs)

    plt.tight_layout()
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

def heatmap_kpr_error_c(crange=CTILDERANGE, koffrange=ZRANGE, scale_factor=alpha, label_style=0):

    def kpr_error_c_full(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        x = c * KON / koff
        g = koff / KF
        factor_A = (1 + x) / ((1 + g) * koff**2 * KP * T**2 * x)
        term1 = 2 * g * (KP + koff * KP * T - koff**2 * T * x)
        term2 = koff * T * (KP + koff * x**2)
        term3 = g**2 * (koff**2 * T + KP * (2 * (1 + x) + koff * T * (2 + x * (2 + x))))
        val = factor_A * (term1 + term2 + term3)
        return scale_factor * val

    def kpr_error_c_high_g_t(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        x = c * KON / koff
        left = (1 + koff / KF) * (1 + x) / (koff * T * x)
        right = koff / KP + 2 + 2 * x + x ** 2
        val = left * right
        return scale_factor * val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = kpr_error_c_full(cval, koffval)
    if label_style == 0:
        label = r'$\alpha \langle\delta c^{2}\rangle$/$c^{2}$'
    elif label_style == 1:
        label = r'$\langle\delta c^{2}\rangle$/$c^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_kpr_heuristic_error_c', label)
    return

def heatmap_kpr_error_koff(crange=CTILDERANGE, koffrange=ZRANGE,
                        scale_factor=alpha, label_style=0):

    def kpr_error_koff(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        # Note the "full" heuristic expression for koff error is the same as the high g high t one
        x = c * KON / koff
        left = (1 + koff / KF) * (1 + x) / (koff * T * x)
        right = koff / KP + 1
        val = left * right
        return scale_factor * val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = kpr_error_koff(cval, koffval)
    if label_style == 0:
        label = r'$\alpha \langle\delta k_{\mathrm{off}}^{2}\rangle$/$k_{\mathrm{off}}^{2}$'
    elif label_style == 1:
        label = r'$\langle\delta k_{\mathrm{off}}^{2}\rangle$/$k_{\mathrm{off}}^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_kpr_heuristic_error_koff', label)
    return

def heatmap_kpr2_error_c(crange=CTILDERANGE, koffrange=ZRANGE,
                        scale_factor=alpha, label_style=0):

    def kpr2_error_c_full(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        kon = KON
        kf = KF
        kp = KP
        t = T
        g = koff / KF
        x = c * kon / koff

        val = ((kf + koff + kf*x + koff*x)**4 * (
                                                (koff**2 * kp**2 * t**3 * x**2 * (1 + x)*((kf**2 * kp**2 * t**2 * x**2)/(kf + koff + kf*x + koff*x)**2 +
                                                                                          (2*kf*koff*kp**2 * t**2 * x**2)/(kf + koff + kf*x + koff*x)**2 +
                                                                                          (koff**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 -
                                                                                          (koff*kp**2*t**2*x)/(kf + koff + kf*x + koff*x))**2 *
                                                                                        ((1 + g)**2*koff*(1 + x)**2 + 2*kp*(1 + g**2*(1 + x)**2 + g*(2 + x))))/
                                                ((1 + g)**3*(kf + koff + kf*x + koff*x)**2) +
                                                (g*kf**2 * kp**2 * t**3 * x**2 * (1 + x)*((1 + g)**2*koff*(1 + x)**2 + 2*g*kp*(1 + g + x + x**2)) *
                                                                                                                       ((koff**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 +
                                                                                                                        (kf*kp*t*x*(-(kp*t) + (kf*kp*t*x)/(kf + koff + kf*x + koff*x)))/(kf + koff + kf*x + koff*x) +
                                                                                                                        (2*koff*kp*t*x*(-(kp*t) + (kf*kp*t*x)/(kf + koff + kf*x + koff*x)))/(kf + koff + kf*x + koff*x))**2) /
                                                ((1 + g)**3*(kf + koff + kf*x + koff*x)**2) +
                                                (2*kf**2 * koff * kp**3 * t**2 * x**2*(-((kf**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2) -
                                                                                       (2*kf*koff*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 -
                                                                                       (koff**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 +
                                                                                       (koff*kp**2*t**2*x)/(kf + koff + kf*x + koff*x))*(-((kf**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2) -
                                                                                                                                         (2*kf*koff*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 -
                                                                                                                                         (koff**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 +
                                                                                                                                         (kf*kp**2*t**2*x)/(kf + koff + kf*x + koff*x) +
                                                                                                                                         (2*koff*kp**2*t**2*x)/(kf + koff + kf*x + koff*x))*(kf**2*(-1 + x)*(-1 + koff*t*(1 + x)) +
                                                                                                                                                                                             kf*koff*(3 - 2*x - x**2 + 2*koff*t*(-1 + x + 3*x**2 + x**3)) +
                                                                                                                                                                                             koff**2*(2 - 3*x - 5*x**2 - 2*x**3 + koff*t*(-1 + 2*x + 5*x**2 + 2*x**3)))) /
                                                ((kf + koff)**4*(kf + koff + kf*x + koff*x)**2))) / \
              (c**2*kf**2*koff*kon**2*kp**3*t**4*x**3*(1 + x)**4*(-(kp*t) + (kf*kp*t*x)/(kf + koff + kf*x + koff*x) + (koff*kp*t*x)/(kf + koff + kf*x + koff*x))**4)
        return scale_factor * val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = kpr2_error_c_full(cval, koffval)
    if label_style == 0:
        label = r'$\alpha \langle\delta c^{2}\rangle$/$c^{2}$'
    elif label_style == 1:
        label = r'$\langle\delta c^{2}\rangle$/$c^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_kpr2_heuristic_error_c', label)
    return

def heatmap_kpr2_error_koff(crange=CTILDERANGE, koffrange=ZRANGE,
                        scale_factor=alpha, label_style=0):

    def kpr2_error_koff(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        # Note the "full" heuristic expression for koff error is the same as the high g high t one
        kon = KON
        kf = KF
        kp = KP
        t = T
        g = koff / KF
        x = c * kon / koff

        val = ((kf + koff + kf*x + koff*x)**4 * ((-2*kf**2 * koff * kp**3 * t**2 * x**2*((kf + koff)*(-kf - 2*koff + koff*(kf + koff)*t) +
                                                                                         (kf**2 + 2*kf*koff + 3*koff**2 - 2*koff**2*(kf + koff)*t)*x -
                                                                                         koff*(kf + 5*koff)*(-1 + (kf + koff)*t)*x**2 -
                                                                                         2*koff**2*(-1 + (kf + koff)*t)*x**3)) /
                                                 ((kf + koff)**4*(kf + koff + kf*x + koff*x)**2) +
                                                 (g*kf**2 * kp**2 * t**3 * x**2*(1 + x)*((1 + g)**2*koff*(1 + x)**2 + 2*g*kp*(1 + g + x + x**2))) /
                                                 ((1 + g)**3*(kf + koff + kf*x + koff*x)**2) +
                                                 (koff**2*kp**2*t**3*x**2*(1 + x)*((1 + g)**2*koff*(1 + x)**2 +
                                                                                   2*kp*(1 + g*(2 + x + g*(1 + x)**2)))) /
                                                 ((1 + g)**3*(kf + koff + kf*x + koff*x)**2))) / \
              (kf**2*koff**3*kp**3*t**4*x**3*(1 + x)**4)

        return scale_factor * val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = kpr2_error_koff(cval, koffval)
    if label_style == 0:
        label = r'$\alpha \langle\delta k_{\mathrm{off}}^{2}\rangle$/$k_{\mathrm{off}}^{2}$'
    elif label_style == 1:
        label = r'$\langle\delta k_{\mathrm{off}}^{2}\rangle$/$k_{\mathrm{off}}^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_kpr2_heuristic_error_koff', label)
    return

def heatmap_ratios(crange=CTILDERANGE, koffrange=ZRANGE,
                        scale_factor=alpha, label_style=0):

    def kpr2_error_koff(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        # Note the "full" heuristic expression for koff error is the same as the high g high t one
        kon = KON
        kf = KF
        kp = KP
        t = T
        g = koff / KF
        x = c * kon / koff

        val = ((kf + koff + kf*x + koff*x)**4 * ((-2*kf**2 * koff * kp**3 * t**2 * x**2*((kf + koff)*(-kf - 2*koff + koff*(kf + koff)*t) +
                                                                                         (kf**2 + 2*kf*koff + 3*koff**2 - 2*koff**2*(kf + koff)*t)*x -
                                                                                         koff*(kf + 5*koff)*(-1 + (kf + koff)*t)*x**2 -
                                                                                         2*koff**2*(-1 + (kf + koff)*t)*x**3)) /
                                                 ((kf + koff)**4*(kf + koff + kf*x + koff*x)**2) +
                                                 (g*kf**2 * kp**2 * t**3 * x**2*(1 + x)*((1 + g)**2*koff*(1 + x)**2 + 2*g*kp*(1 + g + x + x**2))) /
                                                 ((1 + g)**3*(kf + koff + kf*x + koff*x)**2) +
                                                 (koff**2*kp**2*t**3*x**2*(1 + x)*((1 + g)**2*koff*(1 + x)**2 +
                                                                                   2*kp*(1 + g*(2 + x + g*(1 + x)**2)))) /
                                                 ((1 + g)**3*(kf + koff + kf*x + koff*x)**2))) / \
              (kf**2*koff**3*kp**3*t**4*x**3*(1 + x)**4)

        return scale_factor * val
    def kpr2_error_c_full(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        kon = KON
        kf = KF
        kp = KP
        t = T
        g = koff / KF
        x = c * kon / koff

        val = ((kf + koff + kf*x + koff*x)**4 * (
                                                (koff**2 * kp**2 * t**3 * x**2 * (1 + x)*((kf**2 * kp**2 * t**2 * x**2)/(kf + koff + kf*x + koff*x)**2 +
                                                                                          (2*kf*koff*kp**2 * t**2 * x**2)/(kf + koff + kf*x + koff*x)**2 +
                                                                                          (koff**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 -
                                                                                          (koff*kp**2*t**2*x)/(kf + koff + kf*x + koff*x))**2 *
                                                                                        ((1 + g)**2*koff*(1 + x)**2 + 2*kp*(1 + g**2*(1 + x)**2 + g*(2 + x))))/
                                                ((1 + g)**3*(kf + koff + kf*x + koff*x)**2) +
                                                (g*kf**2 * kp**2 * t**3 * x**2 * (1 + x)*((1 + g)**2*koff*(1 + x)**2 + 2*g*kp*(1 + g + x + x**2)) *
                                                                                                                       ((koff**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 +
                                                                                                                        (kf*kp*t*x*(-(kp*t) + (kf*kp*t*x)/(kf + koff + kf*x + koff*x)))/(kf + koff + kf*x + koff*x) +
                                                                                                                        (2*koff*kp*t*x*(-(kp*t) + (kf*kp*t*x)/(kf + koff + kf*x + koff*x)))/(kf + koff + kf*x + koff*x))**2) /
                                                ((1 + g)**3*(kf + koff + kf*x + koff*x)**2) +
                                                (2*kf**2 * koff * kp**3 * t**2 * x**2*(-((kf**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2) -
                                                                                       (2*kf*koff*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 -
                                                                                       (koff**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 +
                                                                                       (koff*kp**2*t**2*x)/(kf + koff + kf*x + koff*x))*(-((kf**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2) -
                                                                                                                                         (2*kf*koff*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 -
                                                                                                                                         (koff**2*kp**2*t**2*x**2)/(kf + koff + kf*x + koff*x)**2 +
                                                                                                                                         (kf*kp**2*t**2*x)/(kf + koff + kf*x + koff*x) +
                                                                                                                                         (2*koff*kp**2*t**2*x)/(kf + koff + kf*x + koff*x))*(kf**2*(-1 + x)*(-1 + koff*t*(1 + x)) +
                                                                                                                                                                                             kf*koff*(3 - 2*x - x**2 + 2*koff*t*(-1 + x + 3*x**2 + x**3)) +
                                                                                                                                                                                             koff**2*(2 - 3*x - 5*x**2 - 2*x**3 + koff*t*(-1 + 2*x + 5*x**2 + 2*x**3)))) /
                                                ((kf + koff)**4*(kf + koff + kf*x + koff*x)**2))) / \
              (c**2*kf**2*koff*kon**2*kp**3*t**4*x**3*(1 + x)**4*(-(kp*t) + (kf*kp*t*x)/(kf + koff + kf*x + koff*x) + (koff*kp*t*x)/(kf + koff + kf*x + koff*x))**4)
        return scale_factor * val
    def kpr_error_koff(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        # Note the "full" heuristic expression for koff error is the same as the high g high t one
        x = c * KON / koff
        left = (1 + koff / KF) * (1 + x) / (koff * T * x)
        right = koff / KP + 1
        val = left * right
        return scale_factor * val
    def kpr_error_c_full(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        x = c * KON / koff
        g = koff / KF
        factor_A = (1 + x) / ((1 + g) * koff**2 * KP * T**2 * x)
        term1 = 2 * g * (KP + koff * KP * T - koff**2 * T * x)
        term2 = koff * T * (KP + koff * x**2)
        term3 = g**2 * (koff**2 * T + KP * (2 * (1 + x) + koff * T * (2 + x * (2 + x))))
        val = factor_A * (term1 + term2 + term3)
        return scale_factor * val
    def combined_error_c(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        x = c * KON / koff
        num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff ** 2 * T * x ** 2 * (1 + x) ** 3
        den = koff ** 2 * KP * T ** 2 * x * (1 + x) ** 2
        val = num / den
        return scale_factor * val
    def combined_error_koff(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        x = c * KON / koff
        num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff**2 * T * (1 + x) ** 3
        den = koff ** 2 * KP * T**2 * x * (1 + x) ** 2
        val = num / den
        return scale_factor * val
    def mode1_error_x(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        x = c * KON / koff
        val = (1 + x) / (KP * T * x) * ((1 + x) ** 2 + 2 * KP / koff)
        return scale_factor * val

    def compute_heatmap(error_func, ctilde, z, scale_factor=scale_factor):
        arr = np.zeros((len(z), len(ctilde)))
        for i, koffval in enumerate(z):
            for j, cval in enumerate(ctilde):
                arr[i, j] = error_func(cval, koffval)
        return arr

    model1_error = compute_heatmap(mode1_error_x, crange, koffrange)
    model2_divided_by_model1 = np.divide(compute_heatmap(combined_error_c, crange, koffrange), model1_error)
    model3_divided_by_model1 = np.divide(compute_heatmap(kpr_error_c_full, crange, koffrange), model1_error)

    #print model2_divided_by_model1

    plot_heatmap(model2_divided_by_model1, crange, koffrange, 'heatmap_ratio_2_over_1',
                 r'Model 2 $\langle \sigma^{2}_{c^{*}}\rangle$/ Model 1 $\langle \sigma_{c^{*}}\rangle$',
                 levels=[0.1, 0.5, 0.99], save=False, show=True)
    plot_heatmap(model3_divided_by_model1, crange, koffrange, 'heatmap_ratio_2_over_1',
                 r'Model 3 $\langle \sigma^{2}_{c^{*}}\rangle$/ Model 1 $\langle \sigma_{c^{*}}\rangle$',
                 levels=[0.1, 0.5, 0.99], save=False, show=True)

def truncate(f, n):
    '''
    Truncates/pads a float f to n decimal places without rounding
    Source: https://stackoverflow.com/questions/783897/truncating-floats-in-python
    '''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return float('.'.join([i, (d+'0'*n)[:n]]))

def heatmap_figure_4():
    nobs = 0.1 * KP * T # 100
    mobs = 0.15 * KP * T # 150
    def mode1_plot(crange=[truncate(f, 3) for f in list(np.arange(0.0 * KON / KP, 5.0 * KON / KP + 0.005, 0.005))[1:]],
                   koffrange=[truncate(f, 2) for f in list(np.arange(0.0 / KP, 50.0 / KP + 0.05, 0.05))[1:]]):
        figname = 'heatmap_log_posterior_mode1'
        def log_posterior_x(c, koff, n):
            if c==0:
                c=5E-4
            if koff==0:
                koff==5E-4
            x = KON * c / koff
            mu = KP * T * x / (1 + x)
            var = (KP * T * x)/(1 + x) + (2 * KP**2 * T * x)/(koff * (1 + x)**3) *(1 + (np.exp(-T * koff * (1 + x)) - 1)/(T * koff * (1 + x)))
            post = -(n - mu)**2 / (2 * var) - 0.5 * np.log(2 * 3.14159 * var)
            return post

        arr = np.zeros((len(koffrange), len(crange)))
        for i, koffval in enumerate(koffrange):
            for j, cval in enumerate(crange):
                arr[i, j] = log_posterior_x(cval * KP/KON, koffval * KP, nobs)

        label = r'ln $P(c, k_{\mathrm{off}}|n)$'
        #log_contours = [-i for i in list(np.logspace(-3, 4,num=20))[::-1]]
        linear_contours = list(range(-500, 0, 50))

        fig, ax = plot_heatmap(arr, crange, koffrange, figname, label, save=False, log_norm=False,
                               levels=linear_contours, vmin=-500, vmax=0, contour_color='w', contour_linewidths=0.5)

        ax.set_xlabel(r'$k_{\mathrm{on}}c/k_{p}$', fontsize=FS)
        ax.set_ylabel(r'$k_{\mathrm{off}}/k_{p}$', fontsize=FS)

        # Superimpose heuristic estimate
        def heuristic_estimate(c, n):
            koff_est = ((KP * T - n) * KON * c) / n
            return koff_est

        # rescale onto weird axes scale which do not run from min to max of koffrange and crange anymore
        xpts = np.linspace(ax.get_xticks()[0], ax.get_xticks()[-1], len(crange))

        koff_star = [ax.get_ylim()[1] * (heuristic_estimate((ctilde - (crange[1]-crange[0])) * KP/KON, nobs)/KP -
                                         (koffrange[1]-koffrange[0]))/max(koffrange) for ctilde in crange]
        ax.plot([i+ax.get_xlim()[0] for i in xpts], [i-ax.get_ylim()[0] for i in koff_star], 'k--')

        # save figure
        fig.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
        fig.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    # ------------------
    # Repeat for Model 2
    # ------------------
    def mode2_plot(crange=[truncate(f, 3) for f in list(np.arange(0.0 * KON / KP, 5 * KON / KP + 0.005, 0.005))[1:]],
                   koffrange=[truncate(f, 2) for f in list(np.arange(0.0 / KP, 50 / KP + 0.05, 0.05))[1:]]):
        figname = 'heatmap_log_posterior_combined'

        def log_posterior_c_koff(c, koff, n, m):
            x = KON * c / koff
            meanN = KP * T * x / (1 + x)
            meanM = koff * T * x / (1 + x)
            varN = (KP * T * x) / (1 + x) + (2 * KP ** 2 * T * x) / (koff * (1 + x) ** 3) * (
                        1 + (np.exp(-T * koff * (1 + x)) - 1) / (T * koff * (1 + x)))
            varM = (koff * T * x) * (1 + x**2) / (1 + x)**3 + 2 * x**2 / (1 + x)**4 * (1 - np.exp(-koff * T * (1 + x)))
            covNM = (KP * T * x * (1-x)) / (1+x)**3 * (1 + (np.exp(-koff * T * (1+x)) - 1) / (T * koff * (1 + x)))
            sigma = np.array([[varN, covNM], [covNM, varM]])
            sigma_inv = np.linalg.inv(sigma)
            q = np.array([n - meanN, m - meanM])
            post = -0.5 * np.transpose(q).dot(sigma_inv).dot(q) - 0.5 * np.log((2 * 3.14159)**2 * np.linalg.det(sigma))
            return post

        arr = np.zeros((len(koffrange), len(crange)))
        for i, koffval in enumerate(koffrange):
            for j, cval in enumerate(crange):
                arr[i, j] = log_posterior_c_koff(cval * KP/KON, koffval * KP, nobs, mobs)

        label = r'ln $P(c, k_{\mathrm{off}}|n, m)$'
        fig, ax = plot_heatmap(arr, crange, koffrange, figname, label, save=False, log_norm=False,
                               levels=list(range(-1000, 0, 100)), vmin=-1000, vmax=0, contour_color='w', contour_linewidths=0.5)
        ax.set_xlabel(r'$k_{\mathrm{on}}c/k_{p}$', fontsize=FS)
        ax.set_ylabel(r'$k_{\mathrm{off}}/k_{p}$', fontsize=FS)

        # Superimpose heuristic estimate
        def heuristic_estimate_c(n, m):
            c_est = (KP/KON) * m / (KP * T - n)
            return c_est
        def heuristic_estimate_koff(n, m):
            koff_est = KP * m / n
            return koff_est

        # rescale onto weird axes scale which do not run from min to max of koffrange and crange anymore
        ypts = np.linspace(float(ax.get_yticklabels()[0].get_text()),
                           float(ax.get_yticklabels()[-1].get_text()), len(koffrange))
        koff_star = [ax.get_ylim()[-1] * (heuristic_estimate_koff(nobs, mobs) / KP) / max(koffrange) for _ in crange]
        ax.plot(list(np.linspace(ax.get_xlim()[0], ax.get_xlim()[-1], len(koffrange)))[:-1],
                [i + ax.get_ylim()[0] for i in koff_star][:-1], 'k--')

        c_star = [ax.get_xlim()[-1] * (heuristic_estimate_c(nobs, mobs) * KON/KP) / max(crange) for _ in koffrange]
        ax.plot([i + ax.get_xlim()[0] for i in c_star][:-1],
                list(np.linspace(ax.get_ylim()[0], ax.get_ylim()[-1], len(crange)))[:-1],
                'k--')

        # save figure
        fig.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
        fig.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    # ------------------
    # Repeat for Model 3
    # ------------------
    def KPR1_plot(crange=[], koffrange=[]):
        figname = 'heatmap_log_posterior_KPR1'

        def log_posterior_c_koff(c, koff, n, m):
            x = KON * c / koff
            meanN = (KF / (KF + koff)) * KP * T * x / (1 + x)
            meanM = (KF / (KF + koff)) * koff * T * x / (1 + x)
            varN = -((KF * KP * T * x)**2 /((KF + koff)**2 * (1 + x)**2)) + \
                   (KF*KP*T*x)/((KF + koff)*(1 + x)) - \
                   (KF * KP**2 * x* ((KF - koff*x)* (2*koff**3*(-1 + koff*T)*(1 + x)**3 + \
                                                KF**3*(-2 + 2*koff*T*(1 + x) + koff**2*T**2*x*(1 + x)**2) + \
                                                2*KF**2*koff*(-3 - x + koff**2*T**2*x*(1 + x)**2 +\
                                                              koff*T*(3 + 4*x + x**2)) +\
                                                KF*koff**2*(koff**2*T**2*x*(1 + x)**2 - 2*(3 + 3*x + x**2) +\
                                                            2*koff*T*(3 + 6*x + 4*x**2 + x**3))) -\
                                2*koff**4*x*(1 + x)**3*np.exp(-((KF + koff)*T)) +\
                                2*KF*(KF + koff)**3*np.exp(-(koff*T*(1 + x))))) / (koff**2 * (KF + koff)**4 * (1 + x)**4 * (-KF + koff*x))
            varM = -((KF**2*koff**2*T**2*x**2)/((KF + koff)**2*(1 + x)**2)) + (KF*koff*T*x)/((KF + koff)*(1 + x)) + \
                   (KF**2*x**2*((KF - koff*x)*(KF**2*(2 - 2*koff*T*(1 + x) + koff**2*T**2*(1 + x)**2) +\
                                               koff**2*(koff**2*T**2*(1 + x)**2 - 2*koff*T*(2 + 3*x + x**2) +\
                                                        2*(3 + 3*x + x**2)) +\
                                               2*KF*koff*(3 + x + koff**2*T**2*(1 + x)**2 -\
                                                          koff*T*(3 + 4*x + x**2))) +\
                                2*koff**3*(1 + x)**3*np.exp(-((KF + koff)*T)) -\
                                2*(KF + koff)**3*np.exp(-(koff*T*(1 + x))))) / ((KF + koff)**4*(1 + x)**4*(KF - koff*x))
            covNM = (KF*KP*x*(-(KF*koff**2*(KF + koff)**2*T**2*x*(1 + x)**2) - \
                              ((KF - koff*x)*(koff**3*(-1 + koff*T)*(1 + x)**3 +\
                                              KF*koff**2*(-3 + 2*x**2 + x**3 + koff**2*T**2*x*(1 + x)**2 +\
                                                          koff*T*(3 + 4*x + x**2)) +\
                                              KF**3*(-1 + x + koff**2*T**2*x*(1 + x)**2 + koff*(T - T*x**2)) +\
                                              KF**2*koff*(-3 + 2*x + x**2 + 2*koff**2*T**2*x*(1 + x)**2 +\
                                                          koff*T*(3 + x - 3*x**2 - x**3))) -\
                               koff**3*(-KF + koff)*x*(1 + x)**3*np.exp(-((KF + koff)*T)) -\
                               KF*(KF + koff)**3*(-1 + x)*np.exp(-(koff*T*(1 + x)))) / (-KF + koff*x)))/(koff*(KF + koff)**4*(1 + x)**4)
            sigma = np.array([[varN, covNM], [covNM, varM]])
            sigma_inv = np.linalg.inv(sigma)
            q = np.array([n - meanN, m - meanM])
            post = -0.5 * np.transpose(q).dot(sigma_inv).dot(q) - 0.5 * np.log(
                (2 * 3.14159) ** 2 * np.linalg.det(sigma))
            return post

        arr = np.zeros((len(koffrange), len(crange)))
        for i, koffval in enumerate(koffrange):
            for j, cval in enumerate(crange):
                arr[i, j] = log_posterior_c_koff(cval, koffval, nobs, mobs)

        label = r'$ln(P(c, k_{\mathrm{off}}|n, m))$'
        fig, ax = plot_heatmap(arr, crange, koffrange, figname, label, save=False, log_norm=False,
                               levels=list(range(-200, 0, 10)), vmin=-200, vmax=0, contour_color='w',
                               contour_linewidths=0.5)
        ax.set_xlabel(r'$k_{\mathrm{on}}c/k_{p}$', fontsize=FS)
        ax.set_ylabel(r'$k_{\mathrm{off}}/k_{p}$', fontsize=FS)

        # Superimpose heuristic estimate
        def heuristic_estimate_c(n, m):
            c_est = - KP * m * (KP * m + KF * n) / (KON * n * (KP * m + KF * n - KF * KP * T))
            return c_est

        def heuristic_estimate_koff(n, m):
            koff_est = KP * m / n
            return koff_est

        # rescale onto weird axes scale which do not run from min to max of koffrange and crange anymore
        xpts = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], len(crange))
        koff_star = [ax.get_ylim()[1] / max(koffrange) * heuristic_estimate_koff(nobs, mobs) + ax.get_ylim()[0]
                     for _ in range(len(xpts))]
        ax.plot(xpts, koff_star, 'k--')

        ypts = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], len(koffrange))
        c_star = [ax.get_xlim()[1] / max(crange) * heuristic_estimate_c(nobs, mobs) + ax.get_xlim()[0]
                  for _ in range(len(ypts))]

        ax.plot(c_star, ypts, 'k--')
        # save figure
        fig.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
        #fig.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    # ------------------
    # Repeat for KPR 2
    # ------------------
    n1obs = nobs * 0.1
    n2obs = mobs * 0.5
    def KPR2_plot(crange=[], koffrange=[]):
        figname = 'heatmap_log_posterior_KPR2'

        def log_posterior_c_koff(c, koff, n1, n2):
            x = KON * c / koff
            meanN1 = koff * KP * T * x / ((1 + x) * (KF + koff))
            meanN2 = KF * KP * T * x / ((1 + x) * (KF + koff))
            varN1 = -((koff**2*KP**2*T**2*x**2)/(KF + koff + KF*x + koff*x)**2) + \
                    (koff*KP*T*x)/(KF + koff + KF*x + koff*x) + \
                    (KP**2*x*((KF - koff*x)*(koff**2*(-2 + 2*koff*T*(1 + x) + koff**2*T**2*x*(1 + x)**2) +\
                                             KF**2*(2*x + koff**2*T**2*x*(1 + x)**2 +\
                                                    2*koff*T*(1 + 2*x + 2*x**2 + x**3)) +\
                                             2*KF*koff*(-1 - 2*x**2 - x**3 + koff**2*T**2*x*(1 + x)**2 +\
                                                        koff*T*(2 + 3*x + 2*x**2 + x**3))) -\
                              2*KF*koff*(-KF + koff*(-1 + x))*(1 + x)**3*np.exp(-((KF + koff)*T)) -\
                              2*(KF + koff)**3*x*np.exp(-(koff*T*(1 + x)))))/((KF + koff)**4*(1 + x)**4*(KF - koff*x))

            varN2 = -((KF**2*KP**2*T**2*x**2)/(KF + koff + KF*x + koff*x)**2) + \
                    (KF*KP*T*x)/(KF + koff + KF*x + koff*x) -\
                    (KF*KP**2*x*((KF - koff*x) * (2*koff**3*(-1 + koff*T)*(1 + x)**3 +\
                                                  KF**3*(-2 + 2*koff*T*(1 + x) + koff**2*T**2*x*(1 + x)**2) +\
                                                  2*KF**2*koff*(-3 - x + koff**2*T**2*x*(1 + x)**2 +\
                                                                koff*T*(3 + 4*x + x**2)) +\
                                                  KF*koff**2*(koff**2*T**2*x*(1 + x)**2 - 2*(3 + 3*x + x**2) +\
                                                              2*koff*T*(3 + 6*x + 4*x**2 + x**3))) -\
                                 2*koff**4*x*(1 + x)**3*np.exp(-((KF + koff)*T)) +\
                                 2*KF*(KF + koff)**3*np.exp(-(koff*T*(1 + x)))))/(koff**2*(KF + koff)**4*(1 + x)**4*(-KF + koff*x))

            covN1N2 = (KF*KP**2*x*(-(koff**2*(KF + koff)**2*T**2*x*(1 + x)**2) + \
                                 (-((KF - koff*x)*(KF**2*(-1 + x + koff**2*T**2*x*(1 + x)**2 + koff*(T - T*x**2)) + \
                                                   KF*koff*(-3 + 2*x + x**2 + 2*koff**2*T**2*x*(1 + x)**2 -\
                                                            2*koff*T*(-1 + x + 3*x**2 + x**3)) +\
                                                   koff**2*(-2 + 3*x + 5*x**2 + 2*x**3 +\
                                                            koff**2*T**2*x*(1 + x)**2 -\
                                                            koff*T*(-1 + 2*x + 5*x**2 + 2*x**3)))) -\
                                  koff**2*(1 + x)**3*(-KF + koff*(-1 + 2*x))*\
                                  np.exp(-((KF + koff)*T)) +\
                                  (KF + koff)**3*(-1 + x)*np.exp(-(koff*T*(1 + x))))/(-KF + koff*x)))/\
                    (koff*(KF + koff)**4*(1 + x)**4)
            sigma = np.array([[varN1, covN1N2], [covN1N2, varN2]])
            sigma_inv = np.linalg.inv(sigma)
            q = np.array([n1 - meanN1, n2 - meanN2])
            post = -0.5 * np.transpose(q).dot(sigma_inv).dot(q) - 0.5 * np.log(
                (2 * 3.14159) ** 2 * np.linalg.det(sigma))
            return post

        arr = np.zeros((len(koffrange), len(crange)))
        for i, koffval in enumerate(koffrange):
            for j, cval in enumerate(crange):
                arr[i, j] = log_posterior_c_koff(cval, koffval, n1obs, n2obs)

        label = r'$ln(P(c, k_{\mathrm{off}}|n_1, n_2))$'
        fig, ax = plot_heatmap(arr, crange, koffrange, figname, label, save=False, log_norm=False,
                               levels=list(range(-200, 0, 10)), vmin=-200, vmax=0, contour_color='w',
                               contour_linewidths=0.5)
        ax.set_xlabel(r'$k_{\mathrm{on}}c/k_{p}$', fontsize=FS)
        ax.set_ylabel(r'$k_{\mathrm{off}}/k_{p}$', fontsize=FS)

        # Superimpose heuristic estimate
        def heuristic_estimate_c(n1, n2):
            c_est = -KF * n1 * (n1 + n2) / (KON * n2 * (n1 + n2 - KP * T))
            return c_est

        def heuristic_estimate_koff(n1, n2):
            koff_est = KF * n1 / n2
            return koff_est

        # rescale onto weird axes scale which do not run from min to max of koffrange and crange anymore
        xpts = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], len(crange))
        koff_star = [ax.get_ylim()[1] / max(koffrange) * heuristic_estimate_koff(n1obs, n2obs) + ax.get_ylim()[0]
                     for _ in range(len(xpts))]
        ax.plot(xpts, koff_star, 'k--')

        ypts = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], len(koffrange))
        c_star = [ax.get_xlim()[1] / max(crange) * heuristic_estimate_c(n1obs, n2obs) + ax.get_xlim()[0]
                  for _ in range(len(ypts))]

        ax.plot(c_star, ypts, 'k--')
        # save figure
        fig.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
        # fig.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    mode1_plot()
    mode2_plot()
    #KPR1_plot()
    #KPR2_plot()
    return 0

def heatmap_ratio(eqn1, eqn2, label, filename, log_norm, crange=CRANGE, koffrange=KOFFRANGE, dedim=False, contour_args=None):
    print("Plotting equations:",eqn1,eqn2)
    # creating arrays to plot
    arr = np.zeros((len(koffrange), len(crange)))
    arr1 = np.zeros((len(koffrange), len(crange)))
    arr2 = np.zeros((len(koffrange), len(crange)))
    # filling arrays with ratio
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = eqn1(cval, koffval)/eqn2(cval, koffval)
            arr1[i, j] = eqn1(cval, koffval) # these are unnecessary, however good for checks
            arr2[i, j] = eqn2(cval, koffval)
    print(np.min(arr),np.max(arr))

    # call heatmap plotting
    plot_heatmap(arr, crange, koffrange, filename, label, log_norm=log_norm, dedim=dedim, **contour_args)
    return 0

def plot_dictionary_ratio(dict_ratio, dedim=False, subdir1='heatmaps', longtime=False, highG=False, contour_args=None):
    # creating subdirectory
    if not os.path.exists(DIR_OUTPUT + os.sep + subdir1 + os.sep + dict_ratio['subdir2']):
        os.makedirs(DIR_OUTPUT + os.sep + subdir1 + os.sep + dict_ratio['subdir2'])

    # calling creation of arrays/plotting for all plots in subdirectory
    for ratio in dict_ratio['plots']:
        heatmap_ratio(dict_ratio['plots'][ratio]['num'], dict_ratio['plots'][ratio]['denom'], dict_ratio['plots'][ratio]['label'], subdir1 + os.sep + dict_ratio['subdir2'] + os.sep + ratio, dict_ratio['log'], dedim=dedim, contour_args=contour_args)
    return 0

def heatmap_one_equation(equation, label, filename, log_norm, crange=CRANGE, koffrange=KOFFRANGE, dedim=False, kon=KON, T=T, KF=KF, KP=KP, contour_args={'None' : None}):
    print("Plotting equation:",equation)
    # create array to plot
    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = equation(cval, koffval)+10**(-10)

    # call heatmap plotting
    plot_heatmap(arr, crange, koffrange, filename, label, dedim=dedim, log_norm=log_norm, **contour_args)
    return 0

def plot_dictionary_one_equation(dict_eqns, subdir1='heatmaps', longtime=False, highG=False, dedim=False, contour_args={'None' : None}):
    # create subdirectory
    if not os.path.exists(DIR_OUTPUT + os.sep + subdir1 + os.sep + dict_eqns['subdir2']):
        os.makedirs(DIR_OUTPUT + os.sep + subdir1 + os.sep + dict_eqns['subdir2'])

    # call plotting function for each plot in dictionary
    for name_eqn in dict_eqns['plots']:
        heatmap_one_equation(dict_eqns['plots'][name_eqn]['eqn'], dict_eqns['plots'][name_eqn]['label'], subdir1 + os.sep + dict_eqns['subdir2'] + os.sep + name_eqn, dict_eqns['log'], dedim=dedim, contour_args=contour_args)
    return 0

def dk_plotting():
    """
    Duncan you can still run any/all of these by running dk_plotting in main
    """
    heatmap_mode1_error_x(make_heatmap=False, make_panel=True)

    #heatmap_mode1_error_x()
    figure_2_combined_cross_sections()

    #heatmap_combined_error_c()
    #heatmap_combined_error_koff()
    #heatmap_kpr_error_c()
    #heatmap_kpr_error_koff()

    #heatmap_kpr2_error_c()
    #heatmap_kpr2_error_koff()

    ctildePosterior = [truncate(f, 3) for f in list(np.arange(0.0 * KON / KP, 5.0 * KON / KP + 0.005, 0.005))[1:]]
    kofftildePosterior = [truncate(f, 2) for f in list(np.arange(0.0 / KP, 50.0 / KP + 0.05, 0.05))[1:]]

    #heatmap_figure_4()

    return 0


def custom_ratio_diagram(subdir1='heatmaps', subdir2='', contour_args=None):
    if not os.path.exists(DIR_OUTPUT + os.sep + subdir1 + os.sep + subdir2):
        os.makedirs(DIR_OUTPUT + os.sep + subdir1 + os.sep + subdir2)

    eps = 10**(-10)
    # params
    fix_C = 1
    fix_KON = 1
    fix_KP = 10
    fix_T = 1000

    ax1label = r'$k_{f}/k_{p}$'
    ax1range = np.logspace(-2, 4, TOTAL_POINTS_KOFF) / KP
    ax2label = r'$k_{\mathrm{off}}/k_{p}$'
    ax2range = ZRANGE
    xy_label_force = [ax1label, ax2label]

    ratio_arr = np.zeros((len(ax2range), len(ax1range)))
    for i, ival in enumerate(ax2range):
        for j, jval in enumerate(ax1range):
            num = eqns.DetSigmacrlb3(fix_C, ival, kon=fix_KON, T=fix_T, KP=fix_KP, KF=jval)
            den = eqns.DetSigmacrlb2(fix_C, ival, kon=fix_KON, T=fix_T, KP=fix_KP, KF=jval)
            ratio_arr[i, j] = num/den + eps  #equation(cval, koffval)
    print('custom_ratio_diagram -- min max:', np.min(ratio_arr), np.max(ratio_arr))
    cbar_label = r'det$(\Sigma^{KPR})$/det$(\Sigma)$'
    plot_heatmap(ratio_arr, ax1range, ax2range, 'custom_ratio_det3det2_notrace', cbar_label,
                 xy_label_force=xy_label_force, show=True, save=True, log_norm=True, dedim=False, **contour_args)
    return

def plot_1E_and_2B():
    figname1 = 'Figure_1E'
    figname2 = 'Figure_2B'
    font_size = FS

    mpl.rcParams['xtick.labelsize'] = font_size
    mpl.rcParams['ytick.labelsize'] = font_size

    # Fig 1E
    koffrange = KOFFRANGE
    c = c0 # = KP/KON
    vecErrorX1 = eqns.dedimRelErrorX1NoTrace(c, koffrange);

    # Fig 2B
    crange = CRANGE*10
    koff = 1.0
    vecRelErrorEst2C = eqns.dedimRelErrC2NoTrace(crange, koff);
    vecRelErrorEst2K = eqns.dedimRelErrK2NoTrace(crange, koff);

    # plot
    fig_width_1E = 2.6
    plt.figure(figsize=(fig_width_1E, fig_width_1E*0.8))
    ax = plt.gca()

    plt.plot(koffrange/KP,vecErrorX1/N, color='purple')
    ax.set_xlabel(r'$k_{\mathrm{off}}/k_{p}$', fontsize=font_size)
    ax.set_xscale('log')
    ax.set_ylabel(r'$\frac{k_{p}t}{N} \frac{\langle\delta x^{2}\rangle} {x^{2}}$',fontsize=font_size)
    ax.set_xlim([1E-1, 1E2])
    plt.ylim([0, 1])
    plt.savefig(DIR_OUTPUT + os.sep + figname1 + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname1 + '.eps')

    #plt.show()
    plt.close()

    # plot
    fig_width_2B = 2.2
    plt.figure(figsize=(fig_width_2B, fig_width_2B*0.8))
    ax = plt.gca()

    plt.plot(crange*KON/KP,vecRelErrorEst2C/N, color='purple', label=r'$c$')
    plt.plot(crange*KON/KP,vecRelErrorEst2K/N, color='orangered', label=r'$k_{\mathrm{off}}$')
    plt.legend()
    ax.set_xlabel(r'$k_{\mathrm{on}}c/k_{p}$', fontsize=font_size)
    ax.set_xscale('log')
    ax.set_ylabel(r'$\frac{k_{p}t}{N} \frac{\langle\delta (\cdot)^{2}\rangle}{(\cdot)^{2}}$',fontsize=font_size)
    #ax.set_xlim([1E-4, 1E1])
    plt.ylim([0, 50])
    plt.savefig(DIR_OUTPUT + os.sep + figname2 + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname2 + '.eps')

    #plt.show()
    plt.close()

    """
    ax.set_ylabel(r'$\frac{k_{p}t}{N} \frac{\langle\delta (\cdot)^{2}\rangle$}{$(\cdot)^{2}}$',fontsize=FS)
    if label_style == 0:
        ln1 = ax1.plot(curve1['xpts'], curve1['ypts'], color=cs['simple_fisher'], label=r'$c$', zorder=1)
        ln2 = ax2.plot(curve2['xpts'], curve2['ypts'], color=cs['heuristic'], label=r'$k_{\mathrm{off}}$', zorder=1)
        #plt.title('Mode 2: MLE relative error comparison\n' + r'($\tilde{c}_0=10$, $\alpha=1 \times 10^4$, $k_{p}=10$)')

    elif label_style == 1:
        ln1 = ax1.plot(curve1['xpts'], curve1['ypts'], color=cs['simple_fisher'], label=r'$\delta c^{2}/c^{2}$', zorder=1)
        ln2 = ax2.plot(curve2['xpts'], curve2['ypts'], color=cs['heuristic'],label=r'$\delta k_{\mathrm{off}}^{2}/k_{\mathrm{off}}^{2}$', zorder=1)
        plt.title('Mode 2: MLE relative error comparison\n' + r'($k_p=10$, $t=100$, $k_{\mathrm{off}}=k_{\mathrm{on}}=1$)')
        plt.ylabel(r'$\langle\delta (\cdot)^{2}\rangle$/$(\cdot)^{2}$')

    # axis
    ax1.set_xlabel(r'$k_{\mathrm{on}}c/k_{p}$')
    ax1.set_ylabel(r'$k_{p}t \langle\delta (\cdot)^{2}\rangle$/$(\cdot)^{2}$')
    #ax2.set_xlabel(r'$k_{\mathrm{off}}$')

    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax1.set_xlim([1E-2, 1E2])
    ax2.set_xlim([1E-2, 1E2])
    plt.ylim([0, 0.01*alpha])

    lns = ln1 + ln2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs)

    plt.tight_layout()
    # save figure
    """

def plot_1F():
    figname = 'heatmap_log_posterior_mode1'

    nobs = 0.1 * KP * T # 10

    LOG_START_KOFF_1F = -1
    LOG_END_KOFF_1F = 3
    TOTAL_POINTS_KOFF_1F = (LOG_END_KOFF_1F - LOG_START_KOFF_1F) * POINTS_BETWEEN_TICKS + 1
    crange = CRANGE
    koffrange = np.logspace(LOG_START_KOFF_1F, LOG_END_KOFF_1F, TOTAL_POINTS_KOFF_1F)

    def log_posterior_x(c, koff, n):
        if c==0:
            c=1E-16
        if koff==0:
            koff=1E-16
        x = KON * c / koff
        mu = KP * T * x / (1 + x)
        var = (KP * T * x)/(1 + x) + (2 * KP**2 * T * x)/(koff * (1 + x)**3) *(1 + (np.exp(-T * koff * (1 + x)) - 1)/(T * koff * (1 + x)))
        post = -(n - mu)**2 / (2 * var) - 0.5 * np.log(2 * 3.14159 * var)
        return post

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = log_posterior_x(cval, koffval, nobs)
    arr = -1.0 * arr
    label = r'-ln $P(c, k_{\mathrm{off}}|n)$'
    #log_contours = [-i for i in list(np.logspace(-3, 4,num=20))[::-1]]
    linear_contours = []#[truncate(np.min(-1.*arr), 1)*1.01]#list(range(-500, 0, 50))

    if np.all(arr) > 0.0:
        fig, ax = plot_heatmap(arr, crange*KON/KP, koffrange/KP, figname, '', save=False, log_norm=True,
                               levels=linear_contours, vmin=1, vmax=500, contour_color='k', contour_linewidths=0.5,
                               cmap_colour='viridis_r',fig_width=2.8)
    else:
        print("Not all posterior points evaluated to non-zero probability")
        return 1

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

    #ax.plot([len(crange)/2 for _ in range(len(koffrange))], range(len(koffrange)), 'k--')


    # save figure
    fig.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    fig.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


if __name__ == '__main__':
    """
    This is how I use the code generally.
    I create a dictionary with a bunch of equations I want to plot in dictionary_plotting.py (for any equations you find in there, they are in equations.py).
    These equations are all taken from mathematica or rescaled versions of these (see equations_txt2python.py for some rescaled equations, trace and eigenvalue, dedimRelErr)
    There are 2 types of dictionaries: ratios or 1 equation. Depending on what type of dictionary you create, you will want to use the respoective plotting function above (ether plot_dictionary_one_equation or plot_dictionary_ratio)
    You have to specify some arguments, check equations above, they should be obvious.
    You can create your own plotting dictionaries and equations! So much fun to be had!
    """
    #plot_1F()
    #plot_1E_and_2B()
    #exit()

    #dk_plotting()
    #exit()
    dictionary = pd.MAIN; dictionary_SI = pd.SI_ALT_RATIO
    want_dedim = True; subdir_2_use = 'heatmaps'

    contour_args = {'cmap_colour' : 'YlGnBu'}; contour_args_SI = {'cmap_colour' : 'PuBu', 'levels' : [1.01, 1.1, 10.0, 100.0, 1000.0], 'fmt' : '%.2f'}

    #plot_dictionary_one_equation(dictionary, subdir1=subdir_2_use, dedim=want_dedim, contour_args=contour_args)
    plot_dictionary_ratio(dictionary_SI, subdir1=subdir_2_use, dedim=want_dedim, contour_args=contour_args_SI)

    #custom_cmap_colour = 'YlGnBu' # 'YlGnBu' or 'pink_r'
    #contour_args_high = {'levels': [1 / (KP * T), 10 / (KP * T), 100 / (KP * T), 1000 / (KP * T), 1E4 / (KP * T)],
    #                     'cmap_colour': custom_cmap_colour,
    #                     'vmin': 1.0}

    #custom_ratio_diagram(contour_args=contour_args_SI)
    """
    Hi
    """
