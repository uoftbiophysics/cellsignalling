import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import plotting_dictionaries as pd

#from load_inputs import DATADICT
from settings import DIR_OUTPUT, DIR_INPUT
from settings import COLOR_SCHEME as cs

plt.style.use('parameters.mplstyle')  # particularIMporting

# unit params
KON = 1.
KP = 10.
T = 1000.
KF = 1.0

# de-dimensionalise params
c0 = KP/KON
alpha = KP*T

# plot params
FS = 20
SHOW = False

# axes
POINTS_BETWEEN_TICKS = 10
LOG_START_C = -1
LOG_END_C = 3
TOTAL_POINTS_C = (LOG_END_C - LOG_START_C) * POINTS_BETWEEN_TICKS + 1
CRANGE = np.logspace(LOG_START_C, LOG_END_C, TOTAL_POINTS_C)
LOG_START_KOFF = -1
LOG_END_KOFF = 3
TOTAL_POINTS_KOFF = (LOG_END_KOFF - LOG_START_KOFF) * POINTS_BETWEEN_TICKS + 1
KOFFRANGE = np.logspace(LOG_START_KOFF, LOG_END_KOFF, TOTAL_POINTS_KOFF)

CTILDERANGE = np.divide(list(CRANGE), c0)
ZRANGE = np.divide(list(KOFFRANGE), KP)

# global vmax and vmin
assert LOG_START_KOFF == -1
assert LOG_END_KOFF == 3
assert LOG_START_C == -1
assert LOG_END_C == 3
vmax_obs = 1000003020
vmin_obs = 0.00022000016528925625


def plot_heatmap(arr, crange, koffrange, fname, label, show=SHOW, save=True, log_norm=True, dedim=False, **kwargs):
    """
    crange: range of values for c
    koffrange: range of koff values
    fname: file saved as pdf and eps format with this name
    show: shows plots
    save: saves plots
    log_norm: heatmap is log, should almost always be the case
    dedim: makes the axis dedimensionalized. scales them (see by how much below)
    kwargs: a variety of keyword args are used below. They are mostly used for contour plot lines if I understand correctly. Using Duncan's default ones mostly, but we can talk about it.
    """
    xy_label = [r'${c}$', r'${k}_{off}$']
    # default parameters
    if 'levels' in kwargs.keys(): levels = kwargs['levels']
    else: levels = [1, 10, 100, 1000, 1E4]

    if 'vmin' in kwargs.keys(): vmin = kwargs['vmin']
    else: vmin = np.min(arr)

    if 'vmax' in kwargs.keys(): vmax = kwargs['vmax']
    else: vmax = np.max(arr)

    if 'contour_linestyle' in kwargs.keys(): contour_linestyle = kwargs['contour_linestyle']
    else: contour_linestyle = 'solid'

    if 'contour_color' in kwargs.keys(): contour_color = kwargs['contour_color']
    else: contour_color = None

    if 'contour_linewidths' in kwargs.keys(): contour_lindewidths = kwargs['contour_linewidths']
    else: contour_lindewidths = None

    if dedim == True:
        # if flag is true, this is how we scale the axis. Simple.
        crange = crange*KON/KP;
        koffrange = koffrange/KP;
        xy_label = [r'$k_{on}c/k_{p}$', r'$k_{off}/k_{p}$'];

    # TODO change colour scheme, see https://matplotlib.org/examples/color/colormaps_reference.html
    # TODO fix ticks randomly disappearing on colourbar + flip colourbar minor ticks or remove?
    """
    Colours viridis, YlGnBu, terrain, plasma
    """
    #print 'arr limits:', np.min(arr), np.max(arr)
    # plot setup
    f = plt.figure()
    if log_norm:
        imshow_kw = {'cmap': 'YlGnBu', 'aspect': None, 'vmin': vmin, 'vmax': vmax, 'norm': mpl.colors.LogNorm()}
        #imshow_kw = {'cmap': 'PuBu', 'aspect': None, 'vmin': vmin, 'vmax': vmax, 'norm': mpl.colors.LogNorm()}
    else:
        imshow_kw = {'cmap': 'YlGnBu', 'aspect': None, 'vmin': vmin, 'vmax': vmax}
        #imshow_kw = {'cmap': 'PuBu', 'aspect': None, 'vmin': vmin, 'vmax': vmax}
    im = plt.imshow(arr, interpolation='spline36', **imshow_kw)

    # axes setup
    fig = plt.gcf(); ax = plt.gca()

    # method 1
    ax.set_xticks([i for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS == 0])
    ax.set_yticks([i for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0])
    ax.set_xticklabels([r'$10^{%d}$' % np.log10(cval) for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
    ax.set_yticklabels([r'$10^{%d}$' % np.log10(kval) for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
    """
    if log_norm:
        ax.set_xticklabels([r'$10^{%d}$' % np.log10(cval) for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
        ax.set_yticklabels([r'$10^{%d}$' % np.log10(kval) for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
    else:
        ax.set_xticklabels(["{:10.2f}".format(cval) for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
        ax.set_yticklabels(["{:10.2f}".format(kval) for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS==0], fontsize=FS)
    """
    ax.invert_yaxis()
    ax.set_xlabel(xy_label[0], fontsize=FS); ax.set_ylabel(xy_label[1], fontsize=FS)

    # create colorbar
    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=FS, labelpad=20); cbar.ax.tick_params(labelsize=FS)
    cbar.ax.minorticks_off(); cbar.update_ticks()

    # TODO IDK why do ticks hide sometimes?
    #for t in cbar.ax.get_yticklabels(): print(t.get_text())
    # contour line for value 1.0
    plt.contour(arr, levels=levels, linestyles=contour_linestyle, colors=contour_color, linewidths=contour_lindewidths)

    # save
    if save == True:
        plt.savefig(DIR_OUTPUT + os.sep + fname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + fname + '.eps')
    if show:
        plt.show()

    plt.close()
    #return fig, ax
    return 0


def heatmap_mode1_error_x(crange=CTILDERANGE, koffrange=ZRANGE, make_heatmap=True, make_panel=False,
                          scale_factor=alpha, label_style=0):
    if make_heatmap == True:
        def mode1_error_x(ctilde, z, scale_factor=scale_factor):
            c = ctilde * c0
            koff = z * KP
            x = c * KON / koff
            val = (1 + x) / (KP * T * x) * ((1 + x) ** 2 + 2 * KP / koff)
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
                return scale_factor * val

            arr = [mode1_error_c(cval, 1) for cval in crange]
            return dict({'xpts': crange, 'ypts':arr})

        figname = 'mode1_error_c_cross_section'
        curve1 = cross_section_mode1_error_c()
        # plot
        plt.figure(figsize=(4, 3))
        #plt.figure()
        plt.plot(curve1['xpts'], curve1['ypts'], color=cs['simple_fisher'], label='Simple Fisher', zorder=1)
        # axis
        if label_style == 0:
            plt.title('Mode 1: MLE relative error comparison \n'+r'($\tilde{c}_0=10$, $\alpha=1 \times 10^4$, $k_{p}=10$)')
            plt.ylabel(r'$\alpha \langle\delta c^{2}\rangle$/$c^{2}$')
        elif label_style == 1:
            plt.title('Mode 1: MLE relative error comparison \n($k_p=10$, $t=100$, $k_{off}=k_{on}=1$)')
            plt.ylabel(r'$\langle\delta c^{2}\rangle$/$c^{2}$')

        plt.xlabel(r'$\tilde{c}_0$')

        plt.gca().set_xscale('log')
        plt.xlim([CTILDERANGE[0], CTILDERANGE[-1]])
        plt.ylim([0, alpha])
        plt.xscale('log')
        #plt.legend()
        # save figure
        plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
        plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    return

def heatmap_combined_error_c(crange=CTILDERANGE, koffrange=ZRANGE,
                             scale_factor=alpha, label_style=0):
    def combined_error_c(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        x = KON * c / koff
        num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff**2 * T * x ** 2 * (1 + x) ** 3
        den = koff ** 2 * KP * T**2 * x * (1 + x) ** 2
        val = num / den
        return scale_factor * val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = combined_error_c(cval, koffval)

    if label_style == 0:
        label = r'$\alpha \langle\delta c^{2}\rangle$/$c^{2}$'
    elif label_style == 1:
        label = r'$\langle\delta c^{2}\rangle$/$c^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_combined_heuristic_error_c', label)
    return

def heatmap_combined_error_koff(crange=CRANGE, koffrange=KOFFRANGE,
                             scale_factor=alpha, label_style=0):

    def combined_error_koff(ctilde, z, scale_factor=scale_factor):
        c = ctilde * c0
        koff = z * KP
        x = c * KON / koff
        num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff**2 * T * (1 + x) ** 3
        den = koff ** 2 * KP * T**2 * x * (1 + x) ** 2
        val = num / den
        return scale_factor * val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = combined_error_koff(cval, koffval)

    if label_style == 0:
        label = r'$\alpha \langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
    elif label_style == 1:
        label = r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
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
        ln1 = ax1.plot(curve1['xpts'], curve1['ypts'], color=cs['simple_fisher'], label=r'$\alpha \delta c^{2}/c^{2}$', zorder=1)
        ln2 = ax2.plot(curve2['xpts'], curve2['ypts'], color=cs['heuristic'], label=r'$\alpha \delta k_{off}^{2}/k_{off}^{2}$', zorder=1)
        plt.title('Mode 2: MLE relative error comparison\n' + r'($\tilde{c}_0=10$, $\alpha=1 \times 10^4$, $k_{p}=10$)')
        plt.ylabel(r'$\alpha \langle\delta (\cdot)^{2}\rangle$/$(\cdot)^{2}$')
    elif label_style == 1:
        ln1 = ax1.plot(curve1['xpts'], curve1['ypts'], color=cs['simple_fisher'], label=r'$\delta c^{2}/c^{2}$', zorder=1)
        ln2 = ax2.plot(curve2['xpts'], curve2['ypts'], color=cs['heuristic'],label=r'$\delta k_{off}^{2}/k_{off}^{2}$', zorder=1)
        plt.title('Mode 2: MLE relative error comparison\n' + r'($k_p=10$, $t=100$, $k_{off}=k_{on}=1$)')
        plt.ylabel(r'$\langle\delta (\cdot)^{2}\rangle$/$(\cdot)^{2}$')

    # axis
    ax1.set_xlabel(r'$c$')
    #ax2.set_xlabel(r'$k_{off}$')

    ax1.set_xscale('log')
    ax2.set_xscale('log')
    #ax1.set_xlim([1E-3, 1E3])
    #ax2.set_xlim([1E-3, 1E3])
    plt.ylim([0, 0.1*alpha])

    lns = ln1 + ln2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs)

    plt.tight_layout()
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

def heatmap_kpr_error_c(crange=CTILDERANGE, koffrange=ZRANGE,
                        scale_factor=alpha, label_style=0):

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
        label = r'$\alpha \langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
    elif label_style == 1:
        label = r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
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
        label = r'$\alpha \langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
    elif label_style == 1:
        label = r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
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

def heatmap_figure_4(crange=CRANGE, koffrange=KOFFRANGE):
    nobs = 0.1 * KP * T # 100
    mobs = 0.15 * KP * T # 150
    def mode1_plot():
        figname = 'heatmap_log_posterior_mode1'
        def log_posterior_x(c, koff, n):
            x = KON * c / koff
            mu = KP * T * x / (1 + x)
            var = (KP * T * x)/(1 + x) + (2 * KP**2 * T * x)/(koff * (1 + x)**3) *(1 + (np.exp(-T * koff * (1 + x)) - 1)/(T * koff * (1 + x)))
            post = -(n - mu)**2 / (2 * var) - 0.5 * np.log(2 * 3.14159 * var)
            return post

        arr = np.zeros((len(koffrange), len(crange)))
        for i, koffval in enumerate(koffrange):
            for j, cval in enumerate(crange):
                arr[i, j] = log_posterior_x(cval, koffval, nobs)

        label = r'$ln(P(x|n))$'
        fig, ax = plot_heatmap(arr, crange, koffrange, figname, label, save=False, log_norm=False,
                               levels=list(range(-50, 0, 5)), vmin=-60, vmax=0, contour_color='w', contour_linewidths=0.5)
        ax.set_xlabel(r'$c$', fontsize=FS)
        ax.set_ylabel(r'$k_{off}$', fontsize=FS)

        # Superimpose heuristic estimate
        def heuristic_estimate(c, n):
            koff_est = ((KP * T - n) * KON * c) / n
            return koff_est

        estimate_line = []
        for c in crange:
            estimate_line.append(heuristic_estimate(c, nobs))
        # rescale onto weird axes scale which do not run from min to max of koffrange and crange anymore
        xpts = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], len(crange))
        ypts = [ax.get_ylim()[1] * koff / max(koffrange) for koff in estimate_line]
        ax.plot(xpts, ypts, 'k--')
        # save figure
        fig.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
        fig.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    # ------------------
    # Repeat for Model 2
    # ------------------
    def mode2_plot():
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
                arr[i, j] = log_posterior_c_koff(cval, koffval, nobs, mobs)

        label = r'$ln(P(c, k_{off}|n, m))$'
        fig, ax = plot_heatmap(arr, crange, koffrange, figname, label, save=False, log_norm=False,
                               levels=list(range(-100, 0, 10)), vmin=-200, vmax=0, contour_color='w', contour_linewidths=0.5)
        ax.set_xlabel(r'$c$', fontsize=FS)
        ax.set_ylabel(r'$k_{off}$', fontsize=FS)

        # Superimpose heuristic estimate
        def heuristic_estimate_c(n, m):
            c_est = KP * m / (KP * T - n)
            return c_est
        def heuristic_estimate_koff(n, m):
            koff_est = KP * m / n
            return koff_est
        # rescale onto weird axes scale which do not run from min to max of koffrange and crange anymore
        xpts = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], len(crange))
        ypts = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], len(koffrange))

        koff_star = [ax.get_ylim()[1]/max(koffrange) * heuristic_estimate_koff(nobs, mobs) for _ in range(len(xpts))]
        c_star = [ax.get_xlim()[1] / max(crange) * heuristic_estimate_c(nobs, mobs) for _ in range(len(ypts))]
        ax.plot(xpts, koff_star, 'k--')
        ax.plot(c_star, ypts, 'k--')
        # save figure
        fig.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
        fig.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    # ------------------
    # Repeat for Model 3
    # ------------------
    def KPR1_plot(crange=crange, koffrange=koffrange):
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

        label = r'$ln(P(c, k_{off}|n, m))$'
        fig, ax = plot_heatmap(arr, crange, koffrange, figname, label, save=False, log_norm=False,
                               levels=list(range(-200, 0, 10)), vmin=-200, vmax=0, contour_color='w',
                               contour_linewidths=0.5)
        ax.set_xlabel(r'$c$', fontsize=FS)
        ax.set_ylabel(r'$k_{off}$', fontsize=FS)

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
    def KPR2_plot(crange=crange, koffrange=koffrange):
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

        label = r'$ln(P(c, k_{off}|n_1, n_2))$'
        fig, ax = plot_heatmap(arr, crange, koffrange, figname, label, save=False, log_norm=False,
                               levels=list(range(-200, 0, 10)), vmin=-200, vmax=0, contour_color='w',
                               contour_linewidths=0.5)
        ax.set_xlabel(r'$c$', fontsize=FS)
        ax.set_ylabel(r'$k_{off}$', fontsize=FS)

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

    #mode1_plot()
    #mode2_plot()
    #KPR1_plot()
    KPR2_plot()
    return

def heatmap_ratio(eqn1, eqn2, label, filename, log_norm, crange=CRANGE, koffrange=KOFFRANGE, dedim=False, contour_args=None):
    print("Plotting equations:",eqn1,eqn2)
    arr = np.zeros((len(koffrange), len(crange)))
    arr1 = np.zeros((len(koffrange), len(crange)))
    arr2 = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = eqn1(cval, koffval)/eqn2(cval, koffval)
            arr1[i, j] = eqn1(cval, koffval)
            arr2[i, j] = eqn2(cval, koffval)

    plot_heatmap(arr, crange, koffrange, filename, label, log_norm=log_norm, dedim=dedim, **contour_args)
    return 0

def plot_dictionary_ratio(dict_ratio, dedim=False, subdir1='heatmaps', longtime=False, highG=False, contour_args=None):

    if not os.path.exists(DIR_OUTPUT + os.sep + subdir1 + os.sep + dict_ratio['subdir2']):
        os.makedirs(DIR_OUTPUT + os.sep + subdir1 + os.sep + dict_ratio['subdir2'])

    for ratio in dict_ratio['plots']:
        heatmap_ratio(dict_ratio['plots'][ratio]['num'], dict_ratio['plots'][ratio]['denom'], dict_ratio['plots'][ratio]['label'], subdir1 + os.sep + dict_ratio['subdir2'] + os.sep + ratio, dict_ratio['log'], dedim=dedim, contour_args=contour_args)
    return 0

def heatmap_one_equation(equation, label, filename, log_norm, crange=CRANGE, koffrange=KOFFRANGE, dedim=False, kon=KON, T=T, KF=KF, KP=KP, contour_args={'None' : None}):
    print("Plotting equation:",equation)
    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = equation(cval, koffval)+10**(-10)

    if np.min(arr[i,j]) < 0:
        print("Cannot plot log, negative values in equation",equation)
        return 0
    plot_heatmap(arr, crange, koffrange, filename, label, dedim=dedim, log_norm=log_norm, **contour_args)
    return 0

def plot_dictionary_one_equation(dict_eqns, subdir1='heatmaps', longtime=False, highG=False, dedim=False, contour_args={'None' : None}):

    if not os.path.exists(DIR_OUTPUT + os.sep + subdir1 + os.sep + dict_eqns['subdir2']):
        os.makedirs(DIR_OUTPUT + os.sep + subdir1 + os.sep + dict_eqns['subdir2'])

    for name_eqn in dict_eqns['plots']:
        heatmap_one_equation(dict_eqns['plots'][name_eqn]['eqn'], dict_eqns['plots'][name_eqn]['label'], subdir1 + os.sep + dict_eqns['subdir2'] + os.sep + name_eqn, dict_eqns['log'], dedim=dedim, contour_args=contour_args)
    return 0

def dk_plotting():
    heatmap_mode1_error_x(make_heatmap=False, make_panel=True)
    heatmap_mode1_error_x()
    figure_2_combined_cross_sections()

    heatmap_combined_error_c()
    heatmap_combined_error_koff()
    heatmap_kpr_error_c()
    heatmap_kpr_error_koff()

    heatmap_kpr2_error_c()
    heatmap_kpr2_error_koff()
    heatmap_figure_4(crange=np.linspace(0.0001, 5, 80), koffrange=np.linspace(0.0001, 50, 80))
    return 0

if __name__ == '__main__':
    """
    This is how I use the code generally.
    I create a dictionary with a bunch of equations I want to plot in dictionary plttoing (see equations.py).
    These equations are all taken from mathematica or rescaled versions of these (see equations_txt2python for some rescaled equations, trace and eigenvalue )
    There are 2 types of dictionaries: ratiosor 1 equation. Depending on what type of dictionary you create, you will want to use the respoective plotting function above (ether plot_dictionary_one_equation or plot_dictionary_ratio)
    You have to specify some arguments, check equations above should be obvious.
    """
    dictionary = pd.MAIN; want_dedim = True; subdir_2_use = 'heatmaps'
    #contour_args = {'levels' : [0.1, 1., 10.], 'contour_linestyle' : ['dashed','solid','dashed'], 'contour_color' : ['b','w','r'], 'contour_linewidths': [2,2,2]}
    counter_args = {'levels' : [1/(KP*T), 10/(KP*T), 100/(KP*T), 1000/(KP*T), 1E4/(KP*T)]}

    plot_dictionary_one_equation(dictionary, subdir1=subdir_2_use, dedim=True, contour_args=contour_args)

#heatmap_ratios()
