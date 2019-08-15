import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

#from load_inputs import DATADICT
from settings import DIR_OUTPUT, DIR_INPUT
from settings import COLOR_SCHEME as cs

plt.style.use('parameters.mplstyle')  # particularIMporting

# unit params
KON = 1
KP = 10
T = 1000
KF = 1.0

# de-dimensionalise params
c0 = KP/KON
alpha = KP*T

# plot params
FS = 16
SHOW = False

# axes
POINTS_BETWEEN_TICKS = 10
LOG_START_C = -3
LOG_END_C = 3
TOTAL_POINTS_C = (LOG_END_C - LOG_START_C) * POINTS_BETWEEN_TICKS + 1
CRANGE = np.logspace(LOG_START_C, LOG_END_C, TOTAL_POINTS_C)
LOG_START_KOFF = -3
LOG_END_KOFF = 3
TOTAL_POINTS_KOFF = (LOG_END_KOFF - LOG_START_KOFF) * POINTS_BETWEEN_TICKS + 1
KOFFRANGE = np.logspace(LOG_START_KOFF, LOG_END_KOFF, TOTAL_POINTS_KOFF)

CTILDERANGE = np.divide(list(CRANGE), c0)
ZRANGE = np.divide(list(KOFFRANGE), KP)

# global vmax and vmin
assert LOG_START_KOFF == -3
assert LOG_END_KOFF == 3
assert LOG_START_C == -3
assert LOG_END_C == 3
vmax_obs = 1000003020
vmin_obs = 0.00022000016528925625


def plot_heatmap(arr, crange, koffrange, fname, label, show=SHOW, save=True, log_norm=True, **kwargs):
    # default parameters
    if 'levels' in kwargs.keys():
        levels = kwargs['levels']
    else:
        levels = [10, 100, 1000, 1E4]
    if 'vmin' in kwargs.keys():
        vmin = kwargs['vmin']
    else:
        vmin = np.min(arr)
    if 'vmax' in kwargs.keys():
        vmax = kwargs['vmax']
    else:
        vmax = np.max(arr)
    if 'contour_linestyle' in kwargs.keys():
        contour_linestyle = 'dashed'
    else:
        contour_linestyle = 'solid'
    if 'contour_color' in kwargs.keys():
        contour_color = kwargs['contour_color']
    else:
        contour_color = None
    if 'contour_linewidths' in kwargs.keys():
        contour_lindewidths = kwargs['contour_linewidths']
    else:
        contour_lindewidths = None
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
    else:
        imshow_kw = {'cmap': 'YlGnBu', 'aspect': None, 'vmin': vmin, 'vmax': vmax}
    #imshow_kw = {'cmap': 'YlGnBu', 'aspect': None, 'vmin': vmin_obs, 'vmax': vmax_obs, 'norm': mpl.colors.LogNorm()}
    im = plt.imshow(arr, **imshow_kw)

    # axes setup
    fig = plt.gcf()
    ax = plt.gca()

    # method 1
    ax.set_xticks([i for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS == 0])
    ax.set_yticks([i for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0])
    #ax.set_xticklabels(['%.3f' % cval for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS == 0], fontsize=FS)
    #ax.set_yticklabels(['%.3f' % kval for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0], fontsize=FS)
    if log_norm:
        ax.set_xticklabels([r'$10^{%d}$' % np.log10(cval) for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS==0],
                           fontsize=FS)
        ax.set_yticklabels([r'$10^{%d}$' % np.log10(kval) for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS==0],
                           fontsize=FS)
    else:
        ax.set_xticklabels(["{:10.2f}".format(cval) for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS==0],
                           fontsize=FS)
        ax.set_yticklabels(["{:10.2f}".format(kval) for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS==0],
                           fontsize=FS)

    ax.invert_yaxis()
    ax.set_xlabel(r'$\tilde{c}$', fontsize=FS)
    ax.set_ylabel('z', fontsize=FS)

    # create colorbar
    cbar = fig.colorbar(im)
    cbar.ax.minorticks_off()
    cbar.update_ticks()
    cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=FS, labelpad=20)
    cbar.ax.tick_params(labelsize=FS)
    # TODO ID why do ticks hide sometimes?
    #for t in cbar.ax.get_yticklabels(): print(t.get_text())

    # contour line for value 1.0
    plt.contour(arr, levels=levels, linestyles=[contour_linestyle], colors=contour_color, linewidths=contour_lindewidths)

    # save
    if save == True:
        plt.savefig(DIR_OUTPUT + os.sep + fname + '.pdf')
        plt.savefig(DIR_OUTPUT + os.sep + fname + '.eps')
    if show:
        plt.show()

    return fig, ax


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


def heatmap_combined_error_koff(crange=CTILDERANGE, koffrange=ZRANGE,
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

T=100
def heatmap_figure_4(crange=CRANGE, koffrange=KOFFRANGE):
    nobs = 0.1 * KP * T
    mobs = 0.15 * KP * T
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
    koff_star = [ax.get_ylim()[1]/max(koffrange) * heuristic_estimate_koff(nobs, mobs) for _ in range(len(xpts))]
    c_star = [ax.get_xlim()[1] / max(crange) *heuristic_estimate_c(nobs, mobs) for _ in range(len(ypts))]
    ax.plot(xpts, koff_star, 'k--')
    ax.plot(c_star, ypts, 'k--')
    # save figure
    fig.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    fig.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    return

if __name__ == '__main__':

    #heatmap_mode1_error_x(make_heatmap=False, make_panel=True)
    #heatmap_mode1_error_x()
    #figure_2_combined_cross_sections()

    #heatmap_combined_error_c()
    #heatmap_combined_error_koff()
    #heatmap_kpr_error_c()
    #heatmap_kpr_error_koff()

    #heatmap_kpr2_error_c()
    #heatmap_kpr2_error_koff()
    heatmap_figure_4(crange=np.linspace(0.0001, 5, 80), koffrange=np.linspace(0.0001, 50, 80))
