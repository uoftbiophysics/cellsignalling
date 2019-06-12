import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

from load_inputs import DATADICT
from settings import DIR_OUTPUT, DIR_INPUT
from settings import COLOR_SCHEME as cs

plt.style.use('parameters.mplstyle')  # particularIMporting

# notebook params
KON = 1
KP = 10
T = 100
KF = 1.0

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

# global vmax and vmin
assert LOG_START_KOFF == -3
assert LOG_END_KOFF == 3
assert LOG_START_C == -3
assert LOG_END_C == 3
vmax_obs = 1000003020
vmin_obs = 0.00022000016528925625


def plot_heatmap(arr, crange, koffrange, fname, label, show=SHOW):
    # TODO change colour scheme, see https://matplotlib.org/examples/color/colormaps_reference.html
    # TODO fix ticks randomly disappearing on colourbar + flip colourbar minor ticks or remove?
    """
    Colours viridis, YlGnBu, terrain, plasma
    """
    print 'arr limits:', np.min(arr), np.max(arr)
    # plot setup
    f = plt.figure()
    imshow_kw = {'cmap': 'YlGnBu', 'aspect': None, 'vmin': np.min(arr), 'vmax': np.max(arr), 'norm': mpl.colors.LogNorm()}
    #imshow_kw = {'cmap': 'YlGnBu', 'aspect': None, 'vmin': vmin_obs, 'vmax': vmax_obs, 'norm': mpl.colors.LogNorm()}
    im = plt.imshow(arr, **imshow_kw)
    # axes setup
    ax = plt.gca()
    # method 1
    ax.set_xticks([i for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS == 0])
    ax.set_yticks([i for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0])
    #ax.set_xticklabels(['%.3f' % cval for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS == 0], fontsize=FS)
    #ax.set_yticklabels(['%.3f' % kval for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS == 0], fontsize=FS)
    ax.set_xticklabels([r'$10^{%d}$' % np.log10(cval) for i, cval in enumerate(crange) if i % POINTS_BETWEEN_TICKS==0],
                       fontsize=FS)
    ax.set_yticklabels([r'$10^{%d}$' % np.log10(kval) for i, kval in enumerate(koffrange) if i % POINTS_BETWEEN_TICKS==0],
                       fontsize=FS)
    ax.invert_yaxis()
    ax.set_xlabel(r'$c$', fontsize=FS)
    ax.set_ylabel(r'$k_{off}$', fontsize=FS)

    # create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=FS, labelpad=20)
    cbar.ax.tick_params(labelsize=FS)
    # TODO ID why do ticks hide sometimes?
    #for t in cbar.ax.get_yticklabels(): print(t.get_text())

    # contour line for value 1.0
    plt.contour(arr, levels=[1.0], linestyles=['dashed'])  # use 'dashed' or 'solid' curve

    # save
    plt.savefig(DIR_OUTPUT + os.sep + fname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + fname + '.eps')
    if show:
        plt.show()
    return


def heatmap_mode1_error_x(crange=CRANGE, koffrange=KOFFRANGE):

    def mode1_error_x(c, koff):
        x = c * KON / koff
        val = (1 + x) / (KP * T * x) * ((1 + x) ** 2 + 2 * KP / koff)
        return val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = mode1_error_x(cval, koffval)

    label = r'$\langle\delta x^{2}\rangle$/$x^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_mode1_heuristic_error_x', label)
    return


def heatmap_combined_error_c(crange=CRANGE, koffrange=KOFFRANGE):

    def combined_error_c(c, koff):
        x = c * KON / koff
        num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff**2 * T * x ** 2 * (1 + x) ** 3
        den = koff ** 2 * KP * T**2 * x * (1 + x) ** 2
        val = num / den
        return val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = combined_error_c(cval, koffval)

    label = r'$\langle\delta c^{2}\rangle$/$c^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_combined_heuristic_error_c', label)
    return


def heatmap_combined_error_koff(crange=CRANGE, koffrange=KOFFRANGE):

    def combined_error_koff(c, koff):
        x = c * KON / koff
        num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff**2 * T * (1 + x) ** 3
        den = koff ** 2 * KP * T**2 * x * (1 + x) ** 2
        val = num / den
        return val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = combined_error_koff(cval, koffval)

    label = r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_combined_heuristic_error_koff', label)
    return


def heatmap_kpr_error_c(crange=CRANGE, koffrange=KOFFRANGE):

    def kpr_error_c_full(c, koff):
        x = c * KON / koff
        g = koff / KF
        factor_A = (1 + x) / ((1 + g) * koff**2 * KP * T**2 * x)
        term1 = 2 * g * (KP + koff * KP * T - koff**2 * T * x)
        term2 = koff * T * (KP + koff * x**2)
        term3 = g**2 * (koff**2 * T + KP * (2 * (1 + x) + koff * T * (2 + x * (2 + x))))
        val = factor_A * (term1 + term2 + term3)
        return val

    def kpr_error_c_high_g_t(c, koff):
        x = c * KON / koff
        left = (1 + koff / KF) * (1 + x) / (koff * T * x)
        right = koff / KP + 2 + 2 * x + x ** 2
        val = left * right
        return val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = kpr_error_c_full(cval, koffval)

    label = r'$\langle\delta c^{2}\rangle$/$c^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_kpr_heuristic_error_c', label)
    return


def heatmap_kpr_error_koff(crange=CRANGE, koffrange=KOFFRANGE):

    def kpr_error_koff(c, koff):
        # Note the "full" heuristic expression for koff error is the same as the high g high t one
        x = c * KON / koff
        left = (1 + koff / KF) * (1 + x) / (koff * T * x)
        right = koff / KP + 1
        val = left * right
        return val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = kpr_error_koff(cval, koffval)

    label = r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_kpr_heuristic_error_koff', label)
    return


def heatmap_kpr2_error_c(crange=CRANGE, koffrange=KOFFRANGE):

    def kpr2_error_c_full(c, koff):
        kon = KON
        kf = KF
        kp = KP
        t = T
        g = koff / KF
        x = c * kon / koff

        val = ((kf + koff + kf * x + koff * x) ** 4 * ((-2 * kf * koff ** 2 * kp ** 2 * t ** 2 * x ** 2 * (
                    1 - x + g ** 2 * (1 + x) ** 2 - g * (-2 + x + x ** 2)) * (-(
                    (koff * kp ** 2 * t ** 2 * x) / (kf + koff + kf * x + koff * x)) + ((kf * kp * t * x) / (
                    kf + koff + kf * x + koff * x) + (koff * kp * t * x) / (kf + koff + kf * x + koff * x)) ** 2) * (((
                                                                                                                                  kf * kp * t * x) / (
                                                                                                                                  kf + koff + kf * x + koff * x) + (
                                                                                                                                  koff * kp * t * x) / (
                                                                                                                                  kf + koff + kf * x + koff * x)) ** 2 - kp * t * (
                                                                                                                                 (
                                                                                                                                             kf * kp * t * x) / (
                                                                                                                                             kf + koff + kf * x + koff * x) + (
                                                                                                                                             2 * koff * kp * t * x) / (
                                                                                                                                             kf + koff + kf * x + koff * x)))) / (
                                                             kf + koff + kf * x + koff * x) ** 2 + (
                                                             g * kf ** 2 * kp ** 2 * t ** 2 * x ** 2 * (
                                                                 (1 + g) ** 2 * koff * (1 + x) ** 2 + 2 * g * kp * (
                                                                     1 + g + x + x ** 2)) * (((kf * kp * t * x) / (
                                                                 kf + koff + kf * x + koff * x) + (
                                                                                                          koff * kp * t * x) / (
                                                                                                          kf + koff + kf * x + koff * x)) ** 2 - kp * t * (
                                                                                                         (
                                                                                                                     kf * kp * t * x) / (
                                                                                                                     kf + koff + kf * x + koff * x) + (
                                                                                                                     2 * koff * kp * t * x) / (
                                                                                                                     kf + koff + kf * x + koff * x))) ** 2) / (
                                                             kf + koff + kf * x + koff * x) ** 2 + (
                                                             koff ** 2 * kp ** 2 * t ** 2 * x ** 2 * (-(
                                                                 (koff * kp ** 2 * t ** 2 * x) / (
                                                                     kf + koff + kf * x + koff * x)) + ((
                                                                                                                    kf * kp * t * x) / (
                                                                                                                    kf + koff + kf * x + koff * x) + (
                                                                                                                    koff * kp * t * x) / (
                                                                                                                    kf + koff + kf * x + koff * x)) ** 2) ** 2 * (
                                                                         (1 + g) ** 2 * koff * (1 + x) ** 2 + 2 * kp * (
                                                                             1 + g * (2 + x + g * (1 + x) ** 2)))) / (
                                                             kf + koff + kf * x + koff * x) ** 2)) / (
                    c ** 2 * (1 + g) ** 3 * kf ** 2 * koff * kon ** 2 * kp ** 3 * t ** 3 * x ** 3 * (1 + x) ** 3 * (
                        -(kp * t) + (kf * kp * t * x) / (kf + koff + kf * x + koff * x) + (koff * kp * t * x) / (
                            kf + koff + kf * x + koff * x)) ** 4)
        return val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = kpr2_error_c_full(cval, koffval)

    label = r'$\langle\delta c^{2}\rangle$/$c^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_kpr2_heuristic_error_c', label)
    return


def heatmap_kpr2_error_koff(crange=CRANGE, koffrange=KOFFRANGE):

    def kpr2_error_koff(c, koff):
        # Note the "full" heuristic expression for koff error is the same as the high g high t one
        kon = KON
        kf = KF
        kp = KP
        t = T
        g = koff / KF
        x = c * kon / koff

        val = ((kf + koff + kf * x + koff * x) ** 4 * (koff * (
                    (g * (1 + g) ** 2 * kf ** 2 * kp ** 2 * t ** 2 * x ** 2 * (1 + x) ** 2) / (
                        kf + koff + kf * x + koff * x) ** 2 + (
                                (1 + g) ** 2 * koff ** 2 * kp ** 2 * t ** 2 * x ** 2 * (1 + x) ** 2) / (
                                kf + koff + kf * x + koff * x) ** 2 - (2 * kf * koff * kp ** 2 * t ** 2 * x ** 2 * (
                        (1 + g) ** 2 + (-1 + g) * (1 + 2 * g) * x + (-1 + g) * g * x ** 2)) / (
                                kf + koff + kf * x + koff * x) ** 2) + 2 * kp * ((
                                                                                             g ** 2 * kf ** 2 * kp ** 2 * t ** 2 * x ** 2 * (
                                                                                                 1 + g + x + x ** 2)) / (
                                                                                             kf + koff + kf * x + koff * x) ** 2 + (
                                                                                             koff ** 2 * kp ** 2 * t ** 2 * x ** 2 * (
                                                                                                 1 + g * (2 + x + g * (
                                                                                                     1 + x) ** 2))) / (
                                                                                             kf + koff + kf * x + koff * x) ** 2))) / (
                    (1 + g) ** 3 * kf ** 2 * koff ** 3 * kp ** 3 * t ** 3 * x ** 3 * (1 + x) ** 3)

        return val

    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = kpr2_error_koff(cval, koffval)

    label = r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_kpr2_heuristic_error_koff', label)
    return


if __name__ == '__main__':
    """
    heatmap_mode1_error_x()
    heatmap_combined_error_c()
    heatmap_combined_error_koff()
    heatmap_kpr_error_c()
    heatmap_kpr_error_koff()
    """
    heatmap_kpr2_error_c()
    heatmap_kpr2_error_koff()
