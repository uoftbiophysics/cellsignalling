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

# plot params
fs = 12


def plot_heatmap(arr, crange, koffrange, fname, label):
    # TODO change colour scheme, see https://matplotlib.org/examples/color/colormaps_reference.html
    # TODO better tick labels and gaps?
    # TODO explore alternate option of https://seaborn.pydata.org/generated/seaborn.heatmap.html

    # plot setup
    f = plt.figure()
    imshow_kw = {'cmap': 'YlGnBu', 'aspect': None, 'vmin': np.min(arr), 'vmax': np.max(arr)}
    im = plt.imshow(arr, **imshow_kw)
    # axes setup
    ax = plt.gca()
    # method 1
    period_x_ticks = 33
    period_y_ticks = 1
    ax.set_xticks([i for i, cval in enumerate(crange) if i % period_x_ticks == 0])
    ax.set_yticks([i for i, kval in enumerate(koffrange) if i % period_y_ticks == 0])
    ax.set_xticklabels(['%.3f' % cval for i, cval in enumerate(crange) if i % period_x_ticks == 0], fontsize=fs)
    ax.set_yticklabels(['%.3f' % kval for i, kval in enumerate(koffrange) if i % period_y_ticks == 0], fontsize=fs)
    ax.invert_yaxis()
    ax.set_xlabel(r'$c$')
    ax.set_ylabel(r'$k_{off}$')
    # create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=fs, labelpad=20)
    # save
    plt.savefig(DIR_OUTPUT + os.sep + fname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + fname + '.eps')
    plt.show()
    return


def heatmap_mode1_error_c():

    def mode1_error_c(c, koff):
        x = c * KON / koff
        val = (1 + x) / (KP * T * x) * ((1 + x) ** 2 + 2 * KP / koff)
        return val

    crange = np.logspace(-2, 1, 100)
    #koffrange = np.logspace(-2, 1, 100)
    koffrange = [1.0]
    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = mode1_error_c(cval, koffval)

    label = r'$\langle\delta c^{2}\rangle$/$c^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_mode1_heuristic_error_c', label)
    return


def heatmap_mode1_error_koff():

    def mode1_error_koff(c, koff):
        x = c * KON / koff
        # TODO derive correct error expression, march 2 note
        assert 1 == 2
        #val = (1 + x) / (KP * T * x) * ((1 + x) ** 2 + 2 * KP / koff)
        return val

    crange = np.arange(0.1, 10.0, 0.1)
    koffrange = np.arange(0.1, 10.0, 0.1)
    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = mode1_error_koff(cval, koffval)

    label = r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_mode1_heuristic_error_koff', label)
    return


def heatmap_combined_error_c():

    def combined_error_c(c, koff):
        x = c * KON / koff
        num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff**2 * T * x ** 2 * (1 + x) ** 3
        den = koff ** 2 * KP * T**2 * x * (1 + x) ** 2
        val = num / den
        return val

    crange = np.arange(0.1, 10.0, 0.1)
    koffrange = np.arange(0.1, 10.0, 0.1)
    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = combined_error_c(cval, koffval)

    label = r'$\langle\delta c^{2}\rangle$/$c^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_combined_heuristic_error_c', label)
    return


def heatmap_combined_error_koff():

    def combined_error_koff(c, koff):
        x = c * KON / koff
        num = 2 * KP * x + koff * KP * T * (1 + x) ** 3 + koff**2 * T * (1 + x) ** 3
        den = koff ** 2 * KP * T**2 * x * (1 + x) ** 2
        val = num / den
        return val

    crange = np.arange(0.1, 10.0, 0.1)
    koffrange = np.arange(0.1, 10.0, 0.1)
    arr = np.zeros((len(koffrange), len(crange)))
    for i, koffval in enumerate(koffrange):
        for j, cval in enumerate(crange):
            arr[i, j] = combined_error_koff(cval, koffval)

    label = r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$'
    plot_heatmap(arr, crange, koffrange, 'heatmap_combined_heuristic_error_koff', label)
    return


if __name__ == '__main__':
    heatmap_mode1_error_c()
    #heatmap_mode1_error_koff()  # TODO derive error expression
    #heatmap_combined_error_c()
    #heatmap_combined_error_koff()
