import matplotlib.pyplot as plt
import numpy as np
import os

from load_inputs import DATADICT
from settings import DIR_OUTPUT, DIR_INPUT
from settings import COLOR_SCHEME as cs

plt.style.use('parameters.mplstyle')  # particularIMporting

GLOB_KON = 10
GLOB_KP = 10
GLOB_T = 100


def heatmap_mode1_error_c():

    def mode1_error_c(c, koff):
        return 1

    crange = np.arange(0.1, 10.0, 0.1)
    koffrange = np.arange(0.1, 10.0, 0.1)
    arr = np.zeros((len(crange), len(koffrange)))
    for i, cval in enumerate(crange):
        for j, koffval in enumerate(koffrange):
            arr[i, j] = mode1_error_c(cval, koffval)

    # plot
    plt.imshow(arr)
    plt.show()
    return


def heatmap_mode1_error_koff():

    def mode1_error_koff(c, koff):
        return 1

    crange = np.arange(0.1, 10.0, 0.1)
    koffrange = np.arange(0.1, 10.0, 0.1)
    arr = np.zeros((len(crange), len(koffrange)))
    for i, cval in enumerate(crange):
        for j, koffval in enumerate(koffrange):
            arr[i, j] = mode1_error_koff(cval, koffval)

    # plot
    plt.imshow(arr)
    plt.show()
    return


def heatmap_combined_error_c():

    def combined_error_c(c, koff):
        return 1

    crange = np.arange(0.1, 10.0, 0.1)
    koffrange = np.arange(0.1, 10.0, 0.1)
    arr = np.zeros((len(crange), len(koffrange)))
    for i, cval in enumerate(crange):
        for j, koffval in enumerate(koffrange):
            arr[i, j] = combined_error_c(cval, koffval)

    # plot
    plt.imshow(arr)
    plt.show()
    return


def heatmap_combined_error_koff():

    def combined_error_koff(c, koff):
        return 1

    crange = np.arange(0.1, 10.0, 0.1)
    koffrange = np.arange(0.1, 10.0, 0.1)
    arr = np.zeros((len(crange), len(koffrange)))
    for i, cval in enumerate(crange):
        for j, koffval in enumerate(koffrange):
            arr[i, j] = combined_error_koff(cval, koffval)

    # plot
    plt.imshow(arr)
    plt.show()
    return


if __name__ == '__main__':
    heatmap_mode1_error_c()
    heatmap_mode1_error_koff()
    heatmap_combined_error_c()
    heatmap_combined_error_koff()
