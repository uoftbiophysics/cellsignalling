from settings import *
import numpy as np
import matplotlib.pyplot as plt
from load_inputs import DATADICT
import os

def make_figure_3():
    """
    Panel A: rel_err_c
    """
    return 0

def make_figure_2():
    """
    Panel A: <n> vs c for two different Kd, showing overlapping "clouds"
    Panel B: rel_err_x vs x for heuristic form of relative error
    """
    fig, axes = plt.subplots(nrows=1, ncols=2)
    fig.set_size_inches(10, 5)
    # Panel A
    axes[0].set_xlabel(DATADICT["mode1_composite_mean_N_Black"]['xlab'])
    axes[0].set_ylabel(DATADICT["mode1_composite_mean_N_Black"]['ylab'])
    #panelA_xticks = ["%.2f" % el for el in DATADICT["mode1_composite_mean_N_Black"]['xpts'][0::int(len(DATADICT["mode1_composite_mean_N_Black"]['xpts'])/8)]]
    #panelA_yticks = [0.2 * r for r in range(6)]
    #axes[0].xticks([float(el) for el in panelA_xticks], panelA_xticks)
    #axes[0].yticks(panelA_yticks, ["%.2f" % el for el in panelA_yticks])

    axes[0].plot(DATADICT["mode1_composite_mean_N_Black"]['xpts'], DATADICT["mode1_composite_mean_N_Black"]['ypts'], 'k')
    axes[0].plot(DATADICT["mode1_composite_Err_N_Black_Upper"]['xpts'], DATADICT["mode1_composite_Err_N_Black_Upper"]['ypts'], 'k--')
    axes[0].plot(DATADICT["mode1_composite_Err_N_Black_Lower"]['xpts'], DATADICT["mode1_composite_Err_N_Black_Lower"]['ypts'], 'k--')
    axes[0].plot(DATADICT["mode1_composite_mean_N_Red"]['xpts'], DATADICT["mode1_composite_mean_N_Red"]['ypts'], 'r')
    axes[0].plot(DATADICT["mode1_composite_Err_N_Red_Upper"]['xpts'], DATADICT["mode1_composite_Err_N_Red_Upper"]['ypts'], 'r--')
    axes[0].plot(DATADICT["mode1_composite_Err_N_Red_Lower"]['xpts'], DATADICT["mode1_composite_Err_N_Red_Lower"]['ypts'], 'r--')

    # Panel B
    axes[1].plot(DATADICT["mode1_error_compare_heuristic"]['xpts'], DATADICT["mode1_error_compare_heuristic"]['ypts'], 'k')
    axes[1].set_xlabel('x')
    axes[1].set_ylabel(r'$\delta x^{2}$/$x^{2}$')
    axes[1].set_ylim([0, 1.5])

    # Save figure
    plt.savefig(os.path.join(DIR_OUTPUT, 'figure2.pdf'))
    plt.savefig(os.path.join(DIR_OUTPUT, 'figure2.png'))


if __name__ == "__main__":
    make_figure_2()

