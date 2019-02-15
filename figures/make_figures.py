import matplotlib.pyplot as plt
import numpy as np
import os

from load_inputs import DATADICT
from settings import DIR_OUTPUT


def make_figure_2():
    """
    Signal specificity and relative error in estimation of x for Mode 1.
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
    plt.savefig(os.path.join(DIR_OUTPUT, 'figure2.eps'))


def make_figure_3():
    """
    Combined model heuristic relative errors in c and koff.
    One panel.
    Left axis: rel_err_c vs x
    Right axis: rel_err_koff vs x
    """
    # Set up figure
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    # Plot
    ax1.plot(DATADICT["combined_error_composite_compare_c_heuristic"]['xpts'],
              DATADICT["combined_error_composite_compare_c_heuristic"]['ypts'], 'k')
    ax2.plot(DATADICT["combined_error_composite_compare_koff_heuristic"]['xpts'],
              DATADICT["combined_error_composite_compare_koff_heuristic"]['ypts'], 'r')
    ax1.set_xlabel(DATADICT["combined_error_composite_compare_c_heuristic"]['xlab'])
    ax1.set_ylabel(DATADICT["combined_error_composite_compare_c_heuristic"]['ylab'])
    ax2.set_ylabel(DATADICT["combined_error_composite_compare_koff_heuristic"]['ylab'])
    # Save figure
    plt.savefig(os.path.join(DIR_OUTPUT, 'figure3.pdf'))
    plt.savefig(os.path.join(DIR_OUTPUT, 'figure3.eps'))


def make_figure_5():
    """
    KPR model heuristic form of relative errors for c and koff.
    One panel.
    Left axis: rel_err_c vs x
    Right axis: rel_err_koff vs x
    """
    # Set up figure
    fig, ax1 = plt.subplots()
    # For two y axes use:
    ax2 = ax1.twinx()

    # Plot
    ax1.plot(DATADICT["KPR_error_composite_compare_c_heuristic"]['xpts'],
              DATADICT["KPR_error_composite_compare_c_heuristic"]['ypts'], 'k')
    ax2.plot(DATADICT["KPR_error_composite_compare_koff_heuristic"]['xpts'],
              DATADICT["KPR_error_composite_compare_koff_heuristic"]['ypts'], 'r')
    ax1.set_xlabel(DATADICT["KPR_error_composite_compare_c_heuristic"]['xlab'])
    ax1.set_ylabel(DATADICT["KPR_error_composite_compare_c_heuristic"]['ylab'])
    ax2.set_ylabel(DATADICT["KPR_error_composite_compare_koff_heuristic"]['ylab'])

    # Save figure
    plt.savefig(os.path.join(DIR_OUTPUT, 'figure5.pdf'))
    plt.savefig(os.path.join(DIR_OUTPUT, 'figure5.eps'))


def make_figure_B1():
    """
    """
    return 0


def make_figure_B2():
    """
    """
    return 0


def make_figure_B3():
    """
    """
    return 0


def make_figure_C1():
    """
    Mode 1 relative error in x as estimated by two different methods.
    One panel.
    """
    # Set up axes
    fig, axes = plt.subplots()
    axes.set_xlabel('x')
    axes.set_ylabel(r'$\delta x^{2}$/$x^{2}$')
    axes.set_ylim([0, 1.5])

    # Plot
    axes.plot(DATADICT["mode1_error_compare_heuristic"]['xpts'], DATADICT["mode1_error_compare_heuristic"]['ypts'], 'k')
    axes.scatter(DATADICT["mode1_error_compare_fullFisher"]['xpts'], DATADICT["mode1_error_compare_fullFisher"]['ypts'], c='b')

    # Save figure
    plt.savefig(os.path.join(DIR_OUTPUT, 'figureC1.pdf'))
    plt.savefig(os.path.join(DIR_OUTPUT, 'figureC1.eps'))



def make_figure_C2():
    """
    Combined model relative error in c and koff, as estimated by three different methods.
    Panel A: relative error in c as estimated by heuristic, saddle point Fisher, and full Fisher
    Panel B: relative error in koff as estimated by heuristic, saddle point Fisher, and full Fisher
    """
    # Set up axes
    fig, axes = plt.subplots(nrows=1, ncols=2)
    axes[0].set_xlabel(DATADICT["combined_error_composite_compare_c_heuristic"]['xlab'])
    axes[1].set_xlabel(DATADICT["combined_error_composite_compare_koff_heuristic"]['xlab'])
    axes[0].set_ylabel(r'$\delta c^{2}$/$c^{2}$')
    axes[1].set_ylabel(r'$\delta koff^{2}$/$koff^{2}$')
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')
    # Plot
    # Panel A:
    axes[0].plot(DATADICT["combined_error_composite_compare_c_heuristic"]['xpts'],
                 DATADICT["combined_error_composite_compare_c_heuristic"]['ypts'], 'k')
    axes[0].scatter(DATADICT["combined_error_composite_compare_c_saddlePointFisher"]['xpts'],
                    DATADICT["combined_error_composite_compare_c_saddlePointFisher"]['ypts'], c='b')
    axes[0].scatter(DATADICT["combined_error_composite_compare_c_fullFisher"]['xpts'],
                    DATADICT["combined_error_composite_compare_c_fullFisher"]['ypts'], c='r')

    # Panel B:
    axes[1].plot(DATADICT["combined_error_composite_compare_koff_heuristic"]['xpts'],
                 DATADICT["combined_error_composite_compare_koff_heuristic"]['ypts'], 'k')
    axes[1].scatter(DATADICT["combined_error_composite_compare_koff_saddlePointFisher"]['xpts'],
                    DATADICT["combined_error_composite_compare_koff_saddlePointFisher"]['ypts'], c='b')
    axes[1].scatter(DATADICT["combined_error_composite_compare_koff_fullFisher"]['xpts'],
                    DATADICT["combined_error_composite_compare_koff_fullFisher"]['ypts'], c='r')

    # Save figure
    plt.savefig(os.path.join(DIR_OUTPUT, 'figureC2.pdf'))
    plt.savefig(os.path.join(DIR_OUTPUT, 'figureC2.eps'))



def make_figure_C3():
    """
    KPR model relative error in c and koff, as estimated by two different methods.
    Panel A: relative error in c as estimated by heuristic, and saddle point Fisher
    Panel B: relative error in koff as estimated by heuristic, and saddle point Fisher
    """
    # Set up axes
    fig, axes = plt.subplots(nrows=1, ncols=2)
    axes[0].set_xlabel(DATADICT["KPR_error_composite_compare_c_heuristic"]['xlab'])
    axes[1].set_xlabel(DATADICT["KPR_error_composite_compare_koff_heuristic"]['xlab'])
    axes[0].set_ylabel(r'$\delta c^{2}$/$c^{2}$')
    axes[1].set_ylabel(r'$\delta koff^{2}$/$koff^{2}$')
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')
    # Plot
    # Panel A:
    axes[0].plot(DATADICT["KPR_error_composite_compare_c_heuristic"]['xpts'],
                 DATADICT["KPR_error_composite_compare_c_heuristic"]['ypts'], 'k')
    axes[0].scatter(DATADICT["KPR_error_composite_compare_c_saddlePointFisher"]['xpts'],
                    DATADICT["KPR_error_composite_compare_c_saddlePointFisher"]['ypts'], c='b')

    # Panel B:
    axes[1].plot(DATADICT["KPR_error_composite_compare_koff_heuristic"]['xpts'],
                 DATADICT["KPR_error_composite_compare_koff_heuristic"]['ypts'], 'k')
    axes[1].scatter(DATADICT["KPR_error_composite_compare_koff_saddlePointFisher"]['xpts'],
                    DATADICT["KPR_error_composite_compare_koff_saddlePointFisher"]['ypts'], c='b')

    # Save figure
    plt.savefig(os.path.join(DIR_OUTPUT, 'figureC3.pdf'))
    plt.savefig(os.path.join(DIR_OUTPUT, 'figureC3.eps'))


def make_figure_D1():
    """
    Mode 1 with a prior, plotting MLE for x vs x and comparing to heuristic (ie. non-prior prediction)
    """
    return 0


if __name__ == "__main__":
    make_figure_2()
    make_figure_3()
    make_figure_5()
    make_figure_B1()
    make_figure_B2()
    make_figure_B3()
    make_figure_C1()
    make_figure_C2()
    make_figure_C3()
    make_figure_D1()
