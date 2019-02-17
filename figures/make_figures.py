import matplotlib.pyplot as plt
import numpy as np
import os

from load_inputs import DATADICT
from settings import DIR_OUTPUT

plt.style.use('parameters.mplstyle') # particularIMporting

def make_figure_2():
    """
    Signal specificity and relative error in estimation of x for Mode 1.
    Panel A: <n> vs c for two different Kd, showing overlapping "clouds"
    Panel B: rel_err_x vs x for heuristic form of relative error
    """
    figname = 'mode1_composite'
    curveLk = DATADICT[figname + '_mean_N_Black']
    curveLkUp = DATADICT[figname + '_Err_N_Black_Upper']
    curveLkLow = DATADICT[figname + '_Err_N_Black_Lower']
    curveLr = DATADICT[figname + '_mean_N_Red']
    curveLrUp = DATADICT[figname + '_Err_N_Red_Upper']
    curveLrLow = DATADICT[figname + '_Err_N_Red_Lower']
    curveR = DATADICT['mode1_error_compare_heuristic']
    # plot
    fig, axarr = plt.subplots(nrows=1, ncols=2)
    #fig.set_size_inches(10, 5)
    plt.suptitle(r'Mode 1: Mean $n$ and relative error in $x$ estimate')
    # left plot
    axarr[0].plot(curveLk['xpts'], curveLk['ypts'], 'k', label=r'$k_{off}=10$')
    axarr[0].plot(curveLkUp['xpts'], curveLkUp['ypts'], 'k--')
    axarr[0].plot(curveLkLow['xpts'], curveLkLow['ypts'], 'k--')
    axarr[0].plot(curveLr['xpts'], curveLr['ypts'], 'r', label=r'$k_{off}=15$')
    axarr[0].plot(curveLrUp['xpts'], curveLrUp['ypts'], 'r--')
    axarr[0].plot(curveLrLow['xpts'], curveLrLow['ypts'], 'r--')
    axarr[0].set_xlabel(r'$c$')
    axarr[0].set_ylabel(r'$\langle n\rangle/k_pt$')
    axarr[0].legend()
    axarr[0].set_title(r'($k_p=10$, $t=120$)')
    # right plot
    axarr[1].plot(curveR['xpts'], curveR['ypts'], 'k')
    axarr[1].set_xlabel(r'$x$')
    axarr[1].set_ylabel(r'$\langle\delta x^{2}\rangle$/$x^{2}$')
    axarr[1].set_ylim([0, 1.5])
    axarr[1].set_title('($k_p=10$, $t=100$, $k_{off}=1$)')
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


def make_figure_3():
    """
    Combined model heuristic relative errors in c and koff.
    One panel.
    Left axis: rel_err_c vs x
    Right axis: rel_err_koff vs x
    """
    figname = 'combined_error_composite_compare'
    curve1 = DATADICT[figname + '_c_heuristic']
    curve2 = DATADICT[figname + '_koff_heuristic']
    # plot
    #plt.figure(figsize=(10, 5))
    plt.figure()
    plt.plot(curve1['xpts'], curve1['ypts'], 'k', label=r'$\langle\delta c^{2}\rangle$/$c^{2}$')
    plt.plot(curve2['xpts'], curve2['ypts'], 'r', label=r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$')
    plt.xlabel(r'$c$')
    plt.ylabel('Relative error')
    plt.title('Combined: MLE Relative error for $c$ and $k_{off}$ ($k_p=10$, $t=100$, $k_{off}=1$)')
    plt.legend()
    # set limits
    plt.ylim(0, 1.3)
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


def make_figure_5():
    """
    KPR model heuristic form of relative errors for c and koff.
    One panel.
    Left axis: rel_err_c vs x
    Right axis: rel_err_koff vs x
    """
    figname = 'KPR_error_composite_compare'
    curve1 = DATADICT[figname + '_c_heuristic']
    curve2 = DATADICT[figname + '_koff_heuristic']
    # plot
    #plt.figure(figsize=(10, 5))
    plt.figure()
    plt.plot(curve1['xpts'], curve1['ypts'], 'k', label=r'$\langle\delta c^{2}\rangle$/$c^{2}$')
    plt.plot(curve2['xpts'], curve2['ypts'], 'r', label=r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$')
    plt.xlabel(r'$c$')
    plt.ylabel('Relative error')
    plt.title('KPR: MLE Relative error for $c$ and $k_{off}$ ($k_p=10$, $t=100$, $k_{off}=1$, $k_f=100$)')
    plt.legend()
    # set limits
    #plt.ylim(0, 0.3)
    #plt.xlim(0, 10)
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


def make_figure_B1():
    """
    'mode1_MLE_compare' figure, 2 curves (numeric vs heuristic)
    """
    figname = 'mode1_MLE_compare'
    curve1 = DATADICT[figname + '_numeric']
    curve2 = DATADICT[figname + '_heuristic']
    # plot
    #plt.figure(figsize=(10, 5))
    plt.figure()
    c2_part1 = plt.plot(curve2['xpts'][0:400], curve2['ypts'][0:400], 'r', label='heuristic')
    c2_part2 = plt.plot(curve2['xpts'][400:], curve2['ypts'][400:], 'r')
    c1 = plt.plot(curve1['xpts'], curve1['ypts'], 'ob', label='numeric')
    plt.title('Mode 1 MLE: Numeric vs Heuristic ($k_p=10$, $t=100$, $k_{off}=1$)')
    plt.xlabel(r'$n_{obs}$')
    plt.ylabel(r'$x_{MLE}$')
    plt.legend()
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


def make_figure_B2():
    """
    'combined_MLE_composite_compare' figure, 4 curves (numeric vs heuristic) for c and koff
    """
    figname = 'combined_MLE_composite_compare'
    curveL1 = DATADICT[figname + '_c_numeric']
    curveL2 = DATADICT[figname + '_c_heuristic']
    curveR1 = DATADICT[figname + '_koff_numeric']
    curveR2 = DATADICT[figname + '_koff_heuristic']
    # plot
    fig, axarr = plt.subplots(nrows=1, ncols=2)
    #fig.set_size_inches(10, 5)
    plt.suptitle('Combined MLE: Numeric vs Heuristic ($k_p=10$, $t=100$, $k_{off}=1$, $m=100$)')
    # left subplot
    cL2 = axarr[0].plot(curveL2['xpts'], curveL2['ypts'], 'r', label='heuristic')
    cL1 = axarr[0].plot(curveL1['xpts'], curveL1['ypts'], 'ob', label='numeric')
    axarr[0].set_title(r'Estimate for $c$')
    axarr[0].set_xlabel(r'$n_{obs}$')
    axarr[0].set_ylabel(r'$c^*$')
    axarr[0].legend()
    # right subplot
    cR1 = axarr[1].plot(curveR1['xpts'], curveR1['ypts'], 'ob', label='numeric')
    cR2 = axarr[1].plot(curveR2['xpts'], curveR2['ypts'], 'r', label='heuristic')
    axarr[1].set_title(r'Estimate for $k_{off}$')
    axarr[1].set_xlabel(r'$n_{obs}$')
    axarr[1].set_ylabel(r'$k_{off}^*$')
    axarr[1].legend()
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


def make_figure_B3():
    """
    'KPR_MLE_composite_compare' figure, 4 curves (numeric vs heuristic) for c and koff
    """
    figname = 'KPR_MLE_composite_compare'
    curveL1 = DATADICT[figname + '_c_numeric']
    curveL2 = DATADICT[figname + '_c_heuristic']
    curveR1 = DATADICT[figname + '_koff_numeric']
    curveR2 = DATADICT[figname + '_koff_heuristic']
    # plot
    fig, axarr = plt.subplots(nrows=1, ncols=2)
    #fig.set_size_inches(10, 5)
    plt.figure()
    plt.suptitle('KPR MLE: Numeric vs Heuristic ($k_p=10$, $t=100$, $k_{off}=1$, $k_f=100$, $m=100$)')
    # left subplot
    cL1 = axarr[0].plot(curveL1['xpts'], curveL1['ypts'], 'ob', label='numeric')
    cL2 = axarr[0].plot(curveL2['xpts'], curveL2['ypts'], 'r', label='heuristic')
    axarr[0].set_title(r'Estimate for $c$')
    axarr[0].set_xlabel(r'$n_{obs}$')
    axarr[0].set_ylabel(r'$c^*$')
    axarr[0].legend()
    # right subplot
    cR1 = axarr[1].plot(curveR1['xpts'], curveR1['ypts'], 'ob', label='numeric')
    cR2 = axarr[1].plot(curveR2['xpts'], curveR2['ypts'], 'r', label='heuristic')
    axarr[1].set_title(r'Estimate for $k_{off}$')
    axarr[1].set_xlabel(r'$n_{obs}$')
    axarr[1].set_ylabel(r'$k_{off}^*$')
    axarr[1].legend()
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


def make_figure_C1():
    """
    Mode 1 relative error in x as estimated by two different methods.
    One panel.
    """
    figname = 'mode1_error_compare'
    curve1 = DATADICT[figname + '_heuristic']
    curve2 = DATADICT[figname + '_saddlePointFisher']
    curve3 = DATADICT[figname + '_fullFisher']
    # plot
    #plt.figure(figsize=(10, 5))
    plt.figure()
    plt.plot(curve1['xpts'], curve1['ypts'], 'k', label='Simple Fisher', zorder=1)
    plt.scatter(curve2['xpts'], curve2['ypts'], c='b', edgecolor='', label='Numeric Fisher (Saddle Point)', zorder=2)
    plt.scatter(curve3['xpts'], curve3['ypts'], c='r', edgecolor='', label='Numeric Fisher (Full)', zorder=3)
    # axis
    plt.title('Mode 1: MLE Relative error comparison ($k_p=10$, $t=100$, $k_{off}=1$)')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\langle\delta x^{2}\rangle$/$x^{2}$')
    plt.gca().set_xscale('log')
    plt.ylim([0, 1.5])
    plt.legend()
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


def make_figure_C2():
    """
    Combined model relative error in c and koff, as estimated by three different methods.
    Panel A: relative error in c as estimated by heuristic, saddle point Fisher, and full Fisher
    Panel B: relative error in koff as estimated by heuristic, saddle point Fisher, and full Fisher
    """
    figname = 'combined_error_composite_compare'
    curveL1 = DATADICT[figname + '_c_heuristic']
    curveL2 = DATADICT[figname + '_c_saddlePointFisher']
    curveL3 = DATADICT[figname + '_c_fullFisher']
    curveR1 = DATADICT[figname + '_koff_heuristic']
    curveR2 = DATADICT[figname + '_koff_saddlePointFisher']
    curveR3 = DATADICT[figname + '_koff_fullFisher']
    # set up axes
    fig, axarr = plt.subplots(nrows=1, ncols=2)
    plt.suptitle('Combined: MLE Relative error comparison ($k_p=10$, $t=100$, $k_{off}=1$)')
    axarr[0].set_xlabel(r'$c$')
    axarr[1].set_xlabel(r'$k_{off}$')
    axarr[0].set_ylabel('Relative error')
    axarr[0].set_xscale('log')
    axarr[1].set_xscale('log')
    # fix to same ylim
    ylow = -0.05
    yhigh = 1.5
    axarr[0].set_ylim(ylow, yhigh)
    axarr[1].set_ylim(ylow, yhigh)
    # left plot
    axarr[0].set_title(r'$\langle\delta c^{2}\rangle$/$c^{2}$')
    axarr[0].plot(curveL1['xpts'], curveL1['ypts'], 'k', label='Simple Fisher', zorder=1)
    axarr[0].scatter(curveL2['xpts'], curveL2['ypts'], c='b', edgecolor='', label='Numeric Fisher (Saddle Point)', zorder=2)
    axarr[0].scatter(curveL3['xpts'], curveL3['ypts'], c='r', edgecolor='', label='Numeric Fisher (Full)', zorder=3)
    axarr[0].legend(fontsize=8)
    # right plot
    axarr[1].set_title(r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$')
    axarr[1].plot(curveR1['xpts'], curveR1['ypts'], 'k', label='Simple Fisher', zorder=1)
    axarr[1].scatter(curveR2['xpts'], curveR2['ypts'], c='b', edgecolor='', label='Numeric Fisher (Saddle Point)', zorder=2)
    axarr[1].scatter(curveR3['xpts'], curveR3['ypts'], c='r', edgecolor='', label='Numeric Fisher (Full)', zorder=3)
    axarr[1].legend(fontsize=8)
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_detail' + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_detail' + '.eps')


def make_figure_C3():
    """
    KPR model relative error in c and koff, as estimated by two different methods.
    Panel A: relative error in c as estimated by heuristic, and saddle point Fisher
    Panel B: relative error in koff as estimated by heuristic, and saddle point Fisher
    """
    figname = 'KPR_error_composite_compare'
    curveL1 = DATADICT[figname + '_c_heuristic']
    curveL2 = DATADICT[figname + '_c_saddlePointFisher']
    curveR1 = DATADICT[figname + '_koff_heuristic']
    curveR2 = DATADICT[figname + '_koff_saddlePointFisher']
    # set up axes
    fig, axarr = plt.subplots(nrows=1, ncols=2)
    plt.suptitle('KPR: MLE Relative error comparison ($k_p=10$, $t=100$, $k_{off}=1$, $k_f=100$)')
    axarr[0].set_xlabel(r'$c$')
    axarr[1].set_xlabel(r'$k_{off}$')
    axarr[0].set_ylabel('Relative error')
    axarr[0].set_xscale('log')
    axarr[1].set_xscale('log')
    # fix to same ylim
    ylow = -0.1
    yhigh = 2.5
    axarr[0].set_ylim(ylow, yhigh)
    axarr[1].set_ylim(ylow, yhigh)
    # left plot
    axarr[0].set_title(r'$\langle\delta c^{2}\rangle$/$c^{2}$')
    axarr[0].plot(curveL1['xpts'], curveL1['ypts'], 'k', label='Simple Fisher', zorder=1)
    axarr[0].scatter(curveL2['xpts'], curveL2['ypts'], c='b', edgecolor='', label='Numeric Fisher (Saddle Point)', zorder=2)
    axarr[0].legend(fontsize=8)
    # right plot
    axarr[1].set_title(r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$')
    axarr[1].plot(curveR1['xpts'], curveR1['ypts'], 'k', label='Simple Fisher', zorder=1)
    axarr[1].scatter(curveR2['xpts'], curveR2['ypts'], c='b', edgecolor='', label='Numeric Fisher (Saddle Point)', zorder=2)
    axarr[1].legend(fontsize=8)
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_detail' + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_detail' + '.eps')


def make_figure_D1():
    """
    Mode 1 with a prior, plotting MLE for x vs x and comparing to heuristic (ie. non-prior prediction)
    One Panel.
    """
    figname = 'mode1_MLE_compare'
    curve1 = DATADICT[figname + '_heuristic']
    curve2 = DATADICT[figname + '_prior']
    # plot setup
    #plt.figure(figsize=(10, 5))
    plt.figure()
    plt.title(r'Mode 1: MLE comparison non-uniform prior ($k_p=10$, $t=100$, $k_{off}=1$, $a=0.001$)')
    plt.plot(curve1['xpts'][0:399], curve1['ypts'][0:399], 'r', label='heuristic', zorder=1)
    plt.plot(curve1['xpts'][400:], curve1['ypts'][400:], 'r', zorder=1)
    plt.scatter(curve2['xpts'], curve2['ypts'], c='b', edgecolor='', label='numeric with prior', zorder=2)
    plt.xlabel(r'$n_{obs}$')
    plt.ylabel(r'$x_{MLE}$')
    plt.legend()
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_prior' + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_prior' + '.eps')


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
