import matplotlib.pyplot as plt
from numpy import divide
import matplotlib.gridspec as gridspec
import os
import seaborn as sns

from load_inputs import DATADICT
from settings import DIR_OUTPUT
from settings import COLOR_SCHEME as cs

plt.style.use('parameters.mplstyle')  # particularIMporting

def make_figure_2():
    """
    Signal specificity and relative error in estimation of x for Mode 1.
    Panel A: <n> vs c for two different Kd, showing overlapping "clouds"
    Panel B: rel_err_x vs x for heuristic form of relative error
    Panel C: over-layed cross sections from heatmaps
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
    fig.set_size_inches(7.2, 3.2)
    #plt.suptitle(r'Mode 1: Mean $n$ and relative error in $x$ estimate')
    # left plot
    colour_palette = sns.color_palette("muted", 2)
    axarr[0].plot(curveLk['xpts'], curveLk['ypts'], c=colour_palette[0], label=r'$k_{off}=10$')
    axarr[0].plot(curveLkUp['xpts'], curveLkUp['ypts'], c=colour_palette[0], linestyle='--')
    axarr[0].plot(curveLkLow['xpts'], curveLkLow['ypts'], c=colour_palette[0], linestyle='--')
    axarr[0].plot(curveLr['xpts'], curveLr['ypts'], c=colour_palette[1], label=r'$k_{off}=15$')
    axarr[0].plot(curveLrUp['xpts'], curveLrUp['ypts'], c=colour_palette[1], linestyle='--')
    axarr[0].plot(curveLrLow['xpts'], curveLrLow['ypts'], c=colour_palette[1], linestyle='--')
    axarr[0].set_xlabel(r'$c$')
    axarr[0].set_ylabel(r'$\langle n\rangle/k_pt$')
    axarr[0].legend()
    ##axarr[0].set_title(r'($k_p=10$, $t=120$)')
    # right plot
    axarr[1].plot(curveR['xpts'], curveR['ypts'], 'k')
    axarr[1].set_xlabel(r'$x$')
    axarr[1].set_ylabel(r'$\langle\delta x^{2}\rangle$/$x^{2}$')
    axarr[1].set_ylim([0, 1.0])
    ##axarr[1].set_title('($k_p=10$, $t=1000$, $k_{off}=1$)')
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
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
    plt.figure(figsize=(4, 3))
    plt.plot(curve1['xpts'], curve1['ypts'], color=cs['c'], label=r'$\langle\delta c^{2}\rangle$/$c^{2}$')
    plt.plot(curve2['xpts'], curve2['ypts'], color=cs['koff'], label=r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$')
    plt.xlabel(r'$c$')
    plt.ylabel('Relative error')
    #plt.title('Combined: MLE Relative error for $c$ and $k_{off}$ ($k_p=10$, $t=1000$, $k_{off}=1$)')
    plt.legend()
    # set limits
    plt.ylim(0, 1.3)
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
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
    plt.figure(figsize=(4, 3))
    plt.plot(curve1['xpts'], curve1['ypts'], 'k', label=r'$\langle\delta c^{2}\rangle$/$c^{2}$')
    plt.plot(curve2['xpts'], curve2['ypts'], 'r', label=r'$\langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$')
    plt.xlabel(r'$c$')
    plt.ylabel('Relative error')
    #plt.title('KPR: MLE Relative error for $c$ and $k_{off}$ ($k_p=10$, $t=1000$, $k_{off}=1$, $k_f=100$)')
    plt.legend()
    # set limits
    #plt.ylim(0, 0.3)
    #plt.xlim(0, 10)
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


def make_figure_B1():
    """
    'mode1_MLE_compare' figure, 2 curves (numeric vs heuristic)
    """
    figname = 'mode1_MLE_compare'
    curve1 = DATADICT[figname + '_numeric']
    curve2 = DATADICT[figname + '_heuristic']
    # plot
    plt.figure(figsize=(6, 4))
    c2_part1 = plt.plot(curve2['xpts'][0:400], curve2['ypts'][0:400], color=cs['heuristic'], label='heuristic')
    #c2_part2 = plt.plot(curve2['xpts'][400:], curve2['ypts'][400:], color=cs['heuristic'])
    c1 = plt.plot(curve1['xpts'], curve1['ypts'], marker='o', linestyle='None', color=cs['numerical_fisher_sp'], label='numeric')
    #plt.title('Mode 1 MLE: Numeric vs Heuristic')# ($k_p=10$, $t=1000$, $k_{off}=1$)')
    plt.xlabel(r'$n$')
    plt.ylabel(r'$x^{*}$')
    plt.legend()
    # save figure
    plt.gca().set_ylim([-5, max(curve1['ypts'])])
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
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
    fig = plt.figure()
    gs1 = gridspec.GridSpec(1, 2)
    axarr = [fig.add_subplot(ss) for ss in gs1]
    fig.set_size_inches(8, 4)
    #plt.suptitle('Combined MLE: Numeric vs Heuristic')# ($k_p=10$, $t=1000$, $k_{off}=1$, $m=100$)')
    # left subplot
    cL2 = axarr[0].plot(curveL2['xpts'], curveL2['ypts'], color=cs['heuristic'], label='heuristic')
    cL1 = axarr[0].plot(curveL1['xpts'], curveL1['ypts'], marker='o', linestyle='None', color=cs['numerical_fisher_sp'], label='numeric')
    ##axarr[0].set_title(r'Estimate for $c$')
    axarr[0].set_xlabel(r'$n$')
    axarr[0].set_ylabel(r'$c^*$')
    axarr[0].legend()
    # right subplot
    cR2 = axarr[1].plot(curveR2['xpts'], curveR2['ypts'], color=cs['heuristic'], label='heuristic')
    cR1 = axarr[1].plot(curveR1['xpts'][1:], curveR1['ypts'][1:], marker='o', linestyle='None', color=cs['numerical_fisher_sp'], label='numeric')
    ##axarr[1].set_title(r'Estimate for $k_{off}$')
    axarr[1].set_xlabel(r'$n$')
    axarr[1].set_ylabel(r'$k_{off}^*$')
    axarr[1].legend()
    # save figure
    axarr[1].set_ylim([0, max(curveR1['ypts'][1:]) * 1.1])
    gs1.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
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
    fig = plt.figure()
    gs1 = gridspec.GridSpec(1, 2)
    axarr = [fig.add_subplot(ss) for ss in gs1]
    fig.set_size_inches(8, 4)
    #plt.suptitle('KPR MLE: Numeric vs Heuristic')# ($k_p=10$, $t=1000$, $k_{off}=1$, $k_f=100$, $m=100$)')
    # left subplot
    cL2 = axarr[0].plot(curveL2['xpts'], curveL2['ypts'], color=cs['heuristic'], label='heuristic')
    cL1 = axarr[0].plot(curveL1['xpts'], curveL1['ypts'], marker='o', linestyle='None', color=cs['numerical_fisher_sp'], label='numeric')
    ##axarr[0].set_title(r'Estimate for $c$')
    axarr[0].set_xlabel(r'$n$')
    axarr[0].set_ylabel(r'$c^*$')
    axarr[0].legend()
    # right subplot
    cR2 = axarr[1].plot(curveR2['xpts'], curveR2['ypts'], color=cs['heuristic'], label='heuristic')
    cR1 = axarr[1].plot(curveR1['xpts'][1:], curveR1['ypts'][1:], marker='o', linestyle='None', color=cs['numerical_fisher_sp'], label='numeric')
    ##axarr[1].set_title(r'Estimate for $k_{off}$')
    axarr[1].set_xlabel(r'$n$')
    axarr[1].set_ylabel(r'$k_{off}^*$')
    axarr[1].legend()
    # save figure
    axarr[0].set_ylim([-5, max(curveL1['ypts']) * 1.1])
    axarr[1].set_ylim([0, max(curveR1['ypts'][1:]) * 1.1])
    gs1.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


def make_figure_B4():
    """
    'KPR2_MLE_composite_compare' figure, 4 curves (numeric vs heuristic) for c and koff
    """
    figname = 'KPR2_MLE_composite_compare'
    curveL1 = DATADICT[figname + '_c_numeric']
    curveL2 = DATADICT[figname + '_c_heuristic']
    curveR1 = DATADICT[figname + '_koff_numeric']
    curveR2 = DATADICT[figname + '_koff_heuristic']
    # plot
    fig = plt.figure()
    gs1 = gridspec.GridSpec(1, 2)
    axarr = [fig.add_subplot(ss) for ss in gs1]
    fig.set_size_inches(8, 4)
    #plt.suptitle('KPR2 MLE: Numeric vs Heuristic')# ($k_p=10$, $t=1000$, $k_{off}=1$, $k_f=100$, $m=100$)')
    # left subplot
    cL2 = axarr[0].plot(curveL2['xpts'], curveL2['ypts'], color=cs['heuristic'], label='heuristic')
    cL1 = axarr[0].plot(curveL1['xpts'], curveL1['ypts'], marker='o', linestyle='None', color=cs['numerical_fisher_sp'], label='numeric')
    ##axarr[0].set_title(r'Estimate for $c$')
    axarr[0].set_xlabel(r'$n$')
    axarr[0].set_ylabel(r'$c^*$')
    axarr[0].legend()
    # right subplot
    cR2 = axarr[1].plot(curveR2['xpts'], curveR2['ypts'], color=cs['heuristic'], label='heuristic')
    cR1 = axarr[1].plot(curveR1['xpts'][1:], curveR1['ypts'][1:], marker='o', linestyle='None', color=cs['numerical_fisher_sp'], label='numeric')
    ##axarr[1].set_title(r'Estimate for $k_{off}$')
    axarr[1].set_xlabel(r'$n$')
    axarr[1].set_ylabel(r'$k_{off}^*$')
    axarr[1].legend()
    # save figure
    axarr[0].set_ylim([-5, max(curveL1['ypts']) * 1.1])
    axarr[1].set_ylim([0, max(curveR1['ypts'][1:]) * 1.1])
    gs1.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
    return


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
    #plt.figure(figsize=(8, 4))
    plt.figure()
    plt.plot(curve1['xpts'], curve1['ypts'], color=cs['simple_fisher'], label='Simple Fisher', zorder=1)
    plt.scatter(curve2['xpts'], curve2['ypts'], color=cs['numerical_fisher_sp'],  edgecolor='', label='Saddle Point Fisher', zorder=2)
    plt.scatter(curve3['xpts'], curve3['ypts'], color=cs['numerical_fisher'], edgecolor='', label='Numeric Fisher', zorder=3)
    # axis
    #plt.title('Mode 1: MLE Relative error comparison')# ($k_p=10$, $t=1000$, $k_{off}=1$)')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\alpha \langle\delta x^{2}\rangle$/$x^{2}$')
    plt.gca().set_xscale('log')
    plt.ylim([0, 1500])
    plt.legend()
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
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
    fig = plt.figure()
    gs1 = gridspec.GridSpec(1, 2)
    axarr = [fig.add_subplot(ss) for ss in gs1]
    fig.set_size_inches(8, 4)
    #plt.suptitle('Combined: MLE Relative error comparison')# ($k_p=10$, $t=1000$, $k_{off}=1$)')
    axarr[0].set_xlabel(r'$c$')
    axarr[1].set_xlabel(r'$c$')
    #axarr[0].set_ylabel('Relative error')
    axarr[0].set_xscale('log')
    axarr[1].set_xscale('log')
    # fix to same ylim
    ylow = -0.05
    yhigh = 1500
    axarr[0].set_ylim(ylow, yhigh)
    axarr[1].set_ylim(ylow, yhigh)
    # left plot
    ##axarr[0].set_title(r'$\alpha \langle\delta c^{2}\rangle$/$c^{2}$')
    axarr[0].plot(curveL1['xpts'], curveL1['ypts'], color=cs['simple_fisher'], label='Simple Fisher', zorder=1)
    axarr[0].scatter(curveL2['xpts'], curveL2['ypts'], marker='o', color=cs['numerical_fisher_sp'], edgecolor='', label='Saddle Point Fisher', zorder=2)
    axarr[0].scatter(curveL3['xpts'], curveL3['ypts'], marker='o', color=cs['numerical_fisher'], edgecolor='', label='Numeric Fisher', zorder=3)
    axarr[0].legend(fontsize=8)
    # right plot
    ##axarr[1].set_title(r'$\alpha \langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$')
    axarr[1].plot(curveR1['xpts'], curveR1['ypts'], color=cs['simple_fisher'], label='Simple Fisher', zorder=1)
    axarr[1].scatter(curveR2['xpts'], curveR2['ypts'], marker='o', color=cs['numerical_fisher_sp'], edgecolor='', label='Saddle Point Fisher', zorder=2)
    axarr[1].scatter(curveR3['xpts'], curveR3['ypts'], marker='o', color=cs['numerical_fisher'], edgecolor='', label='Numeric Fisher', zorder=3)
    axarr[1].legend(fontsize=8)
    # save figure
    gs1.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_detail' + '.pdf', transparent=True)
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
    fig = plt.figure()
    gs1 = gridspec.GridSpec(1, 2)
    axarr = [fig.add_subplot(ss) for ss in gs1]

    fig.set_size_inches(8, 4)
    #plt.suptitle('KPR: MLE Relative error comparison')# ($k_p=10$, $t=1000$, $k_{off}=1$, $k_f=100$)')
    axarr[0].set_xlabel(r'$c$')
    axarr[1].set_xlabel(r'$c$')
    axarr[0].set_ylabel(r'$\alpha \langle \delta c^2 \rangle /c^2$', fontdict={'fontname':'Arial'})
    axarr[1].set_ylabel(r'$\alpha \langle \delta k_{off}^2 \rangle /k_{off}^2$', fontdict={'fontname': 'Arial'})
    axarr[0].set_xscale('log')
    axarr[1].set_xscale('log')
    # fix to same ylim
    ylow = -0.1
    yhigh = 1200
    axarr[0].set_ylim(ylow, yhigh)
    axarr[1].set_ylim(ylow, yhigh)
    # left plot
    #axarr[0].set_title(r'$\alpha \langle\delta c^{2}\rangle$/$c^{2}$')
    axarr[0].plot(curveL1['xpts'], curveL1['ypts'], color=cs['simple_fisher'], label='Simple Fisher', zorder=1)
    axarr[0].scatter(curveL2['xpts'], curveL2['ypts'], marker='o', color=cs['numerical_fisher_sp'], edgecolor='', label='Saddle Point Fisher', zorder=2)
    axarr[0].legend(fontsize=8)
    # right plot
    #axarr[1].set_title(r'$\alpha \langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$')
    axarr[1].plot(curveR1['xpts'], curveR1['ypts'], color=cs['simple_fisher'], label='Simple Fisher', zorder=1)
    axarr[1].scatter(curveR2['xpts'], curveR2['ypts'], marker='o', color=cs['numerical_fisher_sp'], edgecolor='', label='Saddle Point Fisher', zorder=2)
    axarr[1].legend(fontsize=8)
    # save figure
    #fig.set_size_inches(8.0, 7.0)
    gs1.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_detail' + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_detail' + '.eps')


def make_figure_C4():
    """
    Compare HEURISTIC relative error in c and koff, as estimated by different methods.
    Panel A: relative error in c as estimated by heuristic {KPR, KPR2}
    Panel B: relative error in koff as estimated by heuristic {KPR, KPR2}
    """
    figname = 'multimodel_heuristic_error_composite_compare'
    curve1name = 'KPR_error_composite_compare'
    curve2name = 'KPR2_error_composite_compare'
    curveL1 = DATADICT[curve1name + '_c_heuristic']
    curveL2 = DATADICT[curve2name + '_c_heuristic']
    curveR1 = DATADICT[curve1name + '_koff_heuristic']
    curveR2 = DATADICT[curve2name + '_koff_heuristic']

    # set up axes
    fig = plt.figure()
    gs1 = gridspec.GridSpec(1, 2)
    axarr = [fig.add_subplot(ss) for ss in gs1]

    fig.set_size_inches(8, 4)
    #plt.suptitle('KPR vs KPR2 - MLE heuristic relative error')# ($k_p=10$, $t=1000$, $k_{off}=1$, $k_f=100$)')
    axarr[0].set_xlabel(r'$c$')
    axarr[1].set_xlabel(r'$c$')
    axarr[0].set_ylabel('Relative error', fontdict={'fontname':'Arial'})
    axarr[0].set_xscale('log')
    axarr[1].set_xscale('log')
    # fix to same ylim
    ylow = -0.1
    yhigh = 1200
    axarr[0].set_ylim(ylow, yhigh)
    axarr[1].set_ylim(ylow, yhigh)
    # left plot
    #axarr[0].set_title(r'$\alpha \langle\delta c^{2}\rangle$/$c^{2}$')
    axarr[0].plot(curveL1['xpts'], curveL1['ypts'], color=cs['simple_fisher'], marker='o', label='Simple Fisher KPR', zorder=1)
    axarr[0].plot(curveL2['xpts'], curveL2['ypts'], color='blue', marker='^', label='Simple Fisher KPR2', zorder=1)
    axarr[0].legend(fontsize=8)
    # right plot
    #axarr[1].set_title(r'$\alpha \langle\delta k_{off}^{2}\rangle$/$k_{off}^{2}$')
    axarr[1].plot(curveR1['xpts'], curveR1['ypts'], color=cs['simple_fisher'], marker='o', label='Simple Fisher KPR', zorder=1)
    axarr[1].plot(curveR2['xpts'], curveR2['ypts'], color='blue', marker='^', label='Simple Fisher KPR2', zorder=1)
    axarr[1].legend(fontsize=8)
    # save figure
    #fig.set_size_inches(8.0, 7.0)
    gs1.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_detail' + '.pdf', transparent=True)
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
    #plt.figure(figsize=(8, 4))
    plt.figure()
    #plt.title(r'Mode 1: MLE comparison non-uniform prior ')#($k_p=10$, $t=1000$, $k_{off}=1$, $\lambda=0.001$)')
    plt.plot(curve1['xpts'][0:399], curve1['ypts'][0:399], color=cs['heuristic'], label='heuristic', zorder=1)
    plt.plot(curve1['xpts'][400:], curve1['ypts'][400:], color=cs['heuristic'], zorder=1)
    plt.scatter(curve2['xpts'], curve2['ypts'], color=cs['numerical_fisher_sp'], edgecolor='', label='numeric with prior', zorder=2)
    plt.xlabel(r'$n$')
    plt.ylabel(r'$x^{*}$')
    plt.legend()
    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_prior' + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '_prior' + '.eps')


def make_figure_E1():
    """
    Combined estimating one parmaeter vs estimating two parameters vs mode 1 estimate
    One Panel.
    """
    figname = 'Supplementary_E1'
    curve1 = DATADICT[figname + '_mode1']
    curve2 = DATADICT[figname + '_TwoParam']
    curve3 = DATADICT[figname + '_OneParam']
    # plot setup
    plt.figure(figsize=(6, 4))
    #plt.title('Relative error for different receptor schemes')#+'\r\n'+r'($k_p=10$, $t=1000$, $k_{off}=1$, $k_{on}=1$)')
    plt.plot(curve1['xpts'], divide(curve1['ypts'], curve1['ypts']), color=cs['heuristic'], label='Mode 1', zorder=1)
    plt.plot(curve2['xpts'], divide(curve2['ypts'], curve1['ypts']), color=cs['simple_fisher'], label='Combined Estimating Two Parameters', zorder=1)
    plt.plot(curve3['xpts'], divide(curve3['ypts'], curve1['ypts']), color=cs['numerical_fisher_sp'], label='Combined Estimating One Parameter', zorder=1)
    plt.xlabel(r'$c$')
    plt.ylabel(r'$\alpha \langle \delta c^{2}\rangle /c^{2}$')
    plt.legend()
    plt.xlim((0, 20))
    #plt.ylim((0, 0.1))

    # save figure
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')


if __name__ == "__main__":
    #make_figure_2()
    #make_figure_3()
    #make_figure_5()
    #make_figure_B1()
    #make_figure_B2()
    #make_figure_B3()
    make_figure_B4()
    #make_figure_C1()
    #make_figure_C2()
    #make_figure_C3()
    #make_figure_C4()
    #make_figure_D1()
    #make_figure_E1()
