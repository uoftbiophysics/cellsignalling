import equations as eqns
from settings import KON, KP, T, KF

# TODO other dictionary for long time and high G

def trace(diagonal1, diagonal2, c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return (diagonal1(c, koff, kon, T, KF, KP)+diagonal1(c, koff, kon, T, KF, KP))

def eigenvalues(diagonal1, diagonal2, determinant, c, koff, kon=KON, T=T, KF=KF, KP=KP):
    tr = trace(diagonal1, diagonal2, c, koff, kon, T, KF, KP); det = determinant(c, koff, kon, T, KF, KP);
    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) ), ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )

"""
Here there are 2 types of dictionaries: ratios and 1 equation.

ratios have keywords:
    'subdir2': which subdir to create and put figures in
    'log': whether or not it is to be plotted in log-scales
    'plots': all the plots to use
        'name_of_figure': you can change this name
            'num': numerator of ratio
            'denom': denominator of ratio
            'label': what to label the colorbar

1 equation dictionary have keywords:
    'subdir2': which subdir to create and put figures in
    'log': whether or not it is to be plotted in log-scales
    'plots': all the plots to use
        'name_of_figure': you can change this name
            'eqn': the equation to plot
            'label': what to label the colorbar

Essentially very similar, however in one you declare the numerator and denominator equation, in the other 1 equation will suffice.
"""

# Ratios

DIAGANDTRACE_RATIO = { 'subdir2' : 'diagonal_ratio', 'log' : True,
                   'plots' : {
                               'ratioTraceSigma2' : { 'num' : eqns.traceSigmacrlb2NoTrace, 'denom' : eqns.traceSigmacrlb2, 'label' : r'Model 2 tr($\Sigma^A$) / tr($\Sigma^B$)'},
                               'ratioTraceSigma3' : { 'num' : eqns.traceSigmacrlb3NoTrace, 'denom' : eqns.traceSigmacrlb3, 'label' : r'Model 3 tr($\Sigma^A$) / tr($\Sigma^B$)'},
                               'ratioTraceSigma4' : { 'num' : eqns.traceSigmacrlb4NoTrace, 'denom' : eqns.traceSigmacrlb4, 'label' : r'Model 4 tr($\Sigma^A$) / tr($\Sigma^B$)'},
                               'ratioSigmaX1' : { 'num' : eqns.Sigmacrlb1NoTrace, 'denom' : eqns.Sigmacrlb1, 'label' : r'Model 1 $\langle\delta x^{2}\rangle^A$/$\langle\delta x^{2}\rangle^B$'},
                               'ratioSigmaC2' : { 'num' : eqns.SigmacrlbC2NoTrace, 'denom' : eqns.SigmacrlbC2, 'label' : r'Model 2 $\langle\delta c^{2}\rangle^A$/$\langle\delta c^{2}\rangle^B$'},
                               'ratioSigmaC3' : { 'num' : eqns.SigmacrlbC3NoTrace, 'denom' : eqns.SigmacrlbC3, 'label' : r'Model 3 $\langle\delta c^{2}\rangle^A$/$\langle\delta c^{2}\rangle^B$'},
                               'ratioSigmaC4' : { 'num' : eqns.SigmacrlbC4NoTrace, 'denom' : eqns.SigmacrlbC4, 'label' : r'Model 4 $\langle\delta c^{2}\rangle^A$/$\langle\delta c^{2}\rangle^B$'},
                               'ratioSigmaK2' : { 'num' : eqns.SigmacrlbK2NoTrace, 'denom' : eqns.SigmacrlbK2, 'label' : r'Model 2 $\langle\delta k_{\mathrm{off}}^{2}\rangle/^A\langle\delta {k_{\mathrm{off}}^{2}}\rangle^B$'},
                               'ratioSigmaK3' : { 'num' : eqns.SigmacrlbK3NoTrace, 'denom' : eqns.SigmacrlbK3, 'label' : r'Model 3 $\langle\delta k_{\mathrm{off}}^{2}\rangle/^A\langle\delta {k_{\mathrm{off}}^{2}}\rangle^B$'},
                               'ratioSigmaK4' : { 'num' : eqns.SigmacrlbK4NoTrace, 'denom' : eqns.SigmacrlbK4, 'label' : r'Model 4 $\langle\delta k_{\mathrm{off}}^{2}\rangle/^A\langle\delta {k_{\mathrm{off}}^{2}}\rangle^B$'}
                             }
                  }

DETANDEVAL_RATIO = { 'subdir2' : 'determinant_ratio', 'log' : True, #some of the eigenvalues are 0!
                   'plots' : {
                                'ratioEval1Sigma2' : { 'num' : eqns.evalplusSigmacrlb2NoTrace, 'denom' : eqns.evalplusSigmacrlb2, 'label' : r'Model 2 $\lambda_{+}$ / ${\lambda_{+}}_{crlb}$'},
                                'ratioEval2Sigma2' : { 'num' : eqns.evalminusSigmacrlb2NoTrace, 'denom' : eqns.evalminusSigmacrlb2, 'label' : r'Model 2 $\lambda_{-}$ / ${\lambda_{-}}_{crlb}$'},
                                'ratioEval1Sigma3' : { 'num' : eqns.evalplusSigmacrlb3NoTrace, 'denom' : eqns.evalplusSigmacrlb3, 'label' : r'Model 3 $\lambda_{+}$ / ${\lambda_{+}}_{crlb}$'},
                                'ratioEval2Sigma3' : { 'num' : eqns.evalminusSigmacrlb3NoTrace, 'denom' : eqns.evalminusSigmacrlb3, 'label' : r'Model 3 $\lambda_{-}$ / ${\lambda_{-}}_{crlb}$'},
                                'ratioEval1Sigma4' : { 'num' : eqns.evalplusSigmacrlb4NoTrace, 'denom' : eqns.evalplusSigmacrlb4, 'label' : r'Model 4 $\lambda_{+}$ / ${\lambda_{+}}_{crlb}$'},
                                'ratioEval2Sigma4' : { 'num' : eqns.evalminusSigmacrlb4NoTrace, 'denom' : eqns.evalminusSigmacrlb4, 'label' : r'Model 4 $\lambda_{-}$ / ${\lambda_{-}}_{crlb}$'},
                               'ratioDetSigma2': { 'num' : eqns.DetSigmacrlb2NoTrace, 'denom' : eqns.DetSigmacrlb2, 'label' : r'Model 2 $| \Sigma^A |$/$| \Sigma^B |$'},
                               'ratioDetSigma3': { 'num' : eqns.DetSigmacrlb3NoTrace, 'denom' : eqns.DetSigmacrlb3, 'label' : r'Model 3 $| \Sigma^A |$/$| \Sigma^B |$'},
                               'ratioDetSigma4': { 'num' : eqns.DetSigmacrlb4NoTrace, 'denom' : eqns.DetSigmacrlb4, 'label' : r'Model 4 $| \Sigma^A |$/$| \Sigma^B |$'}
                             }
                  }


COV_RATIO = { 'subdir2' : 'cov_ratio', 'log' : False,
                   'plots' : {
                               'ratioSigmacrlbCK2': { 'num' : eqns.SigmacrlbCK2NoTrace, 'denom' : eqns.SigmacrlbCK2, 'label' : r'Model 2 $\langle{\delta c \delta k_{\mathrm{off}}}\rangle$ / $\langle{\delta c \delta k_{\mathrm{off}}}_{crlb}\rangle$'},
                               'ratioSigmacrlbCK3': { 'num' : eqns.SigmacrlbCK3NoTrace, 'denom' : eqns.SigmacrlbCK3, 'label' : r'Model 2 $\langle{\delta c \delta k_{\mathrm{off}}}\rangle$ / $\langle{\delta c \delta k_{\mathrm{off}}}_{crlb}\rangle$'},
                               'ratioSigmacrlbCK4': { 'num' : eqns.SigmacrlbCK4NoTrace, 'denom' : eqns.SigmacrlbCK4, 'label' : r'Model 2 $\langle{\delta c \delta k_{\mathrm{off}}}\rangle$ / $\langle{\delta c \delta k_{\mathrm{off}}}_{crlb}\rangle$'}
                             }
                  }


############################ No ratios No Trace ######################################

RELERRANDTRACE_NOTRACE = { 'subdir2' : 'relerror_notrace', 'log' : True,
                   'plots' : {
                               'relErrX1': { 'eqn' : eqns.RelErrorX1NoTrace, 'label' : r'Model 1 $\langle\delta x^{2}\rangle /x^{2}$'},
                               'relErrC2': { 'eqn' : eqns.RelErrC2NoTrace, 'label' : r'Model 2 $\langle\delta c^{2}\rangle /c^{2}$'},
                               'relErrK2': { 'eqn' : eqns.RelErrK2NoTrace, 'label' : r'Model 2 $\langle{\delta k_{\mathrm{off}}^2}\rangle /k_{\mathrm{off}}^{2}$'},
                               'relErrC3': { 'eqn' : eqns.RelErrC3NoTrace, 'label' : r'Model 3 $\langle\delta c^{2}\rangle /c^{2}$'},
                               'relErrK3': { 'eqn' : eqns.RelErrK3NoTrace, 'label' : r'Model 3 $\langle{\delta k_{\mathrm{off}}^2}\rangle /k_{\mathrm{off}}^{2}$'},
                               'relErrC4': { 'eqn' : eqns.RelErrC4NoTrace, 'label' : r'Model 4 $\langle\delta c^{2}\rangle /c^{2}$'},
                               'relErrK4': { 'eqn' : eqns.RelErrK4NoTrace, 'label' : r'Model 4 $\langle{\delta k_{\mathrm{off}}^2}\rangle /k_{\mathrm{off}}^{2}$'},
                               'SigmacrlbC2': { 'eqn' : eqns.SigmacrlbC2NoTrace, 'label' : r'Model 2 $\langle\delta c^{2}\rangle$'},
                               'SigmacrlbK2': { 'eqn' : eqns.SigmacrlbK2NoTrace, 'label' : r'Model 2 $\langle{\delta k_{\mathrm{off}}^2}\rangle$'},
                               'SigmacrlbC3': { 'eqn' : eqns.SigmacrlbC3NoTrace, 'label' : r'Model 3 $\langle\delta c^{2}\rangle /c^{2}$'},
                               'SigmacrlbK3': { 'eqn' : eqns.SigmacrlbK3NoTrace, 'label' : r'Model 3 $\langle{\delta k_{\mathrm{off}}^2}\rangle$'},
                               'SigmacrlbC4': { 'eqn' : eqns.SigmacrlbC4NoTrace, 'label' : r'Model 4 $\langle\delta c^{2}\rangle$'},
                               'SigmacrlbK4': { 'eqn' : eqns.SigmacrlbK4NoTrace, 'label' : r'Model 4 $\langle{\delta k_{\mathrm{off}}^2}\rangle$'},
                               'traceSigma2': {'eqn' : eqns.traceSigmacrlb2NoTrace, 'label' : r'Model 2 tr($\Sigma$)'},
                               'traceSigma3': {'eqn' : eqns.traceSigmacrlb3NoTrace, 'label' : r'Model 3 tr($\Sigma$)'},
                               'traceSigma4': {'eqn' : eqns.traceSigmacrlb4NoTrace, 'label' : r'Model 4 tr($\Sigma$)'}
                             }
                  }

DETANDEVAL_NOTRACE = { 'subdir2' : 'det_notrace', 'log' : True,
                   'plots' : {
                               'detSigma2': {'eqn' : eqns.DetSigmacrlb2NoTrace, 'label' : r'Model 2 Det($\Sigma^A$)'},
                               'detSigma3': {'eqn' : eqns.DetSigmacrlb3NoTrace, 'label' : r'Model 3 Det($\Sigma^A$)'},
                               'detSigma4': {'eqn' : eqns.DetSigmacrlb4NoTrace, 'label' : r'Model 4 Det($\Sigma^A$)'},
                               'relDetSigma2': {'eqn' : eqns.RelDetSigmacrlb2NoTrace, 'label' : r'Model 2 Det($\Sigma^A$)/$c^2 k_{\mathrm{off}}^2$'},
                               'relDetSigma3': {'eqn' : eqns.RelDetSigmacrlb3NoTrace, 'label' : r'Model 3 Det($\Sigma^A$)/$c^2 k_{\mathrm{off}}^2$'},
                               'relDetSigma4': {'eqn' : eqns.RelDetSigmacrlb4NoTrace, 'label' : r'Model 4 Det($\Sigma^A$)/$c^2 k_{\mathrm{off}}^2$'},
                               'evalplusSigma2': {'eqn' : eqns.evalplusSigmacrlb2NoTrace, 'label' : r'Model 2 $\lambda_{+}$'},
                               'evalminusSigma2': {'eqn' : eqns.evalminusSigmacrlb2NoTrace, 'label' : r'Model 2 $\lambda_{-}$'},
                               'evalplusSigma3': {'eqn' : eqns.evalplusSigmacrlb3NoTrace, 'label' : r'Model 3 $\lambda_{+}$'},
                               'evalminusSigma3': {'eqn' : eqns.evalminusSigmacrlb3NoTrace, 'label' : r'Model 3 $\lambda_{-}$'},
                               'evalplusSigma4': {'eqn' : eqns.evalplusSigmacrlb4NoTrace, 'label' : r'Model 4 $\lambda_{+}$'},
                               'evalminusSigma4': {'eqn' : eqns.evalminusSigmacrlb4NoTrace, 'label' : r'Model 4 $\lambda_{-}$'},
                             }
                  }


COV_NOTRACE = { 'subdir2' : 'cov_notrace', 'log' : False,
                   'plots' : {
                               'SigmaCK2': { 'eqn' : eqns.SigmacrlbCK2, 'label' : r'Model 2 $\langle \delta c \delta k_{\mathrm{off}} \rangle$'},
                               'SigmaCK3': { 'eqn' : eqns.SigmacrlbCK3, 'label' : r'Model 3 $\langle \delta c \delta k_{\mathrm{off}} \rangle$ '},
                               'SigmaCK4': { 'eqn' : eqns.SigmacrlbCK4, 'label' : r'Model 4 $\langle \delta c \delta k_{\mathrm{off}} \rangle$'}
                             }
                  }

EVALS_RATIO = { 'subdir2' : 'notrace', 'log' : True,
                   'plots' : {
                              'ratioEvalsSigma2' : { 'num' : eqns.evalplusSigmacrlb2NoTrace, 'denom' : eqns.evalminusSigmacrlb2NoTrace, 'label' : r'Model 2 $\lambda_{+}$ / ${\lambda_{-}}_{crlb}$'},
                              'ratioEvalsSigma3' : { 'num' : eqns.evalplusSigmacrlb3NoTrace, 'denom' : eqns.evalminusSigmacrlb3NoTrace, 'label' : r'Model 3 $\lambda_{+}$ / ${\lambda_{-}}_{crlb}$'},
                              'ratioEvalsSigma4' : { 'num' : eqns.evalplusSigmacrlb4NoTrace, 'denom' : eqns.evalminusSigmacrlb4NoTrace, 'label' : r'Model 4 $\lambda_{+}$ / ${\lambda_{-}}_{crlb}$'}
                             }
                  }

############################ No ratios Full ######################################

DIAGANDTRACE_FULL = { 'subdir2' : 'full', 'log' : True,
                   'plots' : {
                               #'relErrX1': { 'eqn' : eqns.RelErrorX1, 'label' : r'Model 1 $\langle\delta x^{2}\rangle /x^{2}$'},
                               'SigmacrlbC2': { 'eqn' : eqns.SigmacrlbC2, 'label' : r'Model 2 $\langle \delta c^{2} \rangle$'},
                               'SigmacrlbK2': { 'eqn' : eqns.SigmacrlbK2, 'label' : r'Model 2 $\langle \delta k_{\mathrm{off}}^2 \rangle$'},
                               'SigmacrlbC3': { 'eqn' : eqns.SigmacrlbC3, 'label' : r'Model 3 $\langle\delta c^{2}\rangle$'},
                               'SigmacrlbK3': { 'eqn' : eqns.SigmacrlbK3, 'label' : r'Model 3 $\langle \delta k_{\mathrm{off}}^2 \rangle$'},
                               'SigmacrlbC4': { 'eqn' : eqns.SigmacrlbC4, 'label' : r'Model 4 $\langle \delta c^{2} \rangle$'},
                               'SigmacrlbK4': { 'eqn' : eqns.SigmacrlbK4, 'label' : r'Model 4 $\langle \delta k_{\mathrm{off}}^2 \rangle$'},
                               'traceSigma2': {'eqn' : eqns.traceSigmacrlb2, 'label' : r'Model 2 tr($\Sigma$)'},
                               'traceSigma3': {'eqn' : eqns.traceSigmacrlb3, 'label' : r'Model 3 tr($\Sigma$)'},
                               'traceSigma4': {'eqn' : eqns.traceSigmacrlb4, 'label' : r'Model 4 tr($\Sigma$)'}
                             }
                  }

DETANDEVAL_FULL = { 'subdir2' : 'full', 'log' : True,
                   'plots' : {
                               'detSigma2': {'eqn' : eqns.DetSigmacrlb2, 'label' : r'Model 2 Det($\Sigma$)'},
                               'detSigma3': {'eqn' : eqns.DetSigmacrlb3, 'label' : r'Model 3 Det($\Sigma$)'},
                               'detSigma4': {'eqn' : eqns.DetSigmacrlb4, 'label' : r'Model 4 det($\Sigma$)'},
                               'relDetSigma2': {'eqn' : eqns.RelDetSigmacrlb2, 'label' : r'Model 2 Det($\Sigma$)/$c^2 k_{\mathrm{off}}^2$'},
                               'relDetSigma3': {'eqn' : eqns.RelDetSigmacrlb3, 'label' : r'Model 3 Det($\Sigma$)/$c^2 k_{\mathrm{off}}^2$'},
                               'relDetSigma4': {'eqn' : eqns.RelDetSigmacrlb4, 'label' : r'Model 4 Det($\Sigma$)/$c^2 k_{\mathrm{off}}^2$'},
                               'evalplusSigma2': {'eqn' : eqns.evalplusSigmacrlb2, 'label' : r'Model 2 $\lambda_{+}$'},
                               'evalminusSigma2': {'eqn' : eqns.evalminusSigmacrlb2, 'label' : r'Model 2 $\lambda_{-}$'},
                               'evalplusSigma3': {'eqn' : eqns.evalplusSigmacrlb3, 'label' : r'Model 3 $\lambda_{+}$'},
                               'evalminusSigma3': {'eqn' : eqns.evalminusSigmacrlb3, 'label' : r'Model 3 $\lambda_{-}$'},
                               'evalplusSigma4': {'eqn' : eqns.evalplusSigmacrlb4, 'label' : r'Model 4 $\lambda_{+}$'},
                               'evalminusSigma4': {'eqn' : eqns.evalminusSigmacrlb4, 'label' : r'Model 4 $\lambda_{-}$'},
                             }
                  }


COV_FULL = { 'subdir2' : 'full', 'log' : False,
                   'plots' : {
                               'SigmaCK2': { 'eqn' : eqns.SigmacrlbCK2, 'label' : r'Model 2 $\langle \delta c \delta k_{\mathrm{off}} \rangle$'},
                               'SigmaCK3': { 'eqn' : eqns.SigmacrlbCK3, 'label' : r'Model 3 $\langle \delta c \delta k_{\mathrm{off}} \rangle$ '},
                               'SigmaCK4': { 'eqn' : eqns.SigmacrlbCK4, 'label' : r'Model 4 $\langle \delta c \delta k_{\mathrm{off}} \rangle$'}
                             }
                  }

############################ Supplemen ######################################

SI_RATIO = { 'subdir2' : 'SI_ratios', 'log' : True,
                   'plots' : {
                              'ratioSigmaX1' : { 'num' : eqns.Sigmacrlb1NoTrace, 'denom' : eqns.Sigmacrlb1, 'label' : r'Model 1 $\langle\delta x^{2}\rangle^A$/$\langle\delta x^{2}\rangle^B$'},
                              'ratioSigmaC3' : { 'num' : eqns.SigmacrlbC3NoTrace, 'denom' : eqns.SigmacrlbC3, 'label' : r'Model 3 $\langle\delta c^{2}\rangle^A$/$\langle\delta c^{2}\rangle^B$'},
                              'ratioSigmaC4' : { 'num' : eqns.SigmacrlbC4NoTrace, 'denom' : eqns.SigmacrlbC4, 'label' : r'Model 4 $\langle\delta c^{2}\rangle^A$/$\langle\delta c^{2}\rangle^B$'},
                              'ratioSigmaK3' : { 'num' : eqns.SigmacrlbK3NoTrace, 'denom' : eqns.SigmacrlbK3, 'label' : r'Model 3 $\langle\delta k_{\mathrm{off}}^{2}\rangle/^A\langle\delta {k_{\mathrm{off}}^{2}}\rangle^B$'},
                              'ratioSigmaK4' : { 'num' : eqns.SigmacrlbK4NoTrace, 'denom' : eqns.SigmacrlbK4, 'label' : r'Model 4 $\langle\delta k_{\mathrm{off}}^{2}\rangle/^A\langle\delta {k_{\mathrm{off}}^{2}}\rangle^B$'},
                              'ratioDetSigma3': { 'num' : eqns.DetSigmacrlb3NoTrace, 'denom' : eqns.DetSigmacrlb3, 'label' : r'Model 3 det$ (\Sigma^A)$/det$ (\Sigma^B)$'},
                              'ratioDetSigma4': { 'num' : eqns.DetSigmacrlb4NoTrace, 'denom' : eqns.DetSigmacrlb4, 'label' : r'Model 4 det$ (\Sigma^A)$/det$ (\Sigma^B)$'},
                              'ratioDetSigma2DetSima3_NoTrace': { 'num' : eqns.DetSigmacrlb2NoTrace, 'denom' : eqns.DetSigmacrlb3NoTrace, 'label' : r'Model 2 det$ (\Sigma^A)$ / Model 3 det$ (\Sigma^A )$'},
                              'ratioDetSigma2DetSima3': { 'num' : eqns.DetSigmacrlb2, 'denom' : eqns.DetSigmacrlb3, 'label' : r'Model 2 det$ (\Sigma^B)$ / Model 3 det$ (\Sigma^B )$'},
                              'ratioDetSigma3DetSima4_NoTrace': { 'num' : eqns.DetSigmacrlb3NoTrace, 'denom' : eqns.DetSigmacrlb4NoTrace, 'label' : r'Model 3 det$ (\Sigma^A)$ / Model 4 det$ (\Sigma^A )$'},
                              'ratioDetSigma3DetSima4': { 'num' : eqns.DetSigmacrlb3, 'denom' : eqns.DetSigmacrlb4, 'label' : r'Model 3 det$ (\Sigma^B)$ / Model 4 det$ (\Sigma^B)$'}
                             }
                  }

SI_ALT_RATIO = { 'subdir2' : 'SI_ratios', 'log' : True,
                   'plots' : {
                               'ratioSigmaC2' : { 'num' : eqns.SigmacrlbC2NoTrace, 'denom' : eqns.SigmacrlbC2, 'label' : r'$\langle\delta c^{2}\rangle^A$/$\langle\delta c^{2}\rangle^B$'},
                               'ratioSigmaK2' : { 'num' : eqns.SigmacrlbK2NoTrace, 'denom' : eqns.SigmacrlbK2, 'label' : r'$\langle\delta k_{\mathrm{off}}^{2}\rangle^A/\langle\delta {k_{\mathrm{off}}^{2}}\rangle^B$'},
                               'ratioDetSigma2': { 'num' : eqns.DetSigmacrlb2NoTrace, 'denom' : eqns.DetSigmacrlb2, 'label' : r'det$ (\Sigma^A)$/det$ (\Sigma^B)$'},
                             }
                  }

SI_RATIOS_subset = { 'subdir2' : 'SI_ratios_subset', 'log' : True,
                     'plots' : {
                         'ratioDetSigma3DetSima2_NoTrace': { 'num' : eqns.DetSigmacrlb3NoTrace, 'denom' : eqns.DetSigmacrlb2NoTrace, 'label' : r'Model 3 $| \Sigma^A |$ / Model 2 $| \Sigma^A |$'},
                         'ratioDetSigma3DetSima2': { 'num' : eqns.DetSigmacrlb3, 'denom' : eqns.DetSigmacrlb2, 'label' : r'Model 3 $| \Sigma^B |$ / Model 2 $| \Sigma^B |$'},
                     }}


SI_RATIOS_subset2 = { 'subdir2' : 'SI_ratios_subset2', 'log' : True,
                     'plots' : {
                         'kf_ratioDetSigma3DetSima2_NoTrace': { 'num' : eqns.DetSigmacrlb3NoTrace, 'denom' : eqns.DetSigmacrlb2NoTrace, 'label' : r'Model 3 $| \Sigma^A |$ / Model 2 $| \Sigma^A |$'},
                         'kf_ratioDetSigma3DetSima2': { 'num' : eqns.DetSigmacrlb3, 'denom' : eqns.DetSigmacrlb2, 'label' : r'Model 3 $| \Sigma^B |$ / Model 2 $| \Sigma^B |$'},
                     }}


##################### Figures Main############

MAIN = { 'subdir2' : 'main', 'log' : True,
                   'plots' : {
                   'relErrX1': { 'eqn' : eqns.dedimRelErrorX1NoTrace, 'label' : r'$k_p t \langle\delta x^{2}\rangle /x^{2}$'},
                   'relErrC2': { 'eqn' : eqns.dedimRelErrC2NoTrace, 'label' : r'$k_p t \langle\delta c^{2}\rangle /c^{2}$'},
                   'relErrKoff2': { 'eqn' : eqns.dedimRelErrK2NoTrace, 'label' : r'$k_p t \langle{\delta k_{\mathrm{off}}^2}\rangle /k_{\mathrm{off}}^{2}$'}#,
                   #'relErrC3': { 'eqn' : eqns.dedimRelErrC3NoTrace, 'label' : r'$k_p t \langle\delta c^{2}\rangle /c^{2}$'},
                   #'relErrKoff3': { 'eqn' : eqns.dedimRelErrK3NoTrace, 'label' : r'$k_p t \langle{\delta k_{\mathrm{off}}^2}\rangle /k_{\mathrm{off}}^{2}$'},
                   #'relErrC4': { 'eqn' : eqns.dedimRelErrC4NoTrace, 'label' : r'$k_p t \langle\delta c^{2}\rangle /c^{2}$'},
                   #'relErrKoff4': { 'eqn' : eqns.dedimRelErrK4NoTrace, 'label' : r'$k_p t \langle{\delta k_{\mathrm{off}}^2}\rangle /k_{\mathrm{off}}^{2}$'}
                   }
}
