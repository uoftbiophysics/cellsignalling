import equations as eqns
from heatmaps import KON, KP, T, KF

# TODO other dictionary for long time and high G

def trace(diagonal1, diagonal2, c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return (diagonal1(c, koff, kon, T, KF, KP)+diagonal1(c, koff, kon, T, KF, KP))

def eigenvalues(diagonal1, diagonal2, determinant, c, koff, kon=KON, T=T, KF=KF, KP=KP):
    tr = trace(diagonal1, diagonal2, c, koff, kon, T, KF, KP); det = determinant(c, koff, kon, T, KF, KP);
    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) ), ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )

# Ratios

DIAGANDTRACE_RATIO = { 'subdir2' : 'ratio_notrace_over_full', 'log' : True,
                   'plots' : {
                               'ratioSigmaX1' : { 'num' : eqns.Sigmacrlb1NoTrace, 'denom' : eqns.Sigmacrlb1, 'label' : r'Mode 1 $\langle\delta x^{2}\rangle$ / $\langle\delta x^{2}_{crlb}\rangle$'},
                               'ratioSigmaC2' : { 'num' : eqns.SigmacrlbC2NoTrace, 'denom' : eqns.SigmacrlbC2, 'label' : r'Mode 2 $\langle\delta c^{2}\rangle$ / $\langle\delta c^{2}_{crlb}\rangle$'},
                               'ratioSigmaC3' : { 'num' : eqns.SigmacrlbC3NoTrace, 'denom' : eqns.SigmacrlbC3, 'label' : r'Mode 3 $\langle\delta c^{2}\rangle$ / $\langle\delta c^{2}_{crlb}\rangle$'},
                               'ratioSigmaC4' : { 'num' : eqns.SigmacrlbC4NoTrace, 'denom' : eqns.SigmacrlbC4, 'label' : r'Mode 4 $\langle\delta c^{2}\rangle$ / $\langle\delta c^{2}_{crlb}\rangle$'},
                               'ratioSigmaK2' : { 'num' : eqns.SigmacrlbK2NoTrace, 'denom' : eqns.SigmacrlbK2, 'label' : r'Mode 2 $\langle\delta k_{off}^{2}\rangle / \langle\delta {k_{off}^{2}}_{crlb}\rangle$'},
                               'ratioSigmaK3' : { 'num' : eqns.SigmacrlbK3NoTrace, 'denom' : eqns.SigmacrlbK3, 'label' : r'Mode 3 $\langle\delta k_{off}^{2}\rangle / \langle\delta {k_{off}^{2}}_{crlb}\rangle$'},
                               'ratioSigmaK4' : { 'num' : eqns.SigmacrlbK4NoTrace, 'denom' : eqns.SigmacrlbK4, 'label' : r'Mode 4 $\langle\delta k_{off}^{2}\rangle / \langle\delta {k_{off}^{2}}_{crlb}\rangle$'},
                               'ratioTraceSigma2' : { 'num' : eqns.traceSigmacrlb2NoTrace, 'denom' : eqns.traceSigmacrlb2, 'label' : r'Mode 3 tr($\Sigma$) / tr($\Sigma_{crlb}$)'},
                               'ratioTraceSigma3' : { 'num' : eqns.traceSigmacrlb3NoTrace, 'denom' : eqns.traceSigmacrlb3, 'label' : r'Mode 4 tr($\Sigma$) / tr($\Sigma_{crlb}$)'},
                               'ratioTraceSigma4' : { 'num' : eqns.traceSigmacrlb4NoTrace, 'denom' : eqns.traceSigmacrlb4, 'label' : r'Mode 4 tr($\Sigma$) / tr($\Sigma_{crlb}$)'}
                             }
                  }

DETANDEVAL_RATIO = { 'subdir2' : 'ratio_notrace_over_full', 'log' : False,
                   'plots' : {
                               'ratioDetSigma2': { 'num' : eqns.DetSigmacrlb2NoTrace, 'denom' : eqns.DetSigmacrlb2, 'label' : r'Mode 2 Det($\Sigma_{notrace}$) /    Det($\Sigma_{crlb}$)'},
                               'ratioDetSigma3': { 'num' : eqns.DetSigmacrlb3NoTrace, 'denom' : eqns.DetSigmacrlb3, 'label' : r'Mode 3 Det($\Sigma_{notrace}$) / Det($\Sigma_{crlb}$)'},
                               'ratioDetSigma4': { 'num' : eqns.DetSigmacrlb4NoTrace, 'denom' : eqns.DetSigmacrlb4, 'label' : r'Mode 4 Det($\Sigma_{notrace}$) / Det($\Sigma_{crlb}$)'},
                               'ratioEval1Sigma2' : { 'num' : eqns.evalplusSigmacrlb2NoTrace, 'denom' : eqns.evalplusSigmacrlb2, 'label' : r'Mode 2 $\lambda_{+}$ / ${\lambda_{+}}_{crlb}$'},
                               'ratioEval2Sigma2' : { 'num' : eqns.evalminusSigmacrlb2NoTrace, 'denom' : eqns.evalminusSigmacrlb2, 'label' : r'Mode 2 $\lambda_{-}$ / ${\lambda_{-}}_{crlb}$'},
                               'ratioEval1Sigma3' : { 'num' : eqns.evalplusSigmacrlb3NoTrace, 'denom' : eqns.evalplusSigmacrlb3, 'label' : r'Mode 3 $\lambda_{+}$ / ${\lambda_{+}}_{crlb}$'},
                               'ratioEval2Sigma3' : { 'num' : eqns.evalminusSigmacrlb3NoTrace, 'denom' : eqns.evalminusSigmacrlb3, 'label' : r'Mode 3 $\lambda_{-}$ / ${\lambda_{-}}_{crlb}$'},
                               'ratioEval1Sigma4' : { 'num' : eqns.evalplusSigmacrlb4NoTrace, 'denom' : eqns.evalplusSigmacrlb4, 'label' : r'Mode 4 $\lambda_{+}$ / ${\lambda_{+}}_{crlb}$'},
                               'ratioEval2Sigma4' : { 'num' : eqns.evalminusSigmacrlb4NoTrace, 'denom' : eqns.evalminusSigmacrlb4, 'label' : r'Mode 4 $\lambda_{-}$ / ${\lambda_{-}}_{crlb}$'}
                             }
                  }


COV_RATIO = { 'subdir2' : 'ratio_notrace_over_full', 'log' : False,
                   'plots' : {
                               'ratioSigmacrlbCK2': { 'num' : eqns.SigmacrlbCK2NoTrace, 'denom' : eqns.SigmacrlbCK2, 'label' : r'Mode 2 $\langle{\delta c \delta k_{off}}\rangle$ / $\langle{\delta c \delta k_{off}}_{crlb}\rangle$'},
                               'ratioSigmacrlbCK3': { 'num' : eqns.SigmacrlbCK3NoTrace, 'denom' : eqns.SigmacrlbCK3, 'label' : r'Mode 2 $\langle{\delta c \delta k_{off}}\rangle$ / $\langle{\delta c \delta k_{off}}_{crlb}\rangle$'},
                               'ratioSigmacrlbCK4': { 'num' : eqns.SigmacrlbCK4NoTrace, 'denom' : eqns.SigmacrlbCK4, 'label' : r'Mode 2 $\langle{\delta c \delta k_{off}}\rangle$ / $\langle{\delta c \delta k_{off}}_{crlb}\rangle$'}
                             }
                  }


############################ No ratios No Trace ######################################

RELERRANDTRACE_NOTRACE = { 'subdir2' : 'notrace', 'log' : True,
                   'plots' : {
                               #'relErrX1': { 'eqn' : eqns.RelErrorX1NoTrace, 'label' : r'Mode 1 $\langle\delta x^{2}\rangle /x^{2}$'},
                               'relErrC2': { 'eqn' : eqns.RelErrC2NoTrace, 'label' : r'Mode 2 $\langle\delta c^{2}\rangle /c^{2}$'},
                               'relErrK2': { 'eqn' : eqns.RelErrK2NoTrace, 'label' : r'Mode 2 $\langle{\delta k_{off}^2}\rangle /k_{off}^{2}$'},
                               'relErrC3': { 'eqn' : eqns.RelErrC3NoTrace, 'label' : r'Mode 3 $\langle\delta c^{2}\rangle /c^{2}$'},
                               'relErrK3': { 'eqn' : eqns.RelErrK3NoTrace, 'label' : r'Mode 3 $\langle{\delta k_{off}^2}\rangle /k_{off}^{2}$'},
                               'relErrC4': { 'eqn' : eqns.RelErrC4NoTrace, 'label' : r'Mode 4 $\langle\delta c^{2}\rangle /c^{2}$'},
                               'relErrK4': { 'eqn' : eqns.RelErrK4NoTrace, 'label' : r'Mode 4 $\langle{\delta k_{off}^2}\rangle /k_{off}^{2}$'},
                               'traceSigma2': {'eqn' : eqns.traceSigmacrlb2NoTrace, 'label' : r'Mode 2 tr($\Sigma$)'},
                               'traceSigma3': {'eqn' : eqns.traceSigmacrlb3NoTrace, 'label' : r'Mode 3 tr($\Sigma$)'},
                               'traceSigma4': {'eqn' : eqns.traceSigmacrlb4NoTrace, 'label' : r'Mode 4 tr($\Sigma$)'}
                             }
                  }

DETANDEVAL_NOTRACE = { 'subdir2' : 'notrace', 'log' : False,
                   'plots' : {
                               'detSigma2': {'eqn' : eqns.DetSigmacrlb2NoTrace, 'label' : r'Mode 2 Det($\Sigma$)'},
                               'detSigma3': {'eqn' : eqns.DetSigmacrlb3NoTrace, 'label' : r'Mode 3 Det($\Sigma$)'},
                               'detSigma4': {'eqn' : eqns.DetSigmacrlb4NoTrace, 'label' : r'Mode 4 Det($\Sigma$)'},
                               'evalplusSigma2': {'eqn' : eqns.evalplusSigmacrlb2NoTrace, 'label' : r'Mode 2 $\lambda_{+}$'},
                               'evalminusSigma2': {'eqn' : eqns.evalminusSigmacrlb2NoTrace, 'label' : r'Mode 2 $\lambda_{-}$'},
                               'evalplusSigma3': {'eqn' : eqns.evalplusSigmacrlb3NoTrace, 'label' : r'Mode 3 $\lambda_{+}$'},
                               'evalminusSigma3': {'eqn' : eqns.evalminusSigmacrlb3NoTrace, 'label' : r'Mode 3 $\lambda_{-}$'},
                               'evalplusSigma4': {'eqn' : eqns.evalplusSigmacrlb4NoTrace, 'label' : r'Mode 4 $\lambda_{+}$'},
                               'evalminusSigma4': {'eqn' : eqns.evalminusSigmacrlb4NoTrace, 'label' : r'Mode 4 $\lambda_{-}$'},
                             }
                  }


COV_NOTRACE = { 'subdir2' : 'notrace', 'log' : False,
                   'plots' : {
                               'SigmaCK2': { 'eqn' : eqns.SigmacrlbCK2, 'label' : r'Mode 2 $\langle \delta c \delta k_{off} \rangle$'},
                               'SigmaCK3': { 'eqn' : eqns.SigmacrlbCK3, 'label' : r'Mode 3 $\langle \delta c \delta k_{off} \rangle$ '},
                               'SigmaCK4': { 'eqn' : eqns.SigmacrlbCK4, 'label' : r'Mode 4 $\langle \delta c \delta k_{off} \rangle$'}
                             }
                  }

############################ No ratios Full ######################################

DIAGANDTRACE_FULL = { 'subdir2' : 'notrace', 'log' : True,
                   'plots' : {
                               #'relErrX1': { 'eqn' : eqns.RelErrorX1, 'label' : r'Mode 1 $\langle\delta x^{2}\rangle /x^{2}$'},
                               'SigmacrlbC2': { 'eqn' : eqns.SigmacrlbC2, 'label' : r'Mode 2 $\langle \delta c^{2} \rangle$'},
                               'SigmacrlbK2': { 'eqn' : eqns.SigmacrlbK2, 'label' : r'Mode 2 $\langle \delta k_{off}^2 \rangle$'},
                               'SigmacrlbC3': { 'eqn' : eqns.SigmacrlbC3, 'label' : r'Mode 3 $\langle\delta c^{2}\rangle$'},
                               'SigmacrlbK3': { 'eqn' : eqns.SigmacrlbK3, 'label' : r'Mode 3 $\langle \delta k_{off}^2 \rangle$'},
                               'SigmacrlbC4': { 'eqn' : eqns.SigmacrlbC4, 'label' : r'Mode 4 $\langle \delta c^{2} \rangle$'},
                               'SigmacrlbK4': { 'eqn' : eqns.SigmacrlbK4, 'label' : r'Mode 4 $\langle \delta k_{off}^2 \rangle}$'},
                               'traceSigma2': {'eqn' : eqns.traceSigmacrlb2, 'label' : r'Mode 2 tr($\Sigma$)'},
                               'traceSigma3': {'eqn' : eqns.traceSigmacrlb3, 'label' : r'Mode 3 tr($\Sigma$)'},
                               'traceSigma4': {'eqn' : eqns.traceSigmacrlb4, 'label' : r'Mode 4 tr($\Sigma$)'}
                             }
                  }

DETANDEVAL_FULL = { 'subdir2' : 'notrace', 'log' : False,
                   'plots' : {
                               'detSigma2': {'eqn' : eqns.DetSigmacrlb2, 'label' : r'Mode 2 Det($\Sigma$'},
                               'detSigma3': {'eqn' : eqns.DetSigmacrlb3, 'label' : r'Mode 3 Det($\Sigma$'},
                               'detSigma4': {'eqn' : eqns.DetSigmacrlb4, 'label' : r'Mode 4 det($\Sigma$'},
                               'evalplusSigma2': {'eqn' : eqns.evalplusSigmacrlb2, 'label' : r'Mode 2 $\lambda_{+}$'},
                               'evalminusSigma2': {'eqn' : eqns.evalminusSigmacrlb2, 'label' : r'Mode 2 $\lambda_{-}$'},
                               'evalplusSigma3': {'eqn' : eqns.evalplusSigmacrlb3, 'label' : r'Mode 3 $\lambda_{+}$'},
                               'evalminusSigma3': {'eqn' : eqns.evalminusSigmacrlb3, 'label' : r'Mode 3 $\lambda_{-}$'},
                               'evalplusSigma4': {'eqn' : eqns.evalplusSigmacrlb4, 'label' : r'Mode 4 $\lambda_{+}$'},
                               'evalminusSigma4': {'eqn' : eqns.evalminusSigmacrlb4, 'label' : r'Mode 4 $\lambda_{-}$'},
                             }
                  }


COV_FULL = { 'subdir2' : 'notrace', 'log' : False,
                   'plots' : {
                               'SigmaCK2': { 'eqn' : eqns.SigmacrlbCK2, 'label' : r'Mode 2 \langle \delta c \delta k_{off} \rangle$'},
                               'SigmaCK3': { 'eqn' : eqns.SigmacrlbCK3, 'label' : r'Mode 3 \langle \delta c \delta k_{off} \rangle$ '},
                               'SigmaCK4': { 'eqn' : eqns.SigmacrlbCK4, 'label' : r'Mode 4 \langle \delta c \delta k_{off} \rangle$'}
                             }
                  }
