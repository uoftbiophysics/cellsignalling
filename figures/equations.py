import numpy as np
from numpy import exp as exp
import os
from settings import KON, KP, T, KF, ALPHA

def Sqrt(x):
    return np.sqrt(x)

def DetSigmacrlb2(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return (koff**3*kon**2*KP*T**2)/((koff + c*kon)**2*(koff**2 + c**2*kon**2 + koff*KP)**2) + (c**2*koff*kon**4*KP*T**2)/((koff + c*kon)**2*(koff**2 + c**2*kon**2 + koff*KP)**2) + (koff**2*kon**2*KP**2*T**2)/((koff + c*kon)**2*(koff**2 + c**2*kon**2 + koff*KP)**2)

def DetSigmacrlb2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return 0. + (0.25*koff**10)/(c**2*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (0.5*koff**9*kon)/(c*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (1.25*koff**8*kon**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.*c*koff**7*kon**3)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.5*c**2*koff**6*kon**4)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (3.*c**3*koff**5*kon**5)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.5*c**4*koff**4*kon**6)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.*c**5*koff**3*kon**7)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (1.25*c**6*koff**2*kon**8)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (0.5*c**7*koff*kon**9)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (0.25*c**8*kon**10)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.5*koff**9*KP)/(c**2*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (3.*koff**8*kon*KP)/(c*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (16.5*koff**7*kon**2*KP)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (14.*c*koff**6*kon**3*KP)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (35.*c**2*koff**5*kon**4*KP)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (24.*c**3*koff**4*kon**5*KP)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (31.*c**4*koff**3*kon**6*KP)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (18.*c**5*koff**2*kon**7*KP)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (10.5*c**6*koff*kon**8*KP)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (5.*c**7*kon**9*KP)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (0.5*c**8*kon**10*KP)/(koff*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.25*koff**8*KP**2)/(c**2*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (11.5*koff**7*kon*KP**2)/(c*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (37.25*koff**6*kon**2*KP**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (38.*c*koff**5*kon**3*KP**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (58.25*c**2*koff**4*kon**4*KP**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (41.5*c**3*koff**3*kon**5*KP**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (29.75*c**4*koff**2*kon**6*KP**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (15.*c**5*koff*kon**7*KP**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.5*c**6*kon**8*KP**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.*koff**7*KP**3)/(c**2*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (12.*koff**6*kon*KP**3)/(c*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (30.*koff**5*kon**2*KP**3)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (26.*c*koff**4*kon**3*KP**3)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (28.*c**2*koff**3*kon**4*KP**3)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (14.*c**3*koff**2*kon**5*KP**3)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (4.*c**4*koff*kon**6*KP**3)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.*koff**6*KP**4)/(c**2*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (4.*koff**5*kon*KP**4)/(c*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (8.*koff**4*kon**2*KP**4)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (4.*c*koff**3*kon**3*KP**4)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.*c**2*koff**2*kon**4*KP**4)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (0.5*koff**10*kon*T)/(c*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (1.5*koff**9*kon**2*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (5.*c*koff**8*kon**3*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (7.*c**2*koff**7*kon**4*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (4.*c**3*koff**6*kon**5*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (4.*c**4*koff**5*kon**6*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (3.*c**5*koff**4*kon**7*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (c**6*koff**3*kon**8*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (3.5*c**7*koff**2*kon**9*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.5*c**8*koff*kon**10*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (5.*koff**9*kon*KP*T)/(c*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (4.*koff**8*kon**2*KP*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (23.5*c*koff**7*kon**3*KP*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (18.5*c**2*koff**6*kon**4*KP*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (17.5*c**3*koff**5*kon**5*KP*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (10.5*c**4*koff**4*kon**6*KP*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (0.5*c**5*koff**3*kon**7*KP*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (5.5*c**6*koff**2*kon**8*KP*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (6.5*c**7*koff*kon**9*KP*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (1.5*c**8*kon**10*KP*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (12.5*koff**8*kon*KP**2*T)/(c*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (3.5*koff**7*kon**2*KP**2*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (32.5*c*koff**6*kon**3*KP**2*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (17.5*c**2*koff**5*kon**4*KP**2*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (7.5*c**3*koff**4*kon**5*KP**2*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (12.5*c**4*koff**3*kon**6*KP**2*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (12.5*c**5*koff**2*kon**7*KP**2*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.5*c**6*koff*kon**8*KP**2*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (12.*koff**7*kon*KP**3*T)/(c*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (koff**6*kon**2*KP**3*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (14.*c*koff**5*kon**3*KP**3*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (10.*c**2*koff**4*kon**4*KP**3*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) - (6.*c**3*koff**3*kon**5*KP**3*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (9.*c**4*koff**2*kon**6*KP**3*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (4.*koff**6*kon*KP**4*T)/(c*(koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (4.*c**2*koff**3*kon**4*KP**4*T)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (koff**9*kon**2*KP*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.*c*koff**8*kon**3*KP*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (4.*c**2*koff**7*kon**4*KP*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.*c**3*koff**6*kon**5*KP*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.*c**4*koff**5*kon**6*KP*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.*c**5*koff**4*kon**7*KP*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (4.*c**6*koff**3*kon**8*KP*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.*c**7*koff**2*kon**9*KP*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (c**8*koff*kon**10*KP*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (3.*koff**8*kon**2*KP**2*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.*c*koff**7*kon**3*KP**2*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (9.*c**2*koff**6*kon**4*KP**2*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (12.*c**3*koff**5*kon**5*KP**2*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (9.*c**4*koff**4*kon**6*KP**2*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.*c**5*koff**3*kon**7*KP**2*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (3.*c**6*koff**2*kon**8*KP**2*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (3.*koff**7*kon**2*KP**3*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.*c*koff**6*kon**3*KP**3*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.*c**2*koff**5*kon**4*KP**3*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (6.*c**3*koff**4*kon**5*KP**3*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (3.*c**4*koff**3*kon**6*KP**3*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (koff**6*kon**2*KP**4*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (2.*c*koff**5*kon**3*KP**4*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4) + (c**2*koff**4*kon**4*KP**4*T**2)/((koff + c*kon)**4*(koff**2 + c**2*kon**2 + koff*KP)**4)

def RelErrC2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return ((koff + c*kon)*(c**2*kon**2 + koff*KP))/(c*koff**2*kon*KP*T)

def RelErrK2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return ((koff + c*kon)*(koff + KP))/(c*koff*kon*KP*T)

def RelErrorX1NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return ((koff + c*KON)*(koff**2 + c**2*KON**2 + 2*koff*(c*KON + KP)))/(c*koff**2*KON*KP*T)

def SigmacrlbC2(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return (koff**2*kon*(koff + KP)*T)/(c*(koff + c*kon)*(c**2*kon**2 + koff*(koff + KP)))

def SigmacrlbC2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return (koff*(c**4*kon**4*KP + koff**5*(1. + c*kon*T) + c**2*koff*kon**2*(-6.*c*kon*KP + 3.*KP**2 + c**2*kon**2*(5. + KP*T)) + c*koff**2*kon*(-2.*KP**2 + c**3*kon**3*T + c**2*kon**2*(-6. + KP*T) + c*kon*KP*(9. + KP*T)) + koff**3*(KP**2 + c**3*kon**3*T + c*kon*KP*(-4. + KP*T) + c**2*kon**2*(6. + 2.*KP*T)) + koff**4*(2.*KP + c*kon*(-2. + c*kon*T + 2.*KP*T))))/(c**2*(koff + c*kon)**2*(c**2*kon**2 + koff*(koff + KP))**2)

def SigmacrlbK2(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return ((c**3*kon**3 + c*koff*kon*KP)*T)/(koff*(koff + c*kon)*(c**2*kon**2 + koff*(koff + KP)))

def SigmacrlbK2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return (0.5*koff**6 + 0.5*c**6*kon**6 + koff**5*KP*(3. + c*kon*T) + c**4*koff*kon**4*(KP + c*kon*(-2. + c*kon*T)) + koff**4*(3.*KP**2 + c**3*kon**3*T + c*kon*KP*(-3. + KP*T) + c**2*kon**2*(5.5 + KP*T)) + c*koff**3*kon*(-2.*KP**2 + c**3*kon**3*T + c*kon*KP*(10. + KP*T) + c**2*kon**2*(-6. + 2.*KP*T)) + c**2*koff**2*kon**2*(-5.*c*kon*KP + KP**2 + c**3*kon**3*T + c**2*kon**2*(5.5 + 2.*KP*T)))/(koff**2*(koff + c*kon)**2*(c**2*kon**2 + koff*(koff + KP))**2)

def RelDetSigmacrlb2(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return ( DetSigmacrlb2(c, koff, kon, T, KF, KP) )/( c**2 * koff**2 )

def RelDetSigmacrlb2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return ( DetSigmacrlb2NoTrace(c, koff, kon, T, KF, KP) )/( c**2 * koff**2 )

def traceSigmacrlb2(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return (SigmacrlbC2(c, koff, kon, T, KF, KP)+SigmacrlbK2(c, koff, kon, T, KF, KP))

def evalplusSigmacrlb2(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    tr = traceSigmacrlb2(c, koff, kon, T, KF, KP); det = DetSigmacrlb2(c, koff, kon, T, KF, KP);
    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )

def evalminusSigmacrlb2(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    tr = traceSigmacrlb2(c, koff, kon, T, KF, KP); det = DetSigmacrlb2(c, koff, kon, T, KF, KP);
    return ( 0.5*tr - 0.5*np.sqrt( tr**2-4*det ) )

def dedimRelErrC2(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return KP*T*RelErrC2(c, koff, kon, T, KF, KP)

def dedimRelErrK2(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return KP*T*RelErrK2(c, koff, kon, T, KF, KP)

def traceSigmacrlb2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return (SigmacrlbC2NoTrace(c, koff, kon, T, KF, KP)+SigmacrlbK2NoTrace(c, koff, kon, T, KF, KP))

def evalplusSigmacrlb2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    tr = traceSigmacrlb2NoTrace(c, koff, kon, T, KF, KP); det = DetSigmacrlb2NoTrace(c, koff, kon, T, KF, KP);
    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )

def evalminusSigmacrlb2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    tr = traceSigmacrlb2NoTrace(c, koff, kon, T, KF, KP); det = DetSigmacrlb2NoTrace(c, koff, kon, T, KF, KP);
    return ( 0.5*tr - 0.5*np.sqrt( tr**2-4*det ) )

def dedimRelErrC2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return KP*T*RelErrC2NoTrace(c, koff, kon, T, KF, KP)

def dedimRelErrK2NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return KP*T*RelErrK2NoTrace(c, koff, kon, T, KF, KP)

def dedimRelErrorX1NoTrace(c, koff, kon=KON, T=T, KF=KF, KP=KP):
    return KP*T*RelErrorX1NoTrace(c, koff, kon, T, KF, KP)