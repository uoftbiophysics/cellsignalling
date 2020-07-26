import numpy as np


N = 20000
kp = 10.0; t = 1000.0; koff = 1.0; kon = 1.0;  alpha = 1.0

def prior(alpha, c):
    return alpha*np.exp(-alpha*c)

def mean_n(c, kon, koff, kp, t):
    return ( kp*t*c*kon )/( koff*(1+kon*c/koff) )

def variance_n(c, kon, koff, kp, t):
    return ( kp*t*c*kon )/( koff*( 1 + kon*c/koff ) ) + ( 2*kp**2*t*kon*c ) /( koff**2*( 1+c*kon/koff )**3 )

def likelihood(n, c, koff, kon, kp, t):
    return np.exp(-(n-mean_n(c, kon, koff, kp, t))**2/(2*variance_n(c, kon, koff, kp, t) ) )/np.sqrt( 2*np.pi*variance_n(c, kon, koff, kp, t) )

def mutualInformation(n_vec, c_vec, kon, koff, kp, t, dn, dc, alpha):
    I1 = 0.; I2 = 0.
    i = 0;
    for n in n_vec:
        I1 += np.sum(prior(alpha, c_vec)*likelihood(n, c_vec, koff, kon, kp, t)*np.ma.log2(likelihood(n, c_vec, koff, kon, kp, t))*dc)*dn
        I2 += np.sum(prior(alpha, c_vec)*likelihood(n, c_vec, koff, kon, kp, t)*dc)*np.ma.log2(np.sum(prior(alpha, c_vec)*likelihood(n, c_vec, koff, kon, kp, t)*dc))*dn
        if n % 1000 == 0.:
            print(n)
    print(np.sum(I1),np.sum(I2))
    return np.sum(I1) - np.sum(I2)

n, dn = np.linspace(0., N, N+1, retstep=True);
c, dc = np.linspace(0.001, 10., N+1, retstep=True);
"""
c = np.logspace(0.001, 10., N+1); dc = np.zeros(N+1)
for i, c_i in enumerate(c):
    if i < N:
        dc[i] = c[i+1]-c[i]
    else:
        dc[i] = c[i]-c[i-1]
"""

print(mutualInformation(n, c, kon, koff, kp, t, dn, dc, alpha))
#print(dc, dn)
