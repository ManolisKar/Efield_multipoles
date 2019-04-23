import numpy as np
import math
import scipy.integrate as integrate
from scipy.special import gamma,hyp2f1
from scipy.optimize import leastsq
#import gmpy2
#from bigfloat import BigFloat, precision

#np.seterr(all='ignore')

def olvers(a,b,c,x):
    # Olver's hypergeometric function 
    # See https://dlmf.nist.gov/15.1
    # Note that hyp2f1 does not exist for c=0,-1,-2,...
    if np.abs(c%(-1))<0.1:
        print ('Oops, too close to a pole!')
        # maybe build alternate definition here
    return hyp2f1(a,b,c,x) / gamma(c)

def legendre_P(m,n,x):
    # Legendre/Ferrers function of first kind
    # m,n can be real (non-integers)
    # See https://dlmf.nist.gov/14.3
    if m==0 and n==0:
        value = 1
    elif m==0 and n==1:
        value = x
    else:
        value = ((1+x)/(1-x))**(m/2) * olvers(n+1,-n,1-m,0.5-0.5*x)
    return value


def legendre_Q(m,n,x):
    # Legendre/Ferrers function of second kind
    # m,n can be real (non-integers)
    # See https://dlmf.nist.gov/14.3
    value  = ((1+x)/(1-x))**(m/2) * cos(m*np.pi) * olvers(n+1,-n,1-m,0.5-0.5*x)
    value -= (gamma(n+m+1)/gamma(n-m+1)) * ((1-x)/(1+x))**(m/2) * olvers(n+1,-n,1+m,0.5-0.5*x)
    value *= np.pi/(2*np.sin(m*np.pi))
    return value


def legendre_P_toroidal_integrand(phi, m,n,ksi):
    integrand  = (np.sin(phi))**(2*m) / (np.cosh(ksi)+np.cos(phi)*np.sinh(ksi))**(n+m+0.5)
    return integrand

def legendre_P_toroidal(m,n,ksi):
    # Toroidal Legendre function of first kind
    # Assumes that the order n is really half-integer n-1/2
    # Assumes that the argument ksi is really cosh(ksi)
    # See https://dlmf.nist.gov/14.19(iii)
    value  = gamma(n+m+0.5)*(np.sinh(ksi))**m
    value /= (2**m * np.sqrt(np.pi) * gamma(n-m+0.5) * gamma(m+0.5))
    integral, err = integrate.quad(legendre_P_toroidal_integrand, 0, np.pi, args=(m,n,ksi))
    #print ('Integration error : ', err)
    value *= integral
    return value


def legendre_Q_toroidal_integrand(t, m,n,ksi):
    # this integration goes to infinity, so need bigfloat to handle big numbers in the integrand
    # with precision(25600):
        #integrand  = BigFloat( coshmt / (np.cosh(ksi)+cosht*np.sinh(ksi))**(n+0.5) )
        # eh, this doesn't work either
    if m>n:
        print( 'm>n --> Bad condition! Need m<=n\n')
        exit()
    integrand  = np.cosh(m*t) / (np.cosh(ksi)+np.cosh(t)*np.sinh(ksi))**(n+0.5)
    #print('t=%f , integrand=%f' % (t, integrand))
    if math.isnan(integrand):
        #print('Potential divergence: (m,n,t) = (%d,%d,%f)' % (m,n,t))
        integrand = 0
    return integrand

def legendre_Q_toroidal(m,n,ksi):
    # Toroidal Legendre function of second kind
    # Assumes that the order n is really half-integer n-1/2
    # Assumes that the argument ksi is really cosh(ksi)
    # See https://dlmf.nist.gov/14.19(iii)
    value  = gamma(n+0.5)
    value /= (gamma(n+m+0.5) * gamma(n-m+0.5))
    integral, err = integrate.quad(legendre_Q_toroidal_integrand, 0, np.inf, args=(m,n,ksi))
    value *= integral
    return value
