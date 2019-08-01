#!/usr/bin/env python
import numpy as np
import scipy
from math_functions import *
import matplotlib.pyplot as plt
import sys
import random


__doc__ = '''Extract chi2 values along the parameter space.
The file stem <fstem> is passed as an argument.
<fstem>.dat contains the data, and <fstem>_amias.dat contains the fitted parameters.
For each parameter, there is a line in the parameter file with that parameter's central value and randomization range.
Each parameter is randomized uniformly within that range. Potentially add a switch within the parameter file to randomize eg from a gaussian.
'''


def model_potential(pars, coordinates, n_order, m_order):
    # model function for the potential at a point with coordinates [mu,eta,phi]
    mu = coordinates[:,0]
    eta = coordinates[:,1]
    phi = coordinates[:,2]

    f_mu_eta = np.sqrt(np.cosh(mu) - np.cos(eta))
    model_V = np.zeros(len(mu)) # each element is a multipole sum over n,m at a point
    for i in range(len(model_V)): 
        # calculate each element separately
        # the Legendre functions cannot work with arrays
        mu_i=mu[i]
        eta_i=eta[i]
        phi_i=phi[i]
        for n in range(n_order+1):
            for m in range(n+1): # for now m goes up to n, or up to m_order
                if m>m_order: continue
                nm_id = 4*(sum(range(n+1)) + m) #n-m position in the parameter array
                A_nm=pars[nm_id+0]
                B_nm=pars[nm_id+1]
                eta_nm=pars[nm_id+2]
                phi_nm=pars[nm_id+3]
                model_V[i] += (A_nm*legendre_Q_toroidal(m,n,mu_i) + B_nm*legendre_P_toroidal(m,n,mu_i)) * np.cos(n*(eta_i-eta_nm)) * np.cos(m*(phi_i-phi_nm))
        model_V[i] *= f_mu_eta[i]    
    
    return model_V


def residual(pars, coordinates, V, err_V, n_order, m_order):
    
    model_V = model_potential(pars,coordinates,n_order,m_order)

    norm_residuals = (V - model_V)/err_V
    ndf = len(mu)-len(pars)
    chi2 = sum(np.power(norm_residuals,2.)) / ndf
    
    return chi2



fstem = sys.argv[1] # the stem of the data file

ntrials = 1000
if len(sys.argv)>2:
    ntrials = int(sys.argv[2]) # the stem of the data file

data = np.loadtxt(fstem+ '.dat', comments='!')
x=data[:,0]
y=data[:,1]
z=data[:,2]
r=data[:,3]
# phi=data[:,4] -- unnecessary because not in [0,2pi] - get from y/x instead
# phi = np.pi + np.arctan2(y,x) # this phi will be in [0,2pi]
# this has been created in the file already, just read

V=data[:,8]
err_V=1 ### To do: realistic error estimate for each V value ###

mu=data[:,9]
eta=data[:,10]
phi=data[:,11]
## Read the toroidal coordinates in multi-D array
coordinates = data[:,9:12]

# distance of each "pole" in the bipolar coordinate system from origin.
# defined in cm, same as coordinates in file.
a=5.756 

# Read in fit parameters and their randomization range
pars_data = np.loadtxt(fstem + '_amias.dat', comments='!')
pars=pars_data[:,0]
pars_range=pars_data[:,1]
par_low,par_high = pars-pars_range/2, pars+pars_range/2

n_order,m_order=0,0
nm_order=len(pars)
for i in range(20):
    if 4*(sum(range(i+1))+i+1)==nm_order:
        n_order=i
        m_order=i
        break
print 'n,m order: ', n_order


### Iterate ntrials with randomized pars within their range, and output to data file
rand_pars = np.zeros(len(pars))
for iter in range(ntrials):
    print iter,

    # Create random set of parameters within range
    for i in range(len(pars)):
        rand_pars[i]=random.uniform(par_low[i],par_high[i])
        print rand_pars[i],

    chi2 = residual(rand_pars, coordinates, V, err_V, n_order, m_order)
    print chi2
