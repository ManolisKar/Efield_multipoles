#!/usr/bin/env python
import numpy as np
import scipy
from math_functions import *
import matplotlib.pyplot as plt
import sys


__doc__ = '''Plot residuals from a multiple fit.
The data file <fstem>.dat is passed as an argument.
The parameters (the fit result) are expected to be in a file pars_<fstem>.out, 
and obeying a specific format: Comment lines beginning with '!', '[]' characters removed.
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
    sum_residuals = sum(np.power(norm_residuals,2.))
    print 'normalized residual sum = ', sum_residuals

    return norm_residuals


fstem = sys.argv[1] # the stem of the data file

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

## Read the toroidal coordinates in multi-D array
coordinates = data[:,9:12]

# distance of each "pole" in the bipolar coordinate system from origin.
# defined in cm, same as coordinates in file.
a=5.756 

# Now read in the result for the fit parameters
pars = np.loadtxt('pars_' + fstem + '.dat', comments='!')
nm_order=len(pars)
for i in range(20):
    if 4*(sum(range(i+1))+i+1)==nm_order:
        n_order=i
        m_order=i
        break

print 'n,m order: ', n_order

model_V = np.array(model_potential(pars, coordinates, n_order, m_order))
residuals = np.array(residual(pars, coordinates, V, err_V, n_order, m_order))

mu=np.array(data[:,9])
eta=np.array(data[:,10])
phi=np.array(data[:,11])
r=np.array(data[:,3])
V=np.array(V)


# Plot V in color code vs r,z
fig=plt.figure('V_color', figsize=(18, 7))
plt.scatter(
    r,z,c=V, cmap='coolwarm'
)
plt.xlabel('r [cm]')
plt.ylabel('z [cm]')
plt.colorbar()
plt.show()
fig.savefig('V_color_vs_r_th.png')


# Plot function with fit pars in color code vs r,z
fig=plt.figure('fit_color', figsize=(18, 7))
plt.scatter(
    r,z,c=model_V, cmap='coolwarm'
)
plt.xlabel('r [cm]')
plt.ylabel('z [cm]')
plt.colorbar()
plt.show()
fig.savefig('fit_color_vs_r_th.png')


# Plot residuals in color code vs r,z
fig=plt.figure('res_color', figsize=(18, 7))
plt.scatter(
    r,z,c=residuals, cmap='coolwarm'
)
plt.xlabel('r [cm]')
plt.ylabel('z [cm]')
plt.colorbar()
plt.show()
fig.savefig('res_color_vs_r_th.png')



