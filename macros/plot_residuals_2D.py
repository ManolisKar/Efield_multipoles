#!/usr/bin/env python
import numpy as np
import scipy
from math_functions import *
import matplotlib.pyplot as plt
import sys
import os


__doc__ = '''Plot residuals from a multiple fit.
The data file <fstem>.dat is passed as an argument.
The parameters (the fit result) are expected to be in a file pars_<fstem>.out, 
and obeying a specific format: Comment lines beginning with '!', '[]' characters removed.

argv[1]: file stem
argv[2]: n_order
'''



def model_potential(pars, coordinates, n_order):
    # Model function for the potential at a point with *polar* coordinates (r,theta), with origin at the center of the storage region.
    # pars is an array of the a,b parameters up to order n_order
    r_polar = np.array(coordinates[:,0])
    th_polar = np.array(coordinates[:,1])

    model_V = np.zeros(len(r)) # each element is a multipole sum over n orders
    for i in range(len(model_V)): 
        for n in range(n_order+1):
            an=pars[2*n]
            bn=pars[2*n+1]
            model_V[i] += r_polar[i]**n *(an*np.sin(n*th_polar[i]) + bn*np.cos(n*th_polar[i]))
    return model_V 


def residual(pars, coordinates, V, err_V, n_order):
    model_V = model_potential(pars, coordinates, n_order)
    norm_residuals = (V - model_V)/err_V
    sum_residuals = sum(np.power(norm_residuals,2.))

    print 'pars ::\n', pars
    print 'normalized residual sum = ', sum_residuals
    ndf = len(model_V) - len(pars)
    print 'reduced chi2 = ', sum_residuals/ndf
    sys.stdout.flush()

    return norm_residuals


## Read in file stem and data
if len(sys.argv)>1:# 1st arg given is for file stem
    fstem=str(sys.argv[1])
else:
    fstem='test3000'
print 'Processing file '+fstem+'.dat'
data = np.loadtxt(fstem+'.dat', comments='!')
r=data[:,0]
th=data[:,1]
z=data[:,2]
V=data[:,3]

err_V=1 ### To do: realistic error estimate for each V value ###


## Construct polar coordinates with origin at center of storage region
x = r-711.2 # in cm
th_polar = np.array(np.arctan2(z,x))
r_polar = np.array(x/np.cos(th_polar))

coordinates = np.vstack((r_polar,th_polar)).T
print '\n\nCoordinates:\n', coordinates


## Read in max order in multipole sum from argument
if len(sys.argv)>2:# 2nd arg given is for max n_order
    n_order=int(sys.argv[2])
else:
    n_order=5


## Read in fit parameters
pars = np.zeros(2*(n_order+1))
exists = os.path.isfile('pars_'+fstem+('_%d.dat' % n_order))
if exists:
    pars = np.loadtxt('pars_'+fstem+('_%d.dat' % n_order), comments='!')    
    print '\n\n Fit parameters ::\n', pars
else:
    print '\n\n No fit result found\n'
    exit()


## Calculate model value and residuals based on fit parameters
model_V = np.array(model_potential(pars, coordinates, n_order))
residuals = np.array(residual(pars, coordinates, V, err_V, n_order))

r=np.array(r)
z=np.array(z)
V=np.array(V)


# Plot V in color code vs r,z
fig=plt.figure('V_color', figsize=(18, 7))
plt.scatter(
    r,z,c=V, s=20, alpha=0.5, cmap='coolwarm'
)
plt.title("Voltage data")
plt.xlabel('r [cm]')
plt.ylabel('z [cm]')
plt.colorbar()
plt.show()
fig.savefig('V_color_vs_r_th.png')


# Plot function with fit pars in color code vs r,z
fig=plt.figure('fit_color', figsize=(18, 7))
plt.scatter(
    r,z,c=model_V, s=20, alpha=0.5,cmap='coolwarm'
)
plt.title("Fit result")
plt.xlabel('r [cm]')
plt.ylabel('z [cm]')
plt.colorbar()
plt.show()
fig.savefig('fit_color_vs_r_th.png')


# Plot residuals in color code vs r,z
fig=plt.figure('res_color', figsize=(18, 7))
plt.scatter(
    r,z,c=residuals, s=20, alpha=0.5, cmap='coolwarm'
)
plt.title("Fit residuals")
plt.xlabel('r [cm]')
plt.ylabel('z [cm]')
plt.colorbar()
plt.show()
fig.savefig('res_color_vs_r_th.png')



