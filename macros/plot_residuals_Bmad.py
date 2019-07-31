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
    # Model function for the potential at a point with *cartesian* coordinates (x,y), with origin at the center of the storage region.
    # pars is an array of the a,b parameters up to order n_order
    x = np.array(coordinates[:,0])
    y = np.array(coordinates[:,1])
    r = np.array(x**2 + y**2)

    model_V = np.zeros(len(r)) # each element is a multipole sum over n orders
    for i in range(len(model_V)):
        xi=x[i]
        yi=y[i]
        ri=r[i]
        model_Vi=0
        for n in range(n_order+1):
            an=pars[2*n]
            bn=pars[2*n+1]
            zn=complex(bn,-an)
            xn=complex(xi,yi)
            model_V_n= zn * np.power(xn,n+1) / ((n+1) * np.power(ri,n))
            model_Vi -= model_V_n.real
        model_V[i] = model_Vi
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
r=np.array(data[:,0])
th=np.array(data[:,1])
z=np.array(data[:,2])
V=np.array(data[:,3])

err_V=1 ### To do: realistic error estimate for each V value ###


## Construct polar coordinates with origin at center of storage region
x = np.array(r-711.2) # in cm
th_polar = np.array(np.arctan2(z,x))
r_polar = np.array(x/np.cos(th_polar))

coordinates = np.vstack((x,z)).T
print '\n\nCoordinates:\n', coordinates


## Read in max order in multipole sum from argument
if len(sys.argv)>2:# 2nd arg given is for max n_order
    n_order=int(sys.argv[2])
else:
    n_order=5


## Read in fit parameters
pars = np.zeros(2*(n_order+1))
filename='pars_'+fstem+('_Bmad_%d.dat' % n_order)
exists = os.path.isfile(filename)
if exists:
    pars = np.loadtxt(filename, comments='!')    
    print '\n\n Fit parameters ::\n', pars
else:
    print '\n\n No fit result found\n'
    exit()


## Calculate model value and residuals based on fit parameters
model_V = np.array(model_potential(pars, coordinates, n_order))
residuals = np.array(residual(pars, coordinates, V, err_V, n_order))


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



