#!/usr/bin/env python
import numpy as np
import scipy
from math_functions import *
import matplotlib.pyplot as plt
import sys
import os

plot_results = 0

iter_count=0

__doc__ = '''Fit data of potential from an OPERA map file
to a function of 2D multipoles.
OPERA map file expected in r,z space

1st arg: data file stem
2nd arg: n_order for V model function
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
    global iter_count

    model_V = model_potential(pars, coordinates, n_order)
    norm_residuals = (V - model_V)/err_V
    sum_residuals = sum(np.power(norm_residuals,2.))

    iter_count += 1
    print 'Iteration #', iter_count
    print 'pars ::\n', pars
    print 'normalized residual sum = ', sum_residuals
    ndf = len(model_V) - len(pars)
    print 'reduced chi2 = ', sum_residuals/ndf
    sys.stdout.flush()

    return norm_residuals


#
# Open OPERA map file and read in data.
# In the 2D script, this file only contains r,th,z,V data in a thin slice, 
# but in *cylindrical* coordinates, 
# and with origin at center of ring.
# Need to convert to polar with origin at center of storage region.
#

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
# phi=data[:,4] -- unnecessary because not in [0,2pi] - get from y/x instead
# phi = np.pi + np.arctan2(y,x) # this phi will be in [0,2pi]
# this has been created in the file already, just read

# Construct polar coordinates with origin at center of storage region
x = r-711.2 # in cm
r_polar = np.array(x/np.cos(th))
th_polar = np.array(np.pi + np.arctan2(z,x))

coordinates = np.vstack((r_polar,th_polar)).T
print '\n\nCoordinates:\n', coordinates

err_V=1 ### To do: realistic error estimate for each V value ###

if len(sys.argv)>2:# 2nd arg given is for max n_order
    n_order=int(sys.argv[2])
else:
    n_order=5

# I believe pars is required to be 1-D array of parameters to be optimized.
# pars[2*n] are the A_n, pars[2*n+1] are the B_n.
pars = np.zeros(2*(n_order+1))
# Get initial parameter values from file - if it exists
exists = os.path.isfile('pars_'+fstem+('_%d.dat' % n_order))
if exists:
    pars = np.loadtxt('pars_'+fstem+('_%d.dat' % n_order), comments='!')    
print '\n\nInitial parameters ::\n', pars

# Calculate and minimize the residuals
residual(pars,coordinates,V,err_V,n_order)
fit_result = scipy.optimize.least_squares(residual, pars, args=(coordinates, V, err_V, n_order),ftol=1e-15,xtol=1e-12)
#result_pars, flag = scipy.optimize.least_squares(residual, pars, args=(coordinates, V, err_V, n_order),ftol=1e-15,xtol=1e-12)
print ('Result status = %d\n Fit parameters ::\n'%fit_result.status), fit_result.x

# Output parameters to file
pars_file = open('pars_'+fstem+('_%d.dat' % n_order), 'w+')
np.savetxt(pars_file, fit_result.x)
pars_file.close()

if not plot_results:
    exit()


# Plot results
r=np.array(r)
z=np.array(z)
V=np.array(V)
final_residuals = np.array(residual(result_pars,coordinates,V,err_V,n_order))

# Plot V in color code vs r,z
fig=plt.figure('V_color', figsize=(18, 7))
plt.scatter(
    r,z,c=V, cmap='afmhot'
)
plt.xlabel('r [cm]')
plt.ylabel('z [cm]')
plt.colorbar()
plt.show()
fig.savefig('V_color_vs_r_th.png')

'''
V_color_plot = plt.figure(figsize=(18,7)).add_subplot(1,1,1,
    title='V vs r,z',xlabel='r',ylabel='z'
)
V_color_plot.scatter(
    r,z,c=V, cmap='afmhot'
)
#V_color_plot.figure.savefig('V_color_vs_r_th.png')
#plt.colorbar(V_color_plot)
V_color_plot.figure.show()
'''

# Plot residuals in r,z color code
fig=plt.figure('residual_color', figsize=(18, 7))
plt.scatter(
    r,z,c=final_residuals, cmap='afmhot'
)
plt.xlabel('r [cm]')
plt.ylabel('z [cm]')
plt.colorbar()
plt.show()
fig.savefig('residual_color_vs_r_th.png')
