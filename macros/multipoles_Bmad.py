#!/usr/bin/env python
import numpy as np
import scipy
from math_functions import *
import matplotlib.pyplot as plt
import sys
import os
import time

start = time.time()



## Switches
plot_results = 0
debug = 0

## Global parameters
iter_count=0
best_chi2 = 1e10
best_pars = np.zeros(100)
found_better_solution=0


__doc__ = '''Fit data of potential from an OPERA map file
to a function of 2D multipoles.
The function in this script is taken from the BMad manual eqn.(15.18).
OPERA map file expected in r,z space

1st arg: data file stem
2nd arg: n_order for V model function
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
    global iter_count
    global best_chi2
    global best_pars
    global found_better_solution

    model_V = model_potential(pars, coordinates, n_order)
    norm_residuals = (V - model_V)/err_V
    sum_residuals = sum(np.power(norm_residuals,2.))
    
    ndf = len(model_V) - len(pars)
    reduced_chi2 = sum_residuals/ndf
    if reduced_chi2<best_chi2:
        best_chi2=reduced_chi2
        best_pars=pars
        found_better_solution=1

    if iter_count%100==0:
        print '\n\nIteration #', iter_count
        if found_better_solution:
            print 'Best solution::\n', best_pars
            print '\nnmax   reduced-chi2    max-residual::'
            print n_order, best_chi2, max(norm_residuals)
            print 'Time elapsed :: ', time.time()-start
            sys.stdout.flush()
            found_better_solution=0

    iter_count += 1
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
r=np.array(data[:,0])
th=np.array(data[:,1])
z=np.array(data[:,2])
V=np.array(data[:,3])
# phi=data[:,4] -- unnecessary because not in [0,2pi] - get from y/x instead
# phi = np.pi + np.arctan2(y,x) # this phi will be in [0,2pi]
# this has been created in the file already, just read

# Construct polar coordinates with origin at center of storage region
x = np.array(r-711.2) # in cm
th_polar = np.array(np.arctan2(z,x))
r_polar = np.array(x/np.cos(th_polar))

for i in range(10):
    print('r=%f\tz=%f\nr_p=%f\tth_p=%f\n' %(r[i],z[i],r_polar[i],th_polar[i]))

coordinates = np.vstack((x,z)).T
print '\n\nCoordinates:\n', coordinates

err_V=1 ### To do: realistic error estimate for each V value ###

if len(sys.argv)>2:# 2nd arg given is for max n_order
    n_order=int(sys.argv[2])
else:
    n_order=5

# I believe pars is required to be 1-D array of parameters to be optimized.
# pars[2*n] are the A_n, pars[2*n+1] are the B_n.
pars = np.zeros(2*(n_order+1))
low_limit = np.zeros(2*(n_order+1))
hi_limit = np.zeros(2*(n_order+1))
# Get initial parameter values from file - if it exists
filename='pars_'+fstem+('_Bmad_%d.dat' % n_order)
exists = os.path.isfile(filename)
if exists:
    pars = np.loadtxt(filename, comments='!')  
    print ('Reading starting parameters from file:\n', filename)  
    if debug:
        print pars
filename = 'limits_multipoles_Bmad.dat'
exists = os.path.isfile(filename)
if exists:
    limits_data = np.loadtxt(filename, comments='!')    
    print ('Reading parameter limits from file:\n', filename) 
    if debug:
        print limits_data
    for i in range(2*n_order+2):
        low_limit[i]=limits_data[i,0]
        hi_limit[i]=limits_data[i,1]
else:
    for i in range(2*n_order+2):
        low_limit[i]=-np.inf
        hi_limit[i]=np.inf
print '\n\nInitial parameters ::\n', pars

par_bounds = (low_limit, hi_limit)
if debug:
    print ('Parameter bounds:\n', par_bounds)
print ('Parameter bounds:\n', par_bounds)

for i in range(2*n_order+2):
    if pars[i]<low_limit[i] or pars[i]>hi_limit[i]:
        print(' Uh oh!! Starting par value out of bounds...\n Par #%d=%.5e\tbounds [%.5e,%.5e]' % (i,pars[i],low_limit[i],hi_limit[i]))

# Calculate and minimize the residuals
residual(pars,coordinates,V,err_V,n_order)
fit_result = scipy.optimize.least_squares(residual, pars, method='trf', bounds=par_bounds,max_nfev=1e6,args=(coordinates, V, err_V, n_order),ftol=1e-15,xtol=1e-15)
#result_pars, flag = scipy.optimize.least_squares(residual, pars, args=(coordinates, V, err_V, n_order),ftol=1e-15,xtol=1e-12)
print ('\n\n Number of iterations:: %d / %d\n' % (iter_count, fit_result.nfev))
print ('Result status = ', fit_result.status)
print ('Cost = ', fit_result.cost)
print ('Optimality = ', fit_result.optimality)
print ('Fit parameters ::\n', fit_result.x)

# Output parameters to file
pars_file = open('pars_'+fstem+('_Bmad_%d.dat' % n_order), 'w+')
pars_file.write('!Fitted parameters up to n=%d\n' % n_order)
for i in range(len(pars)):
    pars_file.write('%.9e\n' % fit_result.x[i])
pars_file.close()

if not plot_results:
    exit()


# Plot results
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
