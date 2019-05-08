#!/usr/bin/env python
import numpy as np
import scipy
from math_functions import *
import matplotlib.pyplot as plt
import sys

iter_count=0

__doc__ = '''Fit data of potential from an OPERA map file
to a function of multipoles in toroidal coordinates.
'''
# Should input as argument the number of order of multipoles to go up to.

#def olvers(a,b,c,x):
#    # Olver's hypergeometric function 
#    # See https://dlmf.nist.gov/15.1
#    return hyp2f1(a,b,c,x)/gamma(c)

#def legendre_P(m,n,x):
#    # Legendre/Ferrers function of first kind
#    # m,n can be real (non-integers)
#    # See https://dlmf.nist.gov/14.3
#    return ((1+x)/(1-x))**(m/2) * olvers(n+1,-n,1-m,0.5-0.5*x)

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
                nm_id = 4*( sum(range(n+1)) + m ) #n-m position in the parameter array
                A_nm=pars[0+nm_id]
                B_nm=pars[1+nm_id]
                eta_nm=pars[2+nm_id]
                phi_nm=pars[3+nm_id]
                model_V[i] += (A_nm*legendre_Q_toroidal(m,n,mu_i) + B_nm*legendre_P_toroidal(m,n,mu_i)) * np.cos(n*(eta_i-eta_nm)) * np.cos(m*(phi_i-phi_nm))
        model_V[i] *= f_mu_eta[i]   

    return model_V 

def residual(pars, coordinates, V, err_V, n_order, m_order):
    global iter_count

    model_V = model_potential(pars, coordinates, n_order, m_order)
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
# Open OPERA map file and read in data
# (this file has been parsed to only include measurements within a 9.5 cm diameter circle)
#data = np.loadtxt('short_quad_ideal.dat', comments='!')
#data = np.loadtxt('short_quad_ideal_Vgt10.dat', comments='!')
if len(sys.argv)>1:# 1st arg given is for file stem
    fstem=str(sys.argv[1])
else:
    fstem='test3000'
print 'Processing file '+fstem+'.dat'
data = np.loadtxt(fstem+'.dat', comments='!')
x=data[:,0]
y=data[:,1]
z=data[:,2]
r=data[:,3]
# phi=data[:,4] -- unnecessary because not in [0,2pi] - get from y/x instead
# phi = np.pi + np.arctan2(y,x) # this phi will be in [0,2pi]
# this has been created in the file already, just read

V=data[:,8]

#mu=data[:9]
#eta=data[:10]
#phi=data[:11]
## Read the toroidal coordinates in multi-D array
coordinates = data[:,9:12]

# distance of each "pole" in the bipolar coordinate system from origin.
# defined in cm, same as coordinates in file.
a=5.756 

err_V=1 ### To do: realistic error estimate for each V value ###

if len(sys.argv)>2:# 2nd arg given is for max n_order
    n_order=int(sys.argv[2])
else:
    n_order=5
if len(sys.argv)>3:# 3rd arg given is for max m_order
    m_order=int(sys.argv[3])
else:
    m_order=n_order

pars = np.zeros(4*( sum(range(n_order+1)) + m_order + 1 ))
# I believe pars is required to be 1-D array of parameters to be optimized.
# pars[0+4*n*m] are the A_nm, pars[1+4*n*m] are the B_nm, 
# pars[2+4*n*m] are the eta_nm, pars[3+4*n*m] are the phi_nm.
# Get some initial values from file, up to appropriate order:
starting_pars = np.loadtxt('pars_'+fstem+'.dat', comments='!')
for n in range(n_order+1):
    for m in range(n+1): # for now m goes up to n, or up to m_order
        if m>m_order: continue
        nm_id = 4*( sum(range(n+1)) + m ) #n-m position in the parameter array
        for ipar in range(4):
            if nm_id+ipar>=len(starting_pars): break
            pars[nm_id+ipar] = starting_pars[nm_id+ipar]

'''
n=0
m=0
nm_id = 4*( sum(range(n+1)) + m )
pars[0+nm_id] = 5
pars[1+nm_id] = 10
pars[2+nm_id] = 0.1
pars[3+nm_id] = 0.1
n=1
nm_id = 4*( sum(range(n+1)) + m )
pars[0+nm_id] = 1
pars[1+nm_id] = 1
pars[2+nm_id] = 0.2
pars[3+nm_id] = 2
m=1
nm_id = 4*( sum(range(n+1)) + m )
pars[0+nm_id] = 2
pars[1+nm_id] = 2
pars[2+nm_id] = 0.2
pars[3+nm_id] = 2
'''
print '\n\nInitial parameters ::\n', pars

# Calculate and minimize the residuals
residual(pars,coordinates,V,err_V,n_order,m_order)
result_pars, flag = scipy.optimize.leastsq(residual, pars, args=(coordinates, V, err_V, n_order, m_order),ftol=1e-15,xtol=1e-12)
print ('Result flag = %d\n Fit parameters ::\n'%flag), result_pars

# Output parameters to file
pars_file = open('pars_'+fstem+('.dat' % n_order), 'w+')
for i in range(len(pars)):
    pars_file.write(pars[i])


# Plot results
mu=np.array(data[:,9])
eta=np.array(data[:,10])
phi=np.array(data[:,11])
V=np.array(V)
final_residuals = np.array(residual(result_pars,coordinates,V,err_V,n_order,m_order))

# Plot V vs coordinates
V_plot = plt.figure(figsize=(18,7)).add_subplot(1,1,1,
    title='V vs mu',xlabel='mu'
)
V_plot.plot(
    mu, V, 
    marker='', color='black', 
    linestyle='-', linewidth=0.75, 
    label='V'
)
V_plot.figure.savefig('V_vs_mu.png')
V_plot.figure.show()

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

# Plot residuals vs all coordinates
resid_mu_plot = plt.figure(figsize=(18,7)).add_subplot(1,1,1,
    title='Residuals vs mu',xlabel='mu'
)
resid_mu_plot.plot(
    mu, final_residuals, 
    marker='', color='black', 
    linestyle='-', linewidth=0.75, 
    label='V Residuals'
)
resid_mu_plot.figure.savefig('res_vs_mu.png')
resid_mu_plot.figure.show()
