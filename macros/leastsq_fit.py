import numpy as np

def residual(pars, coordinates, V, err_V):
    pol0 = pars[0]
    pol1 = pars[1]
    pol2 = pars[2]
    pol3 = pars[3]

    x = coordinates[:,0]
    y = coordinates[:,1]
    z = coordinates[:,2]

    model = pol0 + pol1*x + pol2*y + pol3*z
    print 'model ::\n', model 

    norm_residual = (V - model)/err_V
    print 'normalized residual ::\n', norm_residual

    return norm_residual

 #   return (data-model) / eps_data


from scipy.optimize import leastsq


data = np.loadtxt('test.dat',comments='!')
print 'data:\n', data
#x=data[:,0]
#y=data[:,1]
#z=data[:,2]
coordinates = data[:,0:3]
print 'coordinates:\n', coordinates
V=data[:,3]

err_V=0.5

pars = [0.1, 1.2, 1.3, 0.9]
print 'pars = ', pars

residual(pars,coordinates,V,err_V)
result = leastsq(residual, pars, args=(coordinates, V, err_V))
print 'Result :: ', result