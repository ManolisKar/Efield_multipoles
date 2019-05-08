import numpy as np
import math

#
# Open OPERA map file and read in data
#
data = np.loadtxt('/Users/ekargian/Documents/g-2/Magnets/ElectrostaticQuads/E-field/OPERA-simulations/short_quad_ideal.table', comments='!')
x=data[:,0]
y=data[:,1]
z=data[:,2]
r=data[:,3]
# phi=data[:,4] -- unnecessary because not in [0,2pi] - get from y/x instead
phi = np.arctan2(y,x) # this phi will be within [-pi/2,pi/2] -- will be fixed below

V=data[:,8]

# distance of each "pole" in the bipolar coordinate system from origin.
# defined in cm, same as coordinates in file.
a=5.756 

# toroidal coordinates
mu = np.arctanh(2*a*r/(r*r+z*z+a*a))
eta = np.arctan(2*a*z/(r*r+z*z-a*a))
# the above returns eta in [-pi/2,pi/2], fix below

for i in range(len(eta)):
    sineta = np.sin(eta[i])
    taneta = np.tan(eta[i])
    if sineta>=0:
        if taneta>=0: #1st quadrant
            continue # eta value already should be correct
        else: 
            #eta[i]<0, 2nd quadrant
            eta[i] = np.pi + eta[i]
    else:
        if taneta>=0: #eta[i]>0, 3rd quadrant
            eta[i] = np.pi + eta[i]
        else: 
            #eta[i]<0, 4th quadrant
            eta[i] = 2*np.pi + eta[i]
    # Lastly, shift phi values to the [0,2pi] range
    if phi[i]<0:
        phi[i] = 2*np.pi + phi[i]
    
# Open output file
f= open("short_quad_ideal_outererR.dat","w+")

accepted_entries=0
for i in range(len(r)):
    if math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))<5.05 and math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))>4.995 and abs(V[i])>5000:
        for item in data[i]:
            f.write('%.11E  '%item, )
        f.write('%.11E  %.11E  %.11E\n'%(mu[i],eta[i],phi[i]))
        accepted_entries += 1
    
f.close()
print 'Total accepted entries :: ', accepted_entries