import numpy as np
import math

#
# Open OPERA map file and read in data
#
data = np.loadtxt('/Users/ekargian/Documents/g-2/Magnets/ElectrostaticQuads/E-field/OPERA-simulations/short_quad_displacements_demo_added_BC_inner_plate_Q1S.table', comments='!')
r=data[:,0]
th=data[:,1]
z=data[:,2]
V=data[:,3]
    
# Open output file
f= open("short_quad_displacements_demo_added_BC_inner_plate_Q1S.dat","w+")

accepted_entries=0
for i in range(len(r)):
    if math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))<5.1 and abs(V[i])>0.0001:
        for item in data[i]:
            f.write('%.11E  '%item, )
        f.write('\n')
        accepted_entries += 1
    
f.close()
print 'Total accepted entries :: ', accepted_entries