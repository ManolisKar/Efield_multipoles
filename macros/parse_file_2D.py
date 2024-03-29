import numpy as np
import math
import sys

#
# Open OPERA map file and read in data
#

if len(sys.argv)>1:# 1st arg given is for file stem
    filestem=str(sys.argv[1])
else:
    filestem='short_quad_displacements_demo_added_BC_inner_plate_Q1S'

data = np.loadtxt('/Users/ekargian/Documents/g-2/Magnets/ElectrostaticQuads/E-field/OPERA-simulations/map_'+filestem+'.dat', comments='!')
r=data[:,3]
th=data[:,4]
z=data[:,2]
V=data[:,5]
    
# Open output file
f= open('map_'+filestem+'.dat',"w+")
#f= open('test.dat',"w+")

accepted_entries=0
for i in range(len(r)):
    #if (abs(math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))-4.50)<0.1 or abs(math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))-4.00)<0.1 or abs(math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))-3.50)<0.1 or abs(math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))-3.00)<0.1 or abs(math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))-2.50)<0.1 or abs(math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))-2.00)<0.1 or abs(math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))-1.50)<0.1 or abs(math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))-1.00)<0.1 or abs(math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))-0.50)<0.1) and abs(th[i]-6.45)<0.02 and abs(V[i])>0.0001: 
    #if math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))<4.51 and math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))>4.49: 
    if math.sqrt(math.pow(r[i]-711.2,2)+math.pow(z[i],2))<4.8:         
    #if i%20==0:
        for item in data[i]:
            f.write('%.11E  '%item, )
        f.write('\n')
        accepted_entries += 1
    
f.close()
print 'Total accepted entries :: ', accepted_entries
