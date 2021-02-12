#!/usr/bin/env python
import numpy as np
import scipy
from math_functions import *
import matplotlib.pyplot as plt
import sys
import os



__doc__ = '''
    Calculate appropaiet limits of variation for multipole parameters.
    The limits show that each term (other than the normal quadrupole) 
    must contribute less than contribution_limit (~2000V) at r=4.5cm.

    1st arg: switch for 2D/Bmad/toroidal
'''

if len(sys.argv)>1:# 1st arg given is for switch
    switch=str(sys.argv[1])
    if (switch is not '2D') and (switch is not 'Bmad') and (switch is not 'toroidal'):
        print 'Please pass valid switch value\n'
        sys.exit()
else:
    switch='2D'



debug = 1

# Limit to contribution from each term (other than quadrupole) at r=4.5cm, in V:
contribution_limit = 10000 

## Calculate limits up to order 50. 
## The fitting script can read only up to the desired order
nmax=50

## Distance where potential contribution is calculated:
r0=4.5 

# Output file:

outfile = open('limits_multipoles.dat', 'w+')
outfile.write('!! Limits to each multipole term\n')


for n_order in range(nmax):

    if switch is '2D':
        ## Each an and bn multipole must obey abs([an or bn])*r**n < contribution_limit
        ## This formula applied for the eqn used in the NIM paper:
        limit = contribution_limit / np.power(r0,n_order)
        
        ## Output limits for a_n
        if n_order==0:        
            ## For n=0 this contribution is just 0, reagrdless of a_n. Just set limits to zero.
            outfile.write('-1e-15\t1e-15\n')
        else:
            outfile.write('-%.9e\t%.9e\n' % (limit, limit))

        ## Output limits for b_n
        if n_order==2:
            ## If n==2, set the b2 (normal quadrupole) limits to infinity
            outfile.write('-1e20\t1e20\n')
        else:
            ## Otherwise b_n should have same limits as a_n
            outfile.write('-%.9e\t%.9e\n' % (limit, limit))

    elif switch is 'Bmad':
        ## Not sure how to apply limits from the Bmad formula with complex numbers
        continue


outfile.close()


