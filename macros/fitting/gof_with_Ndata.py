#!/usr/bin/env python

__doc__ = '''Naive fits degrade with higher counts.

I fit various combinations of runs between 15923 and 15928 to
see how the three different fits degrade with higher statistical
precision.  The fits take no systematic problems into account.
'''

import math
import numpy as np
import matplotlib.pyplot as plt

import readline,rlcompleter
readline.parse_and_bind('tab:complete')


gof_by_size = np.array([
  [04782813.0, 0.334724,  0.956242, 0.400636 ],
  [05126862.0, 0.344172,  0.629838, 0.186302 ],
  [05181144.0, 0.026289,  0.017196, 0.317741 ],
  [09719796.0, 0.001187,  0.355641, 0.134802 ],
  [09865274.0, 0.036871,  0.486452, 0.148467 ],
  [10094897.0, 0.023669,  0.545283, 0.917546 ],
  [10123994.0, 3.15e-06,  0.374818, 0.240277 ],
  [14846658.0, 3.937e-06, 0.192549, 0.012153 ],
  [15090819.0, 0.005812,  0.474795, 0.646355 ],
  [19989268.0, 0.000171,  0.460516, 0.236411 ],
  [24810615.0, 3.862e-07, 0.260177, 0.535567 ],
  [30084165.0, 5.359e-06, 0.058452, 0.767695 ],
  [34835926.0, 2.05e-10,  0.039935, 0.034444 ], 
  [35050990.0, 2.061e-07, 0.116131, 0.782028 ],
  [40392171.0, 4.527e-11, 0.008007, 0.572863 ],
  [54894780.0, 1.11e-16,  0.000131, 0.639172 ] 
])
gof_by_size[:,0] /= 1e+06 # convert to millions of positrons


plot = plt.figure(figsize=(8,8)).add_subplot(1,1,1,
  title='Goodness-of-Fit', 
  xlabel='Number of Positrons [Millions]', 
  ylabel='P-value ($\chi^2$)'
)
plot.set_yscale('log')
plot.plot(
  gof_by_size[:,0], gof_by_size[:,1], 
  label='Five-param', 
  marker='o', linestyle=''
)
plot.plot(
  gof_by_size[:,0], gof_by_size[:,2], 
  label='Exponential', 
  marker='o', linestyle=''
)
plot.plot(
  gof_by_size[:,0], gof_by_size[:,3], 
  label='Ratio', 
  marker='o', linestyle=''
)

plot.legend(loc=0)
plot.figure.show()
plot.figure.savefig('gof_with_Ndata.png')

