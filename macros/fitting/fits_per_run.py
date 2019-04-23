#!/usr/bin/env python

import os,sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
from pylab import frange

import fitresult

import readline,rlcompleter
readline.parse_and_bind('tab:complete')

cbo = False

if len(sys.argv)<2:
  print 'usage: fits_per_run.py DIRS save_fits_runX [save_fits_runY,...]'
  sys.exit(1)
dirnames = sys.argv[1:]


for dirname in dirnames:
  if not os.path.isdir(dirname):
    raise RuntimeError('directory "%s" not found!'%dirname)

fitresults = []
for dirname in dirnames:
  filepath = os.path.join(dirname,'savedata.pkl')
  try:
    infile = open(filepath,'rb')
    filedata = pickle.load(infile)
    infile.close()
    fitresults += [ [
      int(filepath.split('/')[0].split('_')[-1]), 
      filedata['fivepar_fit_result'], 
      filedata['exp_fit_result'], 
      filedata['ratio_fit_result']
    ] ]
  except:
    print RuntimeWarning('failed loading fit data from %s!'%filepath)
firstrun = fitresults[0][0]
for i in range(len(fitresults)):
  fitresults[i][0] -= firstrun

#param_index = fitresults[0][3].param_names.index('r')
#param_index = fitresults[0][3].param_names.index('omega_cbo')
def grab_fitparams(param_name,fittype):
  fittype = {
    'fivepar':1, 'fiveparam':1, '5par':1, '5param':1, 
    'exp':2, 
    'ratio':3
  }[fittype]
  param_index = fitresults[0][fittype].param_names.index(param_name)
  xlist,ylist,yerrlist = [],[],[]
  for fr in fitresults:
    xvalue,yvalue,yerrvalue = None,None,None
    if type(fr[fittype].params) != type(None):
      xvalue = fr[0] + (0.05,0.2,0.35)[fittype-1]
      yvalue = fr[fittype].params[param_index]
      if type(fr[fittype].param_errors) != type(None):
        yerrvalue = fr[fittype].param_errors[param_index]
        xlist += [ xvalue ]
        ylist += [ yvalue ]
        yerrlist += [ yerrvalue ]
  return xlist,ylist,yerrlist


def plot_param_results(paramname,omit_fittype=tuple()):
  '''
  
  omit_fittype can be one or more of 'ratio', 'fivepar', or 'exp'
  '''
  if not hasattr(omit_fittype,'__iter__'):
    omit_fittype = tuple((omit_fittype,))
  plot = plt.figure(figsize=(12,9)).add_subplot(1,1,1,
    title='Fit Results: %s'%paramname, 
    xlabel='Run-%d'%firstrun, ylabel=paramname
  )
  plot.figure.subplots_adjust(left=0.11,bottom=0.09,right=0.93,top=0.92)
  if 'ratio' not in omit_fittype:
    xlist,ylist,yerrlist = grab_fitparams(paramname,fittype='ratio')
    plot.errorbar(xlist, ylist, yerr=yerrlist, 
      linestyle='', marker='.', capsize=0, 
      label='Ratio Fit'
    )
  if 'fivepar' not in omit_fittype:
    plot.errorbar(
      *grab_fitparams(paramname,fittype='fivepar'), 
      linestyle='', marker='.', capsize=0, 
      label='Fiveparam Fit'
    )
    plotlimits = plot.axis()
  if 'exp' not in omit_fittype:
    plot.errorbar(
      *grab_fitparams(paramname,fittype='exp'), 
      linestyle='', marker='.', capsize=0, 
      label='Exp Fit'
    )
    plot.axis(plotlimits)
  plot.axis(xmax=fitresults[-1][0]+1.)
  plot.legend(loc=0)
  return plot

def save_plot(axis,fit_param_name,filetype='png'):
  outfile_basename = 'run_by_run-%s'%fit_param_name
  print 'writing out %s.%s...'%(outfile_basename,filetype)
  axis.figure.savefig(outfile_basename+'.'+filetype)
  print 'writing out %s.pkl...'%outfile_basename
  pickle.dump(axis.figure,open(outfile_basename+'.pkl','wb'))

A_plot = plot_param_results('A',omit_fittype='exp')
A_plot.figure.show()
save_plot(A_plot,'A')

tau_a_plot = plot_param_results('tau_a',omit_fittype='ratio')
#tau_a_plot.axis(ymin=0,ymax=100000)
tau_a_plot.axis(ymin=64100,ymax=65000)
tau_a_plot.figure.show()
save_plot(tau_a_plot,'tau_a')

omega_a_plot = plt.figure(figsize=(12,9)).add_subplot(1,1,1,
  title='Fit Results: $\omega_a$', 
  xlabel='Run-%d'%firstrun, ylabel='$\omega_a$ [rad*GHz]'
)
omega_a_plot.figure.subplots_adjust(left=0.11,bottom=0.09,right=0.93,top=0.92)
xlist,ylist,yerrlist = grab_fitparams('r',fittype='ratio')
omega_a_plot.errorbar(xlist, ylist, yerr=yerrlist, 
  linestyle='', marker='.', capsize=0, 
  label='Ratio Fit'
)
omega_a_plot.errorbar(
  *grab_fitparams('r',fittype='fivepar'), 
  linestyle='', marker='.', capsize=0, 
  label='Fiveparam Fit'
)
#plotlimits = omega_a_plot.axis()
#omega_a_plot.errorbar(
#  *grab_fitparams('r',fittype='exp'), 
#  linestyle='', marker='.', capsize=0, 
#  label='Exp Fit'
#)
#omega_a_plot.axis(plotlimits)
omega_a_plot.axis(xmax=fitresults[-1][0]+1.)
omega_a_plot.legend(loc=0)
omega_a_plot.figure.show()
save_plot(omega_a_plot,'omega_a')

if cbo:
  t_cbo_plot = plot_param_results('t_cbo')
  t_cbo_plot.axis(ymin=-400000,ymax=600000)
  t_cbo_plot.figure.show()
  save_plot(t_cbo_plot,'t_cbo')

N0_plot = plot_param_results('N0',omit_fittype='ratio')
#N0_plot.axis(ymin=0,ymax=1e+06)
N0_plot.figure.show()
save_plot(N0_plot,'N0')

if cbo:
  omega_cbo_plot = plt.figure(figsize=(12,9)).add_subplot(1,1,1,
    title='Fit Results: $\omega_\mathrm{CBO}$', 
    xlabel='Run-%d'%firstrun, ylabel='$\omega_\mathrm{CBO}$ [rad*GHz]'
  )
  omega_cbo_plot.figure.subplots_adjust(left=0.11,bottom=0.09,right=0.93,top=0.92)
  omega_cbo_plot.errorbar(
    *grab_fitparams('omega_cbo',fittype='ratio'), 
    linestyle='', marker='.', capsize=0, 
    label='Ratio Fit'
  )
  omega_cbo_plot.errorbar(
    *grab_fitparams('omega_cbo',fittype='fiveparam'), 
    linestyle='', marker='.', capsize=0, 
    label='Fiveparam Fit'
  )
  plotlimits = omega_cbo_plot.axis()
  omega_cbo_plot.errorbar(
    *grab_fitparams('omega_cbo',fittype='exp'), 
    linestyle='', marker='.', capsize=0, 
    label='Exp Fit'
  )
  omega_cbo_plot.axis(
    ymin=0.00224,ymax=0.00238,
    xmax=fitresults[-1][0]+1.
  )
  omega_cbo_plot.legend(loc=0)
  omega_cbo_plot.figure.show()
  save_plot(omega_cbo_plot,'omega_cbo')





