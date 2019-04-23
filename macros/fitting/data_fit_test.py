#!/usr/bin/env python

__doc__ = '''Fit data (art files, or dst skim or histogram pickle files).

Accepts multiple files, but only of a single type per invocation.

Creates a directory ./save_fits/ and writes out some results.  This script
REFUSES to overwrite the save_fits directory, so either move it or delete
it with rm -rf.
'''

### some settings ###

verbose = True

# write out plots of all the fits
save_fits = True      # fit info (and plots if enabled) to the 'save_fits_dir'
save_fits_dir = 'save_fits'
make_plots = True     # show the histograms & fits to them
popup_plots = True    # does nothing if make_plots==False
log_console = True    # copy console output to data_fit_test.log

max_art_events = None
max_skimmed_clusters = None

fit_tmin,fit_tmax = 30000,700000
clip_ratio_fit_at_1st_empty_bin = True # False = NotImplementedError
#fit_tmin_override,fit_tmax_override = None,None
fit_tmin_override,fit_tmax_override = 33000,660000

#ft_tmin,ft_tmax = None,None
ft_tmin,ft_tmax = 33000,150000

binwidth = 149.2 # ns


#### load some things ####

# frange is better than arange and linspace, but is tied to matplotlib, which
# is inconvenient if we need to avoid X11 forwarding latency
import matplotlib
if not popup_plots or not make_plots:
  matplotlib.use('svg') # prevent matplotlib from trying to make windows
elif make_plots:
  import matplotlib.pyplot as plt
from pylab import frange

# system imports
import os,sys
import pickle
import math
import random
import numpy as np
import scipy
import scipy.stats as sstats

# our custom stuff
import Blinders
main_verbose = verbose
from util import *     # this replaces with util.verbose
verbose = main_verbose
import fitresult

# tab completion is awesome
import readline,rlcompleter
readline.parse_and_bind('tab:complete')

# copy all output to a logfile
if log_console: tee = Tee()


#### read input file arguments ####

infiles = sys.argv[1:]
if len(infiles)<1:
  raise RuntimeError('Give me a list of input files! (*.root or *.pkl or *.hist)')
input_filetype = infiles[0].split('.')[-1]
if input_filetype=='root' or input_filetype=='art':
  print 'using art files'
  for thing in infiles[1:]:
    assert thing.split('.')[-1]==input_filetype # don't mix file extensions!
elif input_filetype=='pkl': print 'using pickle file(s)'
elif input_filetype=='hist': print 'using pre-histogrammed data in pickle files'
else: raise RuntimeError(
  'input_filetype "%s" not understood! use one of art/root, pkl, or hist...'%input_filetype
)
if input_filetype=='root': input_filetype = 'art' # equivalent
for infile in infiles:
  if not os.path.isfile(infile):
    raise RuntimeError('Input file %s does not exit!'%infile)



#### do some preprocessing ####

# figure out where to save diagnostic plots
if save_fits:
  if os.path.isdir(save_fits_dir):
    raise RuntimeError('ERROR: Fits output directory "%s" already exists!'%save_fits_dir)
  os.mkdir(save_fits_dir)
  if not os.path.isdir(save_fits_dir):
    raise RuntimeError('ERROR: Failed to create output directory "%s"!'%save_fits_dir)
def save_figure(fig,basename,extension='png',pickle_plots=True,):
  '''Save a plot/figure to the specified basefile name.
  
  Intended to save a matplotlib Figure or Axes object.
  
  If pickle=True (default) this also saves a pickled version of the plot.
  
  This function does NOTHING if either save_fits or make_plots is False.
  '''
  if not type(fig)==matplotlib.figure.Figure:
    try:
      fig = fig.figure
      assert type(fig)==matplotlib.figure.Figure
    except: RuntimeError('save_figure() requires a matplotlib.figure.Figure (or at least something that has a data member ".figure" which references a matplotlib Figure)!')
  if save_fits and make_plots:
    filename = os.path.join(save_fits_dir,basename+'.'+extension)
    print 'writing out %s...'%filename
    fig.savefig(filename)
    if pickle_plots:
      filename = os.path.join(save_fits_dir,basename+'.pkl')
      print 'writing out %s...'%filename
      outfile = open(filename,'wb')
      pickle.dump(fig,outfile)
      outfile.close()


#### do stuff! ####


# this has keys for each of the (4?) histograms/binnings that need
# to be built/managed:
# alldata: everything from zero to the end of the fill
# fitdata: all data, but trimmed to fitting range (for wiggle fit)
#    med1: binning for N(t) histogram (fitting domain?)
#    med2: binning for second N(t) histogram
#     low: binning for N(t-T_a0/2) histogram
#    high: binning for N(t+T_a0/2) histogram
#   trial: ALL bin edges (sorted) for the `build_sample_hist` 
#         function (NOT evenly-spaced)
bin_edges,bin_centers,bin_counts,N_bins = {},{},{},{}

# bins for ALL data (zero to end of fill)
bin_edges['alldata'] = frange(0,800000,binwidth)
bin_centers['alldata'] = np.array([0.5*(bin_edges['alldata'][i]+bin_edges['alldata'][i+1]) for i in range(len(bin_edges['alldata'])-1)])

# bins for N(t), N(t-Ta/2), and N(t+Ta/2)
bin_edges['fitdata'] = frange(fit_tmin,fit_tmax,binwidth)
bin_centers['fitdata'] = np.array([0.5*(bin_edges['fitdata'][i]+bin_edges['fitdata'][i+1]) for i in range(len(bin_edges['fitdata'])-1)])
bin_edges['low'] = bin_edges['fitdata'] - half_T_a0
bin_centers['low'] = bin_centers['fitdata'] - half_T_a0
bin_edges['high'] = bin_edges['fitdata'] + half_T_a0
bin_centers['high'] = bin_centers['fitdata'] + half_T_a0
bin_edges['med1'] = bin_edges['fitdata']
bin_centers['med1'] = bin_centers['fitdata']
bin_edges['med2'] = bin_edges['fitdata']
bin_centers['med2'] = bin_centers['fitdata']

## bins for sample histogram builder (and later recombination)
#bin_edges['trial'] = np.append(bin_edges['alldata'],bin_edges['fitdata'])
#bin_edges['trial'] = np.append(bin_edges['trial'],bin_edges['low'])
#bin_edges['trial'] = np.append(bin_edges['trial'],bin_edges['high'])
#bin_edges['trial'].sort()
#bin_centers['trial'] = np.array([0.5*(bin_edges['trial'][i]+bin_edges['trial'][i+1]) for i in range(len(bin_edges['trial'])-1)])

for histtype in ('alldata','low','high'): #,'trial'):
  nbins = len(bin_centers[histtype])
  N_bins[histtype] = nbins
  bin_counts[histtype] = np.zeros(nbins)
N_bins['fitdata'] = N_bins['med1'] = N_bins['med2'] = len(bin_centers['fitdata'])
bin_counts['fitdata'] = np.zeros(N_bins['fitdata'])
bin_counts['med1'] = np.zeros(N_bins['med1'])
bin_counts['med2'] = np.zeros(N_bins['med2'])





fitrecords_dtype = [ 
  ('time','f8'), 
  ('energy','f8'), 
  ('calorimeterIndex','i1'), 
  ('valid','i4'), 
  ('flagged','i1'), 
  ('chiSquared','f8'), 
  ('fillIndex','i4'), 
  ('scale','f8'), 
  ('run','i4'), 
  ('subrun','i4'), 
  ('fill','i4'),
  ('abstime_s','i4'), 
  ('abstime_ns','i4')
]

just_skim = True
if input_filetype=='art':
  import heist
  gps_tag = heist.InputTag(quicktag='gm2common::GPSArtRecord_GPSUnpacker_GPSUnpacker0_offline')
  cluster_tag = heist.InputTag(
    # it's one of these
    #quicktag='gm2reconeast::GlobalFitArtRecords_globalFit_fitter_offline'
    #quicktag='gm2reconeast::GlobalFitArtRecords_energyPartition_partition_offline'
    quicktag='gm2reconeast::GlobalFitArtRecords_energyCalibrator_calibrator_offline'
  )
  ## old test (single art file):
  #infile = '/media/sf_hosthome/data/v9_09_scratch/'
  ##infile += 'run1_5032B/runs_16000/16461/'  # <== comment out to use local file
  #infile += 'gm2offline_full_12521572_16461.00409.root'
  artreader = heist.ArtFileReader(infiles)
  fitrecords = []
  for evt in artreader.event_loop(nmax=max_art_events):
    clusters = evt.get_record(cluster_tag)
    if not clusters: continue
    gps_record = evt.get_record(gps_tag)
    
    if verbose: print 'Found %d clusters in %s'%(len(clusters),evt.get_label())
    run,subrun,fill = evt.get_ID()
    abstime_s,abstime_ns = gps_record.unixTimeFE,gps_record.unixTimeFEFraction
    #print 'fillIndex:',clusters[0].fillIndex
    for c in clusters:
      row = ( 
        1.25*c.time, c.energy, c.calorimeterIndex, 
        c.valid, c.flagged, c.chiSquared, c.fillIndex, c.scale, 
        run, subrun, fill, 
        abstime_s, abstime_ns
      )
      fitrecords += [ row ]
  fitrecords = np.array(fitrecords,dtype=fitrecords_dtype)

  if verbose: print 'sorting records by time...'
  fitrecords.sort(order='time')

  outfile = 'clusters.pkl'
  print 'writing records to %s...'%outfile
  pickle.dump(fitrecords,open(outfile,'wb'))
  print '...done.'
  if just_skim:
    print 'stopping after writing fitrecords!'
    sys.exit(0)
elif input_filetype=='pkl':
  # read first file to start, then just verify that all later files
  # contain numpy arrays with the same field labels (and merge them)
  first_filename = infiles[0]
  outstr = '%s: '%first_filename
  try:
    fitrecords = pickle.load(open(first_filename,'rb'))
    outstr += '%d fit records'%len(fitrecords)
  except:
    outstr += 'FAILED!'
  print outstr
  for filename in infiles[1:]:
    #fitrecords = pickle.load(open(filename,'rb'))
    outstr = '%s: '%filename
    try:
      pickle_input_file = open(filename,'rb')
      filedata = pickle.load(pickle_input_file)
      pickle_input_file.close()
      fitrecords = np.append(fitrecords,filedata) # TODO: is np.concatenate((,,),out=out) faster?
      outstr += '%d fit records'%len(filedata)
    except:
      outstr += 'FAILED!'
    if max_skimmed_clusters and len(fitrecords)>max_skimmed_clusters:
      print 'Reached >=%d clusters!'%max_skimmed_clusters
      fitrecords = fitrecords[:max_skimmed_clusters]
      break
    print outstr
elif input_filetype=='hist':
  # read first file to start, then just verify that all later files
  # contain dictionaries with the same keys (and merge them)
  outstr = '%s: '%infiles[0]
  try:
    first_file = open(infiles[0],'rb')
    cluster_hist_data = pickle.load(first_file)
    assert type(cluster_hist_data)==dict
    assert 'bin_centers' in cluster_hist_data.keys()
    first_file.close()
    outstr += 'histograms with counts\n  low: %d\n  med1: %d\n  med2: %d\n  high: %d\n'%(
      np.sum(cluster_hist_data['bin_counts']['low']), 
      np.sum(cluster_hist_data['bin_counts']['med1']), 
      np.sum(cluster_hist_data['bin_counts']['med2']), 
      np.sum(cluster_hist_data['bin_counts']['high'])
    )
  except: outstr += 'FAILED!'
  print outstr
  for filename in infiles[1:]:
    outstr = '%s: '%filename
    try:
      pickle_input_file = open(filename,'rb')
      filedata = pickle.load(pickle_input_file)
      pickle_input_file.close()
      # do some checks?
      # NOTE: checking for identical N_bins when comparing histograms, but 
      # currently NOT checking bin edges and bin centers are the same!
      paranoid_checks = True
      if paranoid_checks==True:
        keys1 = cluster_hist_data.keys()
        keys1.sort()
        keys2 = filedata.keys()
        keys2.sort()
        # check for identical set of keys
        if keys1!=keys2: raise RuntimeError(
          'keys do not match! ( %s vs %s )'%(str(keys1),str(keys2))
        )
        for key in ('binwidth','half_T_a0','omega_a0','fit_tmin','fit_tmax'):
          if cluster_hist_data[key]!=filedata[key]:
            raise RuntimeError(
              'basic configuration mismatch!\n  %s=%6.4g vs %s=%6.4h'%(
                key, cluster_hist_data[key], key, filedata[key]
              )
            )
        for key in cluster_hist_data['N_bins'].keys(): # check N_bins
          n1 = cluster_hist_data['N_bins'][key]
          n2 = filedata['N_bins'][key]
          if n1!=n2:
            raise RuntimeError('unequal number of bins for %s! (%d vs %d)'%(key,n1,n2))
        # check for identical lower & upper binning range
        for key in cluster_hist_data['bin_edges'].keys():
          low1 = cluster_hist_data['bin_edges'][key][0]
          low2 = filedata['bin_edges'][key][0]
          assert low1==low2  # lower bin edges must match!
          high1 = filedata['bin_edges'][key][-1]
          high2 = cluster_hist_data['bin_edges'][key][-1]
          assert high1==high2  # upper bin edges must match!
      # report
      outstr += 'histograms with counts\n  low: %d\n  med1: %d\n  med2: %d\n  high: %d\n'%(
        np.sum(cluster_hist_data['bin_counts']['low']), 
        np.sum(cluster_hist_data['bin_counts']['med1']), 
        np.sum(cluster_hist_data['bin_counts']['med2']), 
        np.sum(cluster_hist_data['bin_counts']['high'])
      )
      # if acceptable, merge new data into existing structure
      for counts_key in cluster_hist_data['bin_counts'].keys():
        #cluster_hist_data['bin_counts'][counts_key] = np.append(
        #    cluster_hist_data['bin_counts'][counts_key], 
        #    filedata['bin_counts'][counts_key]
        #)
        cluster_hist_data['bin_counts'][counts_key] \
          += filedata['bin_counts'][counts_key]
    except:
      outstr += 'FAILED!'
    print outstr
else: raise RuntimeError('input_filetype "%s" not recognized!'%input_filetype)

#raise NotImplementedError('...here be drag0ns...')

if input_filetype in ('art','pkl'):
  if verbose: print 'performing 4-way binning...'
  #bin_counts['alldata'] = np.histogram(np.ravel(trial_samples),bins=bin_edges['alldata'])[0]
  #bin_counts['low'] = np.histogram(trial_samples[0],bins=bin_edges['low'])[0]
  #bin_counts['med1'] = np.histogram(trial_samples[1],bins=bin_edges['med'])[0]
  #bin_counts['med2'] = np.histogram(trial_samples[2],bins=bin_edges['med'])[0]
  #bin_counts['high'] = np.histogram(trial_samples[3],bins=bin_edges['high'])[0]
  ## bin_counts['med'] is a misnomer in this case...
  #bin_counts['med'] = bin_counts['low']+bin_counts['med1']+bin_counts['med2']+bin_counts['high']
  def do_4way_binning(tlist,bins1,bins2,bins3,bins4):
    '''Bin a list of points four ways (preserving orders of magnitude).
    
    This version is the simplest, dumbest, and fastest thing to write.
    
    (Leaves off up to 3 points at the end.)'''
    sorted_tlist = np.sort(tlist)
    t_list_1,t_list_2,t_list_3,t_list_4 = [],[],[],[]
    for i in range(len(tlist)/4):
      fouri = 4*i
      t_list_1 += [ tlist[fouri] ]
      t_list_2 += [ tlist[fouri+1] ]
      t_list_3 += [ tlist[fouri+2] ]
      t_list_4 += [ tlist[fouri+3] ]
    #for i in range(4*(len(tlist)/4),len(tlist)):
    #remainder = len(tlist)%4
    #if remainder>1: t_list_1 += [ 
    # leave off remainder for now...
    hist1 = np.histogram(t_list_1,bins=bins1)
    hist2 = np.histogram(t_list_2,bins=bins2)
    hist3 = np.histogram(t_list_3,bins=bins3)
    hist4 = np.histogram(t_list_4,bins=bins4)
    return hist1,hist2,hist3,hist4
  def do_4way_binning(tlist,bins1,bins2,bins3,bins4):
    '''Bin a list of points four ways.
    
    This version does it randomly.'''
    sublist_length = int(np.floor(float(len(tlist))/4.))
    tlist1,tlist2,tlist3,tlist4 = np.random.choice(tlist,size=(4,sublist_length),replace=False)
    remainder = len(tlist)-4*sublist_length
    if remainder>0: print 'INFO: dropping %d values in 4-way binning!'%remainder
    hist1 = np.histogram(tlist1,bins=bins1)
    hist2 = np.histogram(tlist2,bins=bins2)
    hist3 = np.histogram(tlist3,bins=bins3)
    hist4 = np.histogram(tlist4,bins=bins4)
    return hist1,hist2,hist3,hist4
  tlist = fitrecords['time']
  quad_hists = do_4way_binning(
    tlist, 
    bin_edges['low'], 
    bin_edges['med1'], 
    bin_edges['med2'], 
    bin_edges['high']
  )
  bin_counts['low'] = np.array(quad_hists[0][0])
  bin_counts['med1'] = np.array(quad_hists[1][0])
  bin_counts['med2'] = np.array(quad_hists[2][0])
  bin_counts['high'] = np.array(quad_hists[3][0])
  #bin_counts['fitdata'] = bin_counts['low']+bin_counts['med1']+bin_counts['med2']+bin_counts['high']
  bin_counts['fitdata'] = np.histogram(tlist,bins=bin_edges['fitdata'])[0]
  bin_counts['alldata'] = np.histogram(tlist,bins=bin_edges['alldata'])[0]
elif input_filetype=='hist':
  # override previously-computed settings with the ones loaded from file
  binwidth = cluster_hist_data['binwidth']
  half_T_a0 = cluster_hist_data['half_T_a0']
  omega_a0 = cluster_hist_data['omega_a0']
  bin_counts = cluster_hist_data['bin_counts']
  bin_edges = cluster_hist_data['bin_edges']
  fit_tmax = cluster_hist_data['fit_tmax']
  N_bins = cluster_hist_data['N_bins']
  fit_tmin = cluster_hist_data['fit_tmin']
  bin_centers = cluster_hist_data['bin_centers']
else: raise ValueError(input_filetype+' input_filetype?')



if make_plots:
  if verbose: print 'making plots...'
  #trial_counts_plot = plt.figure(figsize=(18,7)).add_subplot(1,1,1,
  #  title='Trial Histogram',xlabel='time [ns]'
  #)
  #trial_counts_plot.set_yscale('log')
  #trial_counts_plot.plot(
  #  bin_centers['trial'], trial_hist_counts, 
  #  marker='', color='black', 
  #  linestyle='-', linewidth=0.75, 
  #  label='Trial Data'
  #)
  #save_figure(trial_counts_plot.figure,'trial_counts')
  #if popup_plots: trial_counts_plot.figure.show()
  data_plot = plt.figure(figsize=(18,7)).add_subplot(1,1,1,title='data',xlabel='t [ns]')
  data_plot.set_yscale('log')
  if any(bin_counts['alldata']):
    data_plot.plot(bin_centers['alldata'],bin_counts['alldata'],
      marker='',linestyle='-',linewidth=0.75,
      label='alldata'
    )
  #data_plot.plot(bin_centers['fitdata'],bin_counts['fitdata'],
  #  marker='',linestyle='-',linewidth=0.75,
  #  label='fitdata'
  #)
  data_plot.plot(bin_centers['med1'],bin_counts['med1'],
    marker='',linestyle='-',linewidth=0.75,
    label='med1'
  )
  data_plot.plot(bin_centers['med2'],bin_counts['med2'],
    marker='',linestyle='-',linewidth=0.75,
    label='med2'
  )
  data_plot.plot(bin_centers['fitdata'],bin_counts['low'],
    marker='',linestyle='-',linewidth=0.75,
    label='low'
  )
  data_plot.plot(bin_centers['fitdata'],bin_counts['high'],
    marker='',linestyle='-',linewidth=0.75,
    label='high'
  )
  data_plot.legend(loc=0)
  save_figure(data_plot,'data_plot')
  if popup_plots: data_plot.figure.show()




#if make_plots:
#  dataset_plot = plt.figure(figsize=(18,7)).add_subplot(1,1,1,title='Data',xlabel='t [ns]')
#  dataset_plot.set_yscale('log')
#  dataset_plot.plot(bin_centers['med'],bin_counts['med1'],label='Medium 1')
#  dataset_plot.plot(bin_centers['med'],bin_counts['med2'],label='Medium 2')
#  dataset_plot.plot(bin_centers['med'],bin_counts['low'],label='Low')
#  dataset_plot.plot(bin_centers['med'],bin_counts['high'],label='High')
#  dataset_plot.legend(loc=0)
#  if popup_plots: dataset_plot.figure.show()


#### compose 'graphs' of U, V, and their linear combinations ###
if verbose: print 'generating U & V graphs...'
graph = {}
graph['t'] = bin_centers['fitdata']
graph['U'] = bin_counts['low'] + bin_counts['high']
#if sample_type in ('induced_gaus','sample1'):
#  graph['V'] = 2.*bin_counts['med']
#elif sample_type in ('sample4',):
#  graph['V'] = bin_counts['med1'] + bin_counts['med2']
graph['V'] = bin_counts['med1'] + bin_counts['med2']
graph['D'] = graph['U'] - graph['V']
graph['S'] = graph['U'] + graph['V']


if make_plots:
  if verbose: print 'plotting U & V graphs...'
  UV_plots = plt.figure(figsize=(18,7)).add_subplot(1,1,1,title='U & V',xlabel='t [ns]')
  UV_plots.set_yscale('log')
  UV_plots.plot(graph['t'],graph['U'],label='U')
  UV_plots.plot(graph['t'],graph['V'],label='V')
  UV_plots.legend(loc=0)
  save_figure(UV_plots,'UV_plots')
  if popup_plots: UV_plots.figure.show()
  
  D_plot = plt.figure(figsize=(18,7)).add_subplot(1,1,1,title='D=U-V',xlabel='t [ns]')
  #D_plot.set_yscale('log')
  D_plot.plot(graph['t'],graph['D'],label='data')
  D_plot.legend(loc=0)
  save_figure(D_plot,'D_plot')
  if popup_plots: D_plot.figure.show()
  
  S_plot = plt.figure(figsize=(18,7)).add_subplot(1,1,1,title='S=U+V',xlabel='t [ns]')
  S_plot.set_yscale('log')
  S_plot.plot(graph['t'],graph['S'],label='data')
  S_plot.legend(loc=0)
  save_figure(S_plot,'S_plot')
  if popup_plots: S_plot.figure.show()



#### compose graph of R and clip it to the fitting domain ###

if fit_tmin_override: fit_tmin = fit_tmin_override
if fit_tmax_override: fit_tmax = fit_tmax_override

R_xlist,R_ylist = [],[] # force copy (instead of reference)
#if use_ROOT==False:
if verbose: print 'composing ratio graph...'
if clip_ratio_fit_at_1st_empty_bin:
  for t,D,S in zip(graph['t'],graph['D'],graph['S']):
    if S>0.:
      R_xlist += [ t ]
      R_ylist += [ float(D)/float(S) ]
    else: break
else:
  raise NotImplementedError('r_i=d_i/0 not handled yet')

# clip if desired..? (TODO: np.array.clip)
if verbose: print 'clipping ratio graph domain (if needed)...'
i_start,i_stop = find_iclip(
  R_xlist, 
  xmin=fit_tmin if fit_tmin>-1 else None, 
  xmax=fit_tmax if fit_tmax>-1 else None
)
R_xlist = R_xlist[i_start:i_stop+1]
R_ylist = R_ylist[i_start:i_stop+1]
graph_clipped = {}
for k in graph.keys():
  graph_clipped[k] = np.array(graph[k][i_start:i_stop+1])
fitdata_bin_centers_clipped = bin_centers['fitdata'][i_start:i_stop+1]
fitdata_bin_counts_clipped = bin_counts['fitdata'][i_start:i_stop+1]

#S_ylist_clipped = graph['S'][i_start:i_stop+1] <-- now graph_clipped['S']

R_xlist = np.array(R_xlist)
R_ylist = np.array(R_ylist)

if make_plots:
  ratio_plot = plt.figure(figsize=(18,7)).add_subplot(1,1,1,title='ratio',xlabel='t [ns]')
  ratio_plot.plot(R_xlist,R_ylist,label='data')
  ratio_plot.legend(loc=0)
  save_figure(ratio_plot,'ratio_plot')
  if popup_plots: ratio_plot.figure.show()


#### estimate initial parameters for fitting ###
if verbose: print 'estimating parameters for fit starting point...'
S_distro_mean = np.sum(graph['S']*graph['t'])/np.sum(graph['S']) # estimates tau_a
S_hist_integral = sum(graph['S'])*binwidth # estimates N0
est_lambda = S_distro_mean - graph['t'][0] # mean estimates decay length for exponential
est_N0 = 2.*S_hist_integral*math.exp(-graph['t'][0]/est_lambda) # 'counts' (+/- scale factors of 2 or 4...)
est_N0 *= 1.25
est_A = 0.5*(max(R_ylist[:len(R_ylist)/10])-min(R_ylist[:len(R_ylist)/10]))
est_phi = estimate_phase(R_xlist[50],R_ylist[50],R_xlist[90],R_ylist[90],omega_a0)

#raise NotImplementedError('...there be drag0ns below...')

#### fitting!###

# fit five-parameter function
fivepar_initial_params = (
  est_N0, 
  est_A, 
  0., 
  est_phi, 
  est_lambda
)
if verbose: print 'performing five-parameter fit...'
fivepar_fit_result = scipy.optimize.leastsq(
  fitres_yerr_from_model, 
  fivepar_initial_params, 
  args=(
    fiveparam_model, 
    (
      #bin_centers['alldata'], 
      #bin_counts['alldata']
      #bin_centers['fitdata'], 
      #bin_counts['fitdata']
      fitdata_bin_centers_clipped, 
      fitdata_bin_counts_clipped
    )
  ), 
  full_output=1
)
if verbose: print '...done.'
#if verbose: print '[debug]',fivepar_fit_result
fivepar_fit_result = fitresult.FitResult(
  fivepar_fit_result,
  fivepar_initial_params,
  ('N0','A','r','phi_a','tau_a'), 
  shrink=False
)
print 'Five-parameter fit:\n'+str(fivepar_fit_result)
#retval += [ fivepar_fit_result ]

if make_plots:
  fiveparam_plot = plt.figure(figsize=(18,7)).add_subplot(1,1,1,title='T histogram',xlabel='t [ns]')
  fiveparam_plot.set_yscale('log')
  #fiveparam_plot.plot(bin_centers['alldata'],bin_counts['alldata'],label='all data')
  #fiveparam_plot.plot(bin_centers['fitdata'],bin_counts['fitdata'],label='fitdata')
  #no: fit_tlist = bin_centers['fitdata']
  #fiveparam_plot.plot(fitdata_bin_centers_clipped,fitdata_bin_counts_clipped,label='fitdata')
  fiveparam_plot.plot(bin_centers['fitdata'],bin_counts['fitdata'],label='data')
  fit_tlist = fitdata_bin_centers_clipped
  fiveparam_plot.plot(
    fit_tlist, [ fiveparam_model(t,fivepar_initial_params) for t in fit_tlist ], 
    color='gray', linestyle=':', 
    label='5-param (initial)'
  )
  fiveparam_plot.plot(
    fit_tlist, [ fiveparam_model(t,fivepar_fit_result.params) for t in fit_tlist ], 
    label='5-param fit'
  )
  fiveparam_plot.legend(loc=0)
  save_figure(fiveparam_plot,'fiveparam_plot')
  if popup_plots: fiveparam_plot.figure.show()


###from ROOT import TGraph,TF1
#import ROOT
##ratio_graph = ROOT.TGraph(N_bins['fitdata'],bin_centers['fitdata'],bin_counts['fitdata'])
#ratio_graph = ROOT.TGraphErrors(
#  N_bins['fitdata'],
#  bin_centers['fitdata'], bin_counts['fitdata'],  # x & y
#  #0., # dx
#  #[0.]*N_bins['fitdata'], # dx
#  #None, 
#  np.zeros(N_bins['fitdata']), 
#  np.sqrt(bin_counts['fitdata']) # dy
#)
#fiveparam_tf1 = ROOT.TF1(
#  'fiveparam_tf1', 
#  '([0]/[4])*TMath::Exp(-x/[4])*(1.+[1]*TMath::Cos(%.16f*(1.+[2]*1e-06)*x+[3]))'%omega_a0, 
#  #'([0]/[4])*TMath::Exp(-x/[4])*(1.+[1]*TMath::Cos(TMath::RadToDeg(%.16f*(1.+[2]*1e-06)*x+[3])))'%omega_a0, 
#  0.,800000.
#)
#fiveparam_tf1.SetNpx(int(30.*N_bins['fitdata'])) # curve looks wrong without this
#for i in range(5): 
#  fiveparam_tf1.SetParName(i,('N0','A','r','phi_a','tau_a')[i])
#  fiveparam_tf1.SetParameter(i,fivepar_initial_params[i])
#fitresult_ptr = ratio_graph.Fit(fiveparam_tf1,'S')
#minuit_result = fitresult.ROOTFitResult(fitresult_ptr,fivepar_initial_params)

#if make_plots:
#  c1 = ROOT.TCanvas('c1','c1')
#  c1.SetLogy()
#  #ratio_graph.Draw('A*')
#  fiveparam_tf1.Draw()
#  ratio_graph.Draw('SAME')
#  c1.SaveAs('fiveparam_root.png')





# fit exponential graph (U+V)
exp_initial_params = (
  #S_hist_integral, # 'counts' (+/- scale factors of 2 or 4...)
  est_N0, 
  est_lambda
)
if verbose: print 'performing exponential fit...'
exp_fit_result = scipy.optimize.leastsq(
  fitres_yerr_from_model, 
  exp_initial_params, 
  args=(exp_model,(graph_clipped['t'],graph_clipped['S'])), 
  full_output=1
)
if verbose: print '...done.'
exp_fit_result = fitresult.FitResult(
  exp_fit_result,exp_initial_params,('N0','tau_a'),shrink=False
)
print 'Exponential fit:\n'+str(exp_fit_result)

if make_plots:
  fit_tlist = graph_clipped['t']
  S_plot.plot(
    fit_tlist, [ exp_model(t,exp_initial_params) for t in fit_tlist ], 
    color='gray', linestyle=':', 
    label='initial params'
  )
  S_plot.plot(
    fit_tlist, [ exp_model(t,exp_fit_result.params) for t in fit_tlist ], 
    label='fit result params'
  )
  S_plot.legend(loc=0)
  save_figure(S_plot,'S_plot')
  if popup_plots: S_plot.figure.show()


# fit ratio sinusoid
ratio_initial_params = (
  est_A, 
  0., 
  est_phi
)
if verbose: print 'performing sinusoidal fit...'
ratio_fit_result = scipy.optimize.leastsq(
  #rmeth_sin_residuals_finite_err, 
  rmeth_sin_residuals_finite_err_coshterms, 
  ratio_initial_params, 
  args=((R_xlist,R_ylist,graph_clipped['S']),), 
  full_output=1
)
#ratio_fit_result = scipy.optimize.leastsq(
#  fitres_yerr_from_model, 
#  ratio_initial_params, 
#  args=(sine_lsq_model,(R_xlist,R_ylist),), 
#  full_output=1
#)
if verbose: print '...done.'
ratio_fit_result = fitresult.FitResult(
  ratio_fit_result,ratio_initial_params,('A','r','phi_a'),shrink=False
)
print 'Ratio fit:\n'+str(ratio_fit_result)

if make_plots:
  fit_tlist = R_xlist
  ratio_plot.plot(
    #fit_tlist, [ sine_lsq_model(t,ratio_fit_result.params) for t in fit_tlist ], 
    fit_tlist, [ sine_lsq_model_coshterms(t,ratio_initial_params) for t in fit_tlist ], 
    color='gray', linestyle=':', 
    label='initial params'
  )
  ratio_plot.plot(
    #fit_tlist, [ sine_lsq_model(t,ratio_fit_result.params) for t in fit_tlist ], 
    fit_tlist, [ sine_lsq_model_coshterms(t,ratio_fit_result.params) for t in fit_tlist ], 
    label='fit result params'
  )
  ratio_plot.legend(loc=0)
  save_figure(ratio_plot,'ratio_plot')
  if popup_plots: ratio_plot.figure.show()


# summary data (in case we want to do anything with it)
runsettings = {}
for thing in (
    'verbose',
    'omega_a0','make_plots','save_fits','save_fits_dir','popup_plots',
    'clip_ratio_fit_at_1st_empty_bin',
    'fit_tmin','fit_tmax','binwidth'
  ):
  runsettings[thing] = eval(thing)
savedata = {}
for thing in (
    'runsettings',
    'bin_centers','bin_counts','N_bins',
    'graph','graph_clipped',
    'R_xlist','R_ylist','i_start','i_stop',
    'fivepar_fit_result','exp_fit_result','ratio_fit_result',
  ):
  savedata[thing] = eval(thing)


if save_fits:
  savedata_outfilename = os.path.join(save_fits_dir,'savedata.pkl')
  print 'saving data and run settings to %s...'%savedata_outfilename
  pickle.dump(savedata,open(savedata_outfilename,'wb'))

#---------------------------------------
if not make_plots: sys.exit(0) # done
#---------------------------------------

fivepar_residuals_figure = plot_scaled_residuals(
  fitdata_bin_centers_clipped, fivepar_fit_result.residuals, 
  title='5 Parameter Fit Residuals', 
  avg_scale=(20,100,500)
)
save_figure(fivepar_residuals_figure,'fivepar_residuals')
if popup_plots: fivepar_residuals_figure.show()


exp_residuals_figure = plot_scaled_residuals(
  graph_clipped['t'], exp_fit_result.residuals, 
  title='Exponential Fit Residuals', 
  avg_scale=(20,100,500)
)
save_figure(exp_residuals_figure,'exp_residuals')
if popup_plots: exp_residuals_figure.show()


ratio_residuals_figure = plot_scaled_residuals(
  R_xlist, ratio_fit_result.residuals,
  title='Ratio Fit Residuals', 
  avg_scale=(20,100,500)
)
save_figure(ratio_residuals_figure,'ratio_residuals')
if popup_plots: ratio_residuals_figure.show()


residuals_comparison = plt.figure(figsize=(14,8)).add_subplot(1,1,1,
  title='Fit Residuals Comparison', 
  xlabel='time [ns]', ylabel='Fit Residuals (Error-normalized)'
)
residuals_comparison.plot(
  fitdata_bin_centers_clipped, 
  running_average(fivepar_fit_result.residuals,100),
  label='5-param fit'
)
residuals_comparison.plot(
  graph_clipped['t'], 
  running_average(exp_fit_result.residuals,100),
  label='exponential fit'
)
residuals_comparison.plot(
  R_xlist, 
  running_average(ratio_fit_result.residuals,100),
  label='ratio fit'
)
residuals_comparison.axhline(0.,color='gray',linestyle=':')
residuals_comparison.legend(loc=0)
save_figure(residuals_comparison,'residuals_comparison')
if popup_plots: residuals_comparison.figure.show()




fivepar_unscaled_residuals_figure = plot_unscaled_residuals(
  fitdata_bin_centers_clipped, 
  fivepar_fit_result.residuals, 
  errors=np.sqrt(fitdata_bin_counts_clipped), 
  title='5 Parameter Fit Residuals', 
  avg_scale=(20,100,500)
)
save_figure(fivepar_unscaled_residuals_figure,'fivepar_unscaled_residuals')
if popup_plots: fivepar_unscaled_residuals_figure.show()


exp_unscaled_residuals_figure = plot_unscaled_residuals(
  graph_clipped['t'], 
  exp_fit_result.residuals, 
  errors=np.sqrt(graph_clipped['S']), 
  title='Exponential Fit Residuals', 
  avg_scale=(20,100,500)
)
save_figure(exp_unscaled_residuals_figure,'exp_unscaled_residuals')
if popup_plots: exp_unscaled_residuals_figure.show()


ratio_unscaled_residuals_figure = plot_unscaled_residuals(
  R_xlist, 
  ratio_fit_result.residuals, 
  errors=np.sqrt( (1. - R_ylist**2) / graph_clipped['S'] ), 
  title='Ratio Fit Residuals', 
  avg_scale=(20,100,500)
)
save_figure(ratio_unscaled_residuals_figure,'ratio_unscaled_residuals')
if popup_plots: ratio_unscaled_residuals_figure.show()



fiveparam_residuals_fft_plot = plot_residuals_fft(
  fitdata_bin_centers_clipped, 
  fivepar_fit_result.residuals*np.sqrt(fitdata_bin_counts_clipped), 
  t_min_ns=ft_tmin, 
  t_max_ns=ft_tmax, 
  title='Five-parameter fit residuals FT'
)
save_figure(fiveparam_residuals_fft_plot,'fivepar_residuals_ft')
if popup_plots: fiveparam_residuals_fft_plot.show()


exp_residuals_fft_plot = plot_residuals_fft(
  graph_clipped['t'], 
  exp_fit_result.residuals*np.sqrt(graph_clipped['S']), 
  t_min_ns=ft_tmin, 
  t_max_ns=ft_tmax, 
  title='Exponential fit residuals FT'
)
save_figure(exp_residuals_fft_plot,'exp_residuals_ft')
if popup_plots: exp_residuals_fft_plot.show()


ratio_residuals_fft_plot = plot_residuals_fft(
  R_xlist, ratio_fit_result.residuals*np.sqrt( (1. - R_ylist**2) / graph_clipped['S'] ), 
  t_min_ns=ft_tmin, 
  t_max_ns=ft_tmax, 
  title='Ratio fit residuals FT'
)
save_figure(ratio_residuals_fft_plot,'ratio_residuals_ft')
if popup_plots: ratio_residuals_fft_plot.show()



# test extraction of decay/pileup/loss rates via (U+V)/exp
# (approximations apply...)
#graph['t'],graph['S']
#rate_extract = 149.2*4*graph['S']/(exp_fit_result.params[0]*np.exp(-graph['t']/exp_fit_result.params[1]))
rate_extract = graph['S']/(exp_fit_result.params[0]*np.exp(-graph['t']/exp_fit_result.params[1])/exp_fit_result.params[1])
rate_extract_plot = plt.figure(figsize=(17,9)).add_subplot(1,1,1,title='Rate Extract',xlabel='time [ns]',ylabel='(U+V)/exp')
rate_extract_plot.figure.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.93)
rate_extract_plot.plot(graph['t'],rate_extract,marker='.',linestyle='',color='black',markersize=1.5)
rate_extract_plot.plot(graph['t'],running_average(rate_extract,100))
rate_extract_plot.plot(graph['t'],running_average(rate_extract,15))
rate_extract_plot.axhline(1.,linestyle=':',color='gray')
rate_extract_plot.axis(ymin=0.997,ymax=1.010)
rate_extract_plot.figure.show()

# for fitting inferred rate with a polynomial
#N = 3
#def make_coeffs_array(N,omega_a):
#  half_T_a = math.pi/omega_a
#  return np.array( 
#    [ 
#      [ 
#        np.nan if k>l 
#        else k*np.log(half_T_a) \
#          + np.log(math.factorial(l)/math.factorial(k)/math.factorial(l-k)) 
#        for k in range(N+1) 
#      ] 
#      for l in range(N+1) 
#    ] 
#  )
#print 'rate extraction polynomial coefficients:'
#fivepar_fit_omega_a = omega_a0*(1.+fivepar_fit_result.params[2]*1e-06)
#exp_fit_tau = exp_fit_result.params[1]
#polynom_coeffs_array = make_coeffs_array(N,fivepar_fit_omega_a)
##for row in polynom_coeffs_array:
##  print '  '.join(['%9.3g'%x for x in row])
#def array_print(arr,fmt='9.3g'):
#  for row in arr: print ' '.join([('%'+fmt)%x for x in row])
#array_print(polynom_coeffs_array)

#def make_tp_array(t,params):
#  '''Array of log(params[k]) + (l-k)*log(t).'''
#  # for a given t and params:
#  log_t = np.log(t)
#  log_params = np.log(params)
#  return np.array([ [ np.nan if k>l else log_params[k] + (l-k)*log_t for k in range(N+1) ] for l in range(N+1) ])
#def make_dcpt_array(t,params,debug=False):
#  '''Array of d_l*(l choose k)*half_period^k*t^(l-k).
#  
#  Elements with k>l are zero.
#  '''
#  dpct_array = make_tp_array(t,params) + polynom_coeffs_array
#  if debug: 
#    for row in dpct_array: print ' '.join(['%9.3g'%x for x in row])
#  return dpct_array

#test_t = 20000.
#test_params = range(1,N+2)
#test_params.reverse()
#print 'test_t: %f\ntest_params: %s'%(test_t,test_params)
#exp_dpct_array = np.exp(make_dcpt_array(test_t,test_params))
##for row in exp_dpct_array: print '  '.join(['%9.3g'%x for x in row])
#array_print(exp_dpct_array)

#dl_sum_even,dl_sum_odd = 0.,0.
#for l in range(0,N+1,2): dl_sum_even += sum(exp_dpct_array[l,:l+1])
#for l in range(1,N+1,2): dl_sum_odd += sum(exp_dpct_array[l,:l+1])

#def exp_cosh_term(t,params,omega_a,tau,polynom_coeffs_array):
#  exp_dpct_array = np.exp(make_dcpt_array(t,params))
#  dl_sum_even,dl_sum_odd = 0.,0.
#  for l in range(0,N+1,2): dl_sum_even += sum(exp_dpct_array[l,:l+1])
#  for l in range(1,N+1,2): dl_sum_odd += sum(exp_dpct_array[l,:l+1])
#  return 2.*np.exp(-dl_sum_even)*np.cosh(dl_sum_odd+math.pi/tau/omega_a)
#term = exp_cosh_term(test_t,test_params,fivepar_fit_omega_a,exp_fit_tau,polynom_coeffs_array)
#print term



def polynom_model(t,params,t_scale):
  '''Polynomial with form p[i]*(t/t_scale)^i.
  
  The fitter should be less confused by the (typically) order-of-magnitude
  differences between the lowest and highest coefficient.
  '''
  return np.sum(
    [ p*(t/t_scale)**l for l,p in enumerate(params) ]
  )
def erd_model(t,params,t_scale,half_Ta,tau):
  '''Exponential rate deviation model (approximately ~ S/exp(-t/tau)).
  
  exp(-d(t-half_Ta)+half_Ta/tau) + exp(-d(t+half_Ta)-half_Ta/tau + 2exp(-d(t))
    ~= (U+V)/(4*N0*exp(-t/tau)/tau)
  
  This neglects the cos() term (which should be small).
  
  d(t) is a polynomial with form p[i]*(t/t_scale)^i.
  
  The fitter should be less confused by the (typically) order-of-magnitude
  differences between the lowest and highest coefficient.
  '''
  d = polynom_model(t,params,t_scale)
  d_minus = polynom_model(t-half_Ta,params,t_scale)
  d_plus = polynom_model(t+half_Ta,params,t_scale)
  retval = np.exp(-d_minus+half_Ta/tau) + np.exp(-d_plus-half_Ta/tau) + 2.*np.exp(-d)
  #print t,d,d_minus,d_plus,retval
  return retval

new_i_stop = 0
for new_i_stop in range(len(bin_counts['fitdata'])):
  if bin_counts['fitdata'][new_i_stop]<=0: break
print 'found first zero in fitdata at %d'%new_i_stop

if False:
  d_initial_params = (1.,0.1,0.001,0.0001)
  d_fit_result = scipy.optimize.leastsq(
    fitres_yerr, 
    d_initial_params, 
    args=(
      #erd_model,
      lambda t,params: 4.*np.exp(-polynom_model(t,params,700000.)), 
      (graph['t'][:new_i_stop],rate_extract[:new_i_stop],np.sqrt(bin_counts['fitdata'][:new_i_stop])),
      #700000., half_T_a0, exp_fit_result.params[1]
    ), 
    full_output=1
  )
  d_fit_result = fitresult.FitResult(d_fit_result,d_initial_params,('d0','d1','d2','d3'))
  if d_fit_result.success:
    rate_extract_plot.plot(
      graph['t'][:new_i_stop], 
      [ 
        4.*np.exp(-polynom_model(t,d_fit_result.params,700000.))
        for t in graph['t'][:new_i_stop] 
      ]
    )
    rate_extract_plot.figure.show()
elif False:
  d_initial_params = (1.,0.1,0.001,0.0001)
  d_fit_result = scipy.optimize.leastsq(
    fitres, 
    d_initial_params, 
    args=(
      erd_model,
      (graph['t'][:new_i_stop],rate_extract[:new_i_stop]),
      700000., half_T_a0, exp_fit_result.params[1]
    ), 
    full_output=1
  )
  d_fit_result = fitresult.FitResult(d_fit_result,d_initial_params,('d0','d1','d2','d3'))
  if d_fit_result.success:
    rate_extract_plot.plot(
      graph['t'][:new_i_stop], 
      [ 
        erd_model(t,d_fit_result.params,700000.,half_T_a0, exp_fit_result.params[1])
        for t in graph['t'][:new_i_stop] 
      ]
    )
    rate_extract_plot.figure.show()
elif False: # WORKING!
  d_initial_params = (1.,0.1,0.001,0.0001,0.00001)
  d_fit_result = scipy.optimize.leastsq(
    fitres_yerr, 
    d_initial_params, 
    args=(
      erd_model,
      (
        graph['t'][:new_i_stop],
        rate_extract[:new_i_stop],
        0.05*60.*np.exp(graph['t'][:new_i_stop]*1e-05)*1e-05
      ),
      700000., half_T_a0, exp_fit_result.params[1]
    ), 
    full_output=1
  )
  d_fit_result = fitresult.FitResult(d_fit_result,d_initial_params,('d0','d1','d2','d3'))
  if d_fit_result.success:
    rate_extract_plot.plot(
      graph['t'][:new_i_stop], 
      [ 
        erd_model(t,d_fit_result.params,700000.,half_T_a0, exp_fit_result.params[1])
        for t in graph['t'][:new_i_stop] 
      ]
    )
    rate_extract_plot.figure.show()
elif True:
  d_order = 5
  d_fit_i_start = 0
  d_fit_i_stop = new_i_stop
  d_initial_params = np.power(10.,range(0,-1-d_order,-1))
  d_fit_result = scipy.optimize.leastsq(
    fitres_yerr, 
    d_initial_params, 
    args=(
      erd_model,
      (
        graph['t'][d_fit_i_start:d_fit_i_stop],
        rate_extract[d_fit_i_start:d_fit_i_stop],
        # first guess at errors (emperical): 
        # 0.05*60.*np.exp(graph['t'][d_fit_i_start:d_fit_i_stop]*1e-05)*1e-05
        # first computation of errors / sqrt(149.2) / 4:
        2.*math.sqrt(exp_fit_result.params[1]/exp_fit_result.params[0]/binwidth)*np.exp(graph['t'][d_fit_i_start:d_fit_i_stop]/2./exp_fit_result.params[1])/4.
      ),
      700000., half_T_a0, exp_fit_result.params[1]
    ), 
    full_output=1
  )
  d_fit_result = fitresult.FitResult(
    d_fit_result, 
    d_initial_params, 
    #('d0','d1','d2','d3')
    tuple([ 'p'+str(x) for x in range(5+1) ])
  )
  if d_fit_result.success:
    rate_extract_plot.plot(
      graph['t'][d_fit_i_start:d_fit_i_stop], 
      [ 
        erd_model(t,d_fit_result.params,700000.,half_T_a0, exp_fit_result.params[1])
        for t in graph['t'][d_fit_i_start:d_fit_i_stop] 
      ]
    )
    rate_extract_plot.axis(ymin=0.997,ymax=1.010)
    rate_extract_plot.figure.show()
    save_figure(rate_extract_plot.figure,'rate_extraction')

N_points = len(graph['t'][d_fit_i_start:d_fit_i_stop])
N_params = len(d_fit_result.params)
ymin_list = np.array([None]*N_points)
ymax_list = np.array([None]*N_points)

#for i in range(N_points):
#  ymin_list[i] = None
#  y


#def make_series(width,depth,x=None,inc=0):
#  if x==None: x = range(width)
#  else: x = [x]*width
#  if inc>=depth: return x
#  else: return make_series(width,depth,x,inc+1)



