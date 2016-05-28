#!/usr/bin/python -u
"""
[Advice: run this script using python with unbuffered output:
`python -u dispersion_curves.py`]

This script implements the two-step frequency-time analysis
(FTAN) on a set of cross-correlations, as described by Levshin &
Ritzwoller, "Automated detection, extraction, and measurement
of regional surface waves", Pure Appl. Geoph. (2001) and Bensen
et al., "Processing seismic ambient noise data to obtain
reliable broad-band surface wave dispersion measurements",
Geophys. J. Int. (2007).

In short, the envelope's amplitude of the (analytic representation
of the) cross-correlation is calculated and displayed after applying
narrow bandpass filters: this gives a 2D image, in which the x-axis
is the nominal period (filter's center period) or the instanteneous
period and the  y-axis is time (equivalent to velocity). Tracking
the time (equivalent to velocity) at which the the amplitude reaches
its maximum for each  period gives the dispersion curve, i.e., the
group velocity function  of period. This 'raw' dispersion curve is
used to set up a phase-matched filter, re-apply a 'clean' FTAN and
extract a 'clean' dispersion curve.

The nominal frequency (or period) -- i.e., the filter's center
frequency -- is replaced with the instantaneous frequency
if parameter *USE_INSTANTANEOUS_FREQ* is set to True, or if
use_inst_freq=True is explicitely passed to function xc.FTANs().

This script takes as input one or several binary files containing a
set of cross-correlations (previously calculated with, e.g., script
crosscorrelation.py), located in folder *CROSSCORR_DIR*. A set of
cross-correlations is an instance of pscrosscorr.CrossCorrelationCollection
exported in binary format with module pickle. Two file per set of
cross-correlations are produced:

- a pdf file illustrating the FTAN procedure on all cross-correlations:
  one page per cross-correlation, containing the original and
  band-passed cross-correlation, the amplitude and dispersion curve
  of the raw FTAN, of the clean FTAN, and a map locating the pair
  of stations.

- a binary file containing the clean dispersion curves exported
  with module pickle
"""

import glob
import os
import math
import itertools as it
try:
    import dill as pickle
except:    
    import pickle

import datetime as dt
#import matplotlib.pyplot as plt
import numpy as np

# set numpy to ignore all errors at the moment!
#np.seterr(all='ignore')

#from matplotlib.backends.backend_pdf import PdfPages

# parsing configuration file to import dir of cross-corr results
from seissuite.ant import (pscrosscorr, psutils)

process_type = "serial"

multiprocess = False
#import multiprocessing as mp
if multiprocess:
    
    try:
        import pathos.multiprocessing as mp
    except:
        import multiprocessing as mp
    
    no_of_cores = int(mp.cpu_count())

    process_type = "parallel"



# import CONFIG class initalised in ./configs/tmp_config.pickle
config_pickle = 'configs/tmp_config.pickle'
f = open(name=config_pickle, mode='rb')
CONFIG = pickle.load(f)
f.close()
    
# import variables from initialised CONFIG class.
MSEED_DIR = CONFIG.MSEED_DIR
DATABASE_DIR = CONFIG.DATABASE_DIR
DATALESS_DIR = CONFIG.DATALESS_DIR
STATIONXML_DIR = CONFIG.STATIONXML_DIR
CROSSCORR_DIR = CONFIG.CROSSCORR_DIR
USE_DATALESSPAZ = CONFIG.USE_DATALESSPAZ
USE_STATIONXML = CONFIG.USE_STATIONXML
CROSSCORR_STATIONS_SUBSET = CONFIG.CROSSCORR_STATIONS_SUBSET
CROSSCORR_SKIPLOCS = CONFIG.CROSSCORR_SKIPLOCS
FIRSTDAY = CONFIG.FIRSTDAY
LASTDAY = CONFIG.LASTDAY
MINFILL = CONFIG.MINFILL
FREQMIN = CONFIG.FREQMIN
FREQMAX = CONFIG.FREQMAX
CORNERS = CONFIG.CORNERS
ZEROPHASE = CONFIG.ZEROPHASE
PERIOD_RESAMPLE = CONFIG.PERIOD_RESAMPLE
ONEBIT_NORM = CONFIG.ONEBIT_NORM
FREQMIN_EARTHQUAKE = CONFIG.FREQMIN_EARTHQUAKE
FREQMAX_EARTHQUAKE = CONFIG.FREQMAX_EARTHQUAKE
WINDOW_TIME = CONFIG.WINDOW_TIME
WINDOW_FREQ = CONFIG.WINDOW_FREQ
XCORR_INTERVAL = CONFIG.XCORR_INTERVAL
CROSSCORR_TMAX = CONFIG.CROSSCORR_TMAX
CROSSCORR_DIR = CONFIG.CROSSCORR_DIR
FTAN_DIR = CONFIG.FTAN_DIR
PERIOD_BANDS = CONFIG.PERIOD_BANDS
CROSSCORR_TMAX = CONFIG.CROSSCORR_TMAX
PERIOD_RESAMPLE = CONFIG.PERIOD_RESAMPLE
SIGNAL_WINDOW_VMIN = CONFIG.SIGNAL_WINDOW_VMIN
SIGNAL_WINDOW_VMAX = CONFIG.SIGNAL_WINDOW_VMAX
SIGNAL2NOISE_TRAIL = CONFIG.SIGNAL2NOISE_TRAIL
NOISE_WINDOW_SIZE = CONFIG.NOISE_WINDOW_SIZE
RAWFTAN_PERIODS = CONFIG.RAWFTAN_PERIODS
CLEANFTAN_PERIODS = CONFIG.CLEANFTAN_PERIODS
FTAN_VELOCITIES = CONFIG.FTAN_VELOCITIES
FTAN_ALPHA = CONFIG.FTAN_ALPHA
STRENGTH_SMOOTHING = CONFIG.STRENGTH_SMOOTHING
USE_INSTANTANEOUS_FREQ = CONFIG.USE_INSTANTANEOUS_FREQ
MAX_RELDIFF_INST_NOMINAL_PERIOD = CONFIG.MAX_RELDIFF_INST_NOMINAL_PERIOD
MIN_INST_PERIOD = CONFIG.MIN_INST_PERIOD
HALFWINDOW_MEDIAN_PERIOD = CONFIG.HALFWINDOW_MEDIAN_PERIOD
MAX_RELDIFF_INST_MEDIAN_PERIOD = CONFIG.MAX_RELDIFF_INST_MEDIAN_PERIOD
BBOX_LARGE = CONFIG.BBOX_LARGE
BBOX_SMALL = CONFIG.BBOX_SMALL
FIRSTDAY = CONFIG.FIRSTDAY
LASTDAY = CONFIG.LASTDAY

# loading cross-correlations (looking for *.pickle files in dir *CROSSCORR_DIR*)
folder_list = sorted(glob.glob(os.path.join(CROSSCORR_DIR, '*')))

pickle_list = []
if len(folder_list) < 1: 
    print("There are no files or folders in the data input folder \
please re-run 02_timeseries_process.py to get some results")

else: 
    for folder in folder_list:
    #check to see if there are any pickle files in the xcorr time folder 
        if len(glob.glob(os.path.join(folder, '*.pickle'))) < 1:
            #print("There are no .pickle files in this folder. Skipping ...")
            continue
        else:
            for file_ in glob.glob(os.path.join(folder, '*.pickle')):
                if 'metadata' not in file_ and '.part' not in file_:            
                    pickle_list.append(file_)
                
if len(pickle_list) < 1: 
    print("\nThere are no pickle files to begin from.")
    raise Exception("No pickle files to process, first run the programme.")
    res = ""
        
else:
    print "\nPlease choose a file to process." 
    #print combinations of partial pickle files available
    print '\n'.join('{} - {}'.format(i + 1, f.split('/')[-2])
        for i, f in enumerate(pickle_list))
        
    #change folder_list to pickle_list if this gives problems
    res = raw_input('\n')
    
del f
#create list of pickle files to process FTAN for
if not res or res == "0":
    pickle_files = [f for f in pickle_list if f[-1] != '~']
else:
    pickle_files = [pickle_list[int(i)-1] for i in res.split()]

#usersuffix = raw_input("\nEnter suffix to append: [none]\n").strip()
usersuffix = ""

# processing each set of cross-correlations
for pickle_file in pickle_files:
    print "\nOpening pickle file ... " #+ os.path.basename(pickle_file)
    file_opent0 = dt.datetime.now()
    global xc
    xc = pscrosscorr.load_pickled_xcorr(pickle_file)
    delta = (dt.datetime.now() - file_opent0).total_seconds()
    print "\nThe file took {:.1f} seconds to open.".format(delta)
    del delta
    # copying the suffix of cross-correlations file
    # (everything between 'xcorr_' and the extension)
    suffix = os.path.splitext(os.path.basename(pickle_file))[0].replace('xcorr_', '')
    if usersuffix:
        suffix = '_'.join([suffix, usersuffix])

    # Performing the two-step FTAN, exporting the figures to a
    # pdf file (one page per cross-correlation) and the clean
    # dispersion curves to a binary file using module pickle.
    #
    # The file are saved in dir *FTAN_DIR* (defined in configuration file) as:
    # <prefix>_<suffix>.pdf and <prefix>_<suffix>.pickle
    # You can specify *prefix* as input argument in FTANs(), or leave
    # the function define a default prefix, which will look like:
    # FTAN[_whitenedxc][_mindist=...][_minsSNR=...][_minspectSNR=...] ...
    # [_month-year_month-year]
    #
    # Set whiten=True to whiten the spectrum of the cross-correlations
    # (default is False)
    # Set normalize_ampl=True to normalize FTAN amplitude in plots, so
    # that the max amplitude = 1 at each period (default is True)
    # Set logscale=True to plot log(amplitude^2) instead of amplitude
    # (default is True)
    # Set use_inst_freq=True to replace nominal freq with instantaneous
    # freq (default is given by parameter *USE_INSTANTANEOUS_FREQ*)
    #
    # See other options in the docstring of the function.
    
    FTAN_t0 = dt.datetime.now()
    
    print "Creating FTAN pairs ... " 
    pairs, outputpath = xc.FTAN_pairs()
        

    
    #cd Desktop/Link\ to\ TEST/seissuite/
    def FTAN_pool(pair):
        cleanvg = xc.FTAN_parallel(pair, prefix=None, suffix='', 
              whiten=False, normalize_ampl=True, logscale=True, mindist=None,
              minSNR=None, minspectSNR=None, monthyears=None,
              vmin=SIGNAL_WINDOW_VMIN, vmax=SIGNAL_WINDOW_VMAX,
              signal2noise_trail=SIGNAL2NOISE_TRAIL,
              noise_window_size=NOISE_WINDOW_SIZE,
              savefigs=True, outputpath=os.path.join(FTAN_DIR, 'FTAN'))
        
        return cleanvg

    #pairs = pairs[:45]
    if multiprocess:
        # need to deal with memory clutter by setting max.  pair list size
        max_pairs = 10 # this will only allow a pairs list of len 100 to parse

        if len(pairs) >= max_pairs:
#            n = 0 # this is the counter for the split list parallised loop
            # number of required partial lists to make up total pairs
            cleanvg_total = []
            loop_no =  math.ceil(len(pairs) / float(max_pairs)) 
            print "Processing station pair's partial list FTANs in parallel ... "

            for n in range(loop_no):
                index_min, index_max = n * max_pairs, n * max_pairs + max_pairs
                if index_max > len(pairs):
                    index_max = len(pairs)
                pairs_partial = pairs[index_min:index_max]
                

                pool = mp.Pool(no_of_cores)
                cleanvgcurves = pool.map(FTAN_pool, pairs_partial)
                pool.close()
                pool.join()
                
                cleanvg_total.append(cleanvgcurves)
                
                del cleanvgcurves
                del pairs_partial
            
            cleanvgcurves = list(it.chain(*cleanvg_total))
            del cleanvg_total
                
        
        else:
            print "Processing FTANs of {} station pairs in parallel ... "\
            .format(len(pairs))
            # multiprocessing turned on: one process per station
            pool = mp.Pool(no_of_cores)
            cleanvgcurves = pool.map(FTAN_pool, pairs)
            pool.close()
            pool.join()
    
    else:
        print "Processing FTANs of {} station pairs in serial ... "\
        .format(len(pairs))
        #cleanvgcurves = []
        #for pair in pairs:
        cleanvgcurves = map(FTAN_pool, pairs)
    
    #print "FTAN outputs 1: ", cleanvgcurves
    
    cleanvgcurves = list(cleanvgcurves)
    
    #print "FTAN outputs 2: ", cleanvgcurves

    # convert FTAN_outputs into numpy array
    cleanvgcurves = np.asarray(cleanvgcurves)
    
    # remove None types from cleanvgcurves
    
    cleanvgcurves = cleanvgcurves[cleanvgcurves != np.array(None)]
    
    print cleanvgcurves

    FTAN_delta = (dt.datetime.now() - FTAN_t0).total_seconds()
    
    print "\nIt took {:.1f} seconds to process {} station pairs' FTANS in {}"\
    .format(FTAN_delta, len(pairs), process_type)
    
    
    print "\nFTAN outputs: ", np.asarray(cleanvgcurves)
    #print "FTAN outputs type: ", type(cleanvgcurves)
    #print "FTAN outputs shape: ", cleanvgcurves.shape
    
    
    print "\nExporting FTANs ... " 
    # exporting vg curves to pickle file
    f = psutils.openandbackup(outputpath + '.pickle', mode='wb')
    pickle.dump(cleanvgcurves, f, protocol=2)
    f.close()


    print "\nFTAN calculations completed and saved.\n"
    # delete xc and cleanvgcurves in an attempt to free up memory
    del xc
    del cleanvgcurves
