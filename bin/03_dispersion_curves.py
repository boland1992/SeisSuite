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
import pickle
import datetime as dt
# parsing configuration file to import dir of cross-corr results
from seissuite.ant import pscrosscorr

# import CONFIG class initalised in ./configs/tmp_config.pickle
config_pickle = 'configs/tmp_config.pickle'
f = open(name=config_pickle, mode='rb')
CONFIG = pickle.load(f)
f.close()

# import variables from initialised CONFIG class.
CROSSCORR_DIR = CONFIG.CROSSCORR_DIR

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
    xc = pscrosscorr.load_pickled_xcorr(pickle_file)
    delta = (dt.datetime.now() - file_opent0).total_seconds()
    print "\nThe file took {:.1f} seconds to open.".format(delta)

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

    xc.FTANs(suffix=suffix, whiten=False, normalize_ampl=True, logscale=True)