# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 16:37:02 2015

@author: boland
"""

#!/usr/bin/python -u
"""
The following script has been developed in order to quickly and easily 
plot the SNR with time for various methods of processing. 
"""

from seissuite.ant import pscrosscorr
from seissuite.ant.psconfig import (CROSSCORR_DIR, XCORR_INTERVAL, FIRSTDAY)
import glob
import os
import matplotlib.pyplot as plt
#import itertools as it
from obspy.core import UTCDateTime
import datetime as dt
import numpy as np

def datetime_list(SNRarray):
    """
    Function that outputs a list of N datettime objects between the
    start and end times of a given stream. N is dictated by the length
    of the SNRarray list. 
    """

    nlist = len(SNRarray)

    starttime = UTCDateTime(FIRSTDAY)
    endtime = starttime + XCORR_INTERVAL * 60. * nlist
    
    time_delta = (endtime - starttime)
    starttime = str(starttime)
    endtime = str(endtime)
    # convert to datetime.datetime object        
    starttime = dt.datetime.strptime(str(starttime), 
                                     '%Y-%m-%dT%H:%M:%S.%fZ')
    
    time_space = np.linspace(0, time_delta, nlist)
    
    return [starttime + dt.timedelta(milliseconds=1000*time) 
                  for time in time_space]

# parsing configuration file to import dir of cross-corr results

# loading cross-correlations (looking for *.pickle files in dir *CROSSCORR_DIR*)
folder_list = sorted(glob.glob(os.path.join(CROSSCORR_DIR, '*')))

if len(folder_list) < 1: 
    print("There are no files or folders in the data input folder \
please re-run cross-correlation.py to get some results")
else:
    print("Please select one of the following cross-correlation data \
sets to be be processed: ")
    print '\n0 - All except backups (*~)'
    print '\n'.join('{} - {}'.format(i + 1, os.path.basename(f))
        for i, f in enumerate(folder_list))
    res = raw_input('\n')  

#create list of pickle file names within 
pickle_list = []

for folder in folder_list:
    
    pickle_files = glob.glob(os.path.join(folder, '*.pickle'))
    #check to see if there are any pickle files in the xcorr time folder 
    if len(pickle_files) < 1:
        print("There are no .pickle files in this folder. Skipping ...")
        continue
    else:
        #append name of pickle file path location string to pickle_list
        for pickle_file in pickle_files:
            if 'metadata.pickle' not in pickle_file:
                pickle_list.append(pickle_file)

#create list of pickle files to process FTAN for
if not res or res == "0":
    pickle_files = [f for f in pickle_list if f[-1] != '~']
else:
    pickle_files = [pickle_list[int(i)-1] for i in res.split()]

usersuffix = raw_input("\nEnter suffix to append: [none]\n").strip()


# processing each set of cross-correlations
for pickle_file in pickle_files:
    print "\nProcessing cross-correlations of file: " + pickle_file
    
    # plot SNR vs. time for all station pair xcorr waveforms in dictionary
    xc = pscrosscorr.load_pickled_xcorr(pickle_file)
    
    plt.figure()
    plt.title('SNR vs. Time For XCORR Stacks')
    plt.xlabel('Time (UTC)')
    plt.ylabel('SNR (Max. Signal Amplitude / Noise Standard Deviation)')
    
    for item in xc.items():
        s1 = item[0]
        for s2 in item[1].keys():
            #print '{}-{}'.format(s1, s2)
            try:
                SNRarray = xc[s1][s2].SNR_lin
                timearray = datetime_list(SNRarray)
                
                plt.plot(timearray, SNRarray, alpha=0.5)
            except Exception as err:
                print err
    plt.show()