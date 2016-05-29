# -*- coding: utf-8 -*-
"""
Created on Fri May 27 08:24:22 2016

@author: iese
"""


import os
import warnings
import datetime as dt
import itertools as it
import pickle
import obspy.signal.cross_correlation
import time
import glob
import sqlite3 as lite
import shutil
import numpy as np
import matplotlib.pyplot as plt
from obspy import read
import datetime as dt
from obspy.core import UTCDateTime as utc

def get_filepaths(directory):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.    
    
    
    
    
# set epoch timestamp 
epoch = dt.datetime(1970, 1, 1)

total_verbose = True
psd = False

# create a list of configuration files that will be iterated over! 
# must be 1 or more
from seissuite.ant.psconfig import (create_config_list, run_config, 
                                    remove_config)

config_list = create_config_list()

total_time0 = dt.datetime.now()
for config_file in config_list:
    # global variables MUST be defined 
    # with the function in the seissuite.ant.psconfig module 
    run_config(config_file)
  
  
    from seissuite.ant import (pscrosscorr, psstation, pspreprocess, pserrors, 
                               psstationSQL)

    # import CONFIG class initalised in ./configs/tmp_config.pickle
    config_pickle = 'configs/tmp_config.pickle'
    f = open(name=config_pickle, mode='rb')
    CONFIG = pickle.load(f)
    f.close()
    
    # import variables from initialised CONFIG class.
    MSEED_DIR = CONFIG.MSEED_DIR
    DATABASE_DIR = CONFIG.DATABASE_DIR
    
    abs_paths = get_filepaths(MSEED_DIR)
    #abs_paths = ['/home/iese/Documents/local_data/SAMOA/INPUT/DATA/Stations/AFI/Data/2005/335/IU_AFI_BHZ_10_2005_335.msd']

    for path in abs_paths:
        
        try:
            st = read(path, headonly=True)
            times = []
            for tr in st:
            	
                start_time = tr.stats.starttime
		# reduce the first day to the beginning time i.e. midnight. 
		start_time = start_time - (start_time.hour * 3600) - (start_time.minute * 60) - start_time.second
		start_time = start_time.timestamp
                end_time = tr.stats.endtime.timestamp #tr.stats.endtime
		times.append(start_time)
		times.append(end_time)
		
            days = int((end_time - start_time)/86400) + 1
		
            time_intervals = [utc(i) for i in np.linspace(min(times), max(times), days)]
		
	    for i in range(1, len(time_intervals)):
		starttime = time_intervals[i-1]
		endtime = time_intervals[i]

		st_partial = read(path, starttime=starttime, endtime=endtime)
		net = st_partial[0].stats.network
		stat = st_partial[0].stats.station
		loc = st_partial[0].stats.location
		channel = st_partial[0].stats.channel
		year = starttime.year
		jul_day = starttime.julday
		write_string = '{}_{}_{}_{}_{}_{}.msd'.format(net, stat, channel, loc, year, jul_day)
		print write_string
		mseed_write = os.path.join(os.path.dirname(path), write_string)
		st_partial.write(mseed_write, format='MSEED')

    



			
				
		
                



                
                
                
        except Exception as error:
            print error

	os.remove(path)
        

          
quit()
