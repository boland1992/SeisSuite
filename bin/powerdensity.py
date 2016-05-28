# -*- coding: utf-8 -*-
"""
Created on Fri July 6 11:04:03 2015

@author: boland
"""

import os
import glob
import scipy
import datetime
import numpy as np
import datetime as dt
import multiprocessing as mp
import matplotlib.pyplot as plt
from numpy.lib.stride_tricks import as_strided
from numpy.fft import rfft, irfft
from obspy import read_inventory
from scipy import signal
from obspy import read
import sqlite3 as lite
import pickle

from obspy.core import UTCDateTime as utc
from pyseis.modules.rdreftekc import rdreftek, reftek2stream
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
DATALESS_DIR = CONFIG.DATALESS_DIR
STATIONXML_DIR = CONFIG.STATIONXML_DIR
CROSSCORR_DIR = CONFIG.CROSSCORR_DIR
USE_DATALESSPAZ = CONFIG.USE_DATALESSPAZ
USE_STATIONXML = CONFIG.USE_STATIONXML
CROSSCORR_STATIONS_SUBSET = CONFIG.CROSSCORR_STATIONS_SUBSET
CROSSCORR_SKIPLOCS = CONFIG.CROSSCORR_SKIPLOCS
FIRSTDAY = CONFIG.FIRSTDAY
LASTDAY = CONFIG.LASTDAY




TIMELINE_DB = os.path.join(DATABASE_DIR, 'timeline.db')

PSD_OUTPUT = os.path.join(CROSSCORR_DIR, 'PSD')
if not os.path.exists(PSD_OUTPUT): os.mkdir(PSD_OUTPUT)
    
t_start, t_end = utc(FIRSTDAY), utc(LASTDAY)

def read_ref(path):
    ref_head, ref_data = rdreftek(path)
    st = reftek2stream(ref_head, ref_data)
    return st
    
#------------------------------------------------------------------------------
# IMPORT PATHS TO MSEED FILES
#------------------------------------------------------------------------------

def getpaths(database_loc, starttime, endtime):
        """
        Gets list of paths to mseed files between certain time series intervals
        using initialised SQL timeline database.
        @type starttime: L{UTCDateTime} or L{datetime} or L{date}
        @type endtime: L{UTCDateTime} or L{datetime} or L{date}

        @rtype: unicode
        """

        starttime, endtime = utc(starttime), utc(endtime)
        
        import_start = starttime.timestamp
        #import_end = endtime.timestamp
        
        #connect SQL database

        if not os.path.exists(database_loc):
            raise Exception("Database doesn't exist")
        
        conn = lite.connect(database_loc)
        #print "conn: ", conn
        
        c = conn.cursor()
        extrema = []
        for row in c.execute('''SELECT * FROM 
                             file_extrema ORDER BY station'''):
            extrema.append(row)
        
        # make sure that this test works! 
        #code = '{}.{}.{}'.format(self.network, self.name, self.channel)
    
        file_paths = c.execute('''SELECT file_path FROM 
                             file_extrema WHERE starttime <= ? 
                             AND endtime >= ?''', 
                             (import_start, import_start))        
        
        output_path = []
        for file_path in file_paths: output_path.append(file_path)
        #if len(output_path) > 0:
        #    return str(output_path[0][0])

        
        # close database
        conn.close() 
        return output_path



def spectrum(tr):
    wave = tr.data #this is how to extract a data array from a mseed file
    fs = tr.stats.sampling_rate
    #hour = str(hour).zfill(2) #create correct format for eqstring
    f, Pxx_spec = signal.welch(wave, fs)#, 'flattop', 
                                     #nperseg=1024, scaling='spectrum')
                    #plt.semilogy(f, np.sqrt(Pxx_spec))
    
    return f, Pxx_spec
    
    #if len(f) >= 256:
    #    return np.column_stack((f[:255], Pxx_spec[:255]))
        #column = np.column_stack((f[:255], np.abs(np.sqrt(Pxx_spec)[:255])))   
        #return column
    #else:
    #    return 0.




fig1 = plt.figure(figsize=(15,10))
ax1 = fig1.add_subplot(111)    
ax1.set_title("Average Seismic Power Spectral Density")

ax1.set_xlabel('Frequency [Hz]')
ax1.set_ylabel('Power Density Spectrum [dB/Hz]')

ax1.set_xlim([0,100])
ax1.grid(True, axis='both', color='gray')
ax1.set_autoscaley_on(True)
#ax1.set_yscale('log')

# search through every day period in time window

print "Processing average power spectral densities between {} and {} ..."\
.format(t_start, t_end)

 


    
    
# get only the miniseed paths that have data between t_start and t_end
stream_paths = getpaths(TIMELINE_DB, t_start, t_end)
    

for stream_path in stream_paths: 
        
        try:

            stream_path = str(stream_path[0])
            print "\nCurrently processing file: {}"\
            .format(stream_path)

            st = read(stream_path, starttime=t_start, endtime=t_end)
            
            for tr in st:
                f, Pxx_spec = spectrum(tr)
                
                # convert Pxx_spec into decibels/Hz
                
                Pxx_spec = 10 * np.log10(Pxx_spec)
            
                plt.plot(f, Pxx_spec, c='k', alpha=0.1)
        
        except Exception as error:
            print error



fig1.savefig(os.path.join(PSD_OUTPUT, 'PDS_avg.svg'), 
             format='svg', dpi=300)  
plt.clf


quit()