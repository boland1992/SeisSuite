# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 12:28:15 2016

@author: boland
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 08:35:56 2016

@author: boland
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 23:34:00 2015

@author: boland
"""
from seissuite.ant import (pscrosscorr)
import matplotlib.pyplot as plt
import glob
import os
import pickle
import numpy as np
from obspy.core import UTCDateTime as utc
import datetime as dt
from matplotlib.colors import LogNorm
from matplotlib import colors
#PICKLE_PATH = '/storage/ANT/PROGRAMS/ANT_OUTPUT/OUTPUT/CROSS/06.05.2015-15:53:28/XCORR-STACK_01.01.2014-31.12.2014_datalesspaz.pickle'
#PICKLE_PATH = '/home/boland/Desktop/XCORR-STACK_01.08.1999-10.06.2000_datalesspaz.part.pickle'
# import CONFIG class initalised in ./configs/tmp_config.pickle
config_pickle = 'configs/tmp_config.pickle'
f = open(name=config_pickle, mode='rb')
CONFIG = pickle.load(f)
f.close()
    
# import variables from initialised CONFIG class.
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
PLOT_CLASSIC = CONFIG.PLOT_CLASSIC
PLOT_DISTANCE = CONFIG.PLOT_DISTANCE
MAX_DISTANCE = CONFIG.MAX_DISTANCE

pickle_list = []
folder_list = sorted(glob.glob(os.path.join(CROSSCORR_DIR, '*')))


print MSEED_DIR
print CROSSCORR_TMAX

# station list to plot SNR for
# can only plot one station combination per plot!

stat_list = [['AUBSH', 'AUTOO'], ['AUBSH', 'AUWSH'], ['AUDHS', 'AUKAT'], 
             ['AUMOU', 'AUNRC'], ['AUMOU', 'AUTOO'], ['AUNRC', 'AUTOO'],
             ['AUTOO', 'AUWSH']]

for folder in folder_list:
    #check to see if there are any pickle files in the xcorr time folder 
    if len(glob.glob(os.path.join(folder, '*.pickle'))) < 1:
        #print("There are no .pickle files in this folder. Skipping ...")
        continue
    else:
        for file_ in glob.glob(os.path.join(folder, '*.pickle')):
            if 'metadata' not in file_ and '.part' not in file_:            
                pickle_list.append(file_)
                
start_plot = dt.datetime(2014, 1, 1, 00, 00)
end_plot = dt.datetime(2014, 2, 1, 00, 00)
dpi = 300
counter = 0
total_time = []
total_SNR = []

for pair in stat_list:        
    
    s1, s2 = pair
    fig = plt.figure(figsize=(21.0, 12.0))

    for PICKLE_PATH in pickle_list:    
        OUTFILESPATH = PICKLE_PATH[:-7]
        out_basename = os.path.basename(OUTFILESPATH)        
        OUTPATH = os.path.dirname(OUTFILESPATH)    
        OUT_SNR = os.path.join(OUTPATH, 'SNR_PLOTS')
        print OUTPATH
        # re-initialising .part.pickle collection of cross-correlations
        xc = pscrosscorr.load_pickled_xcorr(PICKLE_PATH)
        
        dataarray = np.asarray(xc[s1][s2].dataarray)
        timearray = np.asarray(xc[s1][s2].timearray)
        
        dataarray = dataarray / np.max(dataarray)
        

        if len(dataarray) > 0 and len(dataarray) == len(timearray):
            x, y = timearray,  dataarray
            s = '{s1}-{s2}: Plot of 100 randomly stacked combination cross-correlations.'
            title = s.format(s1=s1, s2=s2)
                         
            #plt.xlim([start_plot, end_plot])
            plt.title(title)
            plt.ylabel('Amplitude')
            plt.xlabel('Time (UTC)')
        
            plt.plot(x, y, alpha=0.05, c='k')            


        
            file_name = 'Random_Combination_waveform{}-{}-SNR.png'.format(s1, s2)
                                                                     
                                                                     
            print '{s1}-{s2}'.format(s1=s1, s2=s2)
            
            #plt.show()
    fig.savefig(file_name)
    plt.clf()
        
            
        #fig.savefig(outfile_individual, dpi=dpi)
        #fig.clf()
        