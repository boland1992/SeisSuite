# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:35:27 2015

@author: boland
"""
   
from seissuite.ant import (pscrosscorr)

import glob
import os
import pickle
import numpy as np


#PICKLE_PATH = '/storage/ANT/PROGRAMS/ANT_OUTPUT/OUTPUT/CROSS/06.05.2015-15:53:28/XCORR-STACK_01.01.2014-31.12.2014_datalesspaz.pickle'
#PICKLE_PATH = '/home/boland/Desktop/XCORR-STACK_01.08.1999-10.06.2000_datalesspaz.part.pickle'
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
CROSSCORR_TMAX = CONFIG.CROSSCORR_TMAX

pickle_list = []
folder_list = sorted(glob.glob(os.path.join(CROSSCORR_DIR, '*')))



plot_distance = True
plot_classic = False

print MSEED_DIR
print CROSSCORR_TMAX

central_frequencies = np.arange(0.1, 20.0, 0.5)
central_frequencies = central_frequencies[::-1]


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

    
if not res:
    raise Exception("You must choose one a number betwen {} and {}"\
    .format(1, len(pickle_list)))
    
else:
    PICKLE_PATH = pickle_list[int(res)-1]
    OUTFILESPATH = PICKLE_PATH[:-7]
    out_basename = os.path.basename(OUTFILESPATH)        

    print "\nOpening {} file to process ... ".format(out_basename)
    # re-initialising .part.pickle collection of cross-correlations
    xc = pscrosscorr.load_pickled_xcorr(PICKLE_PATH)
                            
    # optimizing time-scale: max time = max distance / vmin (vmin = 2.5 km/s)
    maxdist = max([xc[s1][s2].dist() for s1, s2 in xc.pairs()])
    maxt = max(CROSSCORR_TMAX, maxdist / 2.5)
    
    if plot_distance:
        for central_freq in central_frequencies:
            #plot distance plot of cross-correlations
            xc.plot(plot_type='distance', xlim=(-maxt, maxt), 
                    outfile="{}~.png".format(central_freq), showplot=False, 
                    norm=True, fill=False, absolute=False, 
                    freq_central=central_freq)



    if plot_classic:
        #plot individual cross-correlations
        xc.plot(plot_type='classic', xlim=(-maxt, maxt), 
                outfile="~", showplot=False)