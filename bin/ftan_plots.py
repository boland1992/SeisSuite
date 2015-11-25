# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:31:46 2015

@author: iese
"""

from seissuite.ant import (pscrosscorr)

import os
import glob
import pickle
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

# import CONFIG class initalised in ./configs/tmp_config.pickle
config_pickle = 'configs/tmp_config.pickle'
f = open(name=config_pickle, mode='rb')
CONFIG = pickle.load(f)
f.close()
    
# import variables from initialised CONFIG class.
FTAN_DIR = CONFIG.FTAN_DIR
TOMO_DIR = CONFIG.TOMO_DIR

RAWFTAN_PERIODS = CONFIG.RAWFTAN_PERIODS

PERIODS = RAWFTAN_PERIODS #np.arange(0.07, 9.0, 0.1)


CLEANFTAN_PERIODS = CONFIG.RAWFTAN_PERIODS
FTAN_VELOCITIES = CONFIG.FTAN_VELOCITIES



# make a check to see that all PERIODS are in RAWFTAN_PERIODS

# selecting dispersion curves
pickle_files = sorted(glob.glob(os.path.join(FTAN_DIR, 'FTAN.pickle')))

# loop on pickled curves
for pickle_file in pickle_files:
    print "\nImporting dispersion curves from file: " + pickle_file

    f = open(pickle_file, 'rb')
    curves = pickle.load(f)
    f.close()
    
    for cleanvg in curves:
        stat1, stat2 = cleanvg.station1.name, cleanvg.station2.name 

        plt.figure(1)
        plt.title("Dispersion curve for stacked cross-correlation between \
#{} and {}".format(stat1, stat2))
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Velocity (km/s)")
        plt.plot(1./cleanvg.periods, cleanvg.v)        
        plt.show()
        plt.clf()
