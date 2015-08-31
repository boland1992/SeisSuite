# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:35:27 2015

@author: boland
"""


from pysismo import pscrosscorr, pserrors, psstation
import os
import sys
import warnings
import datetime as dt
import itertools as it
import pickle
import obspy.signal.cross_correlation
import numpy as np
import time
import glob

from pysismo.psconfig import (CROSSCORR_TMAX)

    

#PICKLE_PATH = '/storage/ANT/PROGRAMS/ANT_OUTPUT/OUTPUT/CROSS/06.05.2015-15:53:28/XCORR-STACK_01.01.2014-31.12.2014_datalesspaz.pickle'
PICKLE_PATH = '/storage/ANT/OUTPUT/CROSS/xcorr_2014-2014_datalesspaz_AU-S-2014-FULL.pickle'

xc = pscrosscorr.load_pickled_xcorr(PICKLE_PATH)

# optimizing time-scale: max time = max distance / vmin (vmin = 2.5 km/s)
maxdist = max([xc[s1][s2].dist() for s1, s2 in xc.pairs()])
maxt = min(CROSSCORR_TMAX, maxdist)
    

#plot distance plot of cross-correlations
xc.plot(plot_type='distance', xlim=(-maxt, maxt), 
            outfile="/home/boland/Desktop/something1342.png", showplot=False)
            
            