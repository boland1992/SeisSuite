#!/usr/bin/python -u
"""
Created on Mon Aug  3 13:10:34 2015

@author: boland
"""

# this is a hack to import modules from a parent directory 
import sys; sys.path.append('../..')
from pysismo import pserrors, pspreprocess
import os
import datetime as dt
import pickle
from pysismo.psconfig import (XCORR_INTERVAL, CROSSCORR_SKIPLOCS, MINFILL,
                              FREQMIN, FREQMAX, FREQMIN_EARTHQUAKE,
                              FREQMAX_EARTHQUAKE, CORNERS, ZEROPHASE,
                              PERIOD_RESAMPLE, ONEBIT_NORM, WINDOW_TIME,
                              WINDOW_FREQ)

# turn on multiprocessing to get one merged trace per station?
# to preprocess trace? to stack cross-correlations?
MULTIPROCESSING = {'merge trace': True,
                   'process trace': True,
                   'cross-corr': True}
                   
# how many concurrent processes? (set None to let multiprocessing module decide)
NB_PROCESSES = None
if any(MULTIPROCESSING.values()):
    import multiprocessing as mp

Preprocess = pspreprocess.Preprocess()

def preprocess_total(date, iterate_stations,
                     dataless_inventories,
                     xml_inventories):
    """    
    The following function is used to pre-process an entire set time-series
    length for use with ambient noise cross-correlation functions. 
    """
    
    print 'b1'
