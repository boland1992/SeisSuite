# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 13:36:39 2015

@author: boland
"""

import pickle
import datetime
from obspy import read
import numpy as np
from obspy.core import UTCDateTime

verbose = False

t0 = datetime.datetime.now()

timeline_pickle = 'tmp/timeline_database.pickle'
f = open(name=timeline_pickle, mode='rb')
timeline = pickle.load(f)
f.close()

t1 = datetime.datetime.now()

if verbose:
    print "Time taken to import timeline database pickle file was: ", t1-t0


t0 = datetime.datetime.now()
# find absolute start and end times for the database and set them as var
time_list = []
for stat in timeline.keys():
    for time_info in timeline[stat]:
        time_list.append(time_info[0].timestamp)
        time_list.append(time_info[1].timestamp)

timeline_start = np.min(time_list)
timeline_end = np.max(time_list)
        
t1 = datetime.datetime.now()

if verbose:
    print "Time taken to find start and end times from databse was: ", t1-t0
