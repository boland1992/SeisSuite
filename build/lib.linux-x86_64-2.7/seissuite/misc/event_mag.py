# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 08:44:28 2015

@author: boland
"""

import numpy as np
import obspy
from obspy import UTCDateTime, read
import datetime
from obspy.fdsn import Client
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# set server name. United States Geological Survey. 
client = Client("IRIS") 

# get current UTC date in "yy-mm-dd" format and set that as the end of the searched time-window
date_end = datetime.datetime.utcnow().date(); 

#convert end date from datetime object to a UTCDateTime object
end = UTCDateTime(date_end)

# set the time period to scan. in this case we're looking at the previous 10 days
no_of_days = 1000.0

# define time difference as a datetime object 
number_days = datetime.timedelta(days=no_of_days)

# set start date for the time-window as the current date minus the number of days set 
date_start = date_end - number_days

#convert start date from datetime object to a UTCDateTime object
start = UTCDateTime(date_start)

# set minimum magnitude threshold to search for.
min_mag = 2.0

minlat = -39.1 #minlatitude
maxlat = -37.0 # maxlatitude 
minlon = 143.5 #minlongitude 
maxlon = 147.0 #maxlongitude 

cat = client.get_events(starttime=date_start, 
                        endtime=date_end, 
                        minlatitude=minlat,
                        maxlatitude=maxlat,
                        minlongitude=minlon,
                        maxlongitude=maxlon,
                        minmagnitude=min_mag)

#print(cat)

cat.plot()
print(cat.__str__(print_all=True))

net = 'AU' 

stat = 'TOO'

date_start = UTCDateTime("2003-10-18T10:29:26.580000Z")

date_end = date_start + 3600

st = client.get_waveforms(net, stat, "00", "*Z", 
                          date_start, date_end,
                          attach_response=True)

st.plot()
#st.write('Gippsland_low.MSEED', format='MSEED') 
