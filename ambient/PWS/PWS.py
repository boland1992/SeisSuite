# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 12:17:56 2015

@author: boland

The following script is being used in order to explore and develop 
python methods for phase-stacking, and phase-weighted stacking between
two seismic waveforms. Input uses one file per station waveform. Needs
a minimum of two channels to stack to work!
"""

from math import atan2
from obspy import read
from obspy.core import Stream
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.trigger import *

path = 'NARR.mseed'
path = 'events/20150130_0016.uom.HOLS.mseed'
path = 'events/20150130_0016.uom.MRDN.mseed'
path = 'events/20150130_0016.uom.NARR.mseed'
path = 'events/20150130_0016.uom.OUTU.mseed'

st = read(path)
#st = read(path, headonly=True)
tr = st[0]
start = tr.stats.starttime
end = tr.stats.starttime + tr.stats.npts*(1/tr.stats.sampling_rate)
delta_t = (end-start)

#running normalise time-windows e.g. 5,4,3 and 2 seconds
time_windows = [4,3,2]
#test time windows with delta_t. We want the largest window of the options. 
moduli = []
time_window = False
for i, time in enumerate(time_windows):
    modulus = delta_t % time
    if modulus == 0:
        time_window = time
        break
    else:
        moduli.append([modulus, time])

if not time_window:    
    moduli = np.sort(moduli, axis=0)
    time_window = moduli[0][1]

#new_start = tr.stats.starttime + 60
#st = read(path, starttime=new_start, endtime=end)
st = st.filter('bandpass', freqmin=1., 
               freqmax=10.0, corners=2, 
               zerophase=True)

#==============================================================================
#LINEAR STACKING
#==============================================================================
LS = 0
prev_shape = st[0].data.shape

for tr in st:
    #normalise the traces to between 0 and 1
    curr_shape = tr.data.shape
    if curr_shape == prev_shape:
        LS = (LS - np.mean(LS))
        #stack the traces linearly
        LS += tr.data / np.max(tr.data)

    prev_shape = tr.data.shape


N = len(st)

#c = np.abs(c) / N
LS = LS / N

print "Linear Stack SNR: ", np.abs(np.mean(LS)) / np.abs(np.std(LS))

plt.figure()

plt.plot(LS, c='g', alpha=0.5, zorder=2)

#==============================================================================
#PHASE STACKING
#==============================================================================


PS = 0
    
prev_shape = st[0].data.shape

for tr in st:
    
    curr_shape = tr.data.shape
    if curr_shape == prev_shape:
        tr_norm = tr.data / np.max(tr.data)
        inst_phase = np.arctan2(tr.data, range(0,len(tr.data)))
        #PS = (PS - np.mean(PS))
        PS += np.real(np.exp(1j*inst_phase))
        prev_shape = tr.data.shape

    
N = len(st)
#PS = PS/N
PS = np.abs(PS) / N


# normalise about zero!
#PS = (PS - np.mean(PS))


print "Phase Stack SNR: ", np.abs(np.mean(PS)) / np.abs(np.std(PS))


plt.plot(PS, zorder=3, alpha=0.5)

plt.plot(PS, c='y', zorder=1)

#==============================================================================
#PHASE-WEIGHTED STACKING
#==============================================================================
sharp_v = 2
PWS = (LS * PS ** sharp_v)
PWS = PWS - np.mean(PWS)

PWS = np.max(LS)/np.max(PWS) * PWS
print "Phase-Weighted Stack SNR: ", np.abs(np.mean(PWS)) / np.abs(np.std(PWS))


plt.plot(PWS, c='r', zorder=1)


path = '20150731_0802.agospublic.all.mseed'
st = read(path)

#==============================================================================
#SMOOTHED PHASE STACK
#==============================================================================

SPS = PS
    
prev_shape = st[0].data.shape

for tr in st:
    no_samples = tr.stats.npts
    sample_rate = tr.stats.sampling_rate
    sample_window = tr.stats.sampling_rate * time_window

    time_constant = sample_rate / (2 * sample_window + sample_rate)  
    
    curr_shape = tr.data.shape
    if curr_shape == prev_shape:
        
        for t in np.arange(sample_window/2, no_samples, 
                           int(sample_window/2))[:-1]:
            
            SPS[t-sample_window/2:t-sample_window/2] = time_constant * \
            PS[t-sample_window/2:t-sample_window/2]        
                 

        #PS = (PS - np.mean(PS))
        prev_shape = tr.data.shape

#==============================================================================
#SMOOTHED PHASE-WEIGHTED STACK
#==============================================================================
sharp_v = 2
SPWS = (LS * SPS ** sharp_v)
SPWS = SPWS - np.mean(SPWS)

SPWS = np.max(LS)/np.max(SPWS) * SPWS
print "Phase-Weighted Stack SNR: ", np.abs(np.mean(SPWS)) / np.abs(np.std(SPWS))

plt.plot(SPWS, c='orange', zorder=5)
plt.xlim([0,len(tr.data)])
plt.show()


quit()
#==============================================================================
# split into one file per station!
#==============================================================================
stations = []
for tr in st:
    stat = tr.stats.station
    if stat not in stations:
        stations.append(stat)


station_streams = {}
for stat in stations:
    station_streams[stat] = []    
    for tr in st:
        if tr.stats.station == stat:
            station_streams[stat].append(tr)
            

for key in station_streams.keys():
    
    station_streams[key]
        
    stat = station_streams[key][0].stats.station        
    STREAM = Stream(station_streams[key])
    STREAM.write('{}.mseed'.format(stat), format='MSEED')

#==============================================================================
#==============================================================================
