# -*- coding: utf-8 -*-
"""
Created on Fri May  8 11:16:30 2015

@author: boland
"""
from obspy import read_inventory
from obspy import UTCDateTime
import numpy as np
import os
from obspy import read
from scipy import signal
import matplotlib.pyplot as plt
from numpy.lib.stride_tricks import as_strided
from numpy.fft import rfft, irfft
import glob
import datetime as dt
import scipy

t0 = dt.datetime(2014, 10, 01, 00, 00)
t1 = t0 + dt.timedelta(days = 31)

XCORR_INTERVAL = 45.0

N = int(((t1 - t0).days)*60*24 / XCORR_INTERVAL)
print(N)


times = [t0 + dt.timedelta(minutes=i) for i in \
         [j*XCORR_INTERVAL for j in range(N)]] 
         
times = [UTCDateTime(i) for i in times]
         
dir1 = '/home/boland/Desktop/Link to ANT_OUTPUT/\
INPUT/DATA/2014-10/UM.HOLS.EHZ.mseed'
dir2 = '/home/boland/Desktop/Link to ANT_OUTPUT/\
INPUT/DATA/2014-10/UM.LOCU.EHZ.mseed'
         
def get_stationxml_inventories(stationxml_dir, verbose=False):
    """
    Reads inventories in all StationXML (*.xml) files
    of specified dir

    @type stationxml_dir: unicode or str
    @type verbose: bool
    @rtype: list of L{obspy.station.inventory.Inventory}
    """
    inventories = []

    # list of *.xml files
    flist = glob.glob(pathname=os.path.join(stationxml_dir, "*.xml"))

    if verbose:
        if flist:
            print "Reading inventory in StationXML file:",
        else:
            s = u"Could not find any StationXML file (*.xml) in dir: {}!"
            print s.format(stationxml_dir)

    for f in flist:
        if verbose:
            print os.path.basename(f),
        inv = read_inventory(f, format='stationxml')
        inventories.append(inv)

    if flist and verbose:
        print

    return inventories

def spectrum(tr):

    wave = tr.data #this is how to extract a data array from a mseed file
    fs = tr.stats.sampling_rate
                
                #hour = str(hour).zfill(2) #create correct format for eqstring
    f, Pxx_spec = signal.welch(wave, fs, 'flattop', 1024, scaling='spectrum')
                    #plt.semilogy(f, np.sqrt(Pxx_spec))
    plt.title("Frequency Density Plot of PNG Earthquake from station PMG.IU")
    plt.plot(f, np.sqrt(Pxx_spec))
    plt.xlim([0, 5])
    #plt.ylim([0, 0.01])
            
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Linear spectrum [V RMS]')
    plt.show()
    
    
def resample(trace, dt_resample):
    """
    Subroutine to resample trace

    @type trace: L{obspy.core.trace.Trace}
    @type dt_resample: float
    @rtype: L{obspy.core.trace.Trace}
    """
    dt = 1.0 / trace.stats.sampling_rate
    factor = dt_resample / dt
    if int(factor) == factor:
        # simple decimation (no filt because it shifts the data)
        trace.decimate(int(factor), no_filter=True)
    else:
        # linear interpolation
        tp = np.arange(0, trace.stats.npts) * trace.stats.delta
        zp = trace.data
        ninterp = int(max(tp) / dt_resample) + 1
        tinterp = np.arange(0, ninterp) * dt_resample

        trace.data = np.interp(tinterp, tp, zp)
        trace.stats.npts = ninterp
        trace.stats.delta = dt_resample
        trace.stats.sampling_rate = 1.0 / dt_resample
        #trace.stats.endtime = trace.stats.endtime + max(tinterp)-max(tp)

    
def moving_avg(a, halfwindow, mask=None):
    """
    Performs a fast n-point moving average of (the last
    dimension of) array *a*, by using stride tricks to roll
    a window on *a*.

    Note that *halfwindow* gives the nb of points on each side,
    so that n = 2*halfwindow + 1.

    If *mask* is provided, values of *a* where mask = False are
    skipped.

    Returns an array of same size as *a* (which means that near
    the edges, the averaging window is actually < *npt*).
    """
    # padding array with zeros on the left and on the right:
    # e.g., if halfwindow = 2:
    # a_padded    = [0 0 a0 a1 ... aN 0 0]
    # mask_padded = [F F ?  ?      ?  F F]

    if mask is None:
        mask = np.ones_like(a, dtype='bool')

    zeros = np.zeros(a.shape[:-1] + (halfwindow,))
    falses = zeros.astype('bool')

    a_padded = np.concatenate((zeros, np.where(mask, a, 0), zeros), axis=-1)
    mask_padded = np.concatenate((falses, mask, falses), axis=-1)

    # rolling window on padded array using stride trick
    #
    # E.g., if halfwindow=2:
    # rolling_a[:, 0] = [0   0 a0 a1 ...    aN]
    # rolling_a[:, 1] = [0  a0 a1 a2 ... aN 0 ]
    # ...
    # rolling_a[:, 4] = [a2 a3 ...    aN  0  0]

    npt = 2 * halfwindow + 1  # total size of the averaging window
    rolling_a = as_strided(a_padded,
                           shape=a.shape + (npt,),
                           strides=a_padded.strides + (a.strides[-1],))
    rolling_mask = as_strided(mask_padded,
                              shape=mask.shape + (npt,),
                              strides=mask_padded.strides + (mask.strides[-1],))

    # moving average
    n = rolling_mask.sum(axis=-1)
    return np.where(n > 0, rolling_a.sum(axis=-1).astype('float') / n, np.nan)


def butterworth(trace):
    #filter
    #print("first filter")

    trace.filter(type="bandpass",
             freqmin =freqmin,
             freqmax = freqmax,
             corners=corners,
             zerophase=zerophase)
    


    return trace


def normal(trace,
           freqmin_earthquake, 
           freqmax_earthquake):

    # normalization of the signal by the running mean
    # in the earthquake frequency band
    trcopy = trace
    #print("normalising filter")

    trcopy.filter(type="bandpass",
              freqmin=freqmin_earthquake,
              freqmax=freqmax_earthquake,
              corners=corners,
              zerophase=zerophase)

    # decimating trace
    resample(trcopy, period_resample)

    # Time-normalization weights from smoothed abs(data)
    # Note that trace's data can be a masked array
    halfwindow = int(round(window_time * trcopy.stats.sampling_rate / 2))
    mask = ~trcopy.data.mask if np.ma.isMA(trcopy.data) else None
    tnorm_w = moving_avg(np.abs(trcopy.data),
                             halfwindow=halfwindow,
                             mask=mask)
    if np.ma.isMA(trcopy.data):
        # turning time-normalization weights into a masked array
        s = "[warning: {}.{} trace's data is a masked array]"
        print s.format(trace.stats.network, trace.stats.station),
        tnorm_w = np.ma.masked_array(tnorm_w, trcopy.data.mask)

    # time-normalization
    trace.data /= tnorm_w
    
    return trace


def whiten(trace, window_freq, freqmin, freqmax, corners, zerophase):
    """
    function that produces a whitened spectrum
    """
    
    fft = rfft(trace.data)  # real fast fourier transform
    deltaf = trace.stats.sampling_rate / trace.stats.npts  # frequency step
    # smoothing amplitude spectrum
    halfwindow = int(round(window_freq / deltaf / 2.0))
    
    weight = moving_avg(abs(fft), halfwindow=halfwindow)
    

    # normalizing spectrum and back to time domain
    trace.data = irfft(fft / weight, n=len(trace.data))
    # re bandpass to avoid low/high freq noise
    #print("Whiten filter")
    trace.filter(type="bandpass",
             freqmin =freqmin,
             freqmax = freqmax,
             corners=corners,
             zerophase=zerophase)
    
    return trace

def preprocess(trace):
    
    #trace.attach_response(inventories=xml_inventories)

    trace = butterworth(trace)

    #trace.remove_response(output="VEL", zero_mean=True)

    #plt.figure()
    #spectrum(trace)
    #trace = normal(trace, freqmin_earthquake, freqmax_earthquake)
    #plt.figure()
    #spectrum(trace)
    #print(trace.stats.sampling_rate)
    trace = whiten(trace, window_freq, freqmin, freqmax, corners, zerophase)
    #plt.figure()
    #spectrum(trace)
    return trace
        
xcorr = 0
freqmin = 1.0/25.0
freqmax = 1.0/1
corners = 1
zerophase = True
freqmin_earthquake = 1/50.0
freqmax_earthquake = 1/25.0
window_time = 0.5 * freqmax_earthquake
window_freq = 0.02
period_resample = 0.45
STATIONXML_DIR = '/storage/ANT/PROGRAMS/ANT_OUTPUT/INPUT/XML'
xml_inventories = []

sample_rate = 250
counts = 0
for time in times:

    st0 = read(dir1, starttime=time,
               endtime=time + dt.timedelta(minutes=XCORR_INTERVAL))
    st1 = read(dir2, starttime=time, 
               endtime=time + dt.timedelta(minutes=XCORR_INTERVAL))
    
    tr0 = st0[0]
    tr1 = st1[0]
    
    tr0.stats.sampling_rate = sample_rate     
    tr1.stats.sampling_rate = sample_rate     

    tr0 = preprocess(tr0)
    tr1 = preprocess(tr1)

    xcorr = scipy.signal.correlate(tr0, tr1, mode='same')

    xcorr += xcorr
    
    plt.figure(1)
    plt.plot(xcorr)
    plt.show()
    print(counts)
    counts +=1