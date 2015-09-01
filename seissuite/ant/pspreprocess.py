#!/usr/bin/env python
"""
Module that contains classes pertaining to key pre-prosessing methods
with respect to ambient noise seismic waveforms. These practically
should then be used to find the cross-correlation function of the waveforms
with time.
"""

import pserrors, psstation, psutils
import obspy.signal
import obspy.xseed
import obspy.signal.cross_correlation
import obspy.signal.filter
from obspy import read
import numpy as np
from numpy.fft import rfft, irfft
from obspy import UTCDateTime
import datetime as dt


# ====================================================
# parsing configuration file to import some parameters
# ====================================================
from psconfig import (CROSSCORR_SKIPLOCS, MINFILL, MSEED_DIR)

# ========================
# Constants and parameters
# ========================

EPS = 1.0e-5
ONESEC = dt.timedelta(seconds=1)

class Preprocess:
    """
    Class for performing all possible preprocess steps on a seismic waveform
    to then allow for smooth cross-correlation.
    
    @param freqmin: low frequency of the band-pass filter
    @param freqmax: high frequency of the band-pass filter
    @param freqmin_earthquake: low frequency of the earthquake band
    @param freqmax_earthquake: high frequency of the earthquake band
    @param corners: nb or corners of the band-pass filter
    @param zerophase: set to True for filter not to shift phase
    @type zerophase: bool
    @param period_resample: resampling period in seconds
    @param onebit_norm: set to True to apply one-bit normalization (else,
                        running mean normalization is applied)
    @type onebit_norm: bool
    @param window_time: width of the window to calculate the running mean
                        in the earthquake band (for the time-normalization)
    @param window_freq: width of the window to calculate the running mean
                        of the amplitude spectrum (for the 
                        spectral whitening)
    """
    
    def __init__(self, freqmin, freqmax, freqmin_earthquake,
                       freqmax_earthquake, corners, zerophase,
                       period_resample, window_time, window_freq,
                       onebit_norm, *args, **kwargs):
        
        # initialising config file inter-module variables 
        self.freqmin = freqmin
        self.freqmax = freqmax
        self.onebit_norm = onebit_norm             
        self.freqmin_earthquake = freqmin_earthquake
        self.freqmax_earthquake = freqmax_earthquake
        self.corners = corners
        self.zerophase = zerophase
        self.period_resample = period_resample
        self.window_time = window_time
        self.window_freq = window_freq
        
    def remove_resp(self, trace, paz=None):
        
        # removing response...
        if paz:
            # ...using paz:
            if trace.stats.sampling_rate > 10.0:
                # decimating large trace, else fft crashes
                factor = int(np.ceil(trace.stats.sampling_rate / 10))
                trace.decimate(factor=factor, no_filter=True)
            trace.simulate(paz_remove=paz,
                           paz_simulate=obspy.signal.cornFreq2Paz(0.01),
                           remove_sensitivity=True,
                           simulate_sensitivity=True,
                           nfft_pow2=True)
        else:
            # ...using StationXML:
            # first band-pass to downsample data before removing response
            # (else remove_response() method is slow or even hangs)
            trace.filter(type="bandpass",
                         freqmin=self.freqmin,
                         freqmax=self.freqmax,
                         corners=self.corners,
                         zerophase=self.zerophase)
            psutils.resample(trace, dt_resample=self.period_resample)
            trace.remove_response(output="VEL", zero_mean=True)
        return trace
        
    def bandpass_filt(self, trace):
        """
        Function to apply a butterworth bandpass-filter to an obspy
        trace input object. Note: MUST only be one trace. No streams.
        """
    
        # band-pass
        return trace.filter(type="bandpass",
                            freqmin=self.freqmin,
                            freqmax=self.freqmax,
                            corners=self.corners,
                            zerophase=self.zerophase)
                     
    def trace_downsample(self, trace):

        # downsampling trace if not already done
        if abs(1.0 / trace.stats.sampling_rate - self.period_resample) > EPS:
            psutils.resample(trace, dt_resample=self.period_resample)
            
        return trace

    def time_norm(self, trace, trcopy):
        # normalization of the signal by the running mean
        # in the earthquake frequency band
        trcopy.filter(type="bandpass",
                      freqmin=self.freqmin_earthquake,
                      freqmax=self.freqmax_earthquake,
                      corners=self.corners,
                      zerophase=self.zerophase)
        # decimating trace
        psutils.resample(trcopy, self.period_resample)

        # Time-normalization weights from smoothed abs(data)
        # Note that trace's data can be a masked array
        halfwindow = int(round(self.window_time*trcopy.stats.sampling_rate/2))
        mask = ~trcopy.data.mask if np.ma.isMA(trcopy.data) else None
        tnorm_w = psutils.moving_avg(np.abs(trcopy.data),
                                     halfwindow=halfwindow,
                                     mask=mask)
        if np.ma.isMA(trcopy.data):
            # turning time-normalization weights into a masked array
            s = "[warning: {}.{} trace's data is a masked array]"
            print s.format(trace.stats.network, trace.stats.station),
            tnorm_w = np.ma.masked_array(tnorm_w, trcopy.data.mask)

        if np.any((tnorm_w == 0.0) | np.isnan(tnorm_w)):
            # illegal normalizing value -> skipping trace
            raise pserrors.CannotPreprocess("Zero or NaN normalisation \
                                            weight")

        # time-normalization
        trace.data /= tnorm_w
        
        return trace
        
    def spectral_whitening(self, trace):
        """
        Function that takes an input obspy trace object that has been
        time-normalised, band-pass filtered and had its repsonse removed,
        delimited, demeaned and detrended. 
        """
        # real FFT
        fft = rfft(trace.data) 
        # frequency step
        deltaf = trace.stats.sampling_rate / trace.stats.npts  
        # smoothing amplitude spectrum
        halfwindow = int(round(self.window_freq / deltaf / 2.0))
        weight = psutils.moving_avg(abs(fft), halfwindow=halfwindow)
        # normalizing spectrum and back to time domain
        trace.data = irfft(fft / weight, n=len(trace.data))
        # re bandpass to avoid low/high freq noise
        trace.filter(type="bandpass",
                     freqmin=self.freqmin,
                     freqmax=self.freqmax,
                     corners=self.corners,
                     zerophase=self.zerophase)
        return trace
    
    def get_or_attach_response(self, trace, 
                                     dataless_inventories=(), 
                                     xml_inventories=()):
        """
        Returns or attach instrumental response, from dataless seed inventories
        (as returned by psstation.get_dataless_inventories()) and/or StationXML
        inventories (as returned by psstation.get_stationxml_inventories()).
        If a response if found in a dataless inventory, then a dict of poles
        and zeros is returned. If a response is found in a StationXML
        inventory, then it is directly attached to the trace and nothing is
        returned.

        Raises CannotPreprocess exception if no instrumental response is found.

        @type trace: L{Trace}
        @param dataless_inventories: inventories from dataless seed files 
                                     (as returned by 
                                      psstation.get_dataless_inventories())
        @type dataless_inventories: list of L{obspy.xseed.parser.Parser}
        @param xml_inventories: inventories from StationXML files 
                                (as returned by 
                                 psstation.get_stationxml_inventories())
        @type xml_inventories: list of L{obspy.station.inventory.Inventory}
        """
        t1 = dt.datetime.now()
        # looking for instrument response...
        try:
            # ...first in dataless seed inventories
            paz = psstation.get_paz(channelid=trace.id,
                                    t=trace.stats.starttime,
                                    inventories=dataless_inventories)
            return paz
        except pserrors.NoPAZFound:
            # ... then in dataless seed inventories, replacing 'BHZ' with 'HHZ'
            # in trace's id (trick to make code work with Diogo's data)
            try:
                paz = psstation.get_paz(channelid=trace.id.replace('BHZ','HHZ'),
                                        t=trace.stats.starttime,
                                        inventories=dataless_inventories)
                return paz
            except pserrors.NoPAZFound:
                # ...finally in StationXML inventories
                try:
                    trace.attach_response(inventories=xml_inventories)
                except:
                    # no response found!
                    raise pserrors.CannotPreprocess("No response found")
    
        delta = (dt.datetime.now() - t1).total_seconds()
        print "\nProcessed response attachment in {:.1f} seconds".format(delta)


    def preprocess_trace(self, trace, paz=None):
        """
        Preprocesses a trace (so that it is ready to be cross-correlated),
        by applying the following steps:
            - removal of instrument response, mean and trend
            - band-pass filtering between *freqmin*-*freqmax*
            - downsampling to *period_resample* secs
            - time-normalization (one-bit normalization or normalization
              by the running mean in the earthquake frequency band)
            - spectral whitening (if running mean normalization)

        Raises CannotPreprocess exception if:
            - trace only contains 0 (happens sometimes...)
            - a normalization weight is 0 or NaN
            - a Nan appeared in trace data
            
        Note that the processing steps are performed in-place.

        @type trace: L{Trace}
        @param paz: poles and zeros of instrumental response
                    (set None if response is directly attached to trace)

        """

        # ============================================
        # Removing instrument response, mean and trend
        # ============================================
        
        if paz is not None:
            trace = self.remove_resp(trace, paz=paz)

        # trimming, demeaning, detrending
        midt = trace.stats.starttime + (trace.stats.endtime -
                                        trace.stats.starttime) / 2.0
        t0 = UTCDateTime(midt.date)  # date of trace, at time 00h00m00s
        trace.trim(starttime=t0, endtime=t0 + dt.timedelta(days=1))
        trace.detrend(type='constant')
        trace.detrend(type='linear')

        if np.all(trace.data == 0.0):   
            # no data -> skipping trace
            raise pserrors.CannotPreprocess("Only zeros")
            
        # take a copy of the trace to calculate weights of time-normalization
        trcopy = trace
    
        # =========
        # Band-pass
        # =========
        trace = self.bandpass_filt(trace)
        # =========
        # Downsample
        # =========
        self.trace_downsample(trace)
        
        # ==================
        # Time normalization
        # ==================
        if self.onebit_norm:
            # one-bit normalization
            trace.data = np.sign(trace.data)
        else:
            trace.data = self.time_norm(trace, trcopy)
            # ==================
            # Spectral whitening
            # ==================
        trace = self.spectral_whitening(trace)
                    
        # Verifying that we don't have nan in trace data
        if np.any(np.isnan(trace.data)):
            raise pserrors.CannotPreprocess("Got NaN in trace data")
            
        return trace

    def get_merged_trace(self, station, date, xcorr_interval, 
                         skiplocs=CROSSCORR_SKIPLOCS, 
                         minfill=MINFILL):
        """
        Returns one trace extracted from selected station, at selected date
        (+/- 1 hour on each side to avoid edge effects during subsequent
        processing).
    
        for 45 minute xcorr interval, change edge effects by 1 minute!
    
    
    
        Traces whose location belongs to *skiplocs* are discarded, then
        if several locations remain, only the first is kept. Finally,
        if several traces (with the same location) remain, they are
        merged, WITH GAPS FILLED USING LINEAR INTERPOLATION.

        Raises CannotPreprocess exception if:
            - no trace remain after discarded the unwanted locations
            - data fill is < *minfill*

        @type station: L{psstation.Station}
        @type date: L{datetime.date}
        @param skiplocs: list of locations to discard in station's data
        @type skiplocs: iterable
        @param minfill: minimum data fill to keep trace
        @rtype: L{Trace}
        """
        #calculate edge addition and subtraction as 1/24th of the overall time interval
        #
        startminutes = (xcorr_interval / 24.0)
        endminutes = xcorr_interval + startminutes
    
        # getting station's stream at selected date
        # (+/- one hour to avoid edge effects when removing response)
        t0 = UTCDateTime(date)  # date at time 00h00m00s
        
        path_start = t0 - dt.timedelta(minutes=startminutes)
        path_end = t0 + dt.timedelta(minutes=endminutes)
        
        
        #station_path_old = station.getpath(date, MSEED_DIR)
        station_path_SQL = station.getpath(t0, t0+dt.timedelta(minutes=xcorr_interval))
        #print "station old path: ", station_path_old
        #print "station SQl path: ", station_path_SQL
        
        st = read(pathname_or_url=station_path_SQL,
                  starttime=path_start, endtime=path_end)
                  
        # removing traces of stream from locations to skip
        for tr in [tr for tr in st if tr.stats.location in skiplocs]:
            st.remove(tr)

        if not st.traces:
            # no remaining trace!
            raise pserrors.CannotPreprocess("No trace")

        # if more than one location, we retain only the first one
        if len(set(tr.id for tr in st)) > 1:
            select_loc = sorted(set(tr.stats.location for tr in st))[0]
            for tr in [tr for tr in st if tr.stats.location != select_loc]:
                st.remove(tr)

        # Data fill for current date
        fill = psutils.get_fill(st, starttime=t0, 
                                endtime=t0 + dt.timedelta(minutes=endminutes))
        if fill < minfill:
            # not enough data
            raise pserrors.CannotPreprocess("{:.0f}% fill".format(fill * 100))

        # Merging traces, FILLING GAPS WITH LINEAR INTERP
        st.merge(fill_value='interpolate')
        trace = st[0]

        return trace
        
if __name__ == '__main__':
    # loading pickled cross-correlations
    a=5