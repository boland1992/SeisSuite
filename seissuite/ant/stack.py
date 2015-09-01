# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 12:17:56 2015

@author: boland

The following script is being used in order to explore and develop 
python methods for phase-stacking, and phase-weighted stacking between
two seismic waveforms. Input uses one file per station waveform. Needs
a minimum of two channels to stack to work!
"""

from obspy import read
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import datetime as dt

class Stack:
    """
    The following class contains necessary functions to perform seismic
    waveform stacking. The stacking improves Signal-to-Noise ratio (SNR)
    and allows for better picks when doing event location. 
    """
    
    def __init__(self, st_path, filter_waveform=False, band_lims=[1,10]):
        # set path to multiplexed event waveform
        self.st_path = st_path
        # initialise input stream variable
        self.st = None
        # initialise linear stack variable | type <list>
        self.LS = None
        # initialise phase stack variable | type <list>
        self.PS = None
        # initialise phase-weighted stack variable | type <list>
        self.PWS = None
        self.filter = filter_waveform
        if self.filter:
            # initialise band-pass filter frequency limits
            self.band_lims = band_lims
        self.times = None
        self.starttime = None
        self.endtime = None
        
    def import_stream(self, st_path=None):
        """
        Function that uses the read function from obspy.core to import t
        the given Stream object from miniseed. 
        """
        if st_path is None:
            st_path = self.st_path
        self.st = read(st_path)
        
        # import filtered waveform if specified!
        if self.filter:
            freq_min, freq_max = self.band_lims
            self.st = self.st.filter('bandpass', freqmin=freq_min, 
                                     freqmax=freq_max, corners=2, 
                                     zerophase=True)
        
        # run the datettime_list function on import!
        self.datetime_list()
        
        return self.st
        
    def datetime_list(self, st=None):
        """
        Function that outputs a list of N datettime objects between the
        start and end times of a given stream. N is dictated by the length
        of the st[0].data list. 
        """
        if st is None:
            st = self.st
        # if st is still None
        if st is None:
            self.import_stream()
            st = self.st
        
        starttime = st[0].stats.starttime
        endtime = starttime + st[0].stats.npts * (1/st[0].stats.sampling_rate)
        time_delta = (endtime - starttime)
        self.starttime = str(starttime)
        self.endtime = str(endtime)
        # convert to datetime.datetime object        
        starttime = dt.datetime.strptime(str(starttime), 
                                         '%Y-%m-%dT%H:%M:%S.%fZ')
        
        nlist = len(st[0].data)
        time_space = np.linspace(0, time_delta, nlist)
        self.times = [starttime + dt.timedelta(milliseconds=1000*time) 
                      for time in time_space]
        

    def lin_stack(self, st=None, norm=True):
        """
        The following function takes an obspy stream object and outputs 
        the linear stack of all traces within the stream. This is only useful
        if all traces are various channels from the same station, OR are from 
        the same cross-correlation station pairs for ambient studies. 
        """
        
        if st is None:
            st = self.st
        # if st is STILL None, then run the import stream!
        if st is None:
            self.import_stream()
            st = self.st
            
        LS = 0
        
        for tr in st:
            #stack the traces linearly
            LS += tr.data - np.mean(tr.data)            
        N = len(st)
        if norm:
            #normalise the traces to between 0 and 1
            self.LS = LS / N / np.max(LS)
        else: 
            self.LS = LS / N 
            
        #print "Linear Stack SNR: ", np.abs(np.mean(LS)) / np.abs(np.std(LS))


    def plot_lin_stack(self, LS=None, show=True, save=False):
        """
        Function to plot the linear stack of two or more traces within 
        the input Stream object for the class: Stack. Show is the default, 
        if save is True, then show is automatically set to False
        """
        
        # could add feature to allow for addition of own output path!        
        # could add feature to allow for higher resolution output images!
        
        if LS is None:
            LS = self.LS
                
        # if LS is STILL None, run the lin_stack function
        if LS is None:
            try:
                self.lin_stack()
                LS = self.LS
            except Exception as error:
                print error
                
        fig = plt.figure(figsize=(15,10))
        plt.xticks(rotation=10)
        plt.plot(self.times, LS, color='k')
        plt.title('Linearly Stacked Waveform From File: {} \n {} - {}'\
        .format(os.path.basename(self.st_path), self.starttime, self.endtime))
        plt.ylabel('Amplitude (Units)')
        plt.xlabel('UTC Time (s)')        
        
        if save:
            show = False
        if show:
            plt.show()
        if save:
            # save the linearly stacked waveform
            ext_len = len(os.path.basename(self.st_path).split('.')[-1])
            outpath = os.path.basename(self.st_path)[:-ext_len]+'lin_stack.jpg'
            #print outpath
            fig.savefig(outpath)

    def phase_stack(self, st=None, norm=True):
        """
        The following function takes an obspy stream object and outputs 
        the phase stack of all traces within the stream. This is only useful
        if all traces are various channels from the same station, OR are from 
        the same cross-correlation station pairs for ambient studies. 
        """
        
        if st is None:
            st = self.st
        # if st is STILL None, then run the import stream!
        if st is None:
            self.import_stream()
            st = self.st
            
        PS = 0
        for tr in st:
            trace = tr.data - np.mean(tr.data)
            inst_phase = np.arctan2(trace, range(0,len(trace)))
            PS += np.real(np.exp(1j*inst_phase))

        N = len(st)
        if norm:
            #normalise the traces to between 0 and 1
            self.PS = np.abs(PS) / N / np.max(PS)
        else: 
            self.PS = np.abs(PS) / N
        
        #print "Phase Stack SNR: ", np.abs(np.mean(PS)) / np.abs(np.std(PS))

    def plot_phase_stack(self, PS=None, show=True, save=False):
        """
        Function to plot the phase stack of two or more traces within 
        the input Stream object for the class: Stack. Show is the default, 
        if save is True, then show is automatically set to False
        """
        
        # could add feature to allow for addition of own output path!        
        # could add feature to allow for higher resolution output images!
        
        if PS is None:
            PS = self.PS
                
        # if PS is STILL None, run the lin_stack function
        if PS is None:
            try:
                self.phase_stack()
                PS = self.PS            
            except Exception as error:
                print error
                fig = plt.figure()
  
        fig = plt.figure(figsize=(15,10))
        plt.xticks(rotation=10)
        plt.plot(self.times, PS, color='k')
        plt.title('Phase Stacked Waveform From File: {} \n {} - {}'\
        .format(os.path.basename(self.st_path), self.starttime, self.endtime))
        plt.xlabel('Nomalised Instantaneous Phase (2pi=1)')
        plt.ylabel('UTC Time')        
        
        if save:
            show = False
        if show:
            plt.show()
        if save:
            # save the linearly stacked waveform
            ext_len = len(os.path.basename(self.st_path).split('.')[-1])
            outpath = os.path.basename(self.st_path)[:-ext_len]+'phase_stack.jpg'
            #print outpath
            fig.savefig(outpath)      
            
            
    def pw_stack(self, st=None, LS=None, PS=None, norm=True, sharp_v=2):
        """
        The following function takes an obspy stream object and outputs 
        the phase-weighted stack of all traces within the stream. 
        This is only useful if all traces are various channels from the 
        same station, OR are from the same cross-correlation station pairs 
        for ambient studies. Default to the second power. Even powers work. 
        """
        if st is None:
            st = self.st
        # if st is STILL None, then run the import stream!
        if st is None:
            self.import_stream()
            st = self.st        
        
        
        if LS is None:
            LS = self.LS
                
        # if LS is STILL None, run the lin_stack function
        if LS is None:
            try:
                self.lin_stack()
                LS = self.LS
            except Exception as error:
                print error        
        
        if PS is None:
            PS = self.PS
                
        # if PS is STILL None, run the lin_stack function
        if PS is None:
            try:
                self.phase_stack()
                PS = self.PS            
            except Exception as error:
                print error
                


        PWS = (LS * PS ** sharp_v)
        PWS = PWS - np.mean(PWS)
        
        if norm:
            PWS = np.max(LS)/np.max(PWS) * PWS
        #print "Phase-Weighted Stack SNR: ", np.abs(np.mean(PWS)) / \
        #np.abs(np.std(PWS))

        self.PWS = PWS


    def plot_pw_stack(self, PWS=None, show=True, save=False):
        """
        Function to plot the phase-weighted stack of two or more traces within 
        the input Stream object for the class: Stack. Show is the default, 
        if save is True, then show is automatically set to False
        """
        
        # could add feature to allow for addition of own output path!        
        # could add feature to allow for higher resolution output images!
        
        if PWS is None:
            PWS = self.PWS
                
        # if PS is STILL None, run the lin_stack function
        if PWS is None:
            try:
                self.pw_stack()
                PWS = self.PWS            
            except Exception as error:
                print error
     
        fig = plt.figure(figsize=(15,10))
        plt.xticks(rotation=10)
        plt.plot(self.times, PWS, color='k')
        plt.title('Phase-Weighted Stacked Waveform From File: {} \n {} - {}'\
        .format(os.path.basename(self.st_path), self.starttime, self.endtime))
        plt.ylabel('Amplitude')
        plt.xlabel('UTC Time (s)')        
        
        if save:
            show = False
        if show:
            plt.show()
        if save:
            # save the linearly stacked waveform
            ext_len = len(os.path.basename(self.st_path).split('.')[-1])
            outpath = os.path.basename(self.st_path)[:-ext_len]+'pw_stack.jpg'
            #print outpath
            fig.savefig(outpath)  

if __name__ == '__main__':
    args = sys.argv
    
    if len(args) < 3:
        raise Exception('Please enter the path to the mseed waveforms files \
you want to stack and the type of operation you wish to perform!')

    operation = args[2]
    
    acceptible_ops = ['lin_stack', 'plot_lin_stack', 'phase_stack', 
                      'plot_phase_stack', 'pw_stack', 'plot_pw_stack']
    
    if operation not in acceptible_ops:
        raise Exception('Please choose an operation from the list: {}'.format(acceptible_ops))
                 
    if len(args) < 2:
        raise Exception('Please enter the path to the mseed waveforms files \
you want to stack.')
    # create arguments examples document
    
    # python stack.py path filter_waveform band_lims    
    
    if len(args) < 4:
        args.append(False)
    
    if len(args) < 5:
        args.append([1, 10])
    
    STACK = Stack(args[1], filter_waveform=args[3], band_lims=args[4])
            

    if operation == 'lin_stack':
        STACK.lin_stack()
    if operation == 'plot_lin_stack':
        STACK.plot_lin_stack()
    if operation == 'phase_stack':
        STACK.phase_stack()
    if operation == 'plot_phase_stack':
        STACK.plot_phase_stack()
    if operation == 'pw_stack':
        STACK.pw_stack()
    if operation == 'plot_pw_stack':
        STACK.plot_pw_stack()
