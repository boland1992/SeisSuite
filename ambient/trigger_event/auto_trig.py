# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 12:10:34 2015

@author: boland
"""

import numpy as np
from obspy import read 
from obspy.signal.trigger import *
import os, os.path
from obspy.station import stationxml
import matplotlib.pyplot as plt
from obspy.core import Stream
import pickle
import glob

class Auto_trigger: 
    
    def __init__(self, tr, paths=None, diff_thres=90.0, eq_len=10.0):
        self.tr = tr
        self.diff_thres = diff_thres
        self.eq_len = eq_len
        self.paths = paths
        
    def area_cond(self, tr=None):
        """
        Function that returns a 95% noise to signal area ratio condition. 
        This sorts the trace by value, returns the area beneath the graph
        before the 95% noise line and then computes the area ratio of this
        over the total rectangular area beneath the 99.9% trace value.
        """
        if tr is None:
            tr = self.tr
        
        #print "area cond"
        #print tr.stats.starttime
        #print tr.stats.starttime + tr.stats.npts*(1/tr.stats.sampling_rate)
        avg_tr = np.average(tr)
        shifted_tr =  tr - avg_tr   
        abs_tr = abs(shifted_tr)
        sort_tr =  np.sort(abs_tr)
        tr_99 = sort_tr[int(0.99*len(sort_tr))]
        area_big = tr_99 * int(0.95*len(sort_tr))
        area_small = np.sum(sort_tr[:int(0.95*len(sort_tr))])
        area_condition = area_small /  area_big
        #print "Noise to signal area condition: %0.3f"%(area_condition)
        return area_condition
    
    def signal_cond(self, tr=None, signal='clean', check=False):
        """
        Function inputs a certain trace (as determined by the area condition)
        a trace's characteristic function, a set trigger on and trigger off
        value as determined by being 40% between the ctf_low and 99% ctf 
        values. Three signal types: 'clean', 'noisy', 'too noisy'. If signal 
        is quite noisy, use recursive STALTA, else use classic STALTA.
        """
        if tr is None:
            tr = self.tr
        
        #print "signal cond"
        #print tr.stats.starttime
        #print tr.stats.starttime + tr.stats.npts*(1/tr.stats.sampling_rate)
        
        df = tr.stats.sampling_rate
        
        if signal == 'clean':
            #print("CLEAN SIGNAL!\n")
            #grad = np.gradient(self.tr.data)

            ctf = abs(1-carlSTATrig(tr.data, int(5*df), 
                                    int(10*df), 1e-4, 1e-4))
            #ctf = abs(1-recSTALTA(self.tr.data, int(5 * df), int(10 * df)))

            low_c = 0.8
            above_mean = 0.5
         
        elif signal == 'noisy':        
            #print("NOISY SIGNAL!\n")
            ctf = abs(1-carlSTATrig(tr.data, int(5*df), 
                                    int(10*df), 1e-4, 1e-4))
            #ctf = abs(1-recSTALTA(self.tr.data, int(5 * df), int(10 * df)))
            low_c = 0.9
            above_mean = 0.5
            
        if check:
            above_mean = 0.9
        
        # take 4th kurtosis gradient!
        ctf = abs(np.gradient(ctf))
        #ctf = abs(np.gradient(ctf))
        #ctf = abs(np.gradient(ctf))
        #ctf = abs(np.gradient(ctf))
        
        # remove anomalies
        ctf[975:1025]  = 0.

        #===========================
        sort_ctf = np.sort(ctf)
        ctf_low = sort_ctf[int(low_c*len(sort_ctf))]
        ctf_99 = sort_ctf[int(0.999*len(sort_ctf))]
        max_trig = ctf_99
        
        
        trig_on = ctf_low + above_mean * (max_trig - ctf_low)
        trig_off = 0.95*trig_on

        if not max_trig > 1.25 * ctf_low:
            trig_on, trig_off, ctf = False, False, False
            
        return trig_on, trig_off, ctf

    def remove_n(self, x, y1, n):
        """
        This function removes every nth instance in a list
        """
        y = []
        new_x = np.linspace(x[0], x[-1], len(x)/n)
        for i in range(0, len(x)-1):
            if i % n == 0:
                y.append(y1[i])
        
        y = np.abs(y)
        
        plt.figure()
        plt.plot(new_x, y)
        plt.show()
        y1 = y
        #self.regress(x, y)

    def regress(self, x, y1):
        
        y = [y1[0], y1[1], y1[2], y1[3]]
        for i in range(0, len(x)-4):
            y.append((y1[i+4] - y1[i])/4)
        
        y = np.abs(y)
        
        plt.figure()
        plt.plot(x, y)
        plt.show()
        y1 = y
        self.regress(x, y)
        
    def trig_conditions(self, tr=None, counter=0, check=False):
    
        """
        Function outputs a trigger value based on maximum of characteristic 
        function for a given waveform series. This means that if there is a 
        larger than normal earthquake that occurs during this waveform, it'll 
        take that as peak and the trigger detections won't be thrown off 
        because every detection is based on the maximum amplitude during 
        that imported wave-series. 
    
        The key is to find the optimum seismic waveform time-series length to 
        detect as many earthquakes as possible.
        """
        if tr is None:
            tr = self.tr
        
        #print "trig conditions"
        #print tr.stats.starttime
        #print tr.stats.starttime + tr.stats.npts*(1/tr.stats.sampling_rate)
        ##############
        #TRACE RELATED
        ##############

        tr = tr.filter('bandpass', freqmin=1.0, freqmax=10.0, 
                       corners=2, zerophase=True)
        
        # normalise the data between zero and one. 
        tr.data = np.asarray([point/np.max(tr.data) for point in tr.data])

        area_condition = self.area_cond(tr=tr)

        #noise_parameter = np.max(tr.data)/(np.std(tr.data)/np.sqrt(len(ctf)))
        #print "noise: ", noise_parameter
        
        #print "SNR: ", abs(np.std(tr.data)/np.mean(tr.data))
        
        
        #for i in ctf:
        
        #x, y = np.linspace(0, len(ctf), len(ctf)), ctf
        
        #self.remove_n(x, y, 2)
        #self.regress(x,y)
        #for i in range(0, len())

        #=====================================================================
        #polyfit

        #x, y = np.linspace(0, len(ctf), len(ctf)), np.asarray(ctf)
        # calculate polynomial
        #z = np.polyfit(x, y, 10)
        
        #f = np.poly1d(z)

        # calculate new x's and y's
        #x_new = np.linspace(x[0], x[-1], len(ctf))
        #y_new = f(x_new)    
        #plt.plot(x_new, abs(y_new))
        #plt.xlim([x[0]-1, x[-1] + 1 ])
        #plt.show() 
        
        #ctf = y_new
        #=====================================================================
       
        if area_condition > 0.15:
           # print("TRACE POTENTIALLY TOO NOISY TO DETECT SIGNAL!")
             

        #     print("ATTEMPTING TO CLEAN SIGNAL ... RUNNING BAND-PASS FILTER")
        #     tr = tr.filter('bandpass', freqmin=0.1, freqmax=20.0, 
        #                              corners=2, zerophase=True)
             # only run two instances of the filter to see if the program 
             # can find a signal!
        #     if counter <= 1:    
        #         counter+=1
        #         self.trig_conditions(tr=tr, counter=counter)
            trig_on, trig_off, ctf = False, False, False
             
        elif 0.1 <= area_condition <= 0.15:
            trig_on, trig_off, ctf = self.signal_cond(tr=tr, 
                                                      signal='noisy',
                                                      check=check)
                                                      
        else:
            trig_on, trig_off, ctf = self.signal_cond(tr=tr, 
                                                      signal='clean',
                                                      check=check)
                                                      

        
        return trig_on, trig_off, ctf


    def UTC_times(self, times, tr=None):
        """
        Function that takes the output of the obspy.signal.trigger function 
        triggerOnset() and computes an average for points that do no differ 
        by a certain difference threshold (diff_thres), and gives back any 
        individual points that differ by more than this threshold on either 
        side.
    
        Finally the function converts this output into seconds and then 
        adds those on to the starting time of the trace. 
    
        diff_thres - set min. number of seconds between signals. 90s default.
        eq_len - set min. length for a signal to be counted. 10s default. 
        """
        # set initial variables
        if tr is None:
            tr = self.tr
            
        #print "UTC_times"
        #print tr.stats.starttime
        #print tr.stats.starttime + tr.stats.npts*(1/tr.stats.sampling_rate)
        
        start_time = tr.stats.starttime
        times = times / tr.stats.sampling_rate
        #remove unwanted parts of times numpy array 
        times = times[:,0]
        event_times = []
        event = [times[0]]

        for i in range(1, len(times)):
            # check if two events in times array have a difference < diff_thres, 
            #if not, run average of those times, if so append that events to a 
            #new events_times list
            time_diff = times[i] - times[i-1]
            #save info until events are far enough apart! 
            if time_diff < self.diff_thres:
                event.append(times[i])
            #raise conditional for if events are far enough apart! 
            else:
                event_start = event[0] - 2 #minus 5 seconds
                event_end = max(event) + 2 #add 5 seconds
                event_times.append([event_start, event_end])    
                event = []      
                event.append(times[i])
        #if event still contains something for any reason, add it to event times
        if len(event) > 0:            
            event_start = event[0] - 2 #minus 5 seconds
            event_end = max(event) + 2 #add 5 seconds
            event_times.append([event_start, event_end])
            event = [] 
        
        UTC_events = []
        # check if the event that is being saved fits the min. length criteria
        for i in event_times:
            estart = start_time + i[0]
            eend = start_time + i[1]
            if eend - estart > self.eq_len:
                UTC_events.append([estart, eend])
    
        return UTC_events


    def trigger_times(self, tr=None, show=False, check=False):
        """
        Function that returns UTCdatetime objects of events from a given input 
        trace. It uses carlSTATrig to find the characteristic function of the 
        trace, then plots the trigger findings using the plotTrigger and finally
        the function UTC_events is used to compute the event timings in UTC 
        format from the ouput of the function triggerOnset.
        """
        
        if tr is None:
            tr = self.tr
        
        #print "trigger_times"
        #print tr.stats.starttime
        #print tr.stats.starttime + tr.stats.npts*(1/tr.stats.sampling_rate)
      
        trig_on, trig_off, ctf = self.trig_conditions(tr=tr, check=check)
        tuple_check, tuple_false = (trig_on, trig_off, ctf), \
                                   (False, False, False)
        if tuple_check is tuple_false:
            return None
        else:
            if show:
                # use plotting the trigger functions using plotTrigger
                plotTrigger(tr, ctf, trig_on, trig_off, show=True) 
                
            times = triggerOnset(ctf, trig_on, trig_off, 
                                 max_len=self.diff_thres)

            if len(times) > 0:
                UTC_events = self.UTC_times(times, tr=tr)
                return UTC_events
        
    def verify_events(self, tr=None):
        """
        Function to take the input from the output of function multiplex_event
        and then determine whether or not there is an event in each of the
        traces in the multiplexed stream. 
        """
        condition = False
        
        if tr is None:
            tr = self.tr
            
        trig_on, trig_off, ctf = self.trig_conditions(tr=tr, check=True)
        #plotTrigger(tr, ctf, trig_on, trig_off, show=False) 

        times = triggerOnset(ctf, trig_on, trig_off, 
                             max_len=self.diff_thres)
                             
        if len(times) > 3:
            condition = True
            
        return condition
        
    def multiplex_event(self, event, paths=None, breadth=60.):
        """
        This function simply takes one breadth (seconds) either side of 
        an event and combines all traces in the database. Event is an
        obspy UTCDateTime object of the start time of a single representative
        event. 
        """
        if self.paths is None: 
            raise Exception('In order to combine traces, \
must import database paths')
        elif paths is None:
            paths = self.paths 
            
        event_start = event - breadth
        event_end = event + breadth
        

        traces = []
        for path in paths:
            # check that both event start and event end are in the databases time window!
            
            st_headers = read(path, headonly=True)
            
            tr_headers = st_headers[0]
            
            path_start = tr_headers.stats.starttime
            
            path_end = tr_headers.stats.starttime + \
                              tr_headers.stats.npts *      \
                              (1/tr_headers.stats.sampling_rate)
            

            # check to see if the time window is within the read-in file.
            if path_start <= event_start < path_end or \
            path_start <= event_end < path_end:

                st = read(path, starttime=event_start, endtime=event_end)
                tr = st[0]
                condition = self.verify_events(tr=tr)
                #append only traces that have a possibility of an event
                if condition:
                    traces.append(tr)

        # set minimum of 2 station occurrences for an event to be counted
        if len(traces) >= 2:
            STREAM =  Stream(traces)
            #STREAM.plot()
            #check to see if there is an events folder already:
            if not os.path.exists('events'): os.mkdir('events')
                
            # scan for previous event files! 
            events_list = sorted(glob.glob('events/event*.mseed'))

            if len(events_list) > 0:
                #find the number of the latest saved event.
                latest_event = events_list[-1].split('.')[0][-1]
                mseed_out = 'events/event{}.mseed'.format(int(latest_event)+1)
            
            else:
                # set unique output for each event 
                mseed_out = 'events/event0.mseed'
            
            # save event to individual mseed file! 
            STREAM.write(mseed_out, format='MSEED')
            

            
            #scan all other station's timelines to remove event time window
            #the following code removes event window from time series database
            #plt.figure()

            for key in event_timeline.keys():

                for i in range(0, len(event_timeline[key])-1):

                    item = event_timeline[key][i]

                    #case 1:
                    if item[0] <= event_start <= item[1] and not\
                    item[0] <= event_end <=  item[1]:
                        # set new end to database window as event_start
                        item[1] = event_start                    
                    
                    #case 2:
                    if item[0] <= event_end <= item[1] and not\
                    item[0] <= event_start <=  item[1]:
                        item[0] = event_end               
                        
                    #case 3:
                    if item[0] <= event_start <= item[1] and\
                    item[0] <= event_end <=  item[1]:
                        # split apart the old item such that the event time is
                        # erased!
                        
                        new_item1 = [item[0], event_start, path]
                        new_item2 = [event_end, item[1], path]
                        
                        #append new items to database                        
                        event_timeline[key] = list(event_timeline[key])
                        event_timeline[key].append(new_item1)
                        event_timeline[key].append(new_item2)
                        # delete ith index of event_timeline dict if case 3                        
                        del event_timeline[key][i]
                        #don't forget to re-sort the database!
                        event_timeline[key] = np.sort(np.asarray(
                                                      event_timeline[key]), 
                                                      axis=0)


                # trying to plot timeline to check for problems            
                #k = np.asarray(event_timeline[key])

                #for i in k:
                #    plot_stuff = np.asarray([[i[0].timestamp, 1], 
                #                             [i[1].timestamp, 1]])
                         
                #    plt.plot(plot_stuff[:,0], plot_stuff[:,1])

            #plt.show()
                                      
            # dump the new updated database to be read in on the next pass 
            with open(event_timeline_path, 'wb') as f:
               #print "\nExporting new event timeline database to: " + f.name
                pickle.dump(event_timeline, f, protocol=2)
        

            return STREAM
        else:
            return None
        
def metadata(read_dir):
    try:
        st = read(read_dir)
        #print(st)
        x = len(st)
    #CORRECTING FOR NETWORK CODE AND LOCATION
        network_change1 = 'UM'
        network_new_name1 = 'UM'#'BM'
        network_change2 = '01'
        network_new_name2 = 'UM'
        network_change3 = 'A'
        network_change4 = ''
        location_blank = ''
        for i in range(0, x):
            tr = st[i]
    # removes LOCATION so it is blank, as listed in the metadata 
    # files (regardless of what it was previously)
            tr.stats["location"] = location_blank
    # Changes BOREHOLE network codes from UM to BM and SURFACE 
    #network codes from 01 to UM 
            net = tr.stats["network"]
            if network_change1 in net:
                tr.stats["network"] = network_new_name1
            elif network_change2 or network_change3 or \
            network_change4 in net:
                tr.stats["network"] = network_new_name2
            else:
                continue
    #CORRECTING BOREHOLE STATION NAMES
        serial_no_1 = 'A346'
        site_name_1 = 'LOYU'
        serial_no_2 = 'BD5E'
        site_name_2 = 'MOSU'
        serial_no_3 = 'BD70'
        site_name_3 = 'SGWU'
        serial_no_4 = 'BD91'
        site_name_4 = 'WILU'

    #Changes station name from serial number to station code

        for i in range(0, x):
            tr = st[i]
            stat = tr.stats["station"] 
            if serial_no_1 in stat:
                tr.stats["station"] = site_name_1
            elif serial_no_2 in stat:
                tr.stats["station"] = site_name_2
            elif serial_no_3 in stat:
                tr.stats["station"] = site_name_3
            elif serial_no_4 in stat:
                tr.stats["station"] = site_name_4
            else:
                continue
        
        # CHANGES TO CHANNEL CODE

    # (this is a bit messy at the moment since the wildcard 
    #feature seemed to be failing)
        
        channel_new_name_E = 'EHE'
        channel_new_name_N = 'EHN'
        channel_new_name_Z = 'EHZ'
        
        for i in range(0, x):
            tr = st[i]
            chan = tr.stats["channel"] 

    # Changes CHANNEL names from '**E', '**N', '**Z', (e.g. BHE, DHZ) 
    # to a consistant format of EHE, EHN, EHZ
    # EXCEPT FOR BOREHOLE STATIONS, which will maintain channel codes
    # BHE, BHN, BHZ

            if 'DHE' in chan:
                tr.stats["channel"] = channel_new_name_E
            elif 'DHN' in chan:
                tr.stats["channel"] = channel_new_name_N
            elif 'DHZ' in chan:
                tr.stats["channel"] = channel_new_name_Z    
            elif 'ENE' in chan:
                tr.stats["channel"] = channel_new_name_E
            elif 'ENN' in chan:
                tr.stats["channel"] = channel_new_name_N
            elif 'ENZ' in chan:
                tr.stats["channel"] = channel_new_name_Z         
            else:
                continue        
        # saves stream as a combination of edited traces       
            st[i] = tr
        # ATTACH METADATA TO STREAM
        #import and read metadata file from GitHub
        metadata_path = 'UOM.xml'
        if not os.path.isfile(metadata_path):
            os.system('wget https://raw.githubusercontent.com/\
unimelb-geophys/metadata/master/UOM.xml')
        metadata = stationxml.read_StationXML(metadata_path)
        #print(metadata)
        #st.write("")
        #attach metadata to stream
        st.attach_response(metadata)
        #st.remove_response()
        #delete downloaded metadata file 
        #os.system('rm -r UOM.xml')
        #st.write(read_dir, format="MSEED")
    except Exception as e:
        print("\nOops, there was a problem with attaching the metadata\n")
        print(e)
    return st
    
min_spread = 2
    
