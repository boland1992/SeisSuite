# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 12:10:34 2015

@author: boland
"""

import numpy as np
from obspy.signal.trigger import *
import os, os.path
import datetime
import pickle
import shutil

def run_time(data_input, min_spread=90):
    """
    Function to use for parallel processing of all time spacings. Input is
    n mins e.g.
    min_spread = 2
    secs = min_spread * 60.
    n_mins = int((endtime - starttime)/60./min_spread)
    """
    
    starttime, time, path = data_input

    start = starttime + time*min_spread*60
    end = starttime + secs + time*min_spread*60
    
    print "Processing between {} and {}".format(start, end)
    st = read(path, starttime=start, endtime=end)

    st.select(component="Z")

    if len(st) > 0:
        tr = st[0]    
        AUTO_TRIGGER = Auto_trigger(tr, paths=path, 
                                    diff_thres=20.0, 
                                    eq_len=10.0)
                                                
        UTC_events = AUTO_TRIGGER.trigger_times(show=False)
                
        all_events = []
        if UTC_events is not None or UTC_events is []:
            if len(UTC_events) > 0: 
                
                            
                for event in UTC_events:
                                
                    event_start = event[0]; event_end = event[1]
                    

                    multi_st = AUTO_TRIGGER.multiplex_event(
                                event_start, paths=all_paths, 
                                breadth=30.)
                    #if multi_st is not None:
                    #    multi_st.plot()
                                
                                
                                    #AUTO_TRIGGER.verify_multi_events(multi_st)
                    event_string = "%s,%s,%s,%s\n" \
                                    %(tr.stats.network, 
                                      tr.stats.station, 
                                      str(event_start), 
                                      str(event_end))
                                    
                    if not os.path.exists("events"): 
                        os.makedirs("events")
                    all_events.append(event_string)
                    # save events to text file
                if len(all_events) > 0:
                    #print "writing to output text file"
                    event_file.writelines(all_events)
                        #print(all_events)
                        


######################
# SET INPUT PARAMETERS
######################

# set events output txt file name
event_fname = "event_outputs.txt"

t0 = datetime.datetime.now()

if not os.path.exists('tmp'): os.mkdir('tmp')
    
timeline_path = 'tmp/timeline_database.pickle'
# set up copy of event database to update and delete timelines as event search
# continues

event_timeline_path = 'tmp/event_timeline.pickle'
if not os.path.exists(event_timeline_path) and os.path.exists(timeline_path):
    # copy database to remove event times for future scans!
    shutil.copy(timeline_path, event_timeline_path)
                
if not os.path.exists(event_timeline_path):
    raise Exception("Please run stat_timeline.py to create the necessary \n\
timeline database for this programme.")


print "Importing event timeline database copy ... " 

f = open(event_timeline_path, mode='rb')
event_timeline = pickle.load(f)
f.close()

t1 = datetime.datetime.now()
print "Time taken to import timeline database pickle file was: ", t1-t0


# create list of all files in database from timeline database pickle file!
all_paths = []
for key in event_timeline.keys():
    for item in event_timeline[key]:
        if item[2] not in all_paths:
            all_paths.append(item[2])


# save all events to event_fname file location set above.
# Note that if the file already exists, the program will append events to it!
with open(event_fname, "a+") as event_file:
    # run for loop on all files!
    for key in event_timeline.keys():
        
        print key
        counter = 0
        for time_window in event_timeline[key]:
            starttime, endtime, path = time_window
            
            if counter == 0:
                print "\nProcessing file: {}".format(os.path.basename(path)) 

            #try:

                #read and attached responce with metadata() function
                #st = metadata(f)
            
                #split into min_spread*minute long scans
            secs = min_spread * 60.
            n_mins = int((endtime - starttime)/60./min_spread)
            
            parallel_run = np.column_stack(((n_mins-1)*[starttime], 
                                                range(0, n_mins-1), 
                                                (n_mins-1)*[path])) 
            
            
            for data_input in parallel_run:
                try:
                    run_time(data_input)
                except Exception as error:
                    print error
                #pool = mp.Pool(None)
                #pool.map(run_time, parallel_run)
                #pool.close()
                #pool.join()
                
            #except Exception as error:
            #    print error
            counter += 1
            
        event_timeline_path = 'tmp/event_timeline.pickle'
            
        if not os.path.exists(event_timeline_path):
            raise Exception("Please make sure that both stat_timeline.py\
has been run, and that event_timeline.pickle is in the local tmp directory.")

        print "Importing new event timeline database \
from previous pass ... " 

        f = open(event_timeline_path, mode='rb')
        event_timeline = pickle.load(f)
        f.close()
        
print("file closed")

#------------------------------------------------------------------------------
# NOTES OF THINGS STILL TO COMPLETE BEFORE USE!
#------------------------------------------------------------------------------

# automate time difference between earthquakes to be larger or smaller based
# on earthquake size! ask Abe about some sort of equation. 

# Make sure that if the signal is in fact too noisy, the attempt a band-pass
# filter and check to see if that can find anything. Set limits based on 
# earthquake literature! 

#------------------------------------------------------------------------------
# DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE 
#------------------------------------------------------------------------------

#MAKE SURE YOU ADD IN A MINIMUM EARTHQUAKE TIME OCCURRENCE!

# set a lower limit on the trigger maximum because the signal will give 
# loads of noise if there isn't a single earthquake in the time series! 

# set event time min and time max +- a certain automated amount of time as list
# then go through and save mseed of each event labelling it 
# 0 through to whatever the event number is! 

#NORMALISE THESE POINTS ABOUT ZERO!

