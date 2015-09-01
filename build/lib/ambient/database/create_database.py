# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 11:10:38 2015

@author: boland
"""


from obspy import read
import numpy as np
import itertools
import os
import datetime
import sqlite3 as lite

try:
    import cPickle as pickle
except:
    import pickle
    print "Caution, database code may run slow because cPickle failed to import"

#from pysismo.psconfig import MSEED_DIR

multiprocess = False

if multiprocess:
    import multiprocessing as mp


t_total0 = datetime.datetime.now()

#SCANNING FUNCTIONS 
def paths_sort(path):
    """
    Function defined for customised sorting of the abs_paths list
    and will be used in conjunction with the sorted() built in python
    function in order to produce file paths in chronological order.
    """
    base_name = os.path.basename(path)
    stat_name, date= base_name.split('.')[0], base_name.split('.')[1]     
    try:
        date = datetime.datetime.strptime(date, '%Y-%m-%d')
        return date, stat_name
    except Exception as e:
        print(e)
        
def paths(folder_path, extension, sort=False):
    """
    Function that returns a list of desired absolute paths called abs_paths
    of files that contains a given extension e.g. .txt should be entered as
    folder_path, txt. This function will run recursively through and find
    any and all files within this folder with that extension!
    """
    abs_paths = []
    for root, dirs, files in os.walk(folder_path):
        for f in files:
            fullpath = os.path.join(root, f)
            if os.path.splitext(fullpath)[1] == '.{}'.format(extension):
                abs_paths.append(fullpath)
    if sort:
        abs_paths = sorted(abs_paths, key=paths_sort)
        
    return abs_paths
    

# set folder with all waveform files in it. Can be recursive! 

#folder_path = 'small_borehole_quakes'
#folder_path = '/storage/ABE/borehole_data'
#folder_path = '/storage/ANT/INPUT/DATA/AGOS-FULL-DATA_LOWSAMPLE'


folder_path = 'small_borehole_quakes'

# set file extension
extensions = ['m', 'mseed', 'miniseed', 'MSEED']
 
# set desired component e.g. E, N or Z

abs_paths = []
for extension in extensions:
    abs_paths.append(paths(folder_path, extension))
#['/home/boland/Dropbox/University/UniMelb/Research/SIMULATIONS/Triggers/chch_earthquake.MSEED']
#flatten the list of absolute paths containing all relevant data files
abs_paths = list(itertools.chain(*abs_paths))
# initialise timeline dictionary, one key per station!
timeline = {}

if not os.path.exists('tmp'): os.mkdir('tmp')

LENGTH = len(abs_paths)
counter = 0
loop_times = []
try: 
    for path in abs_paths:
        
        t0 = datetime.datetime.now()
        headers = read(path, headonly=True)
        headers.select(component='Z')

    
        for trace in headers:
            code ='{}.{}'.format(trace.stats.network, trace.stats.station)
            starttime = trace.stats.starttime
            endtime = trace.stats.starttime + trace.stats.npts * \
                (1/trace.stats.sampling_rate)
            if code not in timeline.keys():
                timeline[code] = []
                timeline[code].append([starttime, endtime, path])
            else:
                timeline[code].append([starttime, endtime, path])
        
        t1 = datetime.datetime.now()
        #print 'time taken to process previous loop: ', t1-t0
        loop_times.append((t1-t0).total_seconds())
        avg_time = np.average(loop_times)        
        counter += 1
        loops_remaining = LENGTH - counter
        print "loops remaining: ", loops_remaining
        print "time for previous loop was: ", t1-t0
        #time_remaining = avg_time * loops_remaining
        #mins_remaining = int(time_remaining / 60)
        #seconds_left = time_remaining % 60
        
        #print "estimated processing time remaining: %d mins %0.2f secs" \
        #%(mins_remaining, seconds_left)
        
        
        
except Exception as error:
    print error
#pool = mp.Pool(4)
#info = pool.map(scan_path, abs_paths)
#pool.close()
#pool.join()

# sort the dictionary database using numpy
for key in timeline.keys():
    timeline[key] = np.sort(np.asarray(timeline[key]), axis=0)

# =============================================================================
# USING SQL
# =============================================================================

# create database if it doesn't exist already, if it does, stop the programme.
database_name = 'timeline_database.db'

if os.path.exists(database_name): 
    raise Exception("The SQL database {} already exists, quitting programme."\
    .format(database_name))

conn = lite.connect(database_name)
c = conn.cursor()

# Create table called timeline for timeline database
# the table will contain the station ID which is net.stat_name
# the start time and end time for each trace and the absolute file path loc
# the goal is to have one line per trace!

try:
    c.execute('''CREATE TABLE timeline
                 (stat_name, starttime, endtime, file_path)''')
except Exception as error:
    print error             
                 
for key in timeline.keys():
    try:
        stat_name = key        
        for stat_info in timeline[key]:            
            start, end, file_path = stat_info
            example_row = '\'{}\',\'{}\',\'{}\',\'{}\''.format(stat_name, 
                                                               start,
                                                               end, 
                                                               file_path)
                                                               
            #print example_row
    
            # Insert a row of data
            c.execute("INSERT INTO timeline VALUES ({})".format(example_row))

            # Save (commit) the changes
            conn.commit()
        
        
    except Exception as error: 
        print error
    
    
for row in c.execute('SELECT * FROM timeline ORDER BY stat_name'):
    print row


# =============================================================================
# USING PICKLE
# =============================================================================

#outfile = 'tmp/timeline_database.pickle'

#if not os.path.exists(outfile):
#    with open(outfile, 'wb') as f:
#        print "\nExporting new timeline database to: " + f.name
#        pickle.dump(timeline, f, protocol=2)
#else:
#    raise Exception("Could not create new pickle database as one already exists")
        
#t_total1 = datetime.datetime.now()
#print 'Total time taken to initialise timeline database for {} was: {}'.format(
#            os.path.basename(folder_path), t_total1-t_total0)