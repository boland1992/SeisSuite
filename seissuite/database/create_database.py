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
from seissuite.misc.path_search import paths
import pickle

RANK = False

try:
    import cPickle as pickle
except:
    import pickle
    print "Caution, database code may run slow due to cPickle failed import"
    
# import CONFIG class initalised in ./configs/tmp_config.pickle
config_pickle = 'configs/tmp_config.pickle'
f = open(name=config_pickle, mode='rb')
CONFIG = pickle.load(f)
f.close()
    
# import variables from initialised CONFIG class.
AUTOMATE = CONFIG.AUTOMATE
MSEED_DIR = CONFIG.MSEED_DIR
DATABASE_DIR = CONFIG.DATABASE_DIR
DATALESS_DIR = CONFIG.DATALESS_DIR
STATIONXML_DIR = CONFIG.STATIONXML_DIR


def check_stat(trace, checklist=[], map_dict={}):
    """
    Function written to map a given station name to a new station name if
    it is in a given checklist. 
    
    Args:
    station (obspy.station.Station): Obspy station object containing 
                                     all header information.
    checklist (list): list of station names that should be changed if the 
                      station name is found to be in there.
    map_dict (dictionary): dictionary with keys that are the strings of the 
                           old station names, which store information related
                           to the new station names. 
    Returns:
    station (obspy.station.Station): Obspy station object containing 
                                     all header information.                    
    """
    
    stat_name, net_name = trace.stats.station, trace.stats.network
    

        
    if stat_name in checklist:
        try:
            trace.stats.station = map_dict[stat_name]
            
            if net_name == 'XX':
                # set new network name to just be the first two letter/numbers
                # of the station names
                trace.stats.network = trace.stats.station[:2]

            return trace

        except Exception as error:
            print error


    
    

UNAM_statlist = ['C0200', 'C0230', 'C0120', 'C01E0', 
                'C0180', 'C01B0', 'C0220', 'C00D0', 
                'C00E0', 'COOA0', 'C0240']

UNAM_statmap = {'C0200':'AM01', 'C0230':'AM03', 'C0120':'AM04', 'C01E0':'AM05', 
                'C0180':'AM06', 'C01B0':'AM10', 'C0220':'AM12', 'C00D0':'AM13', 
                'C00E0':'AM15', 'COOA0':'AM16', 'C0240':'AM17'}



multiprocess = False

if multiprocess:
    import multiprocessing as mp


t_total0 = datetime.datetime.now()


    
folder_path = MSEED_DIR

# set file extension
#extensions = ['m', 'mseed', 'miniseed', 'MSEED']
extensions = ['msd']

 
abs_paths = []
for extension in extensions:
    abs_paths.append(paths(folder_path, extension))

#flatten the list of absolute paths containing all relevant data files
abs_paths = list(itertools.chain(*abs_paths))

# initialise timeline dictionary, one key per station!
global timeline
timeline = {}

if not os.path.exists('tmp'): os.mkdir('tmp')


# =============================================================================
# SQL
# =============================================================================

# create database if it doesn't exist already, if it does, stop the programme.
database_name = os.path.join(DATABASE_DIR, 'timeline.db')

resp_db = os.path.join(DATABASE_DIR, 'response.db')

if not AUTOMATE:
    if os.path.exists(database_name):
        yeses = ['y','Y','yes','Yes','YES']    
        nos = ['n','N','no','No','NO']    
    
        condition = False
        while condition is False:
    
            answer = raw_input('Would you like to remove the existing database\
 and start afresh? (y/n): ')
    
            if answer in yeses:
                os.remove(database_name)
                condition = True
            
            elif answer in nos:
                raise Exception("The SQL database {} already exists, \
quitting programme.".format(database_name))
                condition = True

            else:
                print "The input answer must be of yes or no format."
                condition = False
                
else:
    os.remove(database_name)

# =============================================================================
# =============================================================================

<<<<<<< HEAD
def extract_info(info, check=False):    
    trace, path = info
    
    # check that the station names fit within the checklist, else map them.
    
    if check:
        trace = check_stat(trace, checklist=UNAM_statlist, 
                               map_dict=UNAM_statmap)
    
=======



def extract_info(info):    
    trace, path = info
    
    stats = trace.stats
    code ='{}.{}.{}'.format(stats.network, stats.station, stats.channel)
    
    starttime = trace.stats.starttime.timestamp
    endtime = (trace.stats.starttime + trace.stats.npts * \
              (1/trace.stats.sampling_rate)).timestamp

    return (code, starttime, endtime, path)
>>>>>>> 561db556e4ab402ce7b410117402acd2170b7722

    try:
        stats = trace.stats
        code ='{}.{}.{}'.format(stats.network, stats.station, stats.channel)
        print code,
        starttime = trace.stats.starttime.timestamp
        try:
            endtime = trace.stats.endtime.timestamp
        except:    
            endtime = (trace.stats.starttime + trace.stats.npts * \
            (1.0/trace.stats.sampling_rate)).timestamp

        information = (code, starttime, endtime, path)
        return information
    except:
        a=None
        
def info_from_headers(path):
<<<<<<< HEAD

    #print os.path.basename(path)
    try:
        #t0 = datetime.datetime.now()
        headers = read(path, headonly=True)
        headers.select(component='Z')
        info = []
        
        for trace in headers:
            # double check that we are only dealing with the Z channel 
            if 'Z' in trace.stats.channel:
                info.append([trace, path])
        
        timeline_header = map(extract_info, info)
        return timeline_header
        
        
    except Exception as error:
        #t0 = datetime.datetime.now()
        try:
            headers = read(path)

            headers.select(component='Z')
            info = []
            for trace in headers:
                info.append([trace, path])
        
            timeline_header = map(extract_info, info)
    
            return timeline_header

        except Exception as error:
            print error
=======
    try:
        #t0 = datetime.datetime.now()
        print os.path.basename(path)
        headers = read(path, headonly=True)
        headers.select(component='Z')
        info = []
    
        for trace in headers:
            #if code in ranked_list:
            info.append([trace, path])
        
        timeline_header = map(extract_info, info)
    
        return timeline_header
    
    except Exception as error:
        print error    
>>>>>>> 561db556e4ab402ce7b410117402acd2170b7722

    #t1 = datetime.datetime.now()
    #print 'time taken to process previous loop: ', t1-t0
    #print "time for previous loop was: ", t1-t0
    
t0 = datetime.datetime.now()

print "Initialising timeline database. Please be patient ... " 

if multiprocess:
    pool = mp.Pool(None)
    timeline = pool.map(info_from_headers, abs_paths)
    pool.close()
    pool.join()
else:
    timeline = map(info_from_headers, abs_paths)
<<<<<<< HEAD

#
#timeline = np.array([i for i in timeline if all(i)])

#flatten output list and  remove empty lists from timelines database
timeline = np.asarray(list(itertools.chain(*timeline)))
=======
    
try:
    
    #flatten output list
    timeline = np.asarray(list(itertools.chain(*timeline)))

except:
    print "The timeline array is not the correct type and cannot be flattened"
>>>>>>> 561db556e4ab402ce7b410117402acd2170b7722

# remove None's from timeline array. 
timeline = timeline[timeline != np.array(None)]


t1 = datetime.datetime.now()
print "time taken to read in and process timeline database: ", t1-t0


# =============================================================================
# USING SQL
# =============================================================================

    
print 'creating SQL database ... '
conn = lite.connect(database_name)
c = conn.cursor()
print 'SQL database cursor initialised'

try:
    # Create table called timeline for timeline database
    # the table will contain the station ID which is net.stat_name
    # the start time and end time for each trace and the absolute file path loc
    # the goal is to have one line per trace!
    c.execute('''CREATE TABLE timeline
                 (station text, starttime real, 
                 endtime real, file_path text)''')
    print 'timeline table created'

    # create another table called file extrema
    # this table contains the station code, the file path and the start and
    # end times of that specific table!
    # one line per file path!
    c.execute('''CREATE TABLE file_extrema
                 (station text, starttime real, 
                 endtime real, file_path text)''')
    print 'file extrema created'
    
except Exception as error:
    print error             


# do an error check on the timeline list
print len(timeline)



t0 = datetime.datetime.now()
# commit many rows of information! 
c.executemany('INSERT INTO timeline VALUES (?,?,?,?)', timeline)

t1 = datetime.datetime.now()
conn.commit()

print "Time taken to commit rows to SQL db: ", t1-t0

# create another table called file extrema
# this table contains the station code, the file path and the start and
# end times of that specific table!
# one line per file path!

try:
    c.execute('''CREATE TABLE file_extrema
                 (station text, starttime real, 
                  endtime real, file_path text)''')
    print 'file extrema table created'

except Exception as error:
    print error
    

paths = c.execute('''SELECT DISTINCT file_path FROM timeline''')

# find out how many unique paths there are!
global paths_list
paths_list = []
def append_fast(list_):
    paths_list.append(str(list_[0]))    
map(append_fast, paths)
# convert to tuple
paths_list = tuple(paths_list)

file_extrema = []
for path in paths_list:
    path = (path,)

    row_extrema = c.execute('''SELECT station, MIN(starttime), 
                            MAX(endtime), file_path FROM timeline 
                            WHERE file_path=?''', path)   

    for info in row_extrema:
        #if RANK: 
        code, min_time, maxtime, file_path = info
        net, stat, chan = code.split('.')
        file_extrema.append(info)
  
file_extrema = tuple(file_extrema)

c.executemany('INSERT INTO file_extrema VALUES (?,?,?,?)', file_extrema)

conn.commit()

conn.close()