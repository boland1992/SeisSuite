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
from pyseis.modules.rdreftekc import rdreftek, reftek2stream

def read_ref(path):
    ref_head, ref_data = rdreftek(path)
    st = reftek2stream(ref_head, ref_data)
    return st


# mapping old station names to new station names
statmap = {'A324':'W01', 'A320':'W02', 'AA09':'W03', 'A333':'W04', 
           'A971':'W05', 'A356':'W06', 'A316':'W25', '9994':'W08', 
           'B15A':'W26','B2FA':'W10', 'B205':'W11', 'A276':'W19', 
           'A956':'W24', 'B205':'W28', 'ALRZ':'W80', 'HRRZ':'W82', 
           'GRRZ':'W81', 'RITZ':'W97', 'WATZ':'W98', 'WHTZ':'W99',
           'KRVZ':'W83', 'OTVZ':'W84', 'PRRZ':'W85', 'WPRZ':'W86', 
           'HSRZ':'W87', 'MRHZ':'W88', 'ARAZ':'W90', 'HATZ':'W91', 
           'HITZ':'W92', 'KATZ':'W93','KUTZ':'W94','POIZ':'W95', 
           'RATZ':'W96' }



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


def get_filepaths(directory):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.    
    

#UNAM_statlist = ['C0200', 'C0230', 'C0120', 'C01E0', 
#                'C0180', 'C01B0', 'C0220', 'C00D0', 
#                'C00E0', 'COOA0', 'C0240']

#UNAM_statmap = {'C0200':'AM01', 'C0230':'AM03', 'C0120':'AM04', 'C01E0':'AM05', 
#                'C0180':'AM06', 'C01B0':'AM10', 'C0220':'AM12', 'C00D0':'AM13', 
#                'C00E0':'AM15', 'COOA0':'AM16', 'C0240':'AM17'}



multiprocess = False
global reftek
reftek = False

if multiprocess:
    import multiprocessing as mp


t_total0 = datetime.datetime.now()


    
folder_path = MSEED_DIR

# set file extension
#extensions = ['m', 'mseed', 'miniseed', 'MSEED']
extensions = ['msd']

 
#abs_paths = []
#for extension in extensions:
#    abs_paths.append(paths(folder_path, extension))

#print abs_paths



#flatten the list of absolute paths containing all relevant data files
#abs_paths = list(itertools.chain(*abs_paths))



    
abs_paths = get_filepaths(MSEED_DIR)


# create new path REMOVE LATER
abs_paths = [path for path in abs_paths if 'BHZ' in path]

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

def extract_info(info, check=False):    
    trace, path = info
    
    # check that the station names fit within the checklist, else map them.
    
    #if check:
    #    trace = check_stat(trace, checklist=UNAM_statlist, 
    #                           map_dict=UNAM_statmap)
    
    path_info = path.split('/')
    alt_station = path_info[-3]
    
    try:
        stats = trace.stats
        network = stats.network
        station = stats.station
        channel = stats.channel
        
        if station == '' and len(alt_station) == 4:
             station = statmap[path_info[-3]]
        
        if network == '':
            network = 'XX'
        
        if len(channel) == 1:
            channel = 'DH' + channel
            # potentially raw reftek format. 
            
        #    station = 
        code ='{}.{}.{}'.format(network, station, channel)

    except: 
        a=5


def extract_info(info):    
    trace, path = info
    
    stats = trace.stats
    code ='{}.{}.{}'.format(stats.network, stats.station, stats.channel)
    
    starttime = trace.stats.starttime.timestamp
    endtime = (trace.stats.starttime + trace.stats.npts * \
              (1/trace.stats.sampling_rate)).timestamp

    return (code, starttime, endtime, path)


    starttime = trace.stats.starttime.timestamp
    try:
        endtime = trace.stats.endtime.timestamp
    except:    
        endtime = (trace.stats.starttime + trace.stats.npts * \
        (1.0/trace.stats.sampling_rate)).timestamp

        information = (code, starttime, endtime, path)
    
        return information

        
def info_from_headers(path):


    print os.path.basename(path)
    try:
        #t0 = datetime.datetime.now()
        if reftek:
            headers = read_ref(path)
            print headers
        else:
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



timeline = np.array(timeline)
#print "timeline: ", timeline
#timeline = np.array([i for i in timeline if all(i)])
timeline = timeline[timeline != np.array(None)]

#flatten output list and  remove empty lists from timelines database
timeline = np.asarray(list(itertools.chain(*timeline)))
print "timeline: ", timeline

timeline = [tuple(i) for i in timeline]
print "timeline: ", timeline



#timeline = [tuple(i) for i in timeline]

    
#try:
    
    #flatten output list
#    timeline = np.asarray(list(itertools.chain(*timeline)))

#except:
#    print "The timeline array is not the correct type and cannot be flattened"

# remove None's from timeline array. 
timeline = [tuple(i) for i in timeline if i is not None]#timeline[timeline != np.array(None)]
print "timeline yeah: ", timeline


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
    

paths_ = c.execute('''SELECT DISTINCT file_path FROM timeline''')

# find out how many unique paths there are!
global paths_list
paths_list = []
def append_fast(list_):
    paths_list.append(str(list_[0]))    
map(append_fast, paths_)
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
print "file_extrema: ", file_extrema
c.executemany('INSERT INTO file_extrema VALUES (?,?,?,?)', file_extrema)

conn.commit()

conn.close()
