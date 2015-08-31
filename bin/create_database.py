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

from pysismo.psconfig import MSEED_DIR

multiprocess = True

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

folder_path = MSEED_DIR




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
global timeline
timeline = {}

if not os.path.exists('tmp'): os.mkdir('tmp')


# =============================================================================
# SQL
# =============================================================================

# create database if it doesn't exist already, if it does, stop the programme.
database_name = 'timeline_database.db'

if os.path.exists(database_name):
    yeses = ['y','Y','yes','Yes','YES']    
    nos = ['n','N','no','No','NO']    
    
    condition = False
    while condition is False:
        
        answer = raw_input('Would you like to remove the existing database \
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

# =============================================================================
# =============================================================================

def extract_info(info):    
    trace, path = info
    
    stats = trace.stats
    code ='{}.{}.{}'.format(stats.network, stats.station, stats.channel)
    starttime = trace.stats.starttime.timestamp
    endtime = (trace.stats.starttime + trace.stats.npts * \
              (1/trace.stats.sampling_rate)).timestamp

    return (code, starttime, endtime, path)

def info_from_headers(path):
 
    #t0 = datetime.datetime.now()
    headers = read(path, headonly=True)
    headers.select(component='Z')
    info = []
    for trace in headers:
        info.append([trace, path])

    timeline_header = map(extract_info, info)
    
    return timeline_header
    #t1 = datetime.datetime.now()
    #print 'time taken to process previous loop: ', t1-t0
    #print "time for previous loop was: ", t1-t0
    


t0 = datetime.datetime.now()
pool = mp.Pool(None)
timeline = pool.map(info_from_headers, abs_paths)
pool.close()
pool.join()

#flatten output list
timeline = np.asarray(list(itertools.chain(*timeline)))
print timeline.shape

t1 = datetime.datetime.now()
print "time taken to read in and parallel process timeline database: ", t1-t0


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

t0 = datetime.datetime.now()
# commit many rows of information! 
c.executemany('INSERT INTO timeline VALUES (?,?,?,?)', timeline)



t1 = datetime.datetime.now()
conn.commit()

print "Time taken to commit thousands of rows to SQL db: ", t1-t0




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
   
    for info in row_extrema: file_extrema.append(info)
  


file_extrema = tuple(file_extrema)


c.executemany('INSERT INTO file_extrema VALUES (?,?,?,?)', file_extrema)

conn.commit()

conn.close()


















#def add_row_SQL(stat_info, conn=CONN, c=C):
    #""" 
    #This function takes the input of single row to be added to an SQL db
    #and inserts that row. 
    #"""    
 #   stat_name, start, end, file_path = stat_info
  #  example_row = '\'{}\',\'{}\',\'{}\',\'{}\''.format(stat_name, 
  #                                                     start,
  #                                                     end, 
  #                                                     file_path)
                                                                                    
    #print example_row
    
    # Insert a row of data        
   # c.execute("INSERT INTO timeline VALUES ({})".format(example_row))

    # Save (commit) the changes
    #conn.commit()