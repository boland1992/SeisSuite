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
import multiprocessing as mp
#try:
import cPickle 
#except:
import pickle#try using json vs. pickle and see how long this takes!
   #print "Can't import cpickle for some reason"

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

folder_path = 'small_borehole_quakes'
#folder_path = '/storage/ABE/borehole_data'
#folder_path = '/storage/ANT/INPUT/DATA/AGOS-FULL-DATA_LOWSAMPLE'
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


print abs_paths
try: 
    for path in abs_paths:
        t0 = datetime.datetime.now()
        headers = read(path, headonly=True)
        t1 = datetime.datetime.now()
        print 'time to import headers: ', t1-t0
        t0 = datetime.datetime.now()
        headers.select(component='Z')
        t1 = datetime.datetime.now()
        print 'time to select Z component: ', t1-t0
        path = str(path)
        for trace in headers:
            code ='{}.{}'.format(trace.stats.network, trace.stats.station)
            starttime = trace.stats.starttime.timestamp
            endtime = (trace.stats.starttime + trace.stats.npts * \
                (1/trace.stats.sampling_rate)).timestamp
            
            if code not in timeline.keys():
                timeline[code] = []
                timeline[code].append([starttime, endtime, path])
            else:
                timeline[code].append([starttime, endtime, path])
        #print timeline.keys()
        
        
except Exception as error:
    print error
#pool = mp.Pool(4)
#info = pool.map(scan_path, abs_paths)
#pool.close()
#pool.join()


for key in timeline.keys():
    timeline[key] = list(np.sort(np.asarray(timeline[key]), axis=0))


#==============================================================================
# PICKLE 
#==============================================================================

outfile = 'tmp/timeline_database.pickle'

t0 = datetime.datetime.now()
if not os.path.exists(outfile):
    with open(outfile, 'wb') as f:
        print "\nExporting new timeline database to: " + f.name
        pickle.dump(timeline, f, protocol=2)
else:
    raise Exception("Could not create new pickle database as one already exists")

t1 = datetime.datetime.now()
print "Time taken to dump database to pickle: ", t1-t0


#==============================================================================
# CPICKLE 
#==============================================================================
outfile_cpickle = 'tmp/Ctimeline_database.pickle'


t0 = datetime.datetime.now()
if not os.path.exists(outfile_cpickle):
    with open(outfile_cpickle, 'wb') as f:
        print "\nExporting new timeline database to: " + f.name
        cPickle.dump(timeline, f, protocol=2)
else:
    raise Exception("Could not create new pickle database as one already exists")

t1 = datetime.datetime.now()
print "Time taken to dump database to cPickle: ", t1-t0






#import json
# write to JSON file for increase speed, security and smaller file size!
# write to utf-8 format, rather than ascii for smaller file size still
  
#outfile_json = 'tmp/timeline_database.json'
  
  
#t0 = datetime.datetime.now()
#if not os.path.exists(outfile_json):
#    with open(outfile_json, 'wb') as fp:
#        json.dump(timeline, fp)
        
        
#else:
#    raise Exception("Could not create new JSON database as one already exists")

#t1 = datetime.datetime.now()
#print "Time taken to dump database to JSON: ", t1-t0
  
  


t_total1 = datetime.datetime.now()
print 'Total time taken to initialise timeline database for {} was: {}'.format(
            os.path.basename(folder_path), t_total1-t_total0)