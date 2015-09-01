# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 08:31:35 2015

@author: boland
"""

from obspy import read
import os
import itertools
import datetime
import numpy as np
import multiprocessing as mp
# set folder with all raw waveform files in it. Can be recursive! 
folder_path = '/storage/ANT/INPUT/DATA/S-2014/2014-01'



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
    

# set file extension
extensions = ['m', 'mseed', 'miniseed', 'MSEED']
 
# set desired component e.g. E, N or Z

abs_paths = []
for extension in extensions:
    abs_paths.append(paths(folder_path, extension))

#flatten the list of absolute paths containing all relevant data files
abs_paths = list(itertools.chain(*abs_paths))



def extract_info(info):    
    trace, path = info

    code ='{}.{}'.format(trace.stats.network, trace.stats.station)
    print code
    starttime = trace.stats.starttime.timestamp
    endtime = (trace.stats.starttime + trace.stats.npts * \
              (1/trace.stats.sampling_rate)).timestamp
    #if code not in timeline.keys():
    #    timeline[code] = []
    #    timeline[code].append([starttime, endtime, path])
    #else:
    #    timeline[code].append([starttime, endtime, path])

    #print  [code, starttime, endtime, path]
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
#pool = mp.Pool(None)
timeline = map(info_from_headers, abs_paths)
#pool.close()
#pool.join()

#flatten output list
timeline = np.asarray(list(itertools.chain(*timeline)))
print timeline.shape

t1 = datetime.datetime.now()
print "time taken to read in and process the headers of {} files: ".format\
(len(abs_paths)), t1-t0



# todo: write low level function inspired of obspy.mseed.util._getRecordInformation
