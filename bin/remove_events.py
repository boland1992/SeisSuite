# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 16:10:23 2016

@author: boland
"""

database = '/storage/MASTERS/CONFIGURATIONS/S_NETWORK/INPUT/DATABASES/timeline.db'

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 16:54:23 2016
@author: iese
"""

from obspy.core import UTCDateTime as utc
import sqlite3 as lite
import os
import numpy as np
from obspy import read
from obspy.core import Stream

# import configuration file for DAS to station map



#noise_start = utc('2015-10-11T23:46:32')
#noise_end = utc('2015-10-11T23:47:37')
#signal_start =  utc('2015-10-11T23:47:37')
#signal_end =  utc('2015-10-11T23:50:02')

show_ = True # show the events as they are found in the raw data?

# create database if it doesn't exist already, if it does, stop the programme.
database_loc = '/storage/MASTERS/CONFIGURATIONS/S_NETWORK/INPUT/DATABASES/timeline.db'
#event_cat = os.path.join(folder_path, 'catalogue.csv')



events_array = []

if os.path.exists(database_loc):
    conn = lite.connect(database_loc)
    c = conn.cursor()

    try:
       	c.execute('''CREATE TABLE catalogue (id text, lon real, lat real, timestamp real, depth real, mag real)''')

    except Exception as error:
	print error

    existing_list = c.execute('''SELECT * FROM 
               		              catalogue''')


    for i in existing_list: 
	events_array.append(i[1:])

    conn.commit()
    conn.close()


events_array = np.array(events_array)
#print events_array 




print "events array length: ", len(events_array)

# read in event catalogue csv
#event_cat = 'catalogue.csv'



if not os.path.exists(database_loc):
    raise Exception("There is no timeline database to construct events.")

# state number of seconds before an event for time series window
event_minus = 0

# state number of seconds after an event for time series window
event_plus1 = 0
event_plus2 = 300

# remove header
lons, lats, timestamps, depths, mags = events_array[:,0], events_array[:,1], \
events_array[:,2], events_array[:,3], events_array[:,4]

#print lons, lats, timestamps, depths, mags



#SAM_statmap = {'C0200':'AM01', 'C0230':'AM03', 'C0120':'AM04', 'C01E0':'AM05', 
#                'C0180':'AM06', 'C01B0':'AM10', 'C0220':'AM12', 'C00D0':'AM13', 
#                'C00E0':'AM15', 'COOA0':'AM16', 'C0240':'AM17'}


#SBB_statmap = {'C00A0':'B01', 'C00D0':'B02', 'C01B0':'B03', 'C01F0':'B04', 
#                'C0130':'B05', 'C0240':'B06'}



def getpaths(database_loc, starttime, endtime):
        """
        Gets list of paths to mseed files between certain time series intervals
        using initialised SQL timeline database.
        @type starttime: L{UTCDateTime} or L{datetime} or L{date}
        @type endtime: L{UTCDateTime} or L{datetime} or L{date}
        @rtype: unicode
        """

        starttime, endtime = utc(starttime), utc(endtime)
        
        import_start = starttime.timestamp
        #import_end = endtime.timestamp
        
        #connect SQL database

        if not os.path.exists(database_loc):
            raise Exception("Database doesn't exist")
        
        conn = lite.connect(database_loc)
        #print "conn: ", conn
        
        c = conn.cursor()
        extrema = []
        for row in c.execute('''SELECT * FROM 
                             file_extrema ORDER BY station'''):
            extrema.append(row)
        
        # make sure that this test works! 
        #code = '{}.{}.{}'.format(self.network, self.name, self.channel)
    
        file_paths = c.execute('''SELECT file_path FROM 
                             file_extrema WHERE starttime <= ? 
                             AND endtime >= ?''', 
                             (import_start, import_start))        
        
        output_path = []
        for file_path in file_paths: output_path.append(file_path)
        #if len(output_path) > 0:
        #    return str(output_path[0][0])

        
        # close database
        conn.close() 
        return output_path

counter = 1
print "There are {} events in the catalogue to search.".format(len(timestamps))

# list stating whether or not the network was operational during an event!
operational = []

# search through every event's timestamp
for timestamp in timestamps:
    print "processing information for event: ", counter
    print '\n'


    print "the magnitude of the event displayed is: {}".format(mags[counter-1])
    t_start = timestamp + event_plus1
    t_end = timestamp + event_plus2
    
    
    # get only the miniseed paths that have data between t_start and t_end
    event_paths = getpaths(database_loc, t_start, t_end)

    event_traces = []
    SNRs = []

    # multiplex the traces and only show Z components
    for event_path in event_paths: 
        event_path = str(event_path[0])
        
        try:
            print "processing event path: {}".format(event_path)     
            print utc(t_start), utc(t_end)
            st = read(event_path)
            
            end = utc(t_start)+3600

            st = st.cutout(starttime=utc(t_start), endtime=utc(t_end))
            st.merge()
            for tr in st:
                tr.data = np.ma.filled(tr.data)
            st.write(event_path, format='MSEED')
            
        except Exception as error:
            print error

    counter += 1            