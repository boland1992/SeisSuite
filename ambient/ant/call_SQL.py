# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 09:04:28 2015

@author: boland
"""

import sqlite3 as lite
import datetime
from obspy.core import UTCDateTime

verbose = True
database_name = 'timeline_database.db'


if verbose:    
    print 'linking SQL database ... '

conn = lite.connect(database_name)

if verbose:
    print 'SQL database linked correctly'
    print 'creating cursor object'
    
c = conn.cursor()

if verbose:
    print 'cursor object created successfully'
    print 'executing commands'
    

# Call for the first and last time values to run through. 

t0 = datetime.datetime.now()

#for row in c.execute('SELECT * FROM timeline ORDER BY station'):
#    print row

min_start = c.execute('''SELECT MIN(starttime) FROM timeline ORDER 
                      BY starttime DESC LIMIT 1''')

#unpack tuple
(min_start, ) = min_start
#convert to float, then to obspy UTCDateTime object
min_start = UTCDateTime(min_start[0])

max_end = c.execute('''SELECT MAX(endtime) FROM timeline ORDER 
                    BY endtime DESC LIMIT 1''')

#unpack tuple
(max_end, ) = max_end
#convert to float, then to obspy UTCDateTime object
max_end = UTCDateTime(max_end[0])

t1 = datetime.datetime.now()

if verbose:
    print 'Time taken to execute SELECT command was: ', t1-t0
    
    
print min_start
print max_end

extrema = []
for row in c.execute('''SELECT * FROM 
                     file_extrema ORDER BY station'''):

    extrema.append(row)
    #starttime, endtime = UTCDateTime(row[0]), UTCDateTime(row[1])  
    
# We can also close the connection if we are done with it.
# Just be sure any changes have been committed or they will be lost.
extrema = tuple(extrema)
print extrema
conn.close()