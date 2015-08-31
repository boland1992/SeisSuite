# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:34:52 2015

@author: boland
"""

import sqlite3 as lite
from obspy.core import UTCDateTime
import datetime as dt
import sqlite3 as lite
import sys

        
database_name = 'timeline_database.db'

conn = lite.connect(database_name)


c = conn.cursor()

# Create table called timeline for timeline database
# the table will contain the station ID which is net.stat_name
# the start time and end time for each trace and the absolute file path loc
# the goal is to have one line per trace!

try:
    c.execute('''CREATE TABLE timeline
                 (stat name, start time, end time, file path)''')

except Exception as error:
    print error             
             

#con = lite.connect('post.db')
#postid = '52642'

#cur = con.cursor()
#cur.execute("select id from POSTS where id=?", (postid,))
#data = cur.fetchall()
#if data is None:
#        print ('not found')
#else:
#       print ('found')
                 
# create an example row of data

try:
    stat_name = 'UM.HOLS'             
    start = UTCDateTime(dt.datetime.now()).timestamp
    end = UTCDateTime(dt.datetime.now()).timestamp
    file_path = '/storage/UM.HOLS.MSEED'

    example_row = '\'{}\',\'{}\',\'{}\',\'{}\''.format(stat_name, start,
                                                       end, file_path)


    print example_row
    
    # Insert a row of data
    c.execute("INSERT INTO timeline VALUES ({})".format(example_row))

    # Save (commit) the changes
    conn.commit()
except Exception as error: 
    print error
    
    
for row in c.execute('SELECT * FROM timeline ORDER BY station'):
    print row

