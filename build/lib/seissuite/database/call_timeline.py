# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 13:09:01 2015

@author: iese
"""

import sqlite3 as lite
from obspy.core import UTCDateTime as dt


database_name = '/home/iese/Documents/Ben/UNAM/INPUT/DATABASES/timeline.db'

print 'Initialising SQL database ... '
conn = lite.connect(database_name)
c = conn.cursor()
print 'SQL database cursor initialised.'


# find unique station names
unique_stations = list(c.execute('''SELECT DISTINCT station FROM timeline'''))

starttime = dt(list(c.execute('''SELECT MIN(starttime) 
                              FROM file_extrema'''))[0][0]).date

endtime = dt(list(c.execute('''SELECT MAX(endtime) 
                            FROM file_extrema'''))[0][0]).date



UNAM_statmap = {'C020':'AM01', 'C023':'AM03', 'C012':'AM04', 'C01E':'AM05', 
                'C018':'AM06', 'C01B':'AM10', 'C022':'AM12', 'C00D':'AM13', 
                'C00E':'AM15', 'COOA':'AM16', 'C024':'AM17'}

print unique_stations
quit()
