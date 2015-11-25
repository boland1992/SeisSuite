# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 13:09:01 2015

@author: iese
"""

import sqlite3 as lite
from obspy.core import UTCDateTime as dt
from obspy import read
import itertools as it
import numpy as np
import os, sys
import shutil


database_name = '/home/iese/Documents/Ben/UNAM/INPUT/DATABASES/timeline.db'

print 'Initialising SQL database ... '
conn = lite.connect(database_name)
c = conn.cursor()
print 'SQL database cursor initialised.'


# find unique station names
unique_stations = list(c.execute('''SELECT DISTINCT station FROM timeline'''))

starttime = list(c.execute('''SELECT starttime
                              FROM file_extrema'''))

endtime = list(c.execute('''SELECT endtime 
                            FROM file_extrema'''))

starttime = [dt(i[0]) for i in starttime]
endtime = [dt(i[0]) for i in endtime]



UNAM_statmap = {'C020':'AM01', 'C023':'AM03', 'C012':'AM04', 'C01E':'AM05', 
                'C018':'AM06', 'C01B':'AM10', 'C022':'AM12', 'C00D':'AM13', 
                'C00E':'AM15', 'COOA':'AM16', 'C024':'AM17'}


paths = list(c.execute('''SELECT file_path 
                            FROM file_extrema'''))
                            
paths = np.array(list(it.chain(*paths)))

dest_folder = '/home/iese/Documents/Ben/UNAM2'

if not os.path.exists(dest_folder): os.mkdir(dest_folder)
    
    
for path in paths:
    # copy each file that is known to contain information, to another directory
    #print "Copying {} ... ".format(os.path.basename(path))
    #shutil.copy(path, dest_folder)
    st = read(path)
    st.plot()
    
quit()
