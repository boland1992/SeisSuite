# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 08:53:56 2015

@author: boland
"""

from obspy.fdsn.header import URL_MAPPINGS
from obspy import UTCDateTime
from obspy.fdsn import Client

# create example start and end times for event search
starttime = UTCDateTime('2014-01-01T00:00.000')
endtime = UTCDateTime('2015-01-01T00:00.000')

endtime = UTCDateTime('2014-02-01T00:00.000')

# create list of possible servers to find earthquake events
server_list = []
for key in sorted(URL_MAPPINGS.keys()):
    server_list.append(key)

for server in server_list:
    print server
    client = Client(server)
    try:
        cat = client.get_events(starttime=starttime, endtime=endtime, 
                                minmagnitude=4)#, catalog="ISC")
        print cat
        cat.plot()
    except:
        continue
    
print "done"