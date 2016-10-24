# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:09:17 2016
@author: Benjamin Boland
"""

from obspy.fdsn import Client
from obspy.core import UTCDateTime as utc
from obspy.core.event import Catalog
import sqlite3 as lite
import numpy as np
import os


client_list = ['IRIS', 'USGS', 'USP', 'NIEP', 'ETH', 'GFZ', 'INGV', 'ORFEUS']

# Operation dates for IESE network
t1 = utc("2014-01-01T00:00:00")
t2 = utc("2015-01-01T00:00:00")



reset_catalogue = True
timelimit = utc("2005-01-01T00:00:00").timestamp
# set time limit and all events before this limit will be erased
# if reset_catalogue = True

database_name = '/storage/MASTERS/CONFIGURATIONS/S_NETWORK/INPUT/DATABASES/timeline.db'

# enter the information for catalogue 1 search
minlat, maxlat, minlon, maxlon = (-40.0, -12.5, 113.0, 154.0) 

event_list = []

for c in client_list:
    print "Processing events from the {} catalogue ... \n".format(c)


    try:
        client = Client(c)

        catalogue = client.get_events(starttime=t1, endtime=t2,
                                           minlatitude=minlat,
                                           maxlatitude=maxlat,
                                           minlongitude=minlon, 
                                           maxlongitude=maxlon)

            
        for i in catalogue:
            event_list.append(i)
            
    except Exception as error: 
        print error
        

print event_list
 
 

if len(event_list) > 0:
    final_catalogue = Catalog(events=event_list)

    print final_catalogue
    

lons = []
lats = []
times = []
depths = []
types = []
ids = []
mags = []
dists = []

existing_ids = []

if os.path.exists(database_name):
    conn = lite.connect(database_name)
    c = conn.cursor()
    
    try:
       	c.execute('''CREATE TABLE catalogue (id text, lon real, lat real, timestamp real, depth real, mag real)''')

    except Exception as error:
	print error



    existing_events = c.execute('''SELECT * FROM 
               		              catalogue''')


    for i in existing_events: 
	existing_ids.append(str(i[0]))

    conn.commit()
    conn.close()


for event in final_catalogue:

    resource_id = event.resource_id.id



    if resource_id not in existing_ids:
    
    	try:
    	    mags.append(event.magnitudes[0].mag)
    	    ids.append(resource_id)
            lons.append(event.origins[0].longitude)
    	    lats.append(event.origins[0].latitude)
    	    times.append(event.origins[0].time.timestamp)
    	    depths.append(event.origins[0].depth)
    	    types.append(event.event_type)

    #	dists.append(dist(event.origins[0].longitude, event.origins[0].latitude,
    #    	              -101.140612, 19.922488,))
    
    
    #dists.append(dist(event.origins[0].longitude, event.origins[0].latitude,
    #                  -99.5652864, 20.53126886))

        except Exception as error:
            print error


event_info = np.column_stack((ids, lons, lats, times, depths, mags))



#event_info = event_info[event_info[:,4].argsort()]
#print event_info

#if reset_catalogue:

#csv_info = np.column_stack((lons, lats, times, depths, mags))
#print "Saving csv file to : {}".format(os.path.join(folder_path, "catalogue.csv")) 
#np.savetxt(os.path.join(folder_path, "catalogue.csv"), csv_info, delimiter=",")





event_info = [tuple(i) for i in event_info]

print 'Creating SQL database: {} ...'.format(database_name)
conn = lite.connect(database_name)
c = conn.cursor()

print 'SQL database cursor initialised'

try:
    # Create table called catalogue for earthquake times, locations and so on. 

    c.execute('''CREATE TABLE catalogue (id text, lon real, lat real, timestamp real, depth real, mag real)''')
    print 'Catalogue table created'

except Exception as error:
    print error


try:
    c.executemany('INSERT INTO catalogue VALUES (?,?,?,?,?,?)', event_info)
    print '{} new events added to database'.format(len(event_info))


    
except Exception as error:
    print error


conn.commit()
conn.close()



conn = lite.connect(database_name)
c = conn.cursor()

existing_events = c.execute('''SELECT * FROM 
               		              catalogue''')

existing_list = []



#for i in existing_events: existing_list.append(i)

#print

conn.commit()
conn.close()
quit()

