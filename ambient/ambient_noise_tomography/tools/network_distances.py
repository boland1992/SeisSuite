# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 18:56:35 2015

@author: boland
"""


from __future__ import with_statement
import sys
sys.path.append('/home/boland/Anaconda/lib/python2.7/site-packages')
import math, os, io, glob, shutil, time, obspy, \
       datetime, calendar, random, urllib2, \
       StringIO, csv, numpy as np, \
       matplotlib.pyplot as plt, time,\
       matplotlib.lines as mlines

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.fdsn import Client

from obspy import read
from os import open
import string



network_name = "S"

client = Client()
#get all the station data from a network of a specific network_name
#info = client.get_stations(network = network_name)
#get stations based on a set radius (in degrees) from a centralised point

#info = client.get_stations(latitude=-27.5, longitude=132.5, maxradius=35)




#get stations based on a box of lats and longs

#latitude and longitude definitions for Australia
minlatitude = -10.0
maxlatitude = -45.0
minlongitude = 110.0
maxlongitude = 155.0


#set networks of interest
#network_names = ["AU", "S", "IU", "IM", "II", "G"]

network_names = ["AU","S", "IU", "GE"]#, "XD", "ZN", ]
network_names = []

alphabet = string.ascii_uppercase

for i in range(0,len(alphabet)):
    for j in range(0,len(alphabet)):
        net_string = "{}{}".format(alphabet[i],alphabet[j])
        network_names.append(net_string)
        
print(network_names)

counts2 = 0


colours = [ "black", "blue", "green", "yellow", "purple", "orange", "white", "red", "brown"]
#cols = map(lambda i: colour[i], nums)

minlatitude=-45.0
minlongitude =100.0
maxlatitude=10.0
maxlongitude=160.0

minlatitude=-18.0
minlongitude =142.5
maxlatitude=-1.0
maxlongitude=160.0

minlatitude=-89.0
minlongitude =-179.5
maxlatitude=89.0
maxlongitude=179.0


plt.figure(1)
stat_names = []; channels = []; latitudes = []; longitudes = []; 
end_dates = []; creation_dates = []; termination_dates = []; start_dates = [];

for net in network_names:
    #info = client.get_stations(network = "%s"%(net), minlatitude=-45.0, minlongitude =110.0, maxlatitude=-10.0, maxlongitude=155.0)
    try:
        #set lat long limits below!

        info = client.get_stations(network = "%s"%(net),\
        minlatitude=minlatitude, minlongitude =minlongitude, maxlatitude=maxlatitude, \
        maxlongitude=maxlongitude)

###attributes of get_stations()

#client.get_stations(starttime=None, endtime=None, startbefore=None,
#startafter=None, endbefore=None, endafter=None,
#network=None, station=None, location=None, channel=None,
#minlatitude=None, maxlatitude=None, minlongitude=None,
#maxlongitude=None, latitude=None, longitude=None,
#minradius=None, maxradius=None, level=None,
#includerestricted=None, includeavailability=None,
#updatedafter=None, matchtimeseries=None, filename=None)

        info = info[0][:]

#get lists of information about this network!


        for i in info:
            stat_names.append(i.code)
            channels.append(i.channels)
            latitudes.append(i.latitude)
            longitudes.append(i.longitude)
            creation_dates.append(i.creation_date)
            end_dates.append(i.end_date)
            termination_dates.append(i.termination_date)
            start_dates.append(i.start_date)
    
    
        print(stat_names) 
        print(len(stat_names))
#for i in channels: print(i)   
#for i in latitudes: print(i)   
#for i in longitudes: print(i)
#for i in creation_dates: print(i)
#for i in end_dates: print(i)
#for i in start_dates: print(i)
#for i in termination_dates: print(i)

#http://matplotlib.org/basemap/users/examples.html




#network, stations, latitude, longitude, data = stations()

#LAT = filter(None, latitude) #remove all null strings from list
#LONG = filter(None, longitude) #remove all null strings from list
#station_names = filter(None, stations) #remove all null strings from list

    

#latitudes = [float(i) for i in LAT]   
#longitudes = [float(i) for i in LONG]   

#print(latitudes); print(longitudes)

#print(len(latitudes))

#LAT = filter(None, latitude) #remove all null strings from list
#LONG = filter(None, longitude) #remove all null strings from list
#station_names = filter(None, stations) #remove all null strings from list

    

#latitude = [float(i) for i in LAT]   
#longitude = [float(i) for i in LONG]   



#data_folder_string, dist  = distance(latitude, longitude, station_names)

#counts = 0 
#distances = []

#for i in range(0, len(data_folder_string)):
    
#    print("%s : %s" %(data_folder_string[counts], dist[counts]))
#    distances.append(dist[counts])



# make sure the value of resolution is a lowercase L,
#  for 'low', not a numeral 18


         
    except(obspy.fdsn.header.FDSNException):
        k=2



#plt.figure()

#Custom adjust of the subplots
plt.subplots_adjust(left=0.05,right=0.95,top=0.90,bottom=0.05,wspace=0.15,hspace=0.05)
ax = plt.subplot(111)
#draw blank map

m = Basemap(resolution='i',projection='merc',\
llcrnrlat=minlatitude,urcrnrlat=maxlatitude,llcrnrlon=minlongitude,\
urcrnrlon=maxlongitude,lat_ts=45)
    
    
    
        #mid_lat = (-45.0 + -10.0)/2; print(mid_lat)
        #mid_long = (155.0 + 110.0)/2; print(mid_long)
    
        #m.bluemarble()
        #m.etopo()
        #m.shadedrelief()
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)    
    
m.fillcontinents(color='#cc9955',lake_color='aqua')
m.drawmapboundary(fill_color='aqua')
    #m.drawmapboundary() # draw a line around the map region

    #below gives grids
m.drawparallels(np.arange(-90.,120.,5.),labels=[1,0,0,0]) # draw parallels
m.drawmeridians(np.arange(60.,190.,5.),labels=[0,0,0,1]) # draw meridians 
    
m.drawmapscale(120, -42.5, 50, 180, 1000, barstyle='fancy')

x, y = m(longitudes, latitudes)

    #m.scatter(x, y, c="red"2) # we will scale the dots by 10 time the magnitude
        
        #generate a random colour
colour = "yellow"

counts = 0
counts1 = 0
for i in range(0, len(latitudes)):
        
    lat = latitudes[i]
    lon = longitudes[i]
    
    for j in range(0, len(latitudes)):
        lat1 = latitudes[j]
        lon1 = longitudes[j]
        m.drawgreatcircle(lon, lat, lon1, lat1, linewidth=0.03, color="black")
        #counts += 1; print(counts)
        

    #lines between points
    #m.plot(x, y, '^-', markersize=10)
               
        #sm.plot(x,y,'-', label=cat, color=color)
        
       
        
        
        ###legend attempt 1 FAILURE
        #plt.legend(handles=[blue_line])
        #blue_line = mlines.Line2D([], [], color='blue', marker='*',
        #                  markersize=15, label='%s'%(net))
        
        #colors = {'CCM': 'red', 'SCME': 'white', 'SCM': 'yellow'}


            #Plot the points on the map
    #m.plot(x, y, '^', label ="%s"%(net),  markersize = 10, color = colour)
m.plot(x, y, '^', label ="station",  markersize = 5, color = colour)

plt.title("Plot of station pairs for all PNG seismometers minus temporary networks")
counts2 += 1
        
lg = plt.legend()
lg.get_frame().set_facecolor('grey')
        #plt.hold(False)

        


plt.show() 

 
#try and add network names! and or station names to these plots!
#plt.legend(proxy, [])


#rand_colours = [random.choice(colour) for i in range(50)]





