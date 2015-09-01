# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:12:37 2015

@author: boland
"""
from pysismo import pscrosscorr, pserrors, psstation
import os
import sys
sys.path.append('/home/boland/Anaconda/lib/python2.7/site-packages')
import warnings
import datetime as dt
import itertools as it
import pickle
import obspy.signal.cross_correlation
import numpy as np
import time
import glob
import profile
import matplotlib.pyplot as plt
from obspy.xseed import Parser
from mpl_toolkits.basemap import Basemap



#from pysismo.psconfig import (MSEED_DIR, 
#                              DATALESS_DIR, 
#                              STATIONXML_DIR,
#                              USE_DATALESSPAZ, 
#                              USE_STATIONXML)
                              
draw_paths = False


minlatitude = -45.0
minlongitude = 110.0
maxlatitude = 0.0
maxlongitude = 160.0

plt.figure(1)
colours = ["black", "blue", "green", "yellow", "purple", 
           "orange", "white", "red", "brown"]

dataless_files = glob.glob('/home/boland/Dropbox/University/UniMelb/AGOS/PROGRAMS/dataless/*.dataless')

for index, files in enumerate(dataless_files):
    
    network = os.path.basename(files).split('.')[0]
    sp = Parser(files)
        
    info = sp.getInventory()
        
    coordinates = [(i['longitude'], i['latitude'], i['channel_id']) 
    for i in info['channels'][:] if i['channel_id'][-3:] == "BHZ"]

                              
#dataless_inventories = []
#if USE_DATALESSPAZ:
#    with warnings.catch_warnings():
#        warnings.simplefilter('ignore')
#        dataless_inventories = psstation.get_dataless_inventories(DATALESS_DIR,
#                                                                  verbose=True)
        
#        info = dataless_inventories.getInventory()
        
#        coordinates = [(i['longitude'], i['latitude'], i['channel_id']) 
#        for i in info['channels'][:] if i['channel_id'][-3:] == "BHZ"]

        

#xml_inventories = []
#if USE_STATIONXML:
#    xml_inventories = psstation.get_stationxml_inventories(STATIONXML_DIR,
#                                                           verbose=True)

    longitudes = [float(i[0]) for i in coordinates]
    latitudes = [float(i[1]) for i in coordinates]

                                                    
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
    colour = colours[index]
    
    if draw_paths:
        for i in range(0, len(latitudes)):
        
            lat = latitudes[i]
            lon = longitudes[i]
    
            for j in range(0, len(latitudes)):
                lat1 = latitudes[j]
                lon1 = longitudes[j]
                m.drawgreatcircle(lon, lat, lon1, lat1, linewidth=0.03, 
                                  color=colour)
        
        
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
    m.plot(x, y, '^', label ='{}'.format(network),  markersize = 5, color = colour)

        
lg = plt.legend()
lg.get_frame().set_facecolor('grey')

plt.title("Plot of the Australian Based Seismic Networks.")

plt.show() 

