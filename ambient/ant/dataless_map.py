# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:12:37 2015

@author: boland
"""
import os
import numpy as np
import glob
import matplotlib.pyplot as plt
from obspy.xseed import Parser
from mpl_toolkits.basemap import Basemap
                           
draw_paths = False
minlatitude = -45.0
minlongitude = 110.0
maxlatitude = 0.0
maxlongitude = 160.0

plt.figure(1)
colours = ["black", "blue", "green", "yellow", "purple", 
           "orange", "white", "red", "brown"]

dataless_files = glob.glob('/home/boland/Dropbox/University/UniMelb\
                            /AGOS/PROGRAMS/dataless/*.dataless')

for index, files in enumerate(dataless_files):
    network = os.path.basename(files).split('.')[0]
    sp = Parser(files)
    info = sp.getInventory()
    coordinates = [(i['longitude'], i['latitude'], i['channel_id']) 
    for i in info['channels'][:] if i['channel_id'][-3:] == "BHZ"]
    longitudes = [float(i[0]) for i in coordinates]
    latitudes = [float(i[1]) for i in coordinates]                                                
    m = Basemap(resolution='i',projection='merc',\
    llcrnrlat=minlatitude,urcrnrlat=maxlatitude,llcrnrlon=minlongitude,\
    urcrnrlon=maxlongitude,lat_ts=45)
    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)    
    m.fillcontinents(color='#cc9955',lake_color='aqua')
    m.drawmapboundary(fill_color='aqua')
    m.drawparallels(np.arange(-90.,120.,5.),labels=[1,0,0,0]) # draw parallels
    m.drawmeridians(np.arange(60.,190.,5.),labels=[0,0,0,1]) # draw meridians 
    m.drawmapscale(120, -42.5, 50, 180, 1000, barstyle='fancy')
    x, y = m(longitudes, latitudes)
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
    m.plot(x, y, '^', 
           label ='{}'.format(network), 
           markersize=5, color = colour)

lg = plt.legend()
lg.get_frame().set_facecolor('grey')
plt.title("Plot of the Australian Based Seismic Networks.")
plt.show() 

