# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 11:27:31 2015

@author: boland

CODE DESCRIPTION:
The following python script is used to find the lon-lat locations of 
seismic stations from a given dataless SEED metadata file and then plot 
them. 
"""

from mpl_toolkits.basemap import Basemap
from info_dataless import locs_from_dataless

import matplotlib.pyplot as plt
import numpy as np


dataless_path = '/home/boland/Dropbox/University/UniMelb/AGOS/METADATA/metadata/UOM.dataless'

info = locs_from_dataless(dataless_path)

lats = info[:,1].astype(np.float); lons = info[:,2].astype(np.float)

print(lats); print(lons)

minlatitude=np.min(lats) - 0.5
minlongitude =np.min(lons) - 0.5
maxlatitude=np.max(lats) + 0.5
maxlongitude=np.max(lons) + 0.5


plt.figure()

plt.subplots_adjust(left=0.05,right=0.95,top=0.90,bottom=0.05,wspace=0.15,hspace=0.05)
ax = plt.subplot(111)
#draw blank map

m = Basemap(resolution='i',projection='merc',\
llcrnrlat=minlatitude,urcrnrlat=maxlatitude,llcrnrlon=minlongitude,\
urcrnrlon=maxlongitude,lat_ts=45)
    
    
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

x, y = m(lons, lats)


m.plot(x, y, '^', label ="station",  markersize = 5, color = 'k')

plt.title("Plot of station pairs for all UM Gippsland seismometers")
        
lg = plt.legend()
lg.get_frame().set_facecolor('grey')

plt.show() 