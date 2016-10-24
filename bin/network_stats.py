# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 09:08:01 2016

@author: boland

CODE DESCRIPTION: 
The following python script takes a dataless SEED input file for an entire
seismic network, and pumps out desired statistics based on station spacing, 
such as standard deviations, means and any other desirable statistics that
can be thought of!
"""
#------------------------------------------------------------------------------
# MODULES
#------------------------------------------------------------------------------
from seissuite.spacing.search_station import (InShape, 
                                              InPoly, 
                                              Geodesic, 
                                              Coordinates, 
                                              Density)

import os
import pickle
import pyproj
import datetime
import numpy as np
import datetime as dt
#import pointshape as ps
import multiprocessing as mp
from scipy.cluster.vq import kmeans
from seissuite.misc.dataless import Dataless
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#------------------------------------------------------------------------------
# VARIABLES
#------------------------------------------------------------------------------
show = False

region_name = "Japan and The Korean Pennisula"

# Enter path to dataless file

#dataless_path = './networks/australia.dataless'
dataless_path = './networks/Japan.+.Korea.960525.dataless'


# Enter number new stations desired.
N = 3
# Enter km spacing between path density points.
km_points = 20.0
# Reference elipsoid to calculate distance.
wgs84 = pyproj.Geod(ellps='WGS84')
# Enter number of bins for 2D Histogram density calculation. 
nbins = 200
# Enter estimated average shear wave velocity. 3kms-1 is the default!
velocity = 3.0
# Define your ambient noise period range OR individual period in seconds.
global period_range
period_range = [1,40]



dataless_obj = Dataless(dataless_path=dataless_path)

coords = dataless_obj.locs_from_dataless()

lons, lats, elevs = coords[:,0], coords[:,1], coords[:,2]

if show:
    plt.figure()
    plt.scatter(lons, lats)
    plt.show()

#-----------------------------------------------------------------------------
# GENERATE SECOND SET OF VARIABLES AND STATES
#-----------------------------------------------------------------------------
lonmin, lonmax = np.floor(min(coords[:,0])), np.ceil(max(coords[:,0]))
latmin, latmax = np.floor(min(coords[:,1])), np.ceil(max(coords[:,1]))
print lonmin, lonmax, latmin, latmax

kappa = [np.vstack([[coord1[0],coord1[1],coord2[0],coord2[1]]\
                    for coord2 in coords]) for coord1 in coords]
                        
GEODESIC = Geodesic(km_point=km_points)

def spread_paths(coord_list):
    return GEODESIC.fast_paths(coord_list)
    
t0 = datetime.datetime.now()
pool = mp.Pool()    
paths = pool.map(spread_paths, kappa)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print t1-t0

paths = GEODESIC.combine_paths(paths)
paths = GEODESIC.remove_zeros(paths)

lons, lats = paths[:,0], paths[:,1]

if show:
    plt.figure()
    plt.scatter(lons, lats)
    plt.show()
    
DENSITY = Density(paths=paths)

H, xedges, yedges = DENSITY.hist2d(paths=paths)

H = np.rot90(H)
H = np.flipud(H)
H = np.ma.masked_where(H==0,H)  

H_avg = np.average(H)
H_std = np.std(H)

print "The point density distribution average for {} is: {} ".format(region_name, H_avg)
print "The point density distribution standard deviation for {} is: {} ".format(region_name, H_std)



latmin, latmax, lonmin, lonmax = np.min(lats), np.max(lats), np.min(lons), np.max(lons)
swell = 0.05
# Generate InShape class
# SHAPE = InShape(shape_path)
# Create shapely polygon from imported shapefile 
# UNIQUE_SHAPE = SHAPE.shape_poly()
fig = plt.figure(figsize=(15,10), dpi=100)
        
plt.xlabel('longitude (degrees)')
plt.ylabel('latitude (degrees)')
plt.xlim(lonmin-swell*abs(lonmax-lonmin),\
         lonmax+swell*abs(lonmax-lonmin))
plt.ylim(latmin-swell*abs(latmax-latmin),\
         latmax+swell*abs(latmax-latmin))
         
         
plt.title("Seismic Network Path Density Distribution of {}".format(region_name))

plt.pcolor(xedges, yedges, H, norm=LogNorm(\
            vmin=np.min(H), vmax=np.max(H)), cmap='rainbow',\
            alpha=0.6)
col = plt.colorbar()
col.ax.set_ylabel('Points Per Bin')
fig.savefig("Seismic Network Path Density Distribution of {}.png".format(region_name), dpi=300)

