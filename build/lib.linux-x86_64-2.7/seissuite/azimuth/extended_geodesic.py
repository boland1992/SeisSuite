# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 08:44:50 2015

@author: boland
"""
#------------------------------------------------------------------------------
# MODULES
#------------------------------------------------------------------------------
import os
import fiona
import pickle
import pyproj
import datetime
import itertools
import shapefile
import numpy as np
import datetime as dt
import pointshape as ps
import multiprocessing as mp
import matplotlib.pyplot as plt
from math import sqrt, atan2, asin, degrees, radians, tan, sin, cos
from shapely.geometry import asPolygon, Polygon
from info_dataless import locs_from_dataless
from descartes.patch import PolygonPatch
from matplotlib.colors import LogNorm
from scipy.spatial import ConvexHull
from scipy.cluster.vq import kmeans
from shapely.affinity import scale
from matplotlib.path import Path
from shapely import geometry

#------------------------------------------------------------------------------
# VARIABLES
#------------------------------------------------------------------------------

verbose = False
#Enter path to boundary shape file.
shape_boundary = False
shape_path = "/home/boland/Dropbox/University/UniMelb\
/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"
dataless = True
# Enter path to dataless file
dataless_path = 'ALL_AUSTRALIA.870093.dataless'
#dataless_path = 'UOM.dataless'
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


# Enter path to dataless file
dataless_path = 'ALL_AUSTRALIA.870093.dataless'

coords = locs_from_dataless(dataless_path)

shape_path = "/home/boland/Dropbox/University/UniMelb\
/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"

# Generate InShape class
SHAPE = InShape(shape_path)

# Generate Shape lat-lon boundaries 

lonmin, latmin, lonmax, latmax = SHAPE.shape_bounds()

# Create shapely polygon from imported shapefile 
UNIQUE_SHAPE = SHAPE.shape_poly()
GEODESIC = Geodesic()
INPOLY = InPoly(shape_path)
#-----------------------------------------------------------------------------
# GENERATE SECOND SET OF VARIABLES AND STATES
#-----------------------------------------------------------------------------
ideal_path = 'ideal_coordinates.pickle'
#if no paths have been done before, start afresh!
if dataless:
    coords = locs_from_dataless(dataless_path)
    original_coords = coords
elif os.path.exists(ideal_path):
    f = open(name=ideal_path, mode='rb')
    coords = pickle.load(f)
    f.close()


coord_combs = [[coord1[0],  coord1[1],  coord2[0],  coord2[1]] \
               for coord1 in coords for coord2 in coords]

coord_combs = np.vstack(coord_combs)

print coord_combs



# testing canberra to perth 
#coords = np.asarray([[149.1244, -35.3075],[115.8589, -31.9522]])

def extended_paths(coord_list):
    return waypoint_init(coord_list)
    
t0 = datetime.datetime.now()
pool = mp.Pool()    
paths = pool.map(extended_paths, coord_combs)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print t1-t0

paths = np.asarray(list(itertools.chain(*paths)))

    
#keep all but the repeated coordinates by keeping only unique whole rows!
b = np.ascontiguousarray(paths).view(np.dtype\
((np.void, paths.dtype.itemsize * \
paths.shape[1])))
_, idx = np.unique(b, return_index=True)
        
paths = np.unique(b).view(paths.dtype)\
.reshape(-1, paths.shape[1])


paths = INPOLY.points_in(paths, IN=False)        

nbines = 300

H, xedges, yedges = np.histogram2d(paths[:,0],
                                   paths[:,1],
                                   bins=nbins)


H = np.rot90(H)
H = np.flipud(H)
H_masked = np.ma.masked_where(H==0,H)  


fig = plt.figure(figsize=(15,10))
plt.title('Equidistant Points Along Geodesic Path Between Stations \
Extended Behind First Station')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')


patch = PolygonPatch(UNIQUE_SHAPE, facecolor='white',\
                     edgecolor='k', zorder=1)
ax = fig.add_subplot(111)
ax.add_patch(patch)

#plt.scatter(paths[:,0], paths[:,1], s=15, c='orange',alpha=0.5, zorder=2)         
plt.pcolor(xedges, yedges, H_masked, norm=LogNorm(\
            vmin=np.min(H_masked), vmax=np.max(H_masked)), cmap='rainbow',\
            alpha=0.6, zorder = 3)
col = plt.colorbar()
col.ax.set_ylabel('Points Per Bin')  
ax.set_xlim(lonmin-0.15*abs(lonmax-lonmin), \
            lonmax+0.15*abs(lonmax-lonmin))
ax.set_ylim(latmin-0.15*abs(latmax-latmin), \
            latmax+0.15*abs(latmax-latmin))
            
            
fig.savefig('noise_sources.svg', format='SVG')
