# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 12:28:32 2015

@author: boland
"""

import sys
sys.path.append('/home/boland/Anaconda/lib/python2.7/site-packages')
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans
import multiprocessing as mp
import pyproj
import os
import itertools
import datetime
import pointshape as ps
from math import sqrt, atan2, radians,degrees, cos, tan, sin, asin



shape_path = "/home/boland/Dropbox/University/UniMelb/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"

N = 130
#enter km spacing between path density points
km_points = 20.0

# reference elipsoid to calculate distance
wgs84 = pyproj.Geod(ellps='WGS84')
nbins = 200

def haversine(coordinates):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    lon1, lat1, lon2, lat2= coordinates[0],coordinates[1],\
    coordinates[2],coordinates[3]

    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
    # haversine formula 
    dlon, dlat = lon2 - lon1, lat2 - lat1
    
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    
    return km
    
def haversine2(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
    # haversine formula 
    dlon, dlat = lon2 - lon1, lat2 - lat1
    
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    
    return km

def geodesic(coord1, coord2, npts):
    """
    Returns a list of *npts* points along the geodesic between
    (and including) *coord1* and *coord2*, in an array of
    shape (*npts*, 2).
    @rtype: L{ndarray}
    """
    if npts < 2:
        raise Exception('nb of points must be at least 2')

    path = wgs84.npts(lon1=coord1[0], lat1=coord1[1],
                      lon2=coord2[0], lat2=coord2[1],
                      npts=npts-2)
    return np.array([coord1] + path + [coord2])
    
def new_geodesic(lon1,lat1,lon2,lat2, npts):
    """
    Returns a list of *npts* points along the geodesic between
    (and including) *coord1* and *coord2*, in an array of
    shape (*npts*, 2).
    @rtype: L{ndarray}
    """
    if npts < 2:
        raise Exception('nb of points must be at least 2')

    path = wgs84.npts(lon1=lon1, lat1=lat1,
                      lon2=lon2, lat2=lat2,
                      npts=npts-2)
                      
    return np.array([[lon1,lat1]] + path + [[lon2,lat2]])

def cluster_points(coord_points):
    """
    Function that returns k which is an nx2 matrix of lon-lat vector columns
    containing the optimal cluster centroid spacings within a large set of random
    numbers e.g. those produced by the many_points() function above!
    """

    k = kmeans(coord_points, 130)
    
    return k[0]

def new_paths(path_info, km=km_points):

    lon1, lat1, lon2, lat2 = path_info[0], \
    path_info[1], path_info[2], path_info[3]
    
    #lon1, lat1, lon2, lat2, dist = path_info[0], \
    #path_info[1], path_info[2], path_info[3],    \
    #path_info[4]  
  
    dist = haversine2(lon1, lat1, lon2, lat2)
    # interpoint distance <= 1 km, and nb of points >= 100
    npts = max(int((np.ceil(dist) + 1)/km), 100)

    path = new_geodesic(lon1,lat1,lon2,lat2, npts)
    #print("still going strong\n")
    return path


def HIST2D(nbins,paths, grad=False):

    H, xedges, yedges = np.histogram2d(paths[:,0],paths[:,1],bins=nbins)
    
    #name = "path_density_2Dhist.png"
    

    if grad:
        H = np.abs(np.asarray(np.gradient(H)[0]))#name = "path_density_2Dhist_grad.png"

    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
 
     # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero  
    
    return Hmasked
    #fig = plt.figure()
    #plt.pcolormesh(xedges,yedges,Hmasked)
    #plt.xlabel('longitude (degrees)')
    #plt.ylabel('longitude (degrees)')
    #cbar = plt.colorbar()
    #cbar.ax.set_ylabel('Counts')
    
    #fig.savefig(name)
    

def latitude(dist, sigma01, alpha0, lon0):
    sigma = sigma01 + dist#/R
    lat = degrees(asin(cos(alpha0)*sin(sigma)))
    #alpha = atan2(tan(alpha0),cos(sigma))
    return lat
    
def longitude(dist, sigma01, alpha0, lon0):
    sigma = sigma01 + dist#/R
    lon = degrees(atan2(sin(alpha0)*sin(sigma), cos(sigma))) + degrees(lon0)
    #alpha = atan2(tan(alpha0),cos(sigma))
    return lon

vlat_func = np.vectorize(latitude)
vlon_func = np.vectorize(longitude)

def waypoint_init(path_info, km=km_points):
    R = 6371
    lon1, lat1, lon2, lat2, dist = radians(path_info[0]), \
    radians(path_info[1]), radians(path_info[2]), \
    radians(path_info[3]), radians(path_info[4])
    #lon1, lat1, lon2, lat2, dist = map(radians, [path_info[0],path_info[1],path_info[2],path_info[3],path_info[4]])
    lon_diff = lon2-lon1
    alpha1 = atan2(sin(lon_diff),(cos(lat1)*tan(lat2)-sin(lat1)*cos(lon_diff)))
       #alpha2 = atan2(sin(lon_diff),(-cos(lat2)*tan(lat1)+sin(lat2)*cos(lon_diff)))
    #try:
        #sigma12 = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon_diff))
    #except:
        #return
    sigma01, alpha0 = atan2(tan(lat1), cos(alpha1)), asin(sin(alpha1)*cos(lat1))
    #sigma02 = sigma01+sigma12
    lon01 = atan2(sin(alpha0)*sin(sigma01), cos(sigma01))
    lon0 = lon1 - lon01    
    npts = max(int((np.ceil(dist) + 1)/km), 100)
    all_d = np.linspace(0,dist,npts)/R  
    lons, lats = vlon_func(all_d, sigma01, alpha0, lon0), vlat_func(all_d, sigma01, alpha0, lon0)   
    
    return np.column_stack((lons, lats))

t_total0 = datetime.datetime.now()
number=0
while number < 50:
    t0 = datetime.datetime.now()

    #lat lon coordinates of random points generated within set shape file 
    coords = ps.points_in_shape(shape_path, N)
    #lons1,lats1 = coords[:,0], coords[:,1]
    #lons2,lats2 = lons1,lats1
    
    lonmin = np.floor(min(coords[:,0]))
    latmin = np.floor(min(coords[:,1]))

    coords1 = [coord1 for coord1 in coords for coord2 in coords]
                     
    coords2 = [coord2 for coord1 in coords for coord2 in coords]

    columns = np.column_stack((coords1, coords2))

    #dists = map(haversine, columns)
    
    #path_info = zip(coords1,coords2, dists)
    #path_info = np.column_stack((coords1, coords2, dists))
    #parallise the generation of path points for SPEED!
    
    #path_info = np.column_stack((coords1, coords2, dists))
    

    pool = mp.Pool()
    paths = pool.map(new_paths, columns)
    pool.close()
    pool.join()

    #create a flattened numpy array of size 2xN from the paths created! 
    paths = np.asarray(list(itertools.chain(*paths)))

    #keep all but the repeated coordinates by keeping only unique whole rows!
    #method is slowed without the b contiguous array
    b = np.ascontiguousarray(paths).view(np.dtype((np.void, paths.dtype.itemsize * paths.shape[1])))
    _, idx = np.unique(b, return_index=True)
    
    paths = np.unique(b).view(paths.dtype).reshape(-1, paths.shape[1])
     
     # Estimate the 2D histogram
   
    H, xedges, yedges = np.histogram2d(paths[:,0],paths[:,1],bins=nbins)
    
    #name = "path_density_2Dhist.png"
    
    GRAD = np.abs(np.asarray(np.gradient(H)[0]))#name = "path_density_2Dhist_grad.png"

    # H needs to be rotated and flipped
    H = np.rot90(H)
    GRAD = np.rot90(GRAD)

    H = np.flipud(H)
    GRAD = np.flipud(GRAD)

    # Mask zeros
    H = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero  
    GRAD = np.ma.masked_where(GRAD==0,GRAD) # Mask pixels with a value of zero  

    t1 = datetime.datetime.now()
    print t1-t0

    #fig = plt.figure()
    #plt.pcolormesh(xedges,yedges,H)
    #plt.xlabel('longitude (degrees)')
    #plt.ylabel('longitude (degrees)')
    #cbar = plt.colorbar()
    #cbar.ax.set_ylabel('Counts')
    #fig.savefig("SAVE_{}.png".format(number))

    number+=1
    
t_total1 = datetime.datetime.now()
print t_total1-t_total0