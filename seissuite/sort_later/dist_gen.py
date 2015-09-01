# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 09:12:49 2015

@author: boland
"""

import sys
sys.path.append("/home/boland/Anaconda/lib/python2.7/site-packages")
import pickle
import fiona
import seaborn as sns
from shapely import geometry
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans
import multiprocessing as mp
import pyproj
from matplotlib.colors import ColorConverter
import os
from matplotlib.colors import LinearSegmentedColormap
import itertools
import datetime
import pointshape as ps
from numba import jit
from math import sqrt, atan2, radians,degrees, cos, tan, sin, asin, acos
shape_path = "/home/boland/Dropbox/University/UniMelb/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"

N = 130
#enter km spacing between path density points
km_points = 10.0

# reference elipsoid to calculate distance
wgs84 = pyproj.Geod(ellps='WGS84')
#wgs84 = pyproj.Geod(ellps='sphere')

#lat lon coordinates of random points generated within set shape file 
coords = ps.points_in_shape(shape_path, N)


def dist(lons1, lats1, lons2, lats2):
    """
    Returns an array of geodetic distance(s) in km between
    points (lon1, lat1) and (lon2, lat2)
    """
    _, _, d = wgs84.inv(lons1=lons1, lats1=lats1, lons2=lons2, lats2=lats2)
    return np.array(d) / 1000.0
    
    
def new_dist(coordinates):
    """
    Returns an array of geodetic distance(s) in km between
    points (lon1, lat1) and (lon2, lat2)
    """
    lons1=coordinates[0]
    lats1=coordinates[1]
    lons2=coordinates[2]
    lats2=coordinates[3]
    
    _, _, d = wgs84.inv(lons1=lons1, lats1=lats1, lons2=lons2, lats2=lats2)
    return np.array(d) / 1000.0
    
    
    
from math import radians, cos, sin, asin, sqrt
def haversine(coordinates):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    lon1=coordinates[0]
    lat1=coordinates[1]
    lon2=coordinates[2]
    lat2=coordinates[3]
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
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


lons1 = coords[:,0]; lats1 = coords[:,1] 
lons2 = coords[:,0]; lats2 = coords[:,1]
    
lonmin = np.floor(min(lons1))
latmin = np.floor(min(lats1))

dists = []
coords1 = []
coords2 = []



t0 = datetime.datetime.now()
coords11 = [coord1 for coord1 in coords for coord2 in coords]
                 
coords22 = [coord2 for coord1 in coords for coord2 in coords]

dist_list = [dist(coord1[0], coord1[1],coord2[0],coord2[1]) \
             for coord1 in coords for coord2 in coords]

t1 = datetime.datetime.now()

print "dist list comprehension", t1-t0





t0 = datetime.datetime.now()
coords11 = [coord1 for coord1 in coords for coord2 in coords]
                 
coords22 = [coord2 for coord1 in coords for coord2 in coords]

columns = np.column_stack((coords11, coords22))

t0 = datetime.datetime.now()

dist_list = map(new_dist, columns)


t1 = datetime.datetime.now()


print "map dists", t1-t0

t0 = datetime.datetime.now()
coords11 = [coord1 for coord1 in coords for coord2 in coords]
                 
coords22 = [coord2 for coord1 in coords for coord2 in coords]

columns = np.column_stack((coords11, coords22))

t0 = datetime.datetime.now()

dist_list = map(haversine, columns)


t1 = datetime.datetime.now()


print "haversine", t1-t0

@jit
def distances(lats1, lats2):
    for index, item in enumerate(lats1):
        lat1 = item
        lon1 = lons1[index]
        
        for index2, item2 in enumerate(lats2):
            lat2 = item2
            lon2 = lons2[index2]
        
            dists.append(dist(lon1, lat1, lon2, lat2))
            coords1.append([lon1,lat1])
            coords2.append([lon2,lat2])
            
    return dists, coords1, coords2

t0 = datetime.datetime.now()

dists, coords1, coords2 = distances(lats1, lats2)

t1 = datetime.datetime.now()
print t1-t0


t0 = datetime.datetime.now()
path_info = zip(coords1,coords2, dists)
#path_info = np.column_stack((coords1, coords2, dists))
#del dists; del coords1; del coords2

def create_paths(path_point, km=km_points):
    coord1 = path_point[0]
    coord2 = path_point[1]
    dist = path_point[2]
    # interpoint distance <= 1 km, and nb of points >= 100
    npts = max(int((np.ceil(dist) + 1)/km), 100)
    path = geodesic(coord1, coord2, npts)
    #print("still going strong\n")
    return path

#parallise the generation of path points for SPEED!
pool = mp.Pool()
paths = pool.map(create_paths, path_info)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print 'paths parallel create_paths', t1-t0
del path_info

def new_paths(path_info, km=km_points):
 
    lon1, lat1, lon2, lat2, dist = path_info[0], path_info[1], path_info[2], path_info[3], path_info[4]  
  
    # interpoint distance <= 1 km, and nb of points >= 100
    npts = max(int((np.ceil(dist) + 1)/km), 100)

    path = new_geodesic(lon1,lat1,lon2,lat2, npts)
    #print("still going strong\n")
    return path



#parallise the generation of path points for SPEED!
t0 = datetime.datetime.now()

path_info = np.column_stack((coords1, coords2, dists))

pool = mp.Pool()
paths = pool.map(new_paths, path_info)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print 'paths parallel new_paths', t1-t0
del path_info


#plt.figure()
#plt.scatter(paths[1][:,0], paths[1][:,1])
#plt.show()



def points_curve(path_info, km=km_points):
    """
    Function that returns N equidistance lat lon coordinates
    along the great-circle line between the two lat lon coordinates.
    """
    lon1 = path_info[0]   
    lat1 = path_info[1]   
    lon2 = path_info[2]   
    lat2 = path_info[3]   
    dist = path_info[4]
        
    lat_dif = lat2-lat1
    lon_dif = lon2-lon1
    
    npts = max(int((np.ceil(dist) + 1)/km), 100)
    lat_op = (lat_dif)/npts
    lon_op = (lon_dif)/npts

    nums = np.arange(1,npts+1,1)

    latn = np.add(lat1, np.multiply(lat_op,nums))
    lonn = np.add(lon1, np.multiply(lon_op,nums))

    
    #latn = lat1 + lat_op * nums#[lat1 + n * lat_op for n in nums] 
    #lonn = lon1 + lon_op * nums #[lon1 + n * lon_op for n in nums]
    
    path = np.column_stack((lonn,latn))
    
    

    return path


t0 = datetime.datetime.now()

path_info = np.column_stack((coords1, coords2, dists))
pool = mp.Pool()
paths = pool.map(points_curve, path_info)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print 'paths parallel points_curve', t1-t0

del path_info

#plt.figure()
#plt.scatter(paths[1][:,0], paths[1][:,1])
#plt.show()
    
#from geographiclib.geodesic import Geodesic

#number_points = 10

#gd = Geodesic.WGS84.Inverse(35, 0, 35, 90)
#line = Geodesic.WGS84.Line(gd['lat1'], gd['lon1'], gd['azi1'])

#for i in range(number_points + 1):
#    point = line.Position(gd['s12'] / number_points * i)
#    print((point['lat2'], point['lon2']))



#t0 = datetime.datetime.now()

#paths_info = np.column_stack((coords1, coords2, dists))

#paths = map(new_paths, paths_info)

#t1 = datetime.datetime.now()
#print 'paths map', t1-t0

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

 
t0 = datetime.datetime.now()
path_info = np.column_stack((coords1, coords2, dists))
pool = mp.Pool()
paths = pool.map(waypoint_init, path_info)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print "waypoints", t1-t0

#plt.figure()
#plt.scatter(coords[:,0], coords[:,1])
#plt.show()
#parallise the generation of path points for SPEED!





    