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
import random
import uuid





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

def cluster_points(coord_points, N):
    """
    Function that returns k which is an nx2 matrix of lon-lat vector columns
    containing the optimal cluster centroid spacings within a large set of random
    numbers e.g. those produced by the many_points() function above!
    """

    k = kmeans(coord_points, N)
    
    return k[0]

def paths_func(path_info, km=km_points):

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
    length = len(path)
    lons = [lon1 for i in range(0,length)]
    lats = [lat1 for i in range(0,length)]
    
    path = np.column_stack((path,lons,lats))
    
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


t0 = datetime.datetime.now()

ideal_path = 'ideal_coordinates.pickle'
#if no paths have been done before, start afresh!
if not os.path.exists(ideal_path):
    M = 1e5
    many_points = ps.points_in_shape(shape_path, M)
    coords = cluster_points(many_points,N)
    
#else import already processed coordinates if the program has already done so.
else:
    f = open(name=ideal_path, mode='rb')
    coords = pickle.load(f)
    f.close()

#generate N kmeans cluster points from massive M number of randomly distributed
#points inside the shape file. 

lonmin = np.floor(min(coords[:,0]))
latmin = np.floor(min(coords[:,1]))

#coords1 = [coord1 for coord1 in coords for coord2 in coords]
                     
#coords2 = [coord2 for coord1 in coords for coord2 in coords]

#columns = np.column_stack((coords1, coords2))

kappa = [np.vstack([[coord1[0],coord1[1],coord2[0],coord2[1]]\
                    for coord2 in coords]) for coord1 in coords]

def spread_paths(nets):

    #pool = mp.Pool()    
    #paths = pool.map(new_paths, nets)
    #pool.close()
    #pool.join()
    paths = map(paths_func, nets)
    #create a flattened numpy array of size 2xN from the paths created! 
    #paths = np.asarray(list(itertools.chain(*paths)))

    #keep all but the repeated coordinates by keeping only unique whole rows!
    #method is slowed without the b contiguous array
    #b = np.ascontiguousarray(paths).view(np.dtype((np.void, paths.dtype.itemsize * paths.shape[1])))
    #_, idx = np.unique(b, return_index=True)
    
    #paths = np.unique(b).view(paths.dtype).reshape(-1, paths.shape[1])

    #plt.figure()
    #plt.scatter(paths[:,0],paths[:,1])
    #name =  uuid.uuid4()
    #plt.savefig('{}.png'.format(name))
    
    return paths
    
t0 = datetime.datetime.now()
pool = mp.Pool()    
paths = pool.map(spread_paths, kappa)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print t1-t0
#paths = list(paths)

counter = 0

#cd Desktop/Link\ to\ SIMULATIONS/Network_Tracks/smarter_model/
grad_ideal, grad_check1, grad_check2, H_avg1, H_avg2 = 0, 0, 0, 0, 0
SHAPE = (1,1)

while counter <= 1000:
    
    t0 = datetime.datetime.now()
    
    #the following while loop is a work around fix to a:
    #new paths shape: (130, 100, 4) rather than being (130,) like it should be!
    while SHAPE != (130,):
        #remove a random set of paths associated with a single one of the N coordinates 
        old_path = paths.pop(random.randint(0,len(paths)-1))
    
        print "old path shape:", np.asarray(old_path).shape
    
        #figure out which old coordinate to remove from the coordinates list
        old_coord = [old_path[0][0][0],old_path[0][0][1]]
    
        print "old coord:", old_coord
        #NEED TO REMOVE OLD POINT FROM COORDS!
        #find index of array in nested array to remove!
        itemindex = np.where(coords==old_coord)[0][0]
        coords = list(coords)
        #find index of array in nested array to remove!
        del coords[itemindex]
        #coords = np.asarray(coords)
    
        #generate new point coordinate
        new_coord = ps.points_in_shape(shape_path, 1)[0]
    
    
        #place new coordinate in old set of coordinates
        coords = np.append(coords, [new_coord], axis=0)
    
    
    
        #generate new array of points in conjunction with the new randomly generated point!
        new_coord_set = np.vstack([[new_coord[0],new_coord[1],coord1[0],\
                                coord1[1]] for coord1 in coords])
    
        #generate new random point in place of all 'popped' points!
        new_paths = map(paths_func, new_coord_set)
        
        SHAPE = np.asarray(new_paths).shape
        print "new paths shape:", SHAPE
        
        #paths = np.asarray(paths)
        #if np.asarray(new_paths).shape != (130,):
        #    print("This one's trouble")
        #    print np.asarray(new_paths).shape
        #    new_paths =  np.asarray(new_paths[0]).reshape(130,)
    
    
    SHAPE = (1,1)

    #place new_paths in original path set!    
    paths = np.insert(paths, [1], [new_paths], axis=0)
    
    
    #paths = paths.append(new_paths)
    #paths = np.concatenate((paths, [new_paths]), axis=0)
    
    #paths = np.append(paths, new_paths, axis=0)
    
    #create a flattened numpy array of size 2xN from the paths created! 

    paths_density_check = list(itertools.chain(*paths))
    paths_density_check = np.asarray(list(itertools.chain(*paths_density_check)))
    
    #keep all but the repeated coordinates by keeping only unique whole rows!
    #method is slowed without the b contiguous array
    b = np.ascontiguousarray(paths_density_check).view(np.dtype\
    ((np.void, paths_density_check.dtype.itemsize * paths_density_check.shape[1])))
    _, idx = np.unique(b, return_index=True)
        
    paths_density_check = np.unique(b).view(paths_density_check.dtype)\
    .reshape(-1, paths_density_check.shape[1])

    #plt.figure()
    #plt.scatter(paths_density_check[:,0],paths_density_check[:,1])
    #plt.savefig('{}.png'.format(counter))
    
    paths = list(paths)
    
    # Estimate the 2D histogram
    
    H, xedges, yedges = np.histogram2d(paths_density_check[:,0],paths_density_check[:,1],bins=nbins)
        
    GRAD = np.abs(np.asarray(np.gradient(H)[0]))
    
    # H needs to be rotated and flipped
    H = np.rot90(H)
    GRAD = np.rot90(GRAD)
    
    H = np.flipud(H)
    GRAD = np.flipud(GRAD)
    
    # Mask zeros
    H = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero  
    GRAD = np.ma.masked_where(GRAD==0,GRAD) # Mask pixels with a value of zero
    
    print np.argmin(H)