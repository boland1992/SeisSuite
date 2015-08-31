# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 17:25:07 2015

@author: boland
"""
#-----------------------------------------------------------------------------
# IMPORT MODULES
#-----------------------------------------------------------------------------

import os
import fiona
#import pylab
import pickle
import pyproj
import random
import datetime
import itertools
import numpy as np
#import seaborn as sns
import pointshape as ps
import multiprocessing as mp
import matplotlib.pyplot as plt
from math import sqrt, radians, cos, sin, asin
from info_dataless import locs_from_dataless
from matplotlib.colors import LogNorm
from scipy.spatial import ConvexHull
#from scipy.spatial import Delaunay
#from scipy.cluster.vq import kmeans
from shapely import geometry

#-----------------------------------------------------------------------------
# CONVEX HULL FUNCTIONS
#-----------------------------------------------------------------------------
def convex_hull(points):
    """
    Function to produce a convex hull object that surrounds 
    a set of points. The input must be of the type Nx2 matrix/
    numpy array or equivalent. 
    """
    return ConvexHull(points)
    
def points_in_hull(new_point):
    """
    Function that checks if a singular point is contained within
    a convex hull object and returns True or False depending
    on the results. 
    """

    new_points = np.append(points, [new_point], axis=0)
    new_hull = ConvexHull(new_points)
    new_point = np.asarray(new_point)
    if list(hull.vertices) == list(new_hull.vertices):
        return new_point
    else:
        return np.zeros(new_point.shape)

def many_hull(points, new_points):
    """
    Function that quickly maps to the points_in_hull function
    and swiftly checks whether many points are contained within
    a convex hull object. 
    """
    points_contained = map(points_in_hull, new_points)
    return points_contained
    
def plot_hull(points, plot_points, hull, show_points=False):
    """
    Function that plots the boundaries of a convex hull using 
    matplotlib.pyplot. Input hull must be of type:
    scipy.spatial.qhull.ConvexHull 
    
    points input must be of the original coordinates.
    """
    plt.figure()
    for simplex in hull.simplices:
        plt.plot(points[simplex,0], \
        points[simplex,1], 'k-')
    if show_points:
        plt.scatter(plot_points[:,0], \
        plot_points[:,1], s=10,c='g')
        plt.scatter(points[:,0], \
        points[:,1], s=30,c='orange')
    plt.show()
    
def rand_hull(hull, points, K):
    "Generate K new random points contained within a convex hull"
    
    
    minx, maxx = np.min(points[:,0]), np.max(points[:,0])  
    miny, maxy = np.min(points[:,1]), np.max(points[:,1])
    X = abs(maxx - minx) * np.random.rand(10*K**2,1) + minx
    Y = abs(maxy - miny) * np.random.rand(10*K**2,1) + miny
    new_coords = np.column_stack((X,Y))
    
    pool = mp.Pool()    
    points_contained = pool.map(points_in_hull, new_coords)
    pool.close()
    pool.join()
    
    points_contained = np.asarray(points_contained)

    FIND_ZERO1 = np.where(points_contained[:,0]==0)[0]
    FIND_ZERO2 = np.where(points_contained[:,1]==0)[0]

    if len(FIND_ZERO1) != 0 and len(FIND_ZERO2) != 0:
        points_containedx = np.delete(points_contained[:,0], FIND_ZERO1)
        points_containedy = np.delete(points_contained[:,1], FIND_ZERO2)
    
    return np.column_stack((points_containedx,points_containedy))[:K]
        
    

#-----------------------------------------------------------------------------
# GENERATE FIRST SET OF VARIABLES AND STATES
#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
#SHAPEFILE FUNCTIONS
#-----------------------------------------------------------------------------
def shape_(input_shape):
    with fiona.open(input_shape) as fiona_collection:
        # In this case, we'll assume the shapefile only has one record/layer 
        shapefile_record = fiona_collection.next()
        # Use Shapely to create the polygon
        return geometry.asShape( shapefile_record['geometry'] )

shape = shape_(shape_path)

def point_check(coord, shape=shape):
    point = geometry.Point(coord[0], coord[1]) 
    if shape.contains(point):
        return coord   

def many_points(coords, shape):
    """
    Funtion that returns lat-lon coodinates of many random points within 
    a shapefile shape e.g. boundaries of a state or country.
    """
    
    minx, miny, maxx, maxy = shape.bounds
    X, Y = coords[:,0], coords[:,1]
    coords = np.column_stack((X,Y))
    coord_points = []
    #generate points in parallel for speed, order is not preserved.
    pool = mp.Pool()
    coord_points = pool.map(point_check, coords)
    pool.close()
    pool.join()
    #remove None values from above numpy array! 
    X = [i[0] for i in coord_points if i != None]
    Y = [i[1] for i in coord_points if i != None]    
    #convert np array to kmeans function friendly nx2 matrix!         
    #output is a nx2 vector stacked matrix called coord_points, 
    #coord_points[:,0] are lons and coord_points[:,1] are lats
    return np.column_stack((X,Y))

        
#-----------------------------------------------------------------------------
#PATHS AND DISTANCES ON GREAT-CIRCLE FUNCTIONS
#-----------------------------------------------------------------------------
def remove_distance(period_range, max_dist = 2000):
    """
    Function that returns a given possible resolvable ambient noise
    tomography structure distance range, given the maximum period range
    availabe to the study. The distance returned is in km.
    
    Maximum distance default can be reassigned based on the cut-off found
    by your time-lag plots for your study!
    """
    if type(period_range) == list:
        min_dist = min(period_range) * 9
        return [min_dist, max_dist]
        
    elif type(period_range) == int or float:
        return [period_range*9, max_dist]

global dist_range
dist_range = remove_distance(period_range, max_dist = 2000)

def haversine2(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees).
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon, dlat = lon2 - lon1, lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km

def fast_geodesic(lon1,lat1,lon2,lat2,npts):
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

def paths_func(path_info, km=km_points):

    lon1, lat1, lon2, lat2 = path_info[0], \
    path_info[1], path_info[2], path_info[3]    
    # interpoint distance <= 1 km, and nb of points >= 100
    dist = haversine2(lon1, lat1, lon2, lat2)    
    npts = max(int((np.ceil(dist) + 1)/km), 100)    
    path = fast_geodesic(lon1,lat1,lon2,lat2, npts)
    #lons, lats = [lon1 for i in range(0,len(path))], \
    #[lat1 for i in range(0,len(path))]    
    #path = np.column_stack((path,lons,lats))
    
    if min(dist_range) < dist < max(dist_range):
        #remove the closest points along this line that fall below the distance 
        #find the index of the first point that is above this distance away!
        pts_km = npts / float((np.ceil(dist) + 1)) #this gives pts/km
        #remove all points below this index in the paths list
        dist_index = pts_km * min(dist_range) 
        path = path[dist_index:]
        return path    
    else:
        return np.zeros_like(path)
        
def points_and_distances(path_info):
    lon1, lat1, lon2, lat2 = path_info[0], \
    path_info[1], path_info[2], path_info[3]
    dist = haversine2(lon1, lat1, lon2, lat2)
    return [dist, lon2, lat2]
    
def fast_dists(nets):
    dist_points = map(points_and_distances, nets)
    return np.vstack(dist_points)
    
def spread_paths(nets):
    return map(paths_func, nets)
    
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

if not shape_boundary:
# produce convex hull around initial station spacings if no shape file
    #CONVEX HULL
    global points
    points = coords
    hull = convex_hull(coords)
    minx, maxx = np.min(coords[:,0]), np.max(coords[:,0])  
    miny, maxy = np.min(coords[:,1]), np.max(coords[:,1])
    X = abs(maxx - minx) * np.random.rand(1e4,1) + minx
    Y = abs(maxy - miny) * np.random.rand(1e4,1) + miny
    new_coords = np.column_stack((X,Y))
    
    pool = mp.Pool()    
    points_contained = pool.map(points_in_hull, new_coords)
    pool.close()
    pool.join()
    
    points_contained= np.asarray(points_contained)

    FIND_ZERO1 = np.where(points_contained[:,0]==0)[0]
    FIND_ZERO2 = np.where(points_contained[:,1]==0)[0]

    if len(FIND_ZERO1) != 0 and len(FIND_ZERO2) != 0:
        points_containedx = np.delete(points_contained[:,0], FIND_ZERO1)
        points_containedy = np.delete(points_contained[:,1], FIND_ZERO2)
    
    points_contained = np.column_stack((points_containedx,points_containedy))
    
    plot_hull(coords, points_contained, hull, show_points=True)
    
lonmin, lonmax = np.floor(min(coords[:,0])), np.ceil(max(coords[:,0]))
latmin, latmax = np.floor(min(coords[:,1])), np.ceil(max(coords[:,1]))
print lonmin,lonmax,latmin,latmax

kappa = [np.vstack([[coord1[0],coord1[1],coord2[0],coord2[1]]\
                    for coord2 in coords]) for coord1 in coords]
                        
t0 = datetime.datetime.now()
pool = mp.Pool()    
paths = pool.map(spread_paths, kappa)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print t1-t0


counter, counter2 = 0, 0
#cd Desktop/Link\ to\ SIMULATIONS/Network_Tracks/smarter_model/
grad_ideal, grad_check1, grad_check2, H_avg1, H_avg2 = 0, 0, 0, 0, 0
SHAPE = (1,1)
perc_high = 0.01
low_counter = 0
random_counter = 0
#new_coord = 0
infinite_counter = 0
find_it = []
check_coord = None
use_old_path = False
searches_per_point = 3
print len(paths)

while infinite_counter <= 1:
    t0 = datetime.datetime.now()
    
    #----------------------------------------------------------------------
    # Generate N new point coordinates
    #----------------------------------------------------------------------
    if shape_boundary:
        new_coords = ps.points_in_shape(shape_path, N)
    
    else:
        new_coords = rand_hull(hull, original_coords, N)
        
    
    #new_coords = np.asarray([[128.,-28.]])
    # Place new coordinates in old coordinate set
    coords = np.append(coords, new_coords, axis=0)
    # Generate new array of coordinate combinations for new paths.    
    coord_set = [np.vstack([[coord1[0],coord1[1],coord2[0],coord2[1]]\
                 for coord2 in coords]) for coord1 in coords]

    t0 = datetime.datetime.now()
    pool = mp.Pool()    
    paths = pool.map(spread_paths, coord_set)
    pool.close()
    pool.join()
    t1 = datetime.datetime.now()
    print "time to generate new paths", t1-t0
    
    # Append new set of paths now that old set has been deleted.
    print len(paths)
    
    #create a flattened numpy array of size 2xN from the paths created! 
    paths_density_check = list(itertools.chain(*paths))
    paths_density_check = np.asarray(list(itertools.chain\
                                    (*paths_density_check)))
    
    #keep all but the repeated coordinates by keeping only unique whole rows!
    b = np.ascontiguousarray(paths_density_check).view(np.dtype\
    ((np.void, paths_density_check.dtype.itemsize * \
    paths_density_check.shape[1])))
    _, idx = np.unique(b, return_index=True)
        
    paths_density_check = np.unique(b).view(paths_density_check.dtype)\
    .reshape(-1, paths_density_check.shape[1])
    
    #remove all path points that lay outside the shape file polygon
    #paths_density_check = ps.paths_in_shape(paths_density_check)
                 
    # Reset paths as a list to be able to delete indices on next loop.
    paths = list(paths)

    #-------------------------------------------------------------------------
    # Remove zeroes from path_density_check to ensure all paths that
    # were NOT in the distance threshold are removed from the path density
    # calculation! 
    #-------------------------------------------------------------------------

    path_density_lons, path_density_lats = paths_density_check[:,0], \
                                           paths_density_check[:,1]

    FIND_ZERO1 = np.where(paths_density_check[:,0]==0)[0]
    FIND_ZERO2 = np.where(paths_density_check[:,1]==0)[0]
    if len(FIND_ZERO1) != 0 and len(FIND_ZERO2) != 0:
        path_density_lons = np.delete(path_density_lons, FIND_ZERO1)
        path_density_lats = np.delete(path_density_lats, FIND_ZERO2)
        
    #-------------------------------------------------------------------------
    # Set up path density calculations using numpy's histogram2d function.
    #-------------------------------------------------------------------------
    H, xedges, yedges = np.histogram2d(path_density_lons,
                                       path_density_lats,
                                       bins=nbins)
    # Calculate the gradient field of the path density field.
    GRAD = np.abs(np.asarray(np.gradient(H)[0]))
    
    # H needs to be rotated and flipped.
    H = np.rot90(H)
    GRAD = np.rot90(GRAD)
    
    H = np.flipud(H)
    GRAD = np.flipud(GRAD)
    
    # Mask zeros
    H = np.ma.masked_where(H==0,H)  
    GRAD = np.ma.masked_where(GRAD==0,GRAD)

    H_avg1 = np.average(H)
    grad_check1 = np.std(GRAD)
    
    search = np.where(H > np.average(H))#+ 0.20 * (np.max(H) - np.average(H)))
    # Scale these points with respect to the lat-lon limits!
    Hmaxx, Hmaxy =  search[1], search[0]    
    Hmaxx = (lonmax-lonmin)/(nbins) * Hmaxx + lonmin
    Hmaxy = (latmax-latmin)/(nbins) * Hmaxy + latmin
    # Make sure all low density coordinates ARE within shapefile!
    highest_density_coords = ps.paths_in_shape(np.column_stack((Hmaxx, Hmaxy)))
    
    #-------------------------------------------------------------------------
    # Find coordinates of high density within path density field.
    #-------------------------------------------------------------------------
    
    #CONVEX HULL
    # Test:
    #p = np.asarray([[128.,-28.]])
    #is_p_inside_points_hull(list(highest_density_coords), p)

    # Plot:
    #for simplex in hull.simplices:
    #    pylab.plot(highest_density_coords[simplex,0], highest_density_coords[simplex,1], 'k-')
    #pylab.plot(p[:,0], p[:,1], '^r')
    #pylab.show()
    
    #DELAUNAY TRIANGLES
    #tri = Delaunay(list(highest_density_coords))

    #plt.triplot(highest_density_coords[:,0],\
    #            highest_density_coords[:,1], tri.simplices.copy())
    #plt.plot(new_coords[:,0], new_coords[:,1], 'o')
    #plt.show()
    #print tri.find_simplex(new_coords[0])>=0


    if counter == 0:
        grad_ideal = 1e6
        avg_ideal = 0  
    
    if grad_check1 < grad_ideal and avg_ideal < H_avg1:     

        with open(u'ideal_coordinates.pickle', 'wb') as f:
            print "\nExporting new ideal coordinates."
            pickle.dump(coords, f, protocol=2)
            
        fig2 = plt.figure()
        plt.pcolormesh(xedges,yedges,H)
        
        plt.pcolor(xedges, yedges, H, norm=LogNorm(vmin=np.min(H), \
        vmax=np.max(H)), cmap='PuBu')
        plt.xlabel('longitude (degrees)')
        plt.ylabel('latitude (degrees)')
        plt.xlim(lonmin-0.05*abs(lonmax-lonmin), lonmax+0.05*abs(lonmax-lonmin))
        plt.ylim(latmin-0.05*abs(latmax-latmin),latmax+0.05*abs(latmax-latmin))
        col = plt.colorbar()
        col.ax.set_ylabel('Counts')
        #plt.scatter(highest_density_coords[:,0], \
        #highest_density_coords[:,1],c='g', s=30)
        
        plt.scatter(new_coords[:,0], new_coords[:,1],c='r', s=30)
        fig2.savefig("new_coords{}.png".format(counter))
        fig2.clf()
                
        grad_ideal = grad_check1
        avg_ideal = H_avg1
        #remove new coordinate and replace with old coordinate
        coords = list(coords)
        del coords[-N:]
        coords = np.asarray(coords)
        #remove new path and replace it with the old set!
        paths = list(paths)
        del paths[-N:]

    else:
        #RESET values!
        #remove new coordinate and replace with old coordinate
        coords = list(coords)
        del coords[-N:]
        coords = np.asarray(coords)
        #place new coordinate in old set of coordinates
        #remove new path and replace it with the old set!
        paths = list(paths)
        del paths[-N:]
        
    #----------------------------------------------------------------------
    # Delete old paths points, if the program has already run and not found
    # good new solutions.
    #----------------------------------------------------------------------
    
    find_it = []
    counter+=1
    print "counter2:",counter2

    t1 = datetime.datetime.now()
    print "That loop took: ", t1-t0