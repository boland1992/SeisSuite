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
import pickle
import pyproj
import random
import datetime
import itertools
import numpy as np
import pointshape as ps
import multiprocessing as mp
import matplotlib.pyplot as plt
from math import sqrt, radians, cos, sin, asin
from scipy.cluster.vq import kmeans
from shapely import geometry

#-----------------------------------------------------------------------------
# GENERATE SECOND SET OF VARIABLES AND STATES
#-----------------------------------------------------------------------------

verbose = False
#Enter path to boundary shape file.
shape_path = "/home/boland/Dropbox/University/UniMelb\
/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"
# Enter number of stations required.
N = 130
# Enter km spacing between path density points.
km_points = 100.0
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

def many_points(shape, coords):
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

def cluster_points(coord_points, N):
    """
    Function that returns k which is an nx2 matrix of lon-lat vector columns
    containing the optimal cluster centroid spacings within a large set of random
    numbers e.g. those produced by the many_points() function above!
    """
    return kmeans(coord_points, N)[0]
    
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
if not os.path.exists(ideal_path):
    #generate N kmeans cluster points from massive M number of randomly 
    #distributed points inside the shape file. 
    M = 1e5
    lot_points = ps.points_in_shape(shape_path, M)
    coords = cluster_points(lot_points, N)
#else import already processed coordinates if the program has already done so.
else:
    f = open(name=ideal_path, mode='rb')
    coords = pickle.load(f)
    f.close()

lonmin, lonmax = np.floor(min(coords[:,0])), np.ceil(max(coords[:,0]))
latmin, latmax = np.floor(min(coords[:,1])), np.ceil(max(coords[:,1]))
print lonmin,lonmax,latmin,latmax

t0 = datetime.datetime.now()


kappa = [np.vstack([[coord1[0],coord1[1],coord2[0],coord2[1]]\
                    for coord2 in coords]) for coord1 in coords]
                        
t1 = datetime.datetime.now()
dif1 = t1-t0
print dif1
t0 = datetime.datetime.now()

#kappa = np.asarray(list(itertools.combinations(coords, 2)))

#no_rows, no_cols = kappa.shape[0], kappa.shape[1]*kappa.shape[2]
#kappa = kappa.reshape((no_rows,no_cols))
#t1 = datetime.datetime.now()

#kappa = np.asarray([[[i[0][1],i[0][1],i[1][0],i[1][1]] for i in kappa]])

pool = mp.Pool()    
paths = pool.map(spread_paths, kappa)
pool.close()
pool.join()

quit()

counter, counter2 = 0, 0
#cd Desktop/Link\ to\ SIMULATIONS/Network_Tracks/smarter_model/
grad_ideal, grad_check1, grad_check2, H_avg1, H_avg2 = 0, 0, 0, 0, 0
SHAPE = (1,1)
perc_high = 0.01
low_counter = 0
random_counter = 0
new_coord = 0
infinite_counter = 0
find_it = []
check_coord = None
use_old_path = False
searches_per_point = 3

while infinite_counter <= 1:
    
    rand_indicator = random.randint(1,10)

    if verbose:
        print "check coord", check_coord
    t0 = datetime.datetime.now()
    
    print use_old_path
    #the following while loop is a work around.
    #new paths shape: (130, 100, 4) rather than being (130,) 
    while SHAPE != (130,):
        
        if check_coord is None or rand_indicator < 4:
            #------------------------------------------------------------------
            # Option one: randomly select old coordinate to move.
            #------------------------------------------------------------------
            while len(find_it) == 0:
                # Remove random set of paths associated with the N coordinates. 
                rand_int = random.randint(0,len(paths)-1)
                old_path = paths[rand_int]
                # Determine which old coordinate to remove.
                old_coord = coords[rand_int]#[old_path[0][0][0],old_path[0][0][1]]
                # Find index of array in nested array to remove!
                find_it = np.where(coords==old_coord)[0]
            coord_index = find_it[0]
            counter2 = 0
            print "a1"

        
        elif counter2 == searches_per_point-1 or not use_old_path:
        #------------------------------------------------------------------
        # Option two: select new high density point if too many searches per
        # point OR it is stated that use_old_path = False 
        #------------------------------------------------------------------            
            print counter2 % searches_per_point            
            while len(find_it) == 0:
                old_coord = check_coord
                find_it = np.where(coords==old_coord)[0]

                for paths_index, path in enumerate(paths):
                    for dots in path:
                        find_paths = np.where(dots==old_coord[0])[0]
                        if len(find_paths) != 0:
                            rand_int = False
                            break
            coord_index = find_it[0]
            counter2 = 0
            print "a2"
            
        elif use_old_path:
        #------------------------------------------------------------------
        # Option three: if conditions not met, move same old point.
        #------------------------------------------------------------------
            coord_index = -1
            old_path = paths[coord_index]
            old_coord = coords[coord_index]
            counter2 += 1
            print "a3"
                
        if verbose:
            print "coord index to remove is:", coord_index
        coords = list(coords)
        del coords[coord_index]
        coords = np.asarray(coords)
        new_coord_first = new_coord

        #----------------------------------------------------------------------
        # Generate new point coordinate.
        #----------------------------------------------------------------------
        if counter <= 1:
            #------------------------------------------------------------------
            # Option one: generate random new coordinate within boundary shape.
            #------------------------------------------------------------------
            new_coord = ps.points_in_shape(shape_path, 1)[0]
        else:
            #------------------------------------------------------------------
            # Option two: go with whatever previous calculation had been made.
            #------------------------------------------------------------------
            new_coord = new_coord
        
        # Place new coordinate in old set of coordinates
        coords = np.append(coords, [new_coord], axis=0)
        
        # Generate new array of coordinate combinations for new paths.
        new_coord_set = np.vstack([[new_coord[0],new_coord[1],coord1[0],\
                                  coord1[1]] for coord1 in coords])
        # Generate new path points.
        new_paths = map(paths_func, new_coord_set)
        SHAPE = np.asarray(new_paths).shape
        
        if not SHAPE == (130,):
            #remove substitute back the old coordinate for the new coordinate!
            coords = list(coords)
            #find index of array in nested array to remove!
            del coords[-1]
            coords = np.asarray(coords)
            #place new coordinate in old set of coordinates
            coords = np.append(coords, [old_coord], axis=0)
    
    


    #----------------------------------------------------------------------
    # Delete old paths points.
    #----------------------------------------------------------------------
    if not rand_int:
        if verbose:
            print "path index to remove is:", coord_index
        del paths[paths_index]
    elif use_old_path:
        del paths[-1]
    else:
        if verbose:
            print "paths index to remove is:", coord_index
        del paths[rand_int]
    

    # Reset shape to work around error from above. 
    SHAPE = (1,1)
    
    # Append new set of paths now that old set has been deleted.
    paths = np.append(paths, [new_paths], axis=0)

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
    #GRAD = np.rot90(GRAD)
    
    H = np.flipud(H)
    #GRAD = np.flipud(GRAD)
    
    # Mask zeros
    H = np.ma.masked_where(H==0,H)  
    #GRAD = np.ma.masked_where(GRAD==0,GRAD)

    H_avg1 = np.average(H)
    grad_check1 = np.std(GRAD)

    #-------------------------------------------------------------------------
    # Find coordinates of high density within path density field.
    #-------------------------------------------------------------------------
    
    search = np.where(H > 0.50 * np.max(H))
    # Scale these points with respect to the lat-lon limits!
    Hmaxx, Hmaxy =  search[1], search[0]    
    Hmaxx = (lonmax-lonmin)/(nbins) * Hmaxx + lonmin
    Hmaxy = (latmax-latmin)/(nbins) * Hmaxy + latmin
    # Make sure all low density coordinates ARE within shapefile!
    highest_density_coords = ps.paths_in_shape(np.column_stack((Hmaxx, Hmaxy)))
    
    WHERE = np.where(H < perc_high*H_avg1)
    if 0 < rand_indicator <= 3: 
    #half the time move the coordinates to low density locations.

        #scale these points with respect to the lat-lon limits!
    
        Hminx, Hminy =  WHERE[1], WHERE[0]    
        Hminx = (lonmax-lonmin)/(nbins) * Hminx + lonmin
        Hminy = (latmax-latmin)/(nbins) * Hminy + latmin

        # Make sure all low density coordinates ARE within shapefile!
        low_density_coords = ps.paths_in_shape(np.column_stack((Hminx, Hminy)))
        
        # Low_density_coords = many_points(shape, np.column_stack((Hminx, Hminy)))
        if len(low_density_coords) == 0:
            new_coord = ps.points_in_shape(shape_path, 1)[0]
            #increase percentage of search if no new low density points are created!
            perc_high +=0.05
        elif len(low_density_coords) == 1:
            new_coord = low_density_coords[0]
            perc_high +=0.05
        else:
            new_coord = low_density_coords[random.randint(0,len(low_density_coords)-1)]

    elif 3 < rand_indicator <= 10: 
    # Half the time move coordinates to random locations.
        new_coord = ps.points_in_shape(shape_path, 1)[0]

    if counter == 0:
        grad_ideal = 1e6
        avg_ideal = 0
    
    #fig5 = plt.figure()
    #plt.scatter(coords[:,0], coords[:,1],c='b', s=10)
    #plt.scatter(new_coord[0], new_coord[1],c='r', s=30)
    #plt.scatter(old_coord[0], old_coord[1],c='g', s=30)
    #fig5.savefig("coords{}.png".format(counter))
    #fig5.clf()    
    
    if grad_check1 < grad_ideal and avg_ideal < H_avg1:     

        with open(u'ideal_coordinates.pickle', 'wb') as f:
            print "\nExporting new ideal coordinates."
            pickle.dump(coords, f, protocol=2)
            
        fig2 = plt.figure()
        plt.pcolormesh(xedges,yedges,H)
        plt.xlabel('longitude (degrees)')
        plt.ylabel('latitude (degrees)')
        plt.xlim(lonmin-0.05*abs(lonmax-lonmin), lonmax+0.05*abs(lonmax-lonmin))
        plt.ylim(latmin-0.05*abs(latmax-latmin),latmax+0.05*abs(latmax-latmin))
        col = plt.colorbar()
        col.ax.set_ylabel('Counts')
            
        #plt.scatter(highest_density_coords[:,0],\
        #highest_density_coords[:,1],c='orange', s=10)
        plt.scatter(new_coord[0], new_coord[1],c='r', s=30)
        plt.scatter(old_coord[0], old_coord[1],c='g', s=30)
        fig2.savefig("min_density{}.png".format(counter))
        fig2.clf()
        
        # Assign new values.
        use_old_path = False
        grad_ideal = grad_check1
        avg_ideal = H_avg1

    else:
        #RESET values!
        #remove new coordinate and replace with old coordinate
        coords = list(coords)
        del coords[-1]
        coords = np.asarray(coords)
        #place new coordinate in old set of coordinates
        coords = np.append(coords, [old_coord], axis=0)
        #remove new path and replace it with the old set!
        paths = list(paths)
        del paths[-1]
        paths = list(np.append(paths, [old_path], axis=0))
        use_old_path = True
        
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    if not use_old_path:
        # Calculate the combinations of coordinates from real points, to high points! 
        high_point_coord_combs = [np.vstack([[coord1[0],coord1[1],coord2[0],\
        coord2[1]] for coord1 in highest_density_coords for coord2 in coords])]
        
        #print "high point coordinate combitions:", high_point_coord_combs
        pool = mp.Pool()    
        high_point_coords = pool.map(fast_dists, high_point_coord_combs)
        pool.close()
        pool.join()
        
        high_point_coords =  np.vstack(np.asarray(\
        list(itertools.chain(*high_point_coords))))
        high_point_coords = high_point_coords[high_point_coords[:,0].argsort()]
#       t1 = datetime.datetime.now()
#        print "fancy index sorting time: ", t1-t0
        # SELECT 20% OF THESE LOWEST DISTANCE COORDINATES FROM HIGH DENSITY POINTS! 
        point_index = int(0.2*len(high_point_coords))
        high_point_coords = high_point_coords[:point_index]
        # Create unique list of coordinates high_point_coords and remove distances.
        high_point_coords =  np.column_stack((high_point_coords[:,1], 
                                              high_point_coords[:,2]))
        b = np.ascontiguousarray(high_point_coords).view(np.dtype\
        ((np.void, high_point_coords.dtype.itemsize * \
        high_point_coords.shape[1])))
        _, idx = np.unique(b, return_index=True)
        
        high_point_coords = np.unique(b).view(high_point_coords.dtype)\
        .reshape(-1, high_point_coords.shape[1])
        # Find random high density coord. This is the old coordinate to remove.
        check_coord = high_point_coords[random.randint(0,len(high_point_coords)-1)]
        
    if check_coord in coords:
        old_coord = check_coord
    elif not use_old_path:
        check_coord = None
    

    find_it = []
    counter+=1
    print "counter2:",counter2

    t1 = datetime.datetime.now()
    print "That loop took: ", t1-t0